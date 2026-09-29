// Exercise the production RT driver with known, normalized elastic kernels.
// These are transport/dispatch tests, not a Compton-quadrature accuracy test.
#include "compton_rt.h"
#include "save_results.h"
#include "constants.h"
#include "bezier3_transfer.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <vector>

static void check(bool ok, const char* message)
{
	if (!ok) throw std::runtime_error(message);
}

template<class Cache>
static void elastic_storage(Cache& c, const RTGrids& g, int angular_pairs)
{
	c.NT=1;c.NE=g.NE;
	c.T_grid=new double[1]{1e6};c.theta=new double[1]{1e-3};
	c.x_grid=new double[g.NE];c.glo=new int[g.NE];c.ghi=new int[g.NE];
	c.data_size=long(g.NE)*angular_pairs;
	c.data=new double[c.data_size];c.band_lo=new int[c.data_size];
	c.band_hi=new int[c.data_size];c.band_off=new long[c.data_size];
	for(int e=0;e<g.NE;++e) {
		c.x_grid[e]=g.ene[e]/phys::m_e_c2_eV;c.glo[e]=c.ghi[e]=e;
		for(int a=0;a<angular_pairs;++a) {
			long k=long(e)*angular_pairs+a;
			c.band_lo[k]=c.band_hi[k]=e;c.band_off[k]=k;
		}
	}
}

static std::vector<double> faces(const RadField& r, const RTGrids& g)
{
	std::vector<double> values;
	for(auto f : {r.Inu_top,r.Inu_bottom}) for(int m=0;m<g.NA;++m) for(int e=0;e<g.NE;++e)
		values.push_back(f[m][e]);
	return values;
}

static void verify_budget_and_moments(RadField& r, const RTGrids& g, int incidence)
{
	std::vector<double> previous;
	for(int d=0;d<g.ND_MID;++d) for(int e=0;e<g.NE;++e) previous.push_back(r.J0[d][e]);
	r.compute_moments();
	for(int d=0;d<g.ND_MID;++d) for(int e=0;e<g.NE;++e)
		check(std::abs(r.J0[d][e]-previous[d*g.NE+e])<1e-13*std::max(1.0,r.J0[d][e]),
		      "stored mean intensity does not match center moments");
	// Physical absorption needs the integrated intensity, not its center value.
	// Reconstruct the converged isotropic elastic source for a diagnostic sweep.
	const long size=long(g.ND_MID)*g.NA*g.NE;
	std::vector<double> dt(size), source(size), I(size), mean(size);
	std::vector<double> top(g.NE), bottom(g.NE), ft(g.NA*g.NE), fb(g.NA*g.NE);
	compute_boundary_illumination(g.NE,incidence,g,r,top.data(),bottom.data());
	for(int d=0;d<g.ND_MID;++d) for(int m=0;m<g.NA;++m) for(int e=0;e<g.NE;++e) {
		long k=(long(d)*g.NA+m)*g.NE+e;
		double chi=r.kabs[d][e]+r.ksct[d][e];
		dt[k]=chi*g.dr[d]/std::abs(g.mu[m]);
		source[k]=(r.jnu[d][e]+r.ksct[d][e]*r.J0[d][e])/chi;
	}
	dao_rt::Bezier3FormalSolver solver(g.ND_MID,g.NA,g.NE,dt.data());
	solver.solve(incidence,g.mu,top.data(),bottom.data(),source.data(),I.data(),ft.data(),fb.data(),mean.data());
	for(int e=0;e<g.NE;++e) {
		double incoming=0,outgoing=0,emitted=0,absorbed=0;
		for(int m=0;m<g.NA;++m) {
			double factor=g.wt[m]*std::abs(g.mu[m]);
			incoming+=factor*(g.mu[m]<0 ? r.Inu_top[m][e] : r.Inu_bottom[m][e]);
			outgoing+=factor*(g.mu[m]>0 ? r.Inu_top[m][e] : r.Inu_bottom[m][e]);
		}
		for(int d=0;d<g.ND_MID;++d) {
			emitted+=2*g.dr[d]*r.jnu[d][e];
			for(int m=0;m<g.NA;++m)
				absorbed+=g.wt[m]*g.dr[d]*r.kabs[d][e]*mean[(long(d)*g.NA+m)*g.NE+e];
		}
		double residual=(outgoing+absorbed-incoming-emitted)/(incoming+emitted);
		check(std::isfinite(residual) && std::abs(residual)<2e-4,"elastic production RT energy budget");
		std::printf("  budget E=%g eV residual=%.3e\n",g.ene[e],residual);
	}
}

int main()
{
	char temp[]="/tmp/dao-cell-integration-XXXXXX";
	if(!mkdtemp(temp)) return 1;
	try {
		ModelParams p{};p.nh=std::log10(1.0/(phys::reference_electrons_per_hydrogen*phys::sigma_T));
		RTGrids g;g.init_angle();g.init_energy(100,400,3);
		// Resolve this homogeneous test slab independently of production defaults.
		std::vector<double> mids(200), widths(200,2.0/200);
		for(int d=0;d<200;++d) mids[d]=(d+0.5)*widths[d];
		g.init_depth_from_zones(200,mids.data(),widths.data());
		RadField r(g);r.allocate();
		KernelCache directional;avgKernelCache averaged;
		directional.NA_full=g.NA;directional.n_indep=g.NA*g.NA;
		// A full angular table is allowed here: each pair has its own row.
		directional.canon=new int[g.NA*g.NA];
		for(int a=0;a<g.NA*g.NA;++a) directional.canon[a]=a;
		elastic_storage(directional,g,g.NA*g.NA);elastic_storage(averaged,g,1);
		for(int e=0;e<g.NE;++e) {
			double dx=e==0 ? 0.5*(averaged.x_grid[1]-averaged.x_grid[0]) :
			          e==g.NE-1 ? 0.5*(averaged.x_grid[e]-averaged.x_grid[e-1]) :
			          0.5*(averaged.x_grid[e+1]-averaged.x_grid[e-1]);
			averaged.data[e]=1.0/dx;
			for(int a=0;a<g.NA*g.NA;++a) directional.data[e*g.NA*g.NA+a]=0.5/dx;
			r.illum.I_corona[e]=1+e;r.illum.I_disk[e]=0.1*(e+1);
		}
		for(int d=0;d<g.ND_MID;++d) r.T_K[d]=1e6;
		for(double absorption : {0.0,0.7}) for(int incidence : {-1,1}) {
			p.i_incidence=incidence;
			for(int d=0;d<g.ND_MID;++d) for(int e=0;e<g.NE;++e) {
				r.ksct[d][e]=1.0;r.kabs[d][e]=absorption;
				r.jnu[d][e]=absorption>0 ? 0.02*(e+1) : 0;
			}
			compton_rt_solve(r,g,p,directional,1000);
			verify_budget_and_moments(r,g,incidence);auto dir=faces(r,g);
			compton_rt_solve(r,g,p,averaged,1000);
			verify_budget_and_moments(r,g,incidence);auto avg=faces(r,g);
			for(size_t i=0;i<dir.size();++i)
				check(std::abs(dir[i]-avg[i])<1e-11*std::max(1.0,dir[i]),"elastic angular/averaged dispatch mismatch");
			std::printf("PASS production RT: absorption=%g incidence=%d, both kernel dispatch paths\n",absorption,incidence);
		}
		bool rejected=false;
		try {compton_rt_solve(r,g,p,averaged,1);} catch(const std::runtime_error& e) {
			rejected=std::string(e.what()).find("iteration limit")!=std::string::npos;
		}
		check(rejected,"unconverged RT silently accepted");
		r.jnu[0][0]=std::numeric_limits<double>::quiet_NaN();rejected=false;
		try {compton_rt_solve(r,g,p,averaged,1000);} catch(const std::runtime_error&) {rejected=true;}
		check(rejected,"nonfinite emissivity silently accepted");r.jnu[0][0]=0;

		// The saved emergent spectrum must read the boundary, not cell zero.
		for(int m=0;m<g.NA;++m) for(int e=0;e<g.NE;++e) r.Inu[0][m][e]=123456789;
		std::snprintf(p.run_dir,sizeof(p.run_dir),"%s",temp);save_results(r,g,p,1);
		std::ifstream f(std::string(temp)+"/emergent_iter001.dat");std::string line;int e=0;
		while(std::getline(f,line)) {
			if(line.empty() || line[0]=='#') continue;
			std::istringstream row(line);double energy,corona,disk;row>>energy>>corona>>disk;
			for(int m=0;m<g.NA;++m) {double value;check(bool(row>>value),"bad emergent output");
				check(std::abs(value-r.Inu_top[m][e])<5e-6*std::max(1.0,std::abs(value)),"emergent output used cell intensity");}
			++e;
		}
		check(e==g.NE,"missing emergent rows");
		std::puts("PASS center moments, physical boundary output, and explicit nonconvergence/nonfinite failures");
		directional.free_memory();averaged.free_memory();r.deallocate();
		std::filesystem::remove_all(temp);return 0;
	} catch(const std::exception& e) {
		std::fprintf(stderr,"FAIL: %s (diagnostics in %s)\n",e.what(),temp);return 1;
	}
}
