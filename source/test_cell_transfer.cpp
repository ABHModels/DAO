// Standalone spatial-transfer regressions; no Cloudy or kernel cache required.
#include "cell_transfer.h"
#include "rt_grids.h"
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

static void require(bool ok, const char* message)
{
	if (!ok) throw std::runtime_error(message);
}

static void near(double actual, double expected, double tolerance, const char* message)
{
	require(std::isfinite(actual) && std::abs(actual-expected) <=
	        tolerance*std::max(1.0, std::abs(expected)), message);
}

static void analytic_slab(const RTGrids& g)
{
	const int ND=5, NM=g.NA, NE=1;
	const double widths[ND]={1e-12, 0.05, 0.4, 0.8, 1.06};
	const double total=widths[0]+widths[1]+widths[2]+widths[3]+widths[4];
	std::vector<dao_rt::CellWeights> weights(ND*NM);
	std::vector<double> S(ND*NM,2.0), I(S.size()), top(NM), bottom(NM);
	for(int d=0;d<ND;++d) for(int m=0;m<NM;++m)
		weights[d*NM+m]=dao_rt::cell_weights(widths[d]/std::abs(g.mu[m]));
	for(int incidence : {-1,1}) {
		double input_top=incidence==-1 ? 2.0 : 2.0/g.wt[incidence], input_bottom=0.6;
		dao_rt::formal_solution_cells(ND,NM,NE,incidence,g.mu,&input_top,&input_bottom,
			S.data(),weights.data(),I.data(),top.data(),bottom.data());
		double net=0, production=0;
		for(int m=0;m<NM;++m) {
			const bool down=g.mu[m]<0;
			const double incoming=down ? incident_top_intensity(incidence,m,input_top) : input_bottom;
			const double outgoing=down ? bottom[m] : top[m];
			const double expected=2.0+(incoming-2.0)*std::exp(-total/std::abs(g.mu[m]));
			near(outgoing,expected,2e-14,"homogeneous slab analytic surface intensity");
			near(down ? top[m] : bottom[m],incoming,1e-14,"incoming boundary changed");
			net+=g.wt[m]*std::abs(g.mu[m])*(outgoing-incoming);
			for(int d=0;d<ND;++d)
				production+=g.wt[m]*widths[d]*(2.0-I[d*NM+m]);
		}
		near(net,production,2e-13,"emission minus absorption does not equal boundary flux");
		require(std::abs(I[3]-top[3])>1e-14 || std::abs(I[(ND-1)*NM+4]-bottom[4])>0.01,
		        "cell averages confused with boundaries");
	}
	// LTE radiation stays constant through arbitrarily opaque cells.
	for(double t : {0.0,1e-12,1e-6,1.0,1e4}) {
		double incident=2.0;
		std::fill(weights.begin(),weights.end(),dao_rt::cell_weights(t));
		dao_rt::formal_solution_cells(ND,NM,NE,-1,g.mu,&incident,&incident,S.data(),
			weights.data(),I.data(),top.data(),bottom.data());
		for(double value:I) near(value,incident,2e-14,"LTE/vacuum invariant");
		for(double value:top) near(value,incident,2e-14,"LTE upper boundary");
		for(double value:bottom) near(value,incident,2e-14,"LTE lower boundary");
	}
	std::printf("PASS analytic slab, cell energy balance, LTE, both boundaries and illumination modes\n");
}

struct Reflection { std::vector<double> reflected, transmitted, absorbed; };

static Reflection coherent_slab(const RTGrids& g, int NE,
	const std::vector<double>& absorption, const std::vector<double>& scattering,
	int subdivisions, bool thomson, int incidence=1)
{
	const int NM=g.NA, original=int(absorption.size())/NE, ND=original*subdivisions;
	const long size=long(ND)*NM*NE;
	std::vector<double> S(size), next(size), I(size), top(NM*NE), bottom(NM*NE);
	std::vector<double> input_top(NE,incidence==-1 ? 2.0 : 2.0/g.wt[incidence]), input_bottom(NE);
	std::vector<dao_rt::CellWeights> weights(size);
	for(int d=0;d<ND;++d) for(int m=0;m<NM;++m) for(int e=0;e<NE;++e) {
		int old=(d/subdivisions)*NE+e;
		weights[(long(d)*NM+m)*NE+e]=dao_rt::cell_weights(
			(absorption[old]+scattering[old])/(subdivisions*std::abs(g.mu[m])));
	}
	bool converged=false;
	for(int it=0;it<3000;++it) {
		dao_rt::formal_solution_cells(ND,NM,NE,incidence,g.mu,input_top.data(),input_bottom.data(),
			S.data(),weights.data(),I.data(),top.data(),bottom.data());
		double error=0;
		for(int d=0;d<ND;++d) for(int m=0;m<NM;++m) for(int e=0;e<NE;++e) {
			double j=0;
			for(int n=0;n<NM;++n) {
				const double a=g.mu[m]*g.mu[m], b=g.mu[n]*g.mu[n];
				const double phase=thomson ? 3.0/16.0*(3-a-b+3*a*b) : 0.5;
				j+=g.wt[n]*phase*I[(long(d)*NM+n)*NE+e];
			}
			int old=(d/subdivisions)*NE+e;
			long k=(long(d)*NM+m)*NE+e;
			next[k]=scattering[old]/(absorption[old]+scattering[old])*j;
			error=std::max(error,std::abs(next[k]-S[k]));
		}
		if(error<1e-11) {converged=true;break;}
		S.swap(next);
	}
	require(converged,"coherent reference iteration did not converge");
	Reflection result{std::vector<double>(NE),std::vector<double>(NE),std::vector<double>(NE)};
	for(int e=0;e<NE;++e) {
		double incoming=0;
		for(int m=0;m<NM;++m) {
			if(g.mu[m]>0) result.reflected[e]+=g.wt[m]*g.mu[m]*top[m*NE+e];
			else {
				incoming-=g.wt[m]*g.mu[m]*top[m*NE+e];
				result.transmitted[e]-=g.wt[m]*g.mu[m]*bottom[m*NE+e];
			}
			for(int d=0;d<ND;++d)
				result.absorbed[e]+=g.wt[m]*I[(long(d)*NM+m)*NE+e]*
					absorption[(d/subdivisions)*NE+e]/subdivisions;
		}
		result.reflected[e]/=incoming;
		result.transmitted[e]/=incoming;
		result.absorbed[e]/=incoming;
		near(result.reflected[e]+result.transmitted[e]+result.absorbed[e],1.0,2e-9,
		     "coherent reflection + transmission + absorption budget");
	}
	return result;
}

int main(int argc,char** argv)
{
	try {
		RTGrids g;g.init_angle();
		for(double t : {0.0,1e-12,1e-6,1e-3,0.1,1.0,10.0,1e4}) {
			auto w=dao_rt::cell_weights(t);
			near(w.attenuation+w.emission,1,1e-14,"attenuation/emission weights");
			near(w.mean_incoming+w.mean_source,1,1e-14,"mean weights");
			near(t*w.mean_incoming,w.emission,1e-14,"cell integral identity");
			require(w.mean_source>=0 && w.mean_incoming>=0,"negative formal weight");
		}
		near(dao_rt::cell_weights(1e-12).mean_source/1e-12,0.5,1e-12,"thin-cell cancellation");
		analytic_slab(g);
		for(bool thomson : {false,true}) for(int incidence : {-1,1}) {
			auto r=coherent_slab(g,1,std::vector<double>(20,0),std::vector<double>(20,0.1),1,thomson,incidence);
			require(r.reflected[0]>0 && r.transmitted[0]>0,"empty pure-scattering result");
		}
		std::printf("PASS conservative pure scattering, isotropic/Thomson phase functions and both incidence modes\n");
		std::ifstream f(argc>1 ? argv[1] : "tests/data/soft_front_opacity.dat");
		require(bool(f),"cannot open opacity-front fixture");
		std::vector<double> absorption,scattering;
		std::string line;
		while(std::getline(f,line)) {
			if(line.empty() || line[0]=='#') continue;
			std::istringstream row(line);
			for(int e=0;e<4;++e) {double a,s;require(bool(row>>a>>s),"bad fixture row");absorption.push_back(a);scattering.push_back(s);}
		}
		require(absorption.size()==200,"wrong opacity-front fixture size");
		// Independent 256-subcell reference, fixed original material properties.
		const double reference[]={0.412987766,0.406286617,0.380593917,0.163652705};
		for(int sub : {1,16,64}) {
			auto r=coherent_slab(g,4,absorption,scattering,sub,true);
			for(int e=0;e<4;++e)
				near(r.reflected[e]/reference[e],1.0,sub==1 ? 0.015 : (sub==16 ? 0.001 : 0.0001),
				     "cold opacity front: excess reflection or spatial regression");
			std::printf("PASS opacity front sub=%d R(0.30,0.40,0.80,0.87 keV) = %.9f %.9f %.9f %.9f\n",
				sub,r.reflected[0],r.reflected[1],r.reflected[2],r.reflected[3]);
		}
		return 0;
	} catch(const std::exception& e) {std::fprintf(stderr,"FAIL: %s\n",e.what());return 1;}
}
