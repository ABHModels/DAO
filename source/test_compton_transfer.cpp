#include "compton_rt.h"
#include "incidence_boundary.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <stdexcept>
#include <unistd.h>

static void check(bool ok,const char* message)
{ if(!ok) throw std::runtime_error(message); }

// Follow production illumination all the way through the shared slab solver.
// A known absorbing slab checks individual rays; a scattering slab checks
// both faces and the gas/radiation energy exchange on the 48-cell depth grid.
template<class Cache>
static void isotropic_production_boundary(Cache& cache,RTGrids& g,bool directional,
                                          double temperature)
{
    ModelParams p{};p.angsca=directional;p.maxiter=1000;p.i_incidence=-1;
    p.incidence=-2;p.nh=15;p.zeta=3;p.frac=-1;
    std::strcpy(p.corona,"nthcomp");
    p.Gamma=2;p.kT_e=60;p.kT_bb=.01;p.kT_disk=.35;
    RadField r(g);r.allocate();r.illum.compute(p);
    const int ne=g.NE,nd=g.ND_MID,nm=g.NA;
    const auto w=dao_thermal::energy_weights(ne,g.ene);
    double width=0,hemisphere_flux_weight=0;
    for(int d=0;d<nd;++d) width+=g.dr[d];
    for(int m=0;m<nm;++m) if(g.mu[m]<0)
        hemisphere_flux_weight+=0.5*phys::four_pi*g.wt[m]*std::abs(g.mu[m]);
    const double target=std::pow(10.,p.nh+p.zeta)/phys::four_pi;
    for(bool scatter:{false,true}) {
        for(int d=0;d<nd;++d) {
            r.T_K[d]=temperature;r.n_e[d]=scatter ? 1.21e15 : 0;
            for(int e=0;e<ne;++e) {
                r.kabs[d][e]=scatter ? 0 : .7/width;
                r.jnu[d][e]=r.J0[d][e]=0;
                for(int m=0;m<nm;++m) r.Inu[d][m][e]=0;
            }
        }
        auto ops=make_scattering_column(r,g,cache);
        std::vector<double> entering;
        compton_rt_solve(r,g,p,ops,&entering);
        double fin=0,fout=0,nin=0,nout=0,volume=0,worst=0;
        for(int e=0;e<ne;++e) {
            const double boundary=r.illum.I_corona[e]/hemisphere_flux_weight;
            for(int m=0;m<nm;++m) {
                const double in=g.mu[m]<0 ? r.Inu_top[m][e] : r.Inu_bottom[m][e];
                const double out=g.mu[m]<0 ? r.Inu_bottom[m][e] : r.Inu_top[m][e];
                check(std::isfinite(out) && out>=0,"finite emergent intensity for isotropic incidence");
                check(std::abs(in-(g.mu[m]<0 ? boundary : 0))<=1e-12*boundary,
                      "actual incoming boundary must illuminate every downward ray equally");
                const double factor=0.5*phys::four_pi*g.wt[m]*std::abs(g.mu[m])*w[e];
                fin+=factor*in;fout+=factor*out;nin+=factor*in/g.ene[e];nout+=factor*out/g.ene[e];
                if(!scatter && g.mu[m]<0) {
                    const double exact=boundary*std::exp(-.7/std::abs(g.mu[m]));
                    worst=std::max(worst,std::abs(out/exact-1));
                    double distance=0;
                    for(int d=0;d<nd;++d) {
                        const double face=boundary*std::exp(-.7*distance/(width*std::abs(g.mu[m])));
                        const double t=.7*g.dr[d]/(width*std::abs(g.mu[m]));
                        const double average=face*(-std::expm1(-t))/t;
                        worst=std::max(worst,std::abs(entering[(long(d)*nm+m)*ne+e]/face-1));
                        worst=std::max(worst,std::abs(r.Inu[d][m][e]/average-1));
                        distance+=g.dr[d];
                    }
                }
            }
        }
        std::vector<double> intensity(nm*ne);
        for(int d=0;d<nd;++d) {
            for(int m=0;m<nm;++m) std::copy(r.Inu[d][m],r.Inu[d][m]+ne,intensity.data()+m*ne);
            volume+=dao_thermal::budget(ne,w.data(),r.J0[d],r.kabs[d],r.jnu[d],ops[d],
                                        directional ? intensity.data() : nullptr).residual()*g.dr[d];
        }
        const double flux_error=std::abs(fin/target-1);
        const double identity=std::abs((fout-fin+volume)/fin);
        std::printf("Isotropic %s %s %d cells: Fin=%.12e relative flux error %.3e analytic %.3e energy identity %.3e photons %.3e Fout/Fin %.6f\n",
                    directional ? "directional" : "mean",scatter ? "scattering" : "absorption",nd,
                    fin,flux_error,worst,identity,scatter ? std::abs(nout/nin-1) : 0,fout/fin);
        check(flux_error<1e-12,"flux reaching the slab equals xi*nH/(4*pi)");
        check(worst<1e-12,"isotropic cell and transmitted intensities match the analytic solution");
        // Use production's 1e-4 identity tolerance for the hot, thick slab:
        // its finite source-iteration error is measured relative to Fin even
        // when the externally maintained electrons supply much larger power.
        check(identity<(scatter ? 1e-4 : 1e-12),"isotropic slab boundary/volume energy identity");
        if(scatter) check(std::abs(nout/nin-1)<1e-6,"isotropic slab two-face photon conservation");
        std::puts("PASS isotropic production boundary and transfer");
    }
    r.deallocate();
}

template<class Cache>
static void exercise(Cache& cache,RTGrids& g,bool directional,double temperature)
{
    ModelParams p{};p.angsca=directional;p.maxiter=1000;p.i_incidence=-1;
    RadField r(g);r.allocate();
    const int ne=g.NE,nd=g.ND_MID,nm=g.NA;
    const auto w=dao_thermal::energy_weights(ne,g.ene);
    double width=0;for(int d=0;d<nd;++d) width+=g.dr[d];
    const double k=0.7/width;
    for(int d=0;d<nd;++d) {
        r.T_K[d]=temperature;r.n_e[d]=0;
        for(int e=0;e<ne;++e) { r.kabs[d][e]=k;r.jnu[d][e]=0.3*k; }
    }
    // Store the flux of an isotropic unit-intensity upper boundary.
    double unit_flux=0;
    for(int m=0;m<nm;++m) if(g.mu[m]<0)
        unit_flux+=0.5*phys::four_pi*g.wt[m]*std::abs(g.mu[m]);
    for(int e=0;e<ne;++e) { r.illum.I_corona[e]=unit_flux;r.illum.I_disk[e]=1; }
    auto ops=make_scattering_column(r,g,cache);
    std::vector<double> entering;
    compton_rt_solve(r,g,p,ops,&entering);
    double worst=0;
    for(int m=0;m<nm;++m) for(int e=0;e<ne;++e) {
        const double incoming=g.mu[m]<0 ? 1. : 2.;
        const double exact=0.3+(incoming-0.3)*std::exp(-k*width/std::abs(g.mu[m]));
        const double out=g.mu[m]<0 ? r.Inu_bottom[m][e] : r.Inu_top[m][e];
        worst=std::max(worst,std::abs(out/exact-1));
        double distance=0;
        for(int step=0;step<nd;++step) {
            const int d=g.mu[m]<0 ? step : nd-1-step;
            const double t=k*g.dr[d]/std::abs(g.mu[m]);
            const double face=0.3+(incoming-0.3)*std::exp(-k*distance/std::abs(g.mu[m]));
            const double average=0.3+(face-0.3)*(-std::expm1(-t))/t;
            worst=std::max(worst,std::abs(r.Inu[d][m][e]/average-1));
            check(std::abs(entering[(long(d)*nm+m)*ne+e]/face-1)<1e-12,"incoming face intensity");
            distance+=g.dr[d];
        }
    }
    check(worst<1e-12,"constant-source analytic cell means and boundaries");
    std::printf("PASS %s analytic absorption/emission: %.3e\n",directional ? "directional" : "mean",worst);

    // Pure scattering: prescribed hot electrons can supply photon energy,
    // but may not create photons. Include BOTH faces and the signed gas Q.
    for(int d=0;d<nd;++d) {
        r.n_e[d]=1.21e15*(1+0.1*d); // also exercises nonuniform operator construction
        for(int e=0;e<ne;++e) {
            r.kabs[d][e]=r.jnu[d][e]=r.J0[d][e]=0;
            for(int m=0;m<nm;++m) r.Inu[d][m][e]=0;
        }
    }
    for(int e=0;e<ne;++e) {
        r.illum.I_corona[e]=0;
        r.illum.I_disk[e]=g.ene[e]*g.ene[e]*std::exp(-g.ene[e]/2000);
    }
    ops=make_scattering_column(r,g,cache);
    bool failed=false;p.maxiter=3;
    try { compton_rt_solve(r,g,p,ops); } catch(const std::runtime_error&) { failed=true; }
    check(failed,"unconverged source iteration must fail");p.maxiter=1000;
    compton_rt_solve(r,g,p,ops);
    double fin=0,fout=0,nin=0,nout=0,volume=0;
    for(int m=0;m<nm;++m) for(int e=0;e<ne;++e) {
        const double factor=0.5*phys::four_pi*g.wt[m]*std::abs(g.mu[m])*w[e];
        const double in=g.mu[m]<0 ? r.Inu_top[m][e] : r.Inu_bottom[m][e];
        const double out=g.mu[m]<0 ? r.Inu_bottom[m][e] : r.Inu_top[m][e];
        check(std::isfinite(out) && out>=0,"finite nonnegative emergent spectrum");
        fin+=factor*in;fout+=factor*out;nin+=factor*in/g.ene[e];nout+=factor*out/g.ene[e];
    }
    std::vector<double> intensity(nm*ne);
    for(int d=0;d<nd;++d) {
        for(int m=0;m<nm;++m) std::copy(r.Inu[d][m],r.Inu[d][m]+ne,intensity.data()+m*ne);
        volume+=dao_thermal::budget(ne,w.data(),r.J0[d],r.kabs[d],r.jnu[d],ops[d],
                                    directional ? intensity.data() : nullptr).residual()*g.dr[d];
    }
    const double number_error=std::abs(nout/nin-1),identity=std::abs((fout-fin+volume)/fin);
    check(number_error<1e-7,"two-face photon conservation");
    check(identity<1e-7,"boundary/volume energy identity");
    std::printf("PASS %s scattering: photons %.3e energy identity %.3e Fout/Fin %.6f\n",
                directional ? "directional" : "mean",number_error,identity,fout/fin);

    // Uniform Wien radiation is a stationary solution of the same operator.
    const double kT=phys::k_B*temperature/phys::eV_to_erg;
    for(int d=0;d<nd;++d) r.n_e[d]=1.21e15;
    for(int e=0;e<ne;++e) {
        const double B=std::pow(g.ene[e]/kT,3)*std::exp(-g.ene[e]/kT);
        r.illum.I_corona[e]=unit_flux*B; r.illum.I_disk[e]=0.5*B;
        for(int d=0;d<nd;++d) {
            r.J0[d][e]=B;
            for(int m=0;m<nm;++m) r.Inu[d][m][e]=B;
        }
    }
    ops=make_scattering_column(r,g,cache);compton_rt_solve(r,g,p,ops);
    worst=0;
    for(int d=0;d<nd;++d) for(int m=0;m<nm;++m) for(int e=0;e<ne;++e)
        worst=std::max(worst,std::abs(r.Inu[d][m][e]/(2*r.illum.I_disk[e])-1));
    check(worst<1e-8,"Wien stationary radiation field");
    std::printf("PASS %s uniform Wien equilibrium: %.3e\n",directional ? "directional" : "mean",worst);
    r.deallocate();
}

int main()
{
    try {
        char tmp[]="/tmp/dao-shared-rt-XXXXXX";check(mkdtemp(tmp),"temporary directory");
        setenv("KERNEL_DIR",tmp,1);setenv("COMPTON_CACHE_DIR",tmp,1);
        RTGrids g;g.init_energy(10,1e5,40);
        const double temperature=60e3*phys::eV_to_erg/phys::k_B;
        double dr[]={1e7,5e7,9e7,2e8},tau[4],depth=0;
        for(int d=0;d<4;++d) { tau[d]=(depth+0.5*dr[d])*1.21e15*phys::sigma_T;depth+=dr[d]; }
        for(bool double_gauss:{false,true}) {
            if(double_gauss) g.init_angle_double_gauss(); else g.init_angle();
            // Keep kernels for distinct angular grids in separate caches.
            const std::string cache_dir=std::string(tmp)+(double_gauss ? "/double" : "/full");
            std::filesystem::create_directory(cache_dir);
            setenv("KERNEL_DIR",cache_dir.c_str(),1);
            setenv("COMPTON_CACHE_DIR",cache_dir.c_str(),1);
            std::printf("Angular grid: %s\n",double_gauss ? "double-Gauss" : "production full-Gauss");
            g.init_depth_from_zones(4,tau,dr);
            avgKernelCache mean;mean.init(g.NE,g.ene,g.NA,g.mu,g.wt,1,1,&temperature);
            exercise(mean,g,false,temperature);
            g.init_depth(g.TAU_MIN,g.TAU_MAX,15);
            isotropic_production_boundary(mean,g,false,temperature);mean.free_memory();
            g.init_depth_from_zones(4,tau,dr);
            KernelCache directional;directional.init(g.NE,g.ene,g.NA,g.mu,g.wt,1,1,&temperature);
            exercise(directional,g,true,temperature);
            g.init_depth(g.TAU_MIN,g.TAU_MAX,15);
            isotropic_production_boundary(directional,g,true,temperature);directional.free_memory();
        }
        std::filesystem::remove_all(tmp);
    } catch(const std::exception& e) { std::fprintf(stderr,"FAIL: %s\n",e.what());return 1; }
}
