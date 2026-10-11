#include "compton_rt.h"
#include "cell_transfer.h"
#include "rt_parallel.h"
#include "run_log.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace {
using dao_thermal::Scattering;

Scattering scattering_at(const avgKernelCache& cache, double T, double electron_density, const RTGrids&)
{ return Scattering(cache,T,electron_density); }
Scattering scattering_at(const KernelCache& cache, double T, double electron_density, const RTGrids& g)
{ return Scattering(cache,T,electron_density,g.NA,g.wt); }

template<class Cache>
std::vector<Scattering> build_scattering_column(RadField& rad, const RTGrids& g,
                                          const Cache& cache)
{
    std::vector<Scattering> ops(g.ND_MID);
    bool uniform=g.ND_MID>0;
    for(int d=1;d<g.ND_MID;++d)
        uniform=uniform && rad.T_K[d]==rad.T_K[0] && rad.n_e[d]==rad.n_e[0];
    if(uniform) {
        const auto op=scattering_at(cache,rad.T_K[0],rad.n_e[0],g);
        for(int d=0;d<g.ND_MID;++d) {
            ops[d]=op;
            std::copy(op.opacity().begin(),op.opacity().end(),rad.ksct[d]);
        }
        return ops;
    }
    rt_parallel_depths(g.ND_MID,[&](int d) {
        ops[d]=scattering_at(cache,rad.T_K[d],rad.n_e[d],g);
        std::copy(ops[d].opacity().begin(),ops[d].opacity().end(),rad.ksct[d]);
    });
    return ops;
}

} // namespace

std::vector<dao_thermal::Scattering> make_scattering_column(
    RadField& rad,const RTGrids& g,const avgKernelCache& cache)
{ return build_scattering_column(rad,g,cache); }
std::vector<dao_thermal::Scattering> make_scattering_column(
    RadField& rad,const RTGrids& g,const KernelCache& cache)
{ return build_scattering_column(rad,g,cache); }

using dao_thermal::Scattering;

void compton_rt_solve(RadField& rad, const RTGrids& g, const ModelParams& p,
                      const std::vector<Scattering>& ops,std::vector<double>* cell_incoming)
{
    const int nd=g.ND_MID, nm=g.NA, ne=g.NE;
    if(int(ops.size())!=nd || p.maxiter<3)
        throw std::runtime_error("RT requires one operator per cell and maxiter>=3");
    for(const auto& op:ops) if(int(op.opacity().size())!=ne || op.angles()!=(p.angsca ? nm : 1))
        throw std::runtime_error("RT: scattering operator energy/angular grid mismatch");
    const long size=long(nd)*nm*ne;
    if(cell_incoming) cell_incoming->resize(size);
    std::vector<double> source(size),intensity(size),mean(long(nd)*ne),next(mean.size());
    std::vector<double> angular(p.angsca ? size : 0);
    if(p.angsca) for(int d=0;d<nd;++d) for(int m=0;m<nm;++m)
        std::copy(rad.Inu[d][m],rad.Inu[d][m]+ne,angular.data()+(long(d)*nm+m)*ne);
    std::vector<double> top(ne),bottom(ne),ft(long(nm)*ne),fb(long(nm)*ne);
    std::vector<dao_rt::CellWeights> weights(size);
    compute_boundary_illumination(ne,p.i_incidence,g,rad,top.data(),bottom.data());
    for(int e=0;e<ne;++e)
        if(!std::isfinite(top[e]) || !std::isfinite(bottom[e]) || top[e]<0 || bottom[e]<0)
            throw std::runtime_error("RT: invalid boundary intensity");
    for(int d=0;d<nd;++d) for(int e=0;e<ne;++e) {
        const double k=rad.kabs[d][e]+rad.ksct[d][e];
        if(!std::isfinite(k) || k<=0 || !std::isfinite(rad.jnu[d][e]) || rad.jnu[d][e]<0 ||
           !std::isfinite(rad.J0[d][e]) || rad.J0[d][e]<0)
            throw std::runtime_error("RT: invalid cell coefficients");
        mean[long(d)*ne+e]=rad.J0[d][e];
        for(int m=0;m<nm;++m)
            weights[(long(d)*nm+m)*ne+e]=dao_rt::cell_weights(k*g.dr[d]/std::abs(g.mu[m]));
    }
    int successive=0;
    for(int iteration=1;iteration<=p.maxiter;++iteration) {
        rt_parallel_depths(nd,[&](int d) {
            if(p.angsca) {
                std::vector<double> j(nm*ne);
                ops[d].emissivity(angular.data()+long(d)*nm*ne,j.data());
                for(int m=0;m<nm;++m) for(int e=0;e<ne;++e)
                    source[(long(d)*nm+m)*ne+e]=(rad.jnu[d][e]+j[m*ne+e])/(rad.kabs[d][e]+rad.ksct[d][e]);
            } else {
                std::vector<double> j(ne);
                ops[d].emissivity(mean.data()+long(d)*ne,j.data());
                for(int e=0;e<ne;++e) {
                    const double S=(rad.jnu[d][e]+j[e])/(rad.kabs[d][e]+rad.ksct[d][e]);
                    for(int m=0;m<nm;++m) source[(long(d)*nm+m)*ne+e]=S;
                }
            }
        });
        dao_rt::formal_solution_cells(nd,nm,ne,p.i_incidence,g.mu,top.data(),bottom.data(),
                                      source.data(),weights.data(),intensity.data(),ft.data(),fb.data(),
                                      cell_incoming ? cell_incoming->data() : nullptr);
        double change=0;
        for(int d=0;d<nd;++d) for(int e=0;e<ne;++e) {
            double J=0;
            for(int m=0;m<nm;++m) J+=0.5*g.wt[m]*intensity[(long(d)*nm+m)*ne+e];
            if(!std::isfinite(J) || J<0) throw std::runtime_error("RT: invalid intensity");
            const long i=long(d)*ne+e;
            if(!p.angsca) change=std::max(change,std::abs(J-mean[i])/std::max({J,mean[i],1e-100}));
            next[i]=J;
        }
        if(p.angsca) for(long i=0;i<size;++i) {
            if(!std::isfinite(intensity[i]) || intensity[i]<0)
                throw std::runtime_error("RT: invalid directional intensity");
            change=std::max(change,std::abs(intensity[i]-angular[i])/std::max({intensity[i],angular[i],1e-100}));
            angular[i]=intensity[i];
        }
        mean.swap(next);
        successive=change<1e-7 ? successive+1 : 0;
        if(iteration==1 || iteration%25==0 || successive==3)
            dao_log::detail("  [RT] %d max relative d%s %.6e\n",iteration,p.angsca ? "I" : "J",change);
        if(successive==3) break;
    }
    if(successive!=3) throw std::runtime_error("RT: source iteration did not converge");
    for(int d=0;d<nd;++d) for(int e=0;e<ne;++e) {
        rad.J0[d][e]=mean[long(d)*ne+e];
        for(int m=0;m<nm;++m) rad.Inu[d][m][e]=intensity[(long(d)*nm+m)*ne+e];
    }
    for(int m=0;m<nm;++m) for(int e=0;e<ne;++e) {
        rad.Inu_top[m][e]=ft[long(m)*ne+e]; rad.Inu_bottom[m][e]=fb[long(m)*ne+e];
    }
}
