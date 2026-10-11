#include "thermal_production.h"
#include "cell_transfer.h"
#include "compton_rt.h"
#include "cloudy_depth.h"
#include "save_results.h"
#include "rt_parallel.h"
#include "run_log.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>

namespace {
using dao_thermal::Scattering;
using dao_thermal::Budget;

Scattering scattering_at(const avgKernelCache& cache, double T, double electron_density, const RTGrids&)
{ return Scattering(cache,T,electron_density); }
Scattering scattering_at(const KernelCache& cache, double T, double electron_density, const RTGrids& g)
{ return Scattering(cache,T,electron_density,g.NA,g.wt); }

Budget cell_budget(const RadField& rad,const RTGrids& g,int d,
                   const std::vector<double>& weights,const Scattering& op)
{
    std::vector<double> intensity;
    if(op.angles()>1) {
        intensity.resize(g.NA*g.NE);
        for(int m=0;m<g.NA;++m)
            std::copy(rad.Inu[d][m],rad.Inu[d][m]+g.NE,intensity.data()+m*g.NE);
    }
    return dao_thermal::budget(g.NE,weights.data(),rad.J0[d],rad.kabs[d],rad.jnu[d],op,
                               intensity.empty() ? nullptr : intensity.data());
}

struct ColumnBudget {
    double incoming=0,outgoing=0,volume=0,volume_abs=0,max_local=0;
    double closure() const { return (outgoing-incoming)/incoming; }
    double identity_error() const { return (outgoing-incoming+volume)/incoming; }
};

ColumnBudget column_budget(const RadField& rad, const RTGrids& g,
                           const std::vector<Scattering>& ops,
                           const std::string& path)
{
    auto w=dao_thermal::energy_weights(g.NE,g.ene);
    ColumnBudget total;
    std::ofstream f(path); f<<std::setprecision(15);
    if(!f) throw std::runtime_error("Cannot write thermal budget: "+path);
    f<<"# Cell-volume averages; R=4pi integral[kappa_tot*J-j_material-j_Compton] dE\n"
     <<"# j_Compton is angular averaged; directional mode evaluates it from full I(E,mu)\n"
     <<"# depth tau T absorption emission Q_Compton R relative_R dr\n";
    for(int d=0;d<g.ND_MID;++d) {
        const Budget b=cell_budget(rad,g,d,w,ops[d]);
        total.volume+=b.residual()*g.dr[d];
        total.volume_abs+=std::abs(b.residual())*g.dr[d];
        total.max_local=std::max(total.max_local,std::abs(b.relative()));
        f<<d<<' '<<g.tau_mid[d]<<' '<<rad.T_K[d]<<' '<<b.absorption<<' '<<b.emission
         <<' '<<b.compton<<' '<<b.residual()<<' '<<b.relative()<<' '<<g.dr[d]<<'\n';
    }
    for(int m=0;m<g.NA;++m) for(int e=0;e<g.NE;++e) {
        const double factor=0.5*phys::four_pi*g.wt[m]*std::abs(g.mu[m])*w[e];
        total.incoming+=factor*(g.mu[m]<0 ? rad.Inu_top[m][e] : rad.Inu_bottom[m][e]);
        total.outgoing+=factor*(g.mu[m]>0 ? rad.Inu_top[m][e] : rad.Inu_bottom[m][e]);
    }
    if(!std::isfinite(total.incoming) || !std::isfinite(total.outgoing) || total.incoming<=0)
        throw std::runtime_error("Invalid total boundary flux");
    f<<"# Fin "<<total.incoming<<" Fout "<<total.outgoing<<" closure "<<total.closure()
     <<" volume_R "<<total.volume<<" identity_error "<<total.identity_error()<<'\n';
    return total;
}

// Keep failure state separate from spectra: an unfinished run is never labelled
// converged just because it has an emergent spectrum on disk.
void status(const ModelParams& p, const char* state, int iteration, const std::string& reason="")
{
    std::ofstream f(std::string(p.run_dir)+"/thermal_status.json");
    if(!f) throw std::runtime_error("Cannot write thermal status");
    f<<"{\"mode\":\"rt_energy_balance_v1\",\"state\":\""<<state
     <<"\",\"iteration\":"<<iteration<<",\"reason\":\"";
    for(unsigned char c:reason) {
        if(c=='"' || c=='\\') f<<'\\'<<char(c);
        else if(c<32) {
            char escaped[7]; std::snprintf(escaped,sizeof(escaped),"\\u%04x",unsigned(c)); f<<escaped;
        } else f<<char(c);
    }
    f<<"\"}\n";
}
}

template<class Cache>
void equilibrate_column_impl(RadField& rad, const RTGrids& g, const ModelParams& p,
                                const Cache& cache,
                                std::vector<std::vector<LineRec>>& lines, int iteration,
                                const std::vector<double>* cell_incoming=nullptr)
{
    if(p.test_rt) throw std::runtime_error("Thermal roots require production mode");
    if(cache.NE!=g.NE || cache.NT<2)
        throw std::runtime_error("thermal balance: incompatible kernel energy/temperature grid");

    // save other depths' line opacity when the current depth of line is keeping update.
    const FrozenLineColumn column(g,lines);
    const int nd=g.ND_MID;
    if(cell_incoming && cell_incoming->size()!=size_t(nd)*g.NA*g.NE)
        throw std::runtime_error("thermal roots: incoming ray dimensions mismatch");

    // energy weights for integral
    auto w=dao_thermal::energy_weights(g.NE,g.ene);
    double incoming_flux=0;
    if(cell_incoming) {
        for(int m=0;m<g.NA;++m) for(int e=0;e<g.NE;++e)
            incoming_flux+=0.5*phys::four_pi*g.wt[m]*std::abs(g.mu[m])*w[e]*
                (g.mu[m]<0 ? rad.Inu_top[m][e] : rad.Inu_bottom[m][e]);
        if(!std::isfinite(incoming_flux) || incoming_flux<=0)
            throw std::runtime_error("thermal roots: invalid incoming boundary flux");
    }

    std::vector<dao_thermal::TemperatureRoot> roots;
    std::vector<double> trial(nd);
    std::vector<unsigned char> active(nd,1);
    for(int d=0;d<nd;++d) {
        roots.emplace_back(rad.T_K[d],cache.T_grid[0],cache.T_grid[cache.NT-1]);
        trial[d]=rad.T_K[d];
    }
    // Diagnostic log file for the temperature root-finding process
    std::ofstream trace(std::string(p.run_dir)+"/thermal_roots_iter"+std::to_string(iteration)+".dat");
    if(!trace) throw std::runtime_error("Cannot write temperature root trace");
    trace<<std::setprecision(15)<<"# evaluation depth T residual relative_residual\n";
    trace<<"# context "<<(cell_incoming ? "fixed incoming cell rays; locally solved radiation response" : "fixed J diagnostic")<<'\n';

    // Iterate until all depths have converged to a root of the thermal balance equation.
    for(int evaluation=0;evaluation<60;++evaluation) {

        dao_log::detail("  [thermal roots] outer=%d evaluation=%d active=%d/%d\n",
                        iteration,evaluation,int(std::count(active.begin(),active.end(),1)),nd);

        // run cloudy
        run_cloudy_depths(rad,g,p,iteration,lines,0,trial.data(),active.data());
        bool done=true;

        //
        for(int d=0;d<nd;++d) if(active[d]) {
            // Cloudy has just updated n_e at this trial temperature.
            const Scattering scattering=scattering_at(cache,rad.T_K[d],rad.n_e[d],g);
            std::copy(scattering.opacity().begin(),scattering.opacity().end(),rad.ksct[d]);
            apply_frozen_line_escape(d,rad,g,lines[d],column);

            // Production includes this cell's radiation response to its trial
            // opacity/emissivity. Atomic populations still use the previous
            // global J, which is refreshed by every outer RT solve.
            double residual,relative;
            if(cell_incoming) {
                std::vector<double> angular;
                if(p.angsca) {
                    angular.resize(g.NA*g.NE);
                    for(int m=0;m<g.NA;++m)
                        std::copy(rad.Inu[d][m],rad.Inu[d][m]+g.NE,angular.data()+m*g.NE);
                }
                const auto response=dao_thermal::local_response(g.NE,g.NA,w.data(),g.mu,g.wt,g.dr[d],
                    rad.kabs[d],rad.jnu[d],scattering,cell_incoming->data()+long(d)*g.NA*g.NE,
                    p.angsca ? angular.data() : rad.J0[d],p.maxiter);
                const Budget physical=dao_thermal::budget(g.NE,w.data(),response.mean.data(),rad.kabs[d],
                    rad.jnu[d],scattering,p.angsca ? response.intensity.data() : nullptr);
                // Thick absorption: exchange scale avoids cancellation of
                // trapped terms. Thin scattering: physical thermal scale
                // avoids accepting a small fraction of a large through-flux
                // while the actual heating/cooling remains unbalanced.
                residual=response.residual;
                // Also allocate Fin/ND to each cell's integrated error scale.
                // With the 1e-3 root tolerance, sum |R_trial|*dr < 1e-3*Fin.
                // Local relative tolerances alone need not meet the column
                // tolerance when large heating/cooling terms nearly cancel.
                const double column_scale=incoming_flux/(nd*g.dr[d]);
                relative=residual/std::max(std::min({response.scale,physical.scale(),column_scale}),1e-100);
            } else {
                const Budget b=cell_budget(rad,g,d,w,scattering);
                residual=b.residual(); relative=b.relative();
            }
            trace<<evaluation<<' '<<d<<' '<<rad.T_K[d]<<' '<<residual<<' '<<relative<<'\n';
            trace.flush();
            if(std::abs(relative)<1e-3) active[d]=0;
            else {
                roots[d].record(rad.T_K[d],residual);
                try { trial[d]=roots[d].next(); }
                catch(const std::exception& e) {
                    throw std::runtime_error("Depth "+std::to_string(d)+": "+e.what());
                }
                done=false;
            }
        }
        if(done) return;
    }
    throw std::runtime_error("thermal balance: temperature roots did not converge");
}

template<class Cache>
void run_thermal_production_impl(RadField& rad, const RTGrids& g, const ModelParams& p,
                                 const Cache& cache, bool supplied_initial_state)
{
    int iteration=0;
    const auto start=std::chrono::steady_clock::now();
    status(p,"running",iteration);
    try {
        clear_ion_fraction_output(p);
        dao_log::info("Thermal balance: initializing Cloudy and radiation field...\n");
        dao_log::detail("DAO RT energy balance; volume-averaged constant cells; %s Compton\n",p.angsca ? "angle-dependent" : "angle-mean");
        std::vector<double> top(g.NE),bottom(g.NE);
        compute_boundary_illumination(g.NE,p.i_incidence,g,rad,top.data(),bottom.data());
        for(int d=0;d<g.ND_MID;++d) {
            for(int e=0;e<g.NE;++e) {
                if(!supplied_initial_state) {
                    // I_corona stores F_E. Derive
                    // the initial mean from the actual boundary intensities.
                    rad.J0[d][e]=0;
                    for(int m=0;m<g.NA;++m) {
                        rad.Inu[d][m][e]=g.mu[m]<0 ? incident_top_intensity(p.i_incidence,m,top[e]) : bottom[e];
                        rad.J0[d][e]+=0.5*g.wt[m]*rad.Inu[d][m][e];
                    }
                }
                if(!std::isfinite(rad.J0[d][e]) || rad.J0[d][e]<0)
                    throw std::runtime_error("thermal balance: invalid initial radiation field");
                // Initial moments supply an isotropic guess only; the full
                // angular field is subsequently solved and retained for roots.
                if(supplied_initial_state) for(int m=0;m<g.NA;++m) rad.Inu[d][m][e]=rad.J0[d][e];
            }
            rad.line_heat[d]=0;
        }
        rad.compute_ionization_parameter(p.nh);
        std::vector<std::vector<LineRec>> lines(g.ND_MID);
        // Free-T Cloudy is only an initial guess, never the convergence equation.
        // A resolution test may supply an interpolated initial T/J; its atomic
        // state is still recalculated and must pass all the same final checks.
        std::vector<double> initialT;
        if(supplied_initial_state) initialT.assign(rad.T_K,rad.T_K+g.ND_MID);
        run_cloudy_depths(rad,g,p,0,lines,0,supplied_initial_state ? initialT.data() : nullptr);
        // Build a consistent transfer solution and its incoming face rays
        // before the first temperature update, including diffuse emission.
        apply_line_escape(rad,g,lines,p,0);
        std::fill(rad.line_heat,rad.line_heat+g.ND_MID,0);
        std::vector<double> cell_incoming;
        compton_rt_solve(rad,g,p,make_scattering_column(rad,g,cache),&cell_incoming);
        rad.compute_moments(); rad.compute_ionization_parameter(p.nh);
        auto w=dao_thermal::energy_weights(g.NE,g.ene);
        std::ofstream history(std::string(p.run_dir)+"/thermal_history.dat");
        if(!history) throw std::runtime_error("Cannot write thermal history");
        history<<std::setprecision(15)
               <<"# iteration max_dlogT max_relative_L1_d"<<(p.angsca ? "I" : "J")
               <<" max_local_R abs_volume_R/Fin Fout/Fin identity_error\n";

        // main loop
        int successive=0;
        for(iteration=1;iteration<=100;++iteration) {
            const auto iteration_start=std::chrono::steady_clock::now();
            std::vector<double> oldT(rad.T_K,rad.T_K+g.ND_MID), oldJ(long(g.ND_MID)*g.NE);
            for(int d=0;d<g.ND_MID;++d) std::copy(rad.J0[d],rad.J0[d]+g.NE,oldJ.data()+long(d)*g.NE);
            std::vector<double> oldI;
            if(p.angsca) {
                oldI.resize(long(g.ND_MID)*g.NA*g.NE);
                for(int d=0;d<g.ND_MID;++d) for(int m=0;m<g.NA;++m)
                    std::copy(rad.Inu[d][m],rad.Inu[d][m]+g.NE,oldI.data()+(long(d)*g.NA+m)*g.NE);
            }
            equilibrate_column_impl(rad,g,p,cache,lines,iteration,&cell_incoming);
            std::vector<double> damped(g.ND_MID);
            double dT=0;
            for(int d=0;d<g.ND_MID;++d) {
                const double step=std::clamp(0.5*std::log(rad.T_K[d]/oldT[d]),-0.3,0.3);
                damped[d]=oldT[d]*std::exp(step);
                dT=std::max(dT,std::abs(step)/std::log(10.0));
            }
            // Rebuild atomic data and escape at the ACCEPTED temperature/column.
            // Never pass emissivities from an undamped trial to the RT solver.
            run_cloudy_depths(rad,g,p,iteration,lines,0,damped.data());
            auto ops=make_scattering_column(rad,g,cache);
            apply_line_escape(rad,g,lines,p,iteration);
            std::fill(rad.line_heat,rad.line_heat+g.ND_MID,0);
            compton_rt_solve(rad,g,p,ops,&cell_incoming);
            rad.compute_moments(); rad.compute_ionization_parameter(p.nh);
            double dJ=0;
            for(int d=0;d<g.ND_MID;++d) {
                double change=0,scale=0;
                if(p.angsca) {
                    for(int m=0;m<g.NA;++m) for(int e=0;e<g.NE;++e) {
                        const double old=oldI[(long(d)*g.NA+m)*g.NE+e],now=rad.Inu[d][m][e];
                        const double weight=0.5*g.wt[m]*w[e];
                        change+=weight*std::abs(now-old);scale+=weight*std::max(now,old);
                    }
                } else for(int e=0;e<g.NE;++e) {
                    const double old=oldJ[long(d)*g.NE+e],now=rad.J0[d][e];
                    change+=w[e]*std::abs(now-old); scale+=w[e]*std::max(now,old);
                }
                dJ=std::max(dJ,change/std::max(scale,1e-100));
            }
            const ColumnBudget b=column_budget(rad,g,ops,std::string(p.run_dir)+"/thermal_budget_iter"+std::to_string(iteration)+".dat");
            if(std::abs(b.identity_error())>1e-4)
                throw std::runtime_error("thermal balance: boundary/volume energy identity failed");
            history<<iteration<<' '<<dT<<' '<<dJ<<' '<<b.max_local<<' '<<b.volume_abs/b.incoming
                   <<' '<<b.outgoing/b.incoming<<' '<<b.identity_error()<<std::endl;
            save_results(rad,g,p,iteration);
            // Require all five criteria for three consecutive outer iterations:
            // 1. dT: max_d |log10(T_new/T_old)| < 3e-3 dex, after damping.
            // 2. dJ: maximum cell-wise weighted relative L1 field change < 3e-3.
            //    Compare I over energy and angle for angsca=1; J over energy for angsca=0.
            // 3. max_local: max_d |R_d| / (absorption_d + emission_d + |Q_Compton,d|) < 3e-3.
            // 4. volume_abs/Fin: sum_d |R_d| * dr_d / Fin < 1e-2 (no sign cancellation).
            // 5. |Fout - Fin| / Fin < 1e-2; both fluxes include both slab faces.
            const bool converged=dT<3e-3 && dJ<3e-3 && b.max_local<3e-3 &&
                                  b.volume_abs/b.incoming<ModelParams::thermal_column_tolerance &&
                                  std::abs(b.closure())<ModelParams::thermal_column_tolerance;
            successive=converged ? successive+1 : 0;
            dao_log::info("Iter %03d | dlogT=%.2e d%s=%.2e maxR=%.2e absR/Fin=%.2e Fout/Fin=%.6f | %.1fs | pass=%d/3\n",
                          iteration,dT,p.angsca ? "I" : "J",dJ,b.max_local,b.volume_abs/b.incoming,
                          b.outgoing/b.incoming,
                          std::chrono::duration<double>(std::chrono::steady_clock::now()-iteration_start).count(),successive);
            dao_log::detail("  [energy identity] outer=%d error=%.6e\n",iteration,b.identity_error());
            status(p,"running",iteration);
            if(successive==3) {
                apply_line_escape(rad,g,lines,p,iteration,true);
                save_ion_fractions(rad,g,p,iteration);
                status(p,"converged",iteration);
                dao_log::info("Converged after %d iterations: Fout/Fin=%.9f, energy error=%.4f%% (both faces).\n",
                              iteration,b.outgoing/b.incoming,100*std::abs(b.closure()));
                return;
            }
        }
        throw std::runtime_error("thermal balance: outer iteration limit; no converged spectrum");
    } catch(const std::exception& e) {
        dao_log::error("Thermal balance failed at iteration %d after %.1fs: %s\n",iteration,
                       std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count(),e.what());
        status(p,"failed",iteration,e.what()); throw;
    }
}

void thermal_equilibrate_column(RadField& r,const RTGrids& g,const ModelParams& p,
                                const avgKernelCache& c,std::vector<std::vector<LineRec>>& l,int n)
{ equilibrate_column_impl(r,g,p,c,l,n); }
void thermal_equilibrate_column(RadField& r,const RTGrids& g,const ModelParams& p,
                                const KernelCache& c,std::vector<std::vector<LineRec>>& l,int n)
{ equilibrate_column_impl(r,g,p,c,l,n); }
void run_thermal_production(RadField& r,const RTGrids& g,const ModelParams& p,
                            const avgKernelCache& c,bool initial)
{ run_thermal_production_impl(r,g,p,c,initial); }
void run_thermal_production(RadField& r,const RTGrids& g,const ModelParams& p,
                            const KernelCache& c,bool initial)
{ run_thermal_production_impl(r,g,p,c,initial); }
