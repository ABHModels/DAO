#include "thermal_production.h"
#include "cell_transfer.h"
#include "compton_rt.h"
#include "cloudy_depth.h"
#include "save_results.h"
#include "rt_parallel.h"
#include "run_log.h"
#include "thermal_restart.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <filesystem>
#include <memory>
#include <limits>
#include <fcntl.h>
#include <sys/file.h>
#include <unistd.h>

namespace {
using dao_thermal::Scattering;
using dao_thermal::Budget;

class ResumeLock {
    int descriptor=-1;
public:
    explicit ResumeLock(const ModelParams& p) {
        const auto path=std::string(p.run_dir)+"/thermal_resume.lock";
        descriptor=::open(path.c_str(),O_CREAT|O_RDWR,0644);
        if(descriptor<0) throw std::runtime_error("Cannot open resume lock");
        if(flock(descriptor,LOCK_EX|LOCK_NB)<0) {
            ::close(descriptor); descriptor=-1;
            throw std::runtime_error("Another continuation is already using this run directory");
        }
    }
    ~ResumeLock() { if(descriptor>=0) ::close(descriptor); }
};

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

/*
 * THERMAL ROOT SEARCH
 * Per-depth flow;
 *     Hold incoming rays and the line column fixed.
 *     Initialize trial T and retain each depth's sync mode.
 *                            |
 *                            v
 *     +--------------------------------------------------+
 * [A] | Run Cloudy at trial T with the current J.         |
 *     | Update opacity, emission and Compton scattering. |
 *     | Solve local transfer; obtain J_new and residual. |
 *     +--------------------------------------------------+
 *                            |
 *                  Synchronization enabled?
 *                   /                    \
 *                 NO                     YES
 *                  |                      |
 *       Save the best trial.     Compare J_new with J.
 *                  |            Replace J with J_new.
 *                  |                      |
 *                  |              Field change < 1e-5?
 *                  |               /             \
 *                  |             YES              NO
 *                  |              |                |
 *                  |              |       Repeat [A] at SAME T.
 *                  |              |
 *                  +--------------+
 *                            |
 *                            v
 *     +--------------------------------------------------+
 * [B] | Is |normalized thermal residual| < 1e-3?          |
 *     +--------------------------------------------------+
 *                   /                    \
 *                 YES                     NO
 *                  |                      |
 *       Accept this depth.       Record T and residual.
 *                  |                      |
 *       Return when ALL          Normal mode, and
 *       depths are accepted.     next T outside
 *                                [1e1, 5e9] K?
 *                                  /             \
 *                                YES              NO
 *                                 |                |
 *                       Enable persistent     Choose next T
 *                       synchronization.      from root search.
 *                       Clear root samples.        |
 *                       Restart at best T.         |
 *                                 |                |
 *                                 +----------------+
 *                                          |
 *                                  Repeat [A] at NEW T.
 *
 * Sync test: field change < 1e-5 after at least two Cloudy/transfer updates.
 *            Thermal-root and outer convergence tests remain separate.
 *             These is a
 *             new function of code after version v1.0.1. I
 *             add this because we find the temperature is unstable
 *             when the ionization is very low (logxi<0).
 * Stable bracket: net heating below, net cooling above.
 * Exceeding a search or synchronization limit raises an error.
 */

template<class Cache>
void equilibrate_column_impl(RadField& rad, const RTGrids& g, const ModelParams& p,
                                const Cache& cache,
                                std::vector<std::vector<LineRec>>& lines, int iteration,
                                const std::vector<double>* cell_incoming=nullptr,
                                std::vector<unsigned char>* synchronized_depths=nullptr)
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
    std::vector<unsigned char> local_modes(nd,0);
    auto& coupled=synchronized_depths ? *synchronized_depths : local_modes;
    if(coupled.size()!=size_t(nd)) throw std::runtime_error("thermal roots: invalid synchronization mask");
    // Trigger synchronization near Cloudy's temperature limits; these do not limit gas T.
    constexpr double sync_temperature_lo=1e1, sync_temperature_hi=5e9;
    std::vector<double> best_relative(nd,std::numeric_limits<double>::infinity()),best_T(nd);
    std::vector<dao_thermal::CellResponse> best_response(nd);
    const double minimum_temperature=cloudy_temperature_minimum();
    const double maximum_temperature=cloudy_temperature_maximum();
    for(int d=0;d<nd;++d) {
        // Search Cloudy's gas-temperature range; the kernel freezes at either edge.
        trial[d]=std::clamp(rad.T_K[d],minimum_temperature,maximum_temperature);
        roots.emplace_back(trial[d],minimum_temperature,maximum_temperature,coupled[d]!=0);
    }
    // Diagnostic log file for the temperature root-finding process
    std::ofstream trace(std::string(p.run_dir)+"/thermal_roots_iter"+std::to_string(iteration)+".dat");
    if(!trace) throw std::runtime_error("Cannot write temperature root trace");
    trace<<std::setprecision(15)<<"# evaluation depth T residual relative_residual\n";
    trace<<"# context "<<(cell_incoming ? "fixed incoming cell rays; locally solved radiation response" : "fixed J diagnostic")<<'\n';

    std::ofstream sync_trace;
    if(p.verbose) {
        sync_trace.open(std::string(p.run_dir)+"/thermal_sync_iter"+std::to_string(iteration)+".dat");
        if(!sync_trace) throw std::runtime_error("Cannot write atomic/radiation synchronization trace");
        sync_trace<<std::setprecision(15)<<"# evaluation subiteration depth T relative_dJ relative_residual\n";
    }

    // Find temperatures
    // Allow a full ordinary search followed by a full synchronized restart.
    for(int evaluation=0;evaluation<2*60;++evaluation) {

        dao_log::detail("  [thermal roots] outer=%d evaluation=%d active=%d/%d\n",
                        iteration,evaluation,int(std::count(active.begin(),active.end(),1)),nd);

        // Normal depths retain the original fixed-J search. Only failed depths
        // synchronize their atomic state and local radiation at each trial T.
        // coupled[d] is true when the depth is synchronizing its atomic/radiation state.
        // pending[d] is true when the depth has not yet converged at the current trial T.
        std::vector<unsigned char> pending=active;
        std::vector<double> residuals(nd),relatives(nd);
        for(int sync=0;sync<100;++sync) {
            // Run Cloudy at the trial T.
            run_cloudy_depths(rad,g,p,iteration,lines,0,trial.data(),pending.data());
            for(int d=0;d<nd;++d) if(pending[d]) {
                // Use the same 1.21*nH scattering density in local roots and full-slab RT.
                const Scattering scattering=scattering_at(cache,rad.T_K[d],scattering_target_density(rad,p,d),g);
                std::copy(scattering.opacity().begin(),scattering.opacity().end(),rad.ksct[d]);
                // apply surviving factor correction for every lines.
                apply_frozen_line_escape(d,rad,g,lines[d],column);

                double residual,relative; // Thermal equilibirum residual and its normalized value.
                bool synchronized=true; // default to true for non-synchronizing depths; only failed depths are set false.
                if(cell_incoming) {
                    std::vector<double> angular;
                    if(p.angsca) {
                        angular.resize(g.NA*g.NE);
                        for(int m=0;m<g.NA;++m)
                            std::copy(rad.Inu[d][m],rad.Inu[d][m]+g.NE,angular.data()+m*g.NE);
                    }
                    // solve depth-cell average radiation field with fixed incoming rays.
                    auto response=dao_thermal::local_response(g.NE,g.NA,w.data(),g.mu,g.wt,g.dr[d],
                        rad.kabs[d],rad.jnu[d],scattering,cell_incoming->data()+long(d)*g.NA*g.NE,
                        p.angsca ? angular.data() : rad.J0[d],p.maxiter);
                    // compute thermal equilibrium
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
                    const double scale=std::max(std::min({response.scale,physical.scale(),column_scale}),1e-100);
                    relative=residual/scale;


                    // A new function of code after version v1.0.1.
                    // This help the thermal structure to converge when the ionization is very low (logxi<0).
                    // I'd like to compute neutral ionized disk, why don't?
                    // but it increase the computation time.
                    if(coupled[d]) {
                        double difference=0,field_scale=0;
                        // difference of radiation quantities
                        for(int e=0;e<g.NE;++e) {
                            const double delta=std::abs(response.mean[e]-rad.J0[d][e]);
                            difference+=w[e]*delta;
                            field_scale+=w[e]*std::max(response.mean[e],rad.J0[d][e]);
                        }
                        // Measure the energy-integrated relative change in the radiation field.
                        const double field_change=difference/std::max(field_scale,1e-100);
                        // Stabilize the local field before continuing the temperature search.
                        // Thermal equilibrium is checked separately below and after full-slab RT.
                        synchronized=sync>0 && field_change<1e-5;
                        if(sync_trace.is_open())
                            sync_trace<<evaluation<<' '<<sync<<' '<<d<<' '<<rad.T_K[d]<<' '
                                      <<field_change<<' '<<relative<<'\n';
                        // Store the updated intensity for the next Cloudy call.
                        std::copy(response.mean.begin(),response.mean.end(),rad.J0[d]);
                        for(int m=0;m<g.NA;++m)
                            std::copy(response.intensity.data()+m*g.NE,
                                      response.intensity.data()+(m+1)*g.NE,rad.Inu[d][m]);
                    } else if(std::abs(relative)<best_relative[d]) {
                        best_relative[d]=std::abs(relative);
                        best_T[d]=rad.T_K[d];
                        best_response[d]=std::move(response);
                    }
                } else {
                    const Budget b=cell_budget(rad,g,d,w,scattering);
                    residual=b.residual(); relative=b.relative();
                }
                residuals[d]=residual; relatives[d]=relative;
                if(synchronized) pending[d]=0;
            }
            if(cell_incoming) {
                rad.compute_ionization_parameter(p.nh);
            }
            if(sync_trace.is_open()) sync_trace.flush();
            if(std::none_of(pending.begin(),pending.end(),[](unsigned char v){return v!=0;})) break;
            if(sync==99) {
                const int depth=int(std::find(pending.begin(),pending.end(),1)-pending.begin());
                throw std::runtime_error("thermal balance: atomic/radiation synchronization failed at depth "+std::to_string(depth));
            }
        }
        bool done=true;
        for(int d=0;d<nd;++d) if(active[d]) {
            const double residual=residuals[d],relative=relatives[d];
            trace<<evaluation<<' '<<d<<' '<<rad.T_K[d]<<' '<<residual<<' '<<relative<<'\n';
            trace.flush();
            if(std::abs(relative)<1e-3) active[d]=0;
            else {
                roots[d].record(rad.T_K[d],residual);
                double next_T;
                try { next_T=roots[d].next(); }
                catch(const std::exception& e) {
                    throw std::runtime_error("Depth "+std::to_string(d)+": "+e.what());
                }
                // Intercept an out-of-range trial before sending it to Cloudy.
                if(cell_incoming && !coupled[d] &&
                   (next_T<sync_temperature_lo || next_T>sync_temperature_hi)) {
                    const auto& response=best_response[d];
                    std::copy(response.mean.begin(),response.mean.end(),rad.J0[d]);
                    for(int m=0;m<g.NA;++m)
                        std::copy(response.intensity.data()+m*g.NE,
                                  response.intensity.data()+(m+1)*g.NE,rad.Inu[d][m]);
                    rad.compute_ionization_parameter(p.nh);
                    // New residual definition: discard every fixed-J sample.
                    coupled[d]=1;
                    roots[d]=dao_thermal::TemperatureRoot(best_T[d],minimum_temperature,maximum_temperature,true);
                    trial[d]=best_T[d];
                    trace<<"# activate_sync evaluation="<<evaluation<<" depth="<<d
                         <<" proposed_T="<<next_T<<" restart_T="<<best_T[d]<<'\n'; trace.flush();
                    if(p.verbose) {
                        std::ofstream modes(std::string(p.run_dir)+"/thermal_sync_depths.dat",std::ios::app);
                        if(!modes) throw std::runtime_error("Cannot write synchronization activation trace");
                        modes<<iteration<<' '<<d<<'\n'; modes.flush();
                    }
                    dao_log::info("Thermal search: depth %d, proposed T=%.6g K outside [%.6g, %.6g] K; enable persistent local synchronization at T=%.6g K.\n",
                                  d+1,next_T,sync_temperature_lo,sync_temperature_hi,best_T[d]);
                    done=false;
                    continue;
                }
                trial[d]=next_T;
                done=false;
            }
        }
        if(done) return;
    }
    throw std::runtime_error("thermal balance: temperature roots did not converge");
}

// outer iterations
template<class Cache>
void run_thermal_production_impl(RadField& rad, const RTGrids& g, const ModelParams& p,
                                 const Cache& cache, bool supplied_initial_state)
{
    int iteration=p.resume_iteration;
    const auto start=std::chrono::steady_clock::now();
    std::unique_ptr<ResumeLock> resume_lock;
    std::vector<unsigned char> synchronized_depths(g.ND_MID,0);

    // resume from a saved T/J grid, but still rebuild the atomic state at fixed T and solve the angular RT once.
    if(p.resume_iteration) {
        resume_lock=std::make_unique<ResumeLock>(p);
        load_thermal_resume(rad,g,p);
        const auto modes_path=std::filesystem::path(p.run_dir)/
            ("thermal_sync_modes_iter"+std::to_string(iteration)+".dat");
        // Older checkpoints have no synchronization mask; start in normal mode.
        if(std::filesystem::exists(modes_path)) {
            std::ifstream modes(modes_path);
            if(!modes) throw std::runtime_error("Cannot read saved synchronized depth mask");
            for(int d=0;d<g.ND_MID;++d) {
                int saved_depth=-1,mode=-1;
                if(!(modes>>saved_depth>>mode) || saved_depth!=d || (mode!=0 && mode!=1))
                    throw std::runtime_error("Invalid saved synchronized depth mask");
                synchronized_depths[d]=mode;
            }
            std::string extra;
            if(modes>>extra || modes.bad())
                throw std::runtime_error("Invalid saved synchronized depth mask");
        }
        supplied_initial_state=true;
        const auto directory=std::filesystem::path(p.run_dir);
        const auto previous=directory/("thermal_status_before_resume_iter"+std::to_string(iteration)+".json");
        if(std::filesystem::exists(directory/"thermal_status.json") && !std::filesystem::exists(previous))
            std::filesystem::copy_file(directory/"thermal_status.json",previous);
        std::ofstream record(directory/("resume_iter"+std::to_string(iteration)+".json"));
        if(!record) throw std::runtime_error("Cannot write resume record");
        record<<"{\"resume_iteration\":"<<iteration<<",\"max_outer_iterations\":"<<p.max_outer_iterations
              <<",\"method\":\"saved_T_J_fixed_T_Cloudy_then_RT\",\"thermal_update\":\"full\",\"convergence_streak_reset\":true}\n";
        dao_log::info("Continuing from saved T/J at iteration %d; rebuilding atomic data at fixed T and solving angular RT once.\n",iteration);
    }
    status(p,"running",iteration);

    // solution
    try {
        clear_ion_fraction_output(p);
        dao_log::info("Thermal balance: initializing Cloudy and radiation field...\n");
        dao_log::detail("DAO RT energy balance; volume-averaged constant cells; %s Compton\n",p.angsca ? "angle-dependent" : "angle-mean");
        std::vector<double> top(g.NE),bottom(g.NE);
        /*
            Initialization.
        */
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
            // this quantity is never used now.
            rad.line_heat[d]=0;
        }

        // ionization parameter from the initial radiation field and the hydrogen density.
        rad.compute_ionization_parameter(p.nh);
        std::vector<std::vector<LineRec>> lines(g.ND_MID);
        // Free-T Cloudy is only an initial guess, never the convergence equation.
        // A resolution test may supply an interpolated initial T/J; its atomic
        // state is still recalculated and must pass all the same final checks.
        std::vector<double> initialT;
        if(supplied_initial_state) initialT.assign(rad.T_K,rad.T_K+g.ND_MID);
        // Frist guess of radiation quantities, thermal balance is solved at each cell by cloudy.
        run_cloudy_depths(rad,g,p,iteration,lines,0,supplied_initial_state ? initialT.data() : nullptr);
        // Build a consistent transfer solution and its incoming face rays
        // before the first temperature update, including diffuse emission.
        // This reconstruction pass does not request diagnostic file output.
        apply_line_escape(rad,g,lines,p,iteration);
        std::fill(rad.line_heat,rad.line_heat+g.ND_MID,0);
        std::vector<double> cell_incoming;
        compton_rt_solve(rad,g,p,make_scattering_column(rad,g,cache,p),&cell_incoming);
        rad.compute_moments(); rad.compute_ionization_parameter(p.nh);
        auto w=dao_thermal::energy_weights(g.NE,g.ene);
        std::ofstream history(std::string(p.run_dir)+"/thermal_history.dat",
                              p.resume_iteration ? std::ios::app : std::ios::out);
        if(!history) throw std::runtime_error("Cannot write thermal history");
        history<<std::setprecision(15);
        if(p.resume_iteration)
            history<<"# resumed from iteration "<<p.resume_iteration<<" using saved T/J; fixed-T atomic and RT reconstruction\n";
        else
            history<<"# iteration max_dlogT max_relative_L1_d"<<(p.angsca ? "I" : "J")
                   <<" max_local_R abs_volume_R/Fin Fout/Fin identity_error\n";
        /*
            End of initialization.
        */
        /*
            Main loop.
        */
        int successive=0;
        dao_log::info("Outer temperature update: full target; cap=%d iterations\n",p.max_outer_iterations);
        for(;iteration<p.max_outer_iterations;) {
            ++iteration;
            const auto iteration_start=std::chrono::steady_clock::now();
            std::vector<double> oldT(rad.T_K,rad.T_K+g.ND_MID), oldJ(long(g.ND_MID)*g.NE);
            for(int d=0;d<g.ND_MID;++d) std::copy(rad.J0[d],rad.J0[d]+g.NE,oldJ.data()+long(d)*g.NE);
            std::vector<double> oldI;
            if(p.angsca) {
                oldI.resize(long(g.ND_MID)*g.NA*g.NE);
                for(int d=0;d<g.ND_MID;++d) for(int m=0;m<g.NA;++m)
                    std::copy(rad.Inu[d][m],rad.Inu[d][m]+g.NE,oldI.data()+(long(d)*g.NA+m)*g.NE);
            }
            // inner loop (temperature finding at each depths)
            equilibrate_column_impl(rad,g,p,cache,lines,iteration,&cell_incoming,&synchronized_depths);
            const std::vector<double> targetT(rad.T_K,rad.T_K+g.ND_MID);
            double dT=0;
            for(int d=0;d<g.ND_MID;++d) {
                dT=std::max(dT,std::abs(std::log(targetT[d]/oldT[d]))/std::log(10.0));
            }
            // Rebuild atomic data at the full target temperature before transfer.
            // with our new temperature
            run_cloudy_depths(rad,g,p,iteration,lines,0,targetT.data());
            auto ops=make_scattering_column(rad,g,cache,p);
            apply_line_escape(rad,g,lines,p,iteration);
            std::fill(rad.line_heat,rad.line_heat+g.ND_MID,0);
            // Transfer along the whole slab
            compton_rt_solve(rad,g,p,ops,&cell_incoming);
            rad.compute_moments(); rad.compute_ionization_parameter(p.nh);
            // Compute the maximum relative change in the radiation field, either
            // the angular intensity (angsca=1) or the mean intensity (angsca=0).
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
            std::ofstream modes(std::string(p.run_dir)+"/thermal_sync_modes_iter"+std::to_string(iteration)+".dat");
            if(!modes) throw std::runtime_error("Cannot save synchronized depth mask");
            for(int d=0;d<g.ND_MID;++d) modes<<d<<' '<<int(synchronized_depths[d])<<'\n';
            modes.close();
            if(!modes) throw std::runtime_error("Cannot save synchronized depth mask");
            // Require all five criteria for three consecutive outer iterations:
            // 1. dT: max_d |log10(T_new/T_old)| < 3e-3 dex, for the adopted temperature.
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
