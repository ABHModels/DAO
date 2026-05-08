#include "production.h"
#include "cddefines.h"
#include "cddrive.h"
#include "cloudy_interface.h"
#include "compton_cross_section.h"
#include "compton_rt.h"
#include "save_results.h"
#include "constants.h"
#include <cmath>
#include <cstdio>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include "opacity.h"
#include "dense.h"                                                                                                                                                            
#include "elementnames.h"                                                                                                                                                     
#include "abund.h"     
#include "lines.h"   
#include "yield.h"                                                                                                                                                                             
#include "ionbal.h"

void run_production(RadField& rad, const RTGrids& g,
                    const ModelParams& par,
                    KernelCache& kcache)
{
	// --- Cloudy setup ---
	CloudyInput cl_in(g);
	cl_in.init(par);

	// Runtime-sized scratch arrays for the convergence check.
	// H_trap_* are unused in v2 but required by the updated
	// check_convergence signature.
	double* T_old      = new double[g.ND_MID]();
	double* xi_old     = new double[g.ND_MID];
	for (int id = 0; id < g.ND_MID; ++id)
		xi_old[id] = par.zeta;

	for (int id = 0; id < g.ND_MID; ++id)
	{
		for (int ie = 0; ie < g.NE; ++ie){
			rad.J0[id][ie] = rad.illum.I_corona[ie]+rad.illum.I_disk[ie];
		}
		rad.n_e[id] = pow(10,par.nh)*1.21;
	}
	rad.compute_ionization_parameter(par.nh);
	// --- Outer Cloudy-RT iteration ---
	const double T_tol = 1e-3;
	double mean_dT = 1.0, mean_dXi = 1.0;
	int outer_iter = 0;
	const int max_outer = 15;
	while ((mean_dT > T_tol && mean_dXi > T_tol && outer_iter<max_outer))
	{
		++outer_iter;
		fprintf(stdout, "\n========== OUTER ITERATION %d / %d ==========\n",
		        outer_iter, max_outer);

		// --- Cloudy depth loop ---
		fprintf(stdout, "  [Cloudy] Running %d depth points...\n", g.ND_MID);
		for (int id = 0; id < g.ND_MID; ++id)
		{
			cdInit();
			cdTalk(false);

			cl_in.issue_constant();
			cl_in.issue_depth(id, rad ,g, par);
			
			int rc = cdDrive();
			if (rc)
			{
				fprintf(stderr, "  depth %d/%d: cdDrive failed, skipping.\n",
					id + 1, g.ND_MID);
				continue;
			}

			extract_cloudy_output(id, rad, g, outer_iter);

			compute_compton_opacity(rad.ksct[id], g.NE, g.ene,
			                       rad.T_K[id], rad.n_e[id]);

			if (id == 0 || outer_iter%5==0) {
				save_cloudy_opacity(id, rad, g,par, outer_iter);
				save_line_labels(id, g, par,outer_iter);
			}
				
			fprintf(stdout,
				"  depth %3d/%d  tau=%.3e  logT=%.3f  log(I)=%.3f  "
				"ne/nh=%.3f  H/C=%.3f\n",
				id + 1, g.ND_MID, g.tau_mid[id],
				log10(rad.T_K[id]),
				rad.log_xi[id]+par.nh,
				rad.n_e[id] / pow(10.0, par.nh),
				rad.heating[id] / rad.cooling[id]);
		}

		// --- RT solver ---
		fprintf(stdout, "  [RT] Solving radiative transfer...\n");
		compton_rt_solve(rad, g, par, kcache, par.maxiter);

		rad.compute_moments();
		rad.compute_ionization_parameter(par.nh);
		rad.check_convergence(outer_iter,
		                      T_old, xi_old,
		                      mean_dT, mean_dXi);

		fprintf(stdout, "  Outer iter %d: mean|dT/T|=%.4e  mean|d(log_xi)|=%.4e\n",
			outer_iter, mean_dT, mean_dXi);
		
		save_results(rad, g, par, outer_iter);
	}

	// for (int id = 0; id < g.ND_MID; ++id){
	// 	for (int i = 0; i < g.NE; ++i){
	// 		rad.jnu[id][i] += rad.jnu_line[id][i];
	// 	}
	// }
	// compton_rt_solve(rad, g, par, kcache, par.maxiter);
	// rad.compute_moments();
	// rad.compute_ionization_parameter(par.nh);
	// outer_iter = 0;
	// save_results(rad, g, par, outer_iter);
	
	fprintf(stdout, "\n=== %s at outer iteration %d (mean|dT/T|=%.2e) ===\n",
		(mean_dT <= T_tol) ? "Converged" : "Max iterations reached",
		outer_iter, mean_dT);

	delete[] T_old;
	delete[] xi_old;
}
