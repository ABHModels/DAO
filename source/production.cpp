#include "production.h"
#include "cloudy_depth.h"
#include "cddefines.h"
#include "cddrive.h"
#include "cloudy_interface.h"
#include "compton_cross_section.h"
#include "compton_rt.h"
#include "save_results.h"
#include "constants.h"
#include <cmath>
#include <cstdio>
#include <vector>
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

template<class Cache>
void run_production(RadField& rad, const RTGrids& g,
                    const ModelParams& par,
                    Cache& kcache)
{
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
			rad.J0[id][ie] = (rad.illum.I_corona[ie]+rad.illum.I_disk[ie]);
		}
		rad.n_e[id] = pow(10,par.nh)*1.21;
	}
	rad.compute_ionization_parameter(par.nh);

	// --- Outer Cloudy-RT iteration ---
	const double cir = 3e-3;
	double max_dT = 1.0, max_dXi = 1.0;
	int outer_iter = 0;
	const int max_outer = 100;
	while (((max_dT > cir || max_dXi > cir) && outer_iter<max_outer))
	{
		++outer_iter;
		fprintf(stdout, "\n========== OUTER ITERATION %d / %d ==========\n",
		        outer_iter, max_outer);

		// --- Cloudy depth loop ---
		// Per-cell line records, consumed by apply_line_escape() after the column
		// is complete (the cumulative two-sided line optical depth needs all cells).
		std::vector<std::vector<LineRec>> line_store(g.ND_MID);
		run_cloudy_depths(rad, g, par, outer_iter, line_store);

		// --- Two-sided cumulative line escape: inject emiss*P into jnu ---
		apply_line_escape(rad, g, line_store, par, outer_iter);

		// --- RT solver ---
		fprintf(stdout, "  [RT] Solving radiative transfer...\n");
		compton_rt_solve(rad, g, par, kcache, par.maxiter);

		// --- update moments---
		rad.compute_moments();

		// --- update ionization parameter---
		rad.compute_ionization_parameter(par.nh);

		// --- check if convergence ---
		rad.check_convergence(outer_iter,
		                      T_old, xi_old,
		                      max_dT, max_dXi);

		fprintf(stdout, "  Outer iter %d: max|d log T|=%.4e  max|d log xi|=%.4e\n",
			outer_iter, max_dT, max_dXi);

		// --- save results of current iteration ---
		save_results(rad, g, par, outer_iter);
	}

	// This loop is used for write the ions' fraction.

	// {
	// 	for (int id = 0; id < g.ND_MID; ++id)
	// 	{
	// 		cdInit();
	// 		cdTalk(false);
	// 		cl_in.issue_constant();
	// 		cl_in.issue_depth_lastest(id, rad ,g, par);
	// 		int rc = cdDrive();
	// 		if (rc)
	// 		{
	// 			fprintf(stderr, "Final output");
	// 		}
	// 	}
	// }
	
	bool converged = (max_dT <= cir && max_dXi <= cir);
	fprintf(stdout, "\n=== %s at outer iteration %d "
	                "(max|d log T|=%.2e, max|d log xi|=%.2e) ===\n",
		converged ? "Converged" : "Max iterations reached",
		outer_iter, max_dT, max_dXi);

	delete[] T_old;
	delete[] xi_old;
}

// Explicit instantiations for both kernel-cache types.
template void run_production<KernelCache>(
	RadField&, const RTGrids&, const ModelParams&, KernelCache&);
template void run_production<avgKernelCache>(
	RadField&, const RTGrids&, const ModelParams&, avgKernelCache&);
