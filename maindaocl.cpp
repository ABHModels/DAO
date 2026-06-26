// ============================================================
// maindaocl — Compton scattering RT solver using Cloudy
//
// Author:  Yimin Huang (Fudan University; University of Bristol)
// Contact: huangym23@m.fudan.edu.cn
// License: MIT (see LICENSE). Cloudy and HEASoft/Xspec are separate
//          dependencies under their own licenses.
//
// 1. Read input
// 2. Bootstrap Cloudy to get energy grid
// 3. Calculate grids (angle, depth, energy from Cloudy)
// 4. Calculate initial radiation field
// 5. Precompute kernel
// 6. Dispatch: test_rt or production (Cloudy + RT loop)
// ============================================================

#include "cddefines.h"
#include "cddrive.h"
#include "cloudy_exception.h"
#include "constants.h"
#include "rt_grids.h"
#include "params.h"
#include "radiation.h"
#include "cloudy_interface.h"
#include "compton_kernel.h"
#include "avg_compton_kernel.h"
#include "compton_rt.h"
#include "test_rt.h"
#include "production.h"

int main(int argc, char *argv[])
{
	exit_type exit_status = ES_SUCCESS;

	DEBUG_ENTRY("main()");

	try {
		// --- 1. Read input ---
		ModelParams par = read_params(argc, argv);

		// --- 2. Energy grid: load from file, or bootstrap Cloudy ---
		static const char* EGRID_FILE = "cloudy_energy_grid.dat";
		const bool is_test   = par.test_rt;
		const bool is_angavg = par.angsca;
		const bool is_compps = is_test && strcmp(par.test_mode, "compps") == 0;

		RTGrids g;
		// compps benchmark: use the double-Gauss angle grid so our nodes match
		// Xspec compPS's (QDRGSDO on [0,1]); exact when NA=10. Production keeps
		// the full-sphere Gauss-Legendre grid.
		if (is_compps)
			g.init_angle_double_gauss();
		else
			g.init_angle();
		// compps benchmark mode: total Thomson depth set by -tau so it matches
		// the compps "tau" parameter; otherwise use the default slab depth.
		const double tau_max = is_compps ? par.tau_slab : RTGrids::TAU_MAX;
		g.init_depth(RTGrids::TAU_MIN, tau_max, par.nh);

		// Energy grid: test mode uses a synthetic 1000-bin log grid spanning
		// 0.01-1000 keV (no Cloudy needed); production adopts Cloudy's mesh.
		if (is_test)
			g.init_energy(10.0, 1.0e6, 1000);   // 0.01-1000 keV in eV
		else
			bootstrap_cloudy_energy_grid(g, EGRID_FILE);
		snap_incidence(par, g);

		// --- 4. Calculate initial radiation field ---
		RadField rad(g);
		rad.allocate();
		rad.illum.compute(par);

		// --- 5+6. Precompute kernel and dispatch ---
		// In test mode the slab is isothermal, so the redistribution kernel
		// only needs the single slab temperature (NT=1) rather than the full
		// N_T_CACHE grid used in production.
		if (is_angavg)
		{
			// Angle-resolved (directional) Compton kernel.
			KernelCache kcache;
			if (is_test)
			{
				const double kB_eV  = phys::k_B / phys::eV_to_erg;   // [eV/K]
				double T_slab = par.kT_e * 1.0e3 / kB_eV;   // slab Te [keV] -> [K]
				kcache.init(g.NE, g.ene, g.NA, g.mu, g.wt, par.ktype, 1, &T_slab);
			}
			else
			{
				kcache.init(g.NE, g.ene, g.NA, g.mu, g.wt, par.ktype);
			}

			if (par.test_rt)
				run_test_rt(rad, g, par, kcache);
			else
				run_production(rad, g, par, kcache);

			kcache.free_memory();
		}
		else
		{
			// Angle-mean (isotropic) Compton kernel.
			avgKernelCache kcache;
			if (is_test)
			{
				const double kB_eV  = phys::k_B / phys::eV_to_erg;   // [eV/K]
				double T_slab = par.kT_e * 1.0e3 / kB_eV;   // slab Te [keV] -> [K]
				kcache.init(g.NE, g.ene, g.NA, g.mu, g.wt, par.ktype, 1, &T_slab);
			}
			else
			{
				kcache.init(g.NE, g.ene, g.NA, g.mu, g.wt, par.ktype);
			}

			if (par.test_rt)
				run_test_rt(rad, g, par, kcache);
			else
				run_production(rad, g, par, kcache);

			kcache.free_memory();
		}

		cdEXIT(exit_status);
	}
	CLOUDY_CATCH_ALL(exit_status);

	cdPrepareExit(exit_status);
	return exit_status;
}
