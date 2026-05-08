// ============================================================
// maindaocl — Compton scattering RT solver using Cloudy
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
		RTGrids g;
		g.init_angle();
		g.init_depth(RTGrids::TAU_MIN, RTGrids::TAU_MAX, par.nh);
		bootstrap_cloudy_energy_grid(g, EGRID_FILE);
		snap_incidence(par, g);

		// --- 4. Calculate initial radiation field ---
		RadField rad(g);
		rad.allocate();
		rad.illum.compute(par);

		// --- 5. Precompute kernel ---
		KernelCache kcache;
		kcache.init(g.NE, g.ene, g.NA, g.mu, g.wt, par.ktype);
		
		// --- 6. Dispatch ---
		if (par.test_rt)
			run_test_rt(rad, g, par, kcache);
		else
			run_production(rad, g, par, kcache);

		kcache.free_memory();

		cdEXIT(exit_status);
	}
	CLOUDY_CATCH_ALL(exit_status);

	cdPrepareExit(exit_status);
	return exit_status;
}
