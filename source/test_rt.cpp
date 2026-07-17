#include "test_rt.h"
#include "rt_grids.h"
#include "constants.h"
#include "compton_cross_section.h"
#include "compton_rt.h"
#include "corona_models.h"   // blackbody()

#include <cmath>
#include <cstdio>
#include <cstring>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>

// ============================================================
// Fill data for RT test mode: isothermal pure-scattering slab.
// ============================================================
static void fill_test_data(RadField& rad, const RTGrids& g,
                           double nh, double T_slab)
{
	const double nH  = pow(10.0, nh);
	const double n_e = 1.21 * nH;

	// Compute scattering opacity once (same T at all depths)
	compute_compton_opacity(rad.ksct[0], g.NE, g.ene, T_slab, n_e);

	for (int id = 0; id < g.ND_MID; ++id)
	{
		rad.T_K[id] = T_slab;
		rad.n_e[id] = n_e;

		for (int ie = 0; ie < g.NE; ++ie)
		{
			rad.kabs[id][ie] = 0.0;
			rad.jnu[id][ie]  = 0.0;
		}

		// Copy ksct from depth 0
		if (id > 0)
			memcpy(rad.ksct[id], rad.ksct[0], g.NE * sizeof(double));

		rad.heating[id] = 0.0;
		rad.line_heat[id] = 0.0;
		rad.cooling[id] = 0.0;
	}

	double kb_eV = phys::k_B / phys::eV_to_erg;
	fprintf(stdout, "Test mode: pure scattering atmosphere\n");
	fprintf(stdout, "  T=%.1e K  kT=%.2f eV  n_e=%.3e cm^-3\n",
		T_slab, kb_eV * T_slab, n_e);
}

// ============================================================
// compps benchmark mode
//
// Isothermal, pure-scattering plane-parallel slab illuminated
// from the BOTTOM by an isotropic blackbody seed — the exact
// configuration of Xspec's compps with:
//   geom=1 (slab), cov_frac=1 (clean slab), rel_refl=0 (no reflection),
//   Maxwellian electrons at kTe, seed blackbody at kTbb.
//
// Parameter mapping (compps  ->  our flags):
//   kTe   [keV]  ->  -kT_e   (uniform electron temperature T_K)
//   kTbb  [keV]  ->  -kT_bb  (bottom blackbody seed)
//   tau          ->  -tau    (total Thomson depth, set on the grid in main)
//   cosIncl      ->  -incidence (emergent intensity reported at GL nodes)
// ============================================================
static void save_emergent_compps(const RadField& rad, const RTGrids& g,
                                  const ModelParams& par, double Te_keV)
{
	// Save into the run directory, same logic as save_results() in the
	// production path: results/<hash>/.  Column layout matches the
	// production emergent_iter file so existing readers work.
	const char* dir = par.run_dir[0] ? par.run_dir : "results";
	mkdir("results", 0755);
	if (par.run_dir[0])
		mkdir(dir, 0755);

	char fname[256];
	snprintf(fname, sizeof(fname), "%s/emergent_compps.dat", dir);

	int fd = open(fname, O_WRONLY | O_CREAT | O_TRUNC, 0644);
	FILE* fp = fdopen(fd, "w");

	fprintf(fp, "# compps benchmark slab (isothermal pure scattering)\n");
	fprintf(fp, "# kTe=%.3f keV  kTbb=%.3f keV  tau_T=%.3f  geom=slab\n",
	        Te_keV, par.kT_bb, par.tau_slab);
	fprintf(fp, "# Col 1: E [eV]\n");
	fprintf(fp, "# Col 2: I_corona (incident, =0 here)\n");
	fprintf(fp, "# Col 3: I_disk   (bottom blackbody seed, incident)\n");
	for (int nm = 0; nm < g.NA; ++nm)
		fprintf(fp, "# Col %d: I_emergent(mu=%.6f)\n", nm + 4, g.mu[nm]);

	for (int ie = 0; ie < g.NE; ++ie)
	{
		fprintf(fp, "%.6e  %.6e  %.6e",
		        g.ene[ie], rad.illum.I_corona[ie], rad.illum.I_disk[ie]);
		for (int nm = 0; nm < g.NA; ++nm)
			fprintf(fp, "  %.6e", rad.Inu[0][nm][ie]);
		fprintf(fp, "\n");
	}

	fclose(fp);
	fprintf(stdout, "Saved: %s\n", fname);
}

static void save_emergent_tavg(const RadField& rad, const RTGrids& g,
                                  const ModelParams& par, double Te_keV)
{
	// Save into the run directory, same logic as save_results() in the
	// production path: results/<hash>/.  Column layout matches the
	// production emergent_iter file so existing readers work.
	const char* dir = par.run_dir[0] ? par.run_dir : "results";
	mkdir("results", 0755);
	if (par.run_dir[0])
		mkdir(dir, 0755);

	char fname[256];
	snprintf(fname, sizeof(fname), "%s/emergent_tavg.dat", dir);

	int fd = open(fname, O_WRONLY | O_CREAT | O_TRUNC, 0644);
	FILE* fp = fdopen(fd, "w");

	fprintf(fp, "# compps benchmark slab (isothermal pure scattering)\n");
	fprintf(fp, "# kTe=%.3f keV  kTbb=%.3f keV  tau_T=%.3f  geom=slab\n",
	        Te_keV, par.kT_bb, par.tau_slab);
	fprintf(fp, "# Col 1: E [eV]\n");
	fprintf(fp, "# Col 2: I_corona (incident, =0 here)\n");
	fprintf(fp, "# Col 3: I_disk   (bottom blackbody seed, incident)\n");
	for (int nm = 0; nm < g.NA; ++nm)
		fprintf(fp, "# Col %d: I_emergent(mu=%.6f)\n", nm + 4, g.mu[nm]);

	for (int ie = 0; ie < g.NE; ++ie)
	{
		fprintf(fp, "%.6e  %.6e  %.6e",
		        g.ene[ie], rad.illum.I_corona[ie], rad.illum.I_disk[ie]);
		for (int nm = 0; nm < g.NA; ++nm)
			fprintf(fp, "  %.6e", rad.Inu[0][nm][ie]);
		fprintf(fp, "\n");
	}

	fclose(fp);
	fprintf(stdout, "Saved: %s\n", fname);
}

template<class Cache>
static void run_test_rt_compps(RadField& rad, const RTGrids& g,
                               ModelParams& par, Cache& kcache)
{
	// kTe must be supplied (the slab electron temperature) and, for a
	// Maxwellian, compps itself requires kTe > 1 keV.
	if (par.kT_e <= 0.0)
	{
		fprintf(stderr,
		        "Error: -test_rt compps requires -kT_e (slab Te in keV)\n");
		exit(1);
	}

	// kTe [keV] -> uniform slab temperature [K]
	const double kB_eV = phys::k_B / phys::eV_to_erg;   // [eV/K]
	const double Te_K  = par.kT_e * 1.0e3 / kB_eV;

	fprintf(stdout, "\n=== RT TEST MODE: compps benchmark ===\n");
	fprintf(stdout, "  kTe=%.3f keV (T=%.3e K)  kTbb=%.3f keV  tau_T=%.3f\n",
	        par.kT_e, Te_K, par.kT_bb, par.tau_slab);

	// Isothermal pure-scattering slab at Te.
	fill_test_data(rad, g, par.nh, Te_K);

	// Illumination: no top source; isotropic blackbody seed at kTbb
	// enters from the bottom (compute_boundary_illumination uses
	// ill_top = I_corona, ill_bot = I_disk/2).
	blackbody(rad.illum.I_disk, g, par.kT_bb * 1.0e3);
	for (int ie = 0; ie < g.NE; ++ie)
		rad.illum.I_corona[ie] = 0.0;

	// Normalise the seed flux to Fx = xi*nH/(4pi)^2, exactly as
	// IllumSpec::compute does in the main program (xi = (4pi)^2 J / nH).
	double raw_disk = 0.0;
	for (int ie = 0; ie < g.NE - 1; ++ie)
	{
		double dE = g.ene[ie + 1] - g.ene[ie];
		if (g.ene[ie] > g.E_IN_LO && g.ene[ie] < g.E_IN_HI)
			raw_disk += 0.5 * (rad.illum.I_disk[ie] + rad.illum.I_disk[ie + 1]) * dE;
	}
	double xi = pow(10.0, par.zeta);
	double nH = pow(10.0, par.nh);
	double Fx = xi * nH / pow(phys::four_pi, 2);
	double scale_disk = Fx / raw_disk;
	for (int ie = 0; ie < g.NE; ++ie)
		rad.illum.I_disk[ie] *= scale_disk;

	fprintf(stdout, "  Seed normalised: Fx=%.4e  raw_disk=%.4e\n", Fx, raw_disk);

	compton_rt_solve(rad, g, par, kcache, par.maxiter);
	rad.compute_moments();
	save_emergent_compps(rad, g, par, par.kT_e);
}

// ============================================================
// test_avg mode
//
// Identical to the compps benchmark slab (isothermal, pure scattering,
// same normalisation to Fx = xi*nH/(4pi)^2), EXCEPT the illumination.
// Instead of an isotropic blackbody seed entering from the BOTTOM, the
// incident spectrum is the PRODUCTION corona shape (cutoffpl, nthcomp,
// ...) set via compute_corona_shape() exactly as in the main process,
// entering as a pencil beam from the TOP at the incidence node i_inc.
// There is no bottom source (I_disk = 0).
//
// In compute_boundary_illumination, I_corona enters as a pencil beam at
// the i_inc angle node (ill_top[ne] = 2*I_corona/wt[i_inc]); I_disk
// enters isotropically from the bottom (ill_bot = I_disk/2). So a top
// pencil beam means: seed -> I_corona, and I_disk = 0.
// ============================================================
template<class Cache>
static void run_test_rt_avg(RadField& rad, const RTGrids& g,
                            ModelParams& par, Cache& kcache)
{
	if (par.kT_e <= 0.0)
	{
		fprintf(stderr,
		        "Error: -test_rt test_avg requires -kT_e (slab Te in keV)\n");
		exit(1);
	}

	const double kB_eV = phys::k_B / phys::eV_to_erg;   // [eV/K]
	const double Te_K  = par.kT_e * 1.0e3 / kB_eV;

	fprintf(stdout, "\n=== RT TEST MODE: test_avg (top pencil beam) ===\n");
	fprintf(stdout, "  corona=%s  kTe=%.3f keV (T=%.3e K)  tau_T=%.3f\n",
	        par.corona, par.kT_e, Te_K, par.tau_slab);

	// Isothermal pure-scattering slab at Te (same as compps).
	fill_test_data(rad, g, par.nh, Te_K);

	// Illumination: production corona spectrum (set exactly as in the main
	// process) entering as a pencil beam from the TOP (-> I_corona at the
	// i_inc node); no bottom source (-> I_disk = 0).
	compute_corona_shape(rad.illum.I_corona, g, par);
	for (int ie = 0; ie < g.NE; ++ie)
		rad.illum.I_disk[ie] = 0.0;

	// Normalise the top corona flux to Fx = xi*nH/(4pi)^2 (matches the
	// frac<=0, corona-only branch of IllumSpec::compute).
	double raw_corona = 0.0;
	for (int ie = 0; ie < g.NE - 1; ++ie)
	{
		double dE = g.ene[ie + 1] - g.ene[ie];
		if (g.ene[ie] > g.E_IN_LO && g.ene[ie] < g.E_IN_HI)
			raw_corona += 0.5 * (rad.illum.I_corona[ie] + rad.illum.I_corona[ie + 1]) * dE;
	}
	double xi = pow(10.0, par.zeta);
	double nH = pow(10.0, par.nh);
	double Fx = xi * nH / pow(phys::four_pi, 2);
	double scale_corona = Fx / raw_corona;
	for (int ie = 0; ie < g.NE; ++ie)
		rad.illum.I_corona[ie] *= scale_corona;

	fprintf(stdout, "  Seed normalised: Fx=%.4e  raw_corona=%.4e\n",
	        Fx, raw_corona);

	compton_rt_solve(rad, g, par, kcache, par.maxiter);
	rad.compute_moments();
	save_emergent_tavg(rad, g, par, par.kT_e);
}

// ============================================================
// run_test_rt
// ============================================================
template<class Cache>
void run_test_rt(RadField& rad, const RTGrids& g, ModelParams& par,
                 Cache& kcache)
{	
	if (strcmp(par.test_mode, "compps") == 0)
	{
		run_test_rt_compps(rad, g, par, kcache);
		return;
	}
	if (strcmp(par.test_mode, "test_avg") == 0)
	{
		run_test_rt_avg(rad, g, par, kcache);
		return;
	}
	fprintf(stderr,
	        "Error: unknown test mode '%s' (valid: compps, test_avg)\n",
	        par.test_mode);
	exit(1);
}

// Explicit instantiations for both kernel-cache types.
template void run_test_rt<KernelCache>(
	RadField&, const RTGrids&, ModelParams&, KernelCache&);
template void run_test_rt<avgKernelCache>(
	RadField&, const RTGrids&, ModelParams&, avgKernelCache&);
