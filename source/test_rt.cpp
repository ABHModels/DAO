#include "test_rt.h"
#include "rt_grids.h"
#include "constants.h"
#include "compton_cross_section.h"
#include "compton_rt.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <fcntl.h>
#include <unistd.h>

// ============================================================
// Fill data for RT test mode: pure scattering atmosphere.
// ============================================================
static void fill_test_data(RadField& rad, const RTGrids& g,
                           double nh, double T_test)
{
	const double nH  = pow(10.0, nh);
	const double n_e = 1.21 * nH;

	// Compute scattering opacity once (same T at all depths)
	compute_compton_opacity(rad.ksct[0], g.NE, g.ene, T_test, n_e);

	for (int id = 0; id < g.ND_MID; ++id)
	{
		rad.T_K[id] = T_test;
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
		rad.cooling[id] = 0.0;
	}

	double kb_eV = phys::k_B / phys::eV_to_erg;
	fprintf(stdout, "Test mode: pure scattering atmosphere\n");
	fprintf(stdout, "  T=%.1e K  kT=%.2f eV  n_e=%.3e cm^-3\n",
		T_test, kb_eV * T_test, n_e);
}

// ============================================================
// Save emergent intensity
// ============================================================
static void save_emergent(const RadField& rad, const RTGrids& g,
                          const ModelParams& par)
{
	double kT_keV = 8.617333262e-5 * par.T_test / 1.0e3;
	char fname[128];
	snprintf(fname, sizeof(fname), "./data/emergent_kT%.0fkeV.dat", kT_keV);

	int fd = open(fname, O_WRONLY | O_CREAT | O_TRUNC, 0644);
	FILE* fp = fdopen(fd, "w");

	fprintf(fp, "# Emergent intensity at surface (nd=0, mu>0)\n");
	fprintf(fp, "# T=%.2e K  kT=%.2f keV  frac=%.2f\n",
	        par.T_test, kT_keV, par.frac);
	fprintf(fp, "# Col 1: E [eV]\n");
	fprintf(fp, "# Col 2: I_corona (incident)\n");
	fprintf(fp, "# Col 3: I_disk (incident)\n");
	int col = 4;
	for (int nm = 0; nm < g.NA; ++nm)
		if (g.mu[nm] > 0.0)
			fprintf(fp, "# Col %d: I(mu=%.4f)\n", col++, g.mu[nm]);

	for (int ie = 0; ie < g.NE; ++ie)
	{
		fprintf(fp, "%.6e  %.6e  %.6e",
		        g.ene[ie],
		        rad.illum.I_corona[ie],
		        rad.illum.I_disk[ie]);
		for (int nm = 0; nm < g.NA; ++nm)
			if (g.mu[nm] > 0.0)
				fprintf(fp, "  %.6e", rad.Inu[0][nm][ie]);
		fprintf(fp, "\n");
	}

	fclose(fp);
	fprintf(stdout, "Saved: %s\n", fname);
}

// ============================================================
// run_test_rt
// ============================================================
void run_test_rt(RadField& rad, const RTGrids& g, ModelParams& par,
                 KernelCache& kcache)
{
	fprintf(stdout, "\n=== RT TEST MODE ===\n");
	fill_test_data(rad, g, par.nh, par.T_test);
	compton_rt_solve(rad, g, par, kcache, par.maxiter);
	rad.compute_moments();
	save_emergent(rad, g, par);
}

