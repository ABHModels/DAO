#include "save_results.h"
#include "rt_grids.h"
#include "cddefines.h"
#include "cddrive.h"
#include "lines.h"
#include "rfield.h"
#include "dense.h"
#include "iso.h"
#include "taulines.h"
#include "transition.h"
#include "yield.h"
#include "heavy.h"
#include "ionbal.h"
#include "physconst.h"
#include <cstdio>
#include <cmath>
#include <string>
#include <sys/stat.h>

// Cloudy redefines fopen; undo it for file output
#undef fopen

void save_results(const RadField& rad, const RTGrids& g,
                  const ModelParams& par, int iter)
{
	const char* dir = par.run_dir[0] ? par.run_dir : "results";

	mkdir("results", 0755);
	if (par.run_dir[0])
		mkdir(dir, 0755);

	char fname[256];

	// 1. Emergent angle-dependent intensity at surface (nd=0)
	{
		snprintf(fname, sizeof(fname), "%s/emergent_iter%03d.dat", dir, iter);
		FILE* fp = fopen(fname, "w");
		fprintf(fp, "# Emergent specific intensity at surface (nd=0)  iter=%d\n", iter);
		fprintf(fp, "# Col 1: E [eV]\n");
		fprintf(fp, "# Col 2: I_corona (incident) [erg cm^-2 s^-1 eV^-1 sr^-1]\n");
		fprintf(fp, "# Col 3: I_disk (incident) [erg cm^-2 s^-1 eV^-1 sr^-1]\n");
		for (int nm = 0; nm < g.NA; ++nm)
			fprintf(fp, "# Col %d: I(mu=%.6f) [erg cm^-2 s^-1 eV^-1 sr^-1]\n",
			        nm + 4, g.mu[nm]);
		for (int ie = 0; ie < g.NE; ++ie)
		{
			fprintf(fp, "%.6e  %.6e  %.6e",
			        g.ene[ie], rad.illum.I_corona[ie], rad.illum.I_disk[ie]);
			for (int nm = 0; nm < g.NA; ++nm)
				fprintf(fp, "  %.6e", rad.Inu[0][nm][ie]);
			fprintf(fp, "\n");
		}
		fclose(fp);
		fprintf(stdout, "  Saved: %s\n", fname);
	}

	// 2. Angular moments at all depths
	{
		snprintf(fname, sizeof(fname), "%s/moments_iter%03d.dat", dir, iter);
		FILE* fp = fopen(fname, "w");
		fprintf(fp, "# Angular moments: J0, J2, J3 at all depths  iter=%d\n", iter);
		fprintf(fp, "# Col 1: E [eV]\n");
		fprintf(fp, "# Col 2: depth index\n");
		fprintf(fp, "# Col 3: tau_mid\n");
		fprintf(fp, "# Col 4: T_K [K]\n");
		fprintf(fp, "# Col 5: J0 [erg cm^-2 s^-1 eV^-1]\n");
		fprintf(fp, "# Col 6: J2\n");
		fprintf(fp, "# Col 7: J3\n");
		for (int id = 0; id < g.ND_MID; ++id)
		for (int ie = 0; ie < g.NE; ++ie)
		{
			fprintf(fp, "%.6e  %d  %.6e  %.6e  %.6e  %.6e  %.6e\n",
			        g.ene[ie], id, g.tau_mid[id], rad.T_K[id],
			        rad.J0[id][ie], rad.J2[id][ie], rad.J3[id][ie]);
		}
		fclose(fp);
		fprintf(stdout, "  Saved: %s\n", fname);
	}

	// 3. Temperature and opacity profile
	{
		snprintf(fname, sizeof(fname), "%s/profile_iter%03d.dat", dir, iter);
		FILE* fp = fopen(fname, "w");
		fprintf(fp, "# Depth profile  iter=%d\n", iter);
		fprintf(fp, "# Col 1: depth index\n");
		fprintf(fp, "# Col 2: tau_mid [Thomson]\n");
		fprintf(fp, "# Col 3: T_K [K]\n");
		fprintf(fp, "# Col 4: n_e [cm^-3]\n");
		fprintf(fp, "# Col 5: heating [erg cm^-3 s^-1]\n");
		fprintf(fp, "# Col 6: cooling [erg cm^-3 s^-1]\n");
		fprintf(fp, "# Col 7: log_xi = log10(4piJ) [erg cm^-2 s^-1] / nh\n");
		for (int id = 0; id < g.ND_MID; ++id)
		{
			fprintf(fp, "%d  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e\n",
			        id, g.tau_mid[id], rad.T_K[id], rad.n_e[id],
			        rad.heating[id], rad.cooling[id], rad.log_xi[id]);
		}
		fclose(fp);
		fprintf(stdout, "  Saved: %s\n", fname);
	}
}
