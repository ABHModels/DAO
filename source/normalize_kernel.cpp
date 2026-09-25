// ============================================================
// normalize_kernel — load kernel cache, normalize via A23, re-save
//
// For each (iT, ne_in): compute σ_raw from stored kernel,
// then scale all K(..., ne_in, ...) by σ_exact / σ_raw.
//
// Reads energy grid from results/emergent_iter002.dat.
// Overwrites the kernel cache file in place.
//
// Standalone — no Cloudy dependency.
// ============================================================

#include "rt_grids.h"
#include "compton_kernel.h"
#include "compton_cross_section.h"
#include "constants.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

static int read_energy_grid(const char* fname, std::vector<double>& ene_eV)
{
	FILE* fp = fopen(fname, "r");
	if (!fp) { fprintf(stderr, "Cannot open %s\n", fname); exit(1); }
	char line[1024];
	while (fgets(line, sizeof(line), fp))
	{
		if (line[0] == '#') continue;
		double E;
		if (sscanf(line, "%lf", &E) == 1)
			ene_eV.push_back(E);
	}
	fclose(fp);
	return (int)ene_eV.size();
}

int main(int argc, char* argv[])
{
	// ktype: 0 = approximate (_ap), 1 = exact (default)
	int ktype = 1;
	for (int i = 1; i < argc; ++i) {
		if (strcmp(argv[i], "-ktype") == 0 && i + 1 < argc)
			ktype = atoi(argv[++i]);
	}

	// --- Read energy grid ---
	std::vector<double> ene_vec;
	int NE = read_energy_grid("/Users/ym.huang/cloudy_test/results/8b877447/emergent_iter001.dat", ene_vec);
	double* ene_eV = ene_vec.data();
	fprintf(stdout, "Energy grid: NE=%d  E=[%.3f, %.3f] eV\n",
	        NE, ene_eV[0], ene_eV[NE - 1]);

	// --- Setup angle grid ---
	RTGrids g;
	g.init_angle();
	const int NA = g.NA;

	// --- Dimensionless energy ---
	const double mec2 = phys::m_e_c2_eV;
	double* x_grid = new double[NE];
	for (int i = 0; i < NE; ++i)
		x_grid[i] = ene_eV[i] / mec2;

	// nm_neg mapping
	int nm_neg[RTGrids::NA];
	for (int nm = 0; nm < NA; ++nm)
		nm_neg[nm] = NA - 1 - nm;

	// --- Load kernel ---
	fprintf(stdout, "Loading KernelCache...\n");
	KernelCache kcache;
	kcache.init(NE, ene_eV, NA, g.mu, g.wt, ktype);

	// --- Normalize ---
	fprintf(stdout, "Normalizing via A23 sum rule...\n");

	for (int iT = 0; iT < kcache.NT; ++iT)
	{
		fprintf(stdout, "  T=%.2e K (%d/%d)...",
		        kcache.T_grid[iT], iT + 1, kcache.NT);
		fflush(stdout);

		// Compute norm factor for each input energy
		double* norm = new double[NE];

		for (int ne = 0; ne < NE; ++ne)
		{
			double x = x_grid[ne];

			// A23 integral
			double integral = 0.0;
			for (int ne1 = 1; ne1 < NE; ++ne1)
			{
				double dx1 = x_grid[ne1] - x_grid[ne1 - 1];
				double f_a = 0.0, f_b = 0.0;

				for (int inm = 0; inm < NA; ++inm)
				{
					if (g.mu[inm] <= 0.0) continue;
					for (int inm1 = 0; inm1 < NA; ++inm1)
					{
						if (g.mu[inm1] <= 0.0) continue;
						double ww = g.wt[inm] * g.wt[inm1];
						double Ka = kcache.K(iT, ne1-1, inm1, ne, inm)
						          + kcache.K(iT, ne1-1, inm1, ne, nm_neg[inm]);
						double Kb = kcache.K(iT, ne1, inm1, ne, inm)
						          + kcache.K(iT, ne1, inm1, ne, nm_neg[inm]);
						f_a += ww * Ka;
						f_b += ww * Kb;
					}
				}
				f_a *= x_grid[ne1 - 1];
				f_b *= x_grid[ne1];
				integral += 0.5 * (f_a + f_b) * dx1;
			}

			double sigma_raw = integral / x;
			double sigma_exact = compton_cross_section(ene_eV[ne], kcache.T_grid[iT])
			                   / phys::sigma_T;

			norm[ne] = (sigma_raw > 1e-30) ? sigma_exact / sigma_raw : 1.0;
		}

		// Apply: scale K(iT, ne_out, nm_out, ne_in, nm_in) by norm[ne_in]
		for (int ne_out = 0; ne_out < NE; ++ne_out)
		for (int nm_out = 0; nm_out < NA; ++nm_out)
		{
			long r = kcache.row(iT, ne_out, nm_out);
			int lo = kcache.band_lo[r];
			int hi = kcache.band_hi[r];
			if (hi < lo) continue;

			double* dst = kcache.data + kcache.band_off[r];
			for (int ne_in = lo; ne_in <= hi; ++ne_in)
			{
				double nf = norm[ne_in];
				for (int nm_in = 0; nm_in < NA; ++nm_in)
					dst[long(ne_in - lo) * NA + nm_in] *= nf;
			}
		}

		delete[] norm;
		fprintf(stdout, " done.\n");
	}

	// --- Re-save ---
	const char* cache_dir = getenv("COMPTON_CACHE_DIR");
	if (!cache_dir) cache_dir = ".";

	const char* ksuffix = (ktype == 0) ? "_ap" : "";
	char cache_file[512];
	snprintf(cache_file, sizeof(cache_file),
	         "%s/kernel_norm_NE%d_NA%d_NT%d%s.bin", cache_dir, NE, NA, kcache.NT, ksuffix);
	fprintf(stdout, "Saving normalized kernel to %s...\n", cache_file);
	kcache.save(cache_file);

	kcache.free_memory();
	delete[] x_grid;

	fprintf(stdout, "Done.\n");
	return 0;
}
