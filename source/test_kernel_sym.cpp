// ============================================================
// test_kernel_sym — verify symmetry-reduced kernel against direct computation
//
// For selected (iT, ne, ne1), compute all 36 (nm,nm1) kernel
// elements both ways:
//   (A) kcache.K()              — from symmetry-reduced cache (normalized)
//   (B) compton_kernel_element() — direct computation (raw)
//
// Since the cache is normalized, (A) = (B) × norm_factor.
// We verify:
//   1. The norm factor is the same for all 36 pairs at a given (iT,ne,ne1)
//   2. Symmetry-related pairs give identical kcache.K() values
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
	int ktype = 1;
	const char* egrid_file = nullptr;
	for (int i = 1; i < argc; ++i) {
		if (strcmp(argv[i], "-ktype") == 0 && i + 1 < argc)
			ktype = atoi(argv[++i]);
		else if (strcmp(argv[i], "-egrid") == 0 && i + 1 < argc)
			egrid_file = argv[++i];
	}
	if (!egrid_file) {
		fprintf(stderr, "Usage: %s -egrid <energy_grid_file> [-ktype 0|1]\n", argv[0]);
		return 1;
	}

	// --- Read energy grid ---
	std::vector<double> ene_vec;
	int NE = read_energy_grid(egrid_file, ene_vec);
	double* ene_eV = ene_vec.data();

	// --- Setup angle grid ---
	RTGrids g;
	g.init_angle();
	const int NA = g.NA;

	// --- Load kernel cache ---
	KernelCache kcache;
	kcache.init(NE, ene_eV, NA, g.mu, g.wt, ktype);

	// --- Dimensionless energy ---
	const double mec2_eV = 511.0e3;  // must match init()'s value
	double* x_grid = new double[NE];
	for (int i = 0; i < NE; ++i)
		x_grid[i] = ene_eV[i] / mec2_eV;

	// --- Sample points ---
	int iT_samples[] = {0, kcache.NT / 2, kcache.NT - 1};
	int ne_samples[] = {NE / 4, NE / 2, 3 * NE / 4};
	const int N_iT = 3, N_ne = 3;

	int total_pass = 0, total_fail = 0;

	for (int iTs = 0; iTs < N_iT; ++iTs)
	for (int nes = 0; nes < N_ne; ++nes)
	{
		int iT = iT_samples[iTs];
		int ne = ne_samples[nes];
		double T = kcache.T_grid[iT];

		// Pick ne1 in the middle of the band where kernel is large
		int ne1_lo = kcache.lo(iT, ne);
		int ne1_hi = kcache.hi(iT, ne);
		if (ne1_hi < ne1_lo) continue;
		int ne1 = (ne1_lo + ne1_hi) / 2;

		fprintf(stdout,
			"\n======================================================================\n"
			"  iT=%d  T=%.2e K  |  ne=%d  E=%.1f eV  |  ne1=%d  E1=%.1f eV\n"
			"======================================================================\n",
			iT, T, ne, ene_eV[ne], ne1, ene_eV[ne1]);

		// Compute normalization factor for this (iT, ne1):
		// norm = sigma_exact / sigma_raw, applied uniformly to all angle pairs
		// We extract it from the first non-zero pair.
		double cache_val[6][6], direct_val[6][6];
		for (int nm = 0; nm < NA; ++nm)
		for (int nm1 = 0; nm1 < NA; ++nm1)
		{
			cache_val[nm][nm1]  = kcache.K(iT, ne, nm, ne1, nm1);
			direct_val[nm][nm1] = compton_kernel_element(
				x_grid[ne], g.mu[nm], x_grid[ne1], g.mu[nm1], T, ktype);
		}

		// Find norm factor from first non-zero pair
		double norm_factor = 1.0;
		for (int nm = 0; nm < NA && norm_factor == 1.0; ++nm)
		for (int nm1 = 0; nm1 < NA && norm_factor == 1.0; ++nm1)
			if (direct_val[nm][nm1] > 1e-30 && cache_val[nm][nm1] > 1e-30)
				norm_factor = cache_val[nm][nm1] / direct_val[nm][nm1];

		fprintf(stdout, "  A23 norm factor for ne1=%d: %.15f\n\n", ne1, norm_factor);

		fprintf(stdout,
			"  (nm,nm1)  ia  |  cache K          norm×direct K    |  ratio (should be 1)  \n"
			"  -------  ---  |  ---------------  ---------------  |  ---------------------\n");

		double ratio[6][6];
		for (int nm = 0; nm < NA; ++nm)
		for (int nm1 = 0; nm1 < NA; ++nm1)
		{
			double norm_direct = direct_val[nm][nm1] * norm_factor;
			ratio[nm][nm1] = (norm_direct > 1e-30)
			               ? cache_val[nm][nm1] / norm_direct : 0.0;
			int ia = kcache.canon[nm * NA + nm1];
			fprintf(stdout, "  (%d,%d)    %2d  |  %.9e  %.9e  |  %.15f\n",
			        nm, nm1, ia,
			        cache_val[nm][nm1], norm_direct, ratio[nm][nm1]);
		}

		// --- Check 1: all symmetry-related pairs give same cache value ---
		fprintf(stdout, "\n  Check 1: symmetry-related pairs give same cache value\n");
		int sym_ok = 0, sym_fail = 0;
		for (int nm = 0; nm < NA; ++nm)
		for (int nm1 = 0; nm1 < NA; ++nm1)
		{
			// Compare with transpose
			double diff_t = fabs(cache_val[nm][nm1] - cache_val[nm1][nm]);
			// Compare with reflection
			int nm_r  = NA - 1 - nm;
			int nm1_r = NA - 1 - nm1;
			double diff_r = fabs(cache_val[nm][nm1] - cache_val[nm_r][nm1_r]);

			double denom = fabs(cache_val[nm][nm1]) + 1e-30;
			if (diff_t / denom > 1e-14 || diff_r / denom > 1e-14) {
				fprintf(stdout, "    FAIL (%d,%d): trans_diff=%.2e  refl_diff=%.2e\n",
				        nm, nm1, diff_t / denom, diff_r / denom);
				++sym_fail;
			} else {
				++sym_ok;
			}
		}
		fprintf(stdout, "    Result: %d PASS, %d FAIL\n", sym_ok, sym_fail);

		// --- Check 2: cache == norm × direct  (ratio should be 1.0) ---
		fprintf(stdout, "\n  Check 2: cache == norm × direct (ratio should be 1.0)\n");

		double max_ratio_err = 0.0;
		int ratio_ok = 0, ratio_fail = 0;
		for (int nm = 0; nm < NA; ++nm)
		for (int nm1 = 0; nm1 < NA; ++nm1)
		{
			if (ratio[nm][nm1] == 0.0) { ++ratio_ok; continue; }
			double err = fabs(ratio[nm][nm1] - 1.0);
			if (err > max_ratio_err) max_ratio_err = err;
			if (err > 1e-10)
			{
				fprintf(stdout, "    FAIL (%d,%d): ratio=%.15f  err=%.2e\n",
				        nm, nm1, ratio[nm][nm1], err);
				++ratio_fail;
			} else {
				++ratio_ok;
			}
		}
		fprintf(stdout, "    Max deviation from 1.0: %.2e\n", max_ratio_err);
		fprintf(stdout, "    Result: %d PASS, %d FAIL\n", ratio_ok, ratio_fail);

		total_pass += (sym_fail == 0 && ratio_fail == 0) ? 1 : 0;
		total_fail += (sym_fail > 0 || ratio_fail > 0)   ? 1 : 0;
	}

	fprintf(stdout,
		"\n======================================================================\n"
		"  OVERALL: %d sample points PASS, %d FAIL\n"
		"======================================================================\n",
		total_pass, total_fail);

	kcache.free_memory();
	delete[] x_grid;
	return total_fail > 0 ? 1 : 0;
}
