// ============================================================
// Compare old kernel (CKERN02, full NA×NA) vs new kernel (CKERN03, symmetry-reduced)
//
// Loads both kernel files for NE=2768, NA=4, and compares
// K(iT, ne, nm, ne1, nm1) values for sampled points.
// ============================================================

#include "rt_grids.h"
#include "compton_kernel.h"
#include "constants.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

// ============================================================
// Old kernel loader (CKERN02 format from daov2)
// ============================================================
struct OldKernelCache
{
	int NT, NE, NA_K;
	double* T_grid;
	int*    band_lo;
	int*    band_hi;
	long*   band_off;
	double* data;
	long    data_size;

	OldKernelCache() : NT(0), NE(0), NA_K(0),
	    T_grid(nullptr), band_lo(nullptr), band_hi(nullptr),
	    band_off(nullptr), data(nullptr), data_size(0) {}

	inline long row(int iT, int ne, int nm) const
	{ return (long(iT) * NE + ne) * NA_K + nm; }

	inline double K(int iT, int ne, int nm, int ne1, int nm1) const
	{
		long r = row(iT, ne, nm);
		if (ne1 < band_lo[r] || ne1 > band_hi[r]) return 0.0;
		return data[band_off[r] + long(ne1 - band_lo[r]) * NA_K + nm1];
	}

	inline int lo(int iT, int ne, int nm) const { return band_lo[row(iT, ne, nm)]; }
	inline int hi(int iT, int ne, int nm) const { return band_hi[row(iT, ne, nm)]; }

	bool load(const char* filename)
	{
		FILE* fp = fopen(filename, "rb");
		if (!fp) return false;

		char magic[8];
		if (fread(magic, 1, 8, fp) != 8 || memcmp(magic, "CKERN02", 8) != 0)
		{ fclose(fp); return false; }

		fread(&NT, sizeof(int), 1, fp);
		fread(&NE, sizeof(int), 1, fp);
		fread(&NA_K, sizeof(int), 1, fp);
		fread(&data_size, sizeof(long), 1, fp);

		T_grid = new double[NT];
		fread(T_grid, sizeof(double), NT, fp);

		long n_rows = long(NT) * NE * NA_K;
		band_lo  = new int[n_rows];
		band_hi  = new int[n_rows];
		band_off = new long[n_rows];
		fread(band_lo,  sizeof(int),  n_rows, fp);
		fread(band_hi,  sizeof(int),  n_rows, fp);
		fread(band_off, sizeof(long), n_rows, fp);

		data = new double[data_size];
		size_t nread = fread(data, sizeof(double), data_size, fp);
		fclose(fp);

		fprintf(stdout, "  OldKernel: loaded %s  NT=%d NE=%d NA=%d  data=%.1f MB\n",
		        filename, NT, NE, NA_K,
		        data_size * 8.0 / (1024.0 * 1024.0));
		return (long)nread == data_size;
	}

	void free_memory()
	{
		delete[] T_grid;   delete[] band_lo;
		delete[] band_hi;  delete[] band_off;
		delete[] data;
	}
};

// ============================================================
// Read energy grid
// ============================================================
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

int main()
{
	const char* old_file = "kernel/kernel_norm_NE2768_NA4_NT50.bin";
	const char* new_file = "kernel/kernel_norm_NE2768_NI6_NT50.bin";

	// --- Load old kernel ---
	fprintf(stdout, "Loading old kernel...\n");
	OldKernelCache old_k;
	if (!old_k.load(old_file))
	{ fprintf(stderr, "Failed to load %s\n", old_file); return 1; }

	// --- Load new kernel ---
	fprintf(stdout, "Loading new kernel...\n");

	// Build angle grid (NA=4)
	// We need to match the old kernel's NA=4 angle grid
	// The old kernel used RTGrids with NA=4 at the time
	// We'll use 4-point GL quadrature manually
	const int NA = old_k.NA_K;
	fprintf(stdout, "  NA = %d (from old kernel)\n", NA);

	// Compute 4-pt GL nodes on [-1,1]
	double mu[4], wt[4];
	{
		const int n = 4;
		for (int i = 0; i < n; ++i)
		{
			double x = -cos(M_PI * (i + 0.75) / (n + 0.5));
			for (int iter = 0; iter < 100; ++iter)
			{
				double p0 = 1.0, p1 = x;
				for (int j = 2; j <= n; ++j)
				{
					double p2 = ((2*j - 1) * x * p1 - (j - 1) * p0) / j;
					p0 = p1; p1 = p2;
				}
				double dp = n * (x * p1 - p0) / (x * x - 1.0);
				double dx = -p1 / dp;
				x += dx;
				if (fabs(dx) < 1e-15) break;
			}
			mu[i] = x;
			double p0 = 1.0, p1 = x;
			for (int j = 2; j <= n; ++j)
			{
				double p2 = ((2*j - 1) * x * p1 - (j - 1) * p0) / j;
				p0 = p1; p1 = p2;
			}
			double dp = n * (x * p1 - p0) / (x * x - 1.0);
			wt[i] = 2.0 / ((1.0 - x * x) * dp * dp);
		}
	}

	fprintf(stdout, "  GL angles:\n");
	for (int i = 0; i < NA; ++i)
		fprintf(stdout, "    mu[%d] = %+.15f   wt[%d] = %.15f\n", i, mu[i], i, wt[i]);

	// We need an energy grid to init the new kernel
	// Extract NE and T_grid from the old kernel to build compatible new kernel
	int NE = old_k.NE;
	int NT = old_k.NT;

	// Build a dummy energy grid (we just need NE to match; actual values don't matter
	// for loading since we're loading from file)
	// But init() needs ene_eV for x_grid computation if cache miss.
	// Since we're loading from file, we just need NE to match.
	// Use a simple log-spaced grid as placeholder (won't be used if file loads)
	double* ene_eV = new double[NE];
	for (int i = 0; i < NE; ++i)
		ene_eV[i] = 1.0 * pow(1e6, double(i) / (NE - 1));

	KernelCache new_k;
	new_k.init(NE, ene_eV, NA, mu, wt, 1);

	// Verify dimensions match
	fprintf(stdout, "\n  Old: NT=%d NE=%d NA_K=%d\n", old_k.NT, old_k.NE, old_k.NA_K);
	fprintf(stdout, "  New: NT=%d NE=%d NA_full=%d n_indep=%d\n",
	        new_k.NT, new_k.NE, new_k.NA_full, new_k.n_indep);

	if (old_k.NT != new_k.NT || old_k.NE != new_k.NE || old_k.NA_K != new_k.NA_full)
	{
		fprintf(stderr, "Dimension mismatch!\n");
		return 1;
	}

	// --- Compare K values ---
	fprintf(stdout, "\n=== Comparing kernel values ===\n");

	int iT_samples[] = {0, NT/4, NT/2, 3*NT/4, NT-1};
	int ne_samples[] = {NE/8, NE/4, NE/2, 3*NE/4, 7*NE/8};

	double max_abs_err = 0.0, max_rel_err = 0.0;
	long n_compared = 0, n_nonzero = 0, n_mismatch = 0;
	long n_old_only = 0, n_new_only = 0;

	for (int iTs = 0; iTs < 5; ++iTs)
	for (int nes = 0; nes < 5; ++nes)
	{
		int iT = iT_samples[iTs];
		int ne = ne_samples[nes];

		// Find band limits from both kernels
		// Old: per (iT, ne, nm)
		// New: global per (iT, ne)
		int ne1_lo_new = new_k.lo(iT, ne);
		int ne1_hi_new = new_k.hi(iT, ne);

		// Union of old bands across all nm
		int ne1_lo_old = NE, ne1_hi_old = -1;
		for (int nm = 0; nm < NA; ++nm)
		{
			int lo = old_k.lo(iT, ne, nm);
			int hi = old_k.hi(iT, ne, nm);
			if (lo < ne1_lo_old) ne1_lo_old = lo;
			if (hi > ne1_hi_old) ne1_hi_old = hi;
		}

		int ne1_lo = (ne1_lo_old < ne1_lo_new) ? ne1_lo_old : ne1_lo_new;
		int ne1_hi = (ne1_hi_old > ne1_hi_new) ? ne1_hi_old : ne1_hi_new;
		if (ne1_hi < ne1_lo) continue;

		// Sample ne1 points within band
		int ne1_list[] = {ne1_lo, (ne1_lo + ne1_hi)/2, ne1_hi};

		for (int ne1s = 0; ne1s < 3; ++ne1s)
		{
			int ne1 = ne1_list[ne1s];
			if (ne1 < 0 || ne1 >= NE) continue;

			for (int nm = 0; nm < NA; ++nm)
			for (int nm1 = 0; nm1 < NA; ++nm1)
			{
				double v_old = old_k.K(iT, ne, nm, ne1, nm1);
				double v_new = new_k.K(iT, ne, nm, ne1, nm1);

				double abs_err = fabs(v_old - v_new);
				double denom = fabs(v_old) + fabs(v_new);
				double rel_err = (denom > 1e-30) ? 2.0 * abs_err / denom : 0.0;

				if (abs_err > max_abs_err) max_abs_err = abs_err;
				if (rel_err > max_rel_err) max_rel_err = rel_err;
				++n_compared;

				if (v_old != 0.0 || v_new != 0.0) ++n_nonzero;
				if (v_old != 0.0 && v_new == 0.0) ++n_old_only;
				if (v_old == 0.0 && v_new != 0.0) ++n_new_only;

				if (rel_err > 1e-6)
				{
					++n_mismatch;
					if (n_mismatch <= 20)
					{
						fprintf(stdout,
							"  DIFF iT=%d ne=%d nm=%d ne1=%d nm1=%d: "
							"old=%.6e  new=%.6e  rel=%.2e\n",
							iT, ne, nm, ne1, nm1, v_old, v_new, rel_err);
					}
				}
			}
		}
	}

	fprintf(stdout, "\n=== Summary ===\n");
	fprintf(stdout, "  Compared:   %ld values\n", n_compared);
	fprintf(stdout, "  Non-zero:   %ld\n", n_nonzero);
	fprintf(stdout, "  Old-only:   %ld (old non-zero, new zero — trimmed by narrower band)\n", n_old_only);
	fprintf(stdout, "  New-only:   %ld (new non-zero, old zero)\n", n_new_only);
	fprintf(stdout, "  Mismatch:   %ld (rel_err > 1e-6)\n", n_mismatch);
	fprintf(stdout, "  Max abs err: %.6e\n", max_abs_err);
	fprintf(stdout, "  Max rel err: %.6e\n", max_rel_err);

	// --- Also print a full 4x4 table for one sample point ---
	fprintf(stdout, "\n=== Full table: iT=%d  ne=%d  ne1=mid-band ===\n",
	        NT/2, NE/2);
	{
		int iT = NT/2, ne = NE/2;
		int ne1_lo = new_k.lo(iT, ne);
		int ne1_hi = new_k.hi(iT, ne);
		int ne1 = (ne1_lo + ne1_hi) / 2;

		fprintf(stdout, "  ne1=%d\n\n", ne1);
		fprintf(stdout, "  (nm,nm1)  old K            new K            ratio\n");
		fprintf(stdout, "  ------    ---------------  ---------------  -----\n");

		for (int nm = 0; nm < NA; ++nm)
		for (int nm1 = 0; nm1 < NA; ++nm1)
		{
			double v_old = old_k.K(iT, ne, nm, ne1, nm1);
			double v_new = new_k.K(iT, ne, nm, ne1, nm1);
			double ratio = (fabs(v_old) > 1e-30) ? v_new / v_old : 0.0;
			fprintf(stdout, "  (%d,%d)    %.9e  %.9e  %.6f\n",
			        nm, nm1, v_old, v_new, ratio);
		}
	}

	old_k.free_memory();
	new_k.free_memory();
	delete[] ene_eV;
	return 0;
}
