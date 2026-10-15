// ============================================================
// multiscat — iterate the azimuth-integrated Compton kernel R
// to evolve a monochromatic, mono-directional photon distribution
// over n scatterings:
//
//   I^n(x, mu) = x^2 ∫∫ R(x, mu; x_1, mu_1) I^{n-1}(x_1, mu_1) / x_1^2  dmu_1 dx_1
//
// I^0 is a delta at (E_0, mu_0); the closest grid points are used.
// Snapshots at n = 1, 10, 50, 100, 300 are written as 1D spectra
//   I^n(E) = Σ_mu I^n(E, mu) × wt[mu]
//
// Output: data/multiscat_E0_<E0>keV_T<T>.dat
// ============================================================

#include "rt_grids.h"
#include "compton_kernel.h"
#include "constants.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

static int read_energy_grid(const char* fname,
                            std::vector<double>& ene_eV,
                            std::vector<double>& wid_eV)
{
	FILE* fp = fopen(fname, "r");
	if (!fp) { fprintf(stderr, "Cannot open %s\n", fname); exit(1); }
	char line[1024];
	while (fgets(line, sizeof(line), fp))
	{
		if (line[0] == '#' || line[0] == '\n') continue;
		double E, w;
		if (sscanf(line, "%lf %lf", &E, &w) == 2)
		{
			ene_eV.push_back(E);
			wid_eV.push_back(w);
		}
	}
	fclose(fp);
	return (int)ene_eV.size();
}

static int find_nearest_log(const double* arr, int n, double val)
{
	double lv = log10(val);
	int best = 0;
	double bd = fabs(log10(arr[0]) - lv);
	for (int i = 1; i < n; ++i)
	{
		double d = fabs(log10(arr[i]) - lv);
		if (d < bd) { bd = d; best = i; }
	}
	return best;
}

static int find_nearest(const double* arr, int n, double val)
{
	int best = 0;
	double bd = fabs(arr[0] - val);
	for (int i = 1; i < n; ++i)
	{
		double d = fabs(arr[i] - val);
		if (d < bd) { bd = d; best = i; }
	}
	return best;
}

int main(int argc, char* argv[])
{
	const char* kernel_file = "kernel/kernel_norm_NE3403_NI20_NT50.bin";
	const char* egrid_file  = "cloudy_energy_grid.dat";
	const char* out_path    = nullptr;
	double E0_eV = 6.4e3;
	double mu0   = cos(M_PI / 4.0);   // 0.707
	double T     = 1e8;
	int    N_max = 300;

	for (int i = 1; i < argc; ++i)
	{
		if      (!strcmp(argv[i], "-kernel") && i+1 < argc) kernel_file = argv[++i];
		else if (!strcmp(argv[i], "-egrid")  && i+1 < argc) egrid_file  = argv[++i];
		else if (!strcmp(argv[i], "-out")    && i+1 < argc) out_path    = argv[++i];
		else if (!strcmp(argv[i], "-E0")     && i+1 < argc) E0_eV       = atof(argv[++i]);
		else if (!strcmp(argv[i], "-mu0")    && i+1 < argc) mu0         = atof(argv[++i]);
		else if (!strcmp(argv[i], "-T")      && i+1 < argc) T           = atof(argv[++i]);
		else if (!strcmp(argv[i], "-N")      && i+1 < argc) N_max       = atoi(argv[++i]);
		else { fprintf(stderr, "unknown arg: %s\n", argv[i]); return 1; }
	}

	// --- Read energy grid (Cloudy 3403 bins) ---
	std::vector<double> ene_vec, wid_vec;
	int NE = read_energy_grid(egrid_file, ene_vec, wid_vec);
	const double* ene_eV = ene_vec.data();
	const double* wid_eV = wid_vec.data();
	fprintf(stdout, "Energy grid: NE=%d  E=[%.3f, %.3f] eV\n",
	        NE, ene_eV[0], ene_eV[NE-1]);
	fflush(stdout);

	// --- Angle grid (8-pt GL on [-1,1]) ---
	RTGrids g;
	g.init_angle();
	const int NA = g.NA;

	// --- KernelCache shell + load ---
	KernelCache kc;
	kc.NT      = N_T_CACHE;
	kc.NE      = NE;
	kc.NA_full = NA;

	kc.T_grid = new double[kc.NT];
	double log_lo = log10(T_CACHE_LO);
	double log_hi = log10(T_CACHE_HI);
	double dlog   = (kc.NT > 1) ? (log_hi - log_lo) / (kc.NT - 1) : 0.0;
	for (int i = 0; i < kc.NT; ++i)
		kc.T_grid[i] = pow(10.0, log_lo + i * dlog);

	int p = NA / 2;
	int max_indep = p * (p + 1);
	kc.canon     = new int[NA * NA];
	kc.indep_nm  = new int[max_indep];
	kc.indep_nm1 = new int[max_indep];
	for (int i = 0; i < NA*NA; ++i) kc.canon[i] = -1;
	kc.n_indep = 0;
	for (int nm = 0; nm < NA; ++nm)
	for (int nm1 = 0; nm1 < NA; ++nm1)
	{
		if (kc.canon[nm * NA + nm1] >= 0) continue;
		int ia = kc.n_indep++;
		kc.indep_nm[ia]  = nm;
		kc.indep_nm1[ia] = nm1;
		int nm_r  = NA - 1 - nm;
		int nm1_r = NA - 1 - nm1;
		kc.canon[nm   * NA + nm1 ] = ia;
		kc.canon[nm1  * NA + nm  ] = ia;
		kc.canon[nm_r * NA + nm1_r] = ia;
		kc.canon[nm1_r* NA + nm_r ] = ia;
	}
	if (!kc.load(kernel_file))
	{
		fprintf(stderr, "ERROR: failed to load %s\n", kernel_file);
		return 1;
	}

	// --- Snap initial condition to nearest grid points ---
	const double mec2 = 511.0e3;
	int iT    = kc.find_T(T);
	int ne0   = find_nearest_log(ene_eV, NE, E0_eV);
	int nm0   = find_nearest(g.mu, NA, mu0);
	double T_act    = kc.T_grid[iT];
	double E0_act   = ene_eV[ne0];
	double mu0_act  = g.mu[nm0];

	fprintf(stdout, "Initial photon: E0=%.3e eV (ne0=%d), mu0=%.4f (nm0=%d)\n",
	        E0_act, ne0, mu0_act, nm0);
	fprintf(stdout, "T = %.3e K (iT=%d), n_max = %d\n", T_act, iT, N_max);

	// --- Pre-compute dimensionless x and dx ---
	std::vector<double> xg(NE), dxg(NE);
	for (int j = 0; j < NE; ++j) {
		xg[j]  = ene_eV[j] / mec2;
		dxg[j] = wid_eV[j] / mec2;
	}

	// --- Allocate I[ne][nm] and I_new[ne][nm] ---
	std::vector<double> I_curr(long(NE) * NA, 0.0);
	std::vector<double> I_next(long(NE) * NA, 0.0);

	// I^0: delta at (ne0, nm0), normalized so that ∫∫ I^0 dx dmu = 1
	I_curr[long(ne0) * NA + nm0] = 1.0 / (dxg[ne0] * g.wt[nm0]);

	// Snapshots
	std::vector<int> snap_ns;
	snap_ns.push_back(1);
	snap_ns.push_back(10);
	snap_ns.push_back(50);
	snap_ns.push_back(100);
	snap_ns.push_back(N_max);

	std::vector<std::vector<double>> snap_spec;          // I^n(E) [µ-summed]
	std::vector<std::vector<double>> snap_full;          // I^n[ne*NA + nm]
	std::vector<int>                  snap_taken;

	// Save n=0 (delta) too
	{
		std::vector<double> spec(NE, 0.0);
		for (int j = 0; j < NE; ++j)
			for (int k = 0; k < NA; ++k)
				spec[j] += I_curr[long(j)*NA + k] * g.wt[k];
		snap_spec.push_back(spec);
		snap_full.push_back(I_curr);
		snap_taken.push_back(0);
	}

	// --- Iterate ---
	// Track the active range of ne where I > 0 (grows with n).
	int act_lo = ne0, act_hi = ne0;

	for (int n = 1; n <= N_max; ++n)
	{
		for (int j = 0; j < long(NE) * NA; ++j) I_next[j] = 0.0;

		// Output is reachable only from input bins in [act_lo, act_hi];
		// outgoing range expands by the kernel band on each side.
		int out_lo = NE, out_hi = -1;
		for (int ne_i = act_lo; ne_i <= act_hi; ++ne_i)
		{
			int lo_i = kc.lo(iT, ne_i);
			int hi_i = kc.hi(iT, ne_i);
			if (hi_i < lo_i) continue;
			if (lo_i < out_lo) out_lo = lo_i;
			if (hi_i > out_hi) out_hi = hi_i;
		}
		if (out_hi < out_lo) { out_lo = act_lo; out_hi = act_hi; }

		for (int ne_o = out_lo; ne_o <= out_hi; ++ne_o)
		{
			int lo = kc.lo(iT, ne_o);
			int hi = kc.hi(iT, ne_o);
			if (hi < lo) continue;
			// Restrict the incoming sum to the intersection with the active set
			int lo_use = (lo > act_lo) ? lo : act_lo;
			int hi_use = (hi < act_hi) ? hi : act_hi;
			if (hi_use < lo_use) continue;

			double xo2 = xg[ne_o] * xg[ne_o];

			for (int nm_o = 0; nm_o < NA; ++nm_o)
			{
				double sum = 0.0;
				for (int ne_i = lo_use; ne_i <= hi_use; ++ne_i)
				{
					double inv_xi2 = 1.0 / (xg[ne_i] * xg[ne_i]);
					double dxi = dxg[ne_i];
					for (int nm_i = 0; nm_i < NA; ++nm_i)
					{
						double Iv = I_curr[long(ne_i) * NA + nm_i];
						if (Iv == 0.0) continue;
						double Kv = kc.K(iT, ne_o, nm_o, ne_i, nm_i);
						sum += Kv * Iv * inv_xi2 * g.wt[nm_i] * dxi;
					}
				}
				I_next[long(ne_o) * NA + nm_o] = xo2 * sum;
			}
		}
		std::swap(I_curr, I_next);
		act_lo = out_lo;
		act_hi = out_hi;

		// Progress + diagnostics
		if (n == 1 || n % 10 == 0) {
			double total = 0.0, mean_E = 0.0;
			for (int j = act_lo; j <= act_hi; ++j)
				for (int k = 0; k < NA; ++k) {
					double w = I_curr[long(j)*NA + k] * dxg[j] * g.wt[k];
					total  += w;
					mean_E += w * ene_eV[j];
				}
			fprintf(stdout, "  n=%4d  active=[%d,%d]  total=%.4f  <E>=%.3f keV\n",
			        n, act_lo, act_hi, total, mean_E/total/1e3);
			fflush(stdout);
		}

		// Save snapshot
		bool snap = false;
		for (size_t s = 0; s < snap_ns.size(); ++s)
			if (snap_ns[s] == n) { snap = true; break; }
		if (snap) {
			std::vector<double> spec(NE, 0.0);
			for (int j = 0; j < NE; ++j)
				for (int k = 0; k < NA; ++k)
					spec[j] += I_curr[long(j)*NA + k] * g.wt[k];
			snap_spec.push_back(spec);
			snap_full.push_back(I_curr);
			snap_taken.push_back(n);
		}
	}

	// --- Write output: one column per snapshot ---
	char fname[1024];
	if (out_path) snprintf(fname, sizeof(fname), "%s", out_path);
	else snprintf(fname, sizeof(fname),
	              "data/multiscat_E0_%.1fkeV_T%.0e.dat", E0_act/1e3, T_act);

	FILE* fp = fopen(fname, "w");
	if (!fp) { fprintf(stderr, "cannot write %s\n", fname); return 1; }

	fprintf(fp, "# Multi-scatter spectra I^n(E) = Σ_mu I^n(E,mu) wt[mu]\n");
	fprintf(fp, "# E0   = %.6e eV   (ne0=%d)\n", E0_act, ne0);
	fprintf(fp, "# mu0  = %.6f      (nm0=%d, requested 0.707)\n", mu0_act, nm0);
	fprintf(fp, "# T    = %.6e K    (iT=%d, kT_e = %.3f keV)\n",
	        T_act, iT, 8.617333262e-5*T_act/1e3);
	fprintf(fp, "# n_snapshots =");
	for (size_t s = 0; s < snap_taken.size(); ++s)
		fprintf(fp, " %d", snap_taken[s]);
	fprintf(fp, "\n");
	fprintf(fp, "# mu_grid     =");
	for (int k = 0; k < NA; ++k) fprintf(fp, " %+.8f", g.mu[k]);
	fprintf(fp, "\n");
	fprintf(fp, "# wt_grid     =");
	for (int k = 0; k < NA; ++k) fprintf(fp, " %.8f", g.wt[k]);
	fprintf(fp, "\n");
	fprintf(fp, "# === SPECTRA: col1 = E [eV], cols2..= I^n(E) [mu-summed]\n");
	for (int j = 0; j < NE; ++j)
	{
		fprintf(fp, "%.6e", ene_eV[j]);
		for (size_t k = 0; k < snap_spec.size(); ++k)
			fprintf(fp, " %.6e", snap_spec[k][j]);
		fprintf(fp, "\n");
	}

	// 2D blocks: one block per snapshot, NE rows × NA cols of I^n[ne, nm]
	for (size_t k = 0; k < snap_full.size(); ++k)
	{
		fprintf(fp, "# === FULL2D n=%d: rows=E (NE=%d), cols=mu (NA=%d)\n",
		        snap_taken[k], NE, NA);
		const std::vector<double>& Ifull = snap_full[k];
		for (int j = 0; j < NE; ++j)
		{
			for (int kk = 0; kk < NA; ++kk)
				fprintf(fp, "%.6e ", Ifull[long(j)*NA + kk]);
			fprintf(fp, "\n");
		}
	}
	fclose(fp);
	fprintf(stdout, "wrote %s\n", fname);

	kc.free_memory();
	return 0;
}
