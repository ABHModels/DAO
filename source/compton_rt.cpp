#include "compton_rt.h"
#include "source.h"
#include "constants.h"
#include <cmath>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <chrono>

// ============================================================
// Flat-array indexing helpers
// ============================================================
static int g_ND, g_NM, g_NE;

static inline long idx3(int nd, int nm, int ne)
{
	return (long(nd) * g_NM + nm) * g_NE + ne;
}

static inline long idx2(int nd, int ne)
{
	return long(nd) * g_NE + ne;
}

// ============================================================
// Fritsch-Butland monotonic derivative
//   de la Cruz Rodriguez & Piskunov (2013), Eqs. 11-14
//
//   Sm1, S0, Sp1 : source function at three consecutive points
//   hm, hp       : optical depth intervals (backward, forward)
//   Returns dS/dtau at the middle point S0
// ============================================================
static inline double fritsch_butland(double Sm1, double S0, double Sp1,
                                     double hm, double hp)
{
	double dm = (S0 - Sm1) / hm;
	double dp = (Sp1 - S0) / hp;

	if (dm * dp <= 0.0) return 0.0;

	double a = (1.0 / 3.0) * (1.0 + hp / (hp + hm));
	return dm * dp / (a * dp + (1.0 - a) * dm);
}

// ============================================================
//  CUBIC BEZIER SC
//    de la Cruz Rodriguez & Piskunov (2013), Eq. 20
//
//    I_D = I_U exp(-d) + ahat*S_D + bhat*S_U + ghat*E + ehat*F
//
//    Scalar control points (non-DELO):
//      E = S_D + (d/3) * S'_D   (near downwind)
//      F = S_U - (d/3) * S'_U   (near upwind)
// ============================================================

static inline void bezier3_coeffs(double d,
                                  double& ahat, double& bhat,
                                  double& ghat, double& ehat)
{
	if (d < 1e-3)
	{
		// Taylor expansion (Appendix B)
		double d2 = d * d, d3 = d2 * d, d4 = d3 * d;
		ahat = d / 4.0 - d2 / 20.0  + d3 / 120.0 - d4 / 840.0;
		bhat = d / 4.0 - d2 / 5.0   + d3 / 12.0  - d4 / 42.0;
		ghat = d / 4.0 - d2 / 10.0  + d3 / 40.0  - d4 / 210.0;
		ehat = d / 4.0 - 3.0 * d2 / 20.0 + d3 / 20.0 - d4 / 84.0;
	}
	else
	{
		double ed = exp(-d), d2 = d * d, d3 = d2 * d;
		ahat = (-6.0 + 6.0 * d - 3.0 * d2 + d3 + 6.0 * ed) / d3;
		bhat = (6.0 - (6.0 + 6.0 * d + 3.0 * d2 + d3) * ed) / d3;
		ghat = 3.0 * (6.0 + d2 - 4.0 * d - (2.0 * d + 6.0) * ed) / d3;
		ehat = 3.0 * (2.0 * d - 6.0 + (6.0 + 4.0 * d + d2) * ed) / d3;
	}
}

static void formal_solution_bezier3(
	int ND, int NM, int NE,
	int i_inc,
	const double* mu,
	const double* ill_top,
	const double* ill_bot,
	const double* source,
	const double* dt_ang,
	double* intensity)
{
	memset(intensity, 0, long(ND) * NM * NE * sizeof(double));

	for (int nm = 0; nm < NM; ++nm)
	for (int ne = 0; ne < NE; ++ne)
	{
		// === Downward sweep (mu < 0) ===
		for (int nd = 0; nd < ND; ++nd)
		{
			if (mu[nm] >= 0.0) continue;

			if (nd == 0)
			{
				// orginal one 
				intensity[idx3(0, nm, ne)] =
					(nm == i_inc) ? ill_top[ne] : 0.0;
				continue;

				// this is change for test with pexrav (the isotropic incident)
				// intensity[idx3(0, nm, ne)] = ill_top[ne];
				continue;
			}

			double delta = dt_ang[idx3(nd, nm, ne)];
			double atten = exp(-delta);
			double S_D = source[idx3(nd, nm, ne)];
			double S_U = source[idx3(nd - 1, nm, ne)];

			double ah, bh, gh, eh;
			bezier3_coeffs(delta, ah, bh, gh, eh);

			// Control points E (near downwind) and F (near upwind)
			double E, F;
			bool ok_D = (nd > 0 && nd < ND - 1);
			bool ok_U = (nd - 1 > 0 && nd - 1 < ND - 1);

			if (ok_D && ok_U)
			{
				double hm_D = dt_ang[idx3(nd, nm, ne)];
				double hp_D = dt_ang[idx3(nd + 1, nm, ne)];
				double dSdt_D = fritsch_butland(
					source[idx3(nd - 1, nm, ne)], S_D,
					source[idx3(nd + 1, nm, ne)], hm_D, hp_D);

				double hm_U = dt_ang[idx3(nd - 1, nm, ne)];
				double hp_U = dt_ang[idx3(nd, nm, ne)];
				double dSdt_U = fritsch_butland(
					source[idx3(nd - 2, nm, ne)], S_U,
					source[idx3(nd, nm, ne)], hm_U, hp_U);

				// Incoming: h_k = -delta
				E = S_D - (delta / 3.0) * dSdt_D;
				F = S_U + (delta / 3.0) * dSdt_U;
			}
			else
			{
				// 1st-order fallback: control points at midpoints
				E = (2.0 * S_D + S_U) / 3.0;
				F = (S_D + 2.0 * S_U) / 3.0;
			}

			// Monotonicity clamp (each independently)
			double Smin = std::min(S_D, S_U);
			double Smax = std::max(S_D, S_U);
			E = std::max(Smin, std::min(Smax, E));
			F = std::max(Smin, std::min(Smax, F));

			intensity[idx3(nd, nm, ne)] =
				intensity[idx3(nd - 1, nm, ne)] * atten
				+ ah * S_D + bh * S_U + gh * E + eh * F;
		}

		// === Upward sweep (mu > 0) ===
		for (int nd = ND - 1; nd >= 0; --nd)
		{
			if (mu[nm] <= 0.0) continue;

			if (nd == ND - 1)
			{
				intensity[idx3(nd, nm, ne)] = ill_bot[ne];
				continue;
			}

			double delta = dt_ang[idx3(nd + 1, nm, ne)];
			double atten = exp(-delta);
			double S_D = source[idx3(nd, nm, ne)];
			double S_U = source[idx3(nd + 1, nm, ne)];

			double ah, bh, gh, eh;
			bezier3_coeffs(delta, ah, bh, gh, eh);

			// Control points E (near downwind) and F (near upwind)
			double E, F;
			bool ok_D = (nd > 0 && nd < ND - 1);
			bool ok_U = (nd + 1 > 0 && nd + 1 < ND - 1);

			if (ok_D && ok_U)
			{
				double hm_D = dt_ang[idx3(nd, nm, ne)];
				double hp_D = dt_ang[idx3(nd + 1, nm, ne)];
				double dSdt_D = fritsch_butland(
					source[idx3(nd - 1, nm, ne)], S_D,
					source[idx3(nd + 1, nm, ne)], hm_D, hp_D);

				double hm_U = dt_ang[idx3(nd + 1, nm, ne)];
				double hp_U = dt_ang[idx3(nd + 2, nm, ne)];
				double dSdt_U = fritsch_butland(
					source[idx3(nd, nm, ne)], S_U,
					source[idx3(nd + 2, nm, ne)], hm_U, hp_U);

				// Outgoing: h_k = +delta
				E = S_D + (delta / 3.0) * dSdt_D;
				F = S_U - (delta / 3.0) * dSdt_U;
			}
			else
			{
				E = (2.0 * S_D + S_U) / 3.0;
				F = (S_D + 2.0 * S_U) / 3.0;
			}

			double Smin = std::min(S_D, S_U);
			double Smax = std::max(S_D, S_U);
			E = std::max(Smin, std::min(Smax, E));
			F = std::max(Smin, std::min(Smax, F));

			intensity[idx3(nd, nm, ne)] =
				intensity[idx3(nd + 1, nm, ne)] * atten
				+ ah * S_D + bh * S_U + gh * E + eh * F;
		}
	}
}

// ============================================================
// mean_intensity
// ============================================================
static void mean_intensity(
	int ND, int NM, int NE,
	const double* wt,
	const double* intensity,
	double* meani)
{
	memset(meani, 0, long(ND) * NE * sizeof(double));

	for (int nd = 0; nd < ND; ++nd)
	for (int ne = 0; ne < NE; ++ne)
	{
		double s = 0.0;
		for (int nm = 0; nm < NM; ++nm)
			s += intensity[idx3(nd, nm, ne)] * wt[nm];
		meani[idx2(nd, ne)] = s * 0.5;
	}
}

// ============================================================
// check_rt_convergence
// ============================================================
static void check_rt_convergence(
	int ND, int NM, int NE,
	const double* ene, const double* mu, const double* wt,
	const double* meani, const double* meani_old,
	const double* intensity,
	const double* ill_top, const double* ill_bot,
	int i_inc,
	double& max_dJ)
{
	max_dJ = 0.0;
	long size2 = long(ND) * NE;
	for (long k = 0; k < size2; ++k)
	{
		if (meani[k] > 0.0)
		{
			double change = fabs(meani_old[k] / meani[k] - 1.0);
			if (change > max_dJ) max_dJ = change;
		}
	}
}

// ============================================================
// compute_dt_ang — angle-frequency optical depth
// ============================================================
static void compute_dt_ang(
	int ND, int NM, int NE,
	const RTGrids& g,
	const RadField& rad,
	double* dt_ang)
{
	for (int nd = 0; nd < ND; ++nd)
	for (int nm = 0; nm < NM; ++nm)
	for (int ne = 0; ne < NE; ++ne)
	{
		double dtau  = 0;
		// if (g.mu[nm]<0){
		// 	dtau = (rad.kabs[nd][ne] + rad.ksct[nd][ne]+rad.kabs_line[nd][ne]) * g.dr[nd];
		// } else{
		// 	dtau = (rad.kabs[nd][ne] + rad.ksct[nd][ne]) * g.dr[nd];
		// }
		dtau = (rad.kabs[nd][ne] + rad.ksct[nd][ne]) * g.dr[nd];
		dt_ang[idx3(nd, nm, ne)] = dtau / fabs(g.mu[nm]);
	}
}

// ============================================================
// compute_boundary_illumination
// ============================================================
static void compute_boundary_illumination(
	int NE, int i_inc,
	const RTGrids& g,
	const RadField& rad,
	double* ill_top, double* ill_bot)
{	
	// J (E) = 1/2 $$\int_{-1}^{1} $$ I(mu)delta(mu-mu_inc)dmu [Top]
	for (int ne = 0; ne < NE; ++ne)
	{	
		// original one 
		ill_top[ne] = 2*rad.illum.I_corona[ne]/g.wt[i_inc];

		// this is change for test with pexrav (the isotropic incident)
		// ill_top[ne] = 2*rad.illum.I_corona[ne];
		ill_bot[ne] = rad.illum.I_disk[ne] / 2;
	}
}

// ============================================================
// init_source_and_meani
// ============================================================
static void init_source_and_meani(
	int ND, int NM, int NE,
	const double* ill_top,
	double* source, double* meani, double* meani_old)
{
	for (int nd = 0; nd < ND; ++nd)
	for (int ne = 0; ne < NE; ++ne)
	{
		meani[idx2(nd, ne)]     = ill_top[ne];
		meani_old[idx2(nd, ne)] = ill_top[ne];
		for (int nm = 0; nm < NM; ++nm)
			source[idx3(nd, nm, ne)] = ill_top[ne];
	}
}

// ============================================================
// source_dispatch — route to the right source function by cache type
//   KernelCache    → directional compute_source_function (uses intensity)
//   avgKernelCache → angle-mean avgcompute_source_function (uses meani)
// ============================================================
static inline void source_dispatch(
	const KernelCache& kc, int ND, int NM, int NE,
	const double* intensity, const double* /*meani*/,
	const double* x_grid, const double* wmu, const double* T_K,
	const double* const* jnu, const double* const* kabs,
	const double* const* ksct, const double n_h, double* source)
{
	compute_source_function(ND, NM, NE, intensity, x_grid, wmu,
	                        kc, T_K, jnu, kabs, ksct, n_h, source);
}

static inline void source_dispatch(
	const avgKernelCache& kc, int ND, int NM, int NE,
	const double* /*intensity*/, const double* meani,
	const double* x_grid, const double* /*wmu*/, const double* T_K,
	const double* const* jnu, const double* const* kabs,
	const double* const* ksct, const double n_h, double* source)
{
	avgcompute_source_function(ND, NM, NE, meani, x_grid,
	                           kc, T_K, jnu, kabs, ksct, n_h, source);
}

// ============================================================
// compton_rt_solve — main driver
// ============================================================
template<class Cache>
void compton_rt_solve(RadField& rad, const RTGrids& g,
                      const ModelParams& par,
                      const Cache& kcache, int maxiter)
{
	const int ND = g.ND_MID;
	const int NM = g.NA;
	const int NE = g.NE;
	const double mec2 = phys::m_e_c2_eV;

	g_ND = ND;
	g_NM = NM;
	g_NE = NE;

	fprintf(stdout, "\n=== Compton RT solver (bezier3) ===\n");
	fprintf(stdout, "  ND=%d  NM=%d  NE=%d  maxiter=%d\n", ND, NM, NE, maxiter);

	double* x_grid = new double[NE];
	for (int ie = 0; ie < NE; ++ie)
		x_grid[ie] = g.ene[ie] / mec2;

	long size3 = long(ND) * NM * NE;
	long size2 = long(ND) * NE;

	double* source    = new double[size3]();
	double* intensity = new double[size3]();
	double* meani     = new double[size2]();
	double* meani_old = new double[size2]();
	double* dt_ang    = new double[size3]();

	double* ill_top = new double[NE]();
	double* ill_bot = new double[NE]();

	double mem_mb = (3.0 * size3 + 2.0 * size2 + 2.0 * NE) * 8.0 / (1024.0 * 1024.0);
	fprintf(stdout, "  Work arrays: %.1f MB\n", mem_mb);

	// --- Step 1: angle-frequency optical depth ---
	compute_dt_ang(ND, NM, NE, g, rad, dt_ang);

	// --- Step 2: boundary illumination ---
	int i_inc = par.i_incidence;
	compute_boundary_illumination(NE, i_inc, g, rad, ill_top, ill_bot);

	// --- Step 3: initialise source and mean intensity ---
	init_source_and_meani(ND, NM, NE, ill_top, source, meani, meani_old);

	// --- Step 4: Lambda iteration ---
	fprintf(stdout, "  Starting Lambda iteration (max %d)...\n", maxiter);
	const auto t_start = std::chrono::steady_clock::now();
	double formal_seconds=0.0, source_seconds=0.0;

	int iter = 0;
	int n_converged = 0;
	double max_change = 1.0;
	while ((n_converged < 3) && iter < maxiter)
	{
		const auto t_iter = std::chrono::steady_clock::now();

		// --- Formal solution ---
		formal_solution_bezier3(ND, NM, NE, i_inc, g.mu,
		                        ill_top, ill_bot,
		                        source, dt_ang,
		                        intensity);

		mean_intensity(ND, NM, NE, g.wt, intensity, meani);
		const auto t_source = std::chrono::steady_clock::now();
		formal_seconds += std::chrono::duration<double>(t_source-t_iter).count();

		source_dispatch(kcache, ND, NM, NE,
		                intensity, meani, x_grid, g.wt, rad.T_K,
		                rad.jnu, rad.kabs, rad.ksct, pow(10,par.nh),
		                source);
		source_seconds += std::chrono::duration<double>(std::chrono::steady_clock::now()-t_source).count();

		check_rt_convergence(ND, NM, NE, g.ene, g.mu, g.wt,
		                     meani, meani_old, intensity,
		                     ill_top, ill_bot, i_inc,
		                     max_change);

		if (max_change < 1e-4)
			++n_converged;
		else
			n_converged = 0;

		double dt_sec = std::chrono::duration<double>(std::chrono::steady_clock::now()-t_iter).count();
		fprintf(stdout, "    iter %3d  max|dJ/J|=%.6e  converged=%d/3  time=%.3fs\n",
		        iter + 1, max_change, n_converged, dt_sec);

		memcpy(meani_old, meani, size2 * sizeof(double));
		++iter;
	}

	fprintf(stdout, "  Converged at iteration %d.\n", iter);

	double total_sec = std::chrono::duration<double>(std::chrono::steady_clock::now()-t_start).count();
	fprintf(stdout, "  RT solver total time: %.2f s\n", total_sec);
	fprintf(stdout, "  RT wall-time breakdown: formal+mean=%.3fs source=%.3fs\n",
	        formal_seconds, source_seconds);

	// --- Copy results to RadField ---
	for (int nd = 0; nd < ND; ++nd)
	{
		for (int nm = 0; nm < NM; ++nm)
		for (int ne = 0; ne < NE; ++ne)
			rad.Inu[nd][nm][ne] = intensity[idx3(nd, nm, ne)];

		for (int ne = 0; ne < NE; ++ne)
			rad.J0[nd][ne] = meani[idx2(nd, ne)];
	}

	delete[] source;
	delete[] intensity;
	delete[] meani;
	delete[] meani_old;
	delete[] dt_ang;
	delete[] ill_top;
	delete[] ill_bot;
	delete[] x_grid;

	fprintf(stdout, "=== RT solver complete ===\n\n");
}

// Explicit instantiations for both kernel-cache types.
template void compton_rt_solve<KernelCache>(
	RadField&, const RTGrids&, const ModelParams&, const KernelCache&, int);
template void compton_rt_solve<avgKernelCache>(
	RadField&, const RTGrids&, const ModelParams&, const avgKernelCache&, int);
