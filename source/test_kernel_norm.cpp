// ============================================================
// A23 normalization test for the Compton redistribution kernel
//
// Verifies the photon-conservation identity (Poutanen 1994, Eq. 3.19;
// Madej et al. 2017, Eq. A23):
//
//   sigma(x, T) / sigma_T
//     = (1/x) ∫₀^∞ x1 dx1  ∫₀¹ dμ  ∫₀¹ dμ1
//             [ R(x1, μ1;  x,  μ; T) + R(x1, μ1;  x, -μ; T) ]
//
//   LHS  →  compton_cross_section.cpp   (analytic σ / σ_T)
//   RHS  →  compton_kernel.cpp          (integrated redistribution)
//
// Energy grid: log-spaced from E_LO to E_HI (wide enough at high T).
// Angle grid : Gauss-Legendre on [-1, 1] from RTGrids (positive half).
// Frequency integral: trapezoidal in x1.
//
// Output: data/kernel_norm_T{tag}.dat — one file per temperature.
//   Col 1: E [keV]
//   Col 2: sigma_LHS / sigma_T          (analytic)
//   Col 3: sigma_RHS / sigma_T          (A23 integral over kernel)
//   Col 4: ratio = RHS / LHS
// ============================================================

#include "compton_kernel.h"
#include "compton_cross_section.h"
#include "rt_grids.h"
#include "constants.h"
#include <cmath>
#include <cstdio>

int main()
{
	// ---- Configuration ----
	const int    KTYPE  = 1;        // 1 = exact Madej+2017 profile
	const int    NE     = 2000;     // energy samples
	const double E_LO   = 10.0;     // eV
	const double E_HI   = 1.0e7;    // 10 MeV
	const double T_arr[] = { 1.0e4, 1.0e7, 1.0e9 };
	const char*  T_tag[] = { "1e4", "1e7", "1e9" };
	const int    NT      = sizeof(T_arr) / sizeof(T_arr[0]);

	const double mec2_eV = phys::m_e_c2_eV;

	// ---- Angle grid (Gauss-Legendre on [-1, 1], 8 points) ----
	RTGrids g;
	g.init_angle();
	const int NA = g.NA;

	int mu_neg[RTGrids::NA];
	for (int i = 0; i < NA; ++i) mu_neg[i] = NA - 1 - i;

	// ---- Energy grid (log-spaced) ----
	double ene_eV[NE], x_grid[NE];
	for (int i = 0; i < NE; ++i)
	{
		ene_eV[i] = E_LO * pow(E_HI / E_LO, double(i) / (NE - 1));
		x_grid[i] = ene_eV[i] / mec2_eV;
	}

	// ---- A23 test for each temperature ----
	for (int iT = 0; iT < NT; ++iT)
	{
		const double T = T_arr[iT];

		char fname[256];
		snprintf(fname, sizeof(fname), "data/kernel_norm_T%s.dat", T_tag[iT]);
		FILE* fp = fopen(fname, "w");
		if (!fp) { fprintf(stderr, "Cannot open %s\n", fname); return 1; }

		fprintf(fp, "# A23 normalization test, T = %s K, ktype = %d\n", T_tag[iT], KTYPE);
		fprintf(fp, "# LHS: compton_cross_section / sigma_T\n");
		fprintf(fp, "# RHS: (1/x) ∫ x1 dx1 ∫∫ dμ dμ1 [R(x1,μ1;x,μ) + R(x1,μ1;x,-μ)]\n");
		fprintf(fp, "# %14s %16s %16s %12s\n",
		        "E_keV", "sigma_LHS", "sigma_RHS", "RHS/LHS");

		fprintf(stdout, "T = %s K (%d/%d):\n", T_tag[iT], iT + 1, NT);

		for (int ne = 0; ne < NE; ++ne)
		{
			if (ne % 50 == 0)
				fprintf(stdout, "  E = %.2e keV  (%d/%d)\n",
				        ene_eV[ne] / 1e3, ne + 1, NE);

			const double x = x_grid[ne];

			// --- RHS: triple integral over (x1, μ, μ1) ---
			double integral = 0.0;

			for (int ne1 = 1; ne1 < NE; ++ne1)
			{
				const double dx1 = x_grid[ne1] - x_grid[ne1 - 1];

				// Angle sum at the two x1 endpoints
				double f_lo = 0.0, f_hi = 0.0;

				for (int nm = 0; nm < NA; ++nm)
				{
					if (g.mu[nm] <= 0.0) continue;       // positive μ half
					for (int nm1 = 0; nm1 < NA; ++nm1)
					{
						if (g.mu[nm1] <= 0.0) continue;  // positive μ1 half
						const double w = g.wt[nm] * g.wt[nm1];

						// Pair [+μ, -μ] for the reflected hemisphere
						const double K_lo =
						    compton_kernel_element(x_grid[ne1-1], g.mu[nm1],
						                           x, g.mu[nm], T, KTYPE)
						  + compton_kernel_element(x_grid[ne1-1], g.mu[nm1],
						                           x, g.mu[mu_neg[nm]], T, KTYPE);
						const double K_hi =
						    compton_kernel_element(x_grid[ne1],   g.mu[nm1],
						                           x, g.mu[nm], T, KTYPE)
						  + compton_kernel_element(x_grid[ne1],   g.mu[nm1],
						                           x, g.mu[mu_neg[nm]], T, KTYPE);

						f_lo += w * K_lo;
						f_hi += w * K_hi;
					}
				}

				// Trapezoidal step in x1, with the x1 weight inside the integrand
				integral += 0.5 * (f_lo * x_grid[ne1 - 1] + f_hi * x_grid[ne1]) * dx1;
			}

			const double sigma_rhs = integral / x;                                   // RHS
			const double sigma_lhs =
			    compton_cross_section(ene_eV[ne], T) / phys::sigma_T;               // LHS
			const double ratio     = (sigma_lhs > 1e-30) ? sigma_rhs / sigma_lhs : 0.0;

			fprintf(fp, "%16.8e %16.8e %16.8e %12.6f\n",
			        ene_eV[ne] / 1e3, sigma_lhs, sigma_rhs, ratio);
		}

		fclose(fp);
		fprintf(stdout, "  wrote %s\n", fname);
	}

	return 0;
}
