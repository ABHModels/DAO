// ============================================================
// Test normalization condition (Poutanen 1994, Eq. 3.19):
//
// s(x) = (1/x) int_0^inf x1 dx1 int_0^1 deta int_0^1 deta1
//            [R(x1,eta1; x,eta) + R(x1,eta1; x,-eta)]
//
// This should equal sigma_compton(x,T) / sigma_T.
// Tests both exact (ktype=1) and approximate (ktype=0) kernels.
// Output: data/kernel_norm_T{}.dat for T = 1e4, 1e7, 1e9 K
// ============================================================

#include "compton_kernel.h"
#include "compton_cross_section.h"
#include "rt_grids.h"
#include "constants.h"
#include <cmath>
#include <cstdio>

int main()
{
	const double mec2_eV = phys::m_e_c2_eV;

	// Setup angle grid
	RTGrids g;
	g.init_angle();
	const int NA = g.NA;

	int nm_neg[RTGrids::NA];
	for (int nm = 0; nm < NA; ++nm)
		nm_neg[nm] = NA - 1 - nm;

	// Energy grid: log-spaced from 10 eV to 10 MeV (wide enough for T=1e9)
	const int NE = 2000;
	double ene_eV[NE], x_grid[NE];
	for (int i = 0; i < NE; ++i) {
		ene_eV[i] = 10.0 * pow(1.0e7 / 10.0, (double)i / (NE - 1));  // 10 eV to 10 MeV
		x_grid[i] = ene_eV[i] / mec2_eV;
	}

	const double T_arr[] = {1.0e4, 1.0e7, 1.0e9};
	const char*  T_tag[] = {"1e4", "1e7", "1e9"};
	const int NT = 3;

	for (int iT = 0; iT < NT; ++iT)
	{
		double T = T_arr[iT];

		char fname[256];
		snprintf(fname, sizeof(fname), "data/kernel_norm_T%s.dat", T_tag[iT]);

		FILE* fp = fopen(fname, "w");
		if (!fp) { fprintf(stderr, "Cannot open %s\n", fname); return 1; }

		fprintf(fp, "# A23 normalization test for T = %s K\n", T_tag[iT]);
		fprintf(fp, "# Col 1: E [keV]\n");
		fprintf(fp, "# Col 2: sigma_exact / sigma_T\n");
		fprintf(fp, "# Col 3: sigma_raw_exact / sigma_T  (ktype=1)\n");
		fprintf(fp, "# Col 4: sigma_raw_approx / sigma_T (ktype=0)\n");
		fprintf(fp, "# Col 5: ratio exact  (raw/exact)\n");
		fprintf(fp, "# Col 6: ratio approx (raw/exact)\n");
		fprintf(fp, "#\n");
		fprintf(fp, "# %14s %14s %14s %14s %14s %14s\n",
		        "E_keV", "sig_exact", "sig_raw_ex", "sig_raw_ap", "ratio_ex", "ratio_ap");

		fprintf(stdout, "T = %s K: computing A23 integral for %d energies...\n", T_tag[iT], NE);

		for (int ne = 0; ne < NE; ++ne)
		{
			if (ne % 50 == 0)
				fprintf(stdout, "  E = %.1f keV (%d/%d)\n", ene_eV[ne]/1e3, ne+1, NE);

			double x = x_grid[ne];

			// A23 integral for both ktypes
			double integral_ex = 0.0;
			double integral_ap = 0.0;

			for (int ne1 = 1; ne1 < NE; ++ne1)
			{
				double dx1 = x_grid[ne1] - x_grid[ne1 - 1];

				double f_a_ex = 0.0, f_b_ex = 0.0;
				double f_a_ap = 0.0, f_b_ap = 0.0;

				for (int inm = 0; inm < NA; ++inm)
				{
					if (g.mu[inm] <= 0.0) continue;
					for (int inm1 = 0; inm1 < NA; ++inm1)
					{
						if (g.mu[inm1] <= 0.0) continue;
						double ww = g.wt[inm] * g.wt[inm1];

						// Exact kernel (ktype=1)
						double Ka_ex = compton_kernel_element(x_grid[ne1-1], g.mu[inm1], x, g.mu[inm], T, 1)
						             + compton_kernel_element(x_grid[ne1-1], g.mu[inm1], x, g.mu[nm_neg[inm]], T, 1);
						double Kb_ex = compton_kernel_element(x_grid[ne1], g.mu[inm1], x, g.mu[inm], T, 1)
						             + compton_kernel_element(x_grid[ne1], g.mu[inm1], x, g.mu[nm_neg[inm]], T, 1);
						f_a_ex += ww * Ka_ex;
						f_b_ex += ww * Kb_ex;

						// Approx kernel (ktype=0)
						double Ka_ap = compton_kernel_element(x_grid[ne1-1], g.mu[inm1], x, g.mu[inm], T, 0)
						             + compton_kernel_element(x_grid[ne1-1], g.mu[inm1], x, g.mu[nm_neg[inm]], T, 0);
						double Kb_ap = compton_kernel_element(x_grid[ne1], g.mu[inm1], x, g.mu[inm], T, 0)
						             + compton_kernel_element(x_grid[ne1], g.mu[inm1], x, g.mu[nm_neg[inm]], T, 0);
						f_a_ap += ww * Ka_ap;
						f_b_ap += ww * Kb_ap;
					}
				}

				f_a_ex *= x_grid[ne1 - 1];
				f_b_ex *= x_grid[ne1];
				integral_ex += 0.5 * (f_a_ex + f_b_ex) * dx1;

				f_a_ap *= x_grid[ne1 - 1];
				f_b_ap *= x_grid[ne1];
				integral_ap += 0.5 * (f_a_ap + f_b_ap) * dx1;
			}

			double sigma_raw_ex = integral_ex / x;
			double sigma_raw_ap = integral_ap / x;
			double sigma_exact  = compton_cross_section(ene_eV[ne], T) / phys::sigma_T;

			double ratio_ex = (sigma_exact > 1e-30) ? sigma_raw_ex / sigma_exact : 0.0;
			double ratio_ap = (sigma_exact > 1e-30) ? sigma_raw_ap / sigma_exact : 0.0;

			fprintf(fp, "%16.8e %16.8e %16.8e %16.8e %16.8e %16.8e\n",
			        ene_eV[ne] / 1e3, sigma_exact, sigma_raw_ex, sigma_raw_ap,
			        ratio_ex, ratio_ap);
		}

		fclose(fp);
		fprintf(stdout, "Written %s\n", fname);
	}

	return 0;
}