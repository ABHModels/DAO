// ============================================================
// Compare exact vs approximate Compton kernel element
//
// Computes K(x, mu; x1, mu1, T) for ktype=0 (approx) and ktype=1 (exact)
// at fixed mu=0.5, x1=1keV/mec2, mu1=0.7, for T = 1e9, 1e8, 1e4 K.
// Scans x over [0.01, 100] * x1.
//
// Output: data/kernel_compare.dat
// ============================================================

#include "compton_kernel.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>

int main()
{
	const double mec2_eV = 511.0e3;
	const double mu  = 0.5;
	const double mu1 = 0.7;

	const double T_arr[] = {1.0e4, 1.0e5, 1.0e6};
	const int NT = 3;

	// Two input energies: 1 keV and 511 keV
	const double x1_eV[]  = {1.0e3, 511.0e3};
	const char*  x1_tag[] = {"1keV", "511keV"};
	const int N_X1 = 2;

	const int NX = 500;

	for (int ix1 = 0; ix1 < N_X1; ++ix1)
	{
		double x1 = x1_eV[ix1] / mec2_eV;
		double x_lo = 0.01 * x1;
		double x_hi = 100.0 * x1;

		char fname[256];
		snprintf(fname, sizeof(fname), "data/kernel_compare_%s.dat", x1_tag[ix1]);

		FILE* fp = fopen(fname, "w");
		if (!fp) { fprintf(stderr, "Cannot open %s\n", fname); return 1; }

		fprintf(fp, "# Compton kernel comparison: exact (ktype=1) vs approx (ktype=0)\n");
		fprintf(fp, "# x1 = %.6e (%s / mec2),  mu = %.2f,  mu1 = %.2f\n",
		        x1, x1_tag[ix1], mu, mu1);
		fprintf(fp, "# Col 1: x\n# Col 2: E [keV]\n");
		fprintf(fp, "# Col 3-5: K_exact  (T=1e4, 1e5, 1e6 K)\n");
		fprintf(fp, "# Col 6-8: K_approx (T=1e4, 1e5, 1e6 K)\n#\n");

		fprintf(fp, "# %14s %14s %14s %14s %14s %14s %14s %14s\n",
		        "x", "E_keV",
		        "K_ex_1e4", "K_ex_1e5", "K_ex_1e6",
		        "K_ap_1e4", "K_ap_1e5", "K_ap_1e6");

		for (int ix = 0; ix < NX; ++ix)
		{
			double x = x_lo * pow(x_hi / x_lo, (double)ix / (NX - 1));
			double E_keV = x * mec2_eV / 1.0e3;

			fprintf(fp, "%16.8e %16.8e", x, E_keV);

			for (int iT = 0; iT < NT; ++iT)
				fprintf(fp, " %16.8e",
				        compton_kernel_element(x, mu, x1, mu1, T_arr[iT], 1));

			for (int iT = 0; iT < NT; ++iT)
				fprintf(fp, " %16.8e",
				        compton_kernel_element(x, mu, x1, mu1, T_arr[iT], 0));

			fprintf(fp, "\n");
		}

		fclose(fp);
		fprintf(stdout, "Written %s (%d points, %d temperatures)\n", fname, NX, NT);
	}

	return 0;
}