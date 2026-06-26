// ============================================================
// dump_kernel_slice — emit a 2D angular×energy slice of the
// Compton redistribution kernel K(mu_out, E_out) at fixed
// (E_in, mu_in, T).
//
// Calls compton_kernel_element() directly. At T~10^8 K this
// already satisfies the A23 sum rule to ~0.5%, so no extra
// normalization is applied here.
//
// Output:
//   data/kernel_redist2d_Ein{}_T{}_mu{}.dat
//
// Header lines (all start with "#"):
//   T, E_in, mu_in, NMU, NEO, mu_out_grid, Eout_grid_eV
// then NMU rows × NEO cols of K values (rows = mu_out).
//
// Usage:
//   ./dump_kernel_slice -Ein 50e3 -T 1e8 -mu_in 0.96 \
//                       -Eo_lo 5 -Eo_hi 500 -NMU 181 -NEO 240
//
// Author:      Yimin Huang
// Affiliation: Fudan University
// Email:       huangym23@m.fudan.edu.cn
// ============================================================

#include "compton_kernel.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

int main(int argc, char* argv[])
{
	const char* out_path = nullptr;
	double E_in       = 50e3;
	double T          = 1e8;
	double mu_in      = 0.96;
	double Eo_lo_keV  = 5.0;
	double Eo_hi_keV  = 500.0;
	int    NMU        = 181;
	int    NEO        = 240;

	for (int i = 1; i < argc; ++i)
	{
		if      (!strcmp(argv[i], "-out")    && i+1 < argc) out_path  = argv[++i];
		else if (!strcmp(argv[i], "-Ein")    && i+1 < argc) E_in      = atof(argv[++i]);
		else if (!strcmp(argv[i], "-T")      && i+1 < argc) T         = atof(argv[++i]);
		else if (!strcmp(argv[i], "-mu_in")  && i+1 < argc) mu_in     = atof(argv[++i]);
		else if (!strcmp(argv[i], "-Eo_lo")  && i+1 < argc) Eo_lo_keV = atof(argv[++i]);
		else if (!strcmp(argv[i], "-Eo_hi")  && i+1 < argc) Eo_hi_keV = atof(argv[++i]);
		else if (!strcmp(argv[i], "-NMU")    && i+1 < argc) NMU       = atoi(argv[++i]);
		else if (!strcmp(argv[i], "-NEO")    && i+1 < argc) NEO       = atoi(argv[++i]);
		else { fprintf(stderr, "unknown arg: %s\n", argv[i]); return 1; }
	}

	const double mec2_eV = 511.0e3;
	double x_in = E_in / mec2_eV;

	fprintf(stdout, "E_in=%.3e eV, T=%.3e K, mu_in=%.3f\n", E_in, T, mu_in);
	fprintf(stdout, "E_out grid: [%.1f, %.1f] keV, NEO=%d\n", Eo_lo_keV, Eo_hi_keV, NEO);
	fprintf(stdout, "mu_out grid: [-1, +1], NMU=%d\n", NMU);

	// --- (mu_out, E_out) grids ---
	std::vector<double> mu_out(NMU), Eo(NEO);
	for (int i = 0; i < NMU; ++i)
		mu_out[i] = -1.0 + 2.0 * i / (NMU - 1);
	double l_lo = log10(Eo_lo_keV * 1e3);
	double l_hi = log10(Eo_hi_keV * 1e3);
	for (int j = 0; j < NEO; ++j)
		Eo[j] = pow(10.0, l_lo + (l_hi - l_lo) * j / (NEO - 1));

	// --- Output ---
	char fname[1024];
	if (out_path) {
		snprintf(fname, sizeof(fname), "%s", out_path);
	} else {
		snprintf(fname, sizeof(fname),
		         "data/kernel_redist2d_Ein%.0fkeV_T%.0e_mu%.2f.dat",
		         E_in/1e3, T, mu_in);
	}
	FILE* fp = fopen(fname, "w");
	if (!fp) { fprintf(stderr, "cannot write %s\n", fname); return 1; }

	fprintf(fp, "# Compton redistribution K(mu_out, E_out) at fixed E_in, mu_in, T\n");
	fprintf(fp, "# T            = %.6e K\n",  T);
	fprintf(fp, "# E_in         = %.6e eV\n", E_in);
	fprintf(fp, "# mu_in        = %.6f\n", mu_in);
	fprintf(fp, "# NMU          = %d\n", NMU);
	fprintf(fp, "# NEO          = %d\n", NEO);
	fprintf(fp, "# mu_out_grid  =");
	for (int i = 0; i < NMU; ++i) fprintf(fp, " %+.6f", mu_out[i]);
	fprintf(fp, "\n");
	fprintf(fp, "# Eout_grid_eV =");
	for (int j = 0; j < NEO; ++j) fprintf(fp, " %.6e", Eo[j]);
	fprintf(fp, "\n");
	fprintf(fp, "# data block: %d rows (mu_out)  x  %d cols (E_out)\n", NMU, NEO);

	for (int i = 0; i < NMU; ++i)
	{
		double mu_clip = mu_out[i];
		if (mu_clip > 0.99999)  mu_clip = 0.99999;
		if (mu_clip < -0.99999) mu_clip = -0.99999;
		for (int j = 0; j < NEO; ++j)
		{
			double xo = Eo[j] / mec2_eV;
			double v  = compton_kernel_element(
				x_in, mu_in, xo, mu_clip, T, /*ktype*/1);
			fprintf(fp, "%.6e ", v);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	fprintf(stdout, "wrote %s\n", fname);
	return 0;
}