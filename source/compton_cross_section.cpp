#include "compton_cross_section.h"
#include "constants.h"
#include "rt_grids.h"
#include <cmath>
#include <cstdio>
#include <cstring>

// ============================================================
// Modified Bessel function K_2(x) * exp(x)
//
// Uses polynomial approximations from Abramowitz & Stegun,
// with the recurrence K_{n+1}(x) = (2n/x) K_n(x) + K_{n-1}(x).
// Translated from bk2.f (Madej 1989).
// ============================================================
static double bk2_exp(double x)
{
	double k0, k1;

	if (x <= 2.0)
	{
		double z = (x / 2.0) * (x / 2.0);
		double t = (x / 3.75) * (x / 3.75);

		double i0 = (((((0.0045813 * t + 0.0360768) * t + 0.2659732) * t
		            + 1.2067492) * t + 3.0899424) * t + 3.5156229) * t + 1.0;
		double i1 = (((((.00032411 * t + .00301532) * t + .02658733) * t
		            + .15084934) * t + .51498869) * t + .87890594) * t + 0.5;
		i1 *= x;

		k0 = (((((.00000740 * z + .00010750) * z + .00262698) * z
		     + .03488590) * z + .23069756) * z + .42278420) * z
		     - 0.57721566 - i0 * log(x / 2.0);
		k1 = (((((-.00004686 * z - .00110404) * z - .01919402) * z
		     - .18156897) * z - .67278579) * z + .15443144) * z + 1.0;
		k1 = k1 / x + i1 * log(x / 2.0);

		// Recurrence: K_2(x) = (2/x) K_1(x) + K_0(x)
		double k2 = k0 + 2.0 / x * k1;
		return k2 * exp(x);
	}
	else if (x < 500.0)
	{
		double z = 2.0 / x;
		k0 = (((((0.00053208 * z - 0.00251540) * z + 0.00587872) * z
		     - 0.01062446) * z + 0.02189568) * z - 0.07832358) * z
		     + 1.25331414;
		k1 = (((((-0.00068245 * z + 0.00325614) * z - 0.00780353) * z
		     + 0.01504268) * z - 0.03655620) * z + 0.23498619) * z
		     + 1.25331414;
		// Already scaled by exp(x) in this branch
		return (k0 + z * k1) / sqrt(x);
	}
	else
	{
		double z = 1.0 / x;
		return sqrt(M_PI / 2.0) * (sqrt(z) + 1.875 * pow(z, 1.5)
		     + 0.8203125 * pow(z, 2.5));
	}
}

// ============================================================
// Gauss-Legendre quadrature nodes and weights on [x1, x2]
// Translated from gaulegf.f (Numerical Recipes).
// ============================================================
static void gauleg(double x1, double x2, double* xn, double* wn, int n)
{
	const double eps = 3.0e-14;
	double xm = 0.5 * (x2 + x1);
	double xl = 0.5 * (x2 - x1);

	for (int i = 0; i < (n + 1) / 2; ++i)
	{
		double z = cos(M_PI * (i + 0.75) / (n + 0.5));
		double z1 = 0.0;
		double pp;

		while (fabs(z - z1) > eps)
		{
			double p1 = 1.0, p2 = 0.0;
			for (int j = 1; j <= n; ++j)
			{
				double p3 = p2;
				p2 = p1;
				p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j;
			}
			pp = n * (z * p1 - p2) / (z * z - 1.0);
			z1 = z;
			z = z1 - p1 / pp;
		}

		xn[i]       = xm - xl * z;
		xn[n-1-i]   = xm + xl * z;
		wn[i]       = 2.0 * xl / ((1.0 - z * z) * pp * pp);
		wn[n-1-i]   = wn[i];
	}
}

// ============================================================
// Integrand of the exact relativistic Compton cross section
// for a single Lorentz factor γ.
// Translated from crsexact.f: function crsgm(theta, eps, gm)
// ============================================================
static double crsgm(double theta, double eps, double gm)
{
	const int nksi = 20;
	double ksi[20], wksi[20];

	double z = sqrt(gm * gm - 1.0);

	// Integration limits over ksi = ε(γ ± β)
	double a = eps * (gm - z);
	double b = eps * (gm + z);
	gauleg(a, b, ksi, wksi, nksi);

	// Integral over ksi
	double integ = 0.0;
	for (int i = 0; i < nksi; ++i)
		integ += wksi[i] * log(1.0 + 2.0 * ksi[i]) / ksi[i];

	double result = (eps * gm + 4.5 + 2.0 * gm / eps)
	              * log((1.0 + 2.0 * eps * (gm + z)) / (1.0 + 2.0 * eps * (gm - z)))
	              - 2.0 * eps * z + z * (eps - 2.0 / eps)
	              * log(1.0 + 4.0 * eps * gm + 4.0 * eps * eps)
	              + (4.0 * eps * eps * z * (gm + eps))
	              / (1.0 + 4.0 * eps * gm + 4.0 * eps * eps) - 2.0 * integ;

	// Multiply by relativistic Maxwellian factor exp((1-γ)/θ)
	double x_inv = 1.0 / theta;
	result *= exp((1.0 - gm) * x_inv);

	return result;
}

// ============================================================
// Exact relativistic Compton scattering cross section
// averaged over Maxwellian electrons.
// Poutanen & Svensson (1996).
//
// theta = kT / (m_e c^2),  eps = E / (m_e c^2)
// Returns cross section in [cm^2].
//
// Translated from crsexact.f, 
// ============================================================
static double crsexact(double theta, double eps)
{
	if (eps <= 1.0e-4)
		return phys::sigma_T;

	const int gi = 200;
	double g_m1[200], aj[200];

	double mn = 3.0 * phys::sigma_T / (16.0 * eps * eps * theta);

	// Integration over Lorentz factor γ in [1, 20]
	gauleg(1.0, 20.0, g_m1, aj, gi);

	double totalcrs = 0.0;
	for (int i = 0; i < gi; ++i)
		totalcrs += aj[i] * crsgm(theta, eps, g_m1[i]);

	double x_inv = 1.0 / theta;
	totalcrs = mn * totalcrs / bk2_exp(x_inv);

	return totalcrs;
}

// ============================================================
// Compton scattering cross section σ(E, T) [cm^2]
//
// Three regimes (matching scattxs.f):
//   1. Low T, low E: Thomson
//   2. Low T, high E: Klein-Nishina
//   3. High T: exact relativistic (Poutanen & Svensson 1996)
// ============================================================
double compton_cross_section(double E_eV, double T_K)
{
	const double mec2  = phys::m_e_c2_eV;           // 511 keV in eV
	const double kb_eV = phys::k_B / phys::eV_to_erg;  // Boltzmann in eV/K

	double theta = kb_eV * T_K / mec2;   // dimensionless temperature kT/(m_e c^2)

	// xloc = 2 * E / (m_e c^2)  — convention used in scattxs.f
	// 3.913894e-6 = 2 / 511000
	double xloc  = 2.0 * E_eV / mec2;
	double ixloc = 1.0 / xloc;

	if (theta < 0.0169)  // T < ~10^8 K
	{
		if (xloc < 0.002)   // E < ~0.5 keV: Thomson limit
		{
			return phys::sigma_T;
		}
		else
		{
			// Klein-Nishina total cross section
			// σ_KN = (3/4)σ_T × f(x),  where x = 2E/(m_e c^2)
			return 4.9875e-25 * ((1.0 - 4.0 * ixloc - 8.0 * ixloc * ixloc)
			     * log(1.0 + xloc) + 0.5 + 8.0 * ixloc
			     - 0.5 / ((1.0 + xloc) * (1.0 + xloc))) * ixloc;
		}
	}
	else  // T ≥ ~10^8 K: exact relativistic
	{
		double eps = E_eV / mec2;   // dimensionless photon energy
		return crsexact(theta, eps);
	}
}

// ============================================================
// Compute Compton scattering opacity for all energies at one
// temperature and electron density:
//
//   ksct[ie] = n_e × σ_compton(E[ie], T)   [cm^-1]
//
// where σ_compton accounts for Klein-Nishina and relativistic
// thermal corrections (see compton_cross_section above).
//
// Cross section code translated from:
//   based on Poutanen & Svensson (1996)
// ============================================================
void compute_compton_opacity(double* ksct, int NE, const double* ene_eV,
                             double T_K, double n_e)
{
	for (int ie = 0; ie < NE; ++ie)
		ksct[ie] = n_e * compton_cross_section(ene_eV[ie], T_K);
}



