#include "compton_kernel.h"
#include "compton_cross_section.h"
#include "kernel_row_spool.h"
#include "rt_parallel.h"
#include "constants.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

// ============================================================
// Modified Bessel function K₂(x) × exp(x)
//
// Polynomial approximations from Abramowitz & Stegun, with
// recurrence K_{n+1}(x) = (2n/x) K_n(x) + K_{n-1}(x).
// Translated from bk2.f (Madej 1989).
//
// Identical to the version in compton_cross_section.cpp;
// duplicated here to keep the kernel module self-contained.
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
// Exact differential cross section (Madej et al. 2017)
//
// Translated from fexact() in profile.f (Prof. Jerzy Madej,
// University of Warsaw).
//
// eps, eps1 : photon energies in units of m_e c²
// costh     : cosine of scattering angle
// gamma     : electron Lorentz factor
// ============================================================
static double fexact(double eps, double eps1, double costh, double gamma)
{
	double q_down = eps * eps1 * (1.0 - costh);
	double q_up   = sqrt((eps - eps1) * (eps - eps1) + 2.0 * q_down);

	double cot_factor = (1.0 + costh) / (1.0 - costh);
	double a2_minus = (gamma - eps)  * (gamma - eps)  + cot_factor;
	double a2_plus  = (gamma + eps1) * (gamma + eps1) + cot_factor;

	double a_minus = sqrt(a2_minus);
	double a_plus  = sqrt(a2_plus);

	// Numerically stable form of (a_minus - a_plus)
	double a2_diff = -2.0 * gamma * (eps + eps1) + eps * eps - eps1 * eps1;
	double a_diff  = a2_diff / (a_minus + a_plus);

	// Madej et al. (2017) formula
	double term1 = 2.0 / q_up;
	double term2 = -(q_down - 2.0) / (q_down * a_plus * a_minus) * a_diff;
	double term3_num = (a2_minus + a_minus * a_plus + a2_plus) * q_up * q_up
	                 - a2_diff * a2_diff
	                 - a_minus * a_plus * a_diff * a_diff;
	double term3_den = 2.0 * q_down * q_down
	                 * a_plus * a_plus * a_plus
	                 * a_minus * a_minus * a_minus;
	double term3 = term3_num / term3_den * a_diff;

	return term1 + term2 + term3;
}

// ============================================================
// Redistribution function R(ε, ε₁, cosθ; 1/Θ)
//
// Exact form (nr=1) from profile.f (Prof. Jerzy Madej).
// Uses 32-point Gauss-Laguerre quadrature over the electron
// Lorentz factor γ.
//
// The integration variable is transformed as γ = t/x + γ*,
// where t are the GL nodes and x = m_e c²/(kT).
// The GL weight function exp(-t) maps to the Maxwellian
// exp(-x(γ - γ*)) in the Boltzmann factor.
// eps, eps1 : photon energies in units of m_e c²
// costh     : cosine of scattering angle
// x_inv     : m_e c² / (kT)  (inverse dimensionless temperature)
// ============================================================

// 32-point Gauss-Laguerre nodes and weights
static const int NGL = 32;
static const double gl_node[NGL] = {
	0.044489366e0,  0.23452611e0,   0.576884629e0,  1.07244875e0,
	1.72240878e0,   2.528336706e0,  3.492213273e0,  4.61645677e0,
	5.9039585e0,    7.358126733e0,  8.98294092e0,   10.78301863e0,
	12.76369799e0,  14.93113976e0,  17.29245434e0,  19.85586094e0,
	22.63088901e0,  25.62863602e0,  28.86210182e0,  32.34662915e0,
	36.10049481e0,  40.1457198e0,   44.509208e0,    49.22439499e0,
	54.3337213e0,   59.89250916e0,  65.97537729e0,  72.68762809e0,
	80.18744698e0,  88.73534042e0,  98.82954287e0,  111.7513981e0
};
static const double gl_weight[NGL] = {
	1.09218342e-1,  2.10443108e-1,  2.352132297e-1, 1.95903336e-1,
	1.29983786e-1,  7.05786239e-2,  3.17609125e-2,  1.191821484e-2,
	3.738816295e-3, 9.8080331e-4,   2.14864919e-4,  3.920342e-5,
	5.9345416e-6,   7.4164046e-7,   7.6045679e-8,   6.3506022e-9,
	4.28138297e-10, 2.30589949e-11, 9.7993793e-13,  3.2378017e-14,
	8.17182344e-16, 1.54213383e-17, 2.11979229e-19, 2.05442967e-21,
	1.34698259e-23, 5.6612941e-26,  1.41856055e-28, 1.91337549e-31,
	1.19224876e-34, 2.67151122e-38, 1.33861694e-42, 4.5105362e-48
};

double profil_exact(double eps, double eps1, double costh, double x_inv)
{
	// Minimum Lorentz factor for the scattering kinematics
	double q_down   = eps * eps1 * (1.0 - costh);
	double q_up     = sqrt((eps - eps1) * (eps - eps1) + 2.0 * q_down);
	double gam_star = (eps - eps1 + q_up * sqrt(1.0 + 2.0 / q_down)) / 2.0;

	// Gauss-Laguerre quadrature: ∫₀^∞ f(γ) exp(-t) dt
	// with γ = t/x_inv + γ*
	double integ = 0.0;
	for (int j = 0; j < NGL; ++j)
	{
		double gamma = gl_node[j] / x_inv + gam_star;
		integ += gl_weight[j] * fexact(eps, eps1, costh, gamma);
	}

	// bk2_exp(x) = K₂(x) × exp(x), so:
	//   3/(32π) × integ × exp((1-γ*)x) / [K₂(x)×exp(x)]
	// = 3/(32π) × integ × exp((1-γ*)x) / bk2_exp(x)
	// double result = (3.0 / (32.0 * M_PI)) * integ * eps1/eps
	//               * exp((1.0 - gam_star+eps-eps1) * x_inv) / bk2_exp(x_inv);
	double result = (3.0 / (32.0 * M_PI)) * integ
	              * exp((1.0 - gam_star) * x_inv) / bk2_exp(x_inv);

	return result;
}

/// Approximate Compton scattering redistribution function
/// (Poutanen 1994, PhD thesis, Eq. 2.23).
///
/// R(ε, ε₁, cosθ) = [1 / (8π Q)] × exp(−y γ*) / K₂(y)
///
/// where:
///   ε, ε₁  = photon energies in units of m_e c²
///   cosθ    = scattering angle cosine
///   y       = m_e c² / kT  (inverse dimensionless temperature)
///   Q       = √[(ε − ε₁)² + 2εε₁(1 − cosθ)]
///   γ*      = Q / √[2εε₁(1 − cosθ)]   (minimum Lorentz factor)
///
/// Good approximation for weakly/mildly-relativistic temperatures (kT_e ≤ m_e c²).
double profil_exact_ap(double eps, double eps1, double costh, double x_inv)
{
	double q_down   = eps * eps1 * (1.0 - costh);
	double q_up     = sqrt((eps - eps1) * (eps - eps1) + 2.0 * q_down);  // Q
	double gam_star = q_up / sqrt(2.0 * q_down);                         // γ*

	// Use bk2_exp(y) = K₂(y)·exp(y) to avoid overflow:
	//   exp(−y·γ*) / K₂(y) = exp(−y·γ* + y) / bk2_exp(y)
	double result = 1.0 / (8.0 * M_PI) / q_up
	              * exp(-x_inv * gam_star + x_inv) / bk2_exp(x_inv);

	return result;
}

// ============================================================
// Azimuth-integrated Compton kernel for one pair (x,μ)→(x₁,μ₁)
//
// K(x, μ; x₁, μ₁) = ∫₀²π R(x, x₁, η(μ,μ₁,φ)) dφ
//
// where η = μμ₁ + √(1-μ²)√(1-μ₁²) cos(φ).
// 10-point Gauss-Legendre quadrature on [0, 2π].
//
// Nodes and weights from radiative_transfer_sec.f90 (ckernel).
// ============================================================
static const int NPHI = 6;
static const double phi_node[NPHI] = {
	2.1215327807272738e-01, 1.0643421025827622e+00, 2.3919483715852459e+00,
	3.8912369355943404e+00, 5.2188432045968236e+00, 6.0710320291068589e+00
};
static const double phi_weight[NPHI] = {
	5.3823176663840211e-01, 1.1333659075855300e+00, 1.4699949793658615e+00,
	1.4699949793658615e+00, 1.1333659075855300e+00, 5.3823176663840211e-01
};

double compton_kernel_element(double x, double mu, double x1, double mu1,
                              double T_K, const int ktype)
{
	const double boltz_eV = 8.617333262e-5;   // Boltzmann constant [eV/K]
	const double mec2_eV  = 511.0e3;          // electron rest energy [eV]

	double x_inv = mec2_eV / (boltz_eV * T_K);   // m_e c² / kT

	double sin_mu  = sqrt(1.0 - mu  * mu);
	double sin_mu1 = sqrt(1.0 - mu1 * mu1);

	double intephi = 0.0;
	for (int np = 0; np < NPHI; ++np)
	{
		double costh = mu * mu1 + sin_mu * sin_mu1 * cos(phi_node[np]);
		intephi += profil_exact(x, x1, costh, x_inv) * phi_weight[np];
	}

	return intephi;
}

// ============================================================
// KernelCache implementation — banded storage with symmetry reduction
// + detailed balance (upper triangle only)
//
// Binary cache file format (version 04):
//   magic          [8 bytes]  "CKERN04\0"
//   NT,NE,NA_full,n_indep  [4×4 bytes]
//   data_size      [8 bytes]
//   T_grid         [NT doubles]
//   x_grid         [NE doubles]
//   theta          [NT doubles]
//   band_lo        [NT*NE*n_indep ints]
//   band_hi        [NT*NE*n_indep ints]
//   band_off       [NT*NE*n_indep longs]
//   glo            [NT*NE ints]
//   ghi            [NT*NE ints]
//   data           [data_size doubles]
// ============================================================
static const char CACHE_MAGIC[8] = "CKERN04";

// ------------------------------------------------------------
// build_canon — compute symmetry reduction tables from NA_full
//
// Symmetries:
//   (1) K(nm,nm1) = K(nm1,nm)              [angle transpose]
//   (2) K(nm,nm1) = K(NA-1-nm, NA-1-nm1)   [angle reflection]
//
// Combined group Z₂×Z₂ gives p(p+1) independent pairs,
// p = NA_full/2.
// ------------------------------------------------------------
void KernelCache::build_canon()
{
	int NA = NA_full;
	int max_indep = (NA / 2) * (NA / 2 + 1);

	canon     = new int[NA * NA];
	indep_nm  = new int[max_indep];
	indep_nm1 = new int[max_indep];

	for (int i = 0; i < NA * NA; ++i) canon[i] = -1;

	n_indep = 0;
	for (int nm = 0; nm < NA; ++nm)
	for (int nm1 = 0; nm1 < NA; ++nm1)
	{
		if (canon[nm * NA + nm1] >= 0) continue;   // already assigned

		int ia    = n_indep++;
		indep_nm[ia]  = nm;
		indep_nm1[ia] = nm1;

		int nm_r  = NA - 1 - nm;
		int nm1_r = NA - 1 - nm1;

		// Assign all four symmetry images (duplicates handled by check above)
		canon[nm    * NA + nm1 ] = ia;  // original
		canon[nm1   * NA + nm  ] = ia;  // transpose
		canon[nm_r  * NA + nm1_r] = ia; // reflection
		canon[nm1_r * NA + nm_r ] = ia; // transpose + reflection
	}

	fprintf(stdout, "  KernelCache: NA=%d  independent angle pairs = %d"
	        " (expected %d)\n", NA, n_indep, max_indep);
}

void KernelCache::write_header(FILE* fp) const
{
	long n_rows  = long(NT) * NE * n_indep;
	long n_ge    = long(NT) * NE;

	fwrite(CACHE_MAGIC, 1, 8, fp);
	fwrite(&NT,       sizeof(int),  1, fp);
	fwrite(&NE,       sizeof(int),  1, fp);
	fwrite(&NA_full,  sizeof(int),  1, fp);
	fwrite(&n_indep,  sizeof(int),  1, fp);
	fwrite(&data_size, sizeof(long), 1, fp);
	fwrite(T_grid,   sizeof(double), NT,     fp);
	fwrite(x_grid,   sizeof(double), NE,     fp);
	fwrite(theta,    sizeof(double), NT,     fp);
	fwrite(band_lo,  sizeof(int),    n_rows,  fp);
	fwrite(band_hi,  sizeof(int),    n_rows,  fp);
	fwrite(band_off, sizeof(long),   n_rows,  fp);
	fwrite(glo,      sizeof(int),    n_ge,    fp);
	fwrite(ghi,      sizeof(int),    n_ge,    fp);
}

void KernelCache::save(const char* filename) const
{
	KernelCacheOutput output(filename);
	FILE* fp = output.file();
	write_header(fp);
	fwrite(data,     sizeof(double), data_size, fp);
	output.commit();

	fprintf(stdout, "  KernelCache: saved to %s (%.1f MB)\n",
	        filename, data_size * 8.0 / (1024.0 * 1024.0));
}

bool KernelCache::load(const char* filename)
{
	// NT, NE, NA_full, n_indep, T_grid must be set before calling.
	FILE* fp = fopen(filename, "rb");
	if (!fp) return false;

	char magic[8];
	if (fread(magic, 1, 8, fp) != 8 || memcmp(magic, CACHE_MAGIC, 8) != 0)
	{
		fprintf(stderr,
			"  WARNING: cache file %s has wrong magic (expected \"%s\", got \"%.7s\")\n"
			"           File may be from an older format version. Will recompute.\n",
			filename, CACHE_MAGIC, magic);
		fclose(fp); return false;
	}

	int fNT, fNE, fNA, fNI;
	fread(&fNT, sizeof(int), 1, fp);
	fread(&fNE, sizeof(int), 1, fp);
	fread(&fNA, sizeof(int), 1, fp);
	fread(&fNI, sizeof(int), 1, fp);
	if (fNT != NT || fNE != NE || fNA != NA_full || fNI != n_indep)
	{
		fprintf(stderr,
			"  WARNING: cache file %s has mismatched dimensions:\n"
			"           file: NT=%d NE=%d NA=%d NI=%d\n"
			"           need: NT=%d NE=%d NA=%d NI=%d\n"
			"           Will recompute.\n",
			filename, fNT, fNE, fNA, fNI, NT, NE, NA_full, n_indep);
		fclose(fp); return false;
	}

	long f_data_size;
	fread(&f_data_size, sizeof(long), 1, fp);

	// Verify T_grid
	double* fT = new double[NT];
	fread(fT, sizeof(double), NT, fp);
	bool match = true;
	for (int i = 0; i < NT; ++i)
		if (fabs(fT[i] - T_grid[i]) > 1e-6 * fT[i]) { match = false; break; }
	delete[] fT;
	if (!match)
	{
		fprintf(stderr,
			"  WARNING: cache file %s has mismatched T_grid. Will recompute.\n",
			filename);
		fclose(fp); return false;
	}

	// Read x_grid and theta
	x_grid = new double[NE];
	fread(x_grid, sizeof(double), NE, fp);
	theta = new double[NT];
	fread(theta, sizeof(double), NT, fp);

	long n_rows = long(NT) * NE * n_indep;
	long n_ge   = long(NT) * NE;

	band_lo  = new int[n_rows];
	band_hi  = new int[n_rows];
	band_off = new long[n_rows];
	fread(band_lo,  sizeof(int),  n_rows, fp);
	fread(band_hi,  sizeof(int),  n_rows, fp);
	fread(band_off, sizeof(long), n_rows, fp);

	glo = new int[n_ge];
	ghi = new int[n_ge];
	fread(glo, sizeof(int), n_ge, fp);
	fread(ghi, sizeof(int), n_ge, fp);

	data_size = f_data_size;
	data = payload.map(fp,size_t(data_size));
	if(data) {
		fclose(fp);
		fprintf(stdout,"  KernelCache: mapped %s (%.1f MB virtual payload, demand-paged)\n",
		        filename,data_size*8.0/(1024.0*1024.0));
		return true;
	}
	data = new double[data_size];
	size_t nread = fread(data, sizeof(double), data_size, fp);
	fclose(fp);

	if ((long)nread != data_size) return false;

	double fill = 100.0 * data_size / (double(NT) * NE * n_indep * NE);
	fprintf(stdout, "  KernelCache: loaded from %s (%.1f MB, %.1f%% indep-dense fill)\n",
	        filename, data_size * 8.0 / (1024.0 * 1024.0), fill);
	return true;
}

void KernelCache::init(int n_ene, const double* ene_eV,
                       int n_ang, const double* mu, const double* wt,
                       const int ktype,
                       int nt_user, const double* T_user)
{
	NE      = n_ene;
	NA_full = n_ang;

	// Temperature grid: caller-supplied (e.g. single isothermal value in
	// test mode) or the default log-spaced N_T_CACHE grid for production.
	if (nt_user > 0 && T_user)
	{
		NT     = nt_user;
		T_grid = new double[NT];
		for (int i = 0; i < NT; ++i)
			T_grid[i] = T_user[i];
	}
	else
	{
		NT     = N_T_CACHE;
		T_grid = new double[NT];
		double log_lo = log10(T_CACHE_LO);
		double log_hi = log10(T_CACHE_HI);
		double dlog   = (NT > 1) ? (log_hi - log_lo) / (NT - 1) : 0.0;
		for (int i = 0; i < NT; ++i)
			T_grid[i] = pow(10.0, log_lo + i * dlog);
	}

	// Build symmetry tables (sets n_indep, canon, indep_nm/nm1)
	build_canon();

	// Try loading normalized cache first
	const char* cache_dir = getenv("COMPTON_CACHE_DIR");
	if (!cache_dir) cache_dir = ".";

	const char* ksuffix = (ktype == 0) ? "_ap" : "";
	char norm_file[512], raw_file[512];
	snprintf(norm_file, sizeof(norm_file),
	         "%s/kernel_norm_NE%d_NI%d_NT%d%s.bin",
	         cache_dir, NE, n_indep, NT, ksuffix);
	snprintf(raw_file, sizeof(raw_file),
	         "%s/kernel_NE%d_NI%d_NT%d%s.bin",
	         cache_dir, NE, n_indep, NT, ksuffix);

	if (load(norm_file))
	{
		fprintf(stderr, "\n Load normalized Compton scattering kernel");
		return;
	}
	if (load(raw_file))
	{
		fprintf(stderr,
			"\n  WARNING: Loaded un-normalized kernel %s\n"
			"           Normalized kernel %s not found.\n"
			"           Please run: ./normalize_kernel\n",
			raw_file, norm_file);
		return;
	}

	fprintf(stdout, "  KernelCache: no valid cache, computing "
	        "(banded, symmetry-reduced, upper-triangle)...\n");

	// Convert energy grid to dimensionless x = E/(m_e c²)
	const double mec2_eV  = 511.0e3;
	const double boltz_eV = 8.617333262e-5;
	x_grid = new double[NE];
	for (int ie = 0; ie < NE; ++ie)
		x_grid[ie] = ene_eV[ie] / mec2_eV;

	// Dimensionless temperature θ = kT/(m_e c²)
	theta = new double[NT];
	for (int i = 0; i < NT; ++i)
		theta[i] = boltz_eV * T_grid[i] / mec2_eV;

	// Allocate band structure over independent pairs
	long n_rows = long(NT) * NE * n_indep;
	band_lo  = new int[n_rows];
	band_hi  = new int[n_rows];
	band_off = new long[n_rows];

	KernelRowSpool computed_rows;

	const double KMIN_ABS = 1e-30;
	const double KMIN_REL = 1e-10;

	// --- Pass 1: find band limits (upper triangle only: ne1 >= ne) ---
	fprintf(stdout, "  KernelCache: pass 1 — finding band limits (upper triangle)...\n");

	computed_rows.evaluate(size_t(n_rows),size_t(NE),[&](size_t index,double* row_buf,int& lo_ne1,int& hi_ne1) {
		const long r=long(index);
		const int ia=int(index%size_t(n_indep));
		const int ne=int((index/size_t(n_indep))%size_t(NE));
		const int iT=int(index/(size_t(n_indep)*size_t(NE)));
		const int nm=indep_nm[ia], nm1=indep_nm1[ia];
		lo_ne1=NE; hi_ne1=-1;
			double K_max = 0.0;

			// Only scan ne1 >= ne (upper triangle)
			for (int ne1 = ne; ne1 < NE; ++ne1)
			{
				double val = compton_kernel_element(
					x_grid[ne], mu[nm], x_grid[ne1], mu[nm1], T_grid[iT], ktype);
				row_buf[ne1] = val;
				if (val > KMIN_ABS)
				{
					if (ne1 < lo_ne1) lo_ne1 = ne1;
					hi_ne1 = ne1;
				}
				if (val > K_max) K_max = val;
			}

			if (hi_ne1 < 0) { lo_ne1 = 0; hi_ne1 = -1; }

			// Trim band edges with dynamic threshold
			if (hi_ne1 >= 0)
			{
				double thresh = K_max * KMIN_REL;

				while (lo_ne1 <= hi_ne1)
				{
					if (row_buf[lo_ne1] > thresh) break;
					++lo_ne1;
				}
				while (hi_ne1 >= lo_ne1)
				{
					if (row_buf[hi_ne1] > thresh) break;
					--hi_ne1;
				}
				if (hi_ne1 < lo_ne1) { lo_ne1 = 0; hi_ne1 = -1; }
			}

			band_lo[r] = lo_ne1;
			band_hi[r] = hi_ne1;
	});

	// Compute per-(iT,ne) global band including both triangles.
	// Upper triangle: directly from stored bands.
	// Lower triangle: ne1 < ne where K(ne1, nm, ne, nm1) is stored,
	//   i.e., ne is within band of row (iT, ne1, ia).
	long n_ge = long(NT) * NE;
	glo = new int[n_ge];
	ghi = new int[n_ge];
	for (int iT = 0; iT < NT; ++iT)
	for (int ne = 0; ne < NE; ++ne)
	{
		int lo_g = NE, hi_g = -1;

		// Upper triangle (direct)
		for (int ia = 0; ia < n_indep; ++ia)
		{
			long r = row(iT, ne, ia);
			if (band_hi[r] >= band_lo[r])
			{
				if (band_lo[r] < lo_g) lo_g = band_lo[r];
				if (band_hi[r] > hi_g) hi_g = band_hi[r];
			}
		}

		// Lower triangle (via detailed balance): scan ne1 < ne downward
		for (int ne1 = ne - 1; ne1 >= 0; --ne1)
		{
			bool found = false;
			for (int ia = 0; ia < n_indep; ++ia)
			{
				long r = row(iT, ne1, ia);
				if (ne >= band_lo[r] && ne <= band_hi[r])
				{ found = true; break; }
			}
			if (found) {
				lo_g = ne1;
			} else {
				break;   // bands shrink away from diagonal
			}
		}

		glo[long(iT) * NE + ne] = (hi_g >= 0 || lo_g < NE) ? lo_g : 0;
		ghi[long(iT) * NE + ne] = hi_g;
	}

	// Compute offsets and total data size
	data_size = 0;
	for (long r = 0; r < n_rows; ++r)
	{
		band_off[r] = data_size;
		int bw = band_hi[r] - band_lo[r] + 1;
		if (bw > 0) data_size += long(bw);
	}

	double fill = 100.0 * data_size / (double(NT) * NE * n_indep * NE);
	fprintf(stdout, "  KernelCache: data_size = %ld doubles (%.1f MB, %.1f%% indep-dense fill)\n",
	        data_size, data_size * 8.0 / (1024.0 * 1024.0), fill);

	// Reuse the exact values evaluated during the band scan.
	fprintf(stdout, "  KernelCache: normalizing retained rows one temperature at a time...\n");
	KernelCacheOutput output(norm_file);
	write_header(output.file());
	const off_t payload_offset=::ftello(output.file());


	// --- Normalize via A23 sum rule (Poutanen & Svensson 1996, Eq. A23) ---
	// K() accessor handles both triangles (direct + detailed balance),
	// so the normalization sum covers the full energy range correctly.
	fprintf(stdout, "  Normalizing kernel via A23 sum rule...\n");

	for (int iT = 0; iT < NT; ++iT)
	{
		fprintf(stdout, "    T=%.2e K (%d/%d)...", T_grid[iT], iT + 1, NT);
		fflush(stdout);
		data=computed_rows.view(payload,size_t(data_size));

		double* norm = new double[NE];

		rt_parallel_depths(NE, [&](int ne) {
			double x = x_grid[ne];

			double integral = 0.0;
			for (int ne1 = 1; ne1 < NE; ++ne1)
			{
				double dx1 = x_grid[ne1] - x_grid[ne1 - 1];
				double f_a = 0.0, f_b = 0.0;
				// Reuse identical detailed-balance factors across angular pairs.
				// All products and the angular/energy summation order stay exact.
				const double balance_a=(ne<ne1-1)?exp(-(x_grid[ne1-1]-x_grid[ne])/theta[iT]):1.0;
				const double balance_b=(ne<ne1)?exp(-(x_grid[ne1]-x_grid[ne])/theta[iT]):1.0;
				auto lookup_a=[&](int nm,int nm1) {
					return (ne<ne1-1)?K_lower_with_balance(iT,ne1-1,nm,ne,nm1,balance_a)
					                 :K(iT,ne1-1,nm,ne,nm1);
				};
				auto lookup_b=[&](int nm,int nm1) {
					return (ne<ne1)?K_lower_with_balance(iT,ne1,nm,ne,nm1,balance_b)
					               :K(iT,ne1,nm,ne,nm1);
				};
				for (int inm = 0; inm < NA_full; ++inm)
				{
					if (mu[inm] <= 0.0) continue;
					int inm_neg = NA_full - 1 - inm;
					for (int inm1 = 0; inm1 < NA_full; ++inm1)
					{
						if (mu[inm1] <= 0.0) continue;
						double ww = wt[inm] * wt[inm1];
						double Ka = lookup_a(inm1,inm) + lookup_a(inm1,inm_neg);
						double Kb = lookup_b(inm1,inm) + lookup_b(inm1,inm_neg);
						f_a += ww * Ka;
						f_b += ww * Kb;
					}
				}
				f_a *= x_grid[ne1 - 1];
				f_b *= x_grid[ne1];
				integral += 0.5 * (f_a + f_b) * dx1;
			}

			double sigma_raw   = integral / x;
			double sigma_exact = compton_cross_section(ene_eV[ne], T_grid[iT])
			                   / phys::sigma_T;

			norm[ne] = (sigma_raw > 1e-30) ? sigma_exact / sigma_raw : 1.0;
		}, "DAO_KERNEL_THREADS", 8);

		// Apply: scale stored data by norm[ne_in]
		for (int ne_out = 0; ne_out < NE; ++ne_out)
		for (int ia = 0; ia < n_indep; ++ia)
		{
			long r = row(iT, ne_out, ia);
			int lo = band_lo[r];
			int hi = band_hi[r];
			if (hi < lo) continue;

			double* dst = data + band_off[r];
			for (int ne_in = lo; ne_in <= hi; ++ne_in)
				dst[ne_in - lo] *= norm[ne_in];
		}

		delete[] norm;
		const long begin=band_off[long(iT)*NE*n_indep];
		const long end=(iT+1<NT)?band_off[long(iT+1)*NE*n_indep]:data_size;
		if(end>begin && fwrite(data+begin,sizeof(double),size_t(end-begin),output.file())!=size_t(end-begin))
			throw std::runtime_error("kernel: normalized temperature write failed");
		fprintf(stdout, " done.\n");
	}

	// Save normalized kernel
	data=output.map_payload(payload,size_t(data_size),payload_offset);
	output.commit();
	fprintf(stdout, "  Saved normalized kernel to %s\n", norm_file);
}

void KernelCache::free_memory()
{
	delete[] T_grid;    T_grid    = nullptr;
	delete[] x_grid;    x_grid    = nullptr;
	delete[] theta;     theta     = nullptr;
	delete[] indep_nm;  indep_nm  = nullptr;
	delete[] indep_nm1; indep_nm1 = nullptr;
	delete[] canon;     canon     = nullptr;
	delete[] band_lo;   band_lo   = nullptr;
	delete[] band_hi;   band_hi   = nullptr;
	delete[] band_off;  band_off  = nullptr;
	delete[] glo;       glo       = nullptr;
	delete[] ghi;       ghi       = nullptr;
	payload.release(data);
	NT = NE = NA_full = n_indep = 0;
	data_size = 0;
}

int KernelCache::find_T(double T_K) const
{
	double logT = log10(T_K);
	int lo = 0, hi = NT - 1;
	while (lo < hi - 1)
	{
		int mid = (lo + hi) / 2;
		if (log10(T_grid[mid]) < logT) lo = mid;
		else hi = mid;
	}
	if (fabs(log10(T_grid[lo]) - logT) <= fabs(log10(T_grid[hi]) - logT))
		return lo;
	return hi;
}
