#ifndef AVG_COMPTON_KERNEL_H
#define AVG_COMPTON_KERNEL_H

#include "constants.h"
#include "kernel_payload.h"
#include <cmath>

// ============================================================
// Angle-mean Compton redistribution kernel  K̄(x; x₁, T)
//
// The fully angle-averaged counterpart of KernelCache.  Where the
// directional kernel keeps a (μ, μ₁) dependence, here both the
// incoming and outgoing directions are integrated out, so the kernel
// is a function of the two photon energies and the temperature only.
//
// Starting from the redistribution function R(x, x₁, cosθ; 1/Θ)
// (profil_exact), the scattering is azimuthally symmetric about the
// incoming photon direction, so integrating over the full outgoing
// solid angle is just
//
//     K̄(x; x₁) = ∮ dΩ_out R(x, x₁, cosθ)
//              = 2π ∫₋₁¹ R(x, x₁, cosθ) d(cosθ)
//
// i.e. the integral over the scattering-angle cosine times 2π for the
// azimuth.  Used in the angle-averaged RT:
//
//     S̄(x) = x² ∫ K̄(x; x₁) J(x₁) / x₁² dx₁
//
// Storage reduction:
//   With no angular resolution, the only symmetry available is
//   detailed balance,
//
//     K̄(x; x₁) = K̄(x₁; x) × exp(−(x − x₁)/Θ)
//
//   so only the upper triangle (ne1 ≥ ne) is stored; the lower
//   triangle is reconstructed on the fly.  Each (iT, ne) row is kept
//   in banded form, holding only the non-zero range of ne1 ≥ ne.
//
// Reference:
//   Madej J., Różańska A., Majczyna A., Należyty M., 2017,
//   MNRAS, 469, 2032   (exact differential cross section / profil)
// ============================================================

// Angle-mean (azimuth + scattering-angle integrated) kernel element:
//   K̄(x; x₁, T) = 2π ∫₋₁¹ R(x, x₁, cosθ; 1/Θ) d(cosθ)
double angle_mean_kernel_element(double x, double x1, double T_K,
                                 const int ktype);

// ============================================================
// avgKernelCache — banded storage, detailed balance only
// ============================================================
struct avgKernelCache
{
	int NT;   // number of temperature grid points
	int NE;   // number of energy bins

	double* T_grid;   // [NT] temperatures [K]
	double* x_grid;   // [NE] dimensionless energy E/(m_e c²)
	double* theta;    // [NT] dimensionless temperature kT/(m_e c²)

	// Band structure: [NT * NE]
	// Stores only the upper triangle ne1 ≥ ne for each row (iT, ne).
	int*  band_lo;
	int*  band_hi;
	long* band_off;

	// Global band per (iT, ne): covers both triangles.
	int* glo;  // [NT * NE]
	int* ghi;  // [NT * NE]

	double* data;
	KernelPayload payload;
	long    data_size;

	avgKernelCache() : NT(0), NE(0),
	                   T_grid(nullptr), x_grid(nullptr), theta(nullptr),
	                   band_lo(nullptr), band_hi(nullptr), band_off(nullptr),
	                   glo(nullptr), ghi(nullptr),
	                   data(nullptr), data_size(0) {}

	// The angle arguments (n_ang, mu, wt) are accepted only to keep the
	// call signature identical to KernelCache::init — they are unused,
	// since the angle-mean kernel carries no angular resolution.
	//
	// nt_user > 0 with T_user != nullptr builds the kernel on the supplied
	// temperature grid (e.g. a single isothermal value in test mode)
	// instead of the default log-spaced N_T_CACHE grid.
	void init(int n_ene, const double* ene_eV,
	          int n_ang, const double* mu, const double* wt,
	          const int ktype,
	          int nt_user = 0, const double* T_user = nullptr);

	void save(const char* filename) const;
	void write_header(FILE* fp) const;
	bool load(const char* filename);
	void free_memory();

	int find_T(double T_K) const;

	// Row index for (iT, ne)
	inline long row(int iT, int ne) const
	{
		return long(iT) * NE + ne;
	}

	// 5-argument overload mirroring KernelCache::K(iT, ne, nm, ne1, nm1).
	// The angle-mean kernel carries no angular resolution, so the incoming
	// and outgoing angle indices are ignored and the isotropic value is
	// returned.  This lets angle-resolved callers (e.g. the existing source
	// function) bind to avgKernelCache unchanged.
	inline double K(int iT, int ne, int /*nm*/, int ne1, int /*nm1*/) const
	{
		return K(iT, ne, ne1);
	}

	// K̄(iT, ne, ne1):
	//   ne1 >= ne → direct lookup (upper triangle)
	//   ne1 <  ne → detailed balance from stored K̄(ne1, ne)
	inline double K(int iT, int ne, int ne1) const
	{
		if (ne1 >= ne)
		{
			// Upper triangle: direct lookup
			long r = row(iT, ne);
			if (ne1 < band_lo[r] || ne1 > band_hi[r]) return 0.0;
			return data[band_off[r] + (ne1 - band_lo[r])];
		}
		else
		{
			// Lower triangle: detailed balance
			// K̄(ne, ne1) = K̄(ne1, ne) × exp(-(x[ne]-x[ne1])/θ)
			long r = row(iT, ne1);
			if (ne < band_lo[r] || ne > band_hi[r]) return 0.0;
			double K_upper = data[band_off[r] + (ne - band_lo[r])];
			return K_upper * exp(-(x_grid[ne] - x_grid[ne1]) / theta[iT]);
		}
	}

	// Global ne1 band for (iT, ne) — covers both triangles.
	inline int lo(int iT, int ne) const { return glo[long(iT) * NE + ne]; }
	inline int hi(int iT, int ne) const { return ghi[long(iT) * NE + ne]; }
};

#endif // AVG_COMPTON_KERNEL_H
