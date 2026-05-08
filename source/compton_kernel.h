#ifndef COMPTON_KERNEL_H
#define COMPTON_KERNEL_H

#include "constants.h"
#include <cmath>

// ============================================================
// Compton redistribution kernel K(x, μ; x₁, μ₁, T)
//
// The azimuth-integrated kernel of Nagirner & Poutanen (1993),
// using the exact from Madej et al.(2017).
// The redistribution function profil() was provided
// by Prof. Jerzy Madej (University of Warsaw).
//
// Reference:
//   Madej J., Różańska A., Majczyna A., Należyty M., 2017,
//   MNRAS, 469, 2032
//
// Usage in the radiative transfer:
//   S(x,μ) = x² ∫∫ K(x,μ; x₁,μ₁) I(x₁,μ₁) / x₁² dμ₁ dx₁
//
// Symmetries exploited to reduce storage:
//   (1) K(x,μ;x₁,μ₁)  = K(x,μ₁;x₁,μ)    [angle transpose]
//   (2) K(x,μ;x₁,μ₁)  = K(x,−μ;x₁,−μ₁)  [angle reflection]
//   (3) K(x,μ;x₁,μ₁)  = K(x₁,μ;x,μ₁) × exp(−(x−x₁)/Θ)
//                                            [detailed balance]
//
// (1)+(2): p(p+1) independent angle pairs, p = NA/2.
// (3):     only upper triangle (ne1 ≥ ne) stored.
//          Lower triangle derived via detailed balance.
//
// The kernel is stored in banded format: for each (iT, ne, ia),
// only the non-zero range of ne1 ≥ ne is stored.
// ============================================================

// Redistribution function R(ε, ε₁, cosθ; 1/Θ)
double profil_exact(double eps, double eps1, double costh, double x_inv);
double profil_exact_ap(double eps, double eps1, double costh, double x_inv);

// Azimuth-integrated Compton kernel element
double compton_kernel_element(double x, double mu, double x1, double mu1,
                              double T_K, const int ktype);

// ============================================================
// KernelCache — banded storage with angular symmetry + detailed balance
// ============================================================
struct KernelCache
{
	int NT;       // number of temperature grid points
	int NE;       // number of energy bins
	int NA_full;  // total angle points (full GL grid)
	int n_indep;  // independent (nm,nm1) pairs = (NA_full/2)*(NA_full/2+1)

	double* T_grid;   // [NT] temperatures [K]
	double* x_grid;   // [NE] dimensionless energy E/(m_e c²)
	double* theta;    // [NT] dimensionless temperature kT/(m_e c²)

	// Angle-pair symmetry tables
	int* indep_nm;   // [n_indep]
	int* indep_nm1;  // [n_indep]
	int* canon;      // [NA_full * NA_full]  (nm,nm1) → ia

	// Band structure: [NT * NE * n_indep]
	// Stores only upper triangle: ne1 ≥ ne for each row (iT, ne, ia)
	int*    band_lo;
	int*    band_hi;
	long*   band_off;

	// Global band per (iT, ne): union over all ia, both triangles
	int* glo;  // [NT * NE]
	int* ghi;  // [NT * NE]

	double* data;
	long    data_size;

	KernelCache() : NT(0), NE(0), NA_full(0), n_indep(0),
	                T_grid(nullptr), x_grid(nullptr), theta(nullptr),
	                indep_nm(nullptr), indep_nm1(nullptr), canon(nullptr),
	                band_lo(nullptr), band_hi(nullptr), band_off(nullptr),
	                glo(nullptr), ghi(nullptr),
	                data(nullptr), data_size(0) {}

	void init(int n_ene, const double* ene_eV,
	          int n_ang, const double* mu, const double* wt,
	          const int ktype);

	void save(const char* filename) const;
	bool load(const char* filename);
	void free_memory();

	int find_T(double T_K) const;

	// Row index for (iT, ne, ia)
	inline long row(int iT, int ne, int ia) const
	{
		return (long(iT) * NE + ne) * n_indep + ia;
	}

	// K(iT, ne, nm, ne1, nm1):
	//   ne1 >= ne → direct lookup (upper triangle)
	//   ne1 <  ne → detailed balance from stored K(ne1, nm, ne, nm1)
	inline double K(int iT, int ne, int nm, int ne1, int nm1) const
	{
		int ia = canon[nm * NA_full + nm1];

		if (ne1 >= ne)
		{
			// Upper triangle: direct lookup
			long r = row(iT, ne, ia);
			if (ne1 < band_lo[r] || ne1 > band_hi[r]) return 0.0;
			return data[band_off[r] + (ne1 - band_lo[r])];
		}
		else
		{
			// Lower triangle: detailed balance
			// K(ne, nm, ne1, nm1) = K(ne1, nm, ne, nm1) × exp(-(x[ne]-x[ne1])/θ)
			long r = row(iT, ne1, ia);
			if (ne < band_lo[r] || ne > band_hi[r]) return 0.0;
			double K_upper = data[band_off[r] + (ne - band_lo[r])];
			return K_upper * exp(-(x_grid[ne] - x_grid[ne1]) / theta[iT]);
		}
	}

	// Global ne1 band for (iT, ne) — covers both triangles.
	inline int lo(int iT, int ne) const { return glo[long(iT) * NE + ne]; }
	inline int hi(int iT, int ne) const { return ghi[long(iT) * NE + ne]; }

private:
	void build_canon();
};

#endif // COMPTON_KERNEL_H
