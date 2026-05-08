#ifndef COMPTON_KERNEL_POLAR_H
#define COMPTON_KERNEL_POLAR_H

#include "constants.h"
#include <cmath>

// ============================================================
// Polarized Compton redistribution kernel
//
// Computes all five Stokes components (R, RI, RQ, RU, RV) of the
// Compton scattering matrix for a thermal electron gas.
//
// Single-electron cross section: Nagirner & Poutanen (1993)
// Thermal averaging: 32-point Gauss-Laguerre quadrature
// Azimuth integration: 6-point Gauss-Legendre quadrature
//
// Reference:
//   Nagirner D.I., Poutanen J., 1993, A&A, 275, 325
//   Madej J., Rozanska A., Majczyna A., Nalezyty M., 2017, MNRAS, 469, 2032
//
// Symmetries exploited (same as unpolarized kernel):
//   (1) K(x,mu;x1,mu1)  = K(x,mu1;x1,mu)     [angle transpose]
//   (2) K(x,mu;x1,mu1)  = K(x,-mu;x1,-mu1)    [angle reflection]
//   (3) K(x,mu;x1,mu1)  = K(x1,mu;x,mu1) * exp(-(x-x1)/theta)
//                                                [detailed balance]
// ============================================================

static const int NCOMP = 5;   // R, RI, RQ, RU, RV

struct StokesKernel {
	double R;    // total redistribution
	double RI;   // intensity-intensity
	double RQ;   // Q-Stokes
	double RU;   // U-Stokes
	double RV;   // V-Stokes
};

// Redistribution function for all 5 Stokes components
// eps, eps1: photon energies in units of m_e c^2
// costh:     cosine of scattering angle
// x_inv:     m_e c^2 / (kT)
StokesKernel profil_polar(double eps, double eps1, double costh, double x_inv);

// Azimuth-integrated kernel element for all 5 components
StokesKernel compton_kernel_polar_element(double x, double mu, double x1, double mu1,
                                          double T_K, int ktype);

// ============================================================
// PolarKernelCache — banded storage with angular symmetry + detailed balance
//
// Same structure as KernelCache but stores 5 interleaved components
// per energy entry: [R, RI, RQ, RU, RV].
// ============================================================
struct PolarKernelCache
{
	int NT;       // number of temperature grid points
	int NE;       // number of energy bins
	int NA_full;  // total angle points
	int n_indep;  // independent (nm,nm1) pairs = (NA_full/2)*(NA_full/2+1)

	double* T_grid;   // [NT] temperatures [K]
	double* x_grid;   // [NE] dimensionless energy E/(m_e c^2)
	double* theta;    // [NT] dimensionless temperature kT/(m_e c^2)

	// Angle-pair symmetry tables
	int* indep_nm;   // [n_indep]
	int* indep_nm1;  // [n_indep]
	int* canon;      // [NA_full * NA_full]  (nm,nm1) -> ia

	// Band structure: [NT * NE * n_indep]
	// Stores only upper triangle: ne1 >= ne
	int*    band_lo;
	int*    band_hi;
	long*   band_off;   // offset in entry units (each entry = 5 doubles)

	// Global band per (iT, ne): union over all ia, both triangles
	int* glo;  // [NT * NE]
	int* ghi;  // [NT * NE]

	double* data;       // interleaved [R, RI, RQ, RU, RV] per entry
	long    data_size;  // total number of doubles stored

	PolarKernelCache() : NT(0), NE(0), NA_full(0), n_indep(0),
	                     T_grid(nullptr), x_grid(nullptr), theta(nullptr),
	                     indep_nm(nullptr), indep_nm1(nullptr), canon(nullptr),
	                     band_lo(nullptr), band_hi(nullptr), band_off(nullptr),
	                     glo(nullptr), ghi(nullptr),
	                     data(nullptr), data_size(0) {}

	void init(int n_ene, const double* ene_eV,
	          int n_ang, const double* mu, const double* wt,
	          int ktype);

	void save(const char* filename) const;
	bool load(const char* filename);
	void free_memory();

	int find_T(double T_K) const;

	// Row index for (iT, ne, ia)
	inline long row(int iT, int ne, int ia) const
	{
		return (long(iT) * NE + ne) * n_indep + ia;
	}

	// Retrieve all 5 components for K(iT, ne, nm, ne1, nm1)
	inline StokesKernel K(int iT, int ne, int nm, int ne1, int nm1) const
	{
		int ia = canon[nm * NA_full + nm1];
		StokesKernel sk = {0.0, 0.0, 0.0, 0.0, 0.0};

		if (ne1 >= ne)
		{
			long r = row(iT, ne, ia);
			if (ne1 < band_lo[r] || ne1 > band_hi[r]) return sk;
			long idx = NCOMP * (band_off[r] + (ne1 - band_lo[r]));
			sk.R  = data[idx + 0];
			sk.RI = data[idx + 1];
			sk.RQ = data[idx + 2];
			sk.RU = data[idx + 3];
			sk.RV = data[idx + 4];
		}
		else
		{
			// Lower triangle: detailed balance
			long r = row(iT, ne1, ia);
			if (ne < band_lo[r] || ne > band_hi[r]) return sk;
			long idx = NCOMP * (band_off[r] + (ne - band_lo[r]));
			double db = exp(-(x_grid[ne] - x_grid[ne1]) / theta[iT]);
			sk.R  = data[idx + 0] * db;
			sk.RI = data[idx + 1] * db;
			sk.RQ = data[idx + 2] * db;
			sk.RU = data[idx + 3] * db;
			sk.RV = data[idx + 4] * db;
		}
		return sk;
	}

	// Retrieve single R component (for normalization)
	inline double K_R(int iT, int ne, int nm, int ne1, int nm1) const
	{
		int ia = canon[nm * NA_full + nm1];
		if (ne1 >= ne)
		{
			long r = row(iT, ne, ia);
			if (ne1 < band_lo[r] || ne1 > band_hi[r]) return 0.0;
			return data[NCOMP * (band_off[r] + (ne1 - band_lo[r]))];
		}
		else
		{
			long r = row(iT, ne1, ia);
			if (ne < band_lo[r] || ne > band_hi[r]) return 0.0;
			return data[NCOMP * (band_off[r] + (ne - band_lo[r]))]
			     * exp(-(x_grid[ne] - x_grid[ne1]) / theta[iT]);
		}
	}

	// Global ne1 band for (iT, ne) — covers both triangles
	inline int lo(int iT, int ne) const { return glo[long(iT) * NE + ne]; }
	inline int hi(int iT, int ne) const { return ghi[long(iT) * NE + ne]; }

private:
	void build_canon();
};

#endif // COMPTON_KERNEL_POLAR_H
