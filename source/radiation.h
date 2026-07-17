#ifndef RADIATION_H
#define RADIATION_H

#include "rt_grids.h"
#include "params.h"

// ============================================================
// IllumSpec — incident spectra (corona + disk)
// ============================================================
struct IllumSpec
{
	const RTGrids& g;
	double* I_corona;
	double* I_disk;

	IllumSpec(const RTGrids& gr)
		: g(gr), I_corona(nullptr), I_disk(nullptr) {}

	void allocate();
	void deallocate();
	void compute(const ModelParams& par);
};

// ============================================================
// RadField — radiation quantities on the RT grid
//
// Must be created AFTER RTGrids is initialised.
// ============================================================
struct RadField
{
	const RTGrids& g;

	double** jnu;    // [ND_MID][NE]
	double** kabs;   // [ND_MID][NE]
	double** ksct;   // [ND_MID][NE]
	double*** Inu;   // [ND_MID][NA][NE]
	double** J0;     // [ND_MID][NE]
	double** J2;     // [ND_MID][NE]
	double** J3;     // [ND_MID][NE]
	double** jnu_line;  // [ND_MID][NM][NE] — line emissivity (unscaled)
	double** kabs_line; // [ND_MID][NE] — bin-averaged line opacity

	// Runtime-sized depth arrays (heap-allocated in allocate())
	double* T_K;
	double* log_xi;
	double* log_inte;
	double* n_e;
	double* heating;
	double* line_heat; // Pdest line photons deposited locally as heat [erg cm^-3 s^-1]
	double* cooling;

	IllumSpec illum;

	RadField(const RTGrids& gr)
		: g(gr),
		  jnu(nullptr), kabs(nullptr), ksct(nullptr),
		  Inu(nullptr), J0(nullptr), J2(nullptr), J3(nullptr),
		  jnu_line(nullptr), kabs_line(nullptr),
		  T_K(nullptr), log_xi(nullptr), log_inte(nullptr),
		  n_e(nullptr), heating(nullptr), line_heat(nullptr),
		  cooling(nullptr),
		  illum(gr) {}

	void allocate();
	void deallocate();

	void compute_moments();
	void compute_ionization_parameter(double nh);
	void check_convergence(int outer_iter,
	                       double* T_old, double* xi_old,
	                       double& max_dT, double& max_dXi);
};

#endif // RADIATION_H
