#ifndef RADIATION_H
#define RADIATION_H

#include "rt_grids.h"
#include "params.h"
#include <array>
#include <vector>

// Read-only snapshot of one Cloudy cell. Stage 0 is neutral; the last is bare.
// Fractions are relative to the gas-phase elemental abundance, as in Cloudy.
template<std::size_t Stages> struct IonFractions {
	int iteration = -1;
	double gas_density = 0;
	std::array<double,Stages> fraction{};
};

// ============================================================
// IllumSpec — incident spectra (corona + disk)
// ============================================================
struct IllumSpec
{
	const RTGrids& g;
	double* I_corona; // incident normal flux F_E [erg cm^-2 s^-1 eV^-1]
	double* I_disk;   // lower-hemisphere mean intensity J_E [erg cm^-2 s^-1 eV^-1 sr^-1]

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
	double*** Inu;   // [ND_MID][NA][NE], production: cell-volume-average intensity
	double** Inu_top;    // [NA][NE], actual upper-face intensity (both hemispheres)
	double** Inu_bottom; // [NA][NE], actual lower-face intensity (both hemispheres)
	double** J0;     // [ND_MID][NE], angular mean of Inu (same spatial convention)
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
	double* line_heat; // line-destruction diagnostic; never injected as external heat
	double* cooling;
	// Allocated only for requested diagnostics; carried through worker IPC.
	std::vector<IonFractions<9>> oxygen;
	std::vector<IonFractions<27>> iron;

	IllumSpec illum;

	RadField(const RTGrids& gr)
		: g(gr),
		  jnu(nullptr), kabs(nullptr), ksct(nullptr),
		  Inu(nullptr), Inu_top(nullptr), Inu_bottom(nullptr),
		  J0(nullptr), J2(nullptr), J3(nullptr),
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
