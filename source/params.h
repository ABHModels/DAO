#ifndef PARAMS_H
#define PARAMS_H

// Model parameters
struct ModelParams
{
	// --- Slab physical parameters ---
	double nh;         // log hydrogen column density [cm^-3]
	double zeta;       // ionization parameter exponent: xi = 10^zeta [erg cm s^-1]
	double frac;       // flux ratio: F_corona / F_disk
	double incidence;  // cosine of incidence angle (snapped to nearest GL node)
	int    i_incidence; // index into mu_gl[] for incidence angle

	// --- Corona model selection and parameters ---
	char   corona[32]; // corona model name: powerlaw, cutoffpl, nthcomp, comptt, blackbody
	double Gamma;      // photon index (powerlaw, cutoffpl, nthcomp)
	double E_cut;      // high-energy cutoff [keV] (cutoffpl)
	double E_lo_cut;   // low-energy exp cutoff [eV] (powerlaw, cutoffpl) — fixed, not user input
	double kT_e;       // electron temperature [keV] (nthcomp, comptt)
	double kT_bb;      // seed photon temperature [keV] (nthcomp, comptt, blackbody)
	double taup;       // plasma optical depth (comptt)

	// --- Disk parameters ---
	double kT_disk;    // disk blackbody temperature [eV]

	// --- RT solver ---
	double E_rt_lo;    // RT energy range lower bound [eV]
	double E_rt_hi;    // RT energy range upper bound [eV]
	int    maxiter;    // maximum Lambda iterations

	// --- Test mode ---
	bool   test_rt;    // if true, skip Cloudy, use synthetic atmosphere
	char   test_mode[16]; // required: "compps" (compPS benchmark slab) or "test_avg"
	double tau_slab;   // total vertical Thomson optical depth (compps test mode)

	// --- Abundances ---
	double Afe;        // iron abundance: log10(Fe/Fe_solar), 0 = solar

	// --- Angle dependent or not ---
	bool angsca;

	// --- Run management (auto-computed from params) ---
	char run_hash[12];   // 8-char hex hash of all physics params
	char run_dir[256];   // "results/<hash>/"

	// --- kernel type ---
	int ktype; // 0: approximation QED kernel 1: Exactly
};

// Set defaults and parse command-line overrides
ModelParams read_params(int argc, char *argv[]);

struct RTGrids;  // forward declaration

// Snap incidence to nearest GL node; call after grids.init_all()
void snap_incidence(ModelParams& par, const RTGrids& g);

#endif // PARAMS_H
