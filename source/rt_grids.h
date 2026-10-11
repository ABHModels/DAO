#ifndef RT_GRIDS_H
#define RT_GRIDS_H

// ============================================================
// RTGrids — angle, depth, and energy grids for the RT solver
//
// Depth-grid size is runtime.  ND_EDGE and ND_MID are non-
// static int members; tau_edge, tau_mid and dr are heap-
// allocated double arrays whose size equals ND_EDGE / ND_MID.
// ============================================================

// Optional build-time resolution for reference benchmarks. All builds use
// the same transport implementation; production defaults to eight angles.
#ifndef DAO_RT_ANGLES
#define DAO_RT_ANGLES 8
#endif

struct RTGrids {

	// --- Compile-time constants ---
	static const int NA              = DAO_RT_ANGLES;
	static_assert(NA >= 2 && NA % 2 == 0, "RT requires an even angular grid");
	static const int ND_EDGE_DEFAULT = 49;
	static const int ND_MID_DEFAULT  = ND_EDGE_DEFAULT - 1;
	static const int NE_DEFAULT      = 1000;

	// --- Range of intensity ---
	static constexpr double E_IN_LO = 1;
	static constexpr double E_IN_HI = 1000e3;

	// --- Default grid bounds ---
	static constexpr double TAU_MIN = 1e-4;   // reference Thomson optical depth
	static constexpr double TAU_MAX = 5.0;
	static constexpr double E_LO    = 1;         // eV
	static constexpr double E_HI    = 1000e3;    // eV

	// --- Runtime sizes ---
	int NE;
	int ND_EDGE;   // number of depth edges
	int ND_MID;    // number of depth cells = ND_EDGE - 1

	// --- Angle grid (Gauss-Legendre on [-1, 1]) ---
	double mu[NA];
	double wt[NA];

	// --- Depth grid (heap-allocated, runtime-sized) ---
	double* tau_edge;    // size ND_EDGE
	double* tau_mid;     // size ND_MID
	double* dr;          // size ND_MID, physical thickness [cm]

	// --- Energy grid (eV, dynamically allocated) ---
	double* ene;         // bin centres [eV]
	double* wid;         // bin widths  [eV]

	RTGrids();
	~RTGrids();

	void init_angle();

	// Double-Gauss angle grid: NA/2-point Gauss-Legendre applied to each
	// hemisphere separately (nodes on [0,1] mirrored to [-1,0]). Reproduces
	// Xspec compPS's angular quadrature exactly when NA=10 (5 nodes/hemisphere).
	// Requires NA even.
	void init_angle_double_gauss();

	// Fixed reference depth: d tau_ref = 1.21*nH*sigma_T*dr.
	// Matches the RT scattering prescription in the Thomson limit. The free-electron
	// Thomson depth can differ as Cloudy n_e varies with ionization.
	// Positive edges from tau_min to tau_max, plus the surface edge at zero.
	// log_depth=true: log-uniform; false: tanh in log(tau), k=3.
	// Re-sizes the depth arrays to the default ND_EDGE_DEFAULT.
	void init_depth(double tau_min, double tau_max, double nh, bool log_depth = true);

	// Adopt Cloudy's adaptive depth grid.  tau_mids_thomson[i]
	// is the Thomson optical depth to the midpoint of zone i
	// (i in 0..n_zones-1), dr_cm[i] is its physical thickness.
	// Re-allocates tau_edge/tau_mid/dr to fit n_zones cells.
	void init_depth_from_zones(int n_zones,
	                           const double* tau_mids_thomson,
	                           const double* dr_cm);

	// Log-spaced energy grid (default NE_DEFAULT bins)
	void init_energy(double E_lo, double E_hi, int n_bins = NE_DEFAULT);

	// Energy grid from Cloudy's frequency mesh (bins in [E_LO, E_HI])
	void init_energy_from_cloudy(int nflux, const double* anu_ryd,
	                             const double* widflx_ryd);

	void save_energy(const char* path) const;
	bool load_energy(const char* path);

	void init_all(double nh, bool log_depth = true);

private:
	void free_depth();
	void alloc_depth(int n_edge);
};

#endif // RT_GRIDS_H
