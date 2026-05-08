#ifndef RT_GRIDS_H
#define RT_GRIDS_H

// ============================================================
// RTGrids — angle, depth, and energy grids for the RT solver
//
// Depth-grid size is runtime.  ND_EDGE and ND_MID are non-
// static int members; tau_edge, tau_mid and dr are heap-
// allocated double arrays whose size equals ND_EDGE / ND_MID.
// ============================================================

struct RTGrids {

	// --- Compile-time constants ---
	static const int NA              = 8;
	static const int ND_EDGE_DEFAULT = 51;
	static const int ND_MID_DEFAULT  = ND_EDGE_DEFAULT - 1;
	static const int NE_DEFAULT      = 1000;

	// --- Range of intensity ---
	static constexpr double E_IN_LO = 1;
	static constexpr double E_IN_HI = 1000e3;

	// --- Default grid bounds ---
	static constexpr double TAU_MIN = 1e-4;   // Thomson optical depth
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

	// Symmetric-log tanh grid; re-sizes the depth arrays to
	// the default ND_EDGE_DEFAULT.
	void init_depth(double tau_min, double tau_max, double nh);

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

	void init_all(double nh);

private:
	void free_depth();
	void alloc_depth(int n_edge);
};

#endif // RT_GRIDS_H
