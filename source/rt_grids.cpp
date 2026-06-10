#include "rt_grids.h"
#include "constants.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>

RTGrids::RTGrids()
	: NE(0),
	  ND_EDGE(ND_EDGE_DEFAULT),
	  ND_MID(ND_EDGE_DEFAULT - 1),
	  tau_edge(nullptr), tau_mid(nullptr), dr(nullptr),
	  ene(nullptr), wid(nullptr)
{
	alloc_depth(ND_EDGE);
}

RTGrids::~RTGrids()
{
	free_depth();
	delete[] ene;
	delete[] wid;
}

void RTGrids::alloc_depth(int n_edge)
{
	ND_EDGE  = n_edge;
	ND_MID   = n_edge - 1;
	tau_edge = new double[ND_EDGE]();
	tau_mid  = new double[ND_MID]();
	dr       = new double[ND_MID]();
}

void RTGrids::free_depth()
{
	delete[] tau_edge;  tau_edge = nullptr;
	delete[] tau_mid;   tau_mid  = nullptr;
	delete[] dr;        dr       = nullptr;
}

void RTGrids::init_all(double nh)
{
	init_angle();
	init_depth(TAU_MIN, TAU_MAX, nh);
	init_energy(E_LO, E_HI);
}

// ============================================================
// Angle grid: NA-point Gauss-Legendre on [-1, 1]
// ============================================================
void RTGrids::init_angle()
{
	const int n = NA;

	for (int i = 0; i < n; ++i)
	{
		double x = -cos(M_PI * (i + 0.75) / (n + 0.5));

		for (int iter = 0; iter < 100; ++iter)
		{
			double p0 = 1.0, p1 = x;
			for (int j = 2; j <= n; ++j)
			{
				double p2 = ((2*j - 1) * x * p1 - (j - 1) * p0) / j;
				p0 = p1; p1 = p2;
			}
			double dp = n * (x * p1 - p0) / (x * x - 1.0);
			double dx = -p1 / dp;
			x += dx;
			if (fabs(dx) < 1e-15) break;
		}

		mu[i] = x;

		double p0 = 1.0, p1 = x;
		for (int j = 2; j <= n; ++j)
		{
			double p2 = ((2*j - 1) * x * p1 - (j - 1) * p0) / j;
			p0 = p1; p1 = p2;
		}
		double dp = n * (x * p1 - p0) / (x * x - 1.0);
		wt[i] = 2.0 / ((1.0 - x * x) * dp * dp);
	}
}

// ============================================================
// Double-Gauss angle grid (per-hemisphere Gauss-Legendre)
//
// NA/2-point Gauss-Legendre on [-1,1] gives roots x_k, weights w_k.
// Each root is mapped onto the half-interval mu_h = (1+x_k)/2 in (0,1)
// with weight w_k/2, and placed symmetrically in both hemispheres:
//   downward  mu = -mu_h        upward  mu = +mu_h
// Each hemisphere's weights then sum to 1 (total 2), consistent with
// J = 0.5 * sum(wt*I).  With NA=10 (5 nodes/hemisphere) this is exactly
// the angle grid used by Xspec compPS (QDRGSDO mapped to [0,1]).
//
// Nodes are stored ascending in mu, mirroring init_angle().
// ============================================================
void RTGrids::init_angle_double_gauss()
{
	if (NA % 2 != 0)
	{
		fprintf(stderr,
		        "Error: init_angle_double_gauss requires even NA (got %d)\n",
		        NA);
		exit(1);
	}

	const int n = NA / 2;   // Gauss-Legendre nodes per hemisphere

	for (int i = 0; i < n; ++i)
	{
		double x = -cos(M_PI * (i + 0.75) / (n + 0.5));

		for (int iter = 0; iter < 100; ++iter)
		{
			double p0 = 1.0, p1 = x;
			for (int j = 2; j <= n; ++j)
			{
				double p2 = ((2*j - 1) * x * p1 - (j - 1) * p0) / j;
				p0 = p1; p1 = p2;
			}
			double dp = n * (x * p1 - p0) / (x * x - 1.0);
			double dx = -p1 / dp;
			x += dx;
			if (fabs(dx) < 1e-15) break;
		}

		double p0 = 1.0, p1 = x;
		for (int j = 2; j <= n; ++j)
		{
			double p2 = ((2*j - 1) * x * p1 - (j - 1) * p0) / j;
			p0 = p1; p1 = p2;
		}
		double dp = n * (x * p1 - p0) / (x * x - 1.0);
		double w  = 2.0 / ((1.0 - x * x) * dp * dp);

		// x ascending in [-1,1]  ->  mu_h ascending in (0,1)
		double mu_h = 0.5 * (1.0 + x);
		double wt_h = 0.5 * w;

		// Place so the full mu[] array stays ascending:
		//   index 0..n-1   : downward (-mu_h), largest |mu| first
		//   index n..NA-1  : upward   (+mu_h), ascending
		mu[n - 1 - i] = -mu_h;  wt[n - 1 - i] = wt_h;
		mu[n + i]     =  mu_h;  wt[n + i]     = wt_h;
	}
}

// ============================================================
// Depth grid
// ============================================================
void RTGrids::init_depth(double tau_min, double tau_max, double nh)
{
	// Ensure arrays have the default size (in case they were
	// previously resized by init_depth_from_zones()).
	free_depth();
	alloc_depth(ND_EDGE_DEFAULT);

	tau_edge[0] = 0.0;

	// --- Symmetric logarithmic depth grid (TITAN-style) ---
	// log tau(u) = log tau_mid + A * tanh(k * u), u in [-1/2, +1/2]
	// Refines zones near both tau_min and tau_max.
	// -------------------------------------------------------
	const double k        = 3.0;
	const double log_tmin = log10(tau_min);
	const double log_tmax = log10(tau_max);
	const double log_tmid = 0.5 * (log_tmin + log_tmax);
	const double A        = (log_tmax - log_tmid) / tanh(k / 2.0);

	const int N = ND_EDGE - 1;   // tau_edge[1]..tau_edge[N]
	for (int i = 1; i < ND_EDGE; ++i)
	{
		double u = static_cast<double>(i - 1) / (N - 1) - 0.5;
		double log_tau = log_tmid + A * tanh(k * u);
		tau_edge[i] = pow(10.0, log_tau);
	}

	for (int i = 0; i < ND_MID; ++i)
		tau_mid[i] = 0.5 * (tau_edge[i] + tau_edge[i + 1]);

	const double n_e = 1.2 * pow(10.0, nh);
	const double factor = 1.0 / (n_e * phys::sigma_T);

	for (int i = 0; i < ND_MID; ++i)
		dr[i] = (tau_edge[i + 1] - tau_edge[i]) * factor;
}

// ============================================================
// Adopt Cloudy's adaptive depth grid.
//
// Inputs (one element per Cloudy zone):
//   tau_mids_thomson[i] — Thomson optical depth to zone midpoint
//                         (i.e. struc.depth[i] * n_e * sigma_T
//                         if midpoints; or an equivalent measure)
//   dr_cm[i]           — physical thickness of that zone [cm]
//                         (= struc.drad[i] from Cloudy)
//
// Convention:
//   tau_edge[0]   = 0 (illuminated surface)
//   tau_edge[i+1] = tau_edge[i] + (n_e*sigma_T)*dr_cm[i]
// but since dr and tau here are pre-computed outside, we
// derive tau_edge from tau_mid and dr as follows:
//   tau_edge[0]   = max(0, tau_mid[0] - 0.5 * dtau_equiv[0])
// and enforce consecutive edges to be monotonic.
// ============================================================
void RTGrids::init_depth_from_zones(int n_zones,
                                    const double* tau_mids_thomson,
                                    const double* dr_cm)
{
	free_depth();
	alloc_depth(n_zones + 1);

	// Copy zone midpoints and thicknesses
	for (int i = 0; i < n_zones; ++i)
	{
		tau_mid[i] = tau_mids_thomson[i];
		dr[i]      = dr_cm[i];
	}

	// Reconstruct Thomson-tau edges by accumulating (tau_mid[i+1] - tau_mid[i-1])/2.
	// The illuminated surface is tau_edge[0] = 0.
	tau_edge[0] = 0.0;
	for (int i = 1; i < ND_EDGE; ++i)
	{
		if (i == ND_EDGE - 1)
		{
			// Far edge: extend by half the last mid-to-mid spacing
			double half = 0.5 * (tau_mid[i - 1] - (i >= 2 ? tau_mid[i - 2] : 0.0));
			tau_edge[i] = tau_mid[i - 1] + half;
		}
		else
		{
			tau_edge[i] = 0.5 * (tau_mid[i - 1] + tau_mid[i]);
		}
	}
}

// ============================================================
// Log-spaced energy grid
// ============================================================
void RTGrids::init_energy(double E_lo, double E_hi, int n_bins)
{
	delete[] ene;
	delete[] wid;

	NE  = n_bins;
	ene = new double[NE];
	wid = new double[NE];

	double log_lo = log10(E_lo);
	double log_hi = log10(E_hi);
	double dlog   = (log_hi - log_lo) / NE;

	for (int i = 0; i < NE; ++i)
	{
		double e_lo = pow(10.0, log_lo + i * dlog);
		double e_hi = pow(10.0, log_lo + (i + 1) * dlog);
		ene[i] = 0.5 * (e_lo + e_hi);
		wid[i] = e_hi - e_lo;
	}

	fprintf(stdout, "Energy grid (log-spaced): NE=%d  E=[%.3f, %.3f] eV\n",
	        NE, ene[0], ene[NE - 1]);
}

// ============================================================
// Energy grid from Cloudy's frequency mesh.
// Copies bins within [E_LO, E_HI], converts Ryd -> eV.
// ============================================================
void RTGrids::init_energy_from_cloudy(int nflux, const double* anu_ryd,
                                      const double* widflx_ryd)
{
	delete[] ene;
	delete[] wid;

	double E_lo_ryd = E_LO / phys::eV_per_Ryd;
	double E_hi_ryd = E_HI / phys::eV_per_Ryd;

	// Count bins in range
	int count = 0;
	for (int j = 0; j < nflux; ++j)
		if (anu_ryd[j] >= E_lo_ryd && anu_ryd[j] <= E_hi_ryd)
			++count;

	NE  = count;
	ene = new double[NE];
	wid = new double[NE];

	int k = 0;
	for (int j = 0; j < nflux; ++j)
	{
		if (anu_ryd[j] >= E_lo_ryd && anu_ryd[j] <= E_hi_ryd)
		{
			ene[k] = anu_ryd[j]    * phys::eV_per_Ryd;
			wid[k] = widflx_ryd[j] * phys::eV_per_Ryd;
			++k;
		}
	}

	fprintf(stdout, "Energy grid (from Cloudy): NE=%d  E=[%.3f, %.3f] eV\n",
	        NE, ene[0], ene[NE - 1]);
}

// ============================================================
// Save energy grid (ene + wid) to a text file
// ============================================================
void RTGrids::save_energy(const char* path) const
{
	FILE* fp = fopen(path, "w");
	if (!fp)
	{
		fprintf(stderr, "Error: cannot write energy grid to %s\n", path);
		exit(1);
	}
	fprintf(fp, "# Cloudy energy grid: NE=%d\n", NE);
	fprintf(fp, "# Col 1: ene [eV]  Col 2: wid [eV]\n");
	for (int i = 0; i < NE; ++i)
		fprintf(fp, "%.15e  %.15e\n", ene[i], wid[i]);
	fclose(fp);
	fprintf(stdout, "Energy grid saved to %s (%d bins)\n", path, NE);
}

// ============================================================
// Load energy grid from a text file.  Returns true on success.
// ============================================================
bool RTGrids::load_energy(const char* path)
{
	FILE* fp = fopen(path, "r");
	if (!fp) return false;

	// Count data lines (skip comments)
	int count = 0;
	char line[256];
	while (fgets(line, sizeof(line), fp))
		if (line[0] != '#' && line[0] != '\n') ++count;

	if (count == 0) { fclose(fp); return false; }

	delete[] ene;
	delete[] wid;
	NE  = count;
	ene = new double[NE];
	wid = new double[NE];

	rewind(fp);
	int k = 0;
	while (fgets(line, sizeof(line), fp))
	{
		if (line[0] == '#' || line[0] == '\n') continue;
		sscanf(line, "%lf %lf", &ene[k], &wid[k]);
		++k;
	}
	fclose(fp);

	fprintf(stdout, "Energy grid (from file): NE=%d  E=[%.3f, %.3f] eV  (%s)\n",
	        NE, ene[0], ene[NE - 1], path);
	return true;
}