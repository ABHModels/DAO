#include "radiation.h"
#include "corona_models.h"
#include "constants.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <algorithm>

// ============================================================
// IllumSpec
// ============================================================
void IllumSpec::allocate()
{
	I_corona = new double[g.NE]();
	I_disk   = new double[g.NE]();
}

void IllumSpec::deallocate()
{
	delete[] I_corona;  I_corona = nullptr;
	delete[] I_disk;    I_disk   = nullptr;
}

// ============================================================
// RadField allocation helpers
// ============================================================
static double** alloc_2d(int n1, int n2)
{
	double** p = new double*[n1];
	for (int i = 0; i < n1; ++i)
		p[i] = new double[n2]();
	return p;
}


static void free_2d(double**& p, int n1)
{
	if (!p) return;
	for (int i = 0; i < n1; ++i)
		delete[] p[i];
	delete[] p;
	p = nullptr;
}

static double*** alloc_3d(int n1, int n2, int n3)
{
	double*** p = new double**[n1];
	for (int i = 0; i < n1; ++i)
	{
		p[i] = new double*[n2];
		for (int j = 0; j < n2; ++j)
			p[i][j] = new double[n3]();
	}
	return p;
}

static void free_3d(double***& p, int n1, int n2)
{
	if (!p) return;
	for (int i = 0; i < n1; ++i)
	{
		for (int j = 0; j < n2; ++j)
			delete[] p[i][j];
		delete[] p[i];
	}
	delete[] p;
	p = nullptr;
}

void RadField::allocate()
{
	jnu       = alloc_2d(g.ND_MID, g.NE);
	kabs      = alloc_2d(g.ND_MID, g.NE);
	ksct      = alloc_2d(g.ND_MID, g.NE);
	Inu       = alloc_3d(g.ND_MID, g.NA, g.NE);
	J0        = alloc_2d(g.ND_MID, g.NE);
	J2        = alloc_2d(g.ND_MID, g.NE);
	J3        = alloc_2d(g.ND_MID, g.NE);
	jnu_line  = alloc_2d(g.ND_MID, g.NE);
	kabs_line = alloc_2d(g.ND_MID, g.NE);

	// Per-zone scalar arrays (heap-allocated, zero-initialised)
	T_K      = new double[g.ND_MID]();
	log_xi   = new double[g.ND_MID]();
	log_inte = new double[g.ND_MID]();
	n_e      = new double[g.ND_MID]();
	heating  = new double[g.ND_MID]();
	cooling  = new double[g.ND_MID]();

	illum.allocate();
}

void RadField::deallocate()
{
	free_2d(jnu,       g.ND_MID);
	free_2d(kabs,      g.ND_MID);
	free_2d(ksct,      g.ND_MID);
	free_3d(Inu,       g.ND_MID, g.NA);
	free_2d(J0,        g.ND_MID);
	free_2d(J2,        g.ND_MID);
	free_2d(J3,        g.ND_MID);
	free_2d(jnu_line,  g.ND_MID);
	free_2d(kabs_line, g.ND_MID);

	delete[] T_K;      T_K      = nullptr;
	delete[] log_xi;   log_xi   = nullptr;
	delete[] log_inte; log_inte = nullptr;
	delete[] n_e;      n_e      = nullptr;
	delete[] heating;  heating  = nullptr;
	delete[] cooling;  cooling  = nullptr;

	illum.deallocate();
}

// ============================================================
// IllumSpec::compute
// ============================================================
void IllumSpec::compute(const ModelParams& par)
{
	// I_corona, I_disk: spectral mean intensity [erg cm^-2 s^-1 eV^-1]
	compute_corona_shape(I_corona, g, par);
	blackbody(I_disk, g, par.kT_disk * 1.0e3);

	// Integrate spectral flux over energy to get total flux [erg cm^-2 s^-1]
	double raw_corona = 0.0, raw_disk = 0.0;
	for (int i = 0; i < g.NE - 1; ++i)
	{
		double dE = g.ene[i + 1] - g.ene[i];
		if (g.ene[i] > g.E_IN_LO && g.ene[i] < g.E_IN_HI) {
			raw_corona += 0.5 * (I_corona[i] + I_corona[i + 1]) * dE;
			raw_disk   += 0.5 * (I_disk[i]   + I_disk[i + 1])   * dE;
		}
	}

	// Target flux from ionisation parameter
	// xi = (4pi)^2 J / nh
	// J = 1/2 (I)
	// This function return J; 
	double xi = pow(10.0, par.zeta);
	double nH = pow(10.0, par.nh);
	double Fx = xi * nH / pow(phys::four_pi,2);
	// to match the xillver (*0.5 for elminate the 2 factor in compton_rt)
	// double Fx = xi * nH / (2*M_PI*0.7071);


	// Rescale so that corona + disk = Fx, split by frac = F_corona / F_disk
	if (par.frac > 0) {
		double target_corona = par.frac / (1.0 + par.frac) * Fx;
		double target_disk   = 1.0      / (1.0 + par.frac) * Fx;
		double scale_corona = target_corona / raw_corona;
		double scale_disk   = target_disk   / raw_disk;
		for (int i = 0; i < g.NE; ++i)
		{
			I_corona[i] *= scale_corona;
			I_disk[i]   *= scale_disk;
		}
	} else {
		double scale_corona = Fx / raw_corona;
		for (int i = 0; i < g.NE; ++i)
		{
			I_corona[i] *= scale_corona;
			I_disk[i]    = 0.0;
		}
	}

	printf("IllumSpec: Fx=%.4e  raw_corona=%.4e  raw_disk=%.4e\n",
		Fx, raw_corona, raw_disk);
}

// ============================================================
// Compute angular moments via GL quadrature
// ============================================================
void RadField::compute_moments()
{
	for (int id = 0; id < g.ND_MID; ++id)
	for (int ie = 0; ie < g.NE; ++ie)
	{
		double s0 = 0.0, s2 = 0.0, s3 = 0.0;
		for (int ia = 0; ia < g.NA; ++ia)
		{
			double mu = g.mu[ia];
			double w  = g.wt[ia];
			double I  = Inu[id][ia][ie];
			s0 += w * I;
			s2 += w * mu * I;
			s3 += w * mu * mu * I;
		}
		J0[id][ie] = 0.5*s0;
		J2[id][ie] = 0.5*s2;
		J3[id][ie] = 0.5*s3;
	}
}

// ============================================================
// Compute ionization parameter at each depth
// ============================================================
void RadField::compute_ionization_parameter(double lognh)
{	
	double nh=pow(10,lognh);
	for (int id = 0; id < g.ND_MID; ++id)
	{			
		// Now get the mean photons energy <E>

		double Ftot = 0.0;
		for (int ie = 0; ie < g.NE - 1; ++ie)
		{
			if (g.ene[ie] > g.E_IN_LO && g.ene[ie] < g.E_IN_HI){
				double dE = g.ene[ie + 1] - g.ene[ie];
				Ftot += 0.5 * (J0[id][ie] + J0[id][ie+1])*dE;
			}
		}
		double xi = pow((4.0 * M_PI),2) * Ftot / nh; 
		log_inte[id] = log10(xi*nh/phys::four_pi);
		log_xi[id] = log10(xi);
	}
}

void RadField::check_convergence(int outer_iter,
                                 double* T_old, double* xi_old,
                                 double& max_dT, double& max_dXi)
{
	// Force both to 1.0 on the first outer iteration so the
	// outer loop runs at least twice (the feedback is not active
	// until iteration 2).
	if (outer_iter == 1) 
	{
		max_dT = 1.0; max_dXi = 1.0; 
	} else{
		max_dT  = 0.0;
		max_dXi = 0.0;
		for (int id = 0; id < g.ND_MID; ++id)
		{
			max_dT = std::max(max_dT, fabs(log10(T_old[id]) - log10(T_K[id])));
			max_dXi = std::max(max_dXi, fabs(xi_old[id] - log_xi[id]));
			T_old[id]  = T_K[id];
			xi_old[id] = log_xi[id];
		}
	}
}
