#include "corona_models.h"
#include "rt_grids.h"
#include "constants.h"
#include <cmath>
#include <cstring>
#include <cstdio>
#include <cstdlib>

// Xspec Fortran interface (nthcomp)
extern "C" void donthcomp_(float* ear, int* ne, float* param,
                           int* ifl, float* photar, float* photer);

// Xspec C interface (compTT)
extern "C" void C_compTT(const double* energy, int nFlux, const double* params,
                         int spectrumNumber, double* flux, double* fluxError,
                         const char* initStr);

// ============================================================
// Spectrum shape functions
// ============================================================

void powerlaw(double* spec, const RTGrids& g, const ModelParams& par)
{	
	const double E_lowcut = par.E_lo_cut * 1.0e3;
	for (int i = 0; i < g.NE; ++i)
	{
		double E_erg = g.ene[i] * phys::eV_to_erg;
		spec[i] = pow(E_erg, 1.0 - par.Gamma)
		         * exp(-E_lowcut / g.ene[i]);
	}
}

void cutoffpl(double* spec, const RTGrids& g, const ModelParams& par)
{
	const double E_hi = par.E_cut * 1.0e3;  // keV -> eV
	const double E_lowcut = par.E_lo_cut * 1.0e3;

	for (int i = 0; i < g.NE; ++i)
	{
		double E_erg = g.ene[i] * phys::eV_to_erg;
		spec[i] = pow(E_erg, 1.0 - par.Gamma)
		         * exp(-E_lowcut/ g.ene[i])
		         * exp(-g.ene[i] / E_hi);
	}
}

void comptt(double* spec, const RTGrids& g, const ModelParams& par)
{
	int ne = g.NE;
	double* ear      = new double[ne + 1];
	double* flux     = new double[ne];
	double* fluxerr  = new double[ne];
	double param[5];

	for (int i = 0; i < ne; ++i)
		ear[i] = (g.ene[i] - 0.5 * g.wid[i]) * 1.0e-3;
	ear[ne] = (g.ene[ne - 1] + 0.5 * g.wid[ne - 1]) * 1.0e-3;

	param[0] = 0.0;
	param[1] = par.kT_bb;
	param[2] = par.kT_e;
	param[3] = par.taup;
	param[4] = 1.0;

	C_compTT(ear, ne, param, 1, flux, fluxerr, "");

	for (int i = 0; i < ne; ++i)
	{
		double dE_eV = g.wid[i];
		double E_erg = g.ene[i] * phys::eV_to_erg;
		spec[i] = flux[i] * E_erg / dE_eV;
	}

	delete[] ear;
	delete[] flux;
	delete[] fluxerr;
}

void blackbody(double* spec, const RTGrids& g, double kT_eV)
{
	const double coeff = 2.0 / (phys::h * phys::h * phys::h * phys::c * phys::c);
	for (int i = 0; i < g.NE; ++i)
	{
		double E_erg = g.ene[i] * phys::eV_to_erg;
		double x = g.ene[i] / kT_eV;
		double denom = (x < 500.0) ? exp(x) - 1.0 : exp(x);
		spec[i] = M_PI * coeff * E_erg * E_erg * E_erg / denom;
	}
}

static void corona_blackbody(double* spec, const RTGrids& g, const ModelParams& par)
{
	blackbody(spec, g, par.kT_bb * 1.0e3);
}

void nthcomp(double* spec, const RTGrids& g, const ModelParams& par)
{
	int ne = g.NE;
	int ifl = 1;
	float* ear    = new float[ne + 1];
	float* photar = new float[ne];
	float* photer = new float[ne];
	float param[5];

	for (int i = 0; i < ne; ++i)
		ear[i] = (float)((g.ene[i] - 0.5 * g.wid[i]) * 1.0e-3);
	ear[ne] = (float)((g.ene[ne - 1] + 0.5 * g.wid[ne - 1]) * 1.0e-3);

	param[0] = (float)par.Gamma;
	param[1] = (float)par.kT_e;
	param[2] = (float)par.kT_bb;
	param[3] = 1.0f;
	param[4] = 0.0f;

	donthcomp_(ear, &ne, param, &ifl, photar, photer);

	for (int i = 0; i < ne; ++i)
	{
		double dE_eV = g.wid[i];
		double E_erg = g.ene[i] * phys::eV_to_erg;
		spec[i] = photar[i] * E_erg / dE_eV;
	}

	delete[] ear;
	delete[] photar;
	delete[] photer;
}

// ============================================================
// Dispatch map
// ============================================================
#include <map>
#include <string>

typedef void (*CoronaFunc)(double*, const RTGrids&, const ModelParams&);

static const std::map<std::string, CoronaFunc> corona_map = {
	{ "powerlaw",  powerlaw },
	{ "cutoffpl",  cutoffpl },
	{ "nthcomp",   nthcomp  },
	{ "comptt",    comptt   },
	{ "blackbody", corona_blackbody },
};

void compute_corona_shape(double* I_corona, const RTGrids& g, const ModelParams& par)
{
	auto it = corona_map.find(par.corona);
	if (it == corona_map.end())
	{
		fprintf(stderr, "Error: unknown corona model '%s'\n", par.corona);
		exit(1);
	}
	it->second(I_corona, g, par);
}
