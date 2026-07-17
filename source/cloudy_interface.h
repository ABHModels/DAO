#ifndef CLOUDY_INTERFACE_H
#define CLOUDY_INTERFACE_H

#include "rt_grids.h"
#include "params.h"
#include "radiation.h"
#include <string>
#include <vector>

// One bound-bound transition in one depth cell.
// Cloudy supplies the local atomic rates; apply_line_escape() later uses the
// full column of records to estimate how many line photons actually escape.
struct LineRec
{
	int    bin;     // RT-grid bin index (0..NE-1)
	long   ip;      // Cloudy LineSave index (aligns the same line across cells)
	double emiss;   // optically thin line power, 4pi*j_line [erg cm^-3 s^-1]
	double dE_eV;   // continuum bin width used to convert line power to per-eV jnu
	double E_eV;    // line photon energy [eV]
	double wave;    // Cloudy line wavelength label value
	double kappaL;  // line-center opacity [cm^-1], including stimulated emission
	double kappaLo; // lower-level opacity [cm^-1], used only for electron escape
	double damp;    // Voigt damping parameter a
	double y;       // collisional quench probability relative to A_ul: C_ul/A_ul
	int    redis;   // redistribution type iRedisFun()
	std::string label;   // Cloudy line label, e.g. ion + wavelength
	std::string comment; // Cloudy spectroscopic comment when available
};

// ============================================================
// CloudyInput — manages Cloudy commands for each depth call
// ============================================================
struct CloudyInput
{
	const RTGrids& g;
	double nh;
	double Afe;

	CloudyInput(const RTGrids& gr) : g(gr), nh(0.0), Afe(1.0) {}

	void init(const ModelParams& par);
	void issue_constant();
	void issue_depth(int id, const RadField& rad,const RTGrids& g, const ModelParams& par);
	void issue_depth_lastest(int id, const RadField& rad, const RTGrids& g,const ModelParams& par);
};

// Bootstrap Cloudy to extract the energy grid
void bootstrap_cloudy_energy_grid(RTGrids& g, const char* save_file);

// Pull continuum and fluorescence into jnu/jnu_line, and collect bound-bound
// line records for column-scale escape processing.
void extract_cloudy_output(int id, RadField& rad, const RTGrids& g,
                           int outer_iter, std::vector<LineRec>& recs);

// After the full depth column: apply bound-bound escape to every stored line.
// A trapped photon can escape in the line, Thomson-scatter out of the core,
// collisionally quench, or be absorbed by continuum opacity.
//   P = (beta + Pelec)*(1+y)/(beta + Pelec + y + Pdest)
//   j_line = emiss * P / (dE * 4pi)
// Continuum-destroyed line power is stored in rad.line_heat for the next
// Cloudy thermal-balance pass.
void apply_line_escape(RadField& rad, const RTGrids& g,
                       const std::vector<std::vector<LineRec>>& store,
                       const ModelParams& par, int iter);


#endif // CLOUDY_INTERFACE_H
