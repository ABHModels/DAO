#ifndef CLOUDY_INTERFACE_H
#define CLOUDY_INTERFACE_H

#include "rt_grids.h"
#include "params.h"
#include "radiation.h"

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
};

// Bootstrap Cloudy to extract the energy grid
void bootstrap_cloudy_energy_grid(RTGrids& g, const char* save_file);

void extract_cloudy_output(int id, RadField& rad, const RTGrids& g, int outer_iter);
void save_cloudy_opacity(int id, const RadField& rad, const RTGrids& g,const ModelParams& par, int outer_iter);
void save_line_labels(int id, const RTGrids& g,const ModelParams& par, int outer_iter);

#endif // CLOUDY_INTERFACE_H
