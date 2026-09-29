#ifndef DAO_CELL_TRANSFER_H
#define DAO_CELL_TRANSFER_H

#include "incidence_boundary.h"
#include <cmath>
#include <stdexcept>

namespace dao_rt {

// Exact formal solution for a cell with constant extinction and source.
// Keeping both small complementary weights avoids cancellation in thin cells.
struct CellWeights {
	double attenuation;
	double emission;
	double mean_incoming;
	double mean_source;
};

inline CellWeights cell_weights(double optical_depth)
{
	if (!std::isfinite(optical_depth) || optical_depth < 0.0)
		throw std::runtime_error("RT: invalid cell optical depth");
	const double t = optical_depth;
	if (t == 0.0) return {1.0, 0.0, 1.0, 0.0};
	const double emission = -std::expm1(-t);
	double mean_source, mean_incoming;
	if (t < 1e-3) {
		mean_source = t*(0.5 + t*(-1.0/6.0 + t*(1.0/24.0 +
		              t*(-1.0/120.0 + t/720.0))));
		mean_incoming = 1.0 - mean_source;
	} else {
		mean_incoming = emission/t;
		mean_source = 1.0 - mean_incoming;
	}
	return {std::exp(-t), emission, mean_incoming, mean_source};
}

// Flat cell arrays: [depth][angle][energy]. Surface arrays: [angle][energy].
// Each cell is traversed in full, including the first and last cells. The
// returned volume-averaged intensity is the quantity entering local absorption
// and scattering. Surface intensities are stored separately, including their
// prescribed incoming hemispheres. No source is interpolated across a face.
inline void formal_solution_cells(
	int ND, int NM, int NE, int i_inc, const double* mu,
	const double* ill_top, const double* ill_bot,
	const double* source, const CellWeights* weights,
	double* cell_intensity, double* surface_top, double* surface_bottom)
{
	for (int nm = 0; nm < NM; ++nm)
	for (int ne = 0; ne < NE; ++ne) {
		const long ray = long(nm)*NE + ne;
		const bool downward = mu[nm] < 0.0;
		double incoming = downward ?
			incident_top_intensity(i_inc, nm, ill_top[ne]) : ill_bot[ne];
		if (downward) surface_top[ray] = incoming;
		else surface_bottom[ray] = incoming;
		for (int step = 0; step < ND; ++step) {
			const int nd = downward ? step : ND-1-step;
			const long k = (long(nd)*NM + nm)*NE + ne;
			const CellWeights& w = weights[k];
			cell_intensity[k] = w.mean_incoming*incoming + w.mean_source*source[k];
			incoming = w.attenuation*incoming + w.emission*source[k];
		}
		if (downward) surface_bottom[ray] = incoming;
		else surface_top[ray] = incoming;
	}
}

} // namespace dao_rt
#endif
