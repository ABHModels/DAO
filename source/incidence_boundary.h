#ifndef INCIDENCE_BOUNDARY_H
#define INCIDENCE_BOUNDARY_H

#include "radiation.h"
#include "constants.h"
#include <cmath>

// Fixed reflionx flux normalization: I_corona is F_E and ill_top is specific
// intensity I_E [erg cm^-2 s^-1 eV^-1 sr^-1]. It enforces
// F_E = 2*pi*sum_{mu<0}(wt*abs(mu)*I_E), including on an isotropic grid.
// i_inc == -1 denotes isotropic incidence over the downward hemisphere.
inline void compute_boundary_illumination(
	int NE, int i_inc, const RTGrids& g, const RadField& rad,
	double* ill_top, double* ill_bot)
{
	double incident_weight = 0.0;
	if (i_inc == -1)
	{
		for (int nm = 0; nm < g.NA; ++nm)
			if (g.mu[nm] < 0.0)
				incident_weight += g.wt[nm] * std::abs(g.mu[nm]);
	}
	else
	{
		incident_weight = g.wt[i_inc] * std::abs(g.mu[i_inc]);
	}

	for (int ne = 0; ne < NE; ++ne)
	{
		ill_top[ne] = rad.illum.I_corona[ne] / (0.5 * phys::four_pi * incident_weight);
		// I_disk is the incident mean intensity from the lower hemisphere.
		ill_bot[ne] = 2.0 * rad.illum.I_disk[ne];
	}
}

inline double incident_top_intensity(int i_inc, int nm, double ill_top)
{
	return (i_inc == -1 || nm == i_inc) ? ill_top : 0.0;
}

#endif // INCIDENCE_BOUNDARY_H
