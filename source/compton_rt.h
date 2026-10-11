#ifndef COMPTON_RT_H
#define COMPTON_RT_H

#include "radiation.h"
#include "params.h"
#include "thermal_balance.h"

// The single full-slab RT solver for production and prescribed-temperature
// benchmarks. Exact constant-cell transport returns volume-averaged intensities
// and true slab-face intensities. The scattering source is iterated until the
// maximum relative J (angle-mean) or I (directional) change is <1e-7 three times.
// This solves radiation at fixed material properties; it does not update T.
void compton_rt_solve(RadField&, const RTGrids&, const ModelParams&,
                     const std::vector<dao_thermal::Scattering>&,
                     std::vector<double>* cell_incoming=nullptr);

// Build the same photon-conserving scattering operators for both paths and
// copy their opacities into RadField. Density is always rad.n_e (Cloudy in
// production, prescribed free-electron density in synthetic benchmarks).
std::vector<dao_thermal::Scattering> make_scattering_column(
    RadField&, const RTGrids&, const avgKernelCache&);
std::vector<dao_thermal::Scattering> make_scattering_column(
    RadField&, const RTGrids&, const KernelCache&);

#endif // COMPTON_RT_H
