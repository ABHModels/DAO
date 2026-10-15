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
// copy their opacities into RadField. Production treats 1.21*nH scattering
// targets as free electrons; synthetic benchmarks retain their prescribed n_e.
double scattering_target_density(const RadField&, const ModelParams&, int depth);
std::vector<dao_thermal::Scattering> make_scattering_column(
    RadField&, const RTGrids&, const avgKernelCache&, const ModelParams&);
std::vector<dao_thermal::Scattering> make_scattering_column(
    RadField&, const RTGrids&, const KernelCache&, const ModelParams&);

#endif // COMPTON_RT_H
