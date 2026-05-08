#ifndef COMPTON_RT_H
#define COMPTON_RT_H

// ============================================================
// Compton scattering RT solver
//
// Second-order short characteristics formal solution with
// Lambda iteration.
//
// References:
//   Hubeny I., Mihalas D., 2015, Theory of Stellar Atmospheres,
//     Princeton University Press (Section 12.4)
//   Suleimanov V., Poutanen J., Werner K., 2012, A&A, 545, A120
// ============================================================

#include "radiation.h"
#include "compton_kernel.h"
#include "params.h"

void compton_rt_solve(RadField& rad, const RTGrids& g,
                      const ModelParams& par,
                      const KernelCache& kcache, int maxiter);
#endif // COMPTON_RT_H
