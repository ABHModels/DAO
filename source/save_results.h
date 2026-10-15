#ifndef SAVE_RESULTS_H
#define SAVE_RESULTS_H

#include "radiation.h"
#include "params.h"

void save_results(const RadField& rad, const RTGrids& g,
                  const ModelParams& par, int iter);

// Clear only our named final diagnostic files when reusing a run directory.
void clear_ion_fraction_output(const ModelParams& par);
// Called only after the outer solver has passed all convergence checks.
void save_ion_fractions(const RadField& rad, const RTGrids& g,
                        const ModelParams& par, int iter);

#endif // SAVE_RESULTS_H
