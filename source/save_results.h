#ifndef SAVE_RESULTS_H
#define SAVE_RESULTS_H

#include "radiation.h"
#include "params.h"

void save_results(const RadField& rad, const RTGrids& g,
                  const ModelParams& par, int iter);

#endif // SAVE_RESULTS_H
