#ifndef PRODUCTION_H
#define PRODUCTION_H

#include "rt_grids.h"
#include "params.h"
#include "radiation.h"
#include "compton_kernel.h"
// Run the full Cloudy + RT outer iteration loop.
void run_production(RadField& rad, const RTGrids& g,
                    const ModelParams& par,
                    KernelCache& kcache);

#endif // PRODUCTION_H
