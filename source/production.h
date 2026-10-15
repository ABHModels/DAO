#ifndef PRODUCTION_H
#define PRODUCTION_H

#include "rt_grids.h"
#include "params.h"
#include "radiation.h"
#include "compton_kernel.h"
#include "avg_compton_kernel.h"
// Run the full Cloudy + RT outer iteration loop.
// Templated on the kernel-cache type (KernelCache or avgKernelCache);
// explicitly instantiated for both in production.cpp.
template<class Cache>
void run_production(RadField& rad, const RTGrids& g,
                    const ModelParams& par,
                    Cache& kcache);

#endif // PRODUCTION_H
