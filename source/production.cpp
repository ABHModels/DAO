#include "production.h"
#include "thermal_production.h"

// The production model has one temperature/transfer iteration.
template<class Cache>
void run_production(RadField& rad, const RTGrids& g,
                    const ModelParams& par, Cache& kcache)
{
    run_thermal_production(rad, g, par, kcache);
}

// Explicit instantiations for both kernel-cache types.
template void run_production<KernelCache>(
	RadField&, const RTGrids&, const ModelParams&, KernelCache&);
template void run_production<avgKernelCache>(
	RadField&, const RTGrids&, const ModelParams&, avgKernelCache&);
