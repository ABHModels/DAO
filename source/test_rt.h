#ifndef TEST_RT_H
#define TEST_RT_H

#include "radiation.h"
#include "params.h"
#include "compton_kernel.h"
#include "avg_compton_kernel.h"
#include "compton_cross_section.h"

// Templated on the kernel-cache type (KernelCache or avgKernelCache);
// explicitly instantiated for both in test_rt.cpp.
template<class Cache>
void run_test_rt(RadField& rad, const RTGrids& g, ModelParams& par,
                 Cache& kcache);

#endif // TEST_RT_H
