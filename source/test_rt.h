#ifndef TEST_RT_H
#define TEST_RT_H

#include "radiation.h"
#include "params.h"
#include "compton_kernel.h"
#include "compton_cross_section.h"

void run_test_rt(RadField& rad, const RTGrids& g, ModelParams& par,
                 KernelCache& kcache);

#endif // TEST_RT_H
