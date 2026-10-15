#ifndef DAO_THERMAL_PRODUCTION_H
#define DAO_THERMAL_PRODUCTION_H
#include "cloudy_interface.h"
#include "compton_kernel.h"
#include "thermal_balance.h"

// Temperature iteration uses the shared RT solver in compton_rt.h.
// two branchs for angle-dep and angle-mean
// It's stupid way to do it, do we have a better way to do it?
void thermal_equilibrate_column(RadField&, const RTGrids&, const ModelParams&,
                                const avgKernelCache&,
                                std::vector<std::vector<LineRec>>&,
                                int iteration);
void thermal_equilibrate_column(RadField&, const RTGrids&, const ModelParams&,
                                const KernelCache&,
                                std::vector<std::vector<LineRec>>&,
                                int iteration);
void run_thermal_production(RadField&, const RTGrids&, const ModelParams&,
                            const avgKernelCache&, bool supplied_initial_state=false);
void run_thermal_production(RadField&, const RTGrids&, const ModelParams&,
                            const KernelCache&, bool supplied_initial_state=false);
#endif
