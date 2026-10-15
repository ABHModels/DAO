#ifndef DAO_THERMAL_RESTART_H
#define DAO_THERMAL_RESTART_H

#include "radiation.h"

// Legacy outputs retain T and angular moments, but not the full I(E,mu) or
// Cloudy atomic state. This is a warm continuation: reload T/J, rebuild atomic
// coefficients at fixed saved T, and solve directional transfer before roots.
// These functions only read files and reject incompatible/incomplete states.
void validate_thermal_resume(const ModelParams&);
void load_thermal_resume(RadField&, const RTGrids&, const ModelParams&);

#endif
