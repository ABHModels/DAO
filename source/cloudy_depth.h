#ifndef DAO_CLOUDY_DEPTH_H
#define DAO_CLOUDY_DEPTH_H
#include "cloudy_interface.h"

// Run one independent atomic-physics calculation per cell, then return all
// results before column-wide line escape/RT. workers=0 reads DAO_CLOUDY_WORKERS
// (default min(4, hardware concurrency)); workers=1 uses the serial path.
// Call only from the main thread after all kernel/RT workers have joined.
void run_cloudy_depths(RadField& rad, const RTGrids& g, const ModelParams& par,
                      int iteration, std::vector<std::vector<LineRec>>& lines,
                      unsigned workers = 0);
#endif
