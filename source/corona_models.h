#ifndef CORONA_MODELS_H
#define CORONA_MODELS_H

#include "params.h"

struct RTGrids;  // forward declaration

// ============================================================
// Spectrum shape functions
// Each fills spec[NE] with an unnormalised shape on the RT grid.
// ============================================================

void powerlaw(double* spec, const RTGrids& g, const ModelParams& par);
void cutoffpl(double* spec, const RTGrids& g, const ModelParams& par);
void nthcomp(double* spec, const RTGrids& g, const ModelParams& par);
void comptt(double* spec, const RTGrids& g, const ModelParams& par);
void blackbody(double* spec, const RTGrids& g, double kT_eV);

// Dispatch: compute corona shape based on par.corona string.
void compute_corona_shape(double* spec, const RTGrids& g, const ModelParams& par);

#endif // CORONA_MODELS_H
