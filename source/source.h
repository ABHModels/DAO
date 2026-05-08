#ifndef SOURCE_H
#define SOURCE_H

// ============================================================
// Total source function for the RT equation:
//
//   S = jnu/(kabs+ksct) + ksct/(kabs+ksct) * S_compton
//
// Angle-dependent (KernelCache):
//   S_compton(x,μ) = x² ∫∫ K(x,μ;x₁,μ₁) I(x₁,μ₁)/x₁² dμ₁ dx₁
//
// All flat arrays are [ND][NM][NE] or [ND][NE].
// ============================================================

#include "compton_kernel.h"
// Angle-dependent source function (uses KernelCache v1)
void compute_source_function(
	int ND, int NM, int NE,
	const double* intensity,
	const double* x_grid,
	const double* wmu,
	const KernelCache& kcache,
	const double* T_K,
	const double* const* jnu,
	const double* const* kabs,
	const double* const* ksct,
	const double* n_e,
	double* source);
#endif // SOURCE_H