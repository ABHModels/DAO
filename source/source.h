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
#include "avg_compton_kernel.h"
// Source function templated on the kernel-cache type. Works with both
// KernelCache (angle-dependent) and avgKernelCache (angle-mean): both
// expose find_T/lo/hi and a 5-arg K(iT, ne, nm, ne1, nm1).
// Explicitly instantiated for both types in source.cpp.
template<class Cache>
void compute_source_function(
	int ND, int NM, int NE,
	const double* intensity,
	const double* x_grid,
	const double* wmu,
	const Cache& kcache,
	const double* T_K,
	const double* const* jnu,
	const double* const* kabs,
	const double* const* ksct,
	const double* n_e,
	double* source);

// ============================================================
// Angle-mean source function (avgKernelCache).
//
// Takes the mean intensity J = (1/2)∫ I dμ (meani, shape [ND][NE])
// as input.  The Compton term has no angular integral — only the
// energy integration, which is identical to compute_source_function:
//
//   S_compton(x) = x² ∫ K̄(x; x₁) J(x₁) / x₁² dx₁
//
// The resulting source function is the same for every angle, so the
// scalar value is written to all NM angle slots of source[ND][NM][NE].
// ============================================================
void avgcompute_source_function(
	int ND, int NM, int NE,
	const double* meani,        // [ND][NE] mean intensity J
	const double* x_grid,
	const avgKernelCache& kcache,
	const double* T_K,
	const double* const* jnu,
	const double* const* kabs,
	const double* const* ksct,
	const double* n_e,
	double* source);
#endif // SOURCE_H