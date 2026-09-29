#include "compton_rt.h"
#include "bezier3_transfer.h"
#include "source.h"
#include "constants.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <vector>

static inline void source_dispatch(
	const KernelCache& kc, int ND, int NM, int NE,
	const double* intensity, const double* /*meani*/,
	const double* x_grid, const double* wmu, const double* T_K,
	const double* const* jnu, const double* const* kabs,
	const double* const* ksct, const double n_h, double* source)
{
	compute_source_function(ND, NM, NE, intensity, x_grid, wmu,
	                        kc, T_K, jnu, kabs, ksct, n_h, source);
}

static inline void source_dispatch(
	const avgKernelCache& kc, int ND, int NM, int NE,
	const double* /*intensity*/, const double* meani,
	const double* x_grid, const double* /*wmu*/, const double* T_K,
	const double* const* jnu, const double* const* kabs,
	const double* const* ksct, const double n_h, double* source)
{
	avgcompute_source_function(ND, NM, NE, meani, x_grid,
	                           kc, T_K, jnu, kabs, ksct, n_h, source);
}

template<class Cache>
void compton_rt_solve(RadField& rad, const RTGrids& g,
                      const ModelParams& par,
                      const Cache& kcache, int maxiter)
{
	const int ND = g.ND_MID, NM = g.NA, NE = g.NE;
	if (ND < 1 || NE < 1 || maxiter < 1)
		throw std::runtime_error("RT: empty grid or invalid iteration limit");
	if (kcache.NE != NE || kcache.NT < 1)
		throw std::runtime_error("RT: kernel cache does not match the energy grid");
	if (par.i_incidence < -1 || par.i_incidence >= NM ||
	    (par.i_incidence >= 0 && g.mu[par.i_incidence] >= 0.0))
		throw std::runtime_error("RT: invalid incident angle");

	// Flat arrys for the solver: [ND][NM][NE] and [ND][NE] are flattened to 1D.
	const long size3 = long(ND)*NM*NE, size2 = long(ND)*NE;
	std::vector<double> x(NE), top(NE), bottom(NE);
	std::vector<double> source(size3), intensity(size3);
	std::vector<double> mean(size2), previous_mean(size2);
	std::vector<double> face_top(long(NM)*NE), face_bottom(long(NM)*NE);
	std::vector<double> optical_depth(size3);

	// x = E/m_ec^2
	for (int e = 0; e < NE; ++e) x[e] = g.ene[e]/phys::m_e_c2_eV;

	// Compute the boundary illumination for the top and bottom faces of the slab.
	compute_boundary_illumination(NE, par.i_incidence, g, rad,
	                              top.data(), bottom.data());

	// Check for nonfinite or negative boundary illumination.
	for (int e = 0; e < NE; ++e)
		if (!std::isfinite(top[e]) || top[e] < 0.0 ||
		    !std::isfinite(bottom[e]) || bottom[e] < 0.0)
			throw std::runtime_error("RT: invalid boundary illumination");

	// Store each full physical cell depth. The Bezier solver combines adjacent
	// half cells for center-to-center paths and retains both boundary half cells.
	for (int d = 0; d < ND; ++d) {
		if (!std::isfinite(rad.T_K[d]) || rad.T_K[d] <= 0.0)
			throw std::runtime_error("RT: invalid cell temperature");
		if (!std::isfinite(g.dr[d]) || g.dr[d] <= 0.0)
			throw std::runtime_error("RT: invalid cell width");
		for (int e = 0; e < NE; ++e) {
			const double absorption = rad.kabs[d][e], scattering = rad.ksct[d][e];
			const double extinction = absorption + scattering;
			if (!std::isfinite(extinction) || extinction <= 0.0 ||
			    absorption < 0.0 || scattering < 0.0 ||
			    !std::isfinite(rad.jnu[d][e]) || rad.jnu[d][e] < 0.0)
				throw std::runtime_error("RT: invalid extinction or emissivity");
			for (int m = 0; m < NM; ++m) {
				const long k = (long(d)*NM+m)*NE+e;
				optical_depth[k] = extinction*g.dr[d]/std::abs(g.mu[m]);
				source[k] = rad.jnu[d][e]/extinction;
			}
		}
	}

	dao_rt::Bezier3FormalSolver formal_solver(ND,NM,NE,optical_depth.data());
	std::printf("\n=== Compton RT solver (Bezier3 short characteristics) ===\n");
	std::printf("  ND=%d  NM=%d  NE=%d  maxiter=%d\n", ND, NM, NE, maxiter);

	// Memory usage estimate (MB) for the solver's working arrays.
	const double mem_mb = ((3.0*size3 + 2.0*size2 + (3.0+2.0*NM)*NE)*
	                      sizeof(double)+formal_solver.memory_bytes())/(1024.0*1024.0);
	std::printf("  Work arrays: %.1f MB\n", mem_mb);

	// Solve the RT equation iteratively, updating the source function
	const auto start = std::chrono::steady_clock::now();
	double formal_seconds = 0.0, source_seconds = 0.0;
	double max_change = std::numeric_limits<double>::infinity();
	int successive = 0, iterations = 0;
	for (; iterations < maxiter;) {
		const auto step_start = std::chrono::steady_clock::now();

		// Formal solution: compute center intensities and actual surface output.
		formal_solver.solve(par.i_incidence, g.mu, top.data(), bottom.data(),
			source.data(), intensity.data(), face_top.data(), face_bottom.data());

		// convergence check and the mean intensity in each cell
		max_change = 0.0;
		for (int d = 0; d < ND; ++d)
		for (int e = 0; e < NE; ++e) {
			double J = 0.0;
			for (int m = 0; m < NM; ++m) {
				const double I = intensity[(long(d)*NM+m)*NE+e];
				if (!std::isfinite(I) || I < 0.0)
					throw std::runtime_error("RT: nonfinite or negative intensity");
				J += 0.5*g.wt[m]*I;
			}
			const long k = long(d)*NE+e;
			if (!std::isfinite(J)) throw std::runtime_error("RT: nonfinite mean intensity");
			const double scale = std::max(J, previous_mean[k]);
			if (scale > 0.0)
				max_change = std::max(max_change, std::abs(J-previous_mean[k])/scale);
			mean[k] = J;
		}
		const auto source_start = std::chrono::steady_clock::now();
		formal_seconds += std::chrono::duration<double>(source_start-step_start).count();
		++iterations;
		successive = max_change < 1e-4 ? successive+1 : 0;

		// Update the source function for the next iteration, unless we have already converged
		if (successive < 3 && iterations < maxiter) {
			source_dispatch(kcache, ND, NM, NE, intensity.data(), mean.data(),
				x.data(), g.wt, rad.T_K, rad.jnu, rad.kabs, rad.ksct,
				std::pow(10.0, par.nh), source.data());
			for (double S : source)
				if (!std::isfinite(S) || S < 0.0)
					throw std::runtime_error("RT: nonfinite or negative source function");
		}
		source_seconds += std::chrono::duration<double>(
			std::chrono::steady_clock::now()-source_start).count();
		std::printf("    iter %3d  max|dJ|/max(J,Jold)=%.6e  converged=%d/3  time=%.3fs\n",
			iterations, max_change, successive,
			std::chrono::duration<double>(std::chrono::steady_clock::now()-step_start).count());
		if (successive == 3) break;
		previous_mean = mean;
	}
	std::printf("  RT solver total time: %.2f s\n",
		std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count());
	std::printf("  RT wall-time breakdown: formal+mean=%.3fs source=%.3fs\n",
		formal_seconds, source_seconds);
	if (successive != 3) {
		std::fprintf(stderr, "RT failed to converge after %d iterations (residual %.6e)\n",
		             iterations, max_change);
		throw std::runtime_error("RT iteration limit reached before convergence");
	}
	std::printf("  Converged at iteration %d.\n", iterations);

	// Copy the results back to the RadField structure
	for (int d = 0; d < ND; ++d) {
		for (int m = 0; m < NM; ++m)
		for (int e = 0; e < NE; ++e)
			rad.Inu[d][m][e] = intensity[(long(d)*NM+m)*NE+e];
		for (int e = 0; e < NE; ++e) rad.J0[d][e] = mean[long(d)*NE+e];
	}
	for (int m = 0; m < NM; ++m)
	for (int e = 0; e < NE; ++e) {
		rad.Inu_top[m][e] = face_top[long(m)*NE+e];
		rad.Inu_bottom[m][e] = face_bottom[long(m)*NE+e];
	}
	std::printf("=== RT solver complete ===\n\n");
}

template void compton_rt_solve<KernelCache>(
	RadField&, const RTGrids&, const ModelParams&, const KernelCache&, int);
template void compton_rt_solve<avgKernelCache>(
	RadField&, const RTGrids&, const ModelParams&, const avgKernelCache&, int);
