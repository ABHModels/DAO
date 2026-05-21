#include "source.h"
#include "rt_grids.h"
#include "constants.h"
#include <cstring>
#include <cmath>

// ============================================================
// Flat-array indexing: [nd][nm][ne] -> (nd*NM + nm)*NE + ne
// ============================================================
static int s_NM, s_NE;

static inline long sidx3(int nd, int nm, int ne)
{
	return (long(nd) * s_NM + nm) * s_NE + ne;
}

// ============================================================
// compute_source_function
//
// S(nd,nm,ne) = (jnu[nd][ne] + jnu_line[nd][nm][ne]) / ktot
//             + ksct[nd][ne] / ktot * S_compton(nd,nm,ne)
//
// where ktot = kabs + ksct, and
//   S_compton = x² ∫∫ K(x,μ;x₁,μ₁) I(x₁,μ₁)/x₁² dμ₁ dx₁
// ============================================================
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
	double* source)
{
	s_NM = NM;
	s_NE = NE;

	memset(source, 0, long(ND) * NM * NE * sizeof(double));

	for (int nd = 0; nd < ND; ++nd)
	{	
		int iT = 0;
		if (T_K[nd]>1e6){
			iT = kcache.find_T(T_K[nd]);
		}else{
			
			// these is to avoid the numerical instability for low temperature kernel
			// the problem can be fixed by increase the energy resolution
			double temT = 1e6;
			iT = kcache.find_T(temT);
		}

		for (int nm = 0; nm < NM; ++nm)
		for (int ne = 0; ne < NE; ++ne)
		{	
			double ktot = kabs[nd][ne] + ksct[nd][ne];

			// --- Thermal term: (jnu_cont + jnu_line) / ktot ---
			double S_th = jnu[nd][ne]/ ktot;

			// --- Compton scattering term ---
			double S_sct = 0.0;

			int ne1_lo = kcache.lo(iT, ne);
			int ne1_hi = kcache.hi(iT, ne);

			if ((ne1_hi >= ne1_lo))
			{
				double x = x_grid[ne];

				double trapz  = 0.0;
				double f_prev = 0.0;

				for (int ne1 = ne1_lo; ne1 <= ne1_hi; ++ne1)
				{
					double x1 = x_grid[ne1];
					double inte_mu = 0.0;
					for (int nm1 = 0; nm1 < NM; ++nm1)
						inte_mu += kcache.K(iT, ne, nm, ne1, nm1)
						         * intensity[sidx3(nd, nm1, ne1)]
						         * wmu[nm1];

					double f_curr = inte_mu / (x1 * x1);

					if (ne1 > ne1_lo)
						trapz += 0.5 * (f_curr + f_prev)
						       * (x_grid[ne1] - x_grid[ne1 - 1]);

					f_prev = f_curr;
				}

				S_sct = trapz * x * x;
			}
			source[sidx3(nd, nm, ne)] = S_th + n_e[nd]*phys::sigma_T / ktot * S_sct;
		}
	}
}

