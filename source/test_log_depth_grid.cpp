#include "rt_grids.h"
#include "constants.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>

int main()
{
	RTGrids g;
	g.init_depth(RTGrids::TAU_MIN, RTGrids::TAU_MAX, 15.0);
	assert(g.ND_MID == 100 && g.ND_EDGE == 101);
	assert(g.tau_edge[0] == 0.0);
	assert(g.tau_edge[1] == RTGrids::TAU_MIN);
	assert(g.tau_edge[100] == RTGrids::TAU_MAX);
	const double ratio = std::pow(RTGrids::TAU_MAX / RTGrids::TAU_MIN, 1.0 / 99.0);
	for (int i = 2; i <= 100; ++i) {
		assert(g.tau_edge[i] > g.tau_edge[i-1]);
		assert(std::abs(g.tau_edge[i] / g.tau_edge[i-1] - ratio) < 1e-12);
	}
	for (int i = 0; i < 100; ++i) {
		const double expected = (g.tau_edge[i+1]-g.tau_edge[i]) /
			(phys::reference_electrons_per_hydrogen * std::pow(10.0, 15.0) * phys::sigma_T);
		assert(std::abs(g.dr[i]/expected - 1.0) < 1e-12);
	}
	std::printf("PASS 100 log depth cells, monotonic ratio, and 1.21*nH distance scale\n");
}
