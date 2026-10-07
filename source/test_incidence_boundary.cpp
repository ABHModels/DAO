#include "incidence_boundary.h"
#include "constants.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <vector>

static void close_to(double actual, double expected, const char* label)
{
	if (!std::isfinite(actual) ||
	    std::abs(actual-expected) > 1e-11*std::max(1e-100,std::abs(expected)))
		throw std::runtime_error(label);
}

static double integrate_spectrum(const RTGrids& g, const double* values)
{
	double sum=0;
	for (int e=0; e<g.NE-1; ++e)
		if (g.ene[e]>g.E_IN_LO && g.ene[e]<g.E_IN_HI)
			sum+=0.5*(values[e]+values[e+1])*(g.ene[e+1]-g.ene[e]);
	return sum;
}

int main()
{
	try {
		RTGrids g;
		g.init_angle();
		g.init_energy(10.0,1e6,200);
		RadField rad(g);
		rad.allocate();
		ModelParams p{};
		std::strcpy(p.corona,"cutoffpl");
		p.Gamma=2.0;
		p.E_cut=300.0;
		p.E_lo_cut=0.1;
		p.kT_disk=0.35;
		p.nh=15.0;
		p.zeta=3.0;
		const double target_J=std::pow(10.0,p.zeta+p.nh)/
		                      (phys::four_pi*phys::four_pi);
		for (bool double_gauss : {false,true}) {
			if (double_gauss) g.init_angle_double_gauss();
			for (double frac : {-1.0,2.0}) {
				p.frac=frac;
				rad.illum.compute(p);
				const double corona_J=integrate_spectrum(g,rad.illum.I_corona);
				const double disk_J=integrate_spectrum(g,rad.illum.I_disk);
				close_to(corona_J+disk_J,target_J,"IllumSpec total J");
				close_to(corona_J,target_J*(frac>0 ? frac/(1+frac) : 1),
				         "IllumSpec corona J");
				close_to(disk_J,target_J*(frac>0 ? 1/(1+frac) : 0),
				         "IllumSpec disk J");
				for (int d=0; d<g.ND_MID; ++d)
					for (int e=0; e<g.NE; ++e)
						rad.J0[d][e]=rad.illum.I_corona[e]+rad.illum.I_disk[e];
				rad.compute_ionization_parameter(p.nh);
				close_to(rad.log_xi[0],p.zeta,"recovered ionization parameter");

				for (double incidence : {0.7,-2.0}) {
					p.incidence=incidence;
					snap_incidence(p,g);
					if ((incidence==-2.0)!=(p.i_incidence==-1))
						throw std::runtime_error("incidence mode selection");
					std::vector<double> top(g.NE), bottom(g.NE);
					compute_boundary_illumination(g.NE,p.i_incidence,g,rad,
					                              top.data(),bottom.data());
					for (int e=0; e<g.NE; ++e) {
						double top_J=0, bottom_J=0;
						int lit=0;
						for (int m=0; m<g.NA; ++m) {
							if (g.mu[m]<0) {
								const double I=incident_top_intensity(p.i_incidence,m,top[e]);
								if (I>0) ++lit;
								top_J+=0.5*g.wt[m]*I;
							} else {
								bottom_J+=0.5*g.wt[m]*bottom[e];
							}
						}
						close_to(top_J,rad.illum.I_corona[e],"top boundary J");
						close_to(bottom_J,rad.illum.I_disk[e],"bottom boundary J");
						if (lit != (incidence==-2.0 ? g.NA/2 : 1))
							throw std::runtime_error("number of illuminated downward angles");
					}
					std::printf("%s frac=%g incidence=%g lit=%d J_corona=%.9e J_disk=%.9e\n",
					            double_gauss ? "double-Gauss" : "full-Gauss",frac,incidence,
					            incidence==-2.0 ? g.NA/2 : 1,corona_J,disk_J);
				}
			}
		}
		rad.deallocate();
		return 0;
	} catch(const std::exception& e) {
		std::fprintf(stderr,"FAIL: %s\n",e.what());
		return 1;
	}
}
