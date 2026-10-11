// Independent fully ionized RT test: n_e=1.21*nH must scale extinction
// and redistribution identically. Production density is tested in test_thermal_*.
#include "source.h"
#include "rt_grids.h"
#include "compton_cross_section.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <stdexcept>
#include <vector>

static void close_to(double actual, double expected)
{
    if (!std::isfinite(actual) || !std::isfinite(expected) ||
        std::abs(actual - expected) > 1e-11 * std::max(1e-100, std::abs(expected)))
        throw std::runtime_error("density regression: nonfinite value or mismatch");
}

int main()
{
    char temp[] = "/tmp/dao-density-test-XXXXXX";
    if (!mkdtemp(temp)) return 1;
    setenv("COMPTON_CACHE_DIR", temp, 1);
    try {
        RTGrids g;
        g.init_angle();
        g.init_depth(1e-4, 5.0, 15.0);
        double E = 10, opacity;
        compute_compton_opacity(&opacity, 1, &E, 1e4, 1.21e15);
        double tau = 0;
        for (int d = 0; d < g.ND_MID; ++d) tau += opacity * g.dr[d];
        close_to(tau, 5.0);

        g.init_energy(10, 1e4, 32);
        const int ND = 4, NE = g.NE, NM = g.NA;
        const double T = 1e7, temperatures[ND] = {T,T,T,T};
        const double nh[ND] = {1e15, 1e14, 1e13, 0.0};
        KernelCache directional;
        avgKernelCache averaged;
        directional.init(NE, g.ene, NM, g.mu, g.wt, 1, 1, &T);
        averaged.init(NE, g.ene, NM, g.mu, g.wt, 1, 1, &T);
        std::vector<double> j(ND*NE,0), absorption(ND*NE), scattering(ND*NE);
        std::vector<double> intensity(ND*NM*NE,1), mean(ND*NE,1), result(ND*NM*NE);
        const double *jp[ND], *ap[ND], *sp[ND];
        for (int d=0; d<ND; ++d) {
            jp[d] = j.data()+d*NE;
            ap[d] = absorption.data()+d*NE;
            sp[d] = scattering.data()+d*NE;
            compute_compton_opacity(scattering.data()+d*NE, NE, g.ene, T,
                                    phys::reference_electrons_per_hydrogen*nh[d]);
            for (int e=0;e<NE;++e)
                close_to(sp[d][e], 1.21*nh[d]*compton_cross_section(g.ene[e],T));
        }
        for (bool angle_mean : {false,true}) {
            for (bool pure_scattering : {false,true}) {
                std::fill(absorption.begin(),absorption.end(),pure_scattering ? 0 : 1e-8);
                const int count = pure_scattering ? ND-1 : ND;
                for (int d=0;d<count;++d) {
                    if (angle_mean)
                        avgcompute_source_function(1,NM,NE,mean.data()+d*NE,averaged.x_grid,
                            averaged,temperatures+d,jp+d,ap+d,sp+d,nh[d],result.data()+d*NM*NE);
                    else
                        compute_source_function(1,NM,NE,intensity.data()+d*NM*NE,directional.x_grid,
                            g.wt,directional,temperatures+d,jp+d,ap+d,sp+d,nh[d],result.data()+d*NM*NE);
                }
                bool positive = false;
                for (int m=0;m<NM;++m) for (int e=0;e<NE;++e) {
                    const double base = result[m*NE+e];
                    positive = positive || base>0;
                    if (!std::isfinite(base) || base<0) throw std::runtime_error("invalid source");
                    for (int d=1;d<count;++d) {
                        const double value = result[(d*NM+m)*NE+e];
                        if (pure_scattering) close_to(value,base);
                        else close_to(value*(ap[d][e]+sp[d][e]),
                                      base*(ap[0][e]+sp[0][e])*
                                      nh[d]/nh[0]);
                    }
                }
                if (!positive) throw std::runtime_error("vacuous source test");
            }
        }
        directional.free_memory();
        averaged.free_memory();
        std::filesystem::remove_all(temp);
        puts("PASS: requested Thomson depth, 1.21*nH opacity, both source density scalings, pure-scattering density invariance, zero-hydrogen limit (tolerance 1e-11)");
        return 0;
    } catch (const std::exception& e) {
        std::fprintf(stderr,"FAIL: %s\n",e.what());
        std::filesystem::remove_all(temp);
        return 1;
    }
}
