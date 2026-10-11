#ifndef DAO_TEST_LOCAL_RESPONSE_H
#define DAO_TEST_LOCAL_RESPONSE_H
#include "thermal_balance.h"
#include "cell_transfer.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>

// Independently compare the local response with the physical material/Compton
// budget and the flux through the two cell faces, in either scattering mode.
inline void test_local_response(const dao_thermal::Scattering& op,int ne,int na,
                               const double* E,const double* ew,const double* mu,const double* aw)
{
    const bool directional=op.angles()>1;
    std::vector<double> incoming(na*ne),seed(directional ? na*ne : ne,0),absorption(ne),emission(ne);
    for(int e=0;e<ne;++e) {
        absorption[e]=0.2*op.opacity()[e];
        emission[e]=absorption[e]*std::exp(-E[e]/3000)/E[e];
        for(int m=0;m<na;++m) incoming[m*ne+e]=(1+mu[m])*std::exp(-E[e]/8000)/E[e];
    }
    for(double tau:{1e-6,1.,1e3}) {
        const double width=tau/(absorption[ne/2]+op.opacity()[ne/2]);
        auto response=dao_thermal::local_response(ne,na,ew,mu,aw,width,absorption.data(),emission.data(),op,incoming.data(),seed.data());
        auto b=dao_thermal::budget(ne,ew,response.mean.data(),absorption.data(),emission.data(),op,
                                   directional ? response.intensity.data() : nullptr);
        if(std::abs(b.residual()-response.residual)>2e-8*b.scale())
            throw std::runtime_error("local response disagrees with material/Compton energy budget");
        std::vector<double> j(seed.size());
        op.emissivity(directional ? response.intensity.data() : response.mean.data(),j.data());
        double flux=0;
        for(int m=0;m<na;++m) for(int e=0;e<ne;++e) {
            const int i=m*ne+e;
            const double k=absorption[e]+op.opacity()[e];
            const auto w=dao_rt::cell_weights(k*width/std::abs(mu[m]));
            const double source=(emission[e]+j[directional ? i : e])/k;
            const double outgoing=w.attenuation*incoming[i]+w.emission*source;
            flux+=0.5*phys::four_pi*aw[m]*std::abs(mu[m])*ew[e]*(outgoing-incoming[i]);
        }
        if(std::abs(flux+response.residual*width)>2e-8*response.scale*width)
            throw std::runtime_error("local response violates both-face flux identity");
        std::vector<double> large_seed(seed.size(),1.);
        auto other=dao_thermal::local_response(ne,na,ew,mu,aw,width,absorption.data(),emission.data(),op,incoming.data(),large_seed.data());
        if(std::abs(other.residual-response.residual)>2e-8*response.scale)
            throw std::runtime_error("local response depends on its initial field");
        std::printf("PASS %s local cell tau=%g: physical residual, face identity, seed independence\n",directional ? "directional" : "mean",tau);
    }
}
#endif
