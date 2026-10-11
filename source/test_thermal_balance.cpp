#include "thermal_balance.h"
#include "cell_transfer.h"
#include "compton_cross_section.h"
#include "test_local_response.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <limits>
#include <stdexcept>

static void require(bool ok,const char* message)
{
    if(!ok) throw std::runtime_error(message);
}

int main()
{
    try {
        char directory[]="/tmp/dao-thermal-kernel-XXXXXX";
        require(mkdtemp(directory),"temporary directory");
        setenv("KERNEL_DIR",directory,1);
        setenv("COMPTON_CACHE_DIR",directory,1);
        const int n=96;
        std::vector<double> E(n),J(n),j(n),reference(n),scratch(n);
        for(int i=0;i<n;++i) E[i]=100*std::exp(std::log(100.0)*i/(n-1));
        auto w=dao_thermal::energy_weights(n,E.data());
        const double temperatures[]={2.3e6,2.909326994e6,3.680069e6};
        const double T=std::sqrt(temperatures[0]*temperatures[1]),density=1e15;
        avgKernelCache table,exact;
        table.init(n,E.data(),0,nullptr,nullptr,1,3,temperatures);
        exact.init(n,E.data(),0,nullptr,nullptr,1,1,&T);
        dao_thermal::Scattering op(table,T,density),truth(exact,T,density);
        {
            const double mu[]={-0.5773502691896257,0.5773502691896257},aw[]={1.,1.};
            test_local_response(op,n,2,E.data(),w.data(),mu,aw);
            dao_thermal::Scattering zero(table,T,0);
            std::vector<double> incoming(2*n,1.),absorption(n,1.),emission(n),seed(n,0.);
            for(double tau:{1e-6,1.,1e6,1e12}) for(double S:{0.5,1.,2.}) {
                std::fill(emission.begin(),emission.end(),S);
                const auto response=dao_thermal::local_response(n,2,w.data(),mu,aw,tau,
                    absorption.data(),emission.data(),zero,incoming.data(),seed.data());
                const double phi=dao_rt::cell_weights(tau/std::abs(mu[0])).mean_incoming;
                double expected=0;
                for(int e=0;e<n;++e) {
                    expected+=phys::four_pi*w[e]*phi*(1-S);
                    require(std::abs(response.mean[e]-(phi+(1-phi)*S))<1e-14,"analytic grey cell intensity");
                }
                require(std::abs(response.residual-expected)<1e-13*response.scale,"analytic grey cell thermal residual");
                require(std::abs(response.relative()-(1-S)/(1+S))<1e-13,"thick cell root not hidden by trapped radiation");
            }
            std::printf("PASS absorbing cell analytic response and root through tau=1e12\n");
        }
        for(int i=0;i<n;++i) J[i]=std::exp(-E[i]/8000)/E[i];
        op.emissivity(J.data(),j.data()); truth.emissivity(J.data(),reference.data());
        double ni=0,no=0,difference=0,norm=0,q=0;
        for(int i=0;i<n;++i) {
            require(std::isfinite(j[i]) && j[i]>=0,"negative/nonfinite scattering emission");
            ni+=w[i]*op.opacity()[i]*J[i]/E[i]; no+=w[i]*j[i]/E[i];
            q+=phys::four_pi*w[i]*(op.opacity()[i]*J[i]-j[i]);
            difference+=w[i]*std::abs(j[i]-reference[i]); norm+=w[i]*reference[i];
        }
        require(std::abs(no/ni-1)<3e-10,"photon conservation");
        std::printf("diagnostic Q=%.15e operator Q=%.15e emission=%.15e\n",q,op.energy_gain(J.data()),phys::four_pi*norm);
        require(std::abs(q-op.energy_gain(J.data()))<1e-12*phys::four_pi*norm,"Compton diagnostic equals operator residual");
        require(difference/norm<0.01,"interpolation versus directly computed kernel exceeds 1 percent");
        std::printf("PASS photon number %.3e; interpolated/direct spectral L1 %.3e\n",no/ni-1,difference/norm);
        // Changing ionization must scale extinction, emission and Q together.
        // Include neutral and partially ionized limits, with an absolute
        // n_e*sigma reference to detect an accidental extra factor of 1.21.
        for(double fraction:{0.,0.1,1.21}) {
            const double electron_density=fraction*density;
            dao_thermal::Scattering scaled(table,T,electron_density);
            scaled.emissivity(J.data(),scratch.data());
            for(int i=0;i<n;++i) {
                const double expected=electron_density*compton_cross_section(E[i],T);
                require(std::abs(scaled.opacity()[i]-expected)<1e-12*std::max(expected,1e-100),"free-electron opacity normalization");
                require(std::abs(scratch[i]-fraction*j[i])<1e-12*std::max(fraction*j[i],1e-100),"free-electron emissivity scaling");
            }
            require(std::abs(scaled.energy_gain(J.data())-fraction*q)<1e-12*phys::four_pi*norm,"free-electron Compton exchange scaling");
        }
        for(double bad_density:{-1.,std::numeric_limits<double>::quiet_NaN()}) {
            bool rejected=false;
            try { dao_thermal::Scattering bad(table,T,bad_density); }
            catch(const std::runtime_error&) { rejected=true; }
            require(rejected,"invalid electron density accepted");
        }
        std::printf("PASS free-electron opacity/source/Q scaling, zero density, and invalid density rejection\n");
        // No stimulated scattering in this operator: Wien, not Planck, is its
        // stationary spectrum. Check every significant bin, including boundaries.
        const double kT=phys::k_B*T/phys::eV_to_erg;
        for(int i=0;i<n;++i) J[i]=E[i]*E[i]*E[i]*std::exp(-E[i]/kT);
        op.emissivity(J.data(),j.data());
        double worst=0;
        for(int i=0;i<n;++i) worst=std::max(worst,std::abs(j[i]/(op.opacity()[i]*J[i])-1));
        require(worst<3e-10,"Wien equilibrium / detailed balance");
        std::printf("PASS Wien equilibrium max relative residual %.3e\n",worst);
        // Temperature is continuous across a tabulated node. Use a symmetric
        // synthetic third node to evaluate both sides, including exact endpoint.
        dao_thermal::Scattering at(table,temperatures[1],density);
        dao_thermal::Scattering near(table,temperatures[1]*(1-1e-7),density);
        at.emissivity(J.data(),j.data()); near.emissivity(J.data(),scratch.data());
        difference=norm=0;
        for(int i=0;i<n;++i) { difference+=w[i]*std::abs(j[i]-scratch[i]); norm+=w[i]*j[i]; }
        require(difference/norm<1e-6,"temperature interpolation continuity");
        dao_thermal::Scattering above(table,temperatures[1]*(1+1e-7),density);
        above.emissivity(J.data(),scratch.data()); difference=0;
        for(int i=0;i<n;++i) difference+=w[i]*std::abs(j[i]-scratch[i]);
        require(difference/norm<1e-6,"temperature interpolation continuity above table node");
        bool rejected=false;
        try { dao_thermal::Scattering outside(table,temperatures[0]/2,density); }
        catch(const std::runtime_error&) { rejected=true; }
        require(rejected,"out-of-table temperature silently accepted");

        for(double start:{1.01e4,1e5,8e7}) {
            dao_thermal::TemperatureRoot root(start,1e4,1e8);
            double t=start;
            for(int k=0;k<60;++k) {
                double residual=1-std::pow(t/3e6,1.3);
                if(std::abs(residual)<1e-8) break;
                root.record(t,residual); t=root.next();
            }
            require(std::abs(t/3e6-1)<1e-7,"safeguarded temperature root");
        }
        rejected=false;
        try {
            dao_thermal::TemperatureRoot root(1e6,1e4,1e8);
            double t=1e6;
            for(int k=0;k<100;++k) { root.record(t,1);t=root.next(); }
        } catch(const std::runtime_error&) { rejected=true; }
        require(rejected,"no-root case marked converged");
        std::printf("PASS temperature roots and missing-root/out-of-table rejection\n");

        // Cover the cold end and the cooler temperatures encountered in the
        // benchmark atmosphere using newly integrated, intermediate-T kernels.
        for(double low:{1e4,1e5,1e6,1e7}) {
            const double nodes[]={low,low*std::pow(10.0,5.0/49.0)};
            const double middle=std::sqrt(nodes[0]*nodes[1]);
            avgKernelCache pair,direct;
            pair.init(n,E.data(),0,nullptr,nullptr,1,2,nodes);
            direct.init(n,E.data(),0,nullptr,nullptr,1,1,&middle);
            dao_thermal::Scattering interpolated(pair,middle,density),computed(direct,middle,density);
            for(int i=0;i<n;++i) J[i]=std::exp(-E[i]/8000)/E[i];
            interpolated.emissivity(J.data(),j.data()); computed.emissivity(J.data(),reference.data());
            ni=no=difference=norm=0;
            for(int i=0;i<n;++i) {
                ni+=w[i]*interpolated.opacity()[i]*J[i]/E[i]; no+=w[i]*j[i]/E[i];
                difference+=w[i]*std::abs(j[i]-reference[i]); norm+=w[i]*reference[i];
            }
            require(std::abs(no/ni-1)<3e-10,"cold/warm photon conservation");
            require(difference/norm<0.01,"cold/warm interpolation vs direct kernel");
            std::printf("PASS T=%.6e K: photon error %.3e; source L1 vs direct %.3e\n",middle,no/ni-1,difference/norm);
            pair.free_memory(); direct.free_memory();
        }

        // Independent finite-volume identity, all cells and BOTH slab faces.
        const int nd=5,nm=2,ne=1;
        double mu[]={-0.577350269,0.577350269},top[]={2},bottom[]={0.3};
        double kt[]={0.01,0.2,100,3,0.05},dr[]={0.2,0.7,0.4,1.0,0.5};
        double S[nd*nm],I[nd*nm],ft[nm],fb[nm],incoming[nd*nm];
        dao_rt::CellWeights weights[nd*nm];
        for(int d=0;d<nd;++d) for(int m=0;m<nm;++m) {
            S[d*nm+m]=0.1*(d+1);
            weights[d*nm+m]=dao_rt::cell_weights(kt[d]*dr[d]/std::abs(mu[m]));
        }
        dao_rt::formal_solution_cells(nd,nm,ne,0,mu,top,bottom,S,weights,I,ft,fb,incoming);
        require(incoming[0]==top[0] && incoming[(nd-1)*nm+1]==bottom[0],"captured boundary incoming rays");
        for(int d=0;d<nd;++d) for(int m=0;m<nm;++m) {
            const int next=d+(mu[m]<0 ? 1 : -1),i=d*nm+m;
            if(next>=0 && next<nd)
                require(std::abs(incoming[next*nm+m]-(weights[i].attenuation*incoming[i]+weights[i].emission*S[i]))<1e-14,
                        "captured internal incoming rays");
        }
        double volume=0;
        for(int d=0;d<nd;++d) for(int m=0;m<nm;++m)
            volume+=0.5*phys::four_pi*dr[d]*kt[d]*(I[d*nm+m]-S[d*nm+m]);
        double boundary=0;
        for(int m=0;m<nm;++m) boundary+=0.5*phys::four_pi*mu[m]*(ft[m]-fb[m]);
        require(std::abs(boundary+volume)<1e-12,"both-face finite-volume energy identity");
        std::printf("PASS both-face finite-volume identity %.3e\n",boundary+volume);
        table.free_memory(); exact.free_memory();
        std::filesystem::remove_all(directory);
    } catch(const std::exception& e) { std::fprintf(stderr,"FAIL: %s\n",e.what());return 1; }
}
