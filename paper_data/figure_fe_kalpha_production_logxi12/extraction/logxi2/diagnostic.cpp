// Reproduce accepted Cloudy state only; do not solve temperature or transfer.
// Rename the unchanged production extraction routine and append diagnostics.
#define extract_cloudy_output extract_cloudy_output_original
#include "cloudy_interface_v2.cpp"
#undef extract_cloudy_output
#include "cloudy_depth.h"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <chrono>

namespace {
const std::string out="results/fe_lowion_emissivity_effac594_20261009";
void check(bool ok,const char* msg) { if(!ok) throw std::runtime_error(msg); }
}

void extract_cloudy_output(int d,RadField& rad,const RTGrids& g,int iter,
                           std::vector<LineRec>& recs) {
    extract_cloudy_output_original(d,rad,g,iter,recs);
    std::ofstream f(out+"/cells/fluorescence_"+std::to_string(d)+".dat");
    f<<std::setprecision(17)
     <<"# d tau_ref dz_cm T_K ne_cm-3 Fe_cm-3 yield_index parent_stage emitting_stage shell_zero E_line_eV E_bin_eV n_parent_cm-3 Gamma_shell_s-1 yield photons_cm-3_s-1 power_line_erg_cm-3_s-1 power_grid_erg_cm-3_s-1\n";
    const auto& y=t_yield::Inst();
    for(long k=0;k<y.nlines();++k) {
        if(y.nelem(k)!=ipIRON) continue;
        const double E=y.energy(k)*phys::eV_per_Ryd;
        if(E<6200 || E>8200) continue;
        const long ip=y.ipoint(k)-1;
        check(ip>=0 && ip<rfield.nflux,"Bad fluorescence grid index");
        const double ni=dense.xIonDense[ipIRON][y.ion(k)];
        const double rate=ionbal.PhotoRate_Shell[ipIRON][y.ion(k)][y.nshell(k)][0];
        const double photons=ni*rate*y.yield(k);
        check(std::isfinite(photons) && photons>=0,"Invalid fluorescence production");
        const double Ebin=rfield.anu(ip)*phys::eV_per_Ryd;
        f<<d<<' '<<g.tau_mid[d]<<' '<<g.dr[d]<<' '<<rad.T_K[d]<<' '<<rad.n_e[d]<<' '
         <<dense.gas_phase[ipIRON]<<' '<<k<<' '<<y.ion(k)+1<<' '<<y.ion_emit(k)+1<<' '
         <<y.nshell(k)<<' '<<E<<' '<<Ebin<<' '<<ni<<' '<<rate<<' '<<y.yield(k)<<' '
         <<photons<<' '<<photons*E*phys::eV_to_erg<<' '
         <<photons*Ebin*phys::eV_to_erg<<'\n';
    }
    f.close(); check(bool(f),"Cannot write fluorescence diagnostics");
    std::ofstream b(out+"/cells/bound_bound_"+std::to_string(d)+".dat");
    b<<std::setprecision(17)<<"# d line_ip emitting_stage E_line_eV E_bin_eV emiss_thin_erg_cm-3_s-1\n";
    for(const auto& r:recs) {
        auto tr=LineSave.lines[r.ip].getTransition();
        if(!tr.associated() || (*tr.Hi()).nelem()!=26) continue;
        const double E=tr.EnergyRyd()*phys::eV_per_Ryd;
        if(E<6200 || E>8200) continue;
        b<<d<<' '<<r.ip<<' '<<(*tr.Hi()).IonStg()<<' '<<E<<' '<<r.E_eV<<' '<<r.emiss<<'\n';
    }
    b.close(); check(bool(b),"Cannot write bound-bound diagnostics");
}

int main() {
    try {
        const auto started=std::chrono::steady_clock::now();
        ModelParams p{}; p.nh=15; p.Afe=1; p.zeta=2; p.save_iron=true;
        std::snprintf(p.run_dir,sizeof(p.run_dir),"%s",out.c_str());
        RTGrids g; g.init_angle(); g.init_depth(1e-4,5,p.nh,true);
        bootstrap_cloudy_energy_grid(g,(out+"/energy.dat").c_str());
        RadField rad(g); rad.allocate();
        std::vector<double> fixed(g.ND_MID);
        std::ifstream budget("results/effac594/thermal_budget_iter35.dat");
        check(bool(budget),"Missing accepted-temperature budget");
        std::string line; int count=0;
        while(std::getline(budget,line)) {
            if(line.empty() || line[0]=='#') continue;
            std::istringstream s(line); int d; double tau,T,a,j,c,R,rel,dz;
            check(bool(s>>d>>tau>>T>>a>>j>>c>>R>>rel>>dz),"Malformed budget");
            check(d==count++ && d<g.ND_MID,"Bad depth index");
            check(std::abs(tau/g.tau_mid[d]-1)<1e-12 && std::abs(dz/g.dr[d]-1)<1e-12,"Depth grid mismatch");
            fixed[d]=T;
        }
        check(count==g.ND_MID,"Incomplete temperature input");
        std::ifstream moments("results/effac594/moments_iter034.dat");
        check(bool(moments),"Missing preceding outer mean intensity");
        std::vector<int> counts(g.ND_MID,0);
        while(std::getline(moments,line)) {
            if(line.empty() || line[0]=='#') continue;
            std::istringstream s(line); int d; double E,tau,T,J;
            check(bool(s>>E>>d>>tau>>T>>J),"Malformed moments input");
            check(d>=0 && d<g.ND_MID,"Invalid moments depth");
            int e=counts[d]++;
            check(e<g.NE && std::abs(E/g.ene[e]-1)<1e-6 && J>=0,"Energy grid mismatch");
            rad.J0[d][e]=J;
        }
        for(int n:counts) check(n==g.NE,"Incomplete moments input");
        rad.compute_ionization_parameter(p.nh);
        std::vector<std::vector<LineRec>> recs;
        run_cloudy_depths(rad,g,p,35,recs,4,fixed.data());
        const FrozenLineColumn column(g,recs);
        std::ofstream bb(out+"/bound_bound_survival.dat");
        bb<<std::setprecision(17)<<"# d tau_ref dz_cm line_ip E_bin_eV emiss_thin P emiss_surviving\n";
        for(int d=0;d<g.ND_MID;++d) for(const auto& r:recs[d]) {
            if(r.E_eV<6150 || r.E_eV>8250 || r.label.rfind("Fe",0)!=0) continue;
            const auto& v=column.depths.at(r.ip)[d];
            const double own=0.5*r.kappaL*g.dr[d];
            auto prob=line_probabilities(r,rad,d,std::max(0.0,v[0]+own),std::max(0.0,v[1]+own));
            bb<<d<<' '<<g.tau_mid[d]<<' '<<g.dr[d]<<' '<<r.ip<<' '<<r.E_eV<<' '
              <<r.emiss<<' '<<prob.P<<' '<<r.emiss*prob.P<<'\n';
        }
        bb.close(); check(bool(bb),"Cannot write source survival diagnostics");
        std::ofstream ions(out+"/iron_fractions_replayed.dat");
        ions<<std::setprecision(17)<<"# depth tau_ref T_K ne_cm-3 Fe_cm-3 Fe_I ... Fe_XXVII\n";
        for(int d=0;d<g.ND_MID;++d) {
            ions<<d<<' '<<g.tau_mid[d]<<' '<<rad.T_K[d]<<' '<<rad.n_e[d]<<' '<<rad.iron[d].gas_density;
            for(double f:rad.iron[d].fraction) ions<<' '<<f;
            ions<<'\n';
        }
        ions.close(); check(bool(ions),"Cannot write ion validation");
        std::ofstream done(out+"/completed.json");
        done<<std::setprecision(17)<<"{\"state\":\"completed\",\"depth_cells\":"<<g.ND_MID
            <<",\"temperature_iteration\":false,\"transfer_solve\":false,\"wall_seconds\":"
            <<std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()<<"}\n";
        std::printf("Completed local-emissivity diagnostics. No temperature or transfer iteration.\n");
        rad.deallocate(); return 0;
    } catch(const std::exception& e) {
        std::fprintf(stderr,"Diagnostic failed: %s\n",e.what()); return 1;
    }
}
