// Fixed-atmosphere Fe XXV experiment. No production file is changed.
// Include the unchanged atomic interface to call its exact line_probabilities.
#include "cloudy_interface_v2.cpp"
#include "cloudy_depth.h"
#include "compton_rt.h"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace {
const std::string output="results/fe_k_survival_a518ab75/transfer_control_20261009";
void check(bool ok,const char* why) { if(!ok) throw std::runtime_error(why); }
using Field=std::vector<double>;
Field flatten(double** a,const RTGrids& g) {
    Field v(long(g.ND_MID)*g.NE);
    for(int d=0;d<g.ND_MID;++d) std::copy(a[d],a[d]+g.NE,v.data()+long(d)*g.NE);
    return v;
}
void restore(double** a,const Field& v,const RTGrids& g) {
    for(int d=0;d<g.ND_MID;++d) std::copy(v.data()+long(d)*g.NE,v.data()+long(d+1)*g.NE,a[d]);
}
void load_moments(RadField& rad,const RTGrids& g,const std::string& path) {
    std::ifstream file(path); check(bool(file),"Missing moments input");
    std::string line; std::vector<int> counts(g.ND_MID,0);
    while(std::getline(file,line)) {
        if(line.empty() || line[0]=='#') continue;
        std::istringstream row(line); int d; double E,tau,T,J;
        check(bool(row>>E>>d>>tau>>T>>J),"Malformed moments row");
        check(d>=0 && d<g.ND_MID,"Invalid moments depth"); const int e=counts[d]++;
        check(e<g.NE && std::abs(E/g.ene[e]-1)<1e-6 && J>=0,"Moments grid mismatch");
        rad.J0[d][e]=J;
    }
    for(int n:counts) check(n==g.NE,"Incomplete moments input");
}
void write_spectrum(const std::string& name,const RadField& rad,const RTGrids& g) {
    std::ofstream f(output+"/"+name+".dat"); check(bool(f),"Cannot write spectrum");
    f<<std::setprecision(17)<<"# E_eV F_upper_out F_lower_out F_incoming I_top[8] I_bottom[8]\n"
      <<"# Fluxes: erg cm^-2 s^-1 eV^-1; intensities: erg cm^-2 s^-1 eV^-1 sr^-1\n";
    for(int e=0;e<g.NE;++e) {
        double up=0,down=0,in=0;
        for(int m=0;m<g.NA;++m) {
            const double c=0.5*phys::four_pi*g.wt[m]*std::abs(g.mu[m]);
            if(g.mu[m]>0) { up+=c*rad.Inu_top[m][e]; in+=c*rad.Inu_bottom[m][e]; }
            else { down+=c*rad.Inu_bottom[m][e]; in+=c*rad.Inu_top[m][e]; }
        }
        f<<g.ene[e]<<' '<<up<<' '<<down<<' '<<in;
        for(int m=0;m<g.NA;++m) f<<' '<<rad.Inu_top[m][e];
        for(int m=0;m<g.NA;++m) f<<' '<<rad.Inu_bottom[m][e];
        f<<'\n';
    }
}
void write_budget(const std::string& name,double wall,const RadField& rad,const RTGrids& g,
                  const std::vector<dao_thermal::Scattering>& ops,std::ofstream& summary) {
    const auto weights=dao_thermal::energy_weights(g.NE,g.ene);
    double incoming=0,up=0,down=0,residual=0,emission=0,ph_in=0,ph_out=0,ph_volume=0,ph_emit=0;
    std::ofstream f(output+"/budget_"+name+".dat"); f<<std::setprecision(17);
    f<<"# depth tau T ne absorption emission Q_Compton R dr\n";
    for(int d=0;d<g.ND_MID;++d) {
        const auto b=dao_thermal::budget(g.NE,weights.data(),rad.J0[d],rad.kabs[d],rad.jnu[d],ops[d]);
        residual+=b.residual()*g.dr[d]; emission+=b.emission*g.dr[d];
        f<<d<<' '<<g.tau_mid[d]<<' '<<rad.T_K[d]<<' '<<rad.n_e[d]<<' '
         <<b.absorption<<' '<<b.emission<<' '<<b.compton<<' '<<b.residual()<<' '<<g.dr[d]<<'\n';
        for(int e=0;e<g.NE;++e) {
            const double c=phys::four_pi*g.dr[d]*weights[e]/(g.ene[e]*phys::eV_to_erg);
            ph_volume+=c*(rad.kabs[d][e]*rad.J0[d][e]-rad.jnu[d][e]);
            ph_emit+=c*rad.jnu[d][e];
        }
    }
    for(int m=0;m<g.NA;++m) for(int e=0;e<g.NE;++e) {
        const double c=0.5*phys::four_pi*g.wt[m]*std::abs(g.mu[m])*weights[e];
        const double inf=c*(g.mu[m]<0 ? rad.Inu_top[m][e]:rad.Inu_bottom[m][e]);
        const double outf=c*(g.mu[m]>0 ? rad.Inu_top[m][e]:rad.Inu_bottom[m][e]);
        incoming+=inf; if(g.mu[m]>0) up+=outf; else down+=outf;
        ph_in+=inf/(g.ene[e]*phys::eV_to_erg); ph_out+=outf/(g.ene[e]*phys::eV_to_erg);
    }
    const double identity=(up+down-incoming+residual)/std::max({incoming,emission,1e-100});
    const double photons=(ph_out-ph_in+ph_volume)/std::max({ph_in,ph_emit,1e-100});
    summary<<name<<' '<<wall<<' '<<incoming<<' '<<up<<' '<<down<<' '<<emission<<' '
           <<residual<<' '<<identity<<' '<<photons<<std::endl;
    check(std::isfinite(identity) && std::abs(identity)<1e-6,"Boundary/volume energy identity failed");
    check(std::isfinite(photons) && std::abs(photons)<1e-6,"Photon accounting failed");
    std::printf("Completed %s: %.2fs; energy identity=%.3e, photon identity=%.3e\n",name.c_str(),wall,identity,photons);
    std::fflush(stdout);
}
}

int main() {
    try {
        const auto started=std::chrono::steady_clock::now();
        ModelParams p{}; p.nh=15; p.Afe=1; p.zeta=3; p.frac=-1; p.incidence=0.7071;
        p.Gamma=2; p.kT_e=60; p.kT_bb=0.01; p.kT_disk=0.35; p.E_lo_cut=0.1;
        p.E_cut=-1; p.taup=-1; p.ktype=1; p.reflionx_norm=true; p.angsca=false;
        p.maxiter=2000; // Retain the production 1e-7 scattering convergence criterion.
        std::snprintf(p.corona,sizeof(p.corona),"nthcomp");
        std::snprintf(p.run_dir,sizeof(p.run_dir),"%s",output.c_str());
        RTGrids g; g.init_angle(); g.init_depth(1e-4,5,p.nh,true);
        // Match production's parent-process atomic initialization before
        // forking cells. Loading only the energy mesh skips shared Cloudy
        // database initialization and can alter the non-isoelectronic line list.
        bootstrap_cloudy_energy_grid(g,(output+"/energy.dat").c_str());
        snap_incidence(p,g);
        RadField rad(g); rad.allocate(); rad.illum.compute(p);
        std::ifstream bfile("results/a518ab75/thermal_budget_iter33.dat");
        check(bool(bfile),"Missing accepted-temperature budget");
        Field fixed(g.ND_MID); std::string line; int count=0;
        while(std::getline(bfile,line)) {
            if(line.empty() || line[0]=='#') continue;
            std::istringstream row(line); int d; double tau,T,a,j,c,R,rel,dr;
            check(bool(row>>d>>tau>>T>>a>>j>>c>>R>>rel>>dr),"Malformed original budget");
            check(d==count++ && d<g.ND_MID,"Original depth mismatch");
            check(std::abs(tau/g.tau_mid[d]-1)<1e-12 && std::abs(dr/g.dr[d]-1)<1e-12,"Original width mismatch");
            fixed[d]=T;
        }
        check(count==g.ND_MID,"Incomplete original atmosphere");
        load_moments(rad,g,"results/a518ab75/moments_iter032.dat");
        rad.compute_ionization_parameter(p.nh);
        std::vector<std::vector<LineRec>> records;
        run_cloudy_depths(rad,g,p,33,records,4,fixed.data());
        // Strict read-only load: never silently regenerate or alter a shared cache.
        avgKernelCache cache; cache.NE=g.NE; cache.NT=N_T_CACHE;
        cache.T_grid=new double[cache.NT];
        for(int i=0;i<cache.NT;++i) cache.T_grid[i]=std::pow(10.0,std::log10(T_CACHE_LO)+i*(std::log10(T_CACHE_HI)-std::log10(T_CACHE_LO))/(cache.NT-1));
        check(cache.load("/Users/xe26734/Downloads/Models/DAO_KERNEL_TABLE/avgkernel_norm_NE2768_NT50.bin",g.ene),"Existing kernel does not match; no cache computation authorized here");
        const auto ops=make_scattering_column(rad,g,cache);
        apply_line_escape(rad,g,records,p,33);
        std::fill(rad.line_heat,rad.line_heat+g.ND_MID,0);
        apply_line_escape(rad,g,records,p,33,true);
        const Field material=flatten(rad.jnu,g), absorption=flatten(rad.kabs,g), scattering=flatten(rad.ksct,g);
        const Field temperatures(rad.T_K,rad.T_K+g.ND_MID), density(rad.n_e,rad.n_e+g.ND_MID);
        const FrozenLineColumn column(g,records);
        Field bypass=material;
        std::ofstream source(output+"/target_sources.dat"); source<<std::setprecision(17);
        source<<"# depth tau dr line_ip E_bin_eV dE_eV emiss_thin P emiss_surviving\n";
        for(int d=0;d<g.ND_MID;++d) for(const auto& r:records[d]) {
            if(r.ip!=29773 && r.ip!=29767) continue;
            const auto& v=column.depths.at(r.ip)[d]; const double own=0.5*r.kappaL*g.dr[d];
            const auto q=line_probabilities(r,rad,d,std::max(0.0,v[0]+own),std::max(0.0,v[1]+own));
            const long i=long(d)*g.NE+r.bin; const double j=r.emiss/(phys::four_pi*r.dE_eV);
            bypass[i]+=j*(1-q.P);
            source<<d<<' '<<g.tau_mid[d]<<' '<<g.dr[d]<<' '<<r.ip<<' '<<r.E_eV<<' '
                  <<r.dE_eV<<' '<<r.emiss<<' '<<q.P<<' '<<r.emiss*q.P<<'\n';
        }
        source.close();
        load_moments(rad,g,"results/a518ab75/moments_iter033.dat");
        const Field original_final_J=flatten(rad.J0,g);
        std::ofstream budgets(output+"/case_budgets.dat"); budgets<<std::setprecision(17);
        budgets<<"# case wall_seconds F_in F_upper F_lower volume_material_emission volume_R energy_identity photon_identity\n";
        auto solve=[&](const std::string& name,const Field& emiss) {
            std::printf("Starting %s\n",name.c_str()); std::fflush(stdout);
            restore(rad.jnu,emiss,g);
            restore(rad.J0,original_final_J,g);
            const auto start=std::chrono::steady_clock::now();
            compton_rt_solve(rad,g,p,ops);
            const double wall=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
            check(flatten(rad.kabs,g)==absorption && flatten(rad.ksct,g)==scattering,"Opacity changed in fixed-atmosphere experiment");
            check(Field(rad.T_K,rad.T_K+g.ND_MID)==temperatures && Field(rad.n_e,rad.n_e+g.ND_MID)==density,"Temperature/density changed");
            write_spectrum(name,rad,g); write_budget(name,wall,rad,g,ops,budgets);
        };
        solve("full_survival",material);
        solve("full_wx_P1",bypass);
        budgets.close();
        std::ofstream done(output+"/completed.json");
        done<<"{\"state\":\"completed\",\"cases\":2,\"thermal_iterations\":0,\"maxiter\":2000,\"RT_tolerance\":1e-7,\"consecutive_passes\":3,\"wall_seconds\":"
            <<std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()<<"}\n";
        rad.deallocate(); cache.free_memory();
        return 0;
    } catch(const std::exception& e) {
        std::fprintf(stderr,"Fixed-atmosphere experiment failed: %s\n",e.what()); return 1;
    }
}
