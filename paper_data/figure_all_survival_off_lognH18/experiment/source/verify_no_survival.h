#include <stdexcept>

// Verify line sources in BOTH the inner-temperature and full-column paths,
// including opaque/quench-prone lines, multiple lines per bin, zero-emission
// absorbers, and the unchanged pre-existing fluorescence and continuum.
static void verify_no_survival()
{
    RTGrids g; g.init_angle(); g.init_depth(1e-4,5.0,18,true);
    g.init_energy(100.0,1e4,3);
    RadField rad(g); rad.allocate();
    ModelParams p{};
    std::vector<std::vector<LineRec>> lines(g.ND_MID);
    for(int d=0;d<g.ND_MID;++d) {
        rad.n_e[d]=1e18;
        for(int j=0;j<3;++j) {
            LineRec r{};
            r.ip=100+j; r.bin=j%2; r.emiss=(j+1)*(d+1)*1000.;
            r.dE_eV=10.; r.E_eV=g.ene[r.bin];
            r.kappaL=1e6; r.kappaLo=2e6; r.y=100.+d; r.damp=.01;
            lines[d].push_back(r);
        }
        LineRec absorber=lines[d][0]; absorber.ip=500; absorber.emiss=0;
        lines[d].push_back(absorber);
    }
    const FrozenLineColumn column(g,lines);
    for(int path=0;path<2;++path) {
        for(int d=0;d<g.ND_MID;++d) {
            rad.line_heat[d]=0;
            for(int e=0;e<g.NE;++e) {
                rad.jnu[d][e]=3.; rad.jnu_line[d][e]=.7;
                rad.kabs[d][e]=1e4; rad.ksct[d][e]=1e-6;
            }
        }
        if(path==0) apply_line_escape(rad,g,lines,p,0);
        else for(int d=0;d<g.ND_MID;++d) apply_frozen_line_escape(d,rad,g,lines[d],column);
        for(int d=0;d<g.ND_MID;++d) {
            for(int e=0;e<g.NE;++e) {
                double expected=.7;
                for(const auto& r:lines[d]) if(r.bin==e) expected+=r.emiss/(r.dE_eV*phys::four_pi);
                if(std::abs(rad.jnu_line[d][e]-expected)>1e-12*std::max(1.,expected) ||
                   std::abs(rad.jnu[d][e]-(3.+expected))>1e-12*std::max(1.,expected) ||
                   rad.line_heat[d]!=0 || rad.kabs[d][e]!=1e4 || rad.ksct[d][e]!=1e-6)
                    throw std::runtime_error("P=1 line-source verification failed");
            }
        }
    }
    std::printf("PASS: all-line P=1 in local roots and full-column transfer; continuum, fluorescence and opacities preserved.\n");
}
