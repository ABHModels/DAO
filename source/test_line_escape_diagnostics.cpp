#include "cloudy_interface.h"
#include "constants.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <unistd.h>

static void check(bool ok,const char* message)
{
	if(!ok) throw std::runtime_error(message);
}

static std::vector<double> state(const RadField& r,const RTGrids& g)
{
	std::vector<double> values;
	for(int d=0;d<g.ND_MID;++d) {
		values.push_back(r.line_heat[d]);
		for(int e=0;e<g.NE;++e) {
			values.push_back(r.jnu[d][e]); values.push_back(r.jnu_line[d][e]);
		}
	}
	return values;
}

int main()
{
	char temp[]="/tmp/dao-line-diag-XXXXXX";
	if(!mkdtemp(temp)) return 1;
	try {
		RTGrids g; g.init_angle(); g.init_energy(100,400,3);
		const double mids[]={0.05,0.15},widths[]={1e8,1e8};
		g.init_depth_from_zones(2,mids,widths);
		RadField r(g); r.allocate(); ModelParams p{};
		std::snprintf(p.run_dir,sizeof(p.run_dir),"%s",temp);
		for(int d=0;d<g.ND_MID;++d) {
			r.n_e[d]=1e15;
			for(int e=0;e<g.NE;++e) {
				r.jnu[d][e]=0.25; r.jnu_line[d][e]=0.125;
				r.kabs[d][e]=1e-9; r.ksct[d][e]=1e-9;
			}
		}
		std::vector<std::vector<LineRec>> lines(g.ND_MID);
		for(int d=0;d<g.ND_MID;++d) for(int thin=0;thin<2;++thin) {
			LineRec line{}; line.ip=thin+1; line.bin=thin; line.E_eV=g.ene[thin];
			line.dE_eV=g.wid[thin]; line.emiss=(d+1)*2.0; line.redis=-1; // Cloudy ipCRD
			line.kappaL=line.kappaLo=thin ? 0 : 4e-8;
			line.label=thin ? "thin" : "thick";
			lines[d].push_back(line);
		}
		apply_line_escape(r,g,lines,p,1);
		check(!std::filesystem::exists(std::string(temp)+"/line_escape_lines_iter001.dat"),
		      "ordinary pass unexpectedly emitted diagnostics");
		const auto before=state(r,g);
		apply_line_escape(r,g,lines,p,2,true);
		check(state(r,g)==before,"diagnostics changed emissivity or heating");
		apply_line_escape(r,g,lines,p,2,true);
		check(state(r,g)==before,"repeated diagnostics changed emissivity or heating");
		std::ifstream f(std::string(temp)+"/line_escape_lines_iter002.dat");
		check(bool(f),"missing per-line diagnostics");
		std::string text; int rows=0;
		while(std::getline(f,text)) {
			if(text.empty()||text[0]=='#') continue;
			std::istringstream row(text); int rank; long ip;
			double E,wave,thin,escaped,heat;
			check(bool(row>>rank>>ip>>E>>wave>>thin>>escaped>>heat),"bad diagnostic row");
			check(ip==1||ip==2,"unexpected diagnostic line");
			check(std::isfinite(escaped)&&std::isfinite(heat)&&escaped>=0&&heat>=0,
			      "invalid line power");
			check(std::abs(escaped+heat-thin)<2e-8*thin,"line energy branching budget");
			double actual=0;
			for(int d=0;d<g.ND_MID;++d)
				actual+=(r.jnu_line[d][ip-1]-0.125)*g.wid[ip-1]*phys::four_pi;
			check(std::abs(actual-escaped)<2e-8*thin,"diagnostic differs from injected line power");
			if(ip==1) {
				check(heat>0,"thick line did not exercise destruction");
				check(std::abs(r.line_heat[0]+r.line_heat[1]-heat)<2e-8*thin,
				      "diagnostic differs from injected heating");
			} else check(std::abs(escaped-thin)<2e-8*thin&&heat==0,"thin free-escape limit");
			++rows;
		}
		check(rows==2,"missing line diagnostic rows");
		std::ifstream selected(std::string(temp)+"/line_escape_selected_iter002.dat");
		check(bool(selected),"missing selected depth diagnostics");
		int depth_rows=0;
		while(std::getline(selected,text)) if(!text.empty()&&text[0]!='#') ++depth_rows;
		check(depth_rows==4,"missing selected depth rows");
		std::vector<std::vector<LineRec>> empty(g.ND_MID);
		apply_line_escape(r,g,empty,p,3,true);
		check(state(r,g)==before,"empty-line diagnostics changed emissivity");
		r.deallocate(); std::filesystem::remove_all(temp);
		std::puts("PASS line diagnostic isolation, output, thin limit and energy budget");
		return 0;
	} catch(const std::exception& e) {
		std::fprintf(stderr,"FAIL: %s (files in %s)\n",e.what(),temp); return 1;
	}
}
