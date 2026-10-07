// Access both private benchmark writers without expanding the production API.
#include "test_rt.cpp"
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>

int main()
{
	char temp[]="/tmp/dao-test-output-XXXXXX";
	if(!mkdtemp(temp)) return 1;
	try {
		RTGrids g; g.init_angle(); g.init_energy(100,400,3);
		RadField r(g); r.allocate(); ModelParams p{};
		std::snprintf(p.run_dir,sizeof(p.run_dir),"%s",temp);
		for(int m=0;m<g.NA;++m) for(int e=0;e<g.NE;++e) {
			r.Inu[0][m][e]=123456789;
			r.Inu_top[m][e]=10*(m+1)+(e+1);
		}
		save_emergent_compps(r,g,p,60); save_emergent_tavg(r,g,p,60);
		for(const char* name : {"emergent_compps.dat","emergent_tavg.dat"}) {
			std::ifstream f(std::string(temp)+"/"+name); std::string line; int e=0;
			if(!f) throw std::runtime_error("missing benchmark spectrum");
			while(std::getline(f,line)) {
				if(line.empty()||line[0]=='#') continue;
				if(e>=g.NE) throw std::runtime_error("extra spectrum row");
				std::istringstream row(line); double energy,corona,disk;
				row>>energy>>corona>>disk;
				for(int m=0;m<g.NA;++m) {
					double value;
					if(!(row>>value)||value!=r.Inu_top[m][e])
						throw std::runtime_error("benchmark must save upper-face intensity");
				}
				++e;
			}
			if(e!=g.NE) throw std::runtime_error("missing spectrum rows");
		}
		r.deallocate(); std::filesystem::remove_all(temp);
		std::puts("PASS both benchmark spectra save the true upper-face intensity");
		return 0;
	} catch(const std::exception& e) {
		std::fprintf(stderr,"FAIL: %s (files in %s)\n",e.what(),temp); return 1;
	}
}
