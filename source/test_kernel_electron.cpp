// Include the kernel to access its differential cross section and the frozen
// 32-point rule without adding a reference-only entry point to the public API.
#include "compton_kernel.cpp"
#include <stdexcept>
#include <cstring>

static double reference32(double x, double x1, double c, double inverse)
{
	const double q=x*x1*(1-c);
	const double Q=std::sqrt((x-x1)*(x-x1)+2*q);
	const double gamma_min=(x-x1+Q*std::sqrt(1+2/q))/2;
	double integral=0;
	for(int i=0;i<32;++i)
		integral+=gl_weight[i]*fexact(x,x1,c,gamma_min+gl_node[i]/inverse);
	return 3/(32*M_PI)*integral*std::exp((1-gamma_min)*inverse)/bk2_exp(inverse);
}

// Frozen outputs from the user's local source before this restoration.
// These check reproduction of that implementation, not 32-point accuracy.
static void check_local_baseline()
{
    const double samples[][5]={
        {4.999,1e6,1,-0.999,0.00025821689262864799},
        {5.0,1e6,1,-0.999,0.00025670135653932499},
        {5.001,1e6,1,-0.999,0.0002564513744954868},
        {29.999,1e6,1,-0.999,1.1961241135004444e-16},
        {30.0,1e6,1,-0.999,1.1947197194129635e-16},
        {30.001,1e6,1,-0.999,1.1933169731797012e-16},
        {999.999,1e3,1,0.5,480.05082527429374},
        {1000.0,1e3,1,0.5,480.05106572782665},
        {1000.001,1e3,1,0.5,480.05130618244635},
        {10.0,0.1,1,0.5,422395.30900423176},
        {10.0,1e7,0.1,-0.9,2.9375061118746683e-82},
        {1000.0,0.1,1,0.5,4802804.1945089344}
    };
    for(const auto& s:samples) {
        const double x=s[1]/phys::m_e_c2_eV;
        const double actual=profil_exact(x,x*s[2],s[3],s[0]);
        if(!std::isfinite(actual)||std::abs(actual/s[4]-1)>1e-11)
            throw std::runtime_error("electron integral differs from frozen local baseline");
    }
    std::puts("PASS 12 frozen local electron profiles (relative tolerance 1e-11)");
}

int main(int argc,char** argv)
{
    const bool local_baseline=argc==2 && std::strcmp(argv[1],"--local-baseline")==0;
    if(argc!=1 && !local_baseline) {
        std::fprintf(stderr,"Usage: %s [--local-baseline]\n",argv[0]); return 2;
    }
	try {
		if(local_baseline) check_local_baseline();
		double worst=0, worst_case[5]={}; size_t compared=0, exceeding=0;
		for(double inverse : {1.0,4.999,5.0,6.0,10.0,29.999,30.0,30.001,
		                       100.0,300.0,999.999,1000.0,1000.001,1e4,1e5,1e6})
		for(double E : {0.1,1.0,10.0,100.0,1000.0,6400.0,1e4,1e5,3e5,1e6,1e7})
		for(double ratio : {0.1,0.5,0.9,0.99,1.0,1.01,1.1,2.0,10.0})
		for(double c : {-0.999,-0.9,-0.5,0.0,0.5,0.9,0.99,0.999}) {
			const double x=E/phys::m_e_c2_eV,x1=x*ratio;
			const double actual=profil_exact(x,x1,c,inverse);
			const double reference=reference32(x,x1,c,inverse);
			if(!std::isfinite(actual)||!std::isfinite(reference)||actual<0||reference<0)
				throw std::runtime_error("nonfinite or negative electron profile");
			if(reference<1e-100) continue; // underflow tails have negligible weight
			const double error=std::abs(actual/reference-1);
			if(error>worst) {
                worst=error;
                worst_case[0]=inverse; worst_case[1]=E; worst_case[2]=E*ratio;
                worst_case[3]=c; worst_case[4]=reference;
            }
            ++compared;
			if(error>1e-5) {
				++exceeding;
				if(local_baseline) continue;
				std::fprintf(stderr,"inverse=%g E=%g E1=%g cos=%g error=%g\n",
				             inverse,E,E*ratio,c,error);
				throw std::runtime_error("electron quadrature differs from 32-point baseline");
			}
		}
		if(local_baseline) {
            std::printf("AUDIT 32-point reference: %zu profiles, max relative difference %.9g; %zu exceed 1e-5\n",
                        compared,worst,exceeding);
            std::printf("Largest difference: inverse=%g E=%g eV E1=%g eV cos=%g reference=%.9g\n",
                        worst_case[0],worst_case[1],worst_case[2],worst_case[3],worst_case[4]);
        } else
            std::printf("PASS %zu electron profiles; maximum relative difference %.9g (limit 1e-5)\n",
                        compared,worst);
		return 0;
	} catch(const std::exception& e) {
		std::fprintf(stderr,"FAIL: %s\n",e.what()); return 1;
	}
}
