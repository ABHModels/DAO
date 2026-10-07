// Include the kernel to access its differential cross section and the frozen
// 32-point rule without adding a reference-only entry point to the public API.
#include "compton_kernel.cpp"
#include <stdexcept>

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

int main()
{
	try {
		double worst=0; size_t compared=0;
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
			worst=std::max(worst,error); ++compared;
			if(error>1e-5) {
				std::fprintf(stderr,"inverse=%g E=%g E1=%g cos=%g error=%g\n",
				             inverse,E,E*ratio,c,error);
				throw std::runtime_error("electron quadrature differs from 32-point baseline");
			}
		}
		std::printf("PASS %zu electron profiles; maximum relative difference %.9g (limit 1e-5)\n",
		            compared,worst);
		return 0;
	} catch(const std::exception& e) {
		std::fprintf(stderr,"FAIL: %s\n",e.what()); return 1;
	}
}
