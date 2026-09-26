#include "compton_kernel.h"
#include "avg_compton_kernel.h"
#include "kernel_quadrature.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <stdexcept>

template<int Order> double directional_reference(double x,double mu,double x1,
                                                  double mu1,double T,int ktype)
{
    const double inverse=511000.0/(8.617333262e-5*T);
    const double a=mu*mu1,b=std::sqrt(1-mu*mu)*std::sqrt(1-mu1*mu1);
    const double peak=1.0-std::abs(x-x1)/(x*x1);
    std::array<double,5> cuts{0.0,M_PI}; size_t count=2;
    if(b>0 && peak>=a-b && peak<=a+b) {
        const double phi=std::acos(std::clamp((peak-a)/b,-1.0,1.0));
        const double dc=4*std::sqrt(2.0/inverse)/std::max(x,x1);
        const double width=std::min(0.5,std::max(1e-5,dc/std::max(1e-3,b*std::sin(phi))));
        cuts[count++]=std::max(0.0,phi-width);
        cuts[count++]=phi;
        cuts[count++]=std::min(M_PI,phi+width);
    }
    auto integrand=[&](double phi) {
        const double c=std::clamp(a+b*std::cos(phi),-1.0,std::nextafter(1.0,0.0));
        return ktype==0 ? profil_exact_ap(x,x1,c,inverse)
                        : profil_exact(x,x1,c,inverse);
    };
    return 2*kernel_quadrature::integrate_fixed<Order>(integrand,cuts,count);
}

template<int Order> double mean_reference(double x,double x1,double T,int ktype)
{
    const double inverse=511000.0/(8.617333262e-5*T);
    const double peak=1.0-std::abs(x-x1)/(x*x1);
    std::array<double,5> cuts{-1.0,1.0}; size_t count=2;
    if(peak>=-1 && peak<=1) {
        const double dc=4*std::sqrt(2.0/inverse)/std::max(x,x1);
        cuts[count++]=std::max(-1.0,peak-dc);
        cuts[count++]=peak;
        cuts[count++]=std::min(1.0,peak+dc);
    }
    auto integrand=[&](double c) {
        c=std::clamp(c,-1.0,std::nextafter(1.0,0.0));
        return ktype==0 ? profil_exact_ap(x,x1,c,inverse)
                        : profil_exact(x,x1,c,inverse);
    };
    return 2*M_PI*kernel_quadrature::integrate_fixed<Order>(integrand,cuts,count);
}

int main(int argc,char** argv)
{
    const bool report_only=argc==2 && std::strcmp(argv[1],"--report-only")==0;
    struct Case { const char* name; double E, E1, mu, mu1, T; int ktype; };
    const Case cases[] = {
        {"hot_low", 5000, 5200, 0.5, 0.5, 1e9, 1},
        {"hot_high", 200000, 210000, 0.5, 0.5, 1e9, 1},
        {"hot_cross", 300000, 320000, -0.5, 0.5, 1e9, 1},
        {"cold_recoil", 300000, 303000, 0.5, 0.5, 1e4, 1},
        {"warm_upscatter", 300000, 299700, 0.5, 0.5, 1e6, 1},
        {"cold_approx", 300000, 303000, 0.5, 0.5, 1e4, 0},
        {"warm_approx", 300000, 299700, 0.5, 0.5, 1e6, 0},
    };
    for (const auto& c : cases) {
        const double x=c.E/511000.0, x1=c.E1/511000.0;
        const double directional=compton_kernel_element(x,c.mu,x1,c.mu1,c.T,c.ktype);
        const double mean=angle_mean_kernel_element(x,x1,c.T,c.ktype);
        const double d256=directional_reference<256>(x,c.mu,x1,c.mu1,c.T,c.ktype);
        const double d512=directional_reference<512>(x,c.mu,x1,c.mu1,c.T,c.ktype);
        const double m256=mean_reference<256>(x,x1,c.T,c.ktype);
        const double m512=mean_reference<512>(x,x1,c.T,c.ktype);
        const double directional_error=std::abs(directional/d512-1);
        const double mean_error=std::abs(mean/m512-1);
        const double reference_error=std::max(std::abs(d256/d512-1),
                                               std::abs(m256/m512-1));
        std::printf("%s directional %.17g angle_mean %.17g dense_directional %.17g dense_mean %.17g max_rel_error %.3e reference_convergence %.3e\n",
                    c.name,directional,mean,d512,m512,
                    std::max(directional_error,mean_error),reference_error);
        if(!report_only && (!std::isfinite(directional) || !std::isfinite(mean) ||
           directional_error>2e-3 || mean_error>2e-3 || reference_error>1e-4))
            throw std::runtime_error("kernel quadrature differs from converged angular reference");
    }
}
