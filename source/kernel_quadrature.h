#ifndef DAO_KERNEL_QUADRATURE_H
#define DAO_KERNEL_QUADRATURE_H

#include <algorithm>
#include <array>
#include <cmath>

// The caller supplies breakpoints at the minimum-electron-energy angle and
// the edges of the narrow thermal peak. Production uses a fixed high-order
// rule on those intervals; adaptive quadrature remains available as an
// independent convergence reference.
namespace kernel_quadrature {

template<int N> struct GaussRule {
    double node[N], weight[N];
    GaussRule() {
        for (int i=0; i<N/2; ++i) {
            double z=std::cos(M_PI*(i+0.75)/(N+0.5)), p=0, dp=0;
            for (int k=0; k<30; ++k) {
                double p0=1, p1=z;
                for (int j=2; j<=N; ++j) {
                    double next=((2*j-1)*z*p1-(j-1)*p0)/j;
                    p0=p1; p1=next;
                }
                p=p1;
                dp=N*(z*p1-p0)/(z*z-1);
                double dz=p/dp;
                z-=dz;
                if (std::abs(dz)<1e-15) break;
            }
            double w=2/((1-z*z)*dp*dp);
            node[i]=-z; node[N-1-i]=z;
            weight[i]=w; weight[N-1-i]=w;
        }
    }
};

template<int N, class F> double gauss(F& f, double a, double b) {
    static const GaussRule<N> rule;
    const double mid=(a+b)*0.5, half=(b-a)*0.5;
    double sum=0;
    for (int i=0; i<N; ++i) sum+=rule.weight[i]*f(mid+half*rule.node[i]);
    return half*sum;
}

template<class F> double adaptive(F& f, double a, double b, int depth=0) {
    const double coarse=gauss<8>(f,a,b);
    const double fine=gauss<16>(f,a,b);
    if (depth>=16 || std::abs(fine-coarse)<=1e-4*std::abs(fine))
        return fine;
    const double mid=(a+b)*0.5;
    return adaptive(f,a,mid,depth+1)+adaptive(f,mid,b,depth+1);
}

template<class F, size_t N> double integrate(F& f, std::array<double,N> cuts,
                                             size_t count=N) {
    std::sort(cuts.begin(),cuts.begin()+count);
    double sum=0;
    for (size_t i=1; i<count; ++i)
        if (cuts[i]>cuts[i-1]) sum+=adaptive(f,cuts[i-1],cuts[i]);
    return sum;
}

template<int Order, class F, size_t N> double integrate_fixed(F& f,
                    std::array<double,N> cuts,size_t count=N) {
    std::sort(cuts.begin(),cuts.begin()+count);
    double sum=0;
    for(size_t i=1;i<count;++i)
        if(cuts[i]>cuts[i-1]) sum+=gauss<Order>(f,cuts[i-1],cuts[i]);
    return sum;
}

} // namespace kernel_quadrature
#endif
