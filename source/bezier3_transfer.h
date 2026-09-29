#ifndef DAO_BEZIER3_TRANSFER_H
#define DAO_BEZIER3_TRANSFER_H

#include "cell_transfer.h"
#include <algorithm>
#include <array>
#include <vector>

namespace dao_rt {

// Exact integrals of the four cubic Bernstein basis functions, ordered
// from the incoming to the outgoing end of a ray segment.
struct Bezier3Weights {
	double attenuation, mean_incoming;
	std::array<double,4> outgoing, mean;
};

inline Bezier3Weights bezier3_weights(double t)
{
	const auto constant = cell_weights(t); // also validates optical depth
	Bezier3Weights w{};
	w.attenuation = constant.attenuation;
	w.mean_incoming = constant.mean_incoming;
	if (t < 2.0) {
		// Evaluate complementary mean weights directly, avoiding cancellation
		// for thin segments. The convergent series is also valid at t=0.
		for (int j=0; j<4; ++j) {
			double term=0.25, complement=0.0;
			for (int n=1; n<=40; ++n) {
				term *= -t*(3-j+n)/(double(n)*(4+n));
				complement -= term;
			}
			w.mean[j]=complement;
			w.outgoing[j]=t*(0.25-complement);
		}
	} else {
		const double v=1/t, v2=v*v, v3=v2*v, a=w.attenuation;
		w.outgoing = {{6*v3-a*(1+3*v+6*v2+6*v3),
			6*v2-18*v3+a*(3*v+12*v2+18*v3),
			3*v-12*v2+18*v3-a*(6*v2+18*v3),
			1-3*v+6*v2-6*v3+6*v3*a}};
		for (int j=0; j<4; ++j) w.mean[j]=0.25-w.outgoing[j]/t;
	}
	return w;
}

inline double bezier3_slope(double left, double center, double right,
	                          double h_left, double h_right)
{
	if (h_left==0 || h_right==0) return 0;
	const double a=(center-left)/h_left, b=(right-center)/h_right;
	if (a==0 || b==0 || std::signbit(a)!=std::signbit(b)) return 0;
	const double weight=(1+h_right/(h_left+h_right))/3;
	// Fritsch-Butland weighted harmonic mean, without the product a*b.
	return 1/(weight/a+(1-weight)/b);
}

// Point-value cubic Bezier short characteristics on physical cell centers.
// optical_depth holds FULL cell depths [depth][angle][energy] and must remain
// alive and unchanged for this object's lifetime. Material coefficients are
// constant per cell; center-to-center paths contain TWO adjacent half cells.
class Bezier3FormalSolver {
	int nd_, nm_, ne_;
	const double* optical_depth_;
	std::vector<Bezier3Weights> half_;
	std::vector<double> derivative_;

public:
	Bezier3FormalSolver(int nd, int nm, int ne, const double* optical_depth)
		: nd_(nd), nm_(nm), ne_(ne), optical_depth_(optical_depth)
	{
		if (nd<1 || nm<1 || ne<1)
			throw std::runtime_error("Bezier3: empty grid");
		half_.resize(long(nd)*nm*ne);
		derivative_.resize(nd);
		for (size_t k=0;k<half_.size();++k)
			half_[k]=bezier3_weights(optical_depth_[k]/2);
	}

	size_t memory_bytes() const
	{
		return half_.size()*sizeof(Bezier3Weights)+derivative_.size()*sizeof(double);
	}

	// intensity is evaluated at centers, NOT averaged over cells. Surface
	// arrays are [angle][energy], including prescribed incoming hemispheres.
	// Optional averages integrate the SAME interpolated source and intensity
	// over each full physical cell; they are useful for transport diagnostics.
	// Point-value source iteration need not conserve the physical cell source
	// integral exactly on a coarse grid: depth convergence must be checked.
	void solve(int incidence, const double* mu,
	           const double* top, const double* bottom, const double* source,
	           double* intensity, double* surface_top, double* surface_bottom,
	           double* cell_intensity=nullptr, double* cell_source=nullptr)
	{
		if (cell_intensity) std::fill(cell_intensity,cell_intensity+half_.size(),0.0);
		if (cell_source) std::fill(cell_source,cell_source+half_.size(),0.0);
		for (int m=0;m<nm_;++m) {
			if (!std::isfinite(mu[m]) || mu[m]==0)
				throw std::runtime_error("Bezier3: invalid ray cosine");
			for (int e=0;e<ne_;++e) {
				auto idx=[&](int d){return (long(d)*nm_+m)*ne_+e;};
				for (int d=0;d<nd_;++d)
					derivative_[d]=(d==0 || d==nd_-1) ? 0 : bezier3_slope(
						source[idx(d-1)],source[idx(d)],source[idx(d+1)],
						(optical_depth_[idx(d-1)]+optical_depth_[idx(d)])/2,
						(optical_depth_[idx(d)]+optical_depth_[idx(d+1)])/2);
				const bool down=mu[m]<0;
				const int first=down ? 0 : nd_-1, sign=down ? 1 : -1;
				const long ray=long(m)*ne_+e;
				double incoming=down ? incident_top_intensity(incidence,m,top[e]) : bottom[e];
				(down ? surface_top : surface_bottom)[ray]=incoming;
				auto integrate=[&](int d,const std::array<double,4>& control) {
					const long k=idx(d);
					const auto& w=half_[k];
					double outgoing=w.attenuation*incoming;
					for (int j=0;j<4;++j) outgoing+=w.outgoing[j]*control[j];
					if (cell_intensity) {
						double mean=w.mean_incoming*incoming;
						for (int j=0;j<4;++j) mean+=w.mean[j]*control[j];
						cell_intensity[k]+=0.5*mean;
					}
					if (cell_source)
						for (double value:control) cell_source[k]+=value/8;
					incoming=outgoing;
				};
				// Actual slab face to first center: use the boundary cell source.
				double s=source[idx(first)];
				integrate(first,{{s,s,s,s}});
				intensity[idx(first)]=incoming;
				for (int step=1;step<nd_;++step) {
					const int d=first+sign*step, u=d-sign;
					const double h=(optical_depth_[idx(d)]+optical_depth_[idx(u)])/2;
					const double su=source[idx(u)], sd=source[idx(d)];
					const double low=std::min(su,sd), high=std::max(su,sd);
					const double c1=std::clamp(su+sign*h*derivative_[u]/3,low,high);
					const double c2=std::clamp(sd-sign*h*derivative_[d]/3,low,high);
					const double q=h>0 ? optical_depth_[idx(u)]/(2*h) : 0.5;
					// Split the cubic at the physical face, which is generally
					// NOT halfway in optical depth when extinction changes.
					auto mix=[&](double a,double b){return (1-q)*a+q*b;};
					const double p01=mix(su,c1), p12=mix(c1,c2), p23=mix(c2,sd);
					const double p012=mix(p01,p12), p123=mix(p12,p23), p=mix(p012,p123);
					integrate(u,{{su,p01,p012,p}});
					integrate(d,{{p,p123,p23,sd}});
					intensity[idx(d)]=incoming;
				}
				// Last center to actual slab face: include the final half cell.
				const int last=first+sign*(nd_-1);
				s=source[idx(last)];
				integrate(last,{{s,s,s,s}});
				(down ? surface_bottom : surface_top)[ray]=incoming;
			}
		}
	}
};

} // namespace dao_rt
#endif
