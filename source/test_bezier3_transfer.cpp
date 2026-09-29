#include "bezier3_transfer.h"
#include <cstdio>
#include <numeric>

static void check(bool ok, const char* message)
{
	if (!ok) throw std::runtime_error(message);
}

static void analytic()
{
	const double mu[]={-0.3,0.7}, up=1.7, bot=0.6;
	for (double scale : {0.0,1e-12,1.0,1e4}) for (int nd : {1,4})
	for (double s : {0.0,2.0}) {
		std::vector<double> dt(nd*2), source(nd*2,s), I(nd*2), mean(nd*2), sm(nd*2);
		double top[2], bottom[2];
		for (int d=0;d<nd;++d) for (int m=0;m<2;++m)
			dt[d*2+m]=scale*(0.1+d*d)/std::abs(mu[m]);
		dao_rt::Bezier3FormalSolver solver(nd,2,1,dt.data());
		solver.solve(-1,mu,&up,&bot,source.data(),I.data(),top,bottom,mean.data(),sm.data());
		for (int m=0;m<2;++m) {
			double tau=0, incoming=m==0 ? up : bot;
			for (int step=0;step<nd;++step) {
				int d=m==0 ? step : nd-1-step, k=d*2+m;
				double expected=s+(incoming-s)*std::exp(-tau-dt[k]/2);
				check(std::abs(I[k]-expected)<2e-12,"analytic center intensity");
				double face=s+(incoming-s)*std::exp(-tau);
				auto w=dao_rt::cell_weights(dt[k]);
				check(std::abs(mean[k]-(w.mean_incoming*face+w.mean_source*s))<2e-12,"analytic cell average");
				tau+=dt[k];
			}
			check(std::abs((m==0?bottom[m]:top[m])-(s+(incoming-s)*std::exp(-tau)))<2e-12,"analytic boundary intensity");
		}
	}
	std::puts("PASS analytic centers, true faces, cell averages, vacuum and thin/thick limits");
}

static void weights()
{
	const double controls[]={0.2,2.1,0.4,1.3};
	for (double t : {1e-10,0.1,1.99,2.0,3.0,20.0,100.0}) {
		auto w=dao_rt::bezier3_weights(t);
		double integral=0, outgoing=0, mean=0;
		const int n=20000;
		for (int i=0;i<=n;++i) {
			double u=double(i)/n, v=1-u;
			double s=controls[0]*v*v*v+3*controls[1]*u*v*v+3*controls[2]*u*u*v+controls[3]*u*u*u;
			integral+=(i==0 || i==n ? 1 : (i%2 ? 4 : 2))*t*std::exp(-t*v)*s/(3*n);
		}
		for (int j=0;j<4;++j) {
			check(w.outgoing[j]>=0 && w.mean[j]>=0,"negative Bernstein weight");
			outgoing+=w.outgoing[j]*controls[j];mean+=w.mean[j]*controls[j];
		}
		check(std::abs(outgoing-integral)<1e-10,"independent cubic quadrature");
		check(std::abs(outgoing-t*(1.0-mean))<1e-12,"cubic volume balance");
	}
	std::puts("PASS independent cubic quadrature and formal volume balance");
}

static double smooth(int nd)
{
	const double mu[]={-1.0}, up=0.4, bot=0;
	std::vector<double> dt(nd,2.0/nd), source(nd), I(nd), mean(nd), sm(nd);
	for (int d=0;d<nd;++d) {double t=(d+0.5)*dt[d];source[d]=1+t+t*t;}
	double top,bottom;
	dao_rt::Bezier3FormalSolver solver(nd,1,1,dt.data());
	solver.solve(-1,mu,&up,&bot,source.data(),I.data(),&top,&bottom,mean.data(),sm.data());
	double volume=0;
	for (int d=0;d<nd;++d) {
		check(std::isfinite(I[d]) && I[d]>=0,"nonfinite/negative smooth intensity");
		volume+=dt[d]*(sm[d]-mean[d]);
	}
	check(std::abs(bottom-up-volume)<1e-12,"split cubic volume balance");
	double previous=bottom;
	solver.solve(-1,mu,&up,&bot,source.data(),I.data(),&top,&bottom);
	check(previous==bottom,"diagnostic averages altered formal solution");
	// Exact solution for S(t)=1+t+t^2: I(t)=t^2-t+2+(I(0)-2)exp(-t).
	return std::abs(bottom-(4+(up-2)*std::exp(-2.0)));
}

int main()
{
	try {
		analytic();weights();
		double coarse=smooth(20), medium=smooth(40), fine=smooth(80);
		check(medium<coarse/3 && fine<medium/3 && fine<5e-4,"smooth-source depth convergence");
		std::printf("PASS smooth-source refinement errors: %.6g %.6g %.6g\n",coarse,medium,fine);
		return 0;
	} catch (const std::exception& e) {std::fprintf(stderr,"FAIL: %s\n",e.what());return 1;}
}
