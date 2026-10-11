#include "thermal_balance.h"
#include "cell_transfer.h"
#include "compton_cross_section.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <limits>

namespace dao_thermal {

double CellResponse::relative() const
{ return residual/std::max(scale,1e-100); }

CellResponse local_response(int ne,int na,const double* ew,const double* mu,
                            const double* aw,double width,const double* absorption,
                            const double* emission,const Scattering& op,
                            const double* incoming,const double* seed,int max_iterations)
{
    if(ne<2 || na<2 || int(op.opacity().size())!=ne ||
       (op.angles()!=1 && op.angles()!=na) || !std::isfinite(width) || width<=0 || max_iterations<3)
        throw std::runtime_error("local thermal response: invalid grid");
    const bool directional=op.angles()>1;
    const int size=na*ne,field_size=directional ? size : ne;
    std::vector<double> field(seed,seed+field_size),j(field_size),source(size),k(ne);
    std::vector<dao_rt::CellWeights> weights(size);
    CellResponse result;
    result.mean.resize(ne); result.intensity.resize(size);
    for(double v:field) if(!std::isfinite(v) || v<0)
        throw std::runtime_error("local thermal response: invalid seed");
    for(int e=0;e<ne;++e) {
        k[e]=absorption[e]+op.opacity()[e];
        if(!std::isfinite(k[e]) || k[e]<=0 || !std::isfinite(emission[e]) || emission[e]<0)
            throw std::runtime_error("local thermal response: invalid coefficients");
        for(int m=0;m<na;++m) {
            const int i=m*ne+e;
            if(!std::isfinite(incoming[i]) || incoming[i]<0 || !std::isfinite(mu[m]) || mu[m]==0 ||
               !std::isfinite(aw[m]) || aw[m]<=0)
                throw std::runtime_error("local thermal response: invalid incoming ray");
            weights[i]=dao_rt::cell_weights(k[e]*width/std::abs(mu[m]));
        }
    }
    int successive=0;
    for(int iteration=0;iteration<max_iterations;++iteration) {
        op.emissivity(field.data(),j.data());
        std::fill(result.mean.begin(),result.mean.end(),0);
        for(int m=0;m<na;++m) for(int e=0;e<ne;++e) {
            const int i=m*ne+e;
            source[i]=(emission[e]+j[directional ? i : e])/k[e];
            result.intensity[i]=weights[i].mean_incoming*incoming[i]+weights[i].mean_source*source[i];
            if(!std::isfinite(result.intensity[i]) || result.intensity[i]<0)
                throw std::runtime_error("local thermal response: invalid intensity");
            result.mean[e]+=0.5*aw[m]*result.intensity[i];
        }
        const auto& next=directional ? result.intensity : result.mean;
        double change=0;
        for(int i=0;i<field_size;++i)
            change=std::max(change,std::abs(next[i]-field[i])/std::max({next[i],field[i],1e-100}));
        field=next;
        successive=change<1e-9 ? successive+1 : 0;
        if(successive==3) {
            // Equals 4pi integral(k*J-j_material-j_scattering), but remains
            // accurate when I and S almost cancel in an optically thick cell.
            for(int m=0;m<na;++m) for(int e=0;e<ne;++e) {
                const int i=m*ne+e;
                const double weight=0.5*phys::four_pi*aw[m]*ew[e]*k[e]*weights[i].mean_incoming;
                result.residual+=weight*(incoming[i]-source[i]);
                result.scale+=weight*(incoming[i]+source[i]);
            }
            if(!std::isfinite(result.residual) || !std::isfinite(result.scale))
                throw std::runtime_error("local thermal response: invalid energy residual");
            return result;
        }
    }
    throw std::runtime_error("local thermal response: scattering iteration did not converge");
}

std::vector<double> energy_weights(int n, const double* energy)
{
    if(n<2) throw std::runtime_error("thermal balance: energy grid too short");
    std::vector<double> w(n);
    for(int i=0;i<n;++i) {
        if(!std::isfinite(energy[i]) || energy[i]<=0 ||
           (i && energy[i]<=energy[i-1]))
            throw std::runtime_error("thermal balance: invalid energy grid");
        w[i]=0.5*(energy[std::min(n-1,i+1)]-energy[std::max(0,i-1)]);
    }
    return w;
}

// Retain the upper triangle once per independent angular pair. Phase-space
// scaling factors are separate: this avoids expanding the cache to NA^2 for
// every cell while conserving photons for EACH incoming direction.
class DirectionalScattering {
    struct Row { int lo=0,hi=-1; size_t offset=0; };
    int ne_,na_,ni_;
    std::vector<Row> rows_,reverse_rows_;
    std::vector<double> values_,reverse_,scale_,x_,w_,mu_weights_,exchange_;
    std::vector<std::vector<std::pair<int,int>>> angles_;
    template<class F> void pairs(F&& visit) const {
        for(int i=0;i<ne_;++i) for(int a=0;a<ni_;++a) {
            const auto& row=rows_[i*ni_+a];
            const auto& reverse=reverse_rows_[i];
            for(int j=row.lo;j<=row.hi;++j) {
                const double value=values_[row.offset+j-row.lo];
                if(value==0) continue;
                const double back=reverse_[reverse.offset+j-i];
                for(const auto& angle:angles_[a]) {
                    const int m=angle.first,n=angle.second;
                    if(i==j && m>n) continue; // one unordered phase-space pair
                    const int p=m*ne_+i,q=n*ne_+j;
                    const double down=value*scale_[p]*scale_[q];
                    visit(i,j,m,n,p,q,down,down*back);
                }
            }
        }
    }
public:
    std::vector<double> opacity;
    double normalization_error=0;
    int angles() const { return na_; }
    DirectionalScattering(const KernelCache& cache,double T,double electron_density,
                          int na,const double* angular_weights)
        :ne_(cache.NE),na_(na),ni_(cache.n_indep)
    {
        if(cache.NT<1 || ne_<2 || na!=cache.NA_full || na<2 || !angular_weights ||
           !std::isfinite(T) || T<cache.T_grid[0] || T>cache.T_grid[cache.NT-1])
            throw std::runtime_error("directional thermal balance: invalid grid or temperature outside kernel table");
        if(!std::isfinite(electron_density) || electron_density<0)
            throw std::runtime_error("directional thermal balance: invalid free-electron density");
        constexpr double mec2=511000.0; // same units as the retained cache
        const double theta=phys::k_B*T/(mec2*phys::eV_to_erg);
        x_.assign(cache.x_grid,cache.x_grid+ne_);w_=energy_weights(ne_,x_.data());
        mu_weights_.assign(angular_weights,angular_weights+na_);
        double sum_weights=0;
        for(double w:mu_weights_) {
            if(!std::isfinite(w) || w<=0) throw std::runtime_error("invalid angular weights");
            sum_weights+=w;
        }
        if(std::abs(sum_weights-2)>1e-12) throw std::runtime_error("angular weights must integrate dmu on [-1,1]");
        const int high=int(std::lower_bound(cache.T_grid,cache.T_grid+cache.NT,T)-cache.T_grid);
        const int low=std::max(0,high-1);
        const double fraction=low==high ? 0 : std::log(T/cache.T_grid[low])/std::log(cache.T_grid[high]/cache.T_grid[low]);
        angles_.resize(ni_);
        for(int m=0;m<na_;++m) for(int n=0;n<na_;++n)
            angles_[cache.canon[m*na_+n]].push_back({m,n});
        rows_.resize(ne_*ni_);reverse_rows_.resize(ne_);
        size_t count=0,reverse_count=0;
        for(int i=0;i<ne_;++i) {
            int end=i;
            for(int a=0;a<ni_;++a) {
                const auto l=cache.row(low,i,a),h=cache.row(high,i,a);
                const int lo=std::max(i,std::min(cache.band_lo[l],cache.band_lo[h]));
                const int hi=std::max(cache.band_hi[l],cache.band_hi[h]);
                rows_[i*ni_+a]={lo,hi,count};
                count+=size_t(std::max(0,hi-lo+1));end=std::max(end,hi);
            }
            reverse_rows_[i]={i,end,reverse_count};reverse_count+=end-i+1;
        }
        values_.resize(count);reverse_.resize(reverse_count);
        for(int i=0;i<ne_;++i) {
            const auto& r=reverse_rows_[i];
            for(int j=i;j<=r.hi;++j) reverse_[r.offset+j-i]=std::exp(-(x_[j]-x_[i])/theta);
            for(int a=0;a<ni_;++a) {
                const auto& row=rows_[i*ni_+a];
                const int m=cache.indep_nm[a],n=cache.indep_nm1[a];
                for(int j=row.lo;j<=row.hi;++j) {
                    const double value=(1-fraction)*cache.K(low,i,m,j,n)+fraction*cache.K(high,i,m,j,n);
                    if(!std::isfinite(value) || value<0) throw std::runtime_error("invalid directional kernel");
                    values_[row.offset+j-row.lo]=value;
                }
            }
        }
        std::vector<double> target(ne_),sum(ne_*na_);
        scale_.assign(ne_*na_,1);opacity.resize(ne_);
        for(int i=0;i<ne_;++i) target[i]=compton_cross_section(x_[i]*mec2,T)/phys::sigma_T;
        bool converged=false;
        for(int iteration=0;iteration<200;++iteration) {
            std::fill(sum.begin(),sum.end(),0);
            pairs([&](int i,int j,int m,int n,int p,int q,double down,double up) {
                sum[q]+=w_[i]*mu_weights_[m]*x_[i]/x_[j]*down;
                if(p!=q) sum[p]+=w_[j]*mu_weights_[n]*x_[j]/x_[i]*up;
            });
            normalization_error=0;
            for(int p=0;p<ne_*na_;++p) {
                const double ratio=target[p%ne_]/sum[p];
                if(!std::isfinite(ratio) || ratio<=0) throw std::runtime_error("unsupported directional kernel column");
                normalization_error=std::max(normalization_error,std::abs(sum[p]/target[p%ne_]-1));
            }
            if(normalization_error<2e-10) { converged=true;break; }
            // Symmetric scaling in (energy, direction) preserves reciprocity.
            for(int p=0;p<ne_*na_;++p) scale_[p]*=std::sqrt(target[p%ne_]/sum[p]);
        }
        if(!converged) throw std::runtime_error("directional thermal kernel normalization did not converge");
        const double density=electron_density*phys::sigma_T;
        exchange_.resize(ne_*na_);
        for(int i=0;i<ne_;++i) {
            opacity[i]=density*target[i];
            for(int m=0;m<na_;++m) exchange_[m*ne_+i]=0.5*mu_weights_[m]*w_[i]*mec2*opacity[i];
            // Convert retained K to emission coefficients, excluding angular
            // weights. Reverse factors depend on energy only, not angular pair.
            const auto& r=reverse_rows_[i];
            for(int j=i;j<=r.hi;++j) {
                const double ratio=x_[j]/x_[i];
                reverse_[r.offset+j-i]*=w_[i]/w_[j]*ratio*ratio*ratio*ratio;
            }
            for(int a=0;a<ni_;++a) {
                const auto& row=rows_[i*ni_+a];
                for(int j=row.lo;j<=row.hi;++j)
                    values_[row.offset+j-row.lo]*=density*w_[j]*x_[i]*x_[i]/(x_[j]*x_[j]);
            }
        }
        pairs([&](int i,int j,int m,int n,int p,int q,double down,double up) {
            exchange_[q]-=0.5*mu_weights_[m]*mu_weights_[n]*w_[i]*mec2*down;
            if(p!=q) exchange_[p]-=0.5*mu_weights_[n]*mu_weights_[m]*w_[j]*mec2*up;
        });
    }
    void emissivity(const double* intensity,double* emission) const {
        std::fill(emission,emission+ne_*na_,0);
        pairs([&](int,int,int m,int n,int p,int q,double down,double up) {
            emission[p]+=down*mu_weights_[n]*intensity[q];
            if(p!=q) emission[q]+=up*mu_weights_[m]*intensity[p];
        });
    }
    double energy_gain(const double* intensity) const {
        double q=0;
        for(size_t p=0;p<exchange_.size();++p) q+=exchange_[p]*intensity[p];
        return phys::four_pi*q;
    }
};

Scattering::Scattering(const KernelCache& cache,double T,double electron_density,int na,const double* weights)
{
    auto op=std::make_shared<DirectionalScattering>(cache,T,electron_density,na,weights);
    opacity_=op->opacity;normalization_error=op->normalization_error;
    directional_=std::move(op);
}

int Scattering::angles() const { return directional_ ? directional_->angles() : 1; }

Scattering::Scattering(const avgKernelCache& cache, double T, double electron_density)
{
    if(cache.NT<1 || !std::isfinite(T) || T<cache.T_grid[0] ||
       T>cache.T_grid[cache.NT-1])
        throw std::runtime_error("thermal balance: temperature outside kernel table (no extrapolation)");
    if(!std::isfinite(electron_density) || electron_density<0)
        throw std::runtime_error("thermal balance: invalid free-electron density");
    const int n=cache.NE;
    const double* x=cache.x_grid;
    auto w=energy_weights(n,x);
    // avgKernelCache stores x=E/511000 eV. Use that SAME conversion here;
    // phys::m_e_c2_eV is 510999 eV and would break the discrete energy identity.
    constexpr double cache_mec2_eV=511000.0;
    const double theta=phys::k_B*T/(cache_mec2_eV*phys::eV_to_erg);
    int high=int(std::lower_bound(cache.T_grid,cache.T_grid+cache.NT,T)-cache.T_grid);
    int low=std::max(0,high-1);
    double f=low==high ? 0 : std::log(T/cache.T_grid[low])/std::log(cache.T_grid[high]/cache.T_grid[low]);
    for(int i=0;i<n;++i) {
        const int hi=std::max(cache.band_hi[cache.row(low,i)],cache.band_hi[cache.row(high,i)]);
        for(int j=i;j<=hi;++j) {
            const double a=(1-f)*cache.K(low,i,j)+f*cache.K(high,i,j);
            if(!std::isfinite(a) || a<0) throw std::runtime_error("thermal balance: invalid kernel");
            if(a>0) pairs.push_back({i,j,a,a*std::exp(-(x[j]-x[i])/theta)});
        }
    }
    std::vector<double> target(n),norm(n),sum(n);
    opacity_.resize(n); exchange_.assign(n,0);
    for(int i=0;i<n;++i)
        target[i]=compton_cross_section(x[i]*cache_mec2_eV,T)/phys::sigma_T;
    bool converged=false;
    for(int iteration=0;iteration<200;++iteration) {
        std::fill(sum.begin(),sum.end(),0);
        for(const auto& p:pairs) {
            sum[p.j]+=w[p.i]*x[p.i]/x[p.j]*p.down;
            if(p.i!=p.j) sum[p.i]+=w[p.j]*x[p.j]/x[p.i]*p.up;
        }
        normalization_error=0;
        for(int i=0;i<n;++i) {
            norm[i]=target[i]/sum[i];
            if(!std::isfinite(norm[i]) || norm[i]<=0)
                throw std::runtime_error("thermal balance: unsupported kernel column");
            normalization_error=std::max(normalization_error,std::abs(sum[i]/target[i]-1));
        }
        if(normalization_error<2e-10) { converged=true; break; }
        for(auto& p:pairs) {
            const double s=std::sqrt(norm[p.i]*norm[p.j]);
            p.down*=s; p.up*=s;
        }
    }
    if(!converged) throw std::runtime_error("thermal balance: kernel normalization did not converge");
    const double density=electron_density*phys::sigma_T;
    for(int i=0;i<n;++i) {
        opacity_[i]=density*target[i];
        exchange_[i]=w[i]*cache_mec2_eV*opacity_[i];
    }
    // Store coefficients of the actual emission operator (not S=j/kappa).
    for(auto& p:pairs) {
        p.down*=density*w[p.j]*x[p.i]*x[p.i]/(x[p.j]*x[p.j]);
        p.up*=density*w[p.i]*x[p.j]*x[p.j]/(x[p.i]*x[p.i]);
        exchange_[p.j]-=w[p.i]*cache_mec2_eV*p.down;
        if(p.i!=p.j) exchange_[p.i]-=w[p.j]*cache_mec2_eV*p.up;
    }
}

void Scattering::emissivity(const double* J, double* j) const
{
    if(directional_) { directional_->emissivity(J,j);return; }
    std::fill(j,j+opacity_.size(),0);
    for(const auto& p:pairs) {
        j[p.i]+=p.down*J[p.j];
        if(p.i!=p.j) j[p.j]+=p.up*J[p.i];
    }
}

double Scattering::energy_gain(const double* J) const
{
    if(directional_) return directional_->energy_gain(J);
    double q=0;
    for(size_t i=0;i<exchange_.size();++i) q+=exchange_[i]*J[i];
    return phys::four_pi*q;
}

double Budget::scale() const { return absorption+emission+std::abs(compton); }
double Budget::relative() const { return residual()/std::max(scale(),1e-100); }

Budget budget(int n, const double* w, const double* J,
              const double* a, const double* j, const Scattering& s, const double* intensity)
{
    Budget b;
    for(int i=0;i<n;++i) {
        if(!std::isfinite(J[i]) || !std::isfinite(a[i]) || !std::isfinite(j[i]) ||
           J[i]<0 || a[i]<0 || j[i]<0)
            throw std::runtime_error("thermal balance: invalid local radiation coefficients");
        b.absorption+=phys::four_pi*w[i]*a[i]*J[i];
        b.emission+=phys::four_pi*w[i]*j[i];
    }
    if(s.angles()>1 && !intensity)
        throw std::runtime_error("directional thermal budget requires the full specific intensity");
    b.compton=s.energy_gain(s.angles()>1 ? intensity : J);
    if(!std::isfinite(b.residual()) || !std::isfinite(b.scale()))
        throw std::runtime_error("thermal balance: nonfinite residual");
    return b;
}

TemperatureRoot::TemperatureRoot(double t, double lo, double hi)
    :initial_(std::log(t)),minimum_(std::log(lo)),maximum_(std::log(hi))
{
    if(!std::isfinite(initial_) || lo<=0 || hi<=lo || t<lo || t>hi)
        throw std::runtime_error("thermal balance: invalid temperature bracket bounds");
}

// record a new sample of the residual of the thermal balance equation at a given temperature.
void TemperatureRoot::record(double t, double r)
{
    if(!std::isfinite(t) || !std::isfinite(r) || t<=0)
        throw std::runtime_error("thermal balance: nonfinite root sample");
    samples_[std::log(t)]=r;
}

// Brackets a stable root (positive heating below, negative heating above) in
double TemperatureRoot::next()
{
    if(samples_.empty()) return std::exp(initial_);
    if(samples_.size()>=60) throw std::runtime_error("thermal balance: temperature root iteration limit");
    for(auto it=samples_.begin();it!=samples_.end();++it) {
        auto hi=std::next(it);
        if(hi==samples_.end()) break;
        if(it->second>0 && hi->second<0) {
            const double a=it->first,b=hi->first;
            double t=a+(b-a)*it->second/(it->second-hi->second);
            if(t<a+0.1*(b-a) || t>b-0.1*(b-a)) t=0.5*(a+b);
            if(b-a<1e-10) throw std::runtime_error("thermal balance: discontinuous/unresolved temperature root");
            return std::exp(t);
        }
    }
    // Search both sides of the original temperature, preferring the direction
    // implied by the first residual. No assumption that the whole curve is monotone.
    auto origin=std::min_element(samples_.begin(),samples_.end(),[&](const auto& a,const auto& b) {
        return std::abs(a.first-initial_)<std::abs(b.first-initial_);
    });
    const double sign=origin->second>=0 ? 1 : -1;
    while(expansion_<100) {
        const int k=expansion_++;
        const double step=(k/2+1)*std::log(1.5)*(k%2 ? -sign : sign);
        const double t=std::clamp(initial_+step,minimum_,maximum_);
        bool seen=false;
        for(const auto& sample:samples_) if(std::abs(sample.first-t)<1e-12) seen=true;
        if(!seen) return std::exp(t);
    }
    throw std::runtime_error("thermal balance: no stable temperature root within kernel table");
}
} // namespace dao_thermal
