#include "compton_kernel.h"
#include "avg_compton_kernel.h"
#include "compton_cross_section.h"
#include "source.h"
#include "rt_grids.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>

static void check(bool ok,const char* message) {
    if (!ok) { std::fprintf(stderr,"FAIL: %s\n",message); std::exit(1); }
}

int main() {
    char cache_dir[]="/tmp/dao-kernel-low-temp-XXXXXX";
    check(mkdtemp(cache_dir)!=nullptr,"temporary cache directory");
    setenv("COMPTON_CACHE_DIR",cache_dir,1);
    const double T=1e4, x=300000.0/511000.0;
    const double value=compton_kernel_element(x,0.5,1.01*x,0.5,T,1);
    check(std::isfinite(value) && std::abs(value/4.2635746-1)<2e-3,
          "low-temperature recoil peak is not resolved");
    const double warm_value=compton_kernel_element(100000.0/511000.0,0.5,
                                                   105000.0/511000.0,0.5,1e6,1);
    check(std::isfinite(warm_value) && std::abs(warm_value/8.4064024-1)<2e-3,
          "high-energy warm kernel disagrees with dense angular integral");
    const double up=compton_kernel_element(x,0.5,0.999*x,0.5,1e6,1);
    check(std::isfinite(up) && std::abs(up/0.397281-1)<2e-3,
          "low-temperature high-energy upscattering peak is missed");
    const double up_mean=angle_mean_kernel_element(x,0.999*x,1e6,1);
    check(std::isfinite(up_mean) && std::abs(up_mean/0.0668171-1)<2e-3,
          "angle-mean upscattering peak is missed");
    // Dense 32-point electron and full-order angular references at 1e9 K.
    const double hot_low=compton_kernel_element(5000.0/511000.0,0.5,
                                                5200.0/511000.0,0.5,1e9,1);
    const double hot_high=compton_kernel_element(200000.0/511000.0,0.5,
                                                 210000.0/511000.0,0.5,1e9,1);
    const double hot_cross=compton_kernel_element(300000.0/511000.0,-0.5,
                                                  320000.0/511000.0,0.5,1e9,1);
    check(std::isfinite(hot_low) && std::abs(hot_low/59.044297465938662-1)<2e-3,
          "hot low-energy electron quadrature disagrees with reference");
    check(std::isfinite(hot_high) && std::abs(hot_high/1.3760938574844477-1)<2e-3,
          "hot high-energy quadrature disagrees with reference");
    check(std::isfinite(hot_cross) && std::abs(hot_cross/0.32812526965723343-1)<2e-3,
          "hot cross-angle quadrature disagrees with reference");
    const double ene[]={10000,30000,100000,200000,300000,400000};
    const double mu[]={-0.5,0.5}, wt[]={1.0,1.0};
    KernelCache directional;
    directional.init(6,ene,2,mu,wt,1,1,&T);
    avgKernelCache averaged;
    averaged.init(6,ene,2,mu,wt,1,1,&T);
    for (int i=0; i<6; ++i) {
        check(std::isfinite(directional.K(0,i,1,i,1)),"directional diagonal finite");
        check(std::isfinite(averaged.K(0,i,i)),"angle-mean diagonal finite");
        const double target=compton_cross_section(ene[i],T)/phys::sigma_T;
        double ds=0, as=0;
        for (int j=1; j<6; ++j) {
            const double dx=directional.x_grid[j]-directional.x_grid[j-1];
            auto term=[&](int k) {
                return directional.x_grid[k]*(directional.K(0,k,1,i,1)
                    +directional.K(0,k,1,i,0));
            };
            ds+=0.5*(term(j-1)+term(j))*dx;
            as+=0.5*(averaged.x_grid[j-1]*averaged.K(0,j-1,i)
                    +averaged.x_grid[j]*averaged.K(0,j,i))*dx;
        }
        check(std::abs(ds/directional.x_grid[i]/target-1)<2e-3,
              "directional A23 sum rule");
        check(std::abs(as/averaged.x_grid[i]/target-1)<2e-3,
              "angle-mean A23 sum rule");
        for (int j=i+1; j<6; ++j) {
            const double db=std::exp(-(directional.x_grid[j]-directional.x_grid[i])/
                                      directional.theta[0]);
            const double upper=directional.K(0,i,1,j,1);
            if (upper>0 && db>1e-100)
                check(std::abs(directional.K(0,j,1,i,1)/(upper*db)-1)<1e-10,
                      "directional detailed balance");
        }
    }
    double intensity[12], mean[6], emitted[6]={}, opacity[6]={}, scattering[6];
    double dir_source[12], avg_source[12];
    for (double& v:intensity) v=1.0;
    for (double& v:mean) v=1.0;
    for (int i=0;i<6;++i) scattering[i]=compton_cross_section(ene[i],T)*1.21;
    const double *jptr[]={emitted}, *aptr[]={opacity}, *sptr[]={scattering};
    compute_source_function(1,2,6,intensity,directional.x_grid,wt,directional,
                            &T,jptr,aptr,sptr,1.0,dir_source);
    avgcompute_source_function(1,2,6,mean,averaged.x_grid,averaged,
                                &T,jptr,aptr,sptr,1.0,avg_source);
    bool singleton=false;
    for (int i=0;i<6;++i) if (directional.lo(0,i)==directional.hi(0,i)) {
        singleton=true;
        check(dir_source[6+i]>0,"single-bin directional scattering source");
    }
    check(singleton,"cold kernel contains a single-bin row");
    for (int i=0;i<12;++i) {
        check(std::isfinite(avg_source[i]) && avg_source[i]>=0,
              "finite angle-mean scattering source");
        if (averaged.K(0,i%6,i%6)>0)
            check(avg_source[i]>0,"single-bin angle-mean scattering source");
    }
    // A finer high-energy mesh exercises off-diagonal redistribution and
    // symmetric normalization, unlike the coarse cold mesh above.
    const double warmer=1e6;
    double fine_ene[24];
    for (int i=0; i<24; ++i)
        fine_ene[i]=200000*std::exp(std::log(2.0)*i/23);
    KernelCache fine;
    fine.init(24,fine_ene,2,mu,wt,1,1,&warmer);
    for (int i=0; i<24; ++i) {
        const double target=compton_cross_section(fine_ene[i],warmer)/phys::sigma_T;
        double integral=0;
        for (int j=1; j<24; ++j) {
            auto angular=[&](int k) {
                return fine.x_grid[k]*(fine.K(0,k,1,i,1)+fine.K(0,k,1,i,0));
            };
            integral+=0.5*(angular(j-1)+angular(j))*(fine.x_grid[j]-fine.x_grid[j-1]);
        }
        check(std::abs(integral/fine.x_grid[i]/target-1)<2e-3,
              "fine-grid A23 sum rule");
    }
    const double upper=fine.K(0,10,1,11,1);
    check(upper>0,"fine-grid off-diagonal support");
    const double balance=std::exp(-(fine.x_grid[11]-fine.x_grid[10])/fine.theta[0]);
    check(std::abs(fine.K(0,11,1,10,1)/(upper*balance)-1)<1e-10,
          "fine-grid detailed balance");
    RTGrids angular_grid;
    angular_grid.init_angle();
    KernelCache cold_angles;
    cold_angles.init(24,fine_ene,angular_grid.NA,angular_grid.mu,
                     angular_grid.wt,1,1,&T);
    for (int i=0;i<24;++i) {
        const double target=compton_cross_section(fine_ene[i],T)/phys::sigma_T;
        double integral=0;
        for (int j=1;j<24;++j) {
            auto angular=[&](int k) {
                double sum=0;
                for (int m=0;m<angular_grid.NA;++m) if (angular_grid.mu[m]>0)
                    for (int n=0;n<angular_grid.NA;++n) if (angular_grid.mu[n]>0)
                        sum+=angular_grid.wt[m]*angular_grid.wt[n]*
                             (cold_angles.K(0,k,n,i,m)+
                              cold_angles.K(0,k,n,i,angular_grid.NA-1-m));
                return cold_angles.x_grid[k]*sum;
            };
            integral+=0.5*(angular(j-1)+angular(j))*
                      (cold_angles.x_grid[j]-cold_angles.x_grid[j-1]);
        }
        check(std::abs(integral/cold_angles.x_grid[i]/target-1)<2e-3,
              "eight-angle cold high-energy A23 sum rule");
    }
    double cold_intensity[24*RTGrids::NA], cold_emission[24]={},
           cold_absorption[24]={}, cold_opacity[24], cold_source[24*RTGrids::NA];
    for (double& v:cold_intensity) v=1;
    for (int i=0;i<24;++i)
        cold_opacity[i]=1.21*compton_cross_section(fine_ene[i],T);
    const double *cold_j[]={cold_emission}, *cold_a[]={cold_absorption},
                 *cold_s[]={cold_opacity};
    compute_source_function(1,RTGrids::NA,24,cold_intensity,cold_angles.x_grid,
                            angular_grid.wt,cold_angles,&T,cold_j,cold_a,cold_s,
                            1.0,cold_source);
    for (int i=2;i<24;++i) {
        const double previous=cold_source[(RTGrids::NA-1)*24+i-1];
        const double current=cold_source[(RTGrids::NA-1)*24+i];
        check(std::isfinite(current) && current>0 && current<=1.03*previous,
              "cold high-energy scattering source has zigzags");
    }
    avgKernelCache cold_mean;
    cold_mean.init(24,fine_ene,angular_grid.NA,angular_grid.mu,
                   angular_grid.wt,1,1,&T);
    for(int i=0;i<24;++i) {
        double integral=0;
        for(int j=1;j<24;++j)
            integral+=0.5*(cold_mean.x_grid[j-1]*cold_mean.K(0,j-1,i)
                          +cold_mean.x_grid[j]*cold_mean.K(0,j,i))*
                      (cold_mean.x_grid[j]-cold_mean.x_grid[j-1]);
        const double target=compton_cross_section(fine_ene[i],T)/phys::sigma_T;
        check(std::abs(integral/cold_mean.x_grid[i]/target-1)<2e-3,
              "cold high-energy angle-mean A23 sum rule");
    }
    double flat_mean[24], mean_source[24*RTGrids::NA];
    for(double& v:flat_mean) v=1;
    avgcompute_source_function(1,RTGrids::NA,24,flat_mean,cold_mean.x_grid,
                                cold_mean,&T,cold_j,cold_a,cold_s,1.0,mean_source);
    for(int i=2;i<24;++i) {
        const double previous=mean_source[i-1], current=mean_source[i];
        check(std::isfinite(current) && current>0 && current<=1.03*previous,
              "cold angle-mean scattering source has zigzags");
    }
    double wide_ene[24];
    for(int i=0;i<24;++i)
        wide_ene[i]=10.0*std::pow(1e5,double(i)/23);
    KernelCache wide_cold;
    wide_cold.init(24,wide_ene,angular_grid.NA,angular_grid.mu,
                   angular_grid.wt,1,1,&T);
    for(int i=0;i<24;++i)
        check(std::isfinite(wide_cold.K(0,i,angular_grid.NA-1,i,angular_grid.NA-1)),
              "wide-grid cold diagonal is finite");
    KernelCache reused;
    reused.init(6,ene,2,mu,wt,1,1,&T);
    check(std::abs(reused.K(0,2,1,2,1)/directional.K(0,2,1,2,1)-1)<1e-12,
          "versioned kernel cache reload");
    double shifted[6];
    for (int i=0;i<6;++i) shifted[i]=ene[i];
    shifted[2]*=1.001;
    KernelCache changed_grid;
    changed_grid.init(6,shifted,2,mu,wt,1,1,&T);
    check(std::abs(changed_grid.x_grid[2]-shifted[2]/511000.0)<1e-15,
          "cache energy grid mismatch is rejected");
    changed_grid.free_memory();
    reused.free_memory();
    cold_angles.free_memory();
    cold_mean.free_memory();
    wide_cold.free_memory();
    fine.free_memory();
    averaged.free_memory();
    directional.free_memory();
    std::filesystem::remove_all(cache_dir);
    std::puts("low-temperature kernel checks passed");
}
