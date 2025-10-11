!
!> @brief Computes convolution due to thermal and turbulent motion
!! @details  this routine Computes convolution due to thermal 
!!     and turbulent motion.  Called by gsmooth.
!!     author:  T. Kallman 
!!
!! @author T. Kallman
!!
!! @param[in] lpri print switch
!! @param[in] lun11 logical unit for printing
!! @param[in] xee electron fraction relative to H
!! @param[in] xpx H number density (cm^-3)
!! @param[in] t temperature in 10^4K
!! @param[in] epi(ncn) photon energy grid (ev)
!! @param[in] ncn2 length of epi
!! @param[out] brcems(ncn)  bremsstrahlung emissivities (erg cm^-3 s^-1 erg^-1)/10^38
!! @param[out] opakc(ncn)  continuum opacities with lines binned in (cm^-1)
      subroutine gsmooth2(lpri,lun11,vtherm,epi,ncn2,ydat)
!
!     Name: gsmooth2.f90  
!     Description:  
!     this routine smooths emissivities and opaicities using gaussian
!     author:  T. Kallman (12/1/2024)                                 
!
!     List of Parameters:
!           Input:
!           lpri=print switch
!           lun11=logical unit for printing
!           xee: electron fraction relative to H
!           xpx: H number density (cm^-3)
!           t: temperature in 10^4K
!           epi(ncn): photon energy grid (ev)
!           ncn2: length of epi
!           brcems:  brems emissivity
!           rccems:  recomb emissivity
!           opakc(ncn):  continuum opacities with lines binned in (cm^-1)
!           Output:
!           opakc(ncn):  continuum opacities with lines binned in (cm^-1)
!     Dependencies:  none
!     Called by:  calc_hmc_all

!     author:  T. Kallman                                               
!
      use globaldata
      use constants
      implicit none 
!                                                                       
!                                                                       
      real(8) epi(ncn),ydat(ncn),xpx,t,xee,prod(ncn)
      integer lpri,lun11,ncn2,numcon,kl,kk,nb1,nb2,kl2,nbinc,kl2m,kl2p,klm
      real(8) temp,gam,gau,dele,vtherm,sum1m,sum2m,tmp1p,tmp2p,rkern,     &
     &        tmp1po,tmp2po,earg,ydatsum1,ydatsum2,tmp1m,tmp2m,           &
     &        tmp1mo,tmp2mo,sum1p,sum2p
      logical done1,done2
    
      integer lprisv 
!
      lprisv=lpri
      lpri=0
!
      ydatsum1=0.
      do kl = 1,ncn2
        if (kl.gt.1) then
        ydatsum1=ydatsum1+(ydat(kl)+ydat(kl-1))*(epi(kl)-epi(kl-1))/2.
        endif
        enddo

      do kl = 3,ncn2
!
         lpri=0
!         if (kl.eq.299) lpri=1
!
!        standard quantities
         prod(kl)=ydat(kl)
         dele=epi(kl)*vtherm/3.e+10
!
         if (epi(kl).lt.2.e+4) then
!
!        sum in two directions
         sum1p=0.
         sum2p=0.
         sum1m=0.
         sum2m=0.
         tmp1p=0.
         tmp2p=0.
         tmp2m=0.
         tmp1m=0.
         done1=.false.
         done2=.false.
         klm=1
         do while ((.not.done1).or.(.not.done2))
!
!          plus
           if (.not.done2) then
             kl2p=kl+klm
             if (kl2p.ge.ncn2-1) done2=.true.
             tmp1po=tmp1p
             tmp2po=tmp2p
             earg=((epi(kl2p)-epi(kl))/dele)**2
             if (earg.ge.30.) done2=.true.
             rkern=exp(-earg)
             tmp1p=rkern*ydat(kl2p)
             tmp2p=rkern
             sum1p=sum1p+(tmp1p+tmp1po)*(epi(kl2p)-epi(kl2p-1))/2.
             sum2p=sum2p+(tmp2p+tmp2po)*(epi(kl2p)-epi(kl2p-1))/2.
             if (lpri.ne.0) write (lun11,*)kl,klm,kl2p,epi(kl2p),       &
     &         epi(kl),earg,rkern,tmp1p,tmp2p,sum1p,sum2p
             endif
!
!          minus
           if (.not.done1) then
             kl2m=kl-klm
             if (kl2m.le.3) done1=.true.
             tmp1mo=tmp1m
             tmp2mo=tmp2m
             earg=((epi(kl2m)-epi(kl))/dele)**2
             if (earg.ge.30.) done1=.true.
             rkern=exp(-earg)
             tmp1m=rkern*ydat(kl2m)
             tmp2m=rkern
             sum1m=sum1m+(tmp1m+tmp1mo)*abs(epi(kl2m)-epi(kl2m+1))/2.
             sum2m=sum2m+(tmp2m+tmp2mo)*abs(epi(kl2m)-epi(kl2m+1))/2.
             if (lpri.ne.0) write (lun11,*)kl,klm,kl2p,epi(kl2p),       &
     &         epi(kl),earg,rkern,tmp1p,tmp2p,sum1p,sum2p
             endif
!
           klm=klm+1
!
           enddo
        
         prod(kl)=(sum1p+sum1m)/max(1.e-38,sum2p+sum2m)
!
         if (lpri.ne.0) write (lun11,*)'prod(kl)=',kl,sum1p,sum1m,      &
     &        sum2p,sum2m,prod(kl)
!
         endif
!
         enddo 
!
!
      ydatsum2=0.
      do kl = 3,ncn2
        ydat(kl)=prod(kl)
        if (kl.gt.1) then
        ydatsum2=ydatsum2+(ydat(kl)+ydat(kl-1))*(epi(kl)-epi(kl-1))/2.
        endif
        enddo
!
!      write (lun11,*)'sum1=',ydatsum1,'sum2=',ydatsum2

!                                                                       
      lpri=lprisv 
!                                                                       
      return 
99001 format (' ',i6,9e12.4) 
      END                                           
