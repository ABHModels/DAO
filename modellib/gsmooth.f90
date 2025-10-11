!> @brief Computes convolution due to thermal and turbulent motion
!! @details  this routine Computes convolution due to thermal 
!!     and turbulent motion
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
      subroutine gsmooth(lpri,lun11,t,vturbi,epi,ncn2,opakc,rccemis,brcems) 
!
!     Name: gsmooth.f90  
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
!     Called by:  xstar

!     author:  T. Kallman                                               
!                                                                       
      use globaldata
      use constants
      implicit none 
!                                                                       
!                                                                       
      real(8) epi(ncn),brcems(ncn),xpx,t,xee,opakc(ncn),rccemis(2,ncn)
      integer lpri,lun11,ncn2,numcon,kl,kk 
      real(8) vturbi,vtherm,rccems1(ncn),rccems2(ncn)
      integer lskp,lprisv,mm
      real(8) bbee 
!                                                                       
      lprisv=lpri 
      if (lpri.gt.1) write (lun11,*)'in gsmooth',t 
!                    
      vtherm=((vturbi*1.e+5)**2+(1.29e+6/sqrt(1/t))**2)**(0.5) 
!                                                   
!      write (lun11,*)'call gsmooth on brcems'
      call gsmooth2(lpri,lun11,vtherm,epi,ncn2,brcems)
      do mm=1,ncn2
        rccems1(mm)=rccemis(1,mm)
        rccems2(mm)=rccemis(2,mm)
        enddo
!      write (lun11,*)'call gsmooth on rccemis1'
      call gsmooth2(lpri,lun11,vtherm,epi,ncn2,rccems1)
!      write (lun11,*)'call gsmooth on rccemis2'
      call gsmooth2(lpri,lun11,vtherm,epi,ncn2,rccems2)
      do mm=1,ncn2
        rccemis(1,mm)=rccems1(mm)
        rccemis(2,mm)=rccems2(mm)
        enddo
!      write (lun11,*)'call gsmooth on opakc'
      call gsmooth2(lpri,lun11,vtherm,epi,ncn2,opakc)
!
!
      return
      end
