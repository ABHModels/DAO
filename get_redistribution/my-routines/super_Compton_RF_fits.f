!-------------------------------------------------------------------------------------------------
      subroutine super_Compton_RF_fits(itrans, theta, nmaxp, wp, df,
     & skn,  mgi, smit, agt)
c
c     This routine writes a file with the super redistribution function (SRF) for Compton
c     scatterting. For a given gas temperature T and final photon energy Ef, the SRF is
c     defined for a set of initial photon energies Ei as:
c
c         SRF(T,Ef,Ei) = IRF(Ef,Ei)/N(Ei)*skn(Ei)*dEi/Ei
c
c     where IRF(Ef,Ei) is in fact the inverse redistribution function for the Compton scattering
c     of a photon from initial energy Ei to final energy Ef; N(Ei) is the normalization to
c     ensure photon number conservation; and skn(Ei) is the Klein-Nishina cross section.
c     This routine implements the exact Compton RF from Madej et al. (2017).
c     The SRF contains all the information needed for the convolution of a given spectrum
c     to account for the Compton scattering at the given temperature.
c     Only significant values of the SRF are actually written, i.e., when RF > limit.
c
c     Input arguments:
c         itrans: Total number of temperatures
c         theta: Array (itrans) of temperatures in kT/mec2
c         nmaxp: Total number of energy points
c         wp: Array (nmaxp) of energies in eV
c         df: Array (nmaxp) of delta energies in eV
c         skn: Array (nmaxp) of Kein-Nishina cross sections
c         mgi: Total number of angles 
c         smit: Array (mgi) Legendre ordinates (angles)
c         agt: Array (mgi) weights 
c
c     Output arguments: 
c         None
c         
c     Requires:
c         probab.f: Routine for the RF calculation
c
      implicit none
      integer itrans, nmaxp, mgi, iz, np
      real*8 theta(itrans), wp(nmaxp), df(nmaxp), skn(nmaxp,itrans)
      real*8 smit(mgi), agt(mgi)
      real*8 check1(nmaxp), limit, pmin, pmax(nmaxp)
      real*8 kernel_silce(nmaxp,nmaxp)
      real*8 prob(nmaxp,nmaxp), srf
      real*8 mec2, ecen, isp, ikbol, temp(itrans),check,checka
      integer jj, kk, point(nmaxp), indices(nmaxp,nmaxp),ind_arr(nmaxp)
      integer :: n ,indmax(nmaxp)
      integer :: ni, nf,record_length

! Added by Gullo
      integer          :: unit, status
      double precision :: srf_arr(nmaxp)
! Added by Yimin for speed
      real*8 bk2_value,bk2input
      double precision bk2,skkk,smel,ccc
      parameter(skkk= 1.38062d-16,smel= 9.10956d-28,ccc= 2.99793d+10)
! Added by Yimin for store
      character (len=200) filename
      character (len=200) FILESRF,FILENER,FILESKN,FILETEM

c Initialize status
      status=0
c
      ikbol   = 1.16d4      ! inverse of kbol (K * ev-1)
      mec2  = 5.11d5        ! m_e c^2 (eV)
      isp = 0.5641895835d0  ! 1/sqrt(pi)
      limit = 1.d-3         ! Limit for the redistribution function
c
      do iz=1, itrans
            temp(iz) = theta(iz)*ikbol*mec2
c            print *,temp(iz)
      enddo

      ! FILE TO STORE kernel_silce 
      record_length = 8*nmaxp*nmaxp 
      FILESRF='exact/kernel3000.dat'
      open(unit=1, file=FILESRF, status='unknown', form='unformatted', 
     1    access='direct', recl=record_length)
      do iz = 1, itrans

         bk2input = smel*ccc**2/skkk/temp(iz)
         bk2_value= bk2(bk2input)
         
         do kk=1,nmaxp
            point(kk)=0
            pmax(kk)=0.d0
            check1(kk) = 0.d0
            indmax(kk)=0
            do np = 1, nmaxp
               prob(np,kk) = 0.d0
               indices(np,kk)=0
            enddo
         enddo

         do np = 1, nmaxp
            ecen = wp(np)
c$omp parallel num_threads(32)
c$omp& shared(iz,wp,prob,point,pmax,nmaxp,np,ecen,temp,mgi,smit,
c$omp& agt,mec2,indmax)
c$omp& private(jj)
c$omp do
      do jj=1,nmaxp
        call probab(bk2_value,temp(iz),wp(jj)/mec2,ecen/mec2,mgi,smit
     1  ,agt,prob(jj,np))
      if(prob(jj,np).gt.pmax(np)) then
            pmax(np) = prob(jj,np)
            indmax(np)=jj
      endif
      enddo
c$omp end do
c$omp end parallel
         if(pmax(np).le.0.d0) then
            print *,np
         endif
         enddo
         do np=1,nmaxp
            check=0.d0
            do jj=indmax(np),1,-1
                  kk=point(jj)+1
                  if(prob(jj,np).ge.(pmax(np)*limit))then
                        indices(jj,kk)=np
                        check=check+df(jj)*prob(jj,np)
                        point(jj)=kk
                  else
                        exit
                  endif
            enddo

            do jj=indmax(np)+1,nmaxp
                  kk=point(jj)+1
                  if(prob(jj,np).ge.(pmax(np)*limit))then
                        indices(jj,kk)=np
                        check=check+df(jj)*prob(jj,np)
                        point(jj)=kk
                  else
                        exit
                  endif

            enddo
            check1(np)=check
         enddo

         
!     From here crated by Gullo Dec 2020
!     Once it calculates the srf it writes directly the fits file
!     It needs to differentiate the first call, where it writes the extension from all the other calls 
      do jj = 1, nmaxp
         do kk=1,nmaxp
            srf_arr(kk) = 0.d0
c            ind_arr(kk) = 0
         enddo
         do kk = 1, point(jj)
            np=indices(jj,kk)
            srf = prob(jj,np)*skn(np,iz)*df(np)/wp(np)/check1(np) 
            srf_arr(kk) = srf
            kernel_silce(jj,np) = srf
c            ind_arr(kk) = np
         enddo
            n = n + 1           ! increment the row number
      enddo

!     Add by yimin, normalization
      do ni=1,nmaxp
      checka = 0.0
      do nf=2,nmaxp
      checka = checka + (kernel_silce(nf,ni)+kernel_silce(nf-1,ni))
     &      *0.5*(wp(nf)-wp(nf-1))
      enddo
      do nf=1,nmaxp
            kernel_silce(nf,ni) = kernel_silce(nf,ni)*skn(ni,iz)/checka
      enddo
      enddo
      write(1,rec=iz) kernel_silce
               
      write(9999,*)  'Wrote slice for ik = ', iz, '/', itrans
      flush(9999)
      enddo   
      close(1)
      

      FILENER='exact/ener.dat'
      FILESKN='exact/skn.dat'
      FILETEM='exact/theta.dat'
      open(unit=2,file=FILENER,status='unknown',FORM='UNFORMATTED')
      open(unit=3,file=FILESKN,status='unknown',FORM='UNFORMATTED')
      open(unit=4,file=FILETEM,status='unknown',FORM='UNFORMATTED')
      write(2) wp
      write(3) skn
      write(4) theta
      close(1)
      close(2)
      close(3)
      close(4)

      return
      end subroutine
