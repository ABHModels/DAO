module interface
    use Reflection_Var
    use globaldata
    use ieee_arithmetic
    implicit none
    ! integer::itf_nener,itf_ntau,itf_nmu
    ! real(8)::itf_t         ! temperaure [K] for xstar initialization
    ! real(8)::itf_xpx            ! density
    ! real(8)::itf_abel(nl)
    ! real(8),dimension(:),allocatable::itf_flux,itf_fluxint,itf_ener
    ! real(8)::itf_xpcol
    ! real(8)::itf_zeta
    ! real(8)::itf_nlimd
    ! real(8)::itf_r,itf_r19

    contains
    ! subroutine itfinit(lun11)
    !     implicit none
    !     integer,intent(in)::lun11
    !     itf_nener = nfrt
    !     itf_ntau = ndrt
    !     itf_nmu = nmurt
    !     allocate(itf_flux(itf_nener),itf_fluxint(itf_nener),itf_ener(itf_nener))

    !     itf_t = 0.0d0
    !     itf_xpx = 0.0d0
    !     itf_ener = 0.0d0
    !     itf_flux = 0.0d0
    !     itf_fluxint = 0.0d0
    !     itf_xpcol = 0.0d0
    !     itf_zeta = 0.0d0
    !     itf_nlimd = 0.0d0
    !     itf_abel(:) = 0.0d0
    !     if (rfdebug.eq.2) then
    !         write(lun11,*) "Interface initialized"
    !     endif
    ! end subroutine itfinit

    ! subroutine itfconst(lun11)
    !     ! Set depth invariant parameters
    !     implicit none
    !     integer,intent(in)::lun11
    !     itf_t = rftemp_gas
    !     itf_xpx = rfnh
    !     itf_nlimd = rfit_temp
    !     itf_ener(:) = rfener(:)
    !     itf_abel(:) = rfabel(:)
    !     print*, " Interface initialization done"
    ! end subroutine


    subroutine loopconstant(trad,xlum,lwri,lpri,r,t,xpx,p,lcdd,          &
        & numrec,npass,nlimd,rmax,xpxcol,xi,zeta,lfix,                     &
        & lun11,abel,cfrac,emult,taumax,xeemin,spectype,specfile,specunit, &
        & kmodelname,nloopctl,critf,vturbi,eptmp,zrtmp,numcon2,ncn2,radexp)
        ! read depth invariant parameters
         implicit none                                                           
         real(8) eptmp(ncn),zrtmp(ncn),abel(nl),abel2(30) 
         character(8) stringst,kblnk8 
         character(80) specfile,spectype,stringsl,kblnk80,stringst2 
         character(30) kmodelname 
         integer nloopctl,specunit,ierr,ll,lcdd2,lun13,nenergy,ncn2 
         integer lwri,lpri,lcdd,numrec,npass,nlimd,lfix,lun11,numcon2,mm 
         real(8) trad,xlum,r,t,xpx,p,rmax,xpxcol,xi,zeta,cfrac,emult,taumax,&
        &     xeemin,critf,vturbi,ccc,xlum2,xpcol,r19,radexp                                                       
     
!       density                                                     
        xpx = rfnh

        numcon2 = nfrt
                                                                                                     
!       number of steps                                             
        numrec = ndrt

!       number of iterations                                        
        nlimd = rfit_temp
                                                      
!       abundances                                                  
        do ll=1,nl
        abel2(ll)=rfabel(ll)
        enddo
        do mm=1,nl 
            abel(mm)=abel2(mm)
        enddo 

        ncn2 = numcon2
        print*, " Xstar constant parameters done"
        return 
    end subroutine loopconstant


    subroutine depvariant(lun11,mainstep,depstep,&
        &   ncn,epi,ncn2,r,xcol,r19,xi,zeta,xlum,t,delr,&
        &   bremsa,bremsint,zrems)
        ! Set depth dependent parameters
        use constants
        implicit none
        integer,intent(in)::mainstep,depstep,ncn,ncn2,lun11
        real(8),intent(out):: r,xcol,r19,xi,zeta,xlum,t,delr
        real(8),intent(out):: bremsa(ncn),bremsint(ncn),zrems(5,ncn)
        real(8),intent(in)::epi(ncn)

        integer::i
        real(8)::temp
        real(8)::lumo(nfrt)
        real(8)::skse

        ! ref to interface
        ! do i =1,nfrt
        !     itf_flux(i) = rfflux(depstep,i)       
        !     itf_fluxint(i) = rffluxint(depstep,i)
        ! enddo
        ! itf_xpcol = rfcol(depstep)
        ! itf_r = rfr(depstep)
        ! itf_r19 = rfr19(depstep)
        ! itf_t = rftemp_gas

        temp = 0.0
        ! call trapz(nfrt,rfflux(depstep,:),rfener(:)*ergsev,temp)


        ! to xstar
        r = rfr(depstep)
        delr = 1e-15

        r19 = r*1e-19
        xcol = r*rfnh
        t = rftemp_gas/1e4     
        
        do i = 1,nfrt
            lumo(i) = rfflux(depstep,i)*fourpi*r19*r19      ! [10^38 erg/s/erg]
            bremsa(i) = rfflux(depstep,i)                  
            bremsint(i) = rffluxint(depstep,i)
        enddo

        xlum = 0.0
        do i =2,nfrt
            if (rfener(i).gt.0.1e3.and.rfener(i).lt.1000e3) then
                xlum = xlum + (lumo(i)+lumo(i-1))*(rfener(i)-rfener(i-1))/2.
            endif
        enddo

        xlum = xlum*ergsev

        xi = xlum/(rfnh*r19*r19) 
        zeta = log10(xi)

        return  
    end subroutine depvariant

    subroutine emis_abso(lun11,mainstep,depstep,ncn,&
        &   rccemis,flinel,brcems,opakc,t)
        ! This rouitne get emissivity and absorption from xstar
        use Reflection_Var
        use constants
        implicit none
        integer,intent(in)::mainstep,depstep,lun11,ncn
        real(8),intent(inout)::rccemis(2,ncn),flinel(ncn),brcems(ncn),opakc(ncn)
        real(8),intent(in)::t

        real(8) emiss(nfrt)
        real(8)::temdo,temna
        integer::i

        emiss = 0.0
        ! Norm emissivity
        temdo = 0.0
        temna = 0.0
        
        do i=1,nfrt
            emiss(i) = (flinel(i)+rccemis(2,i)+rccemis(1,i)+brcems(i))*2*pi
        enddo

        do i=2,nfrt
            temna = temna + (emiss(i)+emiss(i-1))/2.*(rfener(i)-rfener(i-1))
            temdo = temdo+(rfzero(depstep,i)*opakc(i)+rfzero(depstep,i-1)*opakc(i-1))/2.*(rfener(i)-rfener(i-1))
        enddo
        
        do i=1,nfrt
            emiss(i) = emiss(i)*temdo/temna
        enddo

        do i = 1,nfrt
            rfemis(depstep,i) = emiss(i)
            rfabso(depstep,i) = opakc(i)
        enddo
        
        rfdtemp(depstep) =t*1e4
        
        return

    end subroutine emis_abso
end module


