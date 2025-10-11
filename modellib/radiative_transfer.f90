!> @brief Main driver for solving the 1D plane-parallel radiative transfer equation.
!> @author Yimin Huang, Fudan University
!>
!> @details
!> This subroutine orchestrates the full solution of the radiative transfer equation (RTE)
!> for a given hydrostatic structure. It performs the following key steps:
!>
!> 1.  **Setup & Initialization:**
!>     - Allocates memory for all necessary arrays (intensities, source functions, opacities).
!>     - Calculates energy- and depth-dependent opacities. This includes Compton scattering
!>       (`rfopsc`) and total opacity (`rfopak` = absorption + scattering).
!>     - Computes the differential optical depth `dtauE`.
!>
!> 2.  **Iterative RTE Solution:**
!>     - Implements an iterative scheme (Lambda Iteration) to find a
!>       converged solution for the radiation field.
!>     - The core of the solver is a nested loop over angle (`m`) and energy (`n`).
!>     - For each angle and energy, it solves the second-order Feautrier equation
!>       using a tridiagonal solver (`Thomas` algorithm).
!>     - After solving for all angles and energies, it calculates the zeroth moment of the
!>       radiation field (`rfzero`, i.e., the mean intensity J).
!>     - It then updates the source function (`rfsource`) using the new mean intensity.
!>     - Convergence is checked at the end of each iteration.
!>
!> 3.  **Post-processing and Output:**
!>     - After convergence, it calculates the emergent intensity and flux at the top boundary.
!>     - It writes key spectral results (e.g., emergent spectrum) to an output file.
!>     - It computes the integrated flux throughout the slab for use in subsequent steps
!>       of a larger simulation.
!>     - Finally, it deallocates all temporary arrays.
!>
!> @note This subroutine relies heavily on global variables and arrays provided by modules
!>       such as `interface` and `Reflection_Var`. Key parameters like grid sizes (`nfrt`, `ndrt`),
!>       physical conditions (`rfdtemp`, `rfnh`), and boundary conditions (`rfI_inc`) are
!>       assumed to be available from these modules.
!>
!> @param[in] lun11    Fortran logical unit number for a logging/output file.
!> @param[in] mainstep An integer identifier for the current main loop step (e.g., in a
!>                     time-dependent or hydrostatic iteration).
    subroutine radiative_transfer(lun11,mainstep)
    use constants
    use interface
    use Reflection_Var
    use feautrier
    use omp_lib
    use moment
    use redistribution
    implicit none
    real(8),parameter::zero_r8 = 0.0d0
    integer,intent(in)::lun11,mainstep
    integer::nener,ntau,nmu,inmuind,sourcetype

    logical::conv
    real(8)::pt,pb,mu
    integer n,m,d,i,j,itestep,ind_temp
    real(8)::rftstart,rftend
    real(8)::theta

    real(8) rftem_flux
    real(8) emerflux(nmurt,nfrt)
    real(8) emff(nfrt),eetot,rfFx,fracFLUX

    real(8),dimension(:),allocatable::topb,botb

    real(8),dimension(:,:),allocatable::fmoment
    real(8),dimension(:),allocatable::pointu
    real(8),dimension(:,:),allocatable::skn
    real(8),dimension(:,:),allocatable::rfopak,rfopsc
    real(8),dimension(:,:),allocatable::dtauE
    real(8),dimension(:,:),allocatable::eme
    real(8),dimension(:),allocatable::emei
    real(8),dimension(:),allocatable::meantop


    real(8) :: tsstart
    nener = nfrt
    ntau = ndrt
    nmu = nmurt
    inmuind = rfinmu_ind
    conv = .False.
    emerflux = 0.0

    if (rfmodel.eq.0) then
        sourcetype = 30
    else if (rfmodel.eq.1) then
        sourcetype = 20
    else
        print*, "Error: rfmodel not supported"
        stop
    endif

    allocate(pointu(ntau)  , skn(ntau,nener),           &
    &   rfopak(ntau,nener) , rfopsc(ntau,nener),        &
    &   topb(nener)        , botb(nener),               &
    &   eme(nmu,nener), emei(nener), & 
    &   meantop(nener)     , dtauE(ntau,nener),         &
    &   fmoment(ntau,nener),&
    &   SOURCE = zero_r8)


    ! Get compton scattering cross section
    do d=1,ntau
        theta = rfdtemp(d)*boltz/511.d3
        call get_index(theta,file_ntemp,file_temper,ind_temp)
        call get_cross_section(nener,rfener,theta,skn(d,:))

        rfopsc(d,:) = skn(d,:)*rfnh*1.21

        rfopak(d,:) = rfabso(d,:)+rfopsc(d,:)

        dtauE(d,:) = rfopak(d,:)*rfr(d)
    enddo
    
    write(lun11,*) "Main step: ",mainstep,"optical depth calculated"
    flush(lun11)

    
    ! ╔═══════════════════════════════════════════════════════════════════════════╗
    ! ║                         RADIATIVE TRANSFER SOLUTION                       ║
    ! ╚═══════════════════════════════════════════════════════════════════════════╝

    do d=1,ntau
        rfsource(d,:)=rfI_inc(:)
    enddo

    call allocate_Feautrier(ntau,nener,nmu)

    do n =1,nener

     topb(n) = 2*rfI_inc(n)/rfmu(rfinmu_ind) 
     botb(n) = rfI_bot(n)
    enddo

    call cpu_time(rftstart)
    do itestep=1,rfnmaxrte

        call cpu_time(tsstart)

        call ite_memory(itestep,ntau,nener,nmu,     &
        &  rfzero,rfu,rfsource,           &
        &  zero_old,u_old,source_old)

        do m=1,nmu
        
        mu = rfmu(m)

        do n=1,nener

            pt = topb(n)
            pb = botb(n)

            call calc_coeff(ntau,nener,nmu,dtauE(:,n),mu)

            call calc_rhs(ntau,inmuind,mu,m,dtauE(:,n),&
            &   pt,pb,source_old(:,n))

            call Thomas(ntau,pointu)

            do d=1,ntau
                rfu(d,m,n) = pointu(d)
            enddo
            
        end do

        end do

        call zero_moment(ity,ntau,nener,nmu,dmu,rfmu,rfu,rfzero)
            
        call source_function(sourcetype,ntau,nener,rfener,rfdtemp,rfzero,rfsource,skn,&
        &   rfemis,rfopsc,rfopak)

        call checkconvergence(mainstep,lun11,itestep,tsstart,ntau,nener,nmu,&
        &   zero_old,rfzero,source_old,rfsource,rfener,conv)

        if (conv) exit
    enddo

    call cpu_time(rftend)
    write(lun11,*) "Main step: ",mainstep,"RTE solution calculated, Time Used,",rftend-rftstart ,"s"
    flush(lun11)

    call deallocate_Feautrier

    meantop = 0.0
    do n=1,nener    

        do m=1,nmu

            eme(m,n) = 2*rfu(1,m,n)
            if (m.eq.rfinmu_ind) eme(m,n) = eme(m,n)-topb(n)
            emerflux(m,n) = (eme(m,n) + 2*rfu(ntau,m,n)-botb(n))*rfmu(m)
        enddo

        select case(ity)

        case(1) 

            do m=1,nmu
                meantop(n) = meantop(n) + eme(m,n)*dmu(m)
            enddo

        case(2)

        call trapz(nmu,eme(:,n),rfmu,meantop(n))

        end select
    enddo

    call trapz(nmu,emerflux,rfmu,emff)
    call trapz(nener,emff,rfener,eetot)

    eetot = eetot*fourpi/2.*ergsev
    rfFx = rfnh*10**rfzeta/fourpi
    fracFLUX = eetot/rfFx

    write(lun11,0926) fracFLUX
    0926 format ("Emergent total flux / Fx :",1E15.7)

    write(io_flux,6666) (rfener(n),topb(n),botb(n),eme(inmuind,n),meantop(n),n=1,nener)
    write(io_flux,*) " "
    6666 format(5E15.7)  ! Energy, Top boundary intensity, bottom boundary intensity, top eme flux, eme intensity
    flush(io_flux)

    rffluxint = 0.0
    do d=1,ntau
        do n = 1, nener
            rfflux(d,n) = rfzero(d,n)*fourpi
        enddo
    enddo


    rftem_flux = 0.0
    do i=1,ndrt
     do j=2,nfrt
      rftem_flux = rftem_flux + (rfflux(i,j)+rfflux(i,j-1))*(rfener(j)-rfener(j-1))/2.
      rffluxint(i,j) = rffluxint(i,j-1) + rftem_flux*ergsev
     enddo
    enddo

    deallocate(pointu,skn,rfopak,rfopsc,topb,botb,eme,emei,meantop,dtauE,fmoment)
    
    print*,"Main", mainstep,"RTE Done"
    end subroutine radiative_transfer
