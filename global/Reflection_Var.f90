module Reflection_Var
    implicit none

    !--------------------BEGIN OF INPUT PARAMETERS -----------------------------

    !-------------------------------------------------------
    !                           Debug
    !-------------------------------------------------------
    integer::rfdebug

    !-------------------------------------------------------
    !                           Grids
    !-------------------------------------------------------

    integer::ndrt,nfrt,nmurt 
    integer::rfnmaxmain,rfnmaxrte                               ! Max iteration steps for RTE, number of iterations for temperature equlibrium
    real(8)::rfppemin,rfppemax,rfppemax2
    real(8)::rftaumin,rftaumaxrt
    real(8)::rfmumin,rfmumax

    integer::rfmodel

    integer::ity            ! 1: GSpoint, 2: linear

    !-------------------------------------------------------
    !                   illumination
    !-------------------------------------------------------
    integer::rfsbot                                             ! Bottom illumination swith
    real(8)::rfgamma,rfhcut,rfinmu,rfFx_frac                    ! Corona illumination Photon Index,High Cut, incident angle (deg)
    character(80)::rfinci_file                                  ! Corona illumination file
    character(80)::rfstinci                                     ! Corona illumination type (blackbody, powlaw, cutoff)


    !-------------------------------------------------------
    !                   Atmosphere
    !-------------------------------------------------------
    real(8)::rfzeta,rfnh                                        ! log10(ionization para), hydrogen density
    real(8)::rfktbb                                             ! Bottom illumination temperature 
    real(8)::rftemp_gas                                         ! Gas temperature
    character(80)::rftemp_gas_unit                              ! Unit of gas temperature
    integer::rfit_temp                                          ! Temperature equilibrium iteration steps
    
    !-------------------------------------------------------
    !                   Scattering
    !-------------------------------------------------------
    character(80)::rfcomp_file                                        ! Compton redistribution file name
    ! Under the path of rfcomp_file, there must be four files: 
    !   1. Temperature grid
    !   2. Energy grid
    !   3. Scattering kernel
    !   4. Scattering cross section

    !-------------------------------------------------------
    !                   Atomic data
    !-------------------------------------------------------
    real(8):: rfh
    real(8):: rfhe
    real(8):: rfli
    real(8):: rfbe
    real(8):: rfba
    real(8):: rfc
    real(8):: rfn
    real(8):: rfo
    real(8):: rff
    real(8):: rfne
    real(8):: rfna
    real(8):: rfmg
    real(8):: rfal
    real(8):: rfsi
    real(8):: rfp
    real(8):: rfs
    real(8):: rfcl
    real(8):: rfar
    real(8):: rfk
    real(8):: rfca
    real(8):: rfsca
    real(8):: rfti
    real(8):: rfva
    real(8):: rfcr
    real(8):: rfmn
    real(8):: rffe
    real(8):: rfco
    real(8):: rfni
    real(8):: rfcu
    real(8):: rfzn
    real(8):: rfabel(30) ! abundance in array  

    !-------------------------------------------------------
    !                   Output
    !-------------------------------------------------------
   
    character(80)::rfop_spec                                  ! output file name
    character(80)::rfop_temp                                  ! output file name
    character(80)::rfop_inte                                  ! output file name
    character(80)::rfop_log                                   ! output file name
    character(80)::rfop_emis                                  ! output file name
    character(80)::rfdataenv                                  ! Data environment
    character(80)::rfabdfits                                  ! Ion, Heating, Cooling ...
    

    !--------------------END OF INPUT PARAMETERS -----------------------------
    

    !--------------------BEGIN OF GLOBAL VARIABLES -----------------------------


    ! Grids
    real(8),dimension(:),allocatable::rfener,wrfener
    integer,dimension(:),allocatable::rfind_ener
    real(8),dimension(:),allocatable::rfmu,dmu
    real(8),dimension(:),allocatable::rftau,wrftau     
    integer::rfinmu_ind

    ! IO
    integer::io_temper
    integer::io_intensity
    integer::io_flux
    integer::io_emis
    
    ! radiation field
    real(8)::rfFxinc,rfFxbot
    real(8),dimension(:),allocatable::rfI_inc,rfI_bot                       ! Incident radiation field [erg/cm^2/s/erg]
    real(8),dimension(:,:),allocatable::rfemis,rfabso                       ! Emission and absorption 
    real(8),dimension(:,:,:),allocatable::rfu                               ! Radiation field [erg/cm^2/s/ster/erg]
    real(8),dimension(:,:),allocatable::rfflux                              ! Flux [erg/cm^2/s/erg]
    real(8),dimension(:,:),allocatable::rffluxint                           ! Flux integrated [erg/cm^2/s]
    real(8),dimension(:,:),allocatable::rfzero                             ! Zero moment of intensity [erg/cm^2/s/erg] 
    real(8),dimension(:,:),allocatable::rfsource                            ! Source function 

    ! Memory for radiative filed
    real(8),dimension(:,:),allocatable::zero_old,source_old
    real(8),dimension(:,:,:),allocatable::u_old

    ! Depth variant parameters
    real(8),dimension(:),allocatable::rfcol                                 ! Column density [cm^2/g]
    real(8),dimension(:),allocatable::rfr                                   ! Radius
    real(8),dimension(:),allocatable::rfdzeta                               ! Redshift
    real(8),dimension(:),allocatable::rfdtemp                               ! Temperature gradient

    contains
    subroutine allocate_rfglobal
        implicit none
        real(8),parameter :: zero_r8 = 0.0d0
        integer,parameter :: zero_int = 0

        allocate(rfener(nfrt),       wrfener(nfrt),      &
                 rfmu(nmurt),        dmu(nmurt),         &
                 rftau(ndrt),        wrftau(ndrt),       &
                 rfI_inc(nfrt),      rfI_bot(nfrt),      &
                 rfemis(ndrt,nfrt),  rfabso(ndrt,nfrt),  & 
                 rfu(ndrt,nmurt,nfrt),                    &
                 rfflux(ndrt,nfrt),  rffluxint(ndrt,nfrt),&
                 rfzero(ndrt,nfrt), rfsource(ndrt,nfrt), &
                 zero_old(ndrt,nfrt),source_old(ndrt,nfrt),&
                 rfcol(ndrt),        rfr(ndrt),          &
                 rfdzeta(ndrt),      rfdtemp(ndrt),      &
                 u_old(ndrt,nmurt,nfrt),&
                SOURCE = zero_r8 ) 

        allocate(rfind_ener(nfrt),SOURCE=zero_int)
        
    end subroutine allocate_rfglobal

    subroutine deallocate_rfglobal
        deallocate(rfener,wrfener,rfind_ener)
        deallocate(rfmu,dmu)
        deallocate(rftau,wrftau)

        deallocate(rfI_inc,rfI_bot)
        deallocate(rfemis,rfabso)
        deallocate(rfu)
        deallocate(rfflux,rffluxint)
        deallocate(rfzero)
        deallocate(zero_old,source_old)
        deallocate(u_old)
        deallocate(rfcol,rfr,rfdzeta,rfdtemp)
    end subroutine deallocate_rfglobal
end module Reflection_Var