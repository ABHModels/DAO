    subroutine rfinit(lun11)
    use Reflection_Var
    use constants
    use globaldata
    implicit none   
    integer,intent(in)::lun11
    real(8)::rfFx,rfxi
    integer::i,j,k,e_unit
    real(8)::rftem_flux,topbb
    real(8)::cosine_in

    ! Input parameters
    call readinput(lun11)

    call allocate_rfglobal

    ! Grids
    call get_grids
    ! cosine_in = cos(pi*rfinmu/180.)
    ! cosine_in = rfmu(1)
    ! rfinmu_ind = 1

    ! Atomic data
    call atomic_data

    ! Incident radiation field
    rfxi = 10**rfzeta
    rfFx = rfnh*rfxi/fourpi            ! [erg/cm^2/s]
    if (rfsbot.eq.1)  then
        rfFxinc = rfFx*rfFx_frac
        rfFxbot = rfFx*(1-rfFx_frac)              ! [erg/cm^2/s]  
        call inc_bot(rfdebug,nfrt,rfener,rfktbb,rfi_bot,rfFxbot)
    else
        rfFxinc = rfFx
        rfi_bot = 0.0
    endif
    
    call inc_top(rfdebug,nfrt,rfener,rfstinci,rfgamma,rfhcut,rfi_inc,rfFxinc,rfinmu,rfinci_file)
    
    ! This illumination define the top illumination by flux of bottom illumination
    if (rfdebug.eq.2) then
        call t_inc_bot(nfrt,rfener,rfI_bot,rfFx)
        call t_inc_top(nfrt,rfener,rfI_bot,rfI_inc)
    endif
    
    ! First guess of radiation field
    do i=1,ndrt
     rfflux(i,:) = rfI_inc(:)+rfI_bot(:)
    enddo

    rftem_flux = 0.0
    do i=1,ndrt
     do j=2,nfrt
      rftem_flux = rftem_flux + (rfflux(i,j)+rfflux(i,j-1))*(rfener(j)-rfener(j-1))/2.
      rffluxint(i,j) = rffluxint(i,j-1) + rftem_flux*ergsev
     enddo
    enddo

    do i=1,ndrt
        rfcol(i) = wrftau(i)/sigth/1.21
        rfr(i) = wrftau(i)/rfnh/1.21/sigth
        rfdtemp(i) = rftemp_gas
        rfdxi(i) = rfzeta
    enddo

    ! First estimation for radiation field U
    do i=1,ndrt
        do j=1,nmurt
            rfu(i,j,:)=rfflux(i,:)
        enddo
        rfzero(i,:) = rfflux(i,:)/cos(pi*rfinmu/180.)
    enddo

    rfemis = 0.0
    rfabso = 0.0
    rfu = 0.0

    call loopcontrol_initial
    
    print*, " Reflected model initialization done"
end subroutine
    

    