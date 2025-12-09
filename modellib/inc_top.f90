    subroutine inc_top(debug,nener,ener,stinci,gamma,ehcut,I_inc,Fx,inangle,&
        &   rfinci_file)
    use constants
    implicit none

    integer :: ios
    integer :: unit,ener_unit
    
    integer::nrows
    character(len=100) :: line
    integer,intent(in)::debug,nener
    real(8),intent(in)::gamma,ehcut,Fx,inangle
    character(80),intent(in)::stinci
    character(80),intent(in)::rfinci_file


    real(8),intent(out)::I_inc(nener)   
    real(8),intent(inout)::ener(nener)
    

    ! Loacl
    integer,parameter::nth_np = 5,comptt_np=5
    real nth_par(nth_np),comptt_par(comptt_np)
    integer::ifl
    real Ear(0:nener),Earout(nener)
    real Photar(nener),Photer(nener)


    integer i,numcon,fileunit
    real(8) sum,const,tp,nth_kte,nth_gam
    real(8) ecut                                ! Low energy cut
    real(8), dimension(:), allocatable :: zremsi
    real(8), dimension(:), allocatable :: col2,col1
    integer::custumener
    character(100)::command
    
    

    ! return I_inc in unit  erg / s / cm^2 / strd / erg 
    ! **Warning:
    !     The unit must to be right to enter XSTAR, the best ways is to use Xstar choosen.

    allocate(zremsi(nener))

    custumener = 1
    numcon=nener
    sum=0.
    ecut = 0.1
    zremsi(:)=0.0
    I_inc(:)=0.0
    const=1.

    select case(stinci)

    case ("file")

    nrows = 99999
    open(newunit=fileunit, file=rfinci_file, status='old', action='read')
    read(fileunit,*) nrows
    read(fileunit,*) unit      ! 0:erg/cm^2/s/erg 1:Photons/cm^2/s/erg 2:erg^2/cm^2/s/erg
    read(fileunit,*) ener_unit ! 0:eV 1:keV 
    allocate(col1(nrows),col2(nrows))
    col1 = 0.0
    col2 = 0.0
    do i = 1, nrows
        read(fileunit,*) col1(i), col2(i)
        if (ener_unit.eq.1) col1(i) = col1(i)*1e3
    enddo
    close(fileunit)
    call interpolate_withex(zremsi,nener,ener,col2,nrows,col1)
    deallocate(col1,col2)

    ! To standard form erg/cm^2/s/erg 
    do i =1,nener
        if (unit.eq.1) zremsi(i) = zremsi(i)*ergsev*ener(i)
        if (unit.eq.2) zremsi(i) = zremsi(i)/ener(i)/ergsev
    enddo

    case ("powlaw")
    tp = 1.0-gamma

    do i=1,numcon 
        zremsi(i)=1.e-24
        zremsi(i)=ener(i)**tp
    enddo

    case ("cutoff")
    tp = 1.0-gamma

    do i=1,numcon
        zremsi(i)=1.e-24
        zremsi(i)=ener(i)**tp * exp(-(ener(i)/ehcut))
    enddo

    case ("blackbody")
    tp = gamma
    do i =1,numcon
        zremsi(i)=1.e-24
        if (ener(i)/tp.gt.10) then
            zremsi(i) = 2.*exp(-ener(i)/tp)*ener(i)**3/spol**2 
        else
            zremsi(i) = 2.*ener(i)**3/spol**2/(exp(ener(i)/tp)-1)
        endif
    enddo

    case ("nthcomp")
    nth_par(1) = gamma 
    nth_par(2) = ehcut/1e3 
    nth_par(3) = ktb_nthcomp/1e3
    nth_par(4) = 0.0
    nth_par(5) = 0.0

    do i=0,nener-1
        Ear(i) = ener(i+1)/1e3
    enddo
    ! TODO: It's a very rought numbe, could we do better?
    Ear(nener) = 1.001d3

    ! Nthcomp model, return in [emodel]
    ifl=1
    call donthcomp(Ear,nener,nth_par,ifl,Photar,Photer)  
    do i = 1,nener
        Earout(i) = (Ear(i-1)+Ear(i))*0.5
        zremsi(i) = Photar(i)
    enddo
    
    case("comptt")
    comptt_par(1) = 0.0         ! Red shift
    comptt_par(2) = gamma/1e3   ! Wien temperature (keV) 
    comptt_par(3) = hcut/1e3    ! plasma temperature (keV)
    comptt_par(4) = 1.0         ! Optical depth
    comptt_par(5) = 1.0         ! Depth geometry

    do i=0,nener-1
        Ear(i) = ener(i+1)/1e3
    enddo
    Ear(nener) = 1.001d3

    ifl=1
    call XSTITG(Ear,nener,comptt_par,ifl,Photar,Photer)  
    do i = 1,nener
        Earout(i) = (Ear(i-1)+Ear(i))*0.5
        zremsi(i) = Photar(i)
    enddo

    print*,"It hasn't been tested."
    end select

    ! _Normalization_

    do i=2,nener
        if (ener(i).gt.0.1e3.and.ener(i).lt.1000e3) then
            sum=sum+(zremsi(i)+zremsi(i-1))*(ener(i)-ener(i-1))/2.
        endif
    enddo
    const=Fx/sum/ergsev         

    do i=1,numcon
        I_inc(i) = I_inc(i)+zremsi(i)*const   
    enddo

    deallocate(zremsi)
    return
    end subroutine inc_top

   


    

