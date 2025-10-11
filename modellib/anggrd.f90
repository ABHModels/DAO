subroutine anggrd
    use Reflection_Var
    implicit none

    ! Local variable
    integer m
    integer mutypp 
    character(100)::command
    integer::mu_unit
    integer::i
    real(8) delmu

    mutypp =  ity
    select case(mutypp)

    case(1)
    WRITE(command, '(A,I2,A,F4.2,A,F4.2)') 'python python-anaylsis/GLpoint.py ', &
       nmurt, ' ', rfmumin, ' ', rfmumax

    print*, command
    CALL EXECUTE_COMMAND_LINE(command)
    open(newunit=mu_unit,file='grids/mu.dat',status='old')
    do i=1,nmurt
        read(mu_unit,*) rfmu(i), dmu(i)
    end do
    close(mu_unit)

    case(2)
    ! Linear Grids (Use trapz to integral) Always from 0.05 to 0.95
    dmu(:) = 1./nmurt
    if (rfmumin.eq.0) rfmumin = 0.05
    if (rfmumax.eq.1) rfmumax = 0.95
    rfmu(1) = rfmumin
    rfmu(nmurt) = rfmumax
    do m=2,nmurt-1
        rfmu(m) = rfmu(m-1)+dmu(m)
    enddo
    end select
    return
    end subroutine anggrd