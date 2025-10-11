subroutine depgrd
    use Reflection_Var
    implicit none    
    ! Local Variable
    integer::d,d_unit
    real(8)::logmin,logmax,deltatau
    logmin = log10(rftaumin)
    logmax = log10(rftaumaxrt)
    deltatau = (logmax-logmin)/(ndrt-1)
    do d =1,ndrt
        rftau(d) = 10.**(logmin+(d-1)*deltatau)
    enddo

    do d=2,ndrt
        wrftau(d) = rftau(d)-rftau(d-1)
    enddo
    wrftau(1) = wrftau(2)

    end subroutine depgrd