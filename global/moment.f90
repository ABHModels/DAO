module moment
    implicit none
    contains
    subroutine zero_moment(ity,ntau,nener,nmu,dmu,rfmu,rfu,rfzero)
        implicit none
        integer,intent(in)::ntau,nener,nmu,ity
        real(8),intent(in)::rfmu(nmu),dmu(nmu)
        real(8),intent(in)::rfu(ntau,nmu,nener)
        real(8),intent(out)::rfzero(ntau,nener)

    integer::d,m,n
    real(8)::inte

    rfzero = 0.0

    select case (ity) 

    case(1)
    ! Change for GL integration, drop trapz
    do d=1,ntau
        do n=1,nener
            inte = 0.0
            do m=1,nmu
                inte = inte + rfu(d,m,n)*dmu(m)
            enddo
            rfzero(d,n) = inte
        enddo
    enddo

    case(2)
    do d=1,ntau
        do n=1,nener
            call trapz(nmu,rfu(d,:,n),rfmu,rfzero(d,n))
        enddo
    enddo

    end select

    return
    end subroutine

    subroutine first_moment(ntau,nener,nmu,rfmu,rfu,rffiri)
    implicit none
    integer,intent(in)::ntau,nener,nmu
    real(8),intent(in)::rfmu(nmu)
    real(8),intent(in)::rfu(ntau,nmu,nener)
    real(8),intent(out)::rffiri(ntau,nener)

    integer::d,m,n
    real(8)::inte(nmu)

    rffiri = 0.0
    do d=1,ntau
        do n=1,nener
            inte = 0.0
            do m=1,nmu
                inte(m) = rfu(d,m,n)*rfmu(m)
            enddo
            call trapz(nmu,inte,rfmu,rffiri(d,n))
        enddo
    enddo
    return
    end subroutine

    subroutine second_moment(ntau,nener,nmu,rfmu,rfu,rfseci)
    implicit none
    integer,intent(in)::ntau,nener,nmu
    real(8),intent(in)::rfmu(nmu)
    real(8),intent(in)::rfu(ntau,nmu,nener)
    real(8),intent(out)::rfseci(ntau,nener)

    integer::d,m,n
    real(8)::inte(nmu)

    rfseci = 0.0
    do d=1,ntau
        do n=1,nener
            inte = 0.0
            do m=1,nmu
                inte(m) = rfu(d,m,n)*rfmu(m)**2
            enddo
            call trapz(nmu,inte,rfmu,rfseci(d,n))
        enddo
    enddo
    return
    end subroutine

    subroutine emergentinte(nener,nmu,rfmu,rfu,emergent)
        implicit none
        integer,intent(in)::nener,nmu
        real(8),intent(in)::rfmu(nmu)
        real(8),intent(in)::rfu(nmu,nener)
        real(8),intent(out)::emergent(nener)
    
        integer::m,n
        real(8)::inte(nmu)
    
        emergent = 0.0
        do n=1,nener
            inte = 0.0
            do m=1,nmu
                inte(m) = rfu(m,n)
            enddo
            call trapz(nmu,inte,rfmu,emergent(n))
        enddo
        return
        end subroutine
end module