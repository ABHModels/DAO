subroutine inc_bot(debug,nener,ener,ktbb,I_bot,Fx)
    use constants
    implicit none
    integer,intent(in)::debug,nener
    real(8),intent(in)::ener(nener),ktbb,Fx
    real(8),intent(out)::I_bot(nener)

    ! Local variabel
    integer n
    real(8) sum,const
    real(8),dimension(:),allocatable::zremsi

    allocate(zremsi(nener))
    
    zremsi = 0.0
    I_bot = 0.0
    sum = 0.0
    do n=1,nener
        if (ener(n)/ktbb.gt.10) then
            zremsi(n) = 2.*exp(-ener(n)/ktbb) * ener(n)**3/spol**2
        else
            zremsi(n)=2.* ener(n)**3 / spol**2/(exp(ener(n)/ktbb)-1)  
        endif
        if ((n.ge.13.6).and.(n.le.1.36e+4)) &
    &       sum=sum+(zremsi(n)+zremsi(n-1))*(ener(n)-ener(n-1))/2.
    enddo

    const=Fx/sum/ergsev 
    do n=1,nener
        I_bot(n) = I_bot(n)+zremsi(n)*const
    enddo


    deallocate(zremsi)
    end subroutine inc_bot