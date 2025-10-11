subroutine pure_bot(debug,nener,ener,ktbb,I_bot,Fx)
    use constants
    implicit none

    integer,intent(in)::debug,nener
    real(8),intent(in)::ener(nener),ktbb,Fx
    real(8),intent(out)::I_bot(nener)

    real(8) tem,scale,kt,ergener(nener)
    integer i,nf
    nf = nener
    kt = ktbb*ergsev
    do i=1,nf
        ergener(i) = ener(i)*ergsev
    enddo

    do i=1,nf
        if(ergener(i)/kt<0.01) then
            I_bot(i) = 2. * kt * (ergener(i)/(herg*ccc))**2
        else if (ergener(i)/kt>10) then
            I_bot(i) = 2. * ergener(i)**3 / (herg*ccc)**2 * exp(-ergener(i)/kt) 
        else
            I_bot(i) = 2. * ergener(i)**3 / (herg*ccc)**2/(exp(ergener(i)/kt)-1.)
        endif
    enddo

    ! tem = 0
    ! call trapz(nf,I_bot,ener,tem)
    ! scale = Fx/tem/ergsev
    ! do i=1,nf
    !     I_bot(i) = scale*I_bot(i)
    ! enddo
    do i=1,nf
        I_bot(i) = I_bot(i)/ergsev/hev
    enddo

    if (debug.eq.1) then
        open(unit=30,file='./debug/Bot.dat',status="unknown")
        ! Out intenristy in unit of kev^2/cm^2/s/kev
        write(30,31) (ener(i),I_bot(i)*ener(i)*1e-3,i=1,nener)
        31 format(1pe12.5,2x,1pe12.5)
        close(30)
    endif

    return
end subroutine pure_bot