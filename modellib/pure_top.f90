subroutine pure_top(debug,nener,ener,stinci,gamma,hcut,I_inc,I_bot)
    use constants
    implicit none

    integer,intent(in)::debug,nener
    real(8),intent(in)::ener(nener),gamma,hcut,I_bot(nener)
    real(8),intent(out)::I_inc(nener)
    character(80),intent(in)::stinci

    real(8) temb,temt,f
    integer i,nf
    nf = nener
    temb = 0.0
    temt = 0.0
    f = 0.0
    
    do i=1,nf
        if (stinci == "powlaw") &
        &   I_inc(i) =(ener(i))**(-gamma+1.)
        if (stinci == "cutoff") &
        &   I_inc(i) =(ener(i))**(-gamma+1.) * exp(-(ener(i)/hcut))  
    enddo

    call trapz(nf,I_inc,ener,temt)
    call trapz(nf,I_bot,ener,temb)
    f = temb /temt *1e-10*8.5
    ! f = temb /temt *10
    do i=1,nf
        I_inc(i) = f*I_inc(i)
    enddo

    if (debug.eq.1) then
        open(unit=20,file='./debug/Top.dat',status="unknown")
        write(20,201) (ener(i),I_inc(i)*ener(i)*1e-3,i=1,nener)
        201 format(1pe12.5,2x,1pe12.5)
        close(20)
    endif   

    return
end subroutine pure_top