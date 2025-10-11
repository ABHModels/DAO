

    subroutine t_inc_top(nener,ener,I_bb,I_inc)
    use constants
    implicit none
    integer,intent(in)::nener
    real(8),intent(in)::ener(nener)
    real(8),intent(in)::I_bb(nener)
    real(8),intent(out)::I_inc(nener)
    real(8) temb,temt,f,gamma
    integer i
    
    gamma  = 2.0
    do i=1,nener
        I_inc(i) =(ener(i))**(-gamma+1.)
    enddo

    temb= 0.0
    do i=2,nener
        temb = temb + 0.5*(I_bb(i)+I_bb(i-1))*(ener(i)-ener(i-1))
        temt = temt + 0.5*(I_inc(i)+I_inc(i-1))*(ener(i)-ener(i-1))
    enddo

    f = temb * 1e-10 
    do i=1,nener
        I_inc(i) = f*I_inc(i)/temt
    enddo
    call out_illumination(ener,nener,I_inc,'./illumination/top.txt')
    return
end subroutine 
