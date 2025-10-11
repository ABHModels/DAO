subroutine t_inc_bot(nener,ener,I_bb,Fx)

    use constants
    implicit none
    integer,intent(in)::nener
    real(8),intent(in)::ener(nener),Fx
    real(8),intent(out)::I_bb(nener)
    ! ener in ev 
    ! kt in K
    ! Fx in erg/cm/s
    ! return I_bb in unit  erg / s / cm^2 / strd / ev (need to check)
    real(8) kt
    real(8) tem,scale,gamma
    integer i

    ! always set 0.35 kev
    kt = 350

    do i=1,nener
        tem = ener(i)/kt
        if (tem.gt.10) then
            I_bb(i) = 2.*exp(-ener(i)/kt) * ener(i)**3/(ccc)**2 
        else
            I_bb(i)=2.* ener(i)**3 / (ccc)**2/(exp(ener(i)/kt)-1)
        endif
    enddo

    tem = 0
    do i=2,nener
        tem = tem + 0.5*(I_bb(i)+I_bb(i-1))*(ener(i)-ener(i-1))
    enddo

    ! above integration return ev / s / cm^2 
    scale = Fx/tem/ergsev
    do i=1,nener
        I_bb(i) = scale*I_bb(i)
    enddo

    ! illumin(2,:) = I_bb(:)
    call out_illumination(ener,nener,I_bb,'./illumination/bottom.txt')
    return
    end subroutine 