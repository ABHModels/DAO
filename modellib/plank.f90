subroutine plank_law(t,nener,ener,spect)
    use constants
    implicit none

    real(8),intent(in)::t
    integer,intent(in)::nener
    real(8),intent(in)::ener(nener)
    real(8),intent(out)::spect(nener)

    integer n
    real(8) ktbb
    ktbb = t*boltz ! in ev

    do n=1,nener
        if (ener(n)/ktbb.gt.10) then
            spect(n) = 2.*exp(-ener(n)/ktbb) * ener(n)**3/spol**2
        else
            spect(n)=2.* ener(n)**3 / spol**2/(exp(ener(n)/ktbb)-1)  
        endif
    enddo

    return
end subroutine 