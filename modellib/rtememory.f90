subroutine rtememory(step,ntau,nener,nmu,&
    &  u,rfzero,rfsource,&
    &  zero_old,source_old)
    implicit none
    integer,intent(in)::step,ntau,nener,nmu
    real(8),intent(inout)::rfzero(ntau,nener)
    real(8),intent(inout)::u(ntau,nmu,nener),rfsource(ntau,nener)
    real(8),intent(inout)::zero_old(ntau,nener),source_old(ntau,nener)
    ! copy the new variables to the old variables
    zero_old = rfzero
    source_old = rfsource

    rfzero = 0.0
    u = 0.0
    rfsource = 0.0
    return
end subroutine rtememory

