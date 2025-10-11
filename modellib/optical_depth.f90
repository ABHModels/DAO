subroutine optical_depth(debug,lun11,mainstep,ntau,nener,dtau,opacity,delr)
    implicit none
    integer,intent(in)::debug,lun11,mainstep,ntau,nener
    real(8),intent(in)::opacity(ntau,nener),delr(ntau)
    real(8),intent(out)::dtau(ntau,nener)
    integer::d,n

    do d=1,ntau
        do n=1,nener
            dtau(d,n) = opacity(d,n)*delr(d)
        enddo
    enddo

    write(lun11,*) "Main step: ",mainstep,"optical depth calculated"
    flush(lun11)    
end subroutine optical_depth

subroutine opticaldepth(dtau_e,skn)
    use Reflection_Var
    use constants
    implicit none 
    integer::n,d
    real(8)::tau(ndrt)
    real(8),intent(inout)::dtau_e(ndrt,nfrt)
    real(8),intent(in)::skn(nfrt)
    tau = rftau

    dtau_e= 0.0
    do n=1,nfrt
        dtau_e(1,n)=(tau(2)-tau(1))*skn(n)/sigth
        dtau_e(ndrt,n)=(tau(ndrt)-tau(ndrt-1))*skn(n)/sigth
    enddo
    do d=2,ndrt-1
        do n=1,nfrt
            dtau_e(d,n)=(tau(d+1)-tau(d))*skn(n)/sigth
        enddo
    enddo
    end subroutine