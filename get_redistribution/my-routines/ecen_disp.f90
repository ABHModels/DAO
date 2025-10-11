subroutine ecen_disp(theta,ncn2,mean,std_dev,epi)
    implicit none
    real(8) mean(ncn2),std_dev(ncn2)
    integer ncn2,n,mm
    real(8) epi(ncn2),me2,theta,epsilon
    mean(:) =0.
    std_dev(:) =0.
    me2 = 511.d3
    do n = 1,ncn2
        epsilon = epi(n)/me2
        mean(n) = epi(n)*(1.+4.*theta-epsilon) 
        std_dev(n) = epi(n) * sqrt(2*theta+0.4*epsilon**2)
    enddo 
    return
    end
