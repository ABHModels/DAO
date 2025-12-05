module feautrier
    implicit none
    real(8),dimension(:),allocatable::A,B,C,rhs,DD,EE

    contains    
    subroutine allocate_Feautrier(ntau,nener,nmu)
        implicit none
        integer,intent(in)::ntau,nener,nmu
        allocate(A(ntau),B(ntau),C(ntau),rhs(ntau))
        allocate(DD(ntau),EE(ntau))
    end subroutine allocate_Feautrier

        
    subroutine deallocate_Feautrier
        implicit none
        deallocate(A,B,C,rhs,DD,EE)
    end subroutine deallocate_Feautrier
    
    subroutine calc_coeff(ntau,nener,nmu,dtau,mu)
        implicit none
        integer,intent(in)::ntau,nener,nmu
        real(8),intent(in)::dtau(ntau),mu

        integer::d
        A=0.0
        B=0.0
        C=0.0
     
        do d=2,ntau-1
            A(d) = 2*mu**2/dtau(d-1)/(dtau(d)+dtau(d-1))
            C(d) = 2*mu**2/dtau(d)/(dtau(d)+dtau(d-1))
            B(d) = 1+A(d)+C(d)
        enddo

        C(1) = 2*(mu/dtau(1))**2
        B(1) = 1+2*mu/dtau(1)+2*(mu/dtau(1))**2
    
        A(ntau) = 2*(mu/dtau(ntau))**2
        B(ntau) = 1+2*mu/dtau(ntau)+2*(mu/dtau(ntau))**2
        return
    end subroutine calc_coeff

    subroutine calc_rhs(ntau,inmu_ind,mu,mu_ind,dtau,I_inc,I_bb,source)
        implicit none
        integer,intent(in)::ntau,inmu_ind,mu_ind
        real(8),intent(in)::I_inc,I_bb,source(ntau),dtau(ntau),mu
        integer::d
        rhs = 0.0
        do d=1,ntau
            rhs(d) = source(d)
        enddo
        if (inmu_ind.eq.mu_ind) rhs(1) = rhs(1) + 2*mu/dtau(1)*I_inc
        rhs(ntau) = rhs(ntau) + 2*mu/dtau(ntau)*I_bb
        return
    end subroutine calc_rhs

    subroutine Thomas(ntau,u)
        implicit none
        integer,intent(in)::ntau
        real(8),intent(out)::u(ntau)
        real(8) temp(ntau),tem
        integer::i
        DD = 0.0
        EE = 0.0
        u = 0.0
        temp = 0.0
        DD(1) = C(1)/B(1)
        EE(1) = rhs(1)/B(1)
        do i=2,ntau
            tem = B(i)-A(i)*DD(i-1)
            DD(i) = C(i)/tem
            EE(i) = (rhs(i)+A(i)*EE(i-1))/tem
        enddo

        temp(ntau) = EE(ntau)
        do i=ntau-1,1,-1
            temp(i) = EE(i)+DD(i)*temp(i+1)
        enddo
        u = temp
        return
    end subroutine Thomas

    subroutine ite_memory(step,ntau,nener,nmu,&
        &  rfzeroi,u,rfsource,&
        &  zero_old,u_old,source_old)
        implicit none
        integer,intent(in)::step,ntau,nener,nmu
        real(8),intent(inout)::rfzeroi(ntau,nener)
        real(8),intent(inout)::u(ntau,nmu,nener),rfsource(ntau,nener)
        real(8),intent(inout)::zero_old(ntau,nener)
        real(8),intent(inout)::u_old(ntau,nmu,nener),source_old(ntau,nener)
        ! copy the new variables to the old variables
        zero_old = rfzeroi
        u_old = u
        source_old = rfsource

        rfzeroi = 0.0
        u = 0.0
        rfsource = 0.0
        return
    end subroutine ite_memory

    subroutine checkconvergence(mainstep,lun11,itestep,tstart,ntau,nener,nmu,&
    &   zero_old,meanJ,source_old,source,ener,conv)
    use ieee_arithmetic 
    implicit none
    integer,intent(in)::ntau,nener,nmu,mainstep,itestep,lun11
    real(8),intent(in)::zero_old(ntau,nener),meanJ(ntau,nener),source_old(ntau,nener),source(ntau,nener)
    real(8),intent(in)::ener(nener)
    real(8),intent(in)::tstart

    logical::conv
    integer::d,n
    real(8),parameter::cirt = 1e-9
    real(8)::emis(nener),abso(nener),difsu(nener),difsd(nener),su,sd,ab,em
    real(8)::cirS(ntau),meancir,tend,ratio(ntau),meanratio

    meanratio = 0.0
    meancir = 0.0

    do d=1,ntau

        do n=1,nener
            difsu(n) = abs(source(d,n)-source_old(d,n))
            difsd(n) = source_old(d,n)
        enddo

        call trapz(nener,difsu,ener,su)
        call trapz(nener,difsd,ener,sd)
        cirS(d) = su/sd
        meancir = meancir+cirS(d)
    enddo

    meancir = meancir/ntau
    if (meancir.lt.cirt) conv = .True.

    call cpu_time(tend)
    write(lun11,*) "Main step: ",mainstep,"Iteration: ",itestep,"E[c]",meancir,"Time",tend-tstart
    print*, "Main step: ",mainstep,"Iteration: ",itestep,"E[c]",meancir,"Time",tend-tstart
    
    flush(lun11)
    end subroutine checkconvergence

end module feautrier