subroutine gaussiankernel(itrans,theta,nmaxp,wp,df,skn)
    use ieee_arithmetic, only: ieee_is_nan
    implicit none
    integer,intent(in)::itrans,nmaxp
    real(8),intent(in)::theta(itrans),wp(nmaxp),df(nmaxp),skn(nmaxp,itrans)

    integer ik,i,j,k,ni,nf
    real(8) center,sdev,checka
    real(8)::check(nmaxp),prob(nmaxp,nmaxp),kmax(nmaxp)
    real(8),allocatable::kernel_silce(:,:)
    real(8),parameter::mec2 = 511.e3,pi=3.14169265
    real(8) sigma(nmaxp),ec(nmaxp),amp,doublepi,Gasufactor
    character (len=200) FILESRF,FILENER,FILESKN,FILETEM
    real(8),parameter::cirtf=1e-10
    integer record_length

    allocate(kernel_silce(nmaxp,nmaxp))

    record_length = nmaxp * nmaxp * 8
    print*, 'record_length = ', record_length

    doublepi = 2*pi
    FILESRF='gas/kernel.dat'
    open(unit=1, file=FILESRF, status='unknown', form='unformatted', &
         access='direct', recl=record_length)
    do ik=1,itrans
        
        call ecen_disp(theta(ik),nmaxp,ec,sigma,wp)

        kmax = 0.0
        do ni =1,nmaxp 
            amp = 1/sqrt(doublepi)/sigma(ni) 
            do nf = 1,nmaxp
                Gasufactor = exp(-((wp(nf)-ec(ni))**2)/2./(sigma(ni)**2))
                prob(ni,nf) = Gasufactor*amp
            enddo
   
            kmax(ni) = maxval(prob(ni,:))

        enddo

        do ni =1,nmaxp
            do nf=1,nmaxp
                if (prob(ni,nf).gt.kmax(ni)*cirtf) then
                    check(ni) = check(ni) + prob(ni,nf)*df(ni)
                else
                    prob(ni,nf)=0.0
                endif
            enddo
        enddo


         do ni =1,nmaxp
            do nf=1,nmaxp
                prob(ni,nf) = prob(ni,nf)*skn(ni,ik)*df(ni)/wp(ni)/check(ni)
            enddo
        enddo
      

        check =0.0
        do ni=1,nmaxp
            do nf=2,nmaxp
                check(ni) = check(ni) + (prob(ni,nf)+prob(ni,nf-1))*0.5*(wp(nf)-wp(nf-1))
            enddo
        enddo



        do ni=1,nmaxp
            do nf=1,nmaxp
                kernel_silce(nf,ni) = prob(ni,nf)*skn(ni,ik)/check(ni)
            enddo
        enddo

        write(1,rec=ik) kernel_silce
        print*, 'Wrote slice for ik = ', ik, '/', itrans


    enddo
    close(1)
    ! end loop of temperature

    FILENER='gas/ener.dat'
    FILESKN='gas/skn.dat'
    FILETEM='gas/theta.dat'
    open(unit=2,file=FILENER,status='unknown',FORM='UNFORMATTED')
    open(unit=3,file=FILESKN,status='unknown',FORM='UNFORMATTED')
    open(unit=4,file=FILETEM,status='unknown',FORM='UNFORMATTED')
    
    write(2) wp
    write(3) skn
    write(4) theta

    close(2)
    close(3)
    close(4)
    deallocate(kernel_silce)
end subroutine

function center(ener,theta,mec2) result(ec)
    real(8) ener,theta,mec2

    real(8) ec

    ec = ener*(1+4*theta-ener/mec2)

    return
end function

function sdev(ener,theta,mec2) result(sigma)
    real(8) ener,theta,mec2
    real(8) sigma

    sigma = ener*sqrt(2*theta+0.4*(ener/mec2)**2)
    return
end function
