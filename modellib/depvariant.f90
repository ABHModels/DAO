    subroutine depvariant(lun11,mainstep,depstep,&
        &   ncn,epi,ncn2,r,xcol,r19,xi,zeta,xlum,t,delr,&
        &   bremsa,bremsint,zrems)
        ! Set depth dependent parameters
        use constants
        use Reflection_Var
        implicit none
        integer,intent(in)::mainstep,depstep,ncn,ncn2,lun11
        real(8),intent(out):: r,xcol,r19,xi,zeta,xlum,t,delr
        real(8),intent(out):: bremsa(ncn),bremsint(ncn),zrems(5,ncn)
        real(8),intent(in)::epi(ncn)

        integer::i
        real(8)::temp
        real(8)::lumo(nfrt)
        real(8)::skse

        temp = 0.0

        ! to xstar
        r = rfr(depstep)
        delr = 1e-15

        r19 = r*1e-19
        xcol = r*rfnh
        t = rftemp_gas/1e4     
        
        do i = 1,nfrt
            lumo(i) = rfflux(depstep,i)*fourpi*r19*r19      ! [10^38 erg/s/erg]
            bremsa(i) = rfflux(depstep,i)                  
            bremsint(i) = rffluxint(depstep,i)
        enddo

        xlum = 0.0
        do i =2,nfrt
            if (rfener(i).gt.0.1e3.and.rfener(i).lt.1000e3) then
                xlum = xlum + (lumo(i)+lumo(i-1))*(rfener(i)-rfener(i-1))/2.
            endif
        enddo

        xlum = xlum*ergsev

        xi = xlum/(rfnh*r19*r19) 
        zeta = log10(xi)
        rfdxi(depstep) = zeta
        
        return  
    end subroutine depvariant