    subroutine emis_abso(lun11,mainstep,depstep,ncn,&
        &   rccemis,flinel,brcems,opakc,t)
        ! This rouitne get emissivity and absorption from xstar
        use Reflection_Var
        use constants
        implicit none
        integer,intent(in)::mainstep,depstep,lun11,ncn
        real(8),intent(inout)::rccemis(2,ncn),flinel(ncn),brcems(ncn),opakc(ncn)
        real(8),intent(in)::t

        real(8) emiss(nfrt)
        real(8)::temdo,temna
        integer::i

        emiss = 0.0
        ! Norm emissivity
        temdo = 0.0
        temna = 0.0
        
        do i=1,nfrt
            emiss(i) = (flinel(i)+rccemis(2,i)+rccemis(1,i)+brcems(i))*2*pi
        enddo

        do i=2,nfrt
            temna = temna + (emiss(i)+emiss(i-1))/2.*(rfener(i)-rfener(i-1))
            temdo = temdo+(rfzero(depstep,i)*opakc(i)+rfzero(depstep,i-1)*opakc(i-1))/2.*(rfener(i)-rfener(i-1))
        enddo
        
        do i=1,nfrt
            emiss(i) = emiss(i)*temdo/temna
        enddo

        do i = 1,nfrt
            rfemis(depstep,i) = emiss(i)
            rfabso(depstep,i) = opakc(i)
        enddo
        
        rfdtemp(depstep) =t*1e4
        
        return

    end subroutine emis_abso