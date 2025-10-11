!> @brief Calculates the source function S for the radiative transfer equation.
!> @author Yimin Huang, Fudan University
!>
!> @details
!> This subroutine computes the source function S, a crucial term in the equation of
!> radiative transfer. It uses the input flag 'sourcetype' to select one of several
!> physical scenarios:
!>
!> 1.  **Isotropic Source:** A simple model where the source function `S` is set equal
!>     to the mean intensity `J`.
!> 2.  **Compton Source:** Accounts for energy redistribution due to Compton scattering.
!>     It computes the scattering integral of the mean intensity `J` over the
!>     Compton kernel.
!> 3.  **Total Source:** A composite model that includes both the Compton scattering
!>     contribution and a local thermal emissivity term (`emis`), all divided by the
!>     total opacity (`opac`).
!>
!> The routine performs numerical integration over the energy grid using an
!> external `trapz` function to calculate the scattering integral.
!>
!> @param[in]  sourcetype An integer flag that selects the calculation model:
!>                       - 1: Isotropic source
!>                       - 2: Pure Compton scattering source
!>                       - 3: Total source (scattering + emission)
!> @param[in]  ntau       Number of spatial/optical depth grid points.
!> @param[in]  nener      Number of energy grid points.
!> @param[in]  ener(nener) Array of energy grid points (e.g., in keV).
!> @param[in]  J(ntau,nener) The mean intensity of the radiation field, J(depth, energy).
!> @param[in]  skn(nener) Normalization factor for the scattering kernel.
!> @param[in]  emis(ntau,nener) The thermal emissivity, emis(depth, energy).
!> @param[in]  sc(ntau,nener)   The scattering opacity or coefficient, sc(depth, energy).
!> @param[in]  opac(ntau,nener) The total opacity (scattering + absorption), opac(depth, energy).
!> @param[out] S(ntau,nener)    The calculated source function, S(depth, energy).
!
subroutine source_function(sourcetype,ntau,nener,ener,t,J,S,skn,emis,sc,opac)
    use redistribution
    use constants
    integer,intent(in)::ntau,nener,sourcetype
    real(8),intent(in)::J(ntau,nener),skn(ntau,nener),ener(nener),t(ntau)
    real(8),intent(out)::S(ntau,nener)
    real(8),intent(in)::emis(ntau,nener),sc(ntau,nener),opac(ntau,nener)

    integer::d,n,nn,ind_temp
    real(8)::kernel(file_nener,file_nener),theta,source(ntau,nener)
    real(8)::integ(nener),inte_result


    SELECT CASE (sourcetype)
    
    case(10)
    ! Isotropic source
    source = J

    
    case(20)
    ! Compton source
    do d=1,ntau

        theta = t(d)*boltz/511.d3
        call get_index(theta,file_ntemp,file_temper,ind_temp)
        read(iokernel,rec=ind_temp) kernel

        do nn=1,nener ! Ei loop

            do n=1,nener ! Ef loop
                integ(n) = kernel(nn,n)*J(d,n)*(ener(nn)/ener(n))
            enddo
            call trapz(nener,integ,ener,inte_result)
        
            inte_result = inte_result/skn(d,nn)

            source(d,nn) = inte_result
    
        enddo
    enddo

    case(30)

    ! Total source
    do d=1,ntau

        theta = t(d)*boltz/511.d3
        call get_index(theta,file_ntemp,file_temper,ind_temp)
        read(iokernel,rec=ind_temp) kernel

        do nn=1,nener ! Ei loop

            do n=1,nener ! Ef loop
                integ(n) = kernel(nn,n)*J(d,n)*(ener(nn)/ener(n))
            enddo
            call trapz(nener,integ,ener,inte_result)
        
            inte_result = inte_result/skn(d,nn)

            source(d,nn) = inte_result*sc(d,nn)+emis(d,nn)
 
            source(d,nn) = source(d,nn)/opac(d,nn)
        enddo
    enddo

    CASE DEFAULT
    print*, "Error: Invalid sourcetype value: ", sourcetype
    stop "Invalid sourcetype"

    end select

    S = source

    return
end subroutine source_function