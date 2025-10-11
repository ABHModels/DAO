subroutine get_grids
    use Reflection_Var
    implicit none
    real(8)::dif_mu(nmurt)
    integer::i


    call enegrd
    call depgrd
    call anggrd
    
    do i=1,nmurt
        dif_mu(i) = abs(rfmu(i) - rfinmu)
    end do

    rfinmu_ind = minloc(dif_mu,1)

    if (rfdebug.eq.1) then
    call out_girds(nfrt,rfener,wrfener,'./debug/ener.dat')
    call out_girds(nmurt,rfmu,dmu,'./debug/mu.dat')
    call out_girds(ndrt,rftau,wrftau,'./debug/tau.dat')
    end if
end subroutine get_grids

subroutine out_girds(number,array,darr,filename)
    implicit none
    integer,intent(in)::number
    real(8),intent(in)::array(number)
    real(8),intent(in)::darr(number)
    character(len=*),intent(in)::filename

    integer::output_unit
    integer::i
    open(newunit=output_unit,file=filename,status='replace',action='write')
    write(output_unit,'(2E16.7)') (array(i),darr(i),i=1,number)
    flush(output_unit)
    close(output_unit)
    
end subroutine out_girds

