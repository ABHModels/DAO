

subroutine out_illumination(ener,nener,spec,filename)
    implicit none
    real(8),intent(in)::ener(nener),spec(nener)
    integer,intent(in)::nener
    character(len=*),intent(in)::filename
    integer::output_unit
    integer::i
    open(newunit=output_unit,file=filename,status='replace',action='write')
    write(output_unit,3) (ener(i),spec(i),i=1,nener)
    3 format (2E16.8)
    close(output_unit)
end subroutine out_illumination
