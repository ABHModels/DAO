subroutine get_index(targetvalue,dim,arr,ind)
    implicit none
    integer,intent(in)::dim
    real(8),intent(in)::targetvalue,arr(dim)
    integer,intent(out)::ind

    real(8) diff(dim)
    integer::i

    do i=1,dim
        diff(i) = abs(targetvalue-arr(i))
    enddo
    ind = minloc(diff,1)
    return
end subroutine