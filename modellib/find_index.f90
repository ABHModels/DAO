subroutine find_index(find,nearest_index,arr,nmaxp)
    implicit none
    integer nmaxp,i,n,nearest_index
    real(8) find,arr(nmaxp)
    real(8) min_diff,diff

    min_diff = abs(arr(1) - find)
    nearest_index = 1
    n = nmaxp

    do i = 2, n
        diff = abs(arr(i) - find)
        if (diff < min_diff) then
            min_diff = diff
            nearest_index = i
        end if
    end do
    return
end subroutine