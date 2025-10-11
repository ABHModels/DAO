
subroutine trapz(n,inte,x,y)
    implicit none
    integer,intent(in)::n
    real(8),intent(in)::inte(n),x(n)
    real(8),intent(inout)::y

    integer::i
    y = 0.0
    do i=2,n
        y = y + (inte(i)+inte(i-1))*(x(i)-x(i-1))*0.5
    enddo
    return
end subroutine trapz