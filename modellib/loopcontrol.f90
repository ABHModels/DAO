subroutine loopcontrol_initial 
    use Reflection_Var
    implicit none
    allocate(told(ndrt),xiold(ndrt))
    mainloopsw =.False.
    told = rfdtemp
    xiold = rfzeta
end subroutine loopcontrol_initial

subroutine loopcontrol_save
    use Reflection_Var
    implicit none
    told = log10(rfdtemp)
    xiold = rfdxi
end subroutine loopcontrol_save


subroutine loopcontrol_conv(step,lun11)
    use Reflection_Var
    implicit none
    integer,intent(in)::lun11,step
    real(8),parameter::cirT=1e-3,cirXi=1e-3
    real(8)::meanTold,meanXiold,meanT,meanXi
    integer::d

    meanTold = 0.0
    meanXiold = 0.0
    meanT= 0.0
    meanXi = 0.0

    do d=1,ndrt
    meanTold = meanTold+told(d)
    meanT = meanT + log10(rfdtemp(d))

    meanXi = meanXi+rfdxi(d)
    meanXiold = meanXiold+xiold(d)
    enddo

    meanTold = meanTold/ndrt
    meanXiold= meanXiold/ndrt
    meanT = meanT/ndrt
    meanXi = meanXi/ndrt

    if (abs(meanXi-meanXiold).le.cirXi &
    & .and. (abs(meanT-meanTold).le.cirT)) mainloopsw = .True.
    write(lun11, '(A, I0, 2X, A, ES14.6, 2X, A, ES14.6)') &
    "Main step", step, "E[xi]", abs(meanXi-meanXiold), "E[T]", abs(meanT-meanTold)

end subroutine loopcontrol_conv