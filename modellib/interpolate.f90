    subroutine interpolate(tararr, tardim, tarener, &
                        & cusarr, cusdim, cusener)
        implicit none
        real(8), intent(out) :: tararr(tardim)
        real(8), intent(in)  :: tarener(tardim), cusarr(cusdim), cusener(cusdim)
        integer, intent(in)  :: tardim, cusdim

        integer  :: kl, jlo
        real(8)  :: x, epmx, epmn, zr1, zr2, ep1, ep2, aly, alx, y
        logical  :: ascnd ! 用于判断 cusener 是升序还是降序
        integer  :: min_idx, max_idx ! 存储最小和最大能量的索引

        ! --- Pre-calculation before the loop ---
        ! 判断 cusener 数组是升序还是降序，并获取边界索引
        ! 这比在循环中反复调用 max/min 更高效
        ascnd = (cusener(cusdim) > cusener(1))
        if (ascnd) then
            min_idx = 1
            max_idx = cusdim
        else
            min_idx = cusdim
            max_idx = 1
        endif
        epmn = cusener(min_idx)
        epmx = cusener(max_idx)

        jlo = 0 
        do kl = 1, tardim 
            x = tarener(kl)
            
            ! --- 新的外插和插值逻辑 ---
            if (x > epmx) then
                ! 1. 外插: 超出上边界
                ! 直接使用最大能量对应的 cusarr 值
                tararr(kl) = cusarr(max_idx)
                
            else if (x < epmn) then
                ! 2. 外插: 超出下边界
                ! 直接使用最小能量对应的 cusarr 值
                tararr(kl) = cusarr(min_idx)

            else
                ! 3. 内插: 在数据范围内 (原始逻辑)
                call itfhunt3(cusener, cusdim, x, jlo) 
                jlo = max(jlo, 1) 
                
                ! 防止 jlo+1 超出数组边界 (当 x 正好等于最后一个点时可能发生)
                if (jlo >= cusdim) jlo = cusdim - 1

                zr1 = log10(max(cusarr(jlo+1), 1.d-49)) 
                zr2 = log10(max(cusarr(jlo), 1.d-49)) 
                ep1 = log10(max(cusener(jlo+1), 1.d-49)) 
                ep2 = log10(max(cusener(jlo), 1.d-49)) 
                
                alx = log10(x) 
                alx = max(alx, ep2) 
                alx = min(alx, ep1) 
                
                aly = (zr1 - zr2) * (alx - ep2) / (ep1 - ep2 + 1.d-49) + zr2 
                y = 10.0_8**aly ! 使用 10**aly 更清晰，等效于 exp(log(10)*aly)
                
                tararr(kl) = y
            endif
        enddo

    end subroutine interpolate

    subroutine itfhunt3(xx,n,x,jlo) 

        implicit none 
!                                                                       
        integer n,jlo,nint,jhi,inc,jm,nintmx 
        real(8) xx(n),x 
        logical ascnd 
!                                                                       
        data nintmx/1000/ 
!                                                                       
        nint=0                                                       
        jlo=1 
        ascnd=.false. 
        if (xx(n).gt.xx(1)) ascnd=.true. 
        if(jlo.le.0.or.jlo.gt.n)then 
        jlo=1 
        jhi=n+1 
        go to 3 
        endif 
        inc=1 
        if(x.ge.xx(jlo).eqv.ascnd)then 
    1   jhi=jlo+inc 
        if(jhi.gt.n)then 
            jhi=n+1 
        else if(x.ge.xx(jhi).eqv.ascnd)then 
            jlo=jhi 
            inc=inc+inc 
            go to 1 
        endif 
        else 
        jhi=jlo 
    2   jlo=jhi-inc                  
        if(jlo.lt.1)then 
            jlo=0 
        else if(x.lt.xx(jlo).eqv.ascnd)then 
            jhi=jlo 
            inc=inc+inc 
            go to 2 
        endif 
        endif 
        jm=n/2 
    3 continue 
        jlo=min(jlo,n) 
        jlo=max(jlo,1) 
        nint=nint+1 
        if ((jhi-jlo.eq.1).or.(nint.gt.nintmx)) return 
        jm=(jhi+jlo)/2 
        if(x.gt.xx(jm).eqv.ascnd)then 
        jlo=jm 
        else 
        jhi=jm 
        endif 
        go to 3                                                                        
    end subroutine itfhunt3


!> @brief Performs 1D log-log interpolation with boundary value extrapolation.
!> @details This subroutine interpolates a target array `tararr` at specified energy points `tarener`
!>          based on a known dataset (`cusarr` vs `cusener`). It uses linear interpolation
!>          in log-log space. For points outside the range of `cusener`, it returns the
!>          corresponding boundary value from `cusarr` (clamping).
!>
!> @param[out] tararr   The resulting interpolated array of y-values.
!> @param[in]  tardim   The dimension of the target arrays.
!> @param[in]  tarener  The array of x-values (energies) where interpolation is desired.
!> @param[in]  cusarr   The known array of y-values (e.g., cross sections).
!> @param[in]  cusdim   The dimension of the known data arrays.
!> @param[in]  cusener  The known array of x-values (energies), sorted either ascending or descending.
    subroutine interpolate_withex(tararr, tardim, tarener, &
        & cusarr, cusdim, cusener)
        implicit none

        ! --- Argument declarations ---
        real(8), intent(out) :: tararr(tardim)
        real(8), intent(in)  :: tarener(tardim), cusarr(cusdim), cusener(cusdim)
        integer, intent(in)  :: tardim, cusdim

        ! --- Local variables ---
        integer  :: kl                 ! Loop counter for the target array.
        integer  :: jlo                ! Lower index of the interpolation interval.
        real(8)  :: x                  ! The current target x-value from tarener.
        real(8)  :: epmx, epmn         ! Maximum and minimum values of the known x-array (cusener).
        real(8)  :: zr1, zr2           ! Log of y-values at the interval boundaries.
        real(8)  :: ep1, ep2           ! Log of x-values at the interval boundaries.
        real(8)  :: aly, alx           ! Log of the target y and x values.
        real(8)  :: y                  ! The calculated interpolated y-value.
        logical  :: ascnd              ! Flag to check if cusener is in ascending order.
        integer  :: min_idx, max_idx   ! Indices of the minimum and maximum energy points.
        real(8), parameter :: TINY = 1.0d-49 ! A small positive number to avoid log(0).

        ! --- Pre-calculate boundaries before the main loop for efficiency ---
        ! Determine if the source energy array is ascending or descending.
        ascnd = (cusener(cusdim) > cusener(1))

        ! Set the indices for the minimum and maximum energy points accordingly.
        if (ascnd) then
            min_idx = 1
            max_idx = cusdim
        else
            min_idx = cusdim
            max_idx = 1
        endif
        epmn = cusener(min_idx) ! Minimum energy value
        epmx = cusener(max_idx) ! Maximum energy value

        ! Initialize the lower bound index for the hunt routine.
        jlo = 0 

        ! --- Main loop over all target points ---
        do kl = 1, tardim 
            x = tarener(kl)

            ! --- Core Logic: Handle extrapolation and interpolation ---
            if (x > epmx) then
                ! Case 1: Extrapolation (above the upper bound).
                ! Return the value at the maximum energy point (clamping).
                tararr(kl) = cusarr(max_idx)

            else if (x < epmn) then
                ! Case 2: Extrapolation (below the lower bound).
                ! Return the value at the minimum energy point (clamping).
                tararr(kl) = cusarr(min_idx)

            else
                ! Case 3: Interpolation (within the data range).
                ! Find the interval containing x using the hunt routine.
                call itfhunt3(cusener, cusdim, x, jlo) 
                jlo = max(jlo, 1) 

                ! Safety check: prevent jlo+1 from going out of bounds.
                ! This can happen if x is exactly the last point in cusener.
                if (jlo >= cusdim) jlo = cusdim - 1

                ! Get the log of the y-values (cross sections) for the interval.
                zr1 = log10(max(cusarr(jlo+1), TINY)) 
                zr2 = log10(max(cusarr(jlo), TINY)) 

                ! Get the log of the x-values (energies) for the interval.
                ep1 = log10(max(cusener(jlo+1), TINY)) 
                ep2 = log10(max(cusener(jlo), TINY)) 

                ! Take the log of the target x-value.
                alx = log10(x) 
                ! Ensure alx is within the log-space interval to prevent floating point errors.
                alx = max(alx, ep2) 
                alx = min(alx, ep1) 

                ! Perform linear interpolation in log-log space.
                aly = (zr1 - zr2) * (alx - ep2) / (ep1 - ep2 + TINY) + zr2 

                ! Convert the result from log space back to linear space.
                ! 10.0_8**aly is equivalent to exp(log(10.0_8) * aly) but more readable.
                y = 10.0_8**aly 

                tararr(kl) = y
            endif
        enddo

    end subroutine interpolate_withex