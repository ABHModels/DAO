! This routine is use fortran to implenment the fgsl 
! Author Yimin Huang
! 17th Dec

!-----------------------------------------------------------------------
! XSPEC/FGSL - GAMMI function
! P(a, x) = (1/Gamma(a)) * integral(0 to x) t^(a-1) e^(-t) dt
!-----------------------------------------------------------------------
      FUNCTION GAMMI(A,X)
      IMPLICIT NONE
      REAL A, X, GAMMI
      
      ! 局部变量
      DOUBLE PRECISION gamser, gammcf, gln
      DOUBLE PRECISION a_d, x_d
      
      a_d = DBLE(A)
      x_d = DBLE(X)

      IF (x_d .LT. 0.0d0 .OR. a_d .LE. 0.0d0) THEN
         PRINT *, 'Error: GAMMI called with x < 0 or a <= 0'
         GAMMI = 0.0
         RETURN
      ENDIF

      ! 根据 x 和 a 的大小关系选择算法
      IF (x_d .LT. (a_d + 1.0d0)) THEN
         ! 使用级数展开 (Series Representation)
         CALL gser(gamser, a_d, x_d, gln)
         GAMMI = SNGL(gamser)
      ELSE
         ! 使用连分数 (Continued Fraction)
         ! 注意: gcf 计算的是上不完全伽马 Q(a,x)，所以 P(a,x) = 1 - Q(a,x)
         CALL gcf(gammcf, a_d, x_d, gln)
         GAMMI = SNGL(1.0d0 - gammcf)
      ENDIF
      
      RETURN
      END FUNCTION GAMMI

!-----------------------------------------------------------------------
! (Series Representation)
!-----------------------------------------------------------------------
      SUBROUTINE gser(gamser, a, x, gln)
      IMPLICIT NONE
      DOUBLE PRECISION gamser, a, x, gln
      INTEGER ITMAX
      DOUBLE PRECISION EPS
      PARAMETER (ITMAX=100, EPS=3.0d-7)
      INTEGER n
      DOUBLE PRECISION sum, del, ap
      
      DOUBLE PRECISION gammaln ! External function
      
      gln = gammaln(a)
      IF (x .LE. 0.0d0) THEN
         IF (x .LT. 0.0d0) stop "x < 0 in gser"
         gamser = 0.0d0
         RETURN
      ENDIF
      ap = a
      sum = 1.0d0 / a
      del = sum
      DO n = 1, ITMAX
         ap = ap + 1.0d0
         del = del * x / ap
         sum = sum + del
         IF (ABS(del) .LT. ABS(sum)*EPS) GOTO 1
      END DO
      stop "a too large, ITMAX too small in gser"
 1    gamser = sum * EXP(-x + a*LOG(x) - gln)
      RETURN
      END SUBROUTINE gser

!-----------------------------------------------------------------------
! (Continued Fraction)
!-----------------------------------------------------------------------
      SUBROUTINE gcf(gammcf, a, x, gln)
      IMPLICIT NONE
      DOUBLE PRECISION gammcf, a, x, gln
      INTEGER ITMAX
      DOUBLE PRECISION EPS, FPMIN
      PARAMETER (ITMAX=100, EPS=3.0d-7, FPMIN=1.0d-30)
      INTEGER i
      DOUBLE PRECISION an, b, c, d, del, h
      
      DOUBLE PRECISION gammaln ! External function

      gln = gammaln(a)
      b = x + 1.0d0 - a
      c = 1.0d0 / FPMIN
      d = 1.0d0 / b
      h = d
      DO i = 1, ITMAX
         an = -DBLE(i) * (DBLE(i) - a)
         b = b + 2.0d0
         d = an * d + b
         IF (ABS(d) .LT. FPMIN) d = FPMIN
         c = b + an / c
         IF (ABS(c) .LT. FPMIN) c = FPMIN
         d = 1.0d0 / d
         del = d * c
         h = h * del
         IF (ABS(del - 1.0d0) .LT. EPS) GOTO 1
      END DO
      stop "a too large, ITMAX too small in gcf"
 1    gammcf = EXP(-x + a*LOG(x) - gln) * h
      RETURN
      END SUBROUTINE gcf

!-----------------------------------------------------------------------
! Log Gamma Function (ln(Gamma(x)))
! Lanczos approximation
!-----------------------------------------------------------------------
      FUNCTION gammaln(xx)
      IMPLICIT NONE
      DOUBLE PRECISION gammaln, xx
      DOUBLE PRECISION x, y, tmp, ser
      DOUBLE PRECISION cof(6)
      INTEGER j
      
      DATA cof /76.18009172947146d0, -86.50532032941677d0, &
     &          24.01409824083091d0, -1.231739572450155d0, &
     &          0.1208650973866179d-2, -0.5395239384953d-5/

      x = xx
      y = x
      tmp = x + 5.5d0
      tmp = tmp - (x + 0.5d0) * LOG(tmp)
      ser = 1.000000000190015d0
      DO j = 1, 6
         y = y + 1.0d0
         ser = ser + cof(j) / y
      END DO
      gammaln = -tmp + LOG(2.5066282746310005d0 * ser / x)
      RETURN
      END FUNCTION gammaln