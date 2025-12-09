**==xstitg.spg  processed by SPAG 4.50J  at 15:47 on 21 Sep 1995
      SUBROUTINE XSTITG(Ear,Npts,Param,Ifl,Photar,Photer)      
      IMPLICIT NONE
 
      INTEGER Npts , Ifl
      REAL Ear(0:Npts) , Param(5) , Photar(Npts), Photer(Npts)
 
c relativistic compton scattering of wien-law spectrum
c CONVERTED FROM TITARCHUK'S STAND-ALONE CODE. HANDLES OPTICALLY
C THIN AND THICK CASES. ADAPTED FOR XSPEC January 1995, T. Yaqoob.
C **       IMPORTANT **
C      You MUST consult Titarchuk 1994 (ApJ, 434, 313)
C       and HUA & Titarchuk 1994 (ApJ xxx,xxx) to fully understand
C      and appreciate the physical assumptions and approximations.
C       This routine will not warn you when the physical assumptions
C      break down.
C---
C---

      DOUBLE PRECISION aa , z , pi , gam0 , YYIT2 , argg , rdel , taue , 
     &                 xr , xt
      DOUBLE PRECISION x , temp , t5 , beta , x0 , ro
      DOUBLE PRECISION fx , alfax0 , algx0 , al3x0 , bol3
      DOUBLE PRECISION taup , tb1 , tb2 , t23 , bol2i , bol7i , bol6
      DOUBLE PRECISION COMPD0 , tal , tal0 , bola , BETAINT
      REAL GAMMI , GAMLN , bol3g , arg1 , dumt
      DOUBLE PRECISION xrfac , tfx , f0theta
      INTEGER i , ii
      REAL zfac , t0 , apprx
      REAL ens, oens, phot, ophot 

c=====>
C 1. Z    ** REDSHIFT
C 2. T0   ** WIEN TEMPERATURE (keV)
C 3. TEMP ** PLASMA TEMPERATURE (KEV)
C 4. TAUP ** PLASMA OPTICAL DEPTH
C 5. APPRX * ABS(APPRX) <= 1.0 DISK
C           ABS(APPRX) >  1.0 SPHERE
C           APPRX < 0 : Get beta as a function of optical depth by
C                   by interpolation of accurately calculated points
C                   in Sunyaev & Titarchuk 1985 (A&A, 143, 374).
C           APPRX >= 0: Get beta as a function of optical depth from
C                   analytic approximation (e.g. Titarchuk 1994).
      DATA pi/3.1415927/

c suppress a warning message from the compiler
      i = ifl

c this model has no errors
      DO i = 1, Npts
         photer(i) = 0.0
      ENDDO

      ro = 1.
      bol2i = 0.0D0
      temp = 0.0D0
      taup = 0.0D0
      zfac = 1. + Param(1)
      t0 = Param(2)
      temp = Param(3)
      dumt = Param(4)
      taup = dumt
      apprx = Param(5)
      t5 = temp/511.
c dimensionless soft photon energy
      x0 = 3.*t0/temp
      t23 = taup + (2./3.)
      taue = taup/(1.0+(temp/39.2)**0.86)
      rdel = 0.0
      beta = 0.0
      IF ( ABS(apprx).GT.1.0 ) THEN
c BETA FOR SPHERE
         IF ( taup.LE.0.1D0 ) THEN
            beta = DLOG(1./0.75/taup)
         ELSEIF ( taup.GE.10.D0 .OR. apprx.GE.0.0 ) THEN
            tb1 = pi*pi*(1.-EXP(-0.7*taup))/3./t23/t23
            tb2 = DEXP(-1.4*taup)*DLOG(4./(3.*taup))
            beta = tb1 + tb2
         ELSEIF ( taup.GT.0.1D0 .AND. taup.LT.10.D0 .AND. apprx.LT.0.0 )
     &            THEN
            beta = BETAINT(taup,apprx)
         ENDIF
c                 write(*,*) 'tau & beta = ',taup,beta
         IF ( taue.LT.0.01D0 ) THEN
            rdel = taue/2.D0
         ELSE
            rdel = 1.0 - 3./taue*(1.-2./taue+2./taue/taue*(1-EXP(-taue))
     &             )
         ENDIF
      ELSEIF ( ABS(apprx).LE.1.0 ) THEN
c BETA FOR DISK
         IF ( taup.LE.0.1D0 ) THEN
            beta = DLOG(1./taup/DLOG(1.53/taup))
         ELSEIF ( taup.GE.10.D0 .OR. apprx.GE.0.0 ) THEN
            tb1 = pi*pi*(1.-EXP(-1.35*taup))/12./t23/t23
            tb2 = 0.45*DEXP(-3.7*taup)*DLOG(10./(3.*taup))
            beta = tb1 + tb2
         ELSEIF ( taup.GT.0.1D0 .AND. taup.LT.10.D0 .AND. apprx.LT.0.0 )
     &            THEN
            beta = BETAINT(taup,apprx)
         ENDIF
c       write(*,*) 'tau beta = ',taup,beta
         IF ( taue.LT.0.01D0 ) THEN
            rdel = taue/2.D0
         ELSE
            rdel = 1.0 - (1.0-EXP(-taue))/taue
         ENDIF
      ENDIF
      f0theta = 2.5*t5 + 1.875*t5*t5*(1.0-t5)
      gam0 = beta/t5
c        write(*,*) 'gam0 = ',gam0
c compute the power-law index
      tal0 = DSQRT(2.25+gam0) - 1.5D0

C     These changes by CM 22 Mar 2001 to improve the convergence of the
C     solution.  Now the loop goes for a maximum of 50 iterations, and
C     the convergence criterium is looser.  Finally, the corrected
C     expression using COMPD0 has been replaced.

      if (tal0 .GT. 10.)   tal0 = 10.

      DO ii = 1, 50
         bola = 1.D0 + (tal0+3.D0)*t5/(1.D0+t5)
     &               + 4.*COMPD0(tal0)**(1./tal0)*t5*t5
         tal = beta/DLOG(bola)
         IF ( DABS(tal-tal0).LE.1.0D-4 ) GOTO 100
         tal0 = tal
      ENDDO

 100  CONTINUE
      alfax0 = tal
      if (alfax0 .GT. 10.)   alfax0 = 10.

c follwoing is outdated
c      IF (ithick.eq.1) then
c        gamx0=gam0/(1.+f0theta)
c      elseif (ithick.eq.0) then
c        gamx0=gam0*(8.+15.*t5)/(8.+19.*t5)
c      endif
      aa = 3./x0
c      alfax0 = SQRT(2.25+Gamx0) - 1.5
      argg = 3. + alfax0
      arg1 = SNGL(argg)
      xrfac = rdel**(1./argg)
      algx0 = 2.0*alfax0 + 3.0
      al3x0 = alfax0*(alfax0+3.)/2./algx0

      ens = 0.
      oens = 0.
      phot = 0.
      ophot = 0.

      DO 200 i = 0, Npts
         ens = Ear(i)*zfac
         IF ( ens.GT.0.0 ) THEN
            xt = ens/temp
            x = xt
            z = xt*t5
c            write(*,*) 'z x =',z,x
            xr = x*xrfac
c            write(*,*) 'xr = ',xr
            bol3 = x*aa
c            write(*,*) 'bol3 = ',bol3
            bol3g = SNGL(bol3)
c             write(*,*) 'bol3g = ',bol3g
c see if fx is going to overflow
            tfx = (-1.-alfax0)*DLOG10(bol3) - x*0.4343 + DLOG10(algx0)
c      write(*,*) 'tfx = ',tfx
            IF ( DABS(tfx).LT.27. ) THEN
               fx = bol3**(-alfax0-1.)*DEXP(-x)*algx0
            ELSE
               fx = 0.0
            ENDIF
c           write(*,*) 'arg1 bol3g gammi = ',arg1,bol3g,gammi(arg1,bol3g)
c         write(*,*) ' gamln exp(gamln)= ',gamln(arg1),exp(gamln(arg1))
            bol2i = GAMMI(arg1,bol3g)*EXP(GAMLN(arg1))
c         write(*,*) 'bol2i =',bol2i
c         write(*,*) 'fx yyit2',fx,yyit2(xr,alfax0,ro)
            bol7i = YYIT2(xr,alfax0,ro)*fx*bol2i
c         write(*,*) 'bol7i = ',bol7i
            bol6 = bol3**2.*DEXP(-bol3)/(alfax0+3.D0)
            phot = SNGL(al3x0*(bol7i+bol6)*aa)
c      if (i.eq.0) then
c        write(*,*) 'beta = ',beta
c        write(*,*) 'f0theta = ',f0theta
c        write(*,*) 'gamx0 = ',gamx0
c        write(*,*) 'alfax0 = ',alfax0
c        write(*,*) 'argg = ',argg
c        write(*,*) 'bol2i = ',bol2i
c        write(*,*) 'x = ',x
c        write(*,*) 'bol3 = ',bol3
c        write(*,*) 'fx = ',fx
c        write(*,*) 'bol7i bol6 ',bol7i,bol6
c        write(*,*) 'gamma(arg) ',exp(gamln(arg1))
c        write(*,*) 'igamma   ',gammi(arg1,bol3g)
c      endif
         ELSE
            phot = 0.0
         ENDIF
c            write(2,*) ens(i),bol7i,bol6,bol6+bol7i
         if (i .GT. 0) then
            photar(i) = 0.5*(phot+ophot)*(ens-oens)/zfac
         endif

         oens = ens
         ophot = phot
 200  CONTINUE

c====>
      END
**==yyit2.spg  processed by SPAG 4.50J  at 15:47 on 21 Sep 1995
c
c************************************************************************
C ====>
      DOUBLE PRECISION FUNCTION YYIT2(X,Alfa,Ro)
      IMPLICIT NONE
C*** Start of declarations inserted by SPAG
      DOUBLE PRECISION a2 , a3 , v , w , X , z , Ro
      INTEGER i
C*** End of declarations inserted by SPAG
      DOUBLE PRECISION Alfa , db , a2n , v1
      REAL GAMLN , argg
      DIMENSION w(10) , z(10)
      DATA z/.1377934705 , .7294545495 , 1.8083429017 , 3.4014336978 , 
     &     5.55249614 , 8.3301527467 , 11.8437858379 , 16.2792578313 , 
     &     21.9965858119 , 29.9206970122/
      DATA w/3.0844111576D-01 , 4.0111992915D-01 , 2.180682876E-01 , 
     &     6.208745609E-02 , 9.501516975E-03 , 7.530083885E-04 , 
     &     2.825923349E-05 , 4.249313984E-07 , 1.839564823E-09 , 
     &     9.911827219E-13/
      a2 = Alfa + 3.
      a2n = Alfa + 2.
      a3 = Alfa - 1.
      argg = SNGL(a2 + a3 + 2.)
      db = GAMLN(argg)
      YYIT2 = 0.D0
      DO 100 i = 1 , 10
         v = w(i)*DEXP(a2n*DLOG(Ro*X+z(i))+Alfa*DLOG(z(i))-db)
c          write(*,*) 'alfa x z a2 ',alfa,x,z(i),a2
c          write(*,*) alfa*((x+z(i)) - a2)
         v1 = v/Alfa*((X+z(i))-a2)
c            v = w(i)*DEXP(a2*DLOG(ro*x+z(i))+a3*Dlog(z(i))-db)
         YYIT2 = YYIT2 + v1
 100  CONTINUE
      YYIT2 = YYIT2
      RETURN
      END
**==gamln.spg  processed by SPAG 4.50J  at 15:47 on 21 Sep 1995
 
C====>
      REAL FUNCTION GAMLN(Az)
      IMPLICIT NONE
      DOUBLE PRECISION z , a(7) , s , pi
      REAL Az
      INTEGER i
      DATA pi/3.1415926536/
      a(1) = 1./12.
      a(2) = 1./30.
      a(3) = 53./210.
      a(4) = 195./371.
      a(5) = 22999./22737.
      a(6) = 29944523./19733142.
      a(7) = 109535241009./48264275462.
      z = 0.D0
      z = Az
c      write(*,*) 'gamln : az  z = ',az,z
      s = z
 
      DO 100 i = 1 , 6
         s = z + a(8-i)/s
 100  CONTINUE
 
      s = a(1)/s
      s = s - z + (z-0.5)*DLOG(z) + 0.5*DLOG(2.*pi)
      GAMLN = SNGL(s)
      RETURN
      END
**==gammi.spg  processed by SPAG 4.50J  at 15:47 on 21 Sep 1995
c***********************************************************
      FUNCTION GAMMI(A,X)
      use fgsl
      IMPLICIT NONE
      REAL A , X , GAMMI
      DOUBLE PRECISION a_doub, x_doub
      IF ( X.LT.0.0 .OR. A.LE.0.0 ) 
     & CALL XWRITE('Inc. Gamma Fn. called with x < 0 or A <= 0', 10)
      a_doub = A
      x_doub = X
      GAMMI = SNGL(fgsl_sf_gamma_inc_p(a_doub,x_doub))
      RETURN
      END
**==d0.spg  processed by SPAG 4.50J  at 15:47 on 21 Sep 1995
c**************************************************************
      DOUBLE PRECISION FUNCTION COMPD0(X)
      IMPLICIT NONE
      DOUBLE PRECISION X , x1 , x2 , x3 , bol , db
      REAL GAMLN , argg
      x1 = X + 1.D0
      x2 = X + 2.D0
      x3 = X + 3.D0
      bol = x3*X + 4.D0
      argg = SNGL(2.D0*x1)
      db = GAMLN(argg)
      IF ( db.LT.-70. ) THEN
         COMPD0 = 0.0
      ELSEIF ( db.GT.70. ) THEN
         COMPD0 = 1.E34
      ELSE
         COMPD0 = 3.D0*bol*DEXP(db)/x3/x2/x2
      ENDIF
      RETURN
      END
**==betaint.spg  processed by SPAG 4.50J  at 15:47 on 21 Sep 1995
c*****************************************************************
      DOUBLE PRECISION FUNCTION BETAINT(Tau,Apprx)
      IMPLICIT NONE
c Get beta parameter as a function of optical depth for sphere or
c disk geometry by interpolatiing accurately calculated points in
c Sunyaev & Titrachuk (1985), A&A 143, 374. Optical depth must be
c inside the range 0.1-10.0.
c DISK: abs(apprx) le 1.0
C SPHERE: abs(apprx) gt 1.0
      INTEGER NPTS , ix
      PARAMETER (NPTS=99)
      REAL Apprx
      DOUBLE PRECISION sbetal(NPTS) , sbetah(NPTS) , taul(NPTS) , 
     &                 tauh(NPTS)
      DOUBLE PRECISION dbetal(NPTS) , dbetah(NPTS) , Tau
      DOUBLE PRECISION tau1 , tau2 , b1 , b2 , b3
c
      DATA taul/0.1D0 , 0.2D0 , 2*0.3D0 , 2*0.5D0 , 3*0.7D0 , 5*1.0D0 , 
     &     5*1.5D0 , 5*2.0D0 , 5*2.5D0 , 10*3.0D0 , 10*4.D0 , 20*5.0D0 , 
     &     30*7.D0/
      DATA tauh/0.2D0 , 0.3D0 , 2*0.5D0 , 2*0.7D0 , 3*1.0D0 , 5*1.5D0 , 
     &     5*2.0D0 , 5*2.5D0 , 5*3.0D0 , 10*4.0D0 , 10*5.D0 , 20*7.0D0 , 
     &     30*10.D0/
      DATA sbetal/2.97D0 , 1.97D0 , 2*1.60D0 , 2*1.18D0 , 3*0.925D0 , 
     &     5*0.688D0 , 5*0.462D0 , 5*0.334D0 , 5*0.252D0 , 10*0.197D0 , 
     &     10*0.130D0 , 20*0.092D0 , 30*0.0522D0/
      DATA sbetah/1.97D0 , 1.60D0 , 2*1.18D0 , 2*0.925D0 , 3*0.688D0 , 
     &     5*0.462D0 , 5*0.334D0 , 5*0.252D0 , 5*0.197D0 , 10*0.130D0 , 
     &     10*0.092D0 , 20*0.0522D0 , 30*0.0278D0/
      DATA dbetal/1.23D0 , 0.84D0 , 2*0.655D0 , 2*0.45D0 , 3*0.335D0 , 
     &     5*0.234D0 , 5*0.147D0 , 5*0.101D0 , 5*0.0735D0 , 
     &     10*0.0559D0 , 10*0.0355D0 , 20*0.0245D0 , 30*0.013393D0/
      DATA dbetah/0.84D0 , 0.655D0 , 2*0.45D0 , 2*0.335D0 , 3*0.234D0 , 
     &     5*0.147D0 , 5*0.101D0 , 5*0.0735D0 , 5*0.0559D0 , 
     &     10*0.0355D0 , 10*0.0245D0 , 20*0.013393D0 , 30*0.00706D0/
 
c get look-up index
      ix = INT(Tau*10.)
      tau1 = taul(ix)
      tau2 = tauh(ix)
      IF ( ABS(Apprx).GT.1.0 ) THEN
         b1 = sbetal(ix)
         b2 = sbetah(ix)
      ELSEIF ( ABS(Apprx).LE.1.0 ) THEN
         b1 = dbetal(ix)
         b2 = dbetah(ix)
      ENDIF
c      write(*,*) 'tau = ',tau
      CALL DINTER(Tau,tau1,tau2,b3,b1,b2)
      BETAINT = b3
c      write(*,*) tau1,tau2,b1,b2,b3
      RETURN
      END
**==dinter.spg  processed by SPAG 4.50J  at 15:47 on 21 Sep 1995
C *******************************************************************
C Interpolate in log- space
 
      SUBROUTINE DINTER(X0,X1,X2,Y0,Y1,Y2)
      IMPLICIT NONE
      DOUBLE PRECISION X0 , X1 , X2 , Y0 , Y1 , Y2 , tn , tn1 , t0
 
      IF ( X1.NE.X2 ) THEN
 
 
         IF ( Y1.GT.0.0D0 ) THEN
 
            tn = -LOG10(Y1)
 
         ELSE
 
            tn = 36.0D0
 
         ENDIF
 
         IF ( Y2.GT.0.0D0 ) THEN
 
            tn1 = -LOG10(Y2)
 
         ELSE
 
            tn1 = 36.0D0
 
         ENDIF
 
         t0 = (tn*(X2-X0)+tn1*(X0-X1))/(X2-X1)
 
         IF ( tn.GE.36.D0 .AND. tn1.GE.36.D0 ) THEN
 
            Y0 = 0.0D0
 
         ELSE
 
            Y0 = 10.0D0**(-t0)
 
         ENDIF
 
      ELSE
 
         Y0 = Y1
 
      ENDIF
 
      END
 
