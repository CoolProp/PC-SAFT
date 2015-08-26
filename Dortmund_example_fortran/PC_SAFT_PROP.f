      PROGRAM STAND_ALONE

      IMPLICIT NONE

c---------------------------------------------------------------------
c
c  This first SOUBROUTINE calculates phase equilibria (VLE and LLE)
c  of binary systems. It calles a SUBROUTINE delivering physical
c  properties of Perturbed-Chain SAFT equation of state (PC_SAFT_PROP).
c  Generally, this fist routine it is the most simple and stupid
c  algorithm, the autor could think of.
c
c  Please bear in mind, what this file is intended for: it is written,
c  to provide a reference with the help of which, other programs can
c  be tested. By no means is this code numerically robust. The
c  SUBROUTINES delivering physical properties of the PC-SAFT equation
c  of state, however, may be used within a more intelligent software
c  environment.
c
c  Clearly, the autor does not take any liability for the SOUBROUTINES
c  included in this file.
c
c  If you have questions, please feel free to contact me at:
c  joachim.gross@alumni.tu-berlin.de
c
c             Joachim Gross.
c
c---------------------------------------------------------------------
c
      INTEGER nc, nph
      PARAMETER (nc=20,nph=2)
      INTEGER ncomp, nphas
      DOUBLE PRECISION phi(nph,nc)
      DOUBLE PRECISION t,p,densys(nph),x_sys(nph,nc)
      DOUBLE PRECISION parame(nc,5),kij(nc,nc)
      DOUBLE PRECISION d_sta(nph)

      INTEGER i,j, iso, number
      DOUBLE PRECISION resid1,resid2, x_new(nph,nc),fvek1,fvek2,x_save,
     &                 df1dx,df2dx,df1dy,df2dy,betrag,dens(nph),end_pt,
     &                 delta

      COMMON /condit/ t,p,densys,x_sys,ncomp,nphas
      COMMON /system/ parame,kij
      COMMON /staval/ d_sta

      ncomp = 2
      nphas = 2
      number = 11

      OPEN (45,FILE='output.txt')
      write (45,*) '                     x(1) ,   t / �C  ,   p / bar'
      write (45,*) ' '

      OPEN (85,FILE='parame.inp')
      READ (85,*) t
      READ (85,*) p
      DO i=1,ncomp
        DO j=1,3
          READ (85,*) parame(i,j)
        END DO
      END DO
      READ (85,*) kij(1,2)
      kij(2,1)=kij(1,2)
      READ (85,*) x_new(1,1),x_new(2,1)
      READ (85,*) d_sta(1),d_sta(2)
      CLOSE (85)
      t = t + 273.15d0
      p = p * 1.d5

      call TEST()
      RETURN

      write (*,*) ' '
      write (*,*) ' press   1 : for calculation of isotherms'
      write (*,*) ' press   2 : for calculation of isobars'
      read (*,*) iso
      write (*,*) ' enter the end point for the calculation'
      IF (iso.EQ.1) THEN
        write (*,*) ' end pressure in [bar]'
        read (*,*) end_pt
        delta = (end_pt*1d5 - p)/DBLE(number-1)
        p = p - delta
      ELSE
        write (*,*) ' end temperature in [�C]'
        read (*,*) end_pt
        delta = (end_pt+273.15d0 - t)/DBLE(number-1)
        t = t - delta
      ENDIF
      write (*,*) ' '


      DO 2 j = 1 , number

      IF (iso.EQ.1) p = p + delta      
      IF (iso.EQ.2) t = t + delta      

 5    CONTINUE

        x_sys(1,1) = x_new(1,1)
        x_sys(1,2) = 1.d0 - x_sys(1,1)
        x_sys(2,1) = x_new(2,1)
        x_sys(2,2) = 1.d0 - x_sys(2,1)

        CALL PC_SAFT_PROP (phi)
        CALL ZIELFKT (phi,x_sys,resid1,resid2)

        dens(1)=densys(1)
        dens(2)=densys(2)

        fvek1= resid1
        fvek2= resid2

c------numerical derivative to x_sys(1,1)-----------------------------
        x_save = x_sys(1,1)
        x_sys(1,1) = x_save*1.0001d0
        x_sys(1,2) = 1.d0 - x_sys(1,1)

        CALL PC_SAFT_PROP (phi)
        CALL ZIELFKT (phi,x_sys,resid1,resid2)

        df1dx = (fvek1-resid1)/( x_save-x_sys(1,1) )
        df2dx = (fvek2-resid2)/( x_save-x_sys(1,1) )
        x_sys(1,1)=x_save
        x_sys(1,2) = 1.d0 - x_sys(1,1)
        


c------numerical derivative to x_sys(2,1)-----------------------------
        x_save = x_sys(2,1)
        x_sys(2,1) = x_save*1.0001d0
        x_sys(2,2) = 1.d0 - x_sys(2,1)

        CALL PC_SAFT_PROP (phi)
        CALL ZIELFKT (phi,x_sys,resid1,resid2)

        df1dy = (fvek1-resid1)/( x_save-x_sys(2,1) )
        df2dy = (fvek2-resid2)/( x_save-x_sys(2,1) )

        x_sys(2,1) = x_save
        x_sys(2,2) = 1.d0 - x_sys(2,1)

        betrag = df1dx*df2dy - df1dy*df2dx
        x_new(1,1)= x_sys(1,1) -1.d0/betrag*(df2dy*fvek1+(-df1dy)*fvek2)
        x_new(2,1)= x_sys(2,1) -1.d0/betrag*((-df2dx)*fvek1+df1dx*fvek2)

      IF ((DABS(fvek1).GT.1.d-7).OR.(DABS(fvek2).GT.1.d-7)) GOTO 5

      x_sys(1,1) = x_new(1,1)
      x_sys(1,2) = 1.d0 - x_sys(1,1)
      x_sys(2,1) = x_new(2,1)
      x_sys(2,2) = 1.d0 - x_sys(2,1)

      IF ((DABS(fvek1)+DABS(fvek2)).LT.1.d-6) THEN
        IF ( DABS(x_sys(1,1)-x_sys(2,1)).LT.1.d-4) THEN
          write (*,*) ' homogeneous phase ??'
          write (*,*) '  at T=',t-273.15d0,' degC  and p=',p/1.d5,' bar'
          write (45,*) ' homogeneous phase ??'
          write (45,'(a9,e12.6,a11,e12.6,a4)')
     &              '    at T=',t-273.15d0,' degC  and p=',p/1.d5,' bar'
        ELSE
          write (*,'(a18,3(1x,f12.8))') ' phase 1, comp. 1 ',
     &                                    x_sys(1,1),t-273.15d0,p/1.d5
          write (*,'(a18,3(1x,f12.8))') ' phase 2, comp. 1 ',
     &                                    x_sys(2,1),t-273.15d0,p/1.d5
          write (45,'(a18,3(2x,e12.6))')
     &              ' phase 1, comp. 1 ',x_sys(1,1),t-273.15d0,p/1.d5
          write (45,'(a18,3(2x,e12.6))')
     &              ' phase 2, comp. 1 ',x_sys(2,1),t-273.15d0,p/1.d5
        ENDIF
      ELSE
c        write (*,*) fvek1,fvek2
        write (*,*) ' no convergence at T=',t-273.15d0,' and p=',p/1.d5
        write (45,'(a21,e12.6,a7,e12.6)')
     &            ' no convergence at T=',t-273.15d0,' and p=',p/1.d5
      ENDIF
      write (*,*) ' '
      write (45,*) ' '

 2    CONTINUE
      CLOSE (45)

      END

      SUBROUTINE TEST()

      INTEGER nc, nph
      PARAMETER (nc=20, nph=2)
      INTEGER ncomp

      DOUBLE PRECISION kij(nc,nc),parame(nc,25)
      DOUBLE PRECISION x(nc),t,p,mseg(nc)
      DOUBLE PRECISION densta,dense,dap_dx(nc,7),dbp_dx(nc,7)
      DOUBLE PRECISION order1,order2,apar(7),bpar(7)
      DOUBLE PRECISION z0t,z1t,z2t,z3t,z0,z1,z2,z3
      DOUBLE PRECISION PI, RGAS, NA, KBOL, TAU
      DOUBLE PRECISION dij_ab(nc,nc),uij(nc,nc),sig_ij(nc,nc),sig_t(nc)
      DOUBLE PRECISION d_sta(nph),x_new(nph,nc)


      ncomp = 2

      OPEN (85,FILE='parame.inp')
      READ (85,*) t
      READ (85,*) p
      DO i=1,ncomp
        DO j=1,3
          READ (85,*) parame(i,j)
        END DO
      END DO
      READ (85,*) kij(1,2)
      kij(2,1)=kij(1,2)
      READ (85,*) x_new(1,1),x_new(2,1)
      READ (85,*) d_sta(1),d_sta(2)
      CLOSE (85)

      t = 300 ![K]
      dense = 0.00099912374 ![atoms/A^3]
      x(1) = 0.4
      x(2) = 0.6

      ! Get parameters
      CALL PERTPAR (kij,parame,
     &     ncomp,x,t,p,mseg,densta,dense,dap_dx,dbp_dx,
     &     order1,order2,apar,bpar,z0t,z1t,z2t,z3t,dij_ab,
     &     PI,RGAS,NA,KBOL,TAU,sig_t,uij,sig_ij)

      ! Calculate the pressure
      CALL P_EOS (pges,pgesdz,gij,
     &     ncomp,x,t,p,mseg,densta,dense,dap_dx,dbp_dx,
     &     order1,order2,apar,bpar,z0t,z1t,z2t,z3t,dij_ab,
     &     PI,RGAS,NA,KBOL,TAU,sig_t,uij,sig_ij)

      END



c*********************************************************************
c*********************************************************************
c*********************************************************************
      SUBROUTINE ZIELFKT (x,phi,resid1,resid2)

      IMPLICIT NONE

c---------------------------------------------------------------------
c
c  This SOUBROUTINE evaluates the objective functions (resid1 and
c  resid2). In the course of the outer iteration-loop, these objective
c  functions are convereged to be zero.
c
      INTEGER nc, nph
      PARAMETER (nc=20,nph=2)
      DOUBLE PRECISION phi(nph,nc),x(nph,nc)
      DOUBLE PRECISION resid1,resid2


      resid1 = x(1,1)*phi(1,1) - x(2,1)*phi(2,1)
      resid2 = x(1,2)*phi(1,2) - x(2,2)*phi(2,2)

      RETURN
      END
c*********************************************************************







c---------------------------------------------------------------------


      SUBROUTINE PC_SAFT_PROP (phi)

      IMPLICIT NONE

c---------------------------------------------------------------------
c
c  The following SOUBROUTINES deliver physical properties as
c  calculated from the Perturbed-Chain (PC-SAFT) equation of state.
c  Obviously, the autor (Joachim Gross) does not take any liability
c  for this code.
c
c  If you encounter difficulties with these routines, or if you require
c  further derivatives of the PC-SAFT equation of state, please feel
c  free to contact me at
c    joachim.gross@alumni.tu-berlin.de
c
c  The equations coded here were published by
c    Joachim Gross and Gabriele Sadowski, Perturbed-Chain SAFT: An
c    Equation of State Based on a Perturbation Theory for Chain Molecules
c    Ind. Eng. Chem. Res. 2001, 40, 1244-1260.
c  One misprint was found so far: The right-hand side of equation (A.11)
c  should be raised to the power of -1. It is then consistent with eqn(13)
c  and eqn(25).
c
c
c  If you would like to insert this file into your phase equilibrium
c  software, please define your variables here:
c
c  required are:
c
c      the input variables:
c
c      ncomp   :  number of components
c      nphas   :  number of phases
c      t       :  temperature
c      p       :  pressure
c      x_sys   :  mole fractions of all components in al phases
c      kij     :  binary interaction correction
c      d_sta   :  starting value for the densities in all phases
c      parame  :  pure component parameters,
c                 where parame(i,1) is the segment number of component i
c                       parame(i,1) is the segment diameter of comp. i
c                       parame(i,1) is the energy parameter of comp. i
c
c      and the output variables:
c
c      phi     :  fugacity coefficient in all phases ( phi = f(t,p,x) )
c      densys  :  iteration result for the density in all phases
c
c
c
c   it may for examle be defined:
c
c
      INTEGER nc, nph
      PARAMETER (nc=20,nph=2)
      INTEGER ncomp, nphas
      DOUBLE PRECISION t,p,densys(nph),x_sys(nph,nc)
      DOUBLE PRECISION parame(nc,5),kij(nc,nc)
      DOUBLE PRECISION enthal(nph),entrop(nph),gibbs(nph)
      DOUBLE PRECISION d_sta(nph)

      COMMON /condit/ t,p,densys,x_sys,ncomp,nphas
      COMMON /system/ parame,kij
      COMMON /calor/  enthal,entrop,gibbs
      COMMON /staval/ d_sta
c
c   If in your program, 'x_mole' is used for mole fractions instead
c   of 'x_sys' as in the examle given above, then all 'x_sys' have
c   to be replaced by 'x_mole' in this SUBROUTINE.
c
c---------------------------------------------------------------------

c-----local variables-------------------------------------------------
      INTEGER i,ph
      DOUBLE PRECISION phi(nph,nc)
      DOUBLE PRECISION xtrans(nc),dstart,parsys(nc,25),fugcoe(nc),den
      DOUBLE PRECISION h_res,s_res,g_res
c---------------------------------------------------------------------

c---------------------------------------------------------------------
      DO 1 ph = 1,nphas
        dstart = d_sta(ph)
        DO 11 i = 1,ncomp
          xtrans(i)   = x_sys(ph,i) 
          parsys(i,1) = parame(i,1)
          parsys(i,2) = parame(i,2)
          parsys(i,3) = parame(i,3)
 11     CONTINUE

        CALL PHIEOS (fugcoe,xtrans,t,p,parsys,kij,ncomp,dstart,den)
        densys(ph)  = den
        DO 12 i = 1,ncomp
          phi(ph,i) = fugcoe(i)
 12     CONTINUE
        enthal(ph) = h_res
        entrop(ph) = s_res
        gibbs(ph)  = g_res
 1    CONTINUE

      RETURN
      END

c*********************************************************************
c*********************************************************************
c*********************************************************************
      SUBROUTINE PHIEOS (phi,x,t,p,parame,kij,ncomp,densta,dense)

      IMPLICIT NONE

c-----variables used in the parameter list of subroutine--------------
      INTEGER nc
      PARAMETER (nc=20)
      INTEGER ncomp

      DOUBLE PRECISION phi(nc)
      DOUBLE PRECISION kij(nc,nc),parame(nc,25)
      DOUBLE PRECISION h_res,s_res,g_res
      DOUBLE PRECISION pges,pgesdz,gij(nc, nc)
      DOUBLE PRECISION fres
      DOUBLE PRECISION x(nc),t,p,mseg(nc)
      DOUBLE PRECISION densta,dense,dap_dx(nc,7),dbp_dx(nc,7)
      DOUBLE PRECISION order1,order2,apar(7),bpar(7)
      DOUBLE PRECISION z0t,z1t,z2t,z3t,z0,z1,z2,z3
      DOUBLE PRECISION PI, RGAS, NA, KBOL, TAU
      DOUBLE PRECISION dij_ab(nc,nc),uij(nc,nc),sig_ij(nc,nc),sig_t(nc)
c---------------------------------------------------------------------

c-----local variables-------------------------------------------------
      INTEGER i, j, k, m
      DOUBLE PRECISION zms, rho, m_mean
      DOUBLE PRECISION mhs(nc), mhc(nc), mdsp(nc), mpart(nc),
     &                 myres(nc), myresq, lnphi(nc)
      DOUBLE PRECISION dgijdx(nc, nc, nc)
      DOUBLE PRECISION zres, zges
      DOUBLE PRECISION fhs_sg, fhs_sx
      DOUBLE PRECISION z0dx,z1dx,z2dx,z3dx,m_mndx(nc)
      DOUBLE PRECISION I1, I2, I1_dx, I2_dx, c1_con,c2_con,c1_dx,
     &                 ord1dx, ord2dx
c---------------------------------------------------------------------

c------obtain parameters and density independent expressions----------
      CALL PERTPAR (kij,parame,
     1     ncomp,x,t,p,mseg,densta,dense,dap_dx,dbp_dx,
     2     order1,order2,apar,bpar,z0t,z1t,z2t,z3t,dij_ab,
     3     PI,RGAS,NA,KBOL,TAU,sig_t,uij,sig_ij)

c------density iteration----------------------------------------------
       CALL DENSITR (pges,pgesdz,gij,
     1     ncomp,x,t,p,mseg,densta,dense,dap_dx,dbp_dx,
     2     order1,order2,apar,bpar,z0t,z1t,z2t,z3t,dij_ab,
     3     PI,RGAS,NA,KBOL,TAU,sig_t,uij,sig_ij)

c-----residual Helmholtz free energy----------------------------------
      CALL F_EOS (fres,gij,
     1     ncomp,x,t,p,mseg,densta,dense,dap_dx,dbp_dx,
     2     order1,order2,apar,bpar,z0t,z1t,z2t,z3t,dij_ab,
     3     PI,RGAS,NA,KBOL,TAU,sig_t,uij,sig_ij)

c-----residual enthalpy and entropy-----------------------------------
c  From a calculation-efficiency point of view, it is not optimal to
c  call the SOUBROUTINE for caloric properties (H_res,S_res and G_res)
c  here. However, this file is not optimized for calculation efficiency.
      CALL H_EOS (h_res,s_res,g_res,fres,pges,pgesdz,gij,parame,
     1     ncomp,x,t,p,mseg,densta,dense,dap_dx,dbp_dx,
     2     order1,order2,apar,bpar,z0t,z1t,z2t,z3t,dij_ab,
     3     PI,RGAS,NA,KBOL,TAU,sig_t,uij,sig_ij)

      rho = dense/z3t
      z0 = z0t*rho
      z1 = z1t*rho
      z2 = z2t*rho
      z3 = z3t*rho

      zms    = 1.d0 - dense
      m_mean = z0t/(PI/6.d0)

c-----compressibility factor z = p/(kT*rho)---------------------------
      zges = (pges * 1.d-30)/(KBOL*t*rho)
      zres = zges - 1.d0

c-----calcul. the derivatives of f to mole fraction x ( d(f)/d(x) )---
      DO 1 k = 1,ncomp

c-----derivative d(zeta(i))/d(x)---------------------------------------
      z0dx = rho*PI/6.d0*mseg(k)
      z1dx = rho*PI/6.d0*mseg(k)*sig_t(k)
      z2dx = rho*PI/6.d0*mseg(k)*sig_t(k)*sig_t(k)
      z3dx = rho*PI/6.d0*mseg(k)*sig_t(k)**3.d0

c-----derivative d(m_mean)/d(x)---------------------------------------
      m_mndx(k) = mseg(k)

c-------d(f)/d(x) : hard sphere contribution--------------------------
      fhs_sg = (  3.d0*z1*z2/zms + z2**3.d0/z3/zms/zms
     &        + (z2**3.d0/z3/z3-z0)*DLOG(zms)  )/z0
      fhs_sx = - z0dx/z0*fhs_sg
     & +(  3.d0*(z1dx*z2+z1*z2dx)/zms + 3.d0*z1*z2*z3dx/zms/zms
     &     + 3.d0*z2*z2*z2dx/z3/zms/zms
     &     + z2**3.d0*z3dx*(3.d0*z3-1.d0)/z3/z3/zms**3.d0
     &     + ((3.d0*z2*z2*z2dx*z3-2.d0*z2**3.d0*z3dx)/z3**3.d0-z0dx)
     &                                                      *DLOG(zms)
     &     +(z0-z2**3.d0/z3/z3)*z3dx/zms  )/z0

      mhs(k) = m_mndx(k)* fhs_sg + m_mean*fhs_sx

c-------d(f)/d(x) : chain term----------------------------------------
        DO i = 1, ncomp
          DO j = 1, ncomp
            dgijdx(i,j,k) = z3dx/zms/zms
     &          +3.d0*dij_ab(i,j)*(z2dx+2.d0*z2*z3dx/zms)/zms/zms
     &          +dij_ab(i,j)**2.d0*z2/zms**3.d0
     &                              *(4.d0*z2dx+6.d0*z2*z3dx/zms)
          END DO
        END DO

        mhc(k) = 0.d0
        DO i = 1, ncomp
          mhc(k) = mhc(k) + x(i) * (1.d0-mseg(i))
     &                * (1.d0/gij(i,i)) * dgijdx(i,i,k)
        END DO
        mhc(k) = mhc(k)+( 1.d0-mseg(k))*DLOG(gij(k,k))

c-------d(f)/d(x) : dispersion contribution---------------------------
        I1    = 0.d0
        I2    = 0.d0
        I1_dx = 0.d0
        I2_dx = 0.d0
        DO m = 0,6
          I1  = I1 + apar(m+1)*dense**DBLE(m)
          I2  = I2 + bpar(m+1)*dense**DBLE(m)
          I1_dx = I1_dx + apar(m+1)*DBLE(m)*dense**DBLE(m-1)*z3dx
     &             + dap_dx(k,m+1)*dense**DBLE(m)
          I2_dx = I2_dx + bpar(m+1)*DBLE(m)*dense**DBLE(m-1)*z3dx
     &             + dbp_dx(k,m+1)*dense**DBLE(m)
        END DO

        ord1dx  = 0.d0
        ord2dx  = 0.d0
        DO 16 i = 1,ncomp
          ord1dx= ord1dx
     &  + 2.d0*mseg(k) *x(i)*mseg(i)*sig_ij(i,k)**3.d0 *uij(i,k)/t
          ord2dx= ord2dx
     &  + 2.d0*mseg(k)*x(i)*mseg(i)*sig_ij(i,k)**3.d0*(uij(i,k)/t)**2.d0
 16     CONTINUE

        c1_con= 1.d0/
     &    (  1.d0 + m_mean*(8.d0*z3-2.d0*z3*z3)/zms**4.d0
     &         + (1.d0 - m_mean)*(20.d0*z3-27.d0*z3*z3
     &                           +12.d0*z3**3.d0-2.d0*z3**4.d0)
     &                           /(zms*(2.d0-z3))**2.d0  )
        c2_con= - c1_con*c1_con
     &          *(m_mean*(-4.d0*z3*z3+20.d0*z3+8.d0)/zms**5.d0
     &        + (1.d0 - m_mean)
     &        *(2.d0*z3**3.d0+12.d0*z3*z3-48.d0*z3+40.d0)
     &                                /(zms*(2.d0-z3))**3.d0 )
        c1_dx= c2_con*z3dx - c1_con*c1_con*mseg(k)
     &        *( (8.d0*z3-2.d0*z3*z3)/zms**4.d0
     &          - (-2.d0*z3**4.d0+12.d0*z3**3.d0-27.d0*z3*z3+20.d0*z3)
     &                                         /(zms*(2.d0-z3))**2.d0 )

        mdsp(k)= -2.d0*PI*rho*(order1*I1_dx+ord1dx*I1)
     &            -     PI*rho*c1_con*m_mean*(order2*I2_dx+ord2dx*I2)
     &            -     PI*rho*(c1_con*mseg(k)+c1_dx*m_mean)*order2*I2

c-----d(f)/d(x) : summation of all contributions----------------------
        mpart(k) = mhs(k) + mhc(k) + mdsp(k)

 1    CONTINUE

      myresq = 0.d0
      DO i = 1, ncomp
        myresq = myresq - x(i)* mpart(i)
      END DO

      DO k = 1, ncomp
        myres(k) = myresq + mpart(k) + fres + zres
        lnphi(k) = myres(k) - DLOG(zges)
        phi(k)   = DEXP(lnphi(k))
      END DO

      RETURN
      END



c********************************************************************
c********************************************************************
c********************************************************************
      SUBROUTINE P_EOS (pges,pgesdz,gij,
     1     ncomp,x,t,p,mseg,densta,dense,dap_dx,dbp_dx,
     2     order1,order2,apar,bpar,z0t,z1t,z2t,z3t,dij_ab,
     3     PI,RGAS,NA,KBOL,TAU,sig_t,uij,sig_ij)

      IMPLICIT NONE

c-----variables used in the parameter list of subroutine--------------
      INTEGER nc
      PARAMETER (nc=20)
      INTEGER ncomp

      DOUBLE PRECISION h_res,s_res,g_res
      DOUBLE PRECISION pges,pgesdz,gij(nc, nc)
      DOUBLE PRECISION x(nc),t,p,mseg(nc)
      DOUBLE PRECISION densta,dense,dap_dx(nc,7),dbp_dx(nc,7)
      DOUBLE PRECISION order1,order2,apar(7),bpar(7)
      DOUBLE PRECISION z0t,z1t,z2t,z3t,z0,z1,z2,z3
      DOUBLE PRECISION PI, RGAS, NA, KBOL, TAU
      DOUBLE PRECISION dij_ab(nc,nc),uij(nc,nc),sig_ij(nc,nc),sig_t(nc)
c---------------------------------------------------------------------

c-----local variables-------------------------------------------------
      INTEGER  i, j, m
      DOUBLE PRECISION rho,zms,m_mean
      DOUBLE PRECISION zges,zhs,zhc,zdsp,zgesdz,zhsdz,zhcdz,zdspdz
 
      DOUBLE PRECISION dgijdz(nc,nc),dgijd2(nc,nc)
      DOUBLE PRECISION I2, edI1dz, edI2dz, edI1d2, edI2d2,
     &                 c1_con,c2_con,c3_con
c-------------------------------------------------------------------

      rho = dense/z3t
      z0 = z0t*rho
      z1 = z1t*rho
      z2 = z2t*rho
      z3 = z3t*rho

      zms    = 1.d0 - dense
      m_mean = z0t/(PI/6.d0)

C-----gij , the derivative dgijdz=d(gij)/d(dense) ------------------
C-----and dgijd2 = dd(gij)/d(dense)**2 -----------------------------
      DO i = 1, ncomp
c        j=i
        DO j = 1, ncomp
          gij(i,j) = 1.d0/zms + 3.d0*dij_ab(i,j)*z2/zms/zms
     &             + 2.d0*(dij_ab(i,j)*z2)**2.d0/zms**3.d0
          dgijdz(i,j)= 1.d0/zms/zms
     &             + 3.d0*dij_ab(i,j)*z2*(1.d0+z3)/z3/zms**3.d0
     &             + (dij_ab(i,j)*z2/zms/zms)**2.d0*(4.d0+2.d0*z3)/z3
          dgijd2(i,j) = 2.d0/zms**3.d0
     &             + 6.d0*dij_ab(i,j)*z2/z3/zms**4.d0*(2.d0+z3)
     &             + (2.d0*dij_ab(i,j)*z2/z3)**2.d0/zms**5.d0
     &                                          *(1.d0+4.d0*z3+z3*z3)
        END DO
      END DO


c-----p : hard sphere contribution-----------------------------------
      zhs = m_mean* ( z3/zms + 3.d0*z1*z2/z0/zms/zms
     &             + z2**3.d0/z0*(3.d0-z3)/zms**3.d0 )
      zhsdz = m_mean*(  1.d0/zms/zms
     &             + 3.d0*z1*z2/z0/z3*(1.d0+z3)/zms**3.d0
     &             + 6.d0*z2**3.d0/z0/z3/zms**4.d0 )

c-----p : chain term-------------------------------------------------
      zhc  = 0.d0
      zhcdz  = 0.d0
      DO i= 1, ncomp
        zhc = zhc + x(i)*(1.d0-mseg(i))*dense/gij(i,i)* dgijdz(i,i)
        zhcdz = zhcdz + x(i)*(1.d0-mseg(i))
     &     *(-dense*(dgijdz(i,i)/gij(i,i))**2.d0 + dgijdz(i,i)/gij(i,i)
     &          + dense/gij(i,i)*dgijd2(i,i))
      END DO

c------p : dispersion contribution-----------------------------------
c------edI1dz is equal to d(dense*I1)/d(dense)-----------------------
c------edI2dz is equal to d(dense*I2)/d(dense)-----------------------
      I2     = 0.d0
      edI1dz = 0.d0
      edI2dz = 0.d0
      edI1d2 = 0.d0
      edI2d2 = 0.d0
      DO m=0,6
        I2     = I2 + bpar(m+1)*dense**DBLE(m)
        edI1dz = edI1dz + apar(m+1)*DBLE(m+1)*dense**DBLE(m)
        edI2dz = edI2dz + bpar(m+1)*DBLE(m+1)*dense**DBLE(m)
        edI1d2 = edI1d2 + apar(m+1)*DBLE(m+1)*DBLE(m)*dense**DBLE(m-1)
        edI2d2 = edI2d2 + bpar(m+1)*DBLE(m+1)*DBLE(m)*dense**DBLE(m-1)
      END DO

      c1_con= 1.d0/
     &    (  1.d0 + m_mean*(8.d0*dense-2.d0*dense**2.d0)/ZMS**4.d0
     &         + (1.d0 - m_mean)*(20.d0*dense-27.d0*dense**2.d0
     &                           +12.d0*dense**3.d0-2.d0*dense**4.d0)
     &                           /(ZMS*(2.d0-dense))**2.d0  )
      c2_con= - c1_con*c1_con
     &          *(m_mean*(-4.d0*dense**2.d0+20.d0*dense+8.d0)/ZMS**5.d0
     &        + (1.d0 - m_mean)
     &        *(2.d0*dense**3.d0+12.d0*dense**2.d0-48.d0*dense+40.d0)
     &                                /(ZMS*(2.d0-dense))**3.d0 )
      c3_con= 2.d0 * c2_con*c2_con/c1_con
     &       - c1_con*c1_con
     &       *( m_mean*(-12.d0*dense**2.d0+72.d0*dense+60.d0)/ZMS**6.d0
     &       + (1.d0 - m_mean)
     &        *(-6.d0*dense**4.d0-48.d0*dense**3.d0+288.d0*dense**2.d0
     &                                            -480.d0*dense+264.d0)
     &                                /(ZMS*(2.d0-dense))**4.d0 )

      zdsp  = - 2.d0*PI*rho*edI1dz*order1
     &        - PI*rho*order2*m_mean*(c2_con*I2*dense + c1_con*edI2dz)
      zdspdz= zdsp/dense - 2.d0*PI*rho*edI1d2*order1
     &        - PI*rho*order2*m_mean*(c3_con*I2*dense
     &                       + 2.d0*c2_con*edI2dz + c1_con*edI2d2)

c-----p summation, p is obtained in unit [Pa] ------------------------
c      zges = 1.d0 + zhs + zhc + zdsp + zhb + zdd + zqq
      zges = 1.d0 + zhs + zhc + zdsp
      pges = zges * rho * (KBOL*t) / 1.d-30

c      zgesdz = zhsdz + zhcdz + zdspdz + zhbdz + zdddz + zqqdz
      zgesdz = zhsdz + zhcdz + zdspdz
      pgesdz = ( zgesdz*rho + zges*rho/z3 )*(kbol*T)/1.d-30

      RETURN
      END



c********************************************************************
c********************************************************************
c********************************************************************
      SUBROUTINE DENSITR (pges,pgesdz,gij,
     1     ncomp,x,t,p,mseg,densta,dense,dap_dx,dbp_dx,
     2     order1,order2,apar,bpar,z0t,z1t,z2t,z3t,dij_ab,
     3     PI,RGAS,NA,KBOL,TAU,sig_t,uij,sig_ij)

      IMPLICIT NONE

c-----variables used in the parameter list of subroutine--------------
      INTEGER nc
      PARAMETER (nc=20)
      INTEGER ncomp

      DOUBLE PRECISION pges,pgesdz,gij(nc, nc)
      DOUBLE PRECISION x(nc),t,p,mseg(nc)
      DOUBLE PRECISION densta,dense,dap_dx(nc,7),dbp_dx(nc,7)
      DOUBLE PRECISION order1,order2,apar(7),bpar(7)
      DOUBLE PRECISION z0t,z1t,z2t,z3t,z0,z1,z2,z3
      DOUBLE PRECISION PI, RGAS, NA, KBOL, TAU
      DOUBLE PRECISION dij_ab(nc,nc),uij(nc,nc),sig_ij(nc,nc),sig_t(nc)
c---------------------------------------------------------------------

c-----local variables-------------------------------------------------
      INTEGER i,start,max_i
      DOUBLE PRECISION x1,y1,dydx,acc_i
c---------------------------------------------------------------------

      acc_i = 1.d-10
      max_i = 200

      i     = 0
      x1    = densta

 1    CONTINUE

        i=i+1
        dense = x1
        CALL P_EOS (pges,pgesdz,gij,
     1     ncomp,x,t,p,mseg,densta,dense,dap_dx,dbp_dx,
     2     order1,order2,apar,bpar,z0t,z1t,z2t,z3t,dij_ab,
     3     PI,RGAS,NA,KBOL,TAU,sig_t,uij,sig_ij)
        y1 = (pges / p ) - 1.d0

        dydx = pgesdz/p
        x1   = x1 - y1/ dydx
        IF (x1.GT.0.9d0) x1 = 0.6d0
        IF (x1.LE.0.d0)  x1 = 1.d-10

        start = 1
        IF (DABS(y1).LT.acc_i) start = 0
        IF (i.GT.max_i) THEN 
          start = 0
          write (*,*) 'density iteration failed'
c          stop
        ENDIF

      IF (start.EQ.1) GOTO 1
      dense = x1

      RETURN
      END    


c**************************************************************************
c**************************************************************************
c**************************************************************************

      SUBROUTINE F_EOS (fres,gij,
     1     ncomp,x,t,p,mseg,densta,dense,dap_dx,dbp_dx,
     2     order1,order2,apar,bpar,z0t,z1t,z2t,z3t,dij_ab,
     3     PI,RGAS,NA,KBOL,TAU,sig_t,uij,sig_ij)

      IMPLICIT NONE

c-----variables used in the parameter list of subroutine--------------
      INTEGER nc
      PARAMETER (nc=20)
      INTEGER ncomp

      DOUBLE PRECISION fres,gij(nc,nc)
      DOUBLE PRECISION x(nc),t,p,mseg(nc)
      DOUBLE PRECISION densta,dense,dap_dx(nc,7),dbp_dx(nc,7)
      DOUBLE PRECISION order1,order2,apar(7),bpar(7)
      DOUBLE PRECISION z0t,z1t,z2t,z3t,z0,z1,z2,z3
      DOUBLE PRECISION PI, RGAS, NA, KBOL, TAU
      DOUBLE PRECISION dij_ab(nc,nc),uij(nc,nc),sig_ij(nc,nc),sig_t(nc)
c---------------------------------------------------------------------

c-----local variables-------------------------------------------------
      INTEGER i,m
      DOUBLE PRECISION zms,rho,m_mean
      DOUBLE PRECISION I1,I2,c1_con
      DOUBLE PRECISION fhs,fhc,fdsp
c---------------------------------------------------------------------

c-----abbreviations---------------------------------------------------
      rho = dense/z3t
      z0 = z0t*rho
      z1 = z1t*rho
      z2 = z2t*rho
      z3 = z3t*rho

      zms    = 1.d0 - dense
      m_mean = z0t/(PI/6.d0)

c-----Helmh. free energy : hard sphere contribution-------------------
      fhs= m_mean*(  3.d0*z1*z2/zms + z2**3.d0/z3/zms/zms
     &        + (z2**3.d0/z3/z3-z0)*DLOG(zms)  )/z0

c-----Helmh. free energy : chain term---------------------------------
      fhc = 0.d0
      DO i = 1,ncomp
        fhc = fhc + x(i) *(1.d0- mseg(i)) *DLOG(gij(i,i))
      END DO

c-----Helmh. free energy : dispersion contribution--------------------
      I1 = 0.d0
      I2 = 0.d0
      DO m=0,6
        I1 = I1 + apar(m+1)*dense**DBLE(m)
        I2 = I2 + bpar(m+1)*dense**DBLE(m)
      END DO
      c1_con= 1.d0/
     &    (  1.d0 + m_mean*(8.d0*dense-2.d0*dense**2.d0)/ZMS**4.d0
     &         + (1.d0 - m_mean)*(20.d0*dense-27.d0*dense**2.d0
     &                           +12.d0*dense**3.d0-2.d0*dense**4.d0)
     &                           /(ZMS*(2.d0-dense))**2.d0  )

      fdsp  = -2.d0*PI*rho*I1*order1 - PI*rho*c1_con*m_mean*I2*order2

c-----resid. Helmholtz free energy------------------------------------
      fres = fhs + fhc + fdsp
      
      RETURN
      END



c*********************************************************************
c*********************************************************************
c*********************************************************************

      SUBROUTINE PERTPAR (kij,parame,
     1     ncomp,x,t,p,mseg,densta,dense,dap_dx,dbp_dx,
     2     order1,order2,apar,bpar,z0t,z1t,z2t,z3t,dij_ab,
     3     PI,RGAS,NA,KBOL,TAU,sig_t,uij,sig_ij)

      IMPLICIT NONE

c-----variables used in the parameter list of subroutine--------------
      INTEGER nc
      PARAMETER (nc=20)
      INTEGER ncomp

      DOUBLE PRECISION kij(nc,nc),parame(nc,25)
      DOUBLE PRECISION x(nc),t,p,mseg(nc)
      DOUBLE PRECISION densta,dense,dap_dx(nc,7),dbp_dx(nc,7)
      DOUBLE PRECISION order1,order2,apar(7),bpar(7)
      DOUBLE PRECISION z0t,z1t,z2t,z3t,z0,z1,z2,z3
      DOUBLE PRECISION PI, RGAS, NA, KBOL, TAU
      DOUBLE PRECISION dij_ab(nc,nc),uij(nc,nc),sig_ij(nc,nc),sig_t(nc)
c---------------------------------------------------------------------

c-----local variables-------------------------------------------------
      INTEGER i,j,k,m
      DOUBLE PRECISION m_mean
      DOUBLE PRECISION ap(7,3),bp(7,3)
      DOUBLE PRECISION sig(nc),u0k(nc)
c---------------------------------------------------------------------

c-----constants-------------------------------------------------------
      PI   = 3.14159265359d0
      RGAS = 8.31441d0 
      NA   = 6.022045d23
      KBOL = RGAS/NA
      TAU  = PI/3.d0/DSQRT(2.d0)

c-----dispersion term constants---------------------------------------
      ap(1,1)=  0.91056314451539d0
      ap(1,2)= -0.30840169182720d0
      ap(1,3)= -0.09061483509767d0
      ap(2,1)=  0.63612814494991d0
      ap(2,2)=  0.18605311591713d0
      ap(2,3)=  0.45278428063920d0
      ap(3,1)=  2.68613478913903d0
      ap(3,2)= -2.50300472586548d0
      ap(3,3)=  0.59627007280101d0
      ap(4,1)= -26.5473624914884d0
      ap(4,2)=  21.4197936296668d0
      ap(4,3)= -1.72418291311787d0
      ap(5,1)=  97.7592087835073d0
      ap(5,2)= -65.2558853303492d0
      ap(5,3)= -4.13021125311661d0
      ap(6,1)= -159.591540865600d0
      ap(6,2)=  83.3186804808856d0
      ap(6,3)=  13.7766318697211d0
      ap(7,1)=  91.2977740839123d0
      ap(7,2)= -33.7469229297323d0
      ap(7,3)= -8.67284703679646d0

      bp(1,1)=  0.72409469413165d0
      bp(1,2)= -0.57554980753450d0
      bp(1,3)=  0.09768831158356d0
      bp(2,1)=  1.11913959304690d0  *2.d0
      bp(2,2)=  0.34975477607218d0  *2.d0
      bp(2,3)= -0.12787874908050d0  *2.d0
      bp(3,1)= -1.33419498282114d0  *3.d0
      bp(3,2)=  1.29752244631769d0  *3.d0
      bp(3,3)= -3.05195205099107d0  *3.d0
      bp(4,1)= -5.25089420371162d0  *4.d0
      bp(4,2)= -4.30386791194303d0  *4.d0
      bp(4,3)=  5.16051899359931d0  *4.d0
      bp(5,1)=  5.37112827253230d0  *5.d0
      bp(5,2)=  38.5344528930499d0  *5.d0
      bp(5,3)= -7.76088601041257d0  *5.d0
      bp(6,1)=  34.4252230677698d0  *6.d0
      bp(6,2)= -26.9710769414608d0  *6.d0
      bp(6,3)=  15.6044623461691d0  *6.d0
      bp(7,1)= -50.8003365888685d0  *7.d0
      bp(7,2)= -23.6010990650801d0  *7.d0
      bp(7,3)= -4.23812936930675d0  *7.d0


c-----pure component parameters---------------------------------------
      DO 1 i = 1,ncomp
         mseg(i) = parame(i,1)
         sig(i)  = parame(i,2)
         u0k(i)  = parame(i,3)
         sig_t(i)= sig(i)*(1.d0-0.12d0*DEXP(-3.d0*parame(i,3)/t))
 1    CONTINUE


c-----combination rules-----------------------------------------------
      DO 2 i = 1, ncomp
        DO 21 j = 1, ncomp
          sig_ij(i,j)=0.5d0*( sig(i) + sig(j) )
          uij(i,j)= (1.d0-kij(i,j))*(u0k(i)*u0k(j))**0.5d0
 21     CONTINUE
 2    CONTINUE


c-----abbreviations---------------------------------------------------
c-----abbreviations---------------------------------------------------
      z0t = 0.d0
      z1t = 0.d0
      z2t = 0.d0
      z3t = 0.d0
      DO i = 1,ncomp
        z0t = z0t + x(i) * mseg(i)
        z1t = z1t + x(i) * mseg(i) * sig_t(i)
        z2t = z2t + x(i) * mseg(i) * sig_t(i)**2.d0
        z3t = z3t + x(i) * mseg(i) * sig_t(i)**3.d0
      END DO

      m_mean = z0t
      z0t = PI/6.d0*z0t
      z1t = PI/6.d0*z1t
      z2t = PI/6.d0*z2t
      z3t = PI/6.d0*z3t

      DO i = 1,ncomp
        DO j=1,ncomp
          dij_ab(i,j)=sig_t(i)*sig_t(j)/(sig_t(i)+sig_t(j))
        ENDDO
      END DO


c-----dispersion term parameters for chain molecules------------------
      DO m=1,7
        apar(m) = ap(m,1) + (1.d0-1.d0/m_mean)*ap(m,2)
     &          + (1.d0-1.d0/m_mean)*(1.d0-2.d0/m_mean)*ap(m,3)
        bpar(m) = bp(m,1) + (1.d0-1.d0/m_mean)*bp(m,2)
     &           + (1.d0-1.d0/m_mean)*(1.d0-2.d0/m_mean)*bp(m,3)
      END DO

c-----derivatives of apar, bpar to mole fraction ( d(apar)/d(x) )-----
      DO k=1,ncomp
        DO m=1,7
          dap_dx(k,m) = mseg(k)/m_mean**2.d0*ap(m,2)
     &                 +(3.d0*mseg(k)/m_mean**2.d0
     &                            - 4.d0*mseg(k)/m_mean**3.d0)*ap(m,3)
          dbp_dx(k,m) = mseg(k)/m_mean**2.d0*bp(m,2)
     &                 +(3.d0*mseg(k)/m_mean**2.d0
     &                            - 4.d0*mseg(k)/m_mean**3.d0)*bp(m,3)
        END DO
      END DO

c-----van der Waals mixing rules for perturbation terms---------------
      order1 = 0.d0
      order2 = 0.d0
      DO i = 1,ncomp
        DO j = 1,ncomp
          order1 = order1 + x(i)*x(j)* mseg(i)*mseg(j)
     &                      *sig_ij(i,j)**3.d0 * uij(i,j)/t
          order2 = order2 + x(i)*x(j)* mseg(i)*mseg(j)
     &                      *sig_ij(i,j)**3.d0 * (uij(i,j)/t)**2.d0
        END DO
      END DO

      RETURN
      END




c*********************************************************************
c*********************************************************************
c*********************************************************************

      SUBROUTINE H_EOS (h_res,s_res,g_res,fres,pges,pgesdz,gij,parame,
     1     ncomp,x,t,p,mseg,densta,dense,dap_dx,dbp_dx,
     2     order1,order2,apar,bpar,z0t,z1t,z2t,z3t,dij_ab,
     3     PI,RGAS,NA,KBOL,TAU,sig_t,uij,sig_ij)

      IMPLICIT NONE
c----------------------formal parameters-----------------------------------
      INTEGER nc
      PARAMETER (nc=20)
      INTEGER ncomp

      DOUBLE PRECISION parame(nc,25)
      DOUBLE PRECISION h_res,s_res,g_res
      DOUBLE PRECISION pges,pgesdz,gij(nc, nc)
      DOUBLE PRECISION fres
      DOUBLE PRECISION x(nc),t,p,mseg(nc)
      DOUBLE PRECISION densta,dense,dap_dx(nc,7),dbp_dx(nc,7)
      DOUBLE PRECISION order1,order2,apar(7),bpar(7)
      DOUBLE PRECISION z0t,z1t,z2t,z3t,z0,z1,z2,z3
      DOUBLE PRECISION PI, RGAS, NA, KBOL, TAU
      DOUBLE PRECISION dij_ab(nc,nc),uij(nc,nc),sig_ij(nc,nc),sig_t(nc)
c---------------------------------------------------------------------

c-----local variables-------------------------------------------------
      DOUBLE PRECISION zres,df_dT,dfdr,ddfdrdr
c     &                 ,cv_res,df_dTdT, df_drdt

      INTEGER i,j,m
      DOUBLE PRECISION sigtdt(nc),rho,dijabdt(nc,nc),
     &                 z1dt,z2dt,z3dt,zms,fhsdt,fhcdt,fdspdt,
     &                 m_mean,I1,I2,I1dt,
     &                 I2dt,c1_con,c2_con,c1_dt,gijdt(nc,nc)

c-----------------------------------------------------------------------

      rho = dense/z3t
      z0 = z0t*rho
      z1 = z1t*rho
      z2 = z2t*rho
      z3 = z3t*rho

      zms    = 1.d0 - dense
      m_mean = z0t/(PI/6.d0)

      zres = (pges * 1.d-30)/(KBOL*T*rho) - 1.d0

c-----first and second derivative of fres/kT to density (dfdr,ddfdrdr)--
      dfdr    = pges  /(dense*rho*(KBOL*T)/1.d-30)
      ddfdrdr = pgesdz/(dense*rho*(KBOL*T)/1.d-30)
     &          - dfdr*2.d0/dense - 1.d0/dense**2.d0


c------derivatives to temperature---------------------------------------
      DO i = 1,ncomp
        sigtdt(i)=parame(i,2)
     &   *(-3.d0*parame(i,3)/T/T)*0.12d0*DEXP(-3.d0*parame(i,3)/T)
      END DO

      DO i = 1,ncomp
        DO j = 1,ncomp
          dijabdt(i,j)=dij_ab(i,j)*(sigtdt(i)/sig_t(i)
     &   +sigtdt(j)/sig_t(j)-(sigtdt(i)+sigtdt(j))/(sig_t(i)+sig_t(j)))
        END DO
      END DO

      z1dt = 0.d0
      z2dt = 0.d0
      z3dt = 0.d0
      DO i = 1,ncomp
        z1dt = z1dt + x(i) * mseg(i) * sigtdt(i)
        z2dt = z2dt + x(i) * mseg(i) * 2.d0*sig_t(i)*sigtdt(i)
        z3dt = z3dt + x(i) * mseg(i) * 3.d0*sig_t(i)*sig_t(i)*sigtdt(i)
      END DO
      z1dt = PI/6.d0*z1dt*rho
      z2dt = PI/6.d0*z2dt*rho
      z3dt = PI/6.d0*z3dt*rho

      DO i = 1, ncomp
        DO j = 1, ncomp
          gijdt(i,j)= z3dt/zms/zms
     &               + 3.d0*(dijabdt(i,j)*z2+dij_ab(i,j)*z2dt)/zms/zms
     &               + 4.d0*dij_ab(i,j)*z2*(1.5d0*z3dt+dijabdt(i,j)*z2
     &                                    +dij_ab(i,j)*z2dt)/zms**3.d0
     &               + 6.d0*(dij_ab(i,j)*z2)**2.d0*z3dt/zms**4.d0
        ENDDO
      ENDDO

c------derivative of fres/kT to T: hard-sphere term (fhsdt)-------------
      fhsdt = 6.D0/PI/rho*(  3.d0*(z1dt*z2+z1*z2dt)/zms
     &     + 3.d0*z1*z2*z3dt/zms/zms
     &     + 3.d0*z2*z2*z2dt/z3/zms/zms
     &     + z2**3.d0*(2.d0*z3*z3dt-z3dt*zms)/(z3*z3*zms**3.d0)
     &     + (3.d0*z2*z2*z2dt*z3-2.d0*z2**3.d0*z3dt)/z3**3.d0*DLOG(zms)
     &     + (PI/6.d0*rho*m_mean -z2**3.d0/z3/z3)*z3dt/zms  )

c------derivative of fres/kT to T: hard-chain term (fhcdt)--------------
      fhcdt  = 0.d0
      DO i = 1, ncomp
        fhcdt = fhcdt + x(i) * (1.d0-mseg(i)) * gijdt(i,i) / gij(i,i)
      ENDDO

c------Ableitung von f/kT fuer den Dispersions Term (fdspdt)------------
      I1 = 0.d0
      I2 = 0.d0
      I1dt = 0.d0
      I2dt = 0.d0
      DO m=0,6
        I1 = I1 + apar(m+1)*z3**DBLE(m)
        I2 = I2 + bpar(m+1)*z3**DBLE(m)
        I1dt = I1dt + apar(m+1)*z3dt*DBLE(m)*z3**DBLE(m-1)
        I2dt = I2dt + bpar(m+1)*z3dt*DBLE(m)*z3**DBLE(m-1)
      END DO

      c1_con= 1.d0/
     &    (  1.d0 + m_mean*(8.d0*z3-2.d0*z3**2.d0)/ZMS**4.d0
     &         + (1.d0 - m_mean)*(20.d0*z3-27.d0*z3**2.d0
     &                           +12.d0*z3**3.d0-2.d0*z3**4.d0)
     &                           /(ZMS*(2.d0-z3))**2.d0  )
      c2_con= - c1_con*c1_con
     &          *(m_mean*(-4.d0*z3**2.d0+20.d0*z3+8.d0)/ZMS**5.d0
     &        + (1.d0 - m_mean)
     &        *(2.d0*z3**3.d0+12.d0*z3**2.d0-48.d0*z3+40.d0)
     &                                /(ZMS*(2.d0-z3))**3.d0 )
      c1_dt= c2_con*z3dt
c      c3_con= 2.d0 * c2_con*c2_con/c1_con
c     &      - c1_con*c1_con
c     &       *( m_mean*(-12.d0*z3**2.d0+72.d0*z3+60.d0)/ZMS**6.d0
c     &       + (1.d0 - m_mean)
c     &        *(-6.d0*z3**4.d0-48.d0*z3**3.d0+288.d0*z3**2.d0
c     &                                            -480.d0*z3+264.d0)
c     &                                /(ZMS*(2.d0-z3))**4.d0 )

      fdspdt  = -2.d0*PI*rho*(I1dt-I1/T)*order1
     & -PI*rho*m_mean*(c1_dt*I2+c1_con*I2dt-2.d0*c1_con*I2/T)*order2

c-------gesamte Ableitung von fres/kT nach der Temperatur---------------

      df_dt= fhsdt + fhcdt + fdspdt



      s_res  =   ( - df_dt - fres/t )*RGAS*t  +  RGAS *DLOG(zres+1.d0)
      h_res  =   ( - t*df_dt  +  zres  ) * RGAS*t
      g_res  = h_res - T * s_res
c      g_res  =   (fres-DLOG(zres+1.d0)+zres)* RGAS*t

      RETURN
      END
