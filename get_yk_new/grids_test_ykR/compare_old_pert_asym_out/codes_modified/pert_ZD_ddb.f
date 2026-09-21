cpn Sep 2007 Add a separate calculation for W- production, including the
c       option  to compute s+g -> W- + c for ido_cbar=1
CCPY Jan. 9, 2007: Add the option to return (c-bar + Z0 ) rate, which
C                   is coded through the flag "IDO_CBAR"  
C Dec. 20, 2006: Add the option to return (c-bar + W+ ) rate, which
C                   is coded through the flag "IDO_CBAR"  
C **************************************************************************
C
C  MODULE: WXINT                                  Version 3.1
C  -------------
C                   by Peter Arnold
C                   Ref:  P. Arnold and H. Reno, Nucl. Phys. B319 (1989) 37
C
C  Changes from Version 1 to 2: routine XPART
C  Changes from Version 2 to 3: The calls to alpha(s) have been changed
C    to work with the 4-flavor alpha(s) routine I took from Papageno.
C    The routine A2 has been changed so that Feynman x is always <= 1.
C    The X and Q**2 range warnings in INTSET have been changed to work
C    correctly.  Most importantly, GETSTR has been rewritten to call
C    the structure function routines directly each time it is called.
C    Previously, for the sake of time efficiency, it would save all the
C    structure functions in a table and do its own linear interpolation
C    when called.  The matching of the perturbative to resummed calculations
C    requires delicate cancellations, however, and it is important that
C    each use *exactly* the same approximations to the strutcure functions.
C  Changes Version 3 to 3.1: Set up to allow p(bar)+p(bar) collisions.
C
C
C      THIS MODULE CALCULATES D(SIGMA)/DY/D(QT**2) FOR GIVEN Y AND QT.
C  THIS IS RETURNED BY THE TOP-LEVEL ROUTINES DYPT1 AND DYPT2 FOR FIRST
C  AND SECOND ORDER CONTRIBUTIONS RESPECTIVELY.  RESULTS ARE RETURNED
C  IN PB/GEV**2.  ALL OVERALL NORMALIZATIONS, LEFT OUT OF LOWER MODULES,
C  HAVE BEEN INCLUDED.  ALL OTHER QUANTITIES ARE EXPRESSED IN POWERS OF
C  GEV.
C      AT THE MOMENT, THE FACTORIZATION SCALE M IS TAKEN EQUAL TO Q.
C
C      RATHER THAN D(X1)D(X2), THE PARTON MOMENTUM-FRACTION INTEGRALS
C  ARE SPLIT INTO D(X1)D(S2) AND D(X2)D(S2) AS IN ELLIS, ET. AL. AND
C  KAJANTIE & LINDFORS.  TO SMOOTH OUT THE INTEGRALS, WE THEN TRANSFORM X,S2
C  TO ZX,ZS2 WHICH RUN FROM 0 TO 1 AND ARE DESCRIBED IN ROUTINES
C  ZTOX AND ZTOS2.
C
C  SUMMARY OF ROUTINES
C  -------------------
C
C      IN ADDITION TO WDYPT1 AND WDYPT2, ANOTHER USEFUL ROUTINE IS
C  YMAX WHICH RETURNS THE MAXIMUM KINEMATICALLY ALLOWED VALUE OF THE
C  RAPIDITY.
C      THE FOLLOWING GLOBAL VARIABLES MUST BE SET BEFORE USING DYPT1
C  AND DYPT2: JWTYPE,KCOLLD,QQ.
C
C  OTHER MODULES REFERENCED: COMMON, CODES, WPART, 2GAUSS.
C
C  REFERENCES
C  ----------
C      R. ELLIS, ET. AL., NUCL. PHYS. B211 (1983) 106.
C      K. KAJANTIE AND J. LINDFORS, NUCL. PHYS. B146 (1978) 465.
C
C **************************************************************************

          
C **  JWTYPE -- INDICATES THE TYPE OF PHOTON PRODUCED BY THE FOLLOWING
C **  CODES         
C **      
C **     JWP,JWM       W+, W- BOSON         
C **     JZ            Z BOSON
C **     JPH           PHOTON (NOT IMPLEMENTED)   
          

C **  A1 -- THE UPPER BOUND FOR S2 IN THE D(X1)D(S2) INTEGRATION.
C **  THE INTEGRATION REGION IS EVERYTHING KINEMATICALLY ALLOWED
C **  FOR X1 > XLIM(1).  [XLIM(1) IS SET UP BY ROUTINE INTSET.]
C **  THE MAXIMUM VALUE OF X2 IS 1, WHICH IS USED HERE TO DETERMINE
C **  THE UPPER BOUND ON S2 AFTER CHANGING VARIABLES FROM D(X2) TO
C **  D(S2).

      FUNCTION A1(X1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /STU/ SS,TT,UU,QQ,QM,QMM,qmu
      A1 = UU + X1*(SS+TT-QQ)
      RETURN
      END


C **  A2 -- THE UPPER BOUND FOR S2 IN THE D(X2)D(S2) INTEGRATION.
C **  THE INTEGRATION REGION IS EVERYTHING KINEMATICALLY ALLOWED
C **  FOR X1 < XLIM(1).  XLIM(1) IS SET UP BY ROUTINE INTSET AND
C **  IS USED HERE TO DETERMINE THE UPPER BOUND
C **  ON S2 AFTER CHANGING VARIABLES FROM D(X1) TO D(S2).  THERE IS
C **  NO NEED TO INTEGRATE X1>1.

      FUNCTION A2(X2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /STU/ SS,TT,UU,QQ,QM,QMM,qmu
      COMMON /XLIM/ XLIM(2),QT,Y
      X1MAX = XLIM(1)
      IF (X1MAX.GT.1.D0) X1MAX = 1.D0
      A2 = X1MAX*(X2*SS+TT-QQ) + QQ + X2*(UU-QQ)
      RETURN
      END


C **  XTOZ, ZTOX, XJACB -- THESE FUNCTIONS IMPLEMENT A CHANGE OF VARIABLES
C **  THAT SMOOTHS OUT THE X INTEGRATION.  THE INTEGRAND DIVERGES
C **  LOGARITHMICALLY WITH A CUT-OFF ORDER SQRT(QT**2/SS) FROM THE
C **  LOWER LIMIT OF THE X INTEGRATION.  WE TRANSFORM TO AN APPROPROATE
C **  VARIABLE ZX THAT RUNS FROM 0 TO 1 AND SMOOTHS OUT THE INTEGRAND.
C **
C **      XTOZ       CONVERTS Z TO X
C **      ZTOX       CONVERTS X TO Z
C **      XJACB      JACOBIAN OF TRANSFORM: D(X) = XJACB*D(XZ)
C **
C **  I = 1 OR 2 SPECIFIES X1 OR X2.  XFUDGE IS CHOSEN BY TRIAL AND
C **  ERROR TO MAKE INTEGRAL CONVERGE AS WELL AS POSSIBLE.

      FUNCTION XTOZ(I,X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        PARAMETER (XFUDGE=0.1)
      COMMON /STU/ SS,TT,UU,QQ,QM,QMM,qmu
      COMMON /XLIM/ XLIM(2),QT,Y
      XPEAK = XFUDGE * SQRT(QT**2/SS)
      XTOZ = LOG((X-XLIM(I)+XPEAK)/XPEAK) / LOG((1-XLIM(I)+XPEAK)/XPEAK)
      RETURN

      ENTRY ZTOX(I,Z)
      XPEAK = XFUDGE * SQRT(QT**2/SS)
      ZTOX = XPEAK * ((1-XLIM(I)+XPEAK)/XPEAK)**Z + XLIM(I) - XPEAK
      RETURN

      ENTRY XJACB(I,X)
      XPEAK = XFUDGE * SQRT(QT**2/SS)
      XJACB = (X-XLIM(I)+XPEAK)*LOG((1-XLIM(I)+XPEAK)/XPEAK)
      RETURN
      END


C **  ZTOS2, S2JACB  --  SIMILAR TO ZTOX ETC. ABOVE BUT TRANSFORMS THE
C **  S2 INTEGRAL.  THE LOGARITHMIC DIVERGENCE IS CUT-OFF AROUND QT**2.
C **  Z IS DEFINED BY
C **
C **      Z = LOG((X+S2PEAK)/S2PEAK) / LOG((A+S2PEAK)/S2PEAK)

      FUNCTION ZTOS2(Z,A)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        PARAMETER (S2FUDG=0.1)
      COMMON /XLIM/ XLIM(2),QT,Y
      S2PEAK = S2FUDG * QT**2
      ZTOS2 = S2PEAK * ((A+S2PEAK)/S2PEAK)**Z - S2PEAK
      RETURN

      ENTRY S2JACB(S2,A)
      S2PEAK = S2FUDG * QT**2
      S2JACB = (S2+S2PEAK)*LOG((A+S2PEAK)/S2PEAK)
      RETURN
      END


C **  GETSTR -- FOR GIVEN MOMENTUM FRACTIONS X1 AND X2, GETSTR GETS THE
C **  STRUCTURE FUNCTIONS FOR THE COLLIDING NUCLEONS AND STORES THEM IN
C **  COMMON.  PARTICLE 1 IS TAKEN TO BE A PROTON.  PARTICLE 2 IS EITHER A
C **  PROTON OR ANTI-PROTON DEPENDING ON KCOLLD.
C **      THE DISTRIBUTION FUNCTIONS ARE ASSUMED TO BE CALLABLE IN THE
C **  PESKIN FORMAT AND THEREFORE RETURN X*F(X) RATHER THAN F(X).
C **  CARE MUST BE TAKEN TO LATER DIVIDE THE INTEGRAND BY X1*X2 TO
C **  COMPENSATE.
C **      AFTER DISCOVERING THAT THE STRUCTURE FUNCTION EVALUATIONS WERE
C **  TAKING BETWEEN 50 TO 75% OF THE COMPUTATION TIME FOR THIS PROGRAM,
C **  I REORGANIZED TO CALCULATE THEM ONLY ONCE.  THE DIEMOZ ET. AL.
C **  ROUTINES CALCULATE LOG(X) AND LOG(Q**2) EACH TIME THEY'RE CALLED AND
C **  DO A 2-DIMENSIONAL LINEAR INTERPOLATION IN X AND Q**2.  NOW, I TABULATE
C **  THE FUNCTIONS FOR A GIVEN Q**2, TAKE LOG(X) ONLY ONCE FOR ALL THE
C **  STRUCTURE FUNCTIONS, AND LINEAR INTERPOLATE IN X.
C **      WARNING: IF I EVER USE EICHTEN'S ROUTINES, REMEMBER THAT THE
C **  ARGUMENT IS QM**2 RATHER THAN QM.

      SUBROUTINE GETSTR(X1,X2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C        SAVE QMLAST
C        DATA QMLAST /-1/
C      PARAMETER (IZDIM=1000)
C      COMMON /STARRY/ ZFU(IZDIM),ZFD(IZDIM),ZFSEA(IZDIM),ZFS(IZDIM)
C     1  ,ZFC(IZDIM),ZFB(IZDIM),ZFG(IZDIM)
      COMMON /STU/ SS,TT,UU,QQ,QM,QMM,qmu
      COMMON /STRUCT/
     1         U1,UB1,D1,DB1,ST1,STB1,C1,CB1,B1,BB1,G1
     2        ,U2,UB2,D2,DB2,ST2,STB2,C2,CB2,B2,BB2,G2

      REAL*8 ECM
      INTEGER IBEAM, ihadron1, ihadron2
      common/collider/ecm,ibeam
cpn2007
      COMMON /JCODES/ JWTYPE,  JWP,JWM,JZ,JPH

      external pdfh

      if (abs(ibeam).gt.1) then
        write (*,*) 'GETSTR: ibeam = ',ibeam,' not implemented'
        call error
      endif                     !abs(ibeam)



cpn2007 I use the following hotfix for the W- case. The cross section
cpn2007 for W- production in a pp collision must be equal to the cross
cpn2007 section for W+ production in a pbar-pbar collision. Hence, the W-
cpn2007 cross section for pp, p-pbar is obtained by exchanging proton and 
cpn2007 antiproton PDF's in the W+ cross section without modifying the
cpn2007 rest of the code
      if (jwtype.eq.jwm) then !W- 
        ihadron1=-1; ihadron2=-ibeam
      else                    !W+, Z, etc. 
        ihadron1=1; ihadron2=ibeam
      endif                     !jwtype

CZL
      u1   = x1*pdfh(ihadron1, 1,x1,qm)
      ub1  = x1*pdfh(ihadron1,-1,x1,qm)
      d1   = x1*pdfh(ihadron1, 2,x1,qm)
      db1  = x1*pdfh(ihadron1,-2,x1,qm)
      st1  = x1*pdfh(ihadron1, 3,x1,qm)
      stb1 = x1*pdfh(ihadron1,-3,x1,qm)
      c1   = x1*pdfh(ihadron1, 4,x1,qm)
      cb1  = x1*pdfh(ihadron1,-4,x1,qm)
      b1   = x1*pdfh(ihadron1, 5,x1,qm)
      bb1  = x1*pdfh(ihadron1,-5,x1,qm)
      g1   = x1*pdfh(ihadron1, 0,x1,qm)

      u2   = x2*pdfh(ihadron2, 1,x2,qm)
      ub2  = x2*pdfh(ihadron2,-1,x2,qm)
      d2   = x2*pdfh(ihadron2, 2,x2,qm)
      db2  = x2*pdfh(ihadron2,-2,x2,qm)
      st2  = x2*pdfh(ihadron2, 3,x2,qm)
      stb2 = x2*pdfh(ihadron2,-3,x2,qm)
      c2   = x2*pdfh(ihadron2, 4,x2,qm)
      cb2  = x2*pdfh(ihadron2,-4,x2,qm)
      b2   = x2*pdfh(ihadron2, 5,x2,qm)
      bb2  = x2*pdfh(ihadron2,-5,x2,qm)
      g2   = x2*pdfh(ihadron2, 0,x2,qm)

      RETURN
      END


C **  STINIT -- TABULATES THE STRUCTURE FUNCTIONS FOR LATER INTERPOLATION.
C **  WE TAKE IZDIM EVENLY SPACED POINTS IN Z = LOG(X) FROM Z=ZMIN TO 0.

      SUBROUTINE STINIT(QM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (IZDIM=1000)
      PARAMETER (ZMIN=-9.9035)
C                                      ! LOG(5.E-5)
      COMMON /STARRY/ ZFU(IZDIM),ZFD(IZDIM),ZFSEA(IZDIM),ZFS(IZDIM)
     1  ,ZFC(IZDIM),ZFB(IZDIM),ZFG(IZDIM)
      DO 10, I=1,IZDIM
        Z = ZMIN + (I-1)*(-ZMIN)/(IZDIM-1)
        X = EXP(Z)
        ZFU(I) = XFU(X,QM)
        ZFD(I) = XFD(X,QM)
        ZFSEA(I) = XFSEA(X,QM)
        ZFS(I) = XFS(X,QM)
        ZFC(I) = XFC(X,QM)
        ZFB(I) = XFB(X,QM)
        ZFG(I) = XFG(X,QM)
10    CONTINUE
      RETURN
      END


C **  LOOKUP -- AIDS IN LINEAR INTERPOLATION OF STRUCTURE FUNCTIONS.
C **  FOR A GIVEN X, RETURNS THE INDEX I OF THE DATA POINT JUST BELOW
C **  X AND THE WEIGHTS DELA AND DELB OF THAT POINT AND THE NEXT.  THAT
C **  IS, IF ZF CONTAINS THE TABULATED VALUES INITIALIZED BY STINIT,
C **  THEN THE STRUCTURE FUNCTION AT X IS
C **
C **               DELA*ZF(I) + DELB*ZF(I+1)

      SUBROUTINE LOOKUP(X,I,DELA,DELB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (IZDIM=1000)
      PARAMETER (ZMIN=-9.9035)
      Z = LOG(X)
      IF ((Z.GT.ZMIN) .AND. (Z.LT.0)) THEN
        ZINT = (Z-ZMIN)/(-ZMIN)*(IZDIM-1) + 1
        I = INT(ZINT)
        DELB = ZINT-I
        DELA = 1-DELB
      ELSE
        I = 1
        DELA = 0
        DELB = 0
      END IF
      RETURN
      END


C **  FSIG -- THIS IS A COMMON ROUTINE FOR ALL THE INTEGRANDS.  FOR A
C **  GIVEN X1,X2,S2 IT WILL RETURN WHAT WOULD BE THE INTEGRAND OF
C **  A D(X1)D(X2) INTEGRAL (WHICH IS THEN MULTIPLIED BY THE APPROPRIATE
C **  JACOBIAN BY THE CALLER).  IT USES XPART TO CALL THE APPROPRIATE
C **  ROUTINE OF MODULE WPART.
C **      THE ARGUMENTS ALSO INCLUDE THE UPPER LIMIT A OF THE S2
C **  INTEGRATION AND S2A WHICH IS ALWAYS THE S2 INTEGRATION VARIABLE.
C **  THE ARGUMENT S2 IS USUALLY THE SAME EXCEPT WHEN EVALUATING THE
C **  S2=0 PART OF [1/S2] AND [LOG(S2)/S2] TERMS WHEN IT MUST BE ZERO.
C **  LTERM2 INDICATES WHICH TERM IS TO BE EVALUATED
C **
C **       LTERM2 = LDS2     DELTA(S2) TERMS
C **       LTERM2 = LS20     COEFFICIENT OF [1/S2] AT S2=0
C **       LTERM2 = LWS20    COEFFICIENT OF [LOG(S2)/S2] AT S2=0
C **       LTERM2 = L2DIM    EVERYTHING ELSE
C **
C **  THE OVERALL FACTOR 1/4/NC/S IS INCORPORATED HERE AS WELL AS THE
C **  1/X1/X2 TO CORRECT FOR THE STRUCTURE FUNCTIONS (SEE GETSTRUCT).

      FUNCTION FSIG(LTERM2,X1,X2,S2,A)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /DEBUG/ DEBUG,DEBG0
        LOGICAL DEBUG,DEBG0
      COMMON /CLEBSH/ CF,CNP,TR,NC,NU,ND,NUD     
      COMMON /STU/ SS,TT,UU,QQ,QM,QMM,qmu
      COMMON /LCODES/ LTERM, L2DIM,LS20,LWS20,LDS2
      SAVE NCALLS
      DATA NCALLS /0/
      LTERM = LTERM2
      CALL GETSTR(X1,X2)
C     ! COMPUTE PARTON VARIABLES.
      S = X1*X2*SS
      T = QQ - X1*(QQ-TT)
      U = QQ - X2*(QQ-UU)

      FSIG = (1./4./NC)/S * XPART(S,T,U,S2,A)/X1/X2
C     ! FACTOR OF 1/X1/X2 'CAUSE DISTRIBUTION FUNCTIONS ARE XF(X).
      IF (DEBUG.OR.DEBG0) THEN
        NCALLS = NCALLS + 1
CsB        IF (MOD(NCALLS,10000).EQ.0)
C     1    WRITE (*,*) 'NCALLS =',NCALLS,'X1,X2,FSIG =',X1,X2,FSIG
      END IF
            
      RETURN
      END


C **  FSIG1A, FSIG2A -- PROVIDE THE INTEGRANDS FOR D(ZX1)D(ZS2) AND
C **  D(ZX2)D(ZS2) INTEGRALS RESPECTIVELY.  FSIG1B AND FSIG2B ARE CALLED
C **  FOR THE INNER D(ZS2) INTEGRALS.  THE FACTORS IN FSIG1 COME FROM
C **
C **     XJACB  :    JACOBIAN OF X1-->ZX1
C **     FACTOR :    JACOBIAN OF X2-->S2
C **     S2JACB :    JACOBIAN OF S2-->ZS2
C **
C **  FOR EACH X, THE S2=0 PIECES OF [1/S2] AND [LOG(S2)/S2] ARE SAVED
C **  IN COMMON FOR THE APPROPRIATE SUBTRACTIONS IN THE INTEGRAND.


      FUNCTION FSIG1A(ZX1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        EXTERNAL FSIG1B
      COMMON /ERR/ ERREL,NINTVL
      COMMON /STU/ SS,TT,UU,QQ,QM,QMM,qmu
      COMMON /LCODES/ LTERM, L2DIM,LS20,LWS20,LDS2
      COMMON /FSIG1X/ X1,AA1,FS20,FWS20
      X1 = ZTOX(1,ZX1)
      AA1 = A1(X1)
      FACTOR = 1/(X1*SS+UU-QQ)
      X2 = FACTOR * (X1*(QQ-TT)-QQ)
      FS20 = FSIG(LS20,X1,X2,0.D0,AA1)
      FWS20 = FSIG(LWS20,X1,X2,0.D0,AA1)
      FSIG1A = XJACB(1,X1) * FACTOR
     1  * BGAUSS(FSIG1B,0.D0,1.D0,NINTVL,ERREL)
      RETURN
      END

      FUNCTION FSIG1B(ZS2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /STU/ SS,TT,UU,QQ,QM,QMM,qmu
      COMMON /LCODES/ LTERM, L2DIM,LS20,LWS20,LDS2
      COMMON /FSIG1X/ X1,AA1,FS20,FWS20
      S2 = ZTOS2(ZS2,AA1)
      FACTOR = 1/(X1*SS+UU-QQ)
      X2 = FACTOR * (S2+X1*(QQ-TT)-QQ)
      FSIG1B = S2JACB(S2,AA1)
     1  * ( FSIG(L2DIM,X1,X2,S2,AA1) - (FS20 + FWS20*LOG(S2/QMM))/S2 )
      RETURN
      END


      FUNCTION FSIG2A(ZX2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        EXTERNAL FSIG2B
      COMMON /ERR/ ERREL,NINTVL
      COMMON /STU/ SS,TT,UU,QQ,QM,QMM,qmu
      COMMON /LCODES/ LTERM, L2DIM,LS20,LWS20,LDS2
      COMMON /FSIG2X/ X2,AA2,FS20,FWS20
      X2 = ZTOX(2,ZX2)
      AA2 = A2(X2)
      FACTOR = 1/(X2*SS+TT-QQ)
      X1 = FACTOR * (X2*(QQ-UU)-QQ)
      FS20 = FSIG(LS20,X1,X2,0.D0,AA2)
      FWS20 = FSIG(LWS20,X1,X2,0.D0,AA2)
      FSIG2A = XJACB(2,X2) * FACTOR
     1  * BGAUSS(FSIG2B,0.D0,1.D0,NINTVL,ERREL)
      RETURN
      END

      FUNCTION FSIG2B(ZS2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /STU/ SS,TT,UU,QQ,QM,QMM,qmu
      COMMON /LCODES/ LTERM, L2DIM,LS20,LWS20,LDS2
      COMMON /FSIG2X/ X2,AA2,FS20,FWS20
      S2 = ZTOS2(ZS2,AA2)
      FACTOR = 1/(X2*SS+TT-QQ)
      X1 = FACTOR * (S2+X2*(QQ-UU)-QQ)
      FSIG2B = S2JACB(S2,AA2)
     1  * ( FSIG(L2DIM,X1,X2,S2,AA1) - (FS20 + FWS20*LOG(S2/QMM))/S2 )
      RETURN
      END


C ** FDEL1, FDEL2 -- ARE SIMILAR TO FSIG1 AND FSIG2 BUT ARE FOR DELTA(S2)
C ** INTEGRALS.


      FUNCTION FDEL1(ZX1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /STU/ SS,TT,UU,QQ,QM,QMM,qmu
      COMMON /LCODES/ LTERM, L2DIM,LS20,LWS20,LDS2
      X1 = ZTOX(1,ZX1)
      FACTOR = 1/(X1*SS+UU-QQ)
      X2 = FACTOR * (X1*(QQ-TT)-QQ)
      FDEL1 = XJACB(1,X1) * FACTOR * FSIG(LDS2,X1,X2,0.D0,A1(X1))
      
      RETURN
      END


      FUNCTION FDEL2(ZX2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /STU/ SS,TT,UU,QQ,QM,QMM,qmu
      COMMON /LCODES/ LTERM, L2DIM,LS20,LWS20,LDS2
      X2 = ZTOX(2,ZX2)
      FACTOR = 1/(X2*SS+TT-QQ)
      X1 = FACTOR * (X2*(QQ-UU)-QQ)
      FDEL2 = XJACB(2,X2) * FACTOR * FSIG(LDS2,X1,X2,0.D0,A2(X2))
      RETURN
      END


C **  YMAX -- RETURNS THE KINEMATIC LIMIT ON THE RAPIDITY FOR GIVEN
C **  HADRONIC S, Q**2, AND TRANSVERSE MOMENTUM QT.

      FUNCTION YMAX(SS,QQ,QT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      ACOSH(X) = LOG(X+SQRT(X**2-1))

      YMAX = ACOSH((SS+QQ)/2./SQRT(SS*(QQ+QT**2)))
      RETURN
      END


C **  INTSET -- A ROUTINE THAT SETS UP PARAMETERS COMMON TO THE TOP-LEVEL
C **  ROUTINE DYPT1 AND DYPT2.  INTSET SETS UP THE NUCLEON MANDELSTAHM
C **  VARIABLES SS,TT,UU; THE FACTORIZATION SCALE; AND THE LOWER LIMITS
C **  ZX1MIN AND ZX2MIN OF THE ZX1 AND ZX2 INTEGRATIONS.  THESE LIMITS WILL
C **  BE 0 UNLESS THE STANDARD SPLICING POINT (X1MIN,X2MIN) OF THE INTEGRALS
C **  DOES NOT LIE IN ([0,1],[0,1]).  THEN AT MOST ONE OF THE INTEGRALS, D(X1)
C **  OR D(X2), IS RELEVANT AND THE VALUE OF THE CORRESPONDING ZXMIN IS
C **  CHANGED TO RUN ONLY OVER THAT PART OF THE INTEGRATION REGION WITH X
C **  INSIDE ([0,1],[0,1]).  IF EITHER THE D(X1) OR D(X2) INTERGALS ARE NOT
C **  TO BE DONE, DOX1 OR DOX2 WILL BE SET TO FALSE ON RETURN.
C **      THE INITIAL EXPRESSION FOR XMIN IS SAVED IN COOMON AS XLIM
C **  FOR USE BY THE ROUTINES A2 AND ZTOX AND ITS RELATIVES.

      SUBROUTINE INTSET(QT1,Y1, ZX1MIN,ZX2MIN,DOX1,DOX2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        LOGICAL DOX1,DOX2
      COMMON /STU/ SS,TT,UU,QQ,QM,QMM,qmu
      COMMON /XLIM/ XLIM(2),QT,Y
      DATA XCUT,QQMAX /1E-5,1.31E6/

C     ! STORE QT AND Y INTO COMMON
      QT = QT1
      Y = Y1
C     ! MANDELSTAHM VARIABLES
      TT = QQ - EXP(-Y) * SQRT(SS*(QT**2+QQ))
      UU = QQ - EXP( Y) * SQRT(SS*(QT**2+QQ))
C     ! FACTORIZATION SCALE -- SAVE SQUARE IN COMMON FOR LATER USE
      QMM = QM**2

      IF (ABS(Y).GT.YMAX(SS,QQ,QT)) THEN
C       ! DON'T DO ANY INTEGRALS: RESULT = 0.
        DOX1 = .FALSE.
        DOX2 = .FALSE.
        RETURN
      END IF

C     ! GET NAIVE X LIMITS
      RTAUP = SQRT((QQ+QT**2)/SS) + SQRT(QT**2/SS)
      XLIM(1) = RTAUP * EXP(Y)
      XLIM(2) = RTAUP * EXP(-Y)
      X1MIN = XLIM(1)
      X2MIN = XLIM(2)
      ZX1MIN = 0
      ZX2MIN = 0
      DOX1 = .TRUE.
      DOX2 = .TRUE.
      
      IF ((XLIM(1).GT.1).AND.(XLIM(2).GT.1)) CALL ERROR
C     ! MODIFY X LIMITS IF OUTSIDE [0,1]
      IF (XLIM(1).GT.1) THEN
C       ! NO X1 INTEGRAL TO DO.
        X2MIN = - TT/(SS+UU-QQ)
        DOX1 = .FALSE.
        ZX2MIN = XTOZ(2,X2MIN)
      ELSE IF (XLIM(2).GT.1) THEN
C       ! NO X2 INTEGRAL TO DO.
        X1MIN = - UU/(SS+TT-QQ)
        ZX1MIN = XTOZ(1,X1MIN)
        DOX2 = .FALSE.
      END IF

      IF ((X1MIN.LT.XCUT) .OR. (X2MIN.LT.XCUT))
     1  WRITE(*,*) 'TROUBLE FOR HMRS ROUTINES: X = ',X1MIN,X2MIN
      IF (QQ.GT.QQMAX)
     1  WRITE(*,*) 'TROUBLE FOR HMRS ROUTINES: Q**2 = ',QQ
      
      RETURN
      END


C **  DYPT1 --  CALCULATES THE FIRST-ORDER CONTRIBUTION TO
C **  D(SIGMA)/D(Y)/D(QT**2) FOR GIVEN SQRT(S), QT, AND Y.
C **  THE RESULT IS RETURNED AS SIG1.  THE PARAMETERS ARE
C **
C **     QT          TRANSVERSE MOMENTUM OF W
C **     Y           RAPIDITY
C **
C **  THE GLOBAL VARIABLES SS, QQ, JWTYPE, QM, AND KCOLLID MUST ALSO BE
C **  SET (SEE MODULES COMMON AND CODES).  THE RESULT IS RETURNED
C **  IN SIG1.

      SUBROUTINE DYPT1(QT,Y, SIG1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        EXTERNAL FDEL1,FDEL2
        LOGICAL DOX1,DOX2
      COMMON /DEBUG/ DEBUG,DEBG0
        LOGICAL DEBUG,DEBG0
      COMMON /ERR/ ERREL,NINTVL
      COMMON /STU/ SS,TT,UU,QQ,QM,QMM,qmu
      Real*8 QQCsB
      Common /CsBQQ/ QQCsB
      COMMON /ICODES/ INSTAT,  IQQB1,IQG,IGG,IQQB,IQQ,IALL
      COMMON /MCODES/ MORDER
      DATA IRULE /1/
      DATA GEVXPB /3.8938568E8/
C                                ! CONVERSION 1/(GEV**4) --> PB/(GEV**2)

      IF ( .NOT. ((INSTAT.EQ.IQQB1).OR.(INSTAT.EQ.IQG)
     1            .OR.(INSTAT.EQ.IALL)) ) THEN
        SIG1 = 0
        RETURN
      END IF
CsB >>>
      QQ = QQCsB
c      Print*, ' DYPT1 >>> QQ = ', QQ
CsB <<<
      CALL INTSET(QT,Y, ZX1MIN,ZX2MIN,DOX1,DOX2)
CCPY      ALPHS = ALPHAS(QMu)
      ALPHS = ALPHAS(QMu**2)
C     ! REMAINING OVERALL FACTORS
      FAC1 = GEVXPB * ALPHS
C     ! SET GLOBAL FLAG FOR 1ST ORDER.
      MORDER = 1
      SIG1 = 0
      IF (DOX1) THEN
C       ! D(Z1) DELTA(S2) INTERGAL
        SINTA = AGAUSS(FDEL1,ZX1MIN,1.D0,NINTVL,ERREL)
        SIG1 = SIG1 + FAC1*SINTA
        IF (DEBUG) WRITE(*,*) 'D(Z1)DELTA(S2): ',FAC1*SINTA
      END IF
      IF (DOX2) THEN
C       ! D(Z2) DELTA(S2) INTERGA;
        SINTB = AGAUSS(FDEL2,ZX2MIN,1.D0,NINTVL,ERREL)
        SIG1 = SIG1 + FAC1*SINTB
        IF (DEBUG) WRITE(*,*) 'D(Z2)DELTA(S2): ',FAC1*SINTB
      END IF
      RETURN
      END


C ** DYPT2 -- SIMILAR TO DYPT1 BUT RETURNS THE SECOND-ORDER CONTRIBUTION.

      SUBROUTINE DYPT2(QT,Y, SIG2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        EXTERNAL FSIG1A,FSIG2A,FDEL1,FDEL2
        LOGICAL DOX1,DOX2
      COMMON /DEBUG/ DEBUG,DEBG0
        LOGICAL DEBUG,DEBG0
      COMMON /ERR/ ERREL,NINTVL
      COMMON /STU/ SS,TT,UU,QQ,QM,QMM,qmu
      COMMON /FLAGS/ DS2FLG,NRMFLG,S2AFLG,FACFLG
      Real*8 QQCsB
      Common /CsBQQ/ QQCsB
        LOGICAL DS2FLG,NRMFLG,S2AFLG,FACFLG
      COMMON /MCODES/ MORDER
      DATA GEVXPB /3.8938568E8/
      DATA PI /3.1415 92653 58979 32384 62643/

C                                ! CONVERSION 1/(GEV**2) --> PB/(GEV**4)
      DATA IRULE /1/

CsB >>>
      QQ = QQCsB
c      Print*, ' DYPT1 >>> QQ = ', QQ
CsB <<<

      CALL INTSET(QT,Y, ZX1MIN,ZX2MIN,DOX1,DOX2)
CCPY      ALPHS = ALPHAS(QMu)
      ALPHS = ALPHAS(QMu**2)
CZL
C     ! REMAINING OVERALL FACTORS.
      FAC2 = GEVXPB * ALPHS**2/2/PI
C     ! SET GLOBAL FLAG FOR 2ND ORDER.
      MORDER = 2
      SIG2 = 0

      IF (DS2FLG.OR.FACFLG) THEN
        IF (DOX1) THEN
C         ! D(Z1) DELTA(S2) INTEGRAL
          SINTC = AGAUSS(FDEL1,ZX1MIN,1.D0,NINTVL,ERREL)
          SIG2 = SIG2 + FAC2*SINTC
          IF (DEBUG) WRITE(*,*) 'D(Z1)DELTA(S2): ',FAC2*SINTC
        END IF
        IF (DOX2) THEN
C         ! D(Z2) DELTA(S2) INTEGRAL
          SINTD = AGAUSS(FDEL2,ZX2MIN,1.D0,NINTVL,ERREL)
          SIG2 = SIG2 + FAC2*SINTD
          IF (DEBUG) WRITE(*,*) 'D(Z2)DELTA(S2): ',FAC2*SINTD
        END IF
      END IF
      IF (NRMFLG.OR.S2AFLG.OR.FACFLG) THEN
        IF (DOX1) THEN
C         ! D(Z1)D(S2) INTEGRAL
          SINTA = AGAUSS(FSIG1A,ZX1MIN,1.D0,NINTVL,ERREL)
          SIG2 = SIG2 + FAC2*SINTA
          IF (DEBUG) WRITE(*,*) 'D(Z1)D(S2): ',FAC2*SINTA
        END IF
        IF (DOX2) THEN
C         ! D(Z2)D(S2) INTEGRAL
          SINTB = AGAUSS(FSIG2A,ZX2MIN,1.D0,NINTVL,ERREL)
          SIG2 = SIG2 + FAC2*SINTB
          IF (DEBUG) WRITE(*,*) 'D(Z2)D(S2): ',FAC2*SINTB
        END IF
      END IF
      RETURN
      END


C **  XPART, XPRTQQ -- ACT AS A SWITCHYARD TO CALL THE APPROPRIATE ROUTINES
C **  IN MODULE WPART BASED ON THE GLOBAL FLAGS INSTAT, MORDER, LTERM,
C **  AND JWTYPE.  (SEE COMMON.H AND CODES.H FOR DESCRIPTIONS.)

      FUNCTION XPART(S,T,U,S2,A)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        EXTERNAL QQBW,QQBZ,QQW,QQZ
      COMMON /ICODES/ INSTAT,  IQQB1,IQG,IGG,IQQB,IQQ,IALL
      COMMON /JCODES/ JWTYPE,  JWP,JWM,JZ,JPH
      COMMON /MCODES/ MORDER
CCPY===========================================
      INTEGER IDO_CBAR
      COMMON /CBAR_ONLY/ IDO_CBAR      
      
      IF (INSTAT.EQ.IQQB1) THEN
        XPART = QQB1(S,T,U,S2,A)
        IF ((MORDER.EQ.2) .AND. (JWTYPE.EQ.JZ))
     1    XPART = XPART + QQBZX(S,T,U,S2,A)
      ELSE IF (INSTAT.EQ.IQG) THEN
        XPART = QG(S,T,U,S2,A)
        IF ((MORDER.EQ.2) .AND. (JWTYPE.EQ.JZ))
     1    XPART = XPART + QGZX(S,T,U,S2,A)
      ELSE IF (INSTAT.EQ.IGG) THEN
        XPART = GG(S,T,U,S2,A)
      ELSE IF (INSTAT.EQ.IQQB) THEN
        XPART = XPRTQQ(S,T,U,S2,A,QQBW,QQBZ)
      ELSE IF (INSTAT.EQ.IQQ) THEN
        XPART = XPRTQQ(S,T,U,S2,A,QQW,QQZ)
      ELSE IF (INSTAT.EQ.IALL) THEN
C****************************************************************      
        IF (MORDER.EQ.1) THEN  
          
          IF(IDO_CBAR.EQ.1) THEN
            XPART = QG(S,T,U,S2,A)
          ELSE
            XPART = QQB1(S,T,U,S2,A) + QG(S,T,U,S2,A)
          ENDIF
            
        ELSE 
        
          IF(IDO_CBAR.EQ.1) THEN
            XPART = QG(S,T,U,S2,A)
     1        + GG(S,T,U,S2,A) 
          ELSE
            XPART = QQB1(S,T,U,S2,A) + QG(S,T,U,S2,A)
     1        + GG(S,T,U,S2,A) + XPRTQQ(S,T,U,S2,A,QQBW,QQBZ)
     2        + XPRTQQ(S,T,U,S2,A,QQW,QQZ)
          ENDIF

          IF (JWTYPE.EQ.JZ) THEN
            IF(IDO_CBAR.EQ.1) THEN
              XPART = XPART + QGZX(S,T,U,S2,A)
            ELSE 
              XPART = XPART + QQBZX(S,T,U,S2,A) + QGZX(S,T,U,S2,A)
            ENDIF
          END IF
          
        ENDIF

C****************************************************************
      ELSE
        CALL ERROR
      END IF
      RETURN
      END


      FUNCTION XPRTQQ(S,T,U,S2,A,QW,QZ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        EXTERNAL QW,QZ
      COMMON /JCODES/ JWTYPE,  JWP,JWM,JZ,JPH
      COMMON /LCODES/ LTERM, L2DIM,LS20,LWS20,LDS2
      COMMON /MCODES/ MORDER
      IF (MORDER.EQ.1) CALL ERROR
      IF (JWTYPE.EQ.JWP.or.JWTYPE.EQ.JWM) THEN
        XPRTQQ = QW(S,T,U,S2,A)
      ELSE IF (JWTYPE.EQ.JZ) THEN
        XPRTQQ = QZ(S,T,U,S2,A)
      ELSE
        CALL ERROR
      END IF
      RETURN
      END
      
CCPY==================================
          
C **************************************************************************    
C         
C  MODULE: WPART                                 Version 3.0
C  -------------    
C                   by Peter Arnold
C                   Ref:  P. Arnold and H. Reno, Nucl. Phys. B319 (1989) 37
C                   Ref:  R.J. Gonsalves, J. Pawlowski and C.-F. Wai,
C                         Phys. Rev. D40 (1989) 2245
C
C  Changes from Version 1 to 2: GETQQZ, QQBZ, QQZ, QQBZX, QGZX changed to
C    include the extra diagrams of Gonsalves et al. for Z's.
C  Changes from Version 2 to 3: none
C
C                                    
C      THIS MODULE COMBINES THE PARTON CROSS-SECTIONS OF MODULE WFORM 
C  THE FACTORIZATION CONVERSIONS OF MODULE WFAC, AND THE APPROPRIATE  
C  STRUCTURE FUNCTIONS.       
C      WE HAVE THROUGHOUT IGNORED TERMS OF ORDER (B QUARK DISTRIBUTION)         
C  TIMES (B MIXING ANGLE).    
C         
C  SUMMARY OF ROUTINES        
C  -------------------        
C         
C    QQB1      Q-QBAR:  [Q+QBAR-->G+G+GAMMA] + [VIRTUAL] + (F1+F2)*(F1+F2)      
C    QQBW,QQBZ Q-QBAR:  EVERYTHING ELSE 
C    QG        Q-G  
C    GG        G-G  
C    QQW,QQZ   Q-Q
C
C    QQBZX     Q-QBAR:  THE VIRTUAL DIAGRMS PECULIAR TO Z'S BETWEEN THRESHOLDS
C    QGZX      Q-G   :  DITTO  
C         
C      IN THE PARAMETERS FOR THESE ROUTINES, A IS THE UPPER LIMIT OF THE        
C  S2 INTEGRATION.  S2A IS EQUAL TO THE VALUE OF THE S2 INTEGRATION VARIABLE.   
C  THE PARAMETER S2 ITSELF IS USUALLY THE SAME, BUT IT IS TO BE SET TO ZERO     
C  WHEN ASKING FOR THE S2=0 PIECE OF [1/S2] OR [LOG(S2)/S2] TERMS.    
C      THESE ROUTINES WILL IN GENERAL RETURN THE FOLLOWING FOR DIFFERENT        
C  VALUES OF GLOBAL FLAGS SET BY THE CALLER       
C         
C    MORDER=1                   1ST-ORDER         
C    MORDER=2, LTERM=LDS2       2ND-ORDER: DELTA(S2) PIECE  
C    MORDER=2, LTERM=LS20       2ND-ORDER: COEFFICIENT OF [1/S2] AT S2=0
C    MORDER=2, LTERM=LWS20      2ND-ORDER: COEFFICIENT OF [LOG(S20/S2] AT S2=0
C    MORDER=2, LTERM=L2DIM      2ND-ORDER: EVERYTHING ELSE  
C
C  ONE SHOULD BE CAREFUL, HOWEVER, AS SOME OF THE ROUTINES ASSUME MORDER=2.
C      OTHER GLOBAL VARIABLES THAT MUST BE SET ARE: JWTYPE,NRMFLG,S2AFLG, 
C  DS2FLG,FACFLG,QQ,QMM; AND THE PARTON STRUCTURE FUNCTIONS MUST HAVE 
C  BEEN STORED IN COMMON.  (SEE MODULES COMMON AND CODES)       
C         
C  OTHER MODULES REFERENCED: COMMON, CODES, WFORM, WFAC 
C         
C  OVERALL NORMALIZATIONS     
C  ----------------------     
C    THIS MODULE MULTIPLIES THE CROSS-SECTIONS OF WFORM AND WFAC BY   
C  THE APPROPRIATE  
C         
C                        G(GAMMA)**2 ,  
C         
C  SUMS OVER FINAL FLAVORS, AND INCLUDES FACTORS OF 1/2 FOR IDENTICAL 
C  FINAL PARTICLES WHEN RELEVANT.  THIS LEAVES FOR THE CALLER THE     
C  RESPONSIBILITY OF
C         
C          1ST ORDER RESULTS         1/4/NC/S * ALPHA(S)    
C          2ND ORDER RESULTS         1/4/NC/S * ALPHA(S)**2/2/PI      
C         
C  NAMING CONVENTIONS         
C  ------------------         
C         
C    COMBINATIONS OF PARTON DISTRIBUTION FUNCTIONS HAVE NAMES LIKE    
C         
C            UDB      U(1,X1)*DBAR(2,X2)
C            DBG      DBAR(1,X1)*G(2,X2)
C         
C  BOTTOM (B) IS DISTINGUISHED FROM BAR (B) BY CAPITILIZATION.  THE   
C  ABOVE EXAMPLES ACTUALLY CONTAIN SUMS OVER ALL UP AND DOWN FLAVORS. 
C  THE SUFFIX 1 INDICATES A RESTRICTION OF THE SUM BY A DELTA-FN.     
C  IN FLAVOR SPACE.  THE SUFFIX K INDICATES THAT EACH TERM SHOULD BE  
C  ACCOMPANIED BY THE SQUARE OF THE CORRESPONDING KM-MATRIX ELEMENT.  
C  FINALLY, UU0 IS USED FOR U(1,X1)*U(2,X2) TO DISTINGUISH IT FROM    
C  THE MANDELSTAM VARIABLE DESIGNATED UU.         
C    PARTON CROSS-SECTIONS USE THE SAME BASIC NAMING SCHEME BUT ARE   
C  PREFIXED BY AN X.  COUPLING CONSTANTS MAY BE INCLUDED EITHER WITH THE        
C  DISTRIBUTION FNS. OR THE CROSS-SECTIONS AS IS CONVENIENT.
C         
C **************************************************************************    
          
          
C **  ERROR -- A ROUTINE CALLED WHENEVER THERE'S A SITUATION THAT     
C **  SHOULDN'T HAPPEN, INDICATING EITHER A BUG OR THAT A ROUTINE     
C **  HAS BEEN IMPROPERLY CALLED.       
          
CCPY      SUBROUTINE ERROR        
c      IMPLICIT DOUBLE PRECISION (A-H,O-Z)         
c      WRITE(*,*) 'PROGRAM ERROR.'       
c      STOP
c      END 
          
          
C **  QQBGET -- CALLS THE APPROPRAITE ROUTINES TO GET THE DESIRED PIECE         
C **  (DETERMINED BY MORDER AND LTYPE) OF THE MS-BAR PARTON CROSS-SECTION       
C **  USED BY QQB1.  IN NORMAL OPERATION, ALL THE "FLG" VARIABLES WILL
C **  BE SET TO TRUE.         
          
      FUNCTION QQBGET(S,T,U,S2,A)   
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)       
        EXTERNAL QQ1X,QQS2A,QQWS2A,QQDELT,QQLEAD  
      COMMON /STU/ SS,TT,UU,QQ,QM,QMM,qmu
      COMMON /FLAGS/ DS2FLG,NRMFLG,S2AFLG,FACFLG  
        LOGICAL DS2FLG,NRMFLG,S2AFLG,FACFLG       
      COMMON /LCODES/ LTERM, L2DIM,LS20,LWS20,LDS2      
      COMMON /MCODES/ MORDER  
      QQBGET = 0    
      IF (MORDER.EQ.1) THEN   
        QQBGET = QQLEAD(S,T,U)
      ELSE IF (LTERM.EQ.L2DIM) THEN     
        CALL INIT(S,T,U,S2)   
        IF (NRMFLG) QQBGET = QQ1X(S,T,U,S2)       
        IF (S2AFLG) QQBGET =  
     1    QQBGET + (QQS2A(S,T,U,S2)+QQWS2A(S,T,U,S2)*LOG(S2/QMM))/S2  
        CALL SWITCH(.TRUE.)   
        IF (NRMFLG) QQBGET = QQBGET + QQ1X(S,U,T,S2)        
      ELSE IF ((LTERM.EQ.LS20) .AND. S2AFLG) THEN 
        CALL INIT(S,T,U,0.D0) 
        QQBGET = QQS2A(S,T,U,0.D0)
      ELSE IF ((LTERM.EQ.LWS20) .AND. S2AFLG) THEN
        CALL INIT(S,T,U,0.D0)
        QQBGET = QQWS2A(S,T,U,0.D0)   
      ELSE IF ((LTERM.EQ.LDS2) .AND. DS2FLG) THEN 
        CALL INIDS2(S,T,U,A)  
        QQBGET = QQDELT(S,T,U)
      END IF        
      RETURN        
      END 
          
          
C **  QQB1 -- THIS ROUTINE SUMS Q+QBAR AND QBAR+Q CONTRIBUTIONS.  EITHER        
C **  W OR Z PRODUCTION IS HANDLED DEPENDING ON JTYPE.      
          
      FUNCTION QQB1(S,T,U,S2,A)     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)         
      COMMON /FLAGS/ DS2FLG,NRMFLG,S2AFLG,FACFLG  
        LOGICAL DS2FLG,NRMFLG,S2AFLG,FACFLG       
      COMMON /STRUCT/         
     1         U1,UB1,D1,DB1,ST1,STB1,C1,CB1,B1,BB1,G1      
     2        ,U2,UB2,D2,DB2,ST2,STB2,C2,CB2,B2,BB2,G2      
      COMMON /WPARM/ GF,ZMASS,WMASS,GGW,GGU,GGD,GVU,GVD,GAU,GAD       
      COMMON /CABIBO/ CC11,CC12,CC21,CC22         
      COMMON /JCODES/ JWTYPE,  JWP,JWM,JZ,JPH
      COMMON /MCODES/ MORDER  
CCPY===========================================
      INTEGER IDO_CBAR
      COMMON /CBAR_ONLY/ IDO_CBAR      
          
      XQQB = QQBGET(S,T,U,S2,A)     
      IF ((MORDER.EQ.2) .AND. FACFLG)   
     1  XQQB = XQQB + QQB1FC(S,T,U,S2,A)      
          
      IF (JWTYPE.EQ.JWP.or.JWTYPE.EQ.JWM) THEN  
        IF(IDO_CBAR.EQ.1) THEN
          PARTFN = 0.0
        ELSE
          UDBK = U1*(CC11*DB2+CC12*STB2)  
     >      + C1*(CC21*DB2+CC22*STB2)  
          DBUK = U2*(CC11*DB1+CC12*STB1)  
     >      + C2*(CC21*DB1+CC22*STB1)  
          PARTFN = GGW*(UDBK+DBUK)        
        ENDIF                   !if (ido_cbar.eq.1)
            
      ELSE IF (JWTYPE.EQ.JZ) THEN       
        IF(IDO_CBAR.EQ.1) THEN
          PARTFN = 0.0
        ELSE
          UUB1 = U1*UB2 + C1*CB2
          UBU1 = U2*UB1 + C2*CB1
          DDB1 = D1*DB2 + ST1*STB2 + B1*BB2         
          DBD1 = D2*DB1 + ST2*STB1 + B1*BB2         
CZL
!          PARTFN = GGU*(UUB1+UBU1) + GGD*(DDB1+DBD1)
          PARTFN =  GGD*(DDB1+DBD1)
        ENDIF
            
      ELSE
CCPY
        write(*,*) 'This JWTYPE is not yet implemented.'
        CALL ERROR  
      END IF        
CZL
      QQB1 = PARTFN * XQQB  
      
      RETURN        
      END 
          
          
C **  QGGET -- SIMILAR TO QQBGET BUT FOR Q-G AND G-Q PROCESSES.  THE  
C **  RELEVANT PIECES OF THE MS-BAR PARTON CROSS-SECTION ARE RETURNED 
C **  AS XQG AND XGQ RESPECTIVELY.      
          
      SUBROUTINE QGGET(XQG,XGQ,S,T,U,S2,A)    
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)         
      COMMON /STU/ SS,TT,UU,QQ,QM,QMM,qmu
      COMMON /FLAGS/ DS2FLG,NRMFLG,S2AFLG,FACFLG  
        LOGICAL DS2FLG,NRMFLG,S2AFLG,FACFLG       
      COMMON /LCODES/ LTERM, L2DIM,LS20,LWS20,LDS2      
      COMMON /MCODES/ MORDER  
      XQG = 0       
      XGQ = 0       
      IF (MORDER.EQ.1) THEN   
        XQG = QGLEAD(S,T,U)   
        XGQ = QGLEAD(S,U,T)   
      ELSE IF (LTERM.EQ.L2DIM) THEN     
        CALL INIT(S,T,U,S2)   
        IF (NRMFLG) XQG = QG1(S,T,U,S2) 
        IF (S2AFLG) XQG =     
     1    XQG + (QGS2A(S,T,U,S2)+QGWS2A(S,T,U,S2)*LOG(S2/QMM))/S2     
        CALL SWITCH(.TRUE.)   
        IF (NRMFLG) XGQ = QG1(S,U,T,S2) 
        IF (S2AFLG) XGQ =     
     1    XGQ + (QGS2A(S,U,T,S2)+QGWS2A(S,U,T,S2)*LOG(S2/QMM))/S2     
      ELSE IF ((LTERM.EQ.LS20) .AND. S2AFLG) THEN 
        CALL INIT(S,T,U,0.D0) 
        XQG = QGS2A(S,T,U,0.D0)
        CALL SWITCH(.TRUE.)   
        XGQ = QGS2A(S,U,T,0.D0)
      ELSE IF ((LTERM.EQ.LWS20) .AND. S2AFLG) THEN
        CALL INIT(S,T,U,0.D0)
        XQG = QGWS2A(S,T,U,0.D0)
        CALL SWITCH(.TRUE.)
        XGQ = QGWS2A(S,U,T,0.D0)
      ELSE IF ((LTERM.EQ.LDS2) .AND. DS2FLG) THEN 
        CALL INIDS2(S,T,U,A)  
        XQG = QGDELT(S,T,U)   
        CALL SWTDS2(.TRUE.)   
        XGQ = QGDELT(S,U,T)   
      END IF        
      RETURN        
      END 
          
          
C ** QG -- THIS ROUTINE SUMS Q+G, G+Q, QBAR+G, AND G+QBAR CONTRIBUTIONS.        
C ** EITHER W OR Z PRODUCTION IS HANDLED DEPENDING ON JTYPE.
          
      FUNCTION QG(S,T,U,S2,A)       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)         
      COMMON /FLAGS/ DS2FLG,NRMFLG,S2AFLG,FACFLG  
        LOGICAL DS2FLG,NRMFLG,S2AFLG,FACFLG       
      COMMON /STRUCT/         
     1         U1,UB1,D1,DB1,ST1,STB1,C1,CB1,B1,BB1,G1      
     2        ,U2,UB2,D2,DB2,ST2,STB2,C2,CB2,B2,BB2,G2      
      COMMON /WPARM/ GF,ZMASS,WMASS,GGW,GGU,GGD,GVU,GVD,GAU,GAD       
      COMMON /JCODES/ JWTYPE,  JWP,JWM,JZ,JPH
      COMMON /MCODES/ MORDER  

CCPY===========================================
      INTEGER IDO_CBAR
      COMMON /CBAR_ONLY/ IDO_CBAR      
          
      CALL QGGET(XQG,XGQ, S,T,U,S2,A)         
      IF ((MORDER.EQ.2) .AND. FACFLG) THEN        
        XQG = XQG + QGFC(S,T,U,S2,A)
        XGQ = XGQ + QGFC(S,U,T,S2,A)
      END IF        
          
      IF (JWTYPE.EQ.JWP.or.JWTYPE.eq.JWM) THEN  
cpn2007 This hotfix works both for W+ and W- (see also the choice of the
cpn2007 PDF's for the W+ and W- case in the subroutine GetStr)

        IF(IDO_CBAR.EQ.1) THEN
          UG =0.0
          GU =0.0
          DBG = STB1*G2   
          GDB = STB2*G1   
        ELSE
          UG  = (U1+C1)*G2        
          GU  = (U2+C2)*G1        
          DBG = (DB1+STB1)*G2   
          GDB = (DB2+STB2)*G1   
        ENDIF                   !ido_cbar
        QG = GGW * ( (UG+DBG)*XQG + (GU+GDB)*XGQ )
      ELSE IF (JWTYPE.EQ.JZ) THEN       
        IF(IDO_CBAR.EQ.1) THEN
          UG =0.0
          GU =0.0
          DBG =0.0   
          GDB =0.0
          DG = 0.0
          GD =0.0
          UBG = CB1*G2
          GUB = CB2*G1   
        ELSE
          UG  = (U1+C1)*G2        
          GU  = (U2+C2)*G1        
          DBG = (DB1+STB1+BB1)*G2         
          GDB = (DB2+STB2+BB2)*G1         
          DG  = (D1+ST1+B1)*G2  
          GD  = (D2+ST2+B2)*G1  
          UBG = (UB1+CB1)*G2    
          GUB = (UB2+CB2)*G1    
        ENDIF
CZL
!        QG = (GGU*(UG+UBG) + GGD*(DG+DBG)) * XQG  
!     1    + (GGU*(GU+GUB) + GGD*(GD+GDB)) * XGQ  
        QG = (GGD*(DG+DBG)) * XQG  
     1    + (GGD*(GD+GDB)) * XGQ  
      END IF        
      RETURN        
      END 
          
          
C ** GG -- EITHER W OR Z PRODUCTION IS HANDLED DEPENDING ON JTYPE.    
          
      FUNCTION GG(S,T,U,S2,A)       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)         
        LOGICAL FIRST         
      COMMON /FLAGS/ DS2FLG,NRMFLG,S2AFLG,FACFLG  
        LOGICAL DS2FLG,NRMFLG,S2AFLG,FACFLG       
      COMMON /CLEBSH/ CF,CNP,TR,NC,NU,ND,NUD     
      COMMON /STRUCT/         
     1         U1,UB1,D1,DB1,ST1,STB1,C1,CB1,B1,BB1,G1      
     2        ,U2,UB2,D2,DB2,ST2,STB2,C2,CB2,B2,BB2,G2      
      COMMON /WPARM/ GF,ZMASS,WMASS,GGW,GGU,GGD,GVU,GVD,GAU,GAD       
      COMMON /JCODES/ JWTYPE,  JWP,JWM,JZ,JPH
      COMMON /LCODES/ LTERM, L2DIM,LS20,LWS20,LDS2      
      COMMON /MCODES/ MORDER  
CCPY===========================================
      INTEGER IDO_CBAR
      COMMON /CBAR_ONLY/ IDO_CBAR      
          
          
      XGG = 0       
      IF ((LTERM.EQ.L2DIM) .AND. NRMFLG) THEN     
        CALL INIT(S,T,U,S2)   
        XGG = GG1X(S,T,U,S2)  
        CALL SWITCH(.TRUE.) 
        XGG = XGG + GG1X(S,U,T,S2)    
      END IF        
          
      IF ((MORDER.EQ.2) .AND. FACFLG)   
     1  XGG = XGG + GGFC(S,T,U,S2,A)
          
      IF (JWTYPE.EQ.JWP.or.JWTYPE.eq.JWM) THEN  
            IF(IDO_CBAR.EQ.1) NUD=1

        GG = NUD*GGW * G1*G2 * XGG      
      ELSE IF (JWTYPE.EQ.JZ) THEN       
        IF(IDO_CBAR.EQ.1) THEN
          NU=1
          ND=0
        ENDIF
        
CZL NNLO
!        GG = (NU*GGU+ND*GGD) * G1*G2 * XGG        
        GG = (ND*GGD) * G1*G2 * XGG        
!	GG=0
      ELSE
        CALL ERROR  
      END IF        
      RETURN        
      END 
          
          
C **  QQBW -- THIS ROUTINE IS SPECIFIC TO W PRODUCTION.     
          
      FUNCTION QQBW(S,T,U,S2,A)     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)         
      COMMON /CLEBSH/ CF,CNP,TR,NC,NU,ND,NUD     
      COMMON /FLAGS/ DS2FLG,NRMFLG,S2AFLG,FACFLG  
        LOGICAL DS2FLG,NRMFLG,S2AFLG,FACFLG       
      COMMON /STRUCT/         
     1         U1,UB1,D1,DB1,ST1,STB1,C1,CB1,B1,BB1,G1      
     2        ,U2,UB2,D2,DB2,ST2,STB2,C2,CB2,B2,BB2,G2      
      COMMON /WPARM/ GF,ZMASS,WMASS,GGW,GGU,GGD,GVU,GVD,GAU,GAD       
      COMMON /CABIBO/ CC11,CC12,CC21,CC22         
      COMMON /LCODES/ LTERM, L2DIM,LS20,LWS20,LDS2      
          
      IF ((LTERM.EQ.L2DIM).AND.NRMFLG) THEN       
        CALL INIT(S,T,U,S2)   
        FF3434 = F3434X(S,T,U,S2)       
        FF5656 = F5656(S,T,U,S2)        
        FF1256 = F1256(S,T,U,S2)        
        FF3456 = F3456(S,T,U,S2)        
        CALL SWITCH(.TRUE.)   
        FF3434 = FF3434 + F3434X(S,U,T,S2)        
        FF7878 = F5656(S,U,T,S2)        
        FF1278 = F1256(S,U,T,S2)        
        FF3478 = F3456(S,U,T,S2)        
      ELSE
        FF3434 = 0  
        FF5656 = 0  
        FF1256 = 0  
        FF3456 = 0  
        FF7878 = 0  
        FF1278 = 0  
        FF3478 = 0  
      END IF        
          
      IF (FACFLG) THEN        
        FF5656 = FF5656 + F56FC(S,T,U,S2,A)   
        FF7878 = FF7878 + F56FC(S,U,T,S2,A)   
      END IF        
cpn2007 Note: this code applies both for W+ and W- case, for which PDF's
cpn2007 are chosen accordingly in the subroutine GetStr
      UDBK = U1*(CC11*DB2+CC12*STB2)    
     1     + C1*(CC21*DB2+CC22*STB2)    
      DBUK = U2*(CC11*DB1+CC12*STB1)    
     1     + C2*(CC21*DB1+CC22*STB1)    
      UDB  = (U1+C1)*(DB2+STB2)         
      DBU  = (U2+C2)*(DB1+STB1)         
      UBB  = (U1+C1)*BB2      
      BBU  = (U2+C2)*BB1      
      UUB1 = U1*UB2 + C1*CB2  
      UBU1 = U2*UB1 + C2*CB1  
      UUB  = (U1+C1)*(UB2+CB2)
      UBU  = (U2+C2)*(UB1+CB1)
      DDB1 = D1*DB2 + ST1*STB2
      DBD1 = D2*DB1 + ST2*STB1
      DDB  = (D1+ST1)*(DB2+STB2)        
      DBD  = (D2+ST2)*(DB1+STB1)        
      BDB  = B1*(DB2+STB2)    
      DBB  = B2*(DB1+STB1)    
          
C ! Q1+QBAR2-->W+ RELATED TO QBAR1+Q2-->W+ BY CP+ROTATION+ISOROTATION 
C ! TAKING U1,D1,U2,D2 --> DB1,UB1,DB2,UB2.       
C ! OR FLIP T<--> U TAKING U1,D1 <--> U2,D2.      
          
      XUDB = UDBK*(FF1256 + FF1278) + UDB*(FF5656 + FF7878) 
     1       + UBB*FF5656     
      XDBU = DBUK*(FF1256 + FF1278) + DBU*(FF5656 + FF7878) 
     1       + BBU*FF7878     
      XUUB = UUB1*(NUD*FF3434 + FF3456) + UUB*FF5656        
      XUBU = UBU1*(NUD*FF3434 + FF3478) + UBU*FF7878        
      XDDB = DDB1*(NUD*FF3434 + FF3478) + DDB*FF7878 + BDB*FF7878     
      XDBD = DBD1*(NUD*FF3434 + FF3456) + DBD*FF5656 + DBB*FF5656     
          
      QQBW = GGW * (XUDB+XDBU+XUUB+XUBU+XDDB+XDBD)
      RETURN        
      END 
          
          
C **  GETQQZ -- THIS ROUTINE COLLECTS AND RETURNS THE TERMS COMMON    
C **  TO Q+QBAR --> Q+QBAR+Z AND Q+Q --> Q+Q+Z.  IFLAG SHOULD BE +1
C **  OR -1 AND GIVES THE SIGN OF THE VECTORIAL CROSS-TERMS: +1 FOR
C **  Q+QBAR AND -1 FOR Q+Q.  THIS ROUTINE HANDLES THE CALL TO INIT.
          
      SUBROUTINE GETQQZ(IFLAG,S,T,U,S2,A,XUD,XDU,XUU,XDD)     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)         
      COMMON /FLAGS/ DS2FLG,NRMFLG,S2AFLG,FACFLG  
        LOGICAL DS2FLG,NRMFLG,S2AFLG,FACFLG       
      COMMON /WPARM/ GF,ZMASS,WMASS,GGW,GGU,GGD,GVU,GVD,GAU,GAD       
      COMMON /LCODES/ LTERM, L2DIM,LS20,LWS20,LDS2      
          
      IF ((LTERM.EQ.L2DIM) .AND. NRMFLG) THEN     
        CALL INIT(S,T,U,S2)   
        FF5656 = F5656(S,T,U,S2)        
        FV5678 = F5678X(S,T,U,S2)       
        FA5678 = F5678A(S,T,U,S2)       
        CALL SWITCH(.TRUE.)   
        FF7878 = F5656(S,U,T,S2)        
        FV5678 = FV5678 + F5678X(S,U,T,S2)        
        FA5678 = FA5678 + F5678A(S,U,T,S2)        
      ELSE
        FF5656 = 0                   
        FV5678 = 0  
        FA5678 = 0  
        FF7878 = 0  
      END IF        
          
      IF (FACFLG) THEN        
        FF5656 = FF5656 + F56FC(S,T,U,S2,A)   
        FF7878 = FF7878 + F56FC(S,U,T,S2,A)   
      END IF

      IF (ABS(IFLAG).NE.1) CALL ERROR        
          
CZL
!      XUD = GGU*FF5656 + GGD*FF7878
!     1        + IFLAG*GVU*GVD*FV5678 + GAU*GAD*FA5678 
!      XDU = GGU*FF7878 + GGD*FF5656
!     2        + IFLAG*GVU*GVD*FV5678 + GAU*GAD*FA5678 
!      XUU  = IFLAG*GVU**2*FV5678 + GAU**2*FA5678 + GGU*(FF5656+FF7878)      
!      XDD  = IFLAG*GVD**2*FV5678 + GAD**2*FA5678 + GGD*(FF5656+FF7878)      
      XUD = GGD*FF7878
     1        + IFLAG*GVU*GVD*FV5678 *(GGD/(GGU+GGD))
     &        + GAU*GAD*FA5678 *(GGD/(GGU+GGD))
      XDU = GGD*FF5656
     2        + IFLAG*GVU*GVD*FV5678 *(GGD/(GGU+GGD))
     &        + GAU*GAD*FA5678 *(GGD/(GGU+GGD))
      XUU  = 0!IFLAG*GVU**2*FV5678 + GAU**2*FA5678 + GGU*(FF5656+FF7878)      
      XDD  = IFLAG*GVD**2*FV5678 + GAD**2*FA5678 + GGD*(FF5656+FF7878)      
          
      RETURN        
      END 
          
          
C ** QQBZ -- THIS ROUTINE IS SPECIFIC TO Z PRODUCTION.      
          
      FUNCTION QQBZ(S,T,U,S2,A)     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)         
      COMMON /STRUCT/         
     1         U1,UB1,D1,DB1,ST1,STB1,C1,CB1,B1,BB1,G1      
     2        ,U2,UB2,D2,DB2,ST2,STB2,C2,CB2,B2,BB2,G2      
      COMMON /FLAGS/ DS2FLG,NRMFLG,S2AFLG,FACFLG  
        LOGICAL DS2FLG,NRMFLG,S2AFLG,FACFLG       
      COMMON /CLEBSH/ CF,CNP,TR,NC,NU,ND,NUD     
      COMMON /WPARM/ GF,ZMASS,WMASS,GGW,GGU,GGD,GVU,GVD,GAU,GAD       
      COMMON /LCODES/ LTERM, L2DIM,LS20,LWS20,LDS2      
          
      CALL GETQQZ(1,S,T,U,S2,A,XUDB,XDBU,XUUB,XDDB)       
      IF ((LTERM.EQ.L2DIM) .AND. NRMFLG) THEN     
        CALL SWITCH(.FALSE.)  
        FF3434 = F3434X(S,T,U,S2)       
        FF1256 = F1256(S,T,U,S2)        
        FF3456 = F3456(S,T,U,S2)
        FA1234 = F1234A(S,T,U,S2)        
        CALL SWITCH(.TRUE.)   
        FF3434 = FF3434 + F3434X(S,U,T,S2)        
        FF1278 = F1256(S,U,T,S2)        
        FF3478 = F3456(S,U,T,S2)
        FA1234 = FA1234 + F1234A(S,U,T,S2)        
      ELSE
        FF3434 = 0  
        FF1256 = 0  
        FF3456 = 0  
        FF1278 = 0  
        FF3478 = 0
        FA1234 = 0  
      END IF        
          
CZL
!      GGSUM = NU*GGU + ND*GGD 
!      XUUB1 = GGU*(FF1256+FF1278+FF3456+FF3478) + GGSUM*FF3434
!     1          + GAU*(NU*GAU+ND*GAD)*FA1234        
!      XDDB1 = GGD*(FF1256+FF1278+FF3456+FF3478) + GGSUM*FF3434
!     1          + GAD*(NU*GAU+ND*GAD)*FA1234        
      GGSUM =  ND*GGD 
      XUUB1 =  GGSUM*FF3434
     1          + GAU*(ND*GAD *(GGD/(GGU+GGD))   )*FA1234        
      XDDB1 = GGD*(FF1256+FF1278+FF3456+FF3478) +  GGSUM*FF3434
     1          + GAD*(NU*GAU *(GGD/(GGU+GGD)) +ND*GAD)*FA1234        
          
      UDB  = (U1+C1)*(DB2+STB2+BB2)     
      DBU  = (U2+C2)*(DB1+STB1+BB1)     
      DUB  = (D1+ST1+B1)*(UB2+CB2)      
      UBD  = (D2+ST2+B2)*(UB1+CB1)      
      UUB1 = U1*UB2 + C1*CB2  
      UBU1 = U2*UB1 + C2*CB1  
      DDB1 = D1*DB2 + ST1*STB2 + B1*BB2 
      DBD1 = D2*DB1 + ST2*STB1 + B2*BB1 
      UUB  = (U1+C1)*(UB2+CB2)
      UBU  = (U2+C2)*(UB1+CB1)
      DDB  = (D1+ST1+B1)*(DB2+STB2+BB2) 
      DBD  = (D2+ST2+B2)*(DB1+STB1+BB1) 
          
C ! Q1+QBAR2-->Z RELATED TO QBAR1+Q2-->W BY CP+ROTATION     
C ! TAKING U1,D1,U2,D2 --> UB1,DB1,UB2,DB2.       
          
CZL NNLO
      QQBZ = (UDB+UBD)*XUDB + (DBU+DUB)*XDBU 
     &     + (UUB1+UBU1)*XUUB1 + (UUB+UBU)*XUUB 
     1     + (DDB1+DBD1)*XDDB1 + (DDB+DBD)*XDDB      
!      QQBZ = (UDB+UBD)*XUDB  + (DBU+DUB)*XDBU 
!     &     + (UUB1+UBU1)*XUUB1 + (UUB+UBU)*XUUB 
!     1     + (DDB1+DBD1)*XDDB1       
!	QQBZ=0
      RETURN        
      END 
          
          
C ** QQBZX,QGZX -- COMPUTES THE EXTRA VIRUAL DIAGRAMS ASSOCIATED WITH Z
C **   PRODUCTION FOR Q+QBAR AND Q+G.      
          
      FUNCTION QQBZX(S,T,U,S2,A)     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)         
      COMMON /STRUCT/         
     1         U1,UB1,D1,DB1,ST1,STB1,C1,CB1,B1,BB1,G1      
     2        ,U2,UB2,D2,DB2,ST2,STB2,C2,CB2,B2,BB2,G2      
      COMMON /FLAGS/ DS2FLG,NRMFLG,S2AFLG,FACFLG  
        LOGICAL DS2FLG,NRMFLG,S2AFLG,FACFLG       
      COMMON /CLEBSH/ CF,CNP,TR,NC,NU,ND,NUD     
      COMMON /WPARM/ GF,ZMASS,WMASS,GGW,GGU,GGD,GVU,GVD,GAU,GAD       
      COMMON /LCODES/ LTERM, L2DIM,LS20,LWS20,LDS2      
          
      IF (.NOT. ((LTERM.EQ.LDS2) .AND. DS2FLG) ) THEN
        QQBZX = 0
        RETURN
      END IF
      CALL INIDS2(S,T,U,A)
      ZX = QQDELZ(S,T,U)     
      UUB1 = U1*UB2 + C1*CB2  
      UBU1 = U2*UB1 + C2*CB1  
      DDB1 = D1*DB2 + ST1*STB2 + B1*BB2 
      DBD1 = D2*DB1 + ST2*STB1 + B2*BB1 
CZL NNLO
!      QQBZX =  (GAU*(UUB1+UBU1)+GAD*(DDB1+DBD1)) * (NU*GAU+ND*GAD) * ZX
      QQBZX =  (GAD*(DDB1+DBD1)) * (ND*GAD) * ZX
     &  + (GAU*(UUB1+UBU1)) * (ND*GAD) * ZX *(GGD/(GGU+GGD))
     &  + (GAD*(DDB1+DBD1)) * (NU*GAU) * ZX *(GGD/(GGU+GGD))
!	QQBZX=0
      RETURN
      END

 
      FUNCTION QGZX(S,T,U,S2,A)     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)         
      COMMON /STRUCT/         
     1         U1,UB1,D1,DB1,ST1,STB1,C1,CB1,B1,BB1,G1      
     2        ,U2,UB2,D2,DB2,ST2,STB2,C2,CB2,B2,BB2,G2      
      COMMON /FLAGS/ DS2FLG,NRMFLG,S2AFLG,FACFLG  
        LOGICAL DS2FLG,NRMFLG,S2AFLG,FACFLG       
      COMMON /CLEBSH/ CF,CNP,TR,NC,NU,ND,NUD     
      COMMON /WPARM/ GF,ZMASS,WMASS,GGW,GGU,GGD,GVU,GVD,GAU,GAD       
      COMMON /LCODES/ LTERM, L2DIM,LS20,LWS20,LDS2      
          
      IF (.NOT. ((LTERM.EQ.LDS2) .AND. DS2FLG) ) THEN
        QGZX = 0
        RETURN
      END IF
      CALL INIDS2(S,T,U,A)
      ZXQG = QGDELZ(S,T,U)
      CALL SWTDS2(.TRUE.)
      ZXGQ = QGDELZ(S,U,T)
      UG  = (U1+C1)*G2        
      GU  = (U2+C2)*G1        
      DBG = (DB1+STB1+BB1)*G2         
      GDB = (DB2+STB2+BB2)*G1         
      DG  = (D1+ST1+B1)*G2  
      GD  = (D2+ST2+B2)*G1  
      UBG = (UB1+CB1)*G2    
      GUB = (UB2+CB2)*G1    
CZL NNLO
!      QGZX = (GAU*(UG+UBG)+GAD*(DG+DBG)) * (GAU*NU+GAD*ND) * ZXQG
!     1     + (GAU*(GU+GUB)+GAD*(GD+GDB)) * (GAU*NU+GAD*ND) * ZXGQ
      QGZX = (GAD*(DG+DBG)) * (GAD*ND) * ZXQG
     1     + (GAD*(GD+GDB)) * (GAD*ND) * ZXGQ
     &     + (GAU*(UG+UBG)) * (GAD*ND) * ZXQG *(GGD/(GGU+GGD))
     &     + (GAU*(GU+GUB)) * (GAD*ND) * ZXGQ *(GGD/(GGU+GGD))
     &     + (GAD*(DG+DBG)) * (GAU*NU) * ZXQG *(GGD/(GGU+GGD))
     &     + (GAD*(GD+GDB)) * (GAU*NU) * ZXGQ *(GGD/(GGU+GGD))
!	QGZX=0
      RETURN
      END
 
          
          
C ** QQW -- THIS ROUTINE IS SPECIFIC TO W PRODUCTION.       
          
      FUNCTION QQW(S,T,U,S2,A)      
C     ! INCLUDES Q-Q AND QBAR-QBAR.     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)         
      COMMON /FLAGS/ DS2FLG,NRMFLG,S2AFLG,FACFLG  
        LOGICAL DS2FLG,NRMFLG,S2AFLG,FACFLG       
      COMMON /STRUCT/         
     1         U1,UB1,D1,DB1,ST1,STB1,C1,CB1,B1,BB1,G1      
     2        ,U2,UB2,D2,DB2,ST2,STB2,C2,CB2,B2,BB2,G2      
      COMMON /WPARM/ GF,ZMASS,WMASS,GGW,GGU,GGD,GVU,GVD,GAU,GAD       
      COMMON /CABIBO/ CC11,CC12,CC21,CC22         
      COMMON /LCODES/ LTERM, L2DIM,LS20,LWS20,LDS2      
          
      IF ((LTERM.EQ.L2DIM) .AND. NRMFLG) THEN     
        CALL INIT(S,T,U,S2)   
        FF5656 = F5656(S,T,U,S2)        
        HH1256 = H1256(S,T,U,S2)        
        HH1278 = H1278(S,T,U,S2)        
        CALL SWITCH(.TRUE.)   
        FF7878 = F5656(S,U,T,S2)        
        HH3478 = H1256(S,U,T,S2)        
      ELSE
        FF5656 = 0  
        HH1256 = 0  
        HH1278 = 0  
        FF7878 = 0  
        HH3478 = 0  
      END IF        
          
      IF (FACFLG) THEN        
        FF5656 = FF5656 + F56FC(S,T,U,S2,A)   
        FF7878 = FF7878 + F56FC(S,U,T,S2,A)   
      END IF        
          
      UU1 = U1*U2 + C1*C2     
      UU0  = (U1+C1)*(U2+C2)  
      UDK = U1*(CC11*D2+CC12*ST2) + C1*(CC21*D2+CC22*ST2)   
      DUK = U2*(CC11*D1+CC12*ST1) + C2*(CC21*D1+CC22*ST1)   
      UD  = (U1+C1)*(D2+ST2)  
      DU  = (U2+C2)*(D1+ST1)  
      UB  = (U1+C1)*B2        
      BU  = (U2+C2)*B1        
          
      DBDB1 = DB1*DB2 + STB1*STB2       
      DBDB  = (DB1+STB1)*(DB2+STB2)     
      DBBB  = (DB1+STB1)*BB2  
      BBDB  = (DB2+STB2)*BB1  
      UBDBK = UB1*(CC11*DB2+CC12*STB2)  
     1      + CB1*(CC21*DB2+CC22*STB2)  
      DBUBK = UB2*(CC11*DB1+CC12*STB1)  
     1      + CB2*(CC21*DB1+CC22*STB1)  
      UBDB  = (UB1+CB1)*(DB2+STB2)      
      DBUB  = (UB2+CB2)*(DB1+STB1)      
                 
C ! Q1+Q2-->W+ RELATED TO QBAR1+QBAR2-->W+ BY CP+ROTATION+ISOROTATION 
C ! TAKING U1,D1,U2,D2 --> DB1,UB1,DB2,UB2.       
          
      XUU  = (UU0+DBDB)*(FF5656 + FF7878) + (UU1+DBDB1)*HH1278        
     1       + DBBB*FF5656 + BBDB*FF7878
      XUD  = (UD+DBUB)*FF5656 + (UDK+DBUBK)*HH1256/2. + UB*FF5656     
      XDU  = (DU+UBDB)*FF7878 + (DUK+UBDBK)*HH3478/2. + BU*FF7878     
          
      QQW = GGW * (XUU + XUD + XDU)     
      RETURN        
      END 
          
          
C ** QQZ -- THIS ROUTINE IS SPECIFIC TO Z PRODUCTION.       
          
      FUNCTION QQZ(S,T,U,S2,A)      
C     ! INCLUDES Q-Q AND QBAR-QBAR.     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)         
      COMMON /FLAGS/ DS2FLG,NRMFLG,S2AFLG,FACFLG  
        LOGICAL DS2FLG,NRMFLG,S2AFLG,FACFLG       
      COMMON /STRUCT/         
     1         U1,UB1,D1,DB1,ST1,STB1,C1,CB1,B1,BB1,G1      
     2        ,U2,UB2,D2,DB2,ST2,STB2,C2,CB2,B2,BB2,G2      
      COMMON /WPARM/ GF,ZMASS,WMASS,GGW,GGU,GGD,GVU,GVD,GAU,GAD       
      COMMON /LCODES/ LTERM, L2DIM,LS20,LWS20,LDS2      
          
      CALL GETQQZ(-1,S,T,U,S2,A,XUD,XDU,XUU,XDD) 
      IF ((LTERM.EQ.L2DIM) .AND. NRMFLG) THEN     
        CALL SWITCH(.FALSE.)  
        HH1256 = H1256(S,T,U,S2)        
        HH1278 = H1278(S,T,U,S2)        
        CALL SWITCH(.TRUE.)   
        HH3478 = H1256(S,U,T,S2)        
        HH3456 = H1278(S,U,T,S2)        
      ELSE
        HH1256 = 0  
        HH1278 = 0  
        HH3478 = 0  
        HH3456 = 0  
      END IF        
          
CZL
      XUU1 = 0!(1./2.)*GGU*(HH1256+HH1278+HH3456+HH3478)      
      XDD1 = (1./2.)*GGD*(HH1256+HH1278+HH3456+HH3478)      
          
      UD  = (U1+C1)*(D2+ST2+B2)         
      DU  = (U2+C2)*(D1+ST1+B1)         
      UU1 = U1*U2 + C1*C2     
      UU0  = (U1+C1)*(U2+C2)  
      DD1 = D1*D2 + ST1*ST2 + B1*B2     
      DD  = (D1+ST1+B1)*(D2+ST2+B2)     
          
      UBDB  = (UB1+CB1)*(DB2+STB2+BB2)  
      DBUB  = (UB2+CB2)*(DB1+STB1+BB1)  
      UBUB1 = UB1*UB2 + CB1*CB2         
      UBUB  = (UB1+CB1)*(UB2+CB2)       
      DBDB1 = DB1*DB2 + STB1*STB2 + BB1*BB2       
      DBDB  = (DB1+STB1+BB1)*(DB2+STB2+BB2)       
          
C ! Q1+Q2-->Z RELATED TO QBAR1+QBAR2-->Z BY CP+ROTATION TAKING        
C !   U1,D1,U2,D2 --> UB1,DB1,UB2,DB2.  
          
CZL NNLO
      QQZ = (UD+UBDB)*XUD + (DU+DBUB)*XDU + (UU1+UBUB1)*XUU1
     1    + (DD1+DBDB1)*XDD1 + (UU0+UBUB)*XUU + (DD+DBDB)*XDD         
!      QQZ = (UD+UBDB)*XUD  + (DU+DBUB)*XDU  
!     &    + (UU1+UBUB1)*XUU1 + (DD1+DBDB1)*XDD1 
!     1    + (UU0+UBUB)*XUU         
!	QQZ=0
      RETURN        
      END 


C======================
CCPY W_FORM.FOR

          
C **************************************************************************    
C         
C  MODULE: WFORM                                        Version 3.0
C  -------------    
C                   by Peter Arnold
C                   Ref:  P. Arnold and H. Reno, Nucl. Phys. B319 (1989) 37
C                   Ref:  R.J. Gonsalves, J. Pawlowski and C.-F. Wai,
C                         Phys. Rev. D40 (1989) 2245
C
C  Changes from Versions 1 to 2: Routines F1234A, QQDELZ, QGDELZ now
C    incorporate the extra diagrams of Gonsalves et al. for Z production.
C  Changes from Version 2 to 2.1: The result of entry GG1X was previously
C    assigned to QQ1X.  This was synonymous to assignment to GG1X in
C    the original version of Vax Fortran that I used and did not affect results.
C    The error has since now corrected.
C  Changes from Version 2.1 to 3.0: none 
C
C         
C      THIS MODULE CONTAINS FUNCTIONS FOR THE EVALUATION OF THE BASIC PARTON    
C  CROSS-SECTIONS S*D(SIGMA)/DT/DU.     
C      BEFORE USING ANY OF THESE FUNCTIONS, INIT OR INITDS2 MUST BE CALLED FOR  
C  EACH NEW SET OF VALUES (S,T,U,S2).  THESE ROUTINES EVALUATE AND STORES THE   
C  LOGARITHMS THAT WILL BE NEEDED FOR THE CROSS-SECTIONS.  THIS IS DONE         
C  SEPARATELY SO THAT IT NEED BE DONE ONLY ONCE FOR EACH (S,T,U,S2).  USE       
C  INITDS2 FOR THE DELTA(S2) PEICES OF SECOND-ORDER CROSS-SECTIONS AND
C  INIT FOR THE OTHER PIECES OF SECOND-ORDER CROSS-SECTIONS.
C      OCCASIONALLY, ONE WILL WANT TO SWITCH T AND U IN A FORMULA.  TO DO       
C  SO, CALL SWITCH(.TRUE.) ANYTIME AFTER INIT.  THIS WILL SWITCH THE  
C  APPROPRIATE LOGARITHMS.  THEN CALL THE CROSS-SECTION ROUTINE YOU WANT        
C  WITH THE ARGUMENTS T AND U SWITCHED.  FOR DELTA(S2) PIECES, USE    
C  SWITCHDS2(.TRUE.).  CALLS TO SWITCH(.FALSE.) AND SWITCHDS2(.FALSE.)
C  WILL ENSURE THAT T AND U ARE *NOT* SWITCHED IN THE LOGARITHMS.     
C  FOR EXAMPLE, ONE PIECE OF THE Q-QBAR CROSS-SECTIONS COMES FROM     
C  ADDING QQ1X WITH U<-->T.  THIS IS DONE BY      
C         
C            CALL INIT(S,T,U,S2)        
C            PIECE = QQ1X(S,T,U,S2)     
C            CALL SWITCH(.TRUE.)        
C            PIECE = PIECE + QQ1X(S,U,T,S2)       
C         
C         
C  SUMMARY OF ROUTINES        
C  -------------------        
C      THE ROUTINES FOR VARIOUS CROSS-SECTIONS ARE AS FOLLOWS         
C         
C         
C               NORMAL        [1/S2]    [LOG(S2)/S2]   DELTA(S2)      
C      -----------------------------------------------------------    
C         
C      Q-QBAR:  [Q+QBAR-->G+G+GAMMA] + [VIRTUAL] + (F1+F2)*(F1+F2)    
C         
C           QQ1X + (U<-->T)   QQS2A       QQWS2A       QQDELT        
C         
C      Q-QBAR:  [EXTRA VIRTUAL GRAPH FOR Z PRODUCTION TWIXT THRESHOLDS]
C
C                ---           ---         ---         QQDELZ
C
C      Q-G   :  [Q+G-->Q+G+GAMMA] + [VIRTUAL]     
C         
C           QG1               QGS2A       QGWS2A       QGDELT        
C         
C      Q-QBAR:  [EXTRA VIRTUAL GRAPH FOR Z PRODUCTION TWIXT THRESHOLDS]
C
C                ---           ---         ---         QGDELZ
C
C      G-G   :  [G+G-->Q+QBAR+GAMMA]    
C         
C           GG1X + (U<-->T)    ---         ---          ---  
C         
C  WHERE THE [1/S2] AND [LOG(S2)/S2] ROUTINES RETURN THE EXPRESSIONS  
C  MULTIPLYING THE 1/S2 AND LOG(S2)/S2 'S THAT ARE TO BE EVALUATED WITH THE     
C  A+ PERSCRIPTION.  TAKE U<-->T IN Q-G TO GET G-Q.         
C      THE REST OF THE CROSS-SECTIONS JUST HAVE STRAIGHTFORWARD PIECES
C         
C      2(F1+F2)*(F5+F6)                       F1256         
C      2(F1+F2)*(F7+F8)                       F1256 USING U<-->T      
C       (F3+F4)*(F3+F4)                       F3434X + (U<-->T)       
C      2(F3+F4)*(F5+F6)                       F3456         
C      2(F3+F4)*(F7+F8)                       F3456 USING U<-->T      
C       (F5+F6)*(F5+F6)                       F5656         
C       (F7+F8)*(F7+F8)                       F5656 USING U<-->T      
C      2(F5+F6)*(F7+F8) VECTOR PHOTONS        F5678X + (U<-->T)       
C      2(F5+F6)*(F7+F8)  AXIAL PHOTONS        F5678A + (U<-->T)
C      2(F1+F2)*(F3+F4)  AXIAL PHOTONS        F1234A + (U<-->T)     
C      2(H1+H2)*(H5+H6)                       H1256         
C      2(H3+H4)*(H7+H8)                       H1256 WITH U<-->T       
C      2(H1+H2)*(H7+H8)                       H1278         
C      2(H3+H4)*(H5+H6)                       H1278 WITH U<-->T       
C         
C  AND THE FIRST-ORDER PIECES ARE GIVEN BY        
C         
C      Q+QBAR-->G+GAMMA                       QQLEAD(S,T,U) 
C      Q+G-->Q+GAMMA                          QGLEAD(S,T,U) 
C      G+Q-->Q+GAMMA                          QGLEAD(S,U,T) 
C         
C         
C  WHEN CALLING ANY OF THESE ROUTINES, THE GLOBAL VARIABLES QQ AND QMM MUST     
C  HAVE BEEN INITIALIZED.  (SEE COMMON.H)         
C         
C  OTHER MODULES REFERENCED: COMMON.H   
C  LIBRARIES USED: CERN (DDILOG)        
C         
C  OVERALL NORMALIZATIONS     
C  ----------------------     
C      THESE ROUTINES DO *NOT* INCLUDE THE FOLLOWING OVERALL FACTORS --         
C  THEY ARE THE RESPONSIBILITY OF THE CALLER      
C         
C      1ST-ORDER            K(QG)/S * ALPHA(S)    
C      2ND-ORDER            K(QG)/S * ALPHA(S)**2/2/PI      
C         
C  WHERE K(QG) = G(GAMMA)**2/4/NC AND G(GAMMA) IS THE APPROPRIATE COUPLING      
C  CONSTANT.        
C         
C      THE ROUTINES DO NOT INCLUDE FACTORS FOR MULTIPLICITY OF FINAL  
C  PARTICLE FLAVORS.  THIS MUST BE INCLUDED EXPLICITLY BY THE CALLER. 
C  NEITHER HAVE ANY POSSIBLE FACTORS OF 1/2 FOR IDENITICAL FINAL      
C  PARTICLES BEEN INCLUDED.   
C         
C  VERIFICATION     
C  ------------     
C    THE FORMULA FOR ALL 2ND-ORDER PIECES EXCEPT THE VIRTUAL GRAPHS HAVE        
C  BEEN TAKEN DIRECTLY FROM SCHOONSCHIP-CHECKED INPUT.  RATIONAL FRACTIONS      
C  SUCH AS 3/4 HAVE THEN BEEN CONVERTED TO 3./4. BY HAND.  THE VIRTUAL
C  AND FIRST-ORDER PIECES WERE ENTERED COMPLETELY BY HAND.  
C    AT THE MOMENT, THE G-G PROCESS HAS NOT BEEN SIMPLIFIED AND IS TAKEN        
C  DIRECTLY FROM SCHOONSCHIP "PUNCHED" OUTPUT.    
C         
C  REFERENCES       
C  ----------       
C  FOR LABELLING CONVENTION OF GRAPHS, SEE        
C      R. ELLIS, ET. AL., NUCL. PHYS. B211 (1981) 106.      
C         
C **************************************************************************    
          
          
          
C **  INIT -- INITIALIZE LOGARITHMS AND OTHER EXPRESSIONS FOR PARTICULAR        
C **  (S,T,U,S2) FOR USE BY PARTON CROSS-SECTION FORMULAS.  INIT MUST BE        
C **  CALLED *EXPLICITLY* BY THE USER.  
          
      SUBROUTINE INIT(S,T,U,S2)         
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)         
        LOGICAL SWFLAG,SW     
        SAVE SWFLAG 
      COMMON /STU/ SS,TT,UU,QQ,QM,QMM,qmu
      COMMON /XINIT/ S2T,S2U,STS2,SUS2,SQS2,UTS2Q,LAM       
     1  ,LOG4,LOG5,LOG6,LOG8,LOG9,LTAU,LTS2T,LUS2U,WS2      
        DOUBLE PRECISION LAM,LOG4,LOG5,LOG6,LOG8,LOG9,LTAU,LTS2T,LUS2U
      S2T = S2-T    
      S2U = S2-U    
      STS2 = QQ-U   
      SUS2 = QQ-T   
      SQS2 = S+QQ-S2
      UTS2Q = U*T-S2*QQ       
      LAM = SQRT((U+T)**2-4*S2*QQ)      
      LTAU = LOG(QQ/S)        
      LTS2T = LOG(-S2T/T)     
      LUS2U = LOG(-S2U/U)     
      LOG4 = LOG((SQS2+LAM)/(SQS2-LAM)) 
      LOG5 = LOG(S*QMM/S2T/S2U)         
C     $$ THE FACTOR OF QM**2 COMES FROM OUR CONVENTION FOR WS2.       
C     !!!! IS IT RIGHT?       
      LOG6 = LOG(UTS2Q/S2T/S2U)         
      LOG8 = LOG((S2*(2*QQ-U)-QQ*T)**2/S2T**2/S/QQ)         
      LOG9 = LOG((S2*(2*QQ-T)-QQ*U)**2/S2U**2/S/QQ)         
      IF (S2.NE.0) WS2 = LOG(S2/QMM)    
      SWFLAG = .FALSE.        
      RETURN        
          
          
C **  SWITCH -- SWITCH, AS NECESSARY, U<-->T IN ALL EXPRESSIONS STORED
C **  BY INIT.  IF SW = .TRUE., SWITCH WILL ENSURE THAT U<-->T ARE    
C **  SWITCHED WITH RESPECT TO THE DEFINITIONS THAT APPEAR IN INIT.   
C **  IF SW = .FALSE., IT WILL ENSURE THAT THEY ARE NOT.  THE STATUS  
C **  OF WHETHER OR NOT THEY ARE CURRENTLY SWITCHED IS KEPT IN SWFLAG.
          
      ENTRY SWITCH(SW)        
      IF (SW.NEQV.SWFLAG) THEN
        TEMP = S2T  
        S2T = S2U   
        S2U = TEMP  
        TEMP = STS2 
        STS2 = SUS2 
        SUS2 = TEMP 
        TEMP = LTS2T
        LTS2T = LUS2U         
        LUS2U = TEMP
        TEMP = LOG8 
        LOG8 = LOG9 
        LOG9 = TEMP 
        SWFLAG = SW 
      END IF        
      RETURN        
      END 
          
          
C **  INIDS2, SWTDS2 -- SIMILAR TO INIT AND SWITCH BUT DESIGNED FOR   
C **  THE DELTA(S2) PIECES OF SECOND-ORDER CONTRIBUTIONS.   
          
      SUBROUTINE INIDS2(S,T,U,A)        
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)         
        LOGICAL SWFLAG,SW     
        SAVE SWFLAG 
      COMMON /STU/ SS,TT,UU,QQ,QM,QMM,qmu
      COMMON /XINDS2/ SPT,SPU,TPU,WS,WT,WU,WQ,WA  
      TPU = T+U     
      SPU = S+U     
      SPT = S+T     
      WS = LOG(S/QQ)
      WT = LOG(ABS(T/QQ))     
      WU = LOG(ABS(U/QQ))     
      WA = LOG(A/QMM)         
      WQ = LOG(QQ/QMM)        
      SWFLAG = .FALSE.        
      RETURN        
          
      ENTRY SWTDS2(SW)        
      IF (SW.NEQV.SWFLAG) THEN
        TEMP = SPU  
        SPU = SPT   
        SPT = TEMP  
        TEMP = WT   
        WT = WU     
        WU = TEMP   
        SWFLAG = SW 
      END IF        
      RETURN        
      END 
          
          
          
C ** ROUTINES FOR SECOND-ORDER PIECES OTHER THAN DELTA(S2)  
          
      FUNCTION F1256(S,T,U,S2)
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)         
        DOUBLE PRECISION LOGAT,LOGAU,LOGB,LOGC,LOG4M        
      COMMON /STU/ SS,TT,UU,QQ,QM,QMM,qmu
      COMMON /CLEBSH/ CF,CNP,TR,NC,NU,ND,NUD     
      DATA PI /3.1415 92653 58979 32384 62643/

      COMMON /XINIT/ S2T,S2U,STS2,SUS2,SQS2,UTS2Q,LAM       
     1  ,LOG4,LOG5,LOG6,LOG8,LOG9,LTAU,LTS2T,LUS2U,WS2      
        DOUBLE PRECISION LAM,LOG4,LOG5,LOG6,LOG8,LOG9,LTAU,LTS2T,LUS2U
          
          
      T0(QQ,U,T) = U/T + T/U + 2.*QQ*(QQ-U-T)/U/T 
          
          
      F1256 = (2.*CF)  *  CNP*(1./U/T)*(
     1  2.*U - (S**2+3.*U*QQ)/T + S2**2*S**2/T/S2T**2       
     2  + 1./S2*((S+U)**2+S**2*T/S2T)   
     3  + 2./LAM**2*(U**2+U*T+S*T)*(2.*QQ-T+U*(S-QQ)/S2)    
     4  + (U+T)/SQS2 * ( 2.*S/LAM**2*(1+12.*S*QQ/LAM**2)*(U*T-S2*QQ)  
     5    + 1./S2*((S2-U)**2 - 4.*S**2*U**2/LAM**2) )       
     6  + (2.*U*QQ/T - (U**2+2*S*U+2*S**2)/S2) * (LTAU+2.*LTS2T)      
     7  + LOG4/LAM/S2 * ( U**2*(S2+QQ) - 2.*U*S2*QQ         
     8    + S**2/4.*(5.*U+17.*T) + 3.*S*T*U       
     9    + S/LAM**2*(U+T)*(U*(U**2-T**2)-S*(3.*T**2-2.*S2*QQ))       
     1    + 3.*S**2/4./LAM**4*(U+T)*(U**2-T**2)**2 )        
     2  ) 
      RETURN        
          
      ENTRY F3456(S,T,U,S2)   
      F3456 = (2.*CF) * CNP*( 
     1  1./LAM**4 * 48.*S*S2*QQ*(U*T-S2*QQ)/T/SQS2
     2  + 2./LAM**2 * ( -2./T/SQS2*(QQ*(S2-U)**2+3.*S*(T-2*QQ)**2)    
     3    - 1./S*(2.*S2*U-S2**2-QQ**2) + 2.*T-5.*S-8.*QQ    
     4    - 1./S/T*(-S2*(S**2+QQ**2)-7.*S*QQ*(S+QQ)+S**3+QQ**3) )     
     5  + 2./S/T * ( U+2.*T+2.*S2+QQ - 2.*S2*QQ/T + S2*(QQ-U)/S2T     
     6    + ((S2+T)**2-4.*QQ*(S+T))/SQS2 )        
     7  + LOG4 * (  1./LAM**5 * (-24.*QQ*S*S2/T)*(U*T-S2*QQ)
     8    + 2./LAM**3 * ( (S-QQ)/S*(S+S2-QQ)*(U-S2)         
     9      + U**2/T*(QQ-S) - U*(3.*S+3.*S2-QQ)   
     1      - S2/T*(T**2-2.*S2*T-6.*S*QQ+2.*QQ**2) )        
     2    + 1./S/LAM * ( 2.*(2.*QQ*T+S2*(4.*S-S2-T))/SQS2   
     3    +  (S2-U)**2/T + 3.*(S2-T)**2/T + 2.*(S**2/T+S-2.*T)  ))    
     4  - (LTAU+2.*LTS2T) * ((S2-T)**2+(S+QQ)**2)/S/T/SQS2  
     5  - LOG8 * (S2**2+(S2-T)**2)/S/T/SQS2       
     6  ) 
      RETURN        
          
      ENTRY H1256(S,T,U,S2)   
      H1256 = (2.*CF)  *  CNP*(-2./T)*( 
     1  2.*QQ/T*(LTAU+2*LTS2T)
     2  + ((S2-QQ)**2+S**2)/S2T/SQS2 * (LTAU+2.*LTS2T+LOG8)  )        
      RETURN        
          
      ENTRY H1278(S,T,U,S2)   
      H1278 = (2.*CF)  *  CNP*(2./S/T/U)*(        
     1  2.*S2*QQ*(T/U+U/T) - (T+U)**2   
     2  + (2.*S2*(QQ-S)-(T+U)**2)*(LTAU+LTS2T+LUS2U+WS2+LOG5)  )      
      RETURN        
          
      ENTRY F5656(S,T,U,S2)   
      F5656 = (2.*CF)  *  (1./2.)*(     
     1  ( 2.*(S2-U)/T**2*(2.*S2/T-1)**2 - 4.*S2*U/T**3      
     2    - 2.*S/T**2*(1+2.*S2**2/T**2) 
     3    - 1./S*(2.*S2**2/T**2-2.*S2/T+1.)       
     4      *(2.*(S2-U)**2/T**2-2.*(S2-U)/T+1.)   
     5  ) * (3. + LTAU+2.*LTS2T-WS2)    
     6  + (S*(U**2-T**2)/T/LAM**3 + (2.*U+3.*S)/T/LAM) * LOG4         
     7  + (S/S2T**2 + 2.*(S*T+S2*QQ)/T**2/S2T) * WS2        
     8  + (U-T)*(S+S2-QQ)/LAM**2/T      
     9    + 4.*S2/S/T*(S2/T-1.) + 6.*(S2-U)/S/T*((S2-U)/T-1.) + 3./S  
     1    + (9.*U+7.*S-5.*S2)/T**2 + (U+T)/T/S2T - S/S2T**2 
     2  ) 
      RETURN        
          
      ENTRY F5678X(S,T,U,S2)  
C                                             ! MUST ADD TO U <--> T  
      F5678X        
     1 = (2.*CF)  *  (1./2.)*(
     2  S*(S2**2-U*T)/U/T/S2U/S2T       
     3  - (U+T)/2./U/T - (S+S2-QQ)*(U-T)**2/2./LAM**2/U/T   
     4  + 1./LAM * ( (3.*S+2.*S2)*(U+T)/2./U/T + 2.*(S-QQ)/SQS2       
     5    - S*(U+T)*(U-T)**2/LAM**2/2./U/T ) * LOG4         
     6  - QQ*(2.*(S2-T)**2+S2**2+S**2)/U/T/STS2/SUS2 * LOG6 
     7  + ((S2-QQ)**2+(S2-T)**2+S2**2+S**2)/T/STS2/SQS2 * LOG8        
     8  - (2.*QQ-T)*(T**2+2.*(S2**2+S**2-S2*T))/U/T/SUS2/SQS2         
     9    * (LTAU+2*LTS2T)    
     1  ) 
      RETURN        
          
      ENTRY F5678A(S,T,U,S2)  
C                                       ! MUST ADD TO U <--> T        
      F5678A        
     1 = (2.*CF) * (1./2.)*(  
     2  1./T*(S/U-1) + 2.*S/T/S2U + (U-2.*S2)*(U-T)/T/LAM**2
     3  - (S+S2)/U/T * (LTAU+LTS2T+LUS2U)         
     4  + 1./LAM * (-1. + (3.*S+2.*S2)/T + 2.*(S-U)/SQS2    
     5    - S*(U**2-T**2)/T/LAM**2) * LOG4        
     6  - 1./T/SQS2 * (2*(QQ-T) + U*(U-2.*S2)/STS2)         
     7    * (LOG8-2.*LOG6-LTAU-2.*LUS2U)
     8  + 1./U/T * (QQ-2.*S - (U**2+T**2+4.*S*(U+T))/SQS2)  
     9    * (LOG6+LTAU+LTS2T+LUS2U)     
     1  ) 
      RETURN        

      ENTRY F1234A(S,T,U,S2)  
C                                       ! MUST ADD TO U <--> T        
      F1234A
     1 = (2.*CF) * (
     2  QQ/LAM/T*LOG4*( (3.*S/LAM**4*(U-T)**2+(U-S-T)/LAM**2)*(U+T+2*S)
     3    - 1. )
     4  + 1./LAM**4 * 6.*QQ*(U-T)**2/T/SQS2 * ((S2-QQ)**2-3.*S*(S2+QQ))
     5  - 3.*QQ*(U-T)**2*(U+T+S2)/LAM**2/U/T/SQS2
     6  + 2./LAM**2 * ( U + QQ - S2 + QQ*(S-U)/T )
     7 )
      RETURN        
          
      ENTRY F3434X(S,T,U,S2)  
C                                        ! MUST ADD TO U <--> T       
      F3434X = (2.*CF)  *  (1./2.)*(    
     1  1./2./SQS2*(U/S-1) - 5./4./S    
     2  + 1./LAM**2 * ( U/SQS2*(-2.*S+3./2.*U/S*(T-U)+4.*U-2.*T)      
     3    + U/2./S*(2.*S+2.*S2+T-U) )   
     4  + 3.*U**2*(U-T)/LAM**4 * ( (2.*S-U-T)/SQS2 - (S+S2)/S )       
     5  + LOG4/LAM * ( 1./SQS2*(2.*U**2/S+5./2.*U+3./2.*S)  
     6    + 1./S*(3./4.*S+U-1./2.*S2) ) 
     7  + LOG4/LAM**3 * ( U**2/SQS2*(3.*U-T-U**2/S+T**2/S-2.*S)       
     8    + U/S*(2.*S2*T-U*T-2.*S2**2+4.*U*S2-3.*U**2+2.*S*S2-U*S) )  
     9  + LOG4/LAM**5 * ( 3.*(U**2-T**2)*U**2/SQS2*(S-QQ)   
     1    + 3.*U**2*QQ/S*(U-T)*(U+T-2.*S2) )      
     2  ) 
      RETURN        
          
          
      ENTRY QQ1X(S,T,U,S2)    
C                                        ! MUST ADD U <--> T
      QQ1X = (2*CF) *(        
     1  CF*( S/S2T**2 + 2./S2T - S/U/T  
     2    + (2.*QQ*U+T*S2-2.*QQ*S2)/T/UTS2Q - 2.*QQ*(T-U)/T/U/S2T )   
     3  + NC*( -11./6.*(S+QQ)/U/T + S**2*(3.*S2-4.*T)/2./T/U/S2T**2   
     4    + 2.*S/U/S2T + QQ/3./T**2 )   
     5  + TR*( 2./3.*(S+QQ)/U/T - 2./3.*QQ/T**2)  
     6  + 2.*(CF-1./2.*NC)*(LOG5-LOG6+WS2)*( (QQ-U)**2/U/S2U/S2T      
     7    - (QQ-U)**2/U/T/S2T )         
     8  + 2.*(CF-1./2.*NC)*(LOG5+WS2)*(QQ+S)/U/T  
     9  + CF*LOG6*( 4.*QQ*(QQ-T)**2/U/T/UTS2Q     
     1    + (2.*(QQ+S)-S2)/UTS2Q + (S2-2.*S)/U/T )
     2  + CF*(-WS2)*(         
     3    1./U/T/UTS2Q*( 4.*U*T*(U-QQ) - 4.*QQ*(U-QQ)**2 - U*T*S2 )   
     4    + (QQ-U)/S2T**2 - (2.*QQ-U)/T/S2T - 2./T + QQ/T**2 )        
     5 )  
      RETURN        
          
      ENTRY QQS2A(S,T,U,S2)   
      QQS2A = (2.*CF) * T0(QQ,U,T)*(    
     1  (2./3.*TR-11./6.*NC) + 2.*CF*LOG6         
     2  + (CF-1./2.*NC)*2.*(2.*LOG5-LOG6) )       
      RETURN        
          
      ENTRY QQWS2A(S,T,U,S2)  
      QQWS2A = (2.*CF) * T0(QQ,U,T)*(   
     1  4.*CF + (CF-1./2.*NC)*4. )      
      RETURN        
          
          
      ENTRY QGS2A(S,T,U,S2)   
      LOGAT = LTAU + 2*LTS2T  
      LOGAU = LTAU + 2*LUS2U  
      QGS2A =       
     1  NC*(        
     2    (LOG5-LOG6) * (U**2+(QQ-T)**2)/S/S2T    
     3    + LOG5 * ( -((S2+QQ)**2+(S-QQ)**2)/S/T + 2. )     
     4    + LOGAT * ( 2.*(S2-T)*QQ/T**2 
     5      - ((QQ-T)**2+(2.*QQ-U)**2)/2./S/S2T   
     6     - (2.*S+T+2.*QQ)/2./S )      
     7    - LOGAU * (U**2+(S-QQ)**2)/2/S/T        
     8    - LOG4 * (S+2.*U)*(U+T)/2./T/LAM  )     
     9  + CF*(      
     1    - LOG6 * 2./S/T/UTS2Q*(U**2+(S+U)**2)*(U*T+S2*QQ) 
     2    + LOGAU * (U**2+(U+T)**2)/S/T 
     3    - LOG4/LAM * (U+T)/S/T*(U**2+(U+T)**2)  
     4    + 1./S * ((U-QQ)**2/2./T+QQ**2/T+4.*U+3.*T+QQ)    
     5    - 2.*U*(U+S)/S/S2T + U/S2T    
     6    - 4.*QQ*U*T/S/UTS2Q  )        
      RETURN        
          
      ENTRY QGWS2A(S,T,U,S2)  
      QGWS2A =      
     1  NC*(        
     2    - (U**2+(QQ-T)**2)/S/S2T      
     3    - (7.*(U+T)**2+9.*(U+S)**2+4.*U**2)/2./S/T - U**2/S/S2T     
     4    - 2.*(S2-T)*QQ/T**2 
     5    + ((QQ-T)**2+(2.*QQ-U)**2)/2./S/S2T + (2.*S+T+2.*QQ)/2./S  )
     6  + CF*(      
     7    - 2./S/T/UTS2Q*(U**2+(S+U)**2)*(U*T+S2*QQ)        
     8    - 2.*(2.*U+T)/S  )  
      RETURN        
          
      ENTRY QG1(S,T,U,S2)     
      LOGAT = LTAU + 2*LTS2T  
      LOGAU = LTAU + 2*LUS2U  
      LOGB  = LOG6 + WS2      
      LOGC  = LOG8 + WS2      
      LOG4M = LOG4 - 2*LAM/SQS2         
      QG1NC =       
     1  (S-T)/T/S2U - 6.*S*T/S2T**3 - (3.*S+2.*U-7.*T)/T/S2T
     2    - (9.*S+4.*U-4.*T)/S2T**2     
     3    - 3./S/T**4*(U*T-S2*QQ)*(T**2-4.*S2*QQ) 
     4    - (T-QQ)**2/S/T**2 - U*(2.*U-QQ)/S/T**2 + S2/S/T - 1./S     
     5    + (S2+S)*(U-T)/T/LAM**2       
     6  + LOG5 * 1./S         
     7  + WS2 * (   
     8    4./S/T/S2T*(9./4.*U**2-2.*U*T+3.*S*U-3.*S*T+S2*T+2.*S**2)   
     9    + 2.*S*T/S2T**3 + S2*(4.*U+6.*S)/T/S2T**2         
     1    + 11./2./S + (S2-4.*QQ)/S/T - 10./T + 2.*QQ/T**2 )
     2  + (LOGAT-WS2) * ( 2.*(S2-T)/S/T 
     3    + 6.*QQ/S/T - 2.*(S**2+S2**2+4.*QQ**2)/S/T**2     
     4    + 4.*S2*QQ/S/T**4*(U*T-S2*QQ) )         
     5  - LOGC * ((S2-QQ)**2+(U-2*QQ)**2)/2./S/T/S2T        
     6  + LOG4 * 1./2./T/LAM*( 4.*S+3.*U+2.*S2 + 3.*QQ/S*(S2+2.*S-QQ) 
     7    - U/S*(S2-QQ) - 1./LAM**2*(S2+S)*(U-T)*(S+QQ-S2) )
      QG1CF =       
     1  LOGB * 2./S*( (U**2+(S+T-U)**2)/T/STS2 - 2.*S/T     
     2    + 1./UTS2Q*( (2.*U-T)**2 + S2*(S2-U) + QQ*(4.*U-T) ) )      
     3  + LOGC * 1./S/T* ( 2.*S2 - U - 5.*QQ - (U**2+(S2-U)**2)/STS2 )
     4  + LOGAU * 1./S*( (U**2+(S2-U)**2)/T/STS2  
     5    + (S2-2.*U)/T + S2*QQ/U**2 + (2.*S2-T-QQ)/U - 4. )
     6  - WS2 * ( (U**2+(S2-U)**2)/S/T/STS2 + 2./S2U        
     7    + 2.*S*T/S2T**3 - 2.*(U+T)/S2T**2       
     8    + 1./S2T*( (U**2+3.*(U+S)**2)/S/T + 2.*(S2-2*U)/S - 4. )    
     9    - (U+9.*QQ)/S/T + (T-3.*QQ)/S/U + S2*QQ/S/U**2 + 1./S       
     1    - 1./T + 2./U + QQ/T**2 )     
     2  + LOG4M/LAM**3 * ( -3.*S*QQ/2./T/LAM**2*(U-T)**2*(U+T-2.*S2)  
     3    + S2*T/S*(-U-T+3.*S2) + S2**2/S/T*(U-2.*S2)*(T-2.*S)        
     4    + S*S2/T*(7.*S2-11.*S) + S*U/2./T*(-25.*S2-S+3.*U)
     5    + U*(13.*S2+9.*S)   
     6    - 7.*T*(U+T) + S2*(19.*T-18.*S2+7./2.*S) + 3./2.*S*(S+T) )  
     7  + LOG4/LAM * 2./S/T*( QQ*(U-2.*S2-3.*S) + 5.*QQ**2 + S2*U     
     8    - 7./2.*S*U - 2.*S**2 + S2**2 )         
      QG1CF = QG1CF 
     9  + 1./LAM**2/S/T/SQS2 * ( 3.*S2*U*QQ*(2.*U-3.*S2-QQ) 
     1    + 2.*U**2*QQ**2 - 2.*S2**2*QQ*(S2+3.*QQ)
     2    + S2*T*(2.*S2*U+6.*QQ*U+3.*QQ*S2+QQ**2) )         
     3  - (3.*S2+QQ)*U**2/LAM**2/S/T    
     4  + 1./S * ( 1./SQS2*(S2/T*(7.*S2-8.*U)+5.*U+7.*T-16.*S2)       
     5    + 1./2./T*(16.*S2-11.*U+13.*QQ) - (2.*S2+QQ)/U + S2*QQ/U**2 
     6    - S*QQ/2./T**2 - 12. )        
     7  - 2.*S/T/S2U + 8.*S2*(QQ-U)/S2T**3 - (S+4.*QQ)/S2T**2         
     8  + 3.*U*(U+S)/S/T/S2T + 1./S2T*((3.*S2-2.*U)/S-S/T)  
     9  + 1./UTS2Q*( S2 - T - 4.*U**2/S + 4.*QQ/S/T*(U*(U+S)+T*(T+S)) )         
      QG1 = NC*QG1NC + CF*QG1CF         
      RETURN        
          
          
          
          
      ENTRY GG1X(S,T,U,S2)    
C                                        ! MUST ADD U <--> T
      DENU = S2*(2*QQ-T) - QQ*U        
      GG1X = 1./(2.*CF) *(
     1  ( 2.*CF/S2U/STS2 * ( S2**2*QQ**2/S/U**2 + S2*QQ/U 
     2      - (U**2+S2**2)/SQS2 + U**2/S - S2 )       
     3    - NC*(U**2+(S+QQ)**2)/S/S2U/SQS2 ) * (-WS2+LTAU+2.*LUS2U)
     4  + (2.*CF-NC) * (S2**2+(S-QQ)**2)/S/S2T/S2U * (2.*WS2+LOG5)  
     5  + (2.*CF/SUS2-NC/S) * (S2**2+(T-2.*QQ)**2)/S2U/SQS2 
     6    * (LOG9+WS2)      
     7  + ( 2.*CF/UTS2Q * ( 2.*S*(S2+2*QQ)/SQS2 - 4.*S2*(S2+QQ)/S2U 
     8     - S2 - 5.*QQ
     9      + 4.*S2**3*(QQ-S2)/S2T/S2U/SQS2  
     1      + 2.*(S2**2+(T-2.*QQ)**2)
     2        *(S2*(2*QQ-T)-QQ*U)/SUS2/S2U/SQS2 )    
     3    + NC*(U**2+T**2)/S/S2T/S2U ) * (WS2+LOG6)   
     4  - 8.*CF*S2/S2U**2 * WS2      
     5  + CF/LAM * 4.*(2.*QQ-S2)/SQS2 * LOG4    
     6  + NC/LAM * ( S2*(U-T)**2/LAM**4 * ( 4.*QQ*(QQ-S)/SQS2      
     7      - 5.*S2*(QQ-S2)/2./S + 2.*S - 8.*QQ - 9./2.*S2 ) 
     8    + 1./LAM**2 * ( 3.*(U**2-T**2)**2/8./S/SQS2    
     9      + (S2+QQ)/SQS2*(2.*U*T-3.*S2*QQ) + U*T/2./S*(17.*S2+3.*QQ)         
     1      + S2**2/2./S*(S2-21.*QQ) + 3./2.*U*T 
     2      - S2/4.*(7.*U+7.*T+2.*S) )
     3    + 2.*(U*T-5.*QQ**2)/S/SQS2 - 9.*(S2+2.*QQ)/4./SQS2 
     4      + 3.*(S2+7.*QQ)/8./S - 27./8. ) * LOG4         
     5  + 2.*CF * ( 1./S/U**2*(S2*QQ-(S2**2+U**2)/S2U**2*(U*T-S2*QQ))         
     6     - 2.*S2/UTS2Q*(S*U/S2U**2+2.*QQ/S) ) 
     7  + NC * ( 8.*S2/S2T**2 - 4.*(S2-U)/S/S2T + 1./S     
     8     - 1./S/SQS2**2*(4.*U*T+S2*S+(4.*S-3.*S2)*QQ) )   
     9  + NC/LAM**2/SQS2**2 * ( 
     1     3.*S2*QQ*(U-T)**2*(4.*(S+QQ)/LAM**2+1./S)      
     2     + (3.*S2+4.*QQ)*(U-T)**2 - 4.*S2*QQ*(S+QQ) )
     3 )
      RETURN         
          
          
      END 
          
          
          
C **  ROUTINES FOR DELTA(S2) PIECES OF 2ND-ORDER STUFF      
          
      FUNCTION QQDELT(S,T,U)  
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)         
        EXTERNAL DDILOG       
      COMMON /STU/ SS,TT,UU,QQ,QM,QMM,qmu
      COMMON /CLEBSH/ CF,CNP,TR,NC,NU,ND,NUD     
      COMMON /XINDS2/ SPT,SPU,TPU,WS,WT,WU,WQ,WA  
      DATA PI /3.1415 92653 58979 32384 62643/
         
      !wmu = log(qmu**2/qmm)+log(QQ/qmu**2)
          
      T0(QQ,U,T) = U/T + T/U + 2.*QQ*(QQ-U-T)/U/T 
          
      R1(S,T) = LOG(ABS(S/QQ))*LOG(T/(QQ-S)) + LOG(ABS(S/QQ))**2/2.   
     1  - LOG((QQ-T)/QQ)**2/2. + DDILOG(QQ/S) - DDILOG(QQ/(QQ-T))     
C     !!!!! IS THE \S/QQ\ RIGHT?        
          
      R2(T,U) = LOG((QQ-T)/QQ)**2/2. + LOG((QQ-U)/QQ)**2/2. 
     1  + DDILOG(QQ/(QQ-T)) + DDILOG(QQ/(QQ-U))   
          
      F(S,T,U,SPT,SPU,TPU,WS,WU) = CF * (S/SPU+S/SPT+SPT/U+SPU/T)     
     1  + WT * (CF*(4.*S**2+2.*S*T+4.*S*U+T*U)/SPU**2 + NC*T/SPU)     
     2  + WU * (CF*(4.*S**2+2.*S*U+4.*S*T+T*U)/SPT**2 + NC*U/SPT)     
     3  + (2.*CF-NC) * (2.*WS*(S**2/TPU**2+2.*S/TPU)        
     4    - QQ*(T**2+U**2)/T/U/TPU)     
C     !!!! CAREFUL -- SCREWUP IF WT EVER GETS SWAPPED.      
          
          
          
      QQDS2 = (2.*CF) * T0(QQ,U,T)*(    
     1  (11./6.*NC-2./3.*TR)*(WQ-WA) + (67./18.*NC-10./9.*TR)         
     2  + NC*(WQ-WA)**2 + 2.*CF*(WU+WT-2.*WA+2.*WQ)*(-WQ)   
     3  + (CF-1./2.*NC)*( 1./3.*PI**2 + (2.*WA+WS-WQ-WU-WT-WQ)**2 )   
     4  - 3.*CF*(-WQ) )       
      QQFP = - (2.*CF-NC) * ( (S**2+SPU**2)/U/T*R1(S,T)     
     1      + (S**2+SPT**2)/U/T*R1(S,U) )         
     2    + NC*T0(QQ,U,T)*R2(T,U)       
      QQVIRT = (2.*CF) * (    
     1  T0(QQ,U,T)*(
     2    2./3.*PI**2*CF - 1./6.*NC*PI**2         
     3    - 8.*CF - CF*WS**2  
     4    + 1./2.*NC*( WS**2 - (WU+WT)**2 )       
CCPY     5    + (11./6.*NC - 2./3.*TR)*(-WQ + log(qmu**2/qmm)) )        
     5    + (11./6.*NC - 2./3.*TR)*(log(qmu**2/QQ)) )        
     6  + F(S,T,U,SPT,SPU,TPU,WS,WU)    
     7  + QQFP )    
      QQDELT = QQDS2 + QQVIRT 
C     CANCELLING POLE TERMS FOR QQ      
C      (2*CF)*(1-EP) * T0(QQ,U,T)*DS2   
C        * (1-EP*WQ+1/2*EP**2*WQ**2) * ( (2*CF+NC)/EP**2    
C          + 1/EP*( 3*CF + 2*CF*(-WS) + NC*(WS-WU-WT) + 11/6*NC - 2/3*TR ))     
      RETURN        
          
          
      ENTRY QGDELT(S,T,U)     
      QGDS2 = - T0(QQ,S,T)*(  
     1  PI**2/6.*NC + 7./2.*CF - 2.*NC*(WQ-WA)*(WA+WS-WT-WU)
     2  - 2.*NC*WQ*WT         
     3  + NC/2.*(WS-WT-WU)**2 + (11./6.*NC+3.*CF-2./3.*TR)*WQ         
     4  + CF*(WA**2-WQ**2-3./2.*WA-2.*WU*WQ) )    
      QGFP = - (2.*CF-NC)     
     1    * ((U**2+SPU**2)/S/T*(WU*WT-PI**2+0.5*PI**2-R2(T,U))        
     2      + (U**2+TPU**2)/S/T*R1(S,U))
     3    - NC*T0(QQ,S,T)*R1(S,T)       
C         !!!!! IS THE OVERALL SIGN RIGHT???      
      QGVIRT=       
     1  - T0(QQ,S,T)*(        
     2    - 8.*CF - CF*WU**2 - 1./3.*(CF-NC)*PI**2
     3    + 0.5*NC*(WU**2-WS**2-WT**2)
CCPY     4    + (11./6.*NC-2./3.*TR)*(-WQ + log(qmu**2/qmm)) ) 
     4    + (11./6.*NC-2./3.*TR)*(log(qmu**2/QQ)) ) 
     5  - F(U,T,S,TPU,SPU,SPT,WU,WS)    
     6  - QGFP      
      QGDELT = QGDS2 + QGVIRT 
C     CANCELLING POLE TERMS   
C       - DS2*(1-EP*WQ+0.5*EP**2*WQ**2)*T0(QQ,S,T)*(        
C         (2*CF+NC)/EP**2     
C           + (3*CF-2*CF*WU+11/6*NC+NC*(WU-WS-WT)-2/3*TR)/EP )        
      RETURN


C !!!!!  CHECK ME  !!!!!
      ENTRY QQDELZ(S,T,U)
      QQDELZ = (2.*CF) * 2.*(S+QQ)/(S-QQ)*(QQ/(S-QQ)*WS - 1)
      RETURN

C !!!!!  CHECK ME  !!!!!
      ENTRY QGDELZ(S,T,U)
      QGDELZ = - 2.*(U+QQ)/(U-QQ)*(QQ/(U-QQ)*WU - 1)
      RETURN
          
      END 
          
          
          
C **  QQLEAD, QGLEAD -- FIRST-ORDER CROSS-SECTIONS.         
          
      FUNCTION QQLEAD(S,T,U)  
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)         
      COMMON /STU/ SS,TT,UU,QQ,QM,QMM,qmu
      COMMON /CLEBSH/ CF,CNP,TR,NC,NU,ND,NUD     
          
      QQLEAD = 2.*CF*( U/T + T/U + 2.*QQ*S/U/T )  
      RETURN        
          
      ENTRY QGLEAD(S,T,U)     
      QGLEAD = - ( S/T + T/S + 2.*QQ*U/S/T )      
      RETURN        
          
      END 
          
          
C **  DDILOG -- IF YOU HAVE A TRUE IMPLEMENTATION OF DOUBLE-PRECISION 
C **  DILOGARITHM, GET RID OF THIS ROUTINE WHICH JUST CALLS THE       
C **  SINGLE-PRECISION ONE.   
          
C      FUNCTION DDILOG(X)      
C      IMPLICIT DOUBLE PRECISION (A-H,O-Z)         
C        REAL DILOG  
C      DDILOG = DILOG(SNGL(X)) 
C      RETURN        
C      END 

C======================
CCPY W_FAC.FOR

C **************************************************************************    
C         
C  MODULE: WFAC                                       Version 3.0
C  ------------     
C                   by Peter Arnold
C                   Ref:  P. Arnold and H. Reno, Nucl. Phys. B319 (1989) 37
C
C         
C      THIS MODULE CONTAINS THE ROUTINES FOR SWITCHING BETWEEN THE MS-BAR       
C  FACTORIZATION SCHEME AND THAT OF DIEMOZ, FERRONI, LONGO, AND       
C  MARTINELLI.      
C         
C  SUMMARY OF ROUTINES        
C  -------------------        
C      TO USE DIEMOZ, ET. AL. STRUCTURE FUNCTIONS WITH MS-BAR FACTORIZED        
C  PARTON CROSS-SECTIONS, ADD THE FOLLOWING TO THE MS-BAR CROSS-SECTIONS        
C         
C      Q-QBAR:  [Q+QBAR-->G+G+GAMMA] + [VIRTUAL] + (F1+F2)*(F1+F2)    
C                                                       ADD QQB1FAC   
C         
C      Q-G   :  [Q+G-->Q+G+GAMMA]                       ADD QGFAC     
C         
C      G-Q   :  [G+Q-->G+Q+GAMMA]                       ADD QGFAC(U<-->T)       
C         
C      G-G   :  [G+G-->Q+QBAR+GAMMA]                    ADD GGFAC     
C         
C      (F5+F6)*(F5+F6)                                  ADD F56FAC    
C         
C      (F7+F8)*(F7+F8)                                  ADD F56FAC(U<-->T)      
C         
C      THESE CONVERSION EXPRESSIONS MAY HAVE NORMAL, DELTA(S2), [1/S2], AND     
C  [LOG(S2)/S2] PIECES.  WHICH PIECE IS RETURNED IS DETERMINED BY THE 
C  STATE OF THE GLOBAL FLAG LTERM WHICH MUST BE SET APPROPRIATELY BY THE        
C  CALLER (SEE MODULES COMMON AND CODES).  THE GLOBAL VARIABLES QQ
C  AND QMM MUST ALSO BE INITIALIZED.    
C      IN THE PARAMETERS FOR THESE ROUTINES, A IS THE UPPER LIMIT OF THE        
C  S2 INTEGRATION.  
C         
C  OTHER MODULES REFERENCED: COMMON, CODES    
C         
C  OVERALL NORMALIZATIONS     
C  ----------------------     
C      THESE ROUTINES DO *NOT* INCLUDE THE FOLLOWING OVERALL FACTOR --
C  THEY ARE THE RESPONSIBILITY OF THE CALLER      
C         
C                      K(QG)/S * ALPHA(S)**2/2/PI 
C         
C  WHERE K(QG) = G(GAMMA)**2/4/NC AND G(GAMMA) IS THE APPROPRIATE COUPLING      
C  CONSTANT.        
C         
C  REFERENCES       
C  ----------       
C      M. DIEMOZ, ET. AL., CERN-TH.4751/87
C      R. HERROD AND S. WADA, PHYS. LETT. 96B (1980) 195.        
C      G. ALTARELLI, ET. AL., NUCL. PHYS. B157 (1979) 461.  
C      R. ELLIS, ET. AL., NUCL. PHYS. B211 (1981) 106.
C
C  CHANGES FROM PREVIOUS VERSIONS
C  ------------------------------
C    THE CORRECT TRANSFORMATION FROM MS-BAR TO DIS FACTORIZATION SCHEMES
C    IS NOW IMPLEMENTED.      
C         
C **************************************************************************    
          
          
          
C **  FACQQ, FACGQ, FACQG, FACGG -- FACIJ RETURNS THE FINITE PIECE    
C **  FOR THE FACTORIZATION ASSOCIATED WITH J-->I.  THEY CORRESPOND TO
C **  THE NOTATION F_IJ USED BY ELLIS, ET.AL., F_2 USED BY ALTARELLI, 
C **  ET.AL., AND C^(2)_IJ USED BY DIEMOZ, ET.AL. 
C **      THE DIEMOZ, ET.AL. STRUCTURE FUNCTIONS ARE RELATED TO MS-BAR
C **  ONES BY       
C **      
C **       VI(DIEMOZ,X,Q**2) = INT(X,1) DZ/Z VJ(MS-BAR,X/Z,Q**2)      
C **                     *[ DELTA(Z-1) + ALPHA(S,Q**2)/2/PI * F_IJ(Z) ]         
C **      
C **  THE FORMULA FOR F_QQ AND F_QG HAVE BEEN TAKEN FROM ALTARELLI, ET.AL.      
C **  AND REWRITTEN WITH DELTA(1-Z), [1/(1-Z)], [LOG(1-Z)/(1-Z)] IN TERMS       
C **  OF DELTA(S2), [1/S2], [LOG(S2)/S2].         
          
      FUNCTION FACQQ(U,S2,A)        
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)         
      COMMON /CLEBSH/ CF,CNP,TR,NC,NU,ND,NUD     
      COMMON /STU/ SS,TT,UU,QQ,QM,QMM,qmu
      DATA PI /3.1415 92653 58979 32384 62643/

      COMMON /LCODES/ LTERM,  L2DIM,LS20,LWS20,LDS2     
          
      FAC1(Z) = CF*( -(1+Z**2)/(1.-Z)*LOG(Z) + 3. + 2.*Z )  
      FACWS2(Z,S2) = CF*(S2-U)*(1.+Z**2)
      FACS2A(Z,S2) = CF*(S2-U)*( -LOG((S2-U)/QMM)*(1.+Z**2) - 3./2. ) 
      FACDS2() = CF*(-U)*( LOG(-A/U)**2 - LOG(-A/U)*3./2.   
     1  - (9./2.+PI**2/3.) )  
          
      IF (LTERM.EQ.L2DIM) THEN
        Z = U/(U-S2)
        FACQQ = FAC1(Z)       
     1        + (FACS2A(Z,S2)+FACWS2(Z,S2)*LOG(S2/QMM))/S2  
      ELSE IF (LTERM.EQ.LS20) THEN      
        FACQQ = FACS2A(1.D0,0.D0)
      ELSE IF (LTERM.EQ.LWS20) THEN
        FACQQ = FACWS2(1.D0,0.D0)
      ELSE IF (LTERM.EQ.LDS2) THEN      
        FACQQ = FACDS2()      
      END IF        
      RETURN        
      END 
          
          
      FUNCTION FACGQ(U,S2,A)        
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)         
      FACGQ = -FACQQ(U,S2,A)        
      RETURN        
      END 
          
          
      FUNCTION FACQG(U,S2,A)        
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)         
      COMMON /LCODES/ LTERM,  L2DIM,LS20,LWS20,LDS2     
      IF (LTERM.EQ.L2DIM) THEN
        Z = U/(U-S2)
C Old conversion ...
C       FACQG = 1./2.*( (Z**2+(1-Z)**2)*LOG((1-Z)/Z) + 6*Z*(1-Z) )      
C Correct conversion ...
        FACQG = 1./2.*( (Z**2+(1-Z)**2)*LOG((1-Z)/Z) - 1 + 8*Z*(1-Z) )
      ELSE
        FACQG = 0   
      END IF        
      RETURN        
      END 
          
          
      FUNCTION FACGG(U,S2,A)        
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)         
      COMMON /CLEBSH/ CF,CNP,TR,NC,NU,ND,NUD     
      FACGG = -4*TR*FACQG(U,S2,A)   
      RETURN        
      END 
          
          
C **  QQB1FC, QGFC, GGFC, F56FC -- THESE ARE THE MAIN ROUTINE TO      
C **  CALUCLATE THE FACTORIZATION CONTRIBUTIONS DESCRIBED EARLIER.    
C **  THEY SHOULD BE IDENTICAL TO THE FACTORIZATION PIECES USED IN    
C **  THE SCHOONSCHIP CALCULATION OF THE MS-BAR CROSS-SECTIONS BUT    
C **  WITH  P(I,J)/EPSILON --> - F(I,J).
          
      FUNCTION QQB1FC(S,T,U,S2,A)   
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)         
      COMMON /STU/ SS,TT,UU,QQ,QM,QMM,qmu
          
      SP1() = S*U/(U-S2)      
      SP2() = S*T/(T-S2)      
      TP1() = -(U*T-S2*QQ)/(S2-U)       
      UP2() = -(U*T-S2*QQ)/(S2-T)       
          
C !!!! RECHECK OVERALL SIGNS BELOW !!!! 
          
      QQB1FC =    1/U*FACQQ(U,S2,A)*QQLEAD(SP1(),TP1(),U)         
     1          + 1/T*FACQQ(T,S2,A)*QQLEAD(SP2(),T,UP2())         
      RETURN        
          
      ENTRY QGFC(S,T,U,S2,A)        
      QGFC   =    1/U*FACQQ(U,S2,A)*QGLEAD(SP1(),TP1(),U)         
     1          + 1/T*FACQG(T,S2,A)*QQLEAD(SP2(),T,UP2())         
     2          + 1/T*FACGG(T,S2,A)*QGLEAD(SP2(),T,UP2())         
      RETURN        
          
      ENTRY GGFC(S,T,U,S2,A)        
      GGFC =  2*( 1/U*FACQG(U,S2,A)*QGLEAD(SP1(),TP1(),U)         
     1          + 1/T*FACQG(T,S2,A)*QGLEAD(SP2(),UP2(),T) )       
      RETURN        
          
      ENTRY F56FC(S,T,U,S2,A)       
      F56FC  =   1/T*FACGQ(T,S2,A)*QGLEAD(SP2(),T,UP2())
      RETURN        
          
      END 


C=====================-=

