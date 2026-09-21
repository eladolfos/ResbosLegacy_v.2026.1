!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%.%%|%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%                                  %%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%              ResBos              %%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%                                  %%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%         Version 10.02.09cpy      %%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%         Version 06.09.08pn       %%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%                                  %%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
C----------------------------------------------------------------------------
CCPY: March 2008: modified to be used with resbos_p
ccpy: Oct. 2006: correct a misiing factor of 2 in wpres_hvv
C --------------------------------------------------------------------------
      SUBROUTINE WPRES_HVV(WT)
C --------------------------------------------------------------------------
      IMPLICIT NONE

      REAL*8 WT
      REAL*8 W,QCDL,X1X2MIN
      INTEGER*4 NJET,N_PRIMARY_PART,NF,I_STR,IBEAM,NPART,NEVENT,NORDER,
     ,          N_SECOND_PART
      INTEGER*4 IW
      REAL*8 PLAB,PCEM,WHAT,X1,X2
      INTEGER*4 IERR_WT
      REAL*8 XMASS,X_DECAYMASS

      COMMON/AA/W,QCDL,X1X2MIN,NJET,N_PRIMARY_PART,NF,I_STR,IBEAM,NPART,
     &  NEVENT,NORDER,N_SECOND_PART
      COMMON/CHARGE/IW
      COMMON/MOM/PLAB(4,12),PCEM(4,12),WHAT
      COMMON/PMF/X1,X2
      COMMON/MYAMP/IERR_WT
      COMMON/XPART_MASS/XMASS(100),X_DECAYMASS(100)

      REAL*8 XMW,XMZ,XMW2,XMZ2,AEM,GAMW,GAMZ,GAMW2,GAMZ2,
     >SWS,CWS,XMBOT,XMTOP,XMB2,XMT2,GWEAK2,GAMT,
     >XLU,XLD,XRU,XRD,QU,QD,GFERMI,XLE,XRE,
     >XMCHARM,XMCHARM2,XMHIGGS,GAMH,GAMH2
      COMMON/STAND0/XMW,XMZ,XMW2,XMZ2,AEM,GAMW,GAMZ,GAMW2,GAMZ2
      COMMON/STAND1/SWS,CWS,XMBOT,XMTOP,XMB2,XMT2,GWEAK2
      COMMON/STAND2/XLU,XLD,XRU,XRD,QU,QD,GFERMI,XLE,XRE
      COMMON/STAND5/XMCHARM,XMCHARM2,XMHIGGS,GAMH,GAMH2
      COMMON/HTOP/GAMT
CCAO!
      REAL*8 GAMW_NLO,GAMZ_NLO,ALFA_EM_MZ,ACOUPL,WCOUPL,
     &ZCOUPL,HCOUPL,SW2_EFF,EW_ALFA
      COMMON/STAND7/ GAMW_NLO,GAMZ_NLO,ALFA_EM_MZ,ACOUPL,WCOUPL,
     >       ZCOUPL,HCOUPL,SW2_EFF,EW_ALFA
CCAO! The readin W-boson mass and width
      Real*8 XMW_IN,XMW_GRID,GAMW_IN
      COMMON/WMASS_IN/XMW_IN,XMW_GRID,GAMW_IN


      REAL*8 Q_V,QT_V,YB_V,THE_STA,PHI_STA
      COMMON/CSFRAME/Q_V,QT_V,YB_V,THE_STA,PHI_STA
      REAL*8 YMW_AMP,YGW_AMP
      COMMON/VMASS/YMW_AMP,YGW_AMP

      REAL*8 CONST,CONST1,CONST2,CONST3,PI,PROP,WTS,WTA

CsB   3DInterp:
      Integer MD1,MD2,MD3, ID1,ID2,ID3,  IY1,IY2,IY3
      Real*8 X1A,X2A,X3A, YAS,YAA, YDAS,YDAA, X1Y,X2Y,X3Y, YAY
      Parameter (MD1 = 200, MD2 = 500, MD3 = 50)
      Dimension X1A(0:MD1), X2A(0:MD2), X3A(0:MD3)
      Dimension YAS(0:MD1,0:MD2,0:MD3), YAA(0:MD1,0:MD2,0:MD3)
      Dimension YDAS(0:MD1,0:MD2,0:MD3), YDAA(0:MD1,0:MD2,0:MD3)
      Dimension X1Y(0:MD1), X2Y(0:MD2), X3Y(0:MD3)
      Dimension YAY(0:MD1,0:MD2,0:MD3,0:6)
      Common / Arrays / X1A,X2A,X3A, YAS,YAA,  X1Y,X2Y,X3Y, YAY
      Common / Dims / iD1,iD2,iD3, iY1,iY2,iY3
      INTEGER N_QT_INTP,N_Y_INTP,N_Q_INTP

      Character*40 Boson
      Common / Boson / Boson
      REAL*8 ALPHAS,ALFAS
      EXTERNAL ALPHAS
CMRENNA+++
c      COMMON/SCOUPLE/ALFAS
CMRENNA---
      CHARACTER*10 PDF_EVL_NAME
      COMMON/PDFFILE/ PDF_EVL_NAME

CCPY USED FOR SEPERATING NLO AND NLO_SIG CONTRIBUTIONS
      REAL*8 QT_V_SEP
      COMMON/NLO/ QT_V_SEP

      Real*8 cA, cV, VLLS, VLLA
      Real*8 SigS,SigA,ErrSgS,ErrSgA
      Real*8 SigY,ErrSgY
      Integer iYPiece, iReadIn, ii, jj
      Common / YPiece / iYPiece, iReadIN
      Dimension SigY (0:4)
      Real*8 fL,fR, fLA,fRA,fLZ,fRZ, gWeak, sWeak, cWeak
      Real*8 S_Q
      Real*8 Wt1, KFac
      Real*8 sH,tH,uH, CJJ1XFJ_1_AA,SIGS_AA, CJJ1XFJ_1_ZZ,SIGS_ZZ,
     &  SigR_ZZ,SigP_ZZ, x,Mz2,Mz4
      External CJJ1XFJ_1_AA, CJJ1XFJ_1_ZZ
      Real*8 XMAT_PHOT
      Real*8 GGAA
      REAL*8 ETOP,BTOP,CTHE,STHE,BCTH

      Real*8 QT_V_Mat
      Common / Matching / QT_V_Mat
      Integer iOneQ
      Common / OneQ / iOneQ
      Integer iYPSw
      Common / YPertSwitch / iYPSw
      Real*8 eEM, PYalEM
      Real*8 FacMatch
      Real*8 FacN
      Integer iProc
      Common / PartonProcess / iProc
      Real*8 pStar, beta, Sign1

      Real*8 ENFS
      External ENFS

      Real*8 KappaNLO,KappaSoft
      Common / KappaHiggs / KappaNLO,KappaSoft

      Character*2 hDecayMode
      Common / hDM / hDecayMode

      Real*8 BrH2AA
      External BrH2AA

      Real*8 WtPrt
      
      INTEGER*4 I_DECAY,I_HADRONIZE
      COMMON/DEC_AY/I_DECAY,I_HADRONIZE
      

      Real*8 propagator
      
CCAO  
      INTEGER*4 IFLAG
      COMMON/RET/ IFLAG

      Data Pi / 3.141592653589793 /
      Data (SigY(ii), ii=0,4) / 0.d0,0.d0,0.d0,0.d0,0.d0 /
      Data FacN / 1.d0 /
CCAO!
      REAL*8 YXP,YXETA,YXPT,YXDELR,YXPP
      COMMON/SECCONV/ YXP(4,30),YXETA(30),YXPT(30),YXDELR(30,30),
     >YXPP(30)
      REAL*8 PCM_L1(4),PCM_L2(4),PCM_L3(4),PCM_L4(4),PCM_V1(4),
     >       PCM_V2(4) 
      REAL*8 WT_DEC,A12,B12,A34,B34,CPLUS,CMINU,TERMB
      REAL*8 R2,Sw1PROP,Sw2PROP,XMASS_VECT, XWID_VECT
      INTEGER*4 I
      REAL*8 SC
      EXTERNAL SC
      REAL*8 DSIGMA_LEGACY,SYM_IDENTICAL
CCHEN! 
      LOGICAL EENUNU 
         
      IFLAG=0
C This routine is right only in the Collins-Soper frame.

CsB___Initialization
      SigS = 0.d0
      SigA = 0.d0
      SigY(0) = 0.d0
      SigY(1) = 0.d0
      SigY(2) = 0.d0
      SigY(3) = 0.d0
      SigY(4) = 0.d0

CCPY ORDER OF POLYNOMIAL USED IN INTERPOLATION
      N_QT_INTP=4
      N_Y_INTP=4
      N_Q_INTP=4

      IF(NPART.EQ.14) THEN
CCPY FOR LEADING ORDER
        QT_V=0.D0
        N_QT_INTP=1
      ELSEIF(NPART.EQ.13) THEN
C FOR _NLO_SIG
        QT_V=QT_V_SEP
CsB QT_Sep can be a non-grid point value, so the line below is incorrect.
C        N_QT_INTP=1
      ENDIF

      If (nPart.NE.11) then
CsB___Interpolate the Asymptotic/Delta-sigma/LO piece
      Call ThreeDInterp ( N_QT_INTP,N_Y_INTP,N_Q_INTP,
     >      YDAS,YDAA, QT_V,yB_V,Q_V,SigS,SigA, ErrSgS,ErrSgA )      

CsB_____Interpolate the Y piece
        If (iYpiece.Eq.1) then
C          Print*, ' Calling ThreeYInterp '
cC          If (QT_V.LT.QT_V_Sep) then
cCCsB_________If iYPSw=1 then send the Pert piece back
cC            iYPSw = 1
cC          Else
cCCsB_________If iYPSw=0 then send the Y piece back
            iYPSw = 0
cC          End If
            Call ThreeYInterp ( N_QT_INTP,N_Y_INTP,N_Q_INTP,
     >        YDAS,YDAA, QT_V,yB_V,Q_V, SigY, ErrSgY )
CsB Checking interpolation
C          If (jj.LT.10 .Or. (jj.LT.10010 .and. jj.GT.10000)) then
C            Print*, Q_V,QT_V,yB_V,SigY(0)
C          End If
C          jj = jj + 1
c          If (Q_V.GT.36.8 .AND.  Q_V.LT.38.4) then
c            Print*, 'Q_V,QT_V,yB_V,SigY(0),SigS:'
c            Print*, Q_V,QT_V,yB_V,SigY(0),SigS
c          End if
        End IF
      Else ! nPart.EQ.11
CsB_____Resummed case
        If (iYpiece.Eq.1) then
CsB_______For the matched resummed case the CSS piece is included in the Y grid
C         This way the CSS piece interpolation does not have to be called.
          SigA = 0.d0

CsB_______The part below is based on the following observations:
C         1. The Y piece is small for the qQ contribution in the low QT
C         2. The fake Y approximates well the Y piece in the low QT
C         3. After the matching point the fake and real Pert curves have
C            the same slope.
C         Using the above we can use resbos to match the qQ contribution.
          SigS = 0.d0
CsB_______Interpolate the Y piece
          iYPSw = 0
          Call ThreeYInterp ( N_QT_INTP,N_Y_INTP,N_Q_INTP,
     >        YDAS,YDAA, QT_V,yB_V,Q_V, SigY, ErrSgY )
        Else
CsB_______Interpolate the CSS piece
          Call ThreeDInterp ( N_QT_INTP,N_Y_INTP,N_Q_INTP,
     >      YDAS,YDAA, QT_V,yB_V,Q_V,SigS,SigA, ErrSgS,ErrSgA )
        End If
cchen! tmp 
CCHEN!      print*,'QT_V, Q_V: ', QT_V, Q_V
CCHEN!      print*,' SigS, SigY(0): ', SigS, SigY(0)
CCHEN!      print*
        
      End If ! nPart.=!=.11

CCAO  

      IF(IFLAG.EQ.-997)THEN
        WT=0.D0
        RETURN
      ENDIF

      R2=SQRT(2.0D0)

CsB___Weak coupling
      gWeak = DSqrt(GWEAK2)
      sWeak = DSqrt(SWS)
      cWeak = DSqrt(CWS)

CsB___Running alpha_ew implemented on 11/26/97
      eEM   = Sqrt(4.d0*Pi*PYALEM(Q_V**2))

C_____Vector boson - lepton - lepton couplings: fL and fR
C     and mass and width of the vector boson: YMW_AMP, YGW_AMP
      YMW_AMP = xMHiggs
      YGW_AMP = GAMH

CsB_____The general formula for fixed order scalar production and decay:
C       sigma = 4 Pi Gamma(production)*Scalar Propagator*Gamma(decay)

CsB_____The following is implemented in RES.FOR:
C        ELSEIF(TYPE_V.EQ.'H0') THEN
CC         THIS IS THE (\sigma_0)*(Q_V)^2 IN MY NOTES
CC         In SUBROUTINE STANDARD: HCOUPL=R2*GMU*PI/576.0D0 =>
CC         Const = Sqrt(2) Gmu Pi/576 * alpha_s^2/pi^2 Q^2 =
CC                 Sqrt(2) Gmu alpha_s^2 Q^2 / (576 Pi)
C          CONST=HCOUPL*(ALPI(C2*S_Q)*Q_V)**2
C
C_______The resummation formula from PLB283('92)395 (2.7)
C       multiplied by the Higgs propagator, |M(H -> Z Z)|^2,
C       2 Q^2 (a  factor which is included in the phases space for the WBP
C       processes), and the factor which takes the c,b,t masses into account.
C       The total rate for the lowest order process is perfect agreement 
C       (better than .5%) with Papageno.

CsB_____To obtain the real Y piece from the fake one we use the K-factor
C       (only to PLOT the matched QT dist'n at the LHC)
      KFac = 1.d0

CCHEN! CHOICE OF ZZ DECAY MODE
c       EENUNU=.TRUE.
       EENUNU=.FALSE.
      
      If(hDecayMode.Eq.'ZZ')then
        XMASS_VECT=xMZ
        XWID_VECT =GAMZ
        A12=(4.D0*SWS-1)/(4.D0*sWeak*cWeak)
        B12=-1.D0/(4.D0*sWeak*cWeak)
        A34=(4.D0*SWS-1)/(4.D0*sWeak*cWeak)
        B34=-1.D0/(4.D0*sWeak*cWeak)
CCHEN! ADD IN H->ZZ->E E NU NU
        IF(EENUNU) THEN
c         PRINT*,'***********************'
c         PRINT*,'DOING  H->ZZ->E E NU NU'
c         PRINT*,'***********************'
         A34=1.D0/(4.D0*sWeak*cWeak)
         B34=1.D0/(4.D0*sWeak*cWeak)
        ENDIF!END CHEN      
      Elseif(hDecayMode.Eq.'WW')then
        XMASS_VECT=xMW
        XWID_VECT =GAMW
        A12=1.D0/(2.D0*R2*sWeak)
        B12=1.D0/(2.D0*R2*sWeak)
        A34=1.D0/(2.D0*R2*sWeak)
        B34=1.D0/(2.D0*R2*sWeak)
      Else
         print*,'Wrong hDecayMode in subroutine wpres_hvv'
         call exit
      EndIf

      DO I=1,4
        PCM_L1(I)=PCEM(I,3)
        PCM_L2(I)=PCEM(I,4)
        PCM_L3(I)=PCEM(I,5)
        PCM_L4(I)=PCEM(I,6)
        PCM_V1(I)=PCEM(I,3)+PCEM(I,4)
        PCM_V2(I)=PCEM(I,5)+PCEM(I,6)
      ENDDO
  

  
      TERMB=XMASS_VECT*XWID_VECT
      Sw1PROP=(SC(PCM_V1,PCM_V1)-XMASS_VECT**2)**2+TERMB**2
      Sw1PROP=1.0/Sw1PROP
      Sw2PROP=(SC(PCM_V2,PCM_V2)-XMASS_VECT**2)**2+TERMB**2
      Sw2PROP=1.0/Sw2PROP
      
C COLOR FACTOR
      CONST=128.0*(GWEAK2*SWS)**2*R2*GFERMI*XMASS_VECT**4

      CPLUS=(A12**2+B12**2)*(A34**2+B34**2)
     >      +4.0*A12*B12*A34*B34
      CMINU=(A12**2+B12**2)*(A34**2+B34**2)
     >      -4.0*A12*B12*A34*B34
      
cchen! this wt is correct      
      WT_DEC=CONST*Sw1PROP*Sw2PROP*
     >     (CPLUS*SC(PCM_L1,PCM_L3)*SC(PCM_L2,PCM_L4)+
     >      CMINU*SC(PCM_L1,PCM_L4)*SC(PCM_L2,PCM_L3))

CCHEN! FOUND TYPO IN NOTE, this wt_dec is incorrect
CCHEN!      WT_DEC=CONST*Sw1PROP*Sw2PROP*
CCHEN!     >     (CMINU*SC(PCM_L1,PCM_L3)*SC(PCM_L2,PCM_L4)+
CCHEN!     >      CPLUS*SC(PCM_L1,PCM_L4)*SC(PCM_L2,PCM_L3))
       
CCPY SEPT 2006: REWRITE (SigS + KFac*SigY(0)) USED IN RESUMMED
C CALCULATION (NPART.EQ.11)
C_________d(sigma[G G -> H])/d(Q^2)d(y)d(Q_T) from Legacy

c      TYPE(11)=' V^+ --> l l` (RESUMMED) '
c      TYPE(12)=' V^+ --> l l` (NLO, WITH asymp and pert) '
c      TYPE(13)=' V^+ --> l l` (NLO_SIG) '
c      TYPE(14)=' V^+ --> l l` (LO) '
      IF(NPART.EQ.11) THEN 
        DSIGMA_LEGACY= SigS + KFac*SigY(0)
      ELSEIF(NPART.EQ.13) THEN 
        DSIGMA_LEGACY= SigS 
      ELSEIF(NPART.EQ.14) THEN 
        DSIGMA_LEGACY= SigS 
      ELSEIF(NPART.EQ.12) THEN 
        IF(QT_V.LE.QT_V_SEP) THEN
          DSIGMA_LEGACY= KFac*SigY(0)    ! PERT - ASYMP
ctemp      
C      print*,'Q_V,QT_V,YB_V,QT_V_SEP,Sigs,sigy(0) ' 
C      print*,Q_V,QT_V,YB_V,QT_V_SEP,Sigs,sigy(0) 

        ELSE
          DSIGMA_LEGACY= SigS + KFac*SigY(0)         ! PERT PART  
        ENDIF
      ENDIF        

      WT =
C_________Effective coupling of the Higgs to gluons in the heavy top quark limit
C         see Eq.(13) hep-ph/9611272 by Kramer-Laenen-Spira
!     *         (1.0d0+KappaSoft-KappaNLO)  *
C_________Factor which takes into account the quark masses in the g g -> H 
C         fermion loop
     *      ENFS(Q_V) *
C_________d(sigma[G G -> H])/d(Q^2)d(y)d(Q_T) from Legacy
     *      (DSIGMA_LEGACY) *
C_________Higgs propagator
     *      1.d0/((Q_V**2-YMW_AMP**2)**2 + (YMW_AMP*YGW_AMP)**2) *
C_________H -> V1 V2 --> 4f amplitude squared and averaged
     *      WT_DEC*
C_________Extra 2*Q^2 factor to balance the one which is in the phase space
     *       2.0*Q_V**2

      If(hDecayMode.Eq.'ZZ')then
CCPY SYMMETRY FACTOR FOR IDENTICAL PARTICLES IN FINAL STATE
        SYM_IDENTICAL=1.0/2.0
        WT=WT*SYM_IDENTICAL
      ENDIF  
      
      RETURN
      END ! WPRES



C --------------------------------------------------------------------------
      SUBROUTINE HVVRESPH(Y,WT)
C --------------------------------------------------------------------------
      IMPLICIT NONE

      REAL*8 W,QCDL,X1X2MIN
      INTEGER NJET,N_PRIMARY_PART,NF,I_STR,IBEAM,NPART,NORDER,
     &  NEVENT,N_SECOND_PART
      REAL*8 HIGGS_MASS,PLAB,PCEM,WHAT,
     &X1,X2,
     &XMW,XMZ,XMW2,XMZ2,AEM,GAMW,GAMZ,GAMW2,GAMZ2,
     &SWS,CWS,XMBOT,XMTOP,XMB2,XMT2,GWEAK2,
     &XMCHARM,XMCHARM2,XMHIGGS,GAMH,GAMH2,
     &XMASS,X_DECAYMASS,
     &Y_M1,Y_M2,Y_EMIN,
     &GAMT

      COMMON/AA/W,QCDL,X1X2MIN,NJET,N_PRIMARY_PART,NF,I_STR,IBEAM,NPART,
     &  NEVENT,NORDER,N_SECOND_PART
      COMMON/HIGGS/HIGGS_MASS
      COMMON/MOM/PLAB(4,12),PCEM(4,12),WHAT
      COMMON/PMF/X1,X2
      COMMON/STAND0/XMW,XMZ,XMW2,XMZ2,AEM,GAMW,GAMZ,GAMW2,GAMZ2
      COMMON/STAND1/SWS,CWS,XMBOT,XMTOP,XMB2,XMT2,GWEAK2
      COMMON/STAND5/XMCHARM,XMCHARM2,XMHIGGS,GAMH,GAMH2
      COMMON/XPART_MASS/XMASS(100),X_DECAYMASS(100)
      COMMON/MYX12/Y_M1,Y_M2,Y_EMIN
      COMMON/HTOP/GAMT
CCAO!
      REAL*8 GAMW_NLO,GAMZ_NLO,ALFA_EM_MZ,ACOUPL,WCOUPL,
     &ZCOUPL,HCOUPL,SW2_EFF,EW_ALFA
      COMMON/STAND7/ GAMW_NLO,GAMZ_NLO,ALFA_EM_MZ,ACOUPL,WCOUPL,
     >       ZCOUPL,HCOUPL,SW2_EFF,EW_ALFA
CCAO! The readin W-boson mass and width
      Real*8 XMW_IN,XMW_GRID,GAMW_IN
      COMMON/WMASS_IN/XMW_IN,XMW_GRID,GAMW_IN
CCAO! 
      REAL*8 SEC_PLAB
      COMMON/DECAY_MOMENTA/SEC_PLAB(4,50)
      
      REAL*8 YM12MIN,YM12MAX,QTMin,QTMax, TMmin, Tmmax, ETmin
      COMMON/MYYM12/YM12MIN,YM12MAX,QTMin,QTMax, TMmin, Tmmax, ETmin

      REAL*8 Y,WT
      DIMENSION Y(50)
      INTEGER I
      REAL*8 RN
      EXTERNAL RN

      INTEGER*4 IMONTE
      REAL*8 B_QT2,C_QT2,Q2_V,WT_Q2,QT2_V,WT_QT2,WT_YB
      REAL*8 WT_THE,WT_PHI,CTHE,STHE,CPHI,SPHI
      REAL*8 PI,ZTEMP,YM1,YM2,EPS
      REAL*8 ONE,YMW,YGW,YRSH,YMREST,YMW_MIN,YMW_MAX,YMW_VAR,
     >YMW_NORM,YMRES2,YMRES,WT_RES,YBMAX,TAU,YBOOST
      REAL*8 TM2_Q,TM_Q,PHI_Q,QPLUS,QMINUS,R2INV
      REAL*8 P_L1,P_L2,P_Q,P_X,P_Y,P_Z,P_NBAR,P_N,R, P_L1CS
      DIMENSION P_L1(5),P_L2(5),P_Q(5),P_X(5),P_Y(5),
     >P_Z(5),P_NBAR(5),P_N(5),R(2), P_L1CS(5)
      REAL*8 Q_V,QT_V,YB_V,THE_STA,PHI_STA
      COMMON/CSFRAME/Q_V,QT_V,YB_V,THE_STA,PHI_STA

CCPY      REAL*8 QT_V_MAX,YB_V_MAX,Q_V_MIN,Q_V_MAX,QT_V_MXY
CCPY      COMMON/HARDCUT/ QT_V_MAX,YB_V_MAX,Q_V_MIN,Q_V_MAX,QT_V_MXY

      REAL*8 QT_V_MIN, QT_V_MAX,YB_V_MAX,Q_V_MIN,Q_V_MAX,QT_V_MXY
      COMMON/HARDCUT/ QT_V_MIN, QT_V_MAX,YB_V_MAX,
     >  Q_V_MIN,Q_V_MAX,QT_V_MXY

CCPY  USED FOR SEPERATING NLO AND NLO_SIG CONTRIBUTIONS
      REAL*8 QT_V_SEP
      COMMON/NLO/ QT_V_SEP

      REAL*8 T_NBAR,T_N
      DIMENSION T_NBAR(4),T_N(4)
      LOGICAL TESTING
CCPY      REAL*8 QT2_V_MIN,QT2_V_MAX,TEMP0,TEMP1,TEMP,QT_V_MIN
      REAL*8 QT2_V_MIN,QT2_V_MAX,TEMP0,TEMP1,TEMP
      REAL*8 Y_EMAX
      Character*40 Boson
      Common / Boson / Boson
      REAL*8 PCM_L1,ECM_L1,ECM_L2,BETA
      Integer ICalls
      Common / Calls / ICalls
      Integer iYPiece, iReadIn
      Common / YPiece / iYPiece, iReadIn
      Data TESTING /.False./
CCPY      Common / pTMin / QT_V_MIN
      INTEGER*4 I_DECAY,I_HADRONIZE
      COMMON/DEC_AY/I_DECAY,I_HADRONIZE

CCAO!cpn March 1, 2001 
CCAO!      common/cymw/ymw, ygw
      
!      Real*8 QT,MT
      
      Character*2 hDecayMode
      Common / hDM / hDecayMode

      Double Precision ScaPro
      External ScaPro
CCAO! FOR H --> V V
      REAL*8 YMW2_MAX,YMW2_MIN,YMW2_NORM,YMW2_VAR,YM34SQ,YM34
      REAL*8 YMW1_MAX,YMW1_MIN,YMW1_NORM,YMW1_VAR,YM56SQ,YM56
      REAL*8 CP1,CP2,CP3,CP4,CP5,CP6
      DIMENSION CP1(5),CP2(5),CP3(5),CP4(5),CP5(5),CP6(5)
      REAL*8 YP1,YP2,YP3,YP4,YP5,YP6
      DIMENSION YP1(5),YP2(5),YP3(5),YP4(5),YP5(5),YP6(5)
      REAL*8 RT(2),P1(5),P2(5),P3(5),P4(5),P5(5),P6(5)
      REAL*8 WT0, WT1,WT2,WT3,WT4      

C      TESTING=.TRUE.
      If (Testing) Print*, ' in HVVRESPH '

      PI=4.0*ATAN(1.0D0)
      EPS=1.D-4
CsB * At LHC:
      If (W.GT.10.d3) EPS=1.D-3
      WT=-1.0

      Q_V_MIN=YM12MIN
      Q_V_MAX=YM12MAX

CCAO  
      QT_V_MIN=0.D0

CsB For processes requiring grids:
C   QT_V_MAX and YB_V_MAX are determined in ReadIn and passed through COMMON
C   QT_V_MAX is the maximal qT for the resummed grid.

C FIND Q2_V
CCAO! Put in real width of the Higgs Boson in order to generate the 
C     Breit-Wigner shape
      YMW=XMHIGGS
      YGW=GAMH
CsB___Choose a scheme to generate PS
CCAO
      IMONTE=2

CCPY__IN THIS CODE, X1 AND X2 SHOULD BE INTERPRETED AS X_A AND X_B
      IF(IMONTE.EQ.0) THEN ! Use linear shape ------------------------------
        Q2_V = Q_V_Min**2 + (Q_V_MAX**2 - Q_V_MIN**2)*Y(1)
        Q_V=DSQRT(Q2_V)
        WT_Q2=(Q_V_MAX**2 - Q_V_MIN**2)
      ELSEIF(IMONTE.EQ.1) THEN ! USE sigma ~ (1/Q2_V) ----------------------
C_______HARD-WIRED CUT
        If (Q_V_MIN.Eq.0.d0) Q_V_MIN = 10.d-8
        Y_EMIN=Q_V_MIN
        Y_EMAX=Q_V_MAX
        Q2_V=Y_EMIN**2*(Y_EMAX**2/Y_EMIN**2)**Y(1)
        Q_V=DSQRT(Q2_V)
        WT_Q2=Q2_V*DLOG(Y_EMAX**2/Y_EMIN**2)
      ELSE IF(IMONTE.EQ.2) THEN ! BREIT-WIGNER SHAPE ------------------------
        ONE=1.D0
        YRSH=W
        YMREST=0.D0
CCPY
C        YMW_MIN=((YM1+YM2)**2-YMW*YMW)/(YMW*YGW)
C        YMW_MAX=((YRSH-YMREST)**2-YMW*YMW)/(YMW*YGW)
C_______HARD-WIRED CUT
        YMW_MIN=((Q_V_MIN)**2-YMW*YMW)/(YMW*YGW)
        YMW_MAX=((Q_V_MAX)**2-YMW*YMW)/(YMW*YGW)
C
        YMW_VAR=DATAN(YMW_MAX)*Y(1)+(ONE-Y(1))*DATAN(YMW_MIN)
        YMW_NORM=DATAN(YMW_MAX)-DATAN(YMW_MIN)
        YMRES2=YMW*YMW+YMW*YGW*DTAN(YMW_VAR)
        YMRES=DSQRT(YMRES2)
        WT_RES=YMW*YGW*YMW_NORM/(DCOS(YMW_VAR)*DCOS(YMW_VAR))
        Q2_V=YMRES2
        Q_V=YMRES
        WT_Q2=WT_RES

ccpy
c      print*,'Q_V_MIN,Q_V_MAX,YMW,YGW,Y(1)'
c      print*,Q_V_MIN,Q_V_MAX,YMW,YGW,Y(1)
      
                  
      ENDIF
C FIND QT2_V
C ASSUME THE SHAPE AS d(sigma)/d(QT2_V) = c* EXP(-b*QT2_V)
C WHERE b IS FITTED WITH RESUMMED RESULT OF UNPOLARIZED W-BOSON CALCULATIONS
C QT2_V_MAX IS THE MAXIMUM QT2 OF V
CCPY INPUT FROM FIT TO W^+ AT 1.8 TEV, P-PBAR

      
      If (iYPiece.Eq.0) then
        B_QT2=(5.9-1.2)/(QT_V_MAX**2-0.1)
        QT2_V_MAX=MIN(QT_V_MAX**2,QTMAX**2)
      Else If (iYPiece.Eq.1) then
        B_QT2=(5.9-1.2)/(QT_V_MXY**2-0.1)
        QT2_V_MAX=MIN(QT_V_MXY**2,QTMAX**2)
      End If
C HARD-WIRED CUT
      QT2_V_MIN=MAX(QT_V_MIN**2,QTMIN**2)
      TEMP0=EXP(-B_QT2*QT2_V_MIN)
      TEMP=TEMP0-EXP(-B_QT2*QT2_V_MAX)
      C_QT2=B_QT2/TEMP
      QT2_V=-DLOG(TEMP0-TEMP*Y(2))/B_QT2
      QT_V=DSQRT(QT2_V)
      WT_QT2=1.0/C_QT2/(TEMP0-TEMP*Y(2))

C GET YBOOST
C      Print*, ' CALL YMAXIMUM '
      CALL YMAXIMUM(W,Q_V,QT_V,YBMAX)
C      Print*, ' CALLED YMAXIMUM '



C HARD-WIRED CUT
      IF(YBMAX.GT.YB_V_MAX) YBMAX=YB_V_MAX

      WHAT=Q_V
      TAU=(WHAT/W)**2
      YBOOST=YBMAX*(2.0D0*Y(3)-1.0D0)
      YB_V=YBOOST
      WT_YB=2.0D0*YBMAX
      X1=DSQRT(TAU)*DEXP(YBOOST)
      X2=TAU/X1

C ADD THE FACTOR FROM d^4(Q): PI/2
      WT=PI/2.0D0*WT_QT2*WT_Q2*WT_YB

      
CCAO GENERATE 4-BODY PHASE SPACE FOR 4-fermion decay
      IF(hDecayMode.Eq.'ZZ')then
        YMW=XMZ
        YGW=GAMZ
      ElseIf(hDecayMode.Eq.'WW')then
        YMW=XMW
        YGW=GAMW
      Endif
      

      YMW1_MIN=(-YMW*YMW)/(YMW*YGW)
      YMW1_MAX=((WHAT)**2-YMW*YMW)/(YMW*YGW)
      YMW1_VAR=DATAN(YMW1_MAX)*Y(4)+(ONE-Y(4))*DATAN(YMW1_MIN)
      YMW1_NORM=DATAN(YMW1_MAX)-DATAN(YMW1_MIN)
      YM34SQ=YMW*YMW+YMW*YGW*DTAN(YMW1_VAR)
      YM34=DSQRT(YM34SQ)
      WT0=YMW*YGW*YMW1_NORM/(DCOS(YMW1_VAR)*DCOS(YMW1_VAR))

      YMW2_MIN=(-YMW*YMW)/(YMW*YGW)
      YMW2_MAX=((WHAT-YM34)**2-YMW*YMW)/(YMW*YGW)
      YMW2_VAR=DATAN(YMW2_MAX)*Y(5)+(ONE-Y(5))*DATAN(YMW2_MIN)
      YMW2_NORM=DATAN(YMW2_MAX)-DATAN(YMW2_MIN)
      YM56SQ=YMW*YMW+YMW*YGW*DTAN(YMW2_VAR)
      YM56=DSQRT(YM56SQ)
      WT1=YMW*YGW*YMW2_NORM/(DCOS(YMW2_VAR)*DCOS(YMW2_VAR))

      IF(WHAT.LT.(YM34+YM56))THEN
        PRINT*,' WHAT < (YM34+YM56) in HVVRESPH'
        PRINT*,' WHAT,YM34,YM56'
        PRINT*, WHAT,YM34,YM56
        WT=0.d0
        RETURN
      ENDIF

C GENERATE HIGGS DECAY INTO W-BOSON PAIR:  H --> W1 + W2
      DO I=1,2
        RT(I)=Y(I+5)
      ENDDO
      CP1(5)=YM34
      CP2(5)=YM56
      CALL X2PHASE(RT,WHAT,CP1(5),CP2(5),YP1,YP2,WT2)
      DO I=1,4
        CP1(I)=YP1(I)
        CP2(I)=YP2(I)
      ENDDO

C GENERATE THE DECAY OF W1 --> L1 + N1
      DO I=1,2
        RT(I)=Y(I+7)
      ENDDO
      CP3(5)=0.D0
      CP4(5)=0.D0
      CALL X2PHASE(RT,YM34,CP3(5),CP4(5),YP3,YP4,WT3)
      CALL BOOST_5(YP3,CP1,CP3)
      DO I=1,4
        CP4(I)=CP1(I)-CP3(I)
      ENDDO

C GENERATE THE DECAY OF W2 --> L2 + N2
      DO I=1,2
        RT(I)=Y(I+9)
      ENDDO
      CP5(5)=0.D0
      CP6(5)=0.D0
      CALL X2PHASE(RT,YM56,CP5(5),CP6(5),YP5,YP6,WT4)
      CALL BOOST_5(YP5,CP2,CP5)
      DO I=1,4
        CP6(I)=CP2(I)-CP5(I)
      ENDDO


CCAO FILL IN MOMENTUM IN CS FRAME
      DO 410 I=1,2
        PCEM(I,1)=0.D0
        PCEM(I,2)=0.D0
410   CONTINUE
      PCEM(3,1)=Q_V/2.D0
      PCEM(4,1)=PCEM(3,1)
      PCEM(3,2)=-PCEM(3,1)
      PCEM(4,2)=PCEM(3,1)
      DO 414 I=1,4
        PCEM(I,3)=CP3(I)  ! L1
        PCEM(I,4)=CP4(I)  ! L2
        PCEM(I,5)=CP5(I)  ! L3
        PCEM(I,6)=CP6(I)  ! L4
414   CONTINUE

          
      WT=WT*            ! PS of g g --> H
     &   WT0*WT1*WT2*WT3*WT4/(2.0*PI)**2* ! PS of H --> 4f
     &   1.D0/(2*Q_V**2)*
     &   1.D0/(2*PI)*   ! Factor compensate the legacy grid 
     &   1.d0/(Pi*qT_V) ! <= d4Q = pi*qT*dqT*dQ^2*dy

      
      
C In the Collins-Soper frame
C GIVE A RANDOM CHOICE OF AZIMUTHAL ANGLE FOR W-BOSON, IN LAB FRAME
      PHI_Q=2.0*PI*RN(2)
      TM2_Q=Q2_V+QT2_V
      TM_Q=SQRT(TM2_Q)
            
CsB 1&2 (vector boson) 4 momentum in lab frame
      P_Q(1)=QT_V*COS(PHI_Q)
      P_Q(2)=QT_V*SIN(PHI_Q)
      P_Q(3)=TM_Q*(EXP(YB_V)-EXP(-YB_V))/2.0D0
      P_Q(4)=TM_Q*(EXP(YB_V)+EXP(-YB_V))/2.0D0

      P_Q(5)=SQRT(-P_Q(1)**2-P_Q(2)**2-P_Q(3)**2+
     >        P_Q(4)**2)
      
CCAO RETURN TO LAB FRAME
      CALL CS2LAB(CP3,P3,P_Q)
      CALL CS2LAB(CP4,P4,P_Q)     
      CALL CS2LAB(CP5,P5,P_Q)
      CALL CS2LAB(CP6,P6,P_Q)
CCHEN!     
      DO I=1,4           !W-MODE !Z-MODE
        PLAB(I,3)=P3(I)  ! Nu    !E
        PLAB(I,4)=P4(I)  ! E+    !E
        PLAB(I,5)=P5(I)  ! E-    !MU OR NU
        PLAB(I,6)=P6(I)  ! Nu    !MU OR NU

        SEC_PLAB(I,1)=P3(I)  ! Nu
        SEC_PLAB(I,2)=P4(I)  ! E+
        SEC_PLAB(I,3)=P5(I)  ! E-
        SEC_PLAB(I,4)=P6(I)  ! Nu
      ENDDO
      RETURN
      END 


C --------------------------------------------------------------------------
      SUBROUTINE hvvLOPH(Y,WT)
C --------------------------------------------------------------------------
      IMPLICIT NONE

      REAL*8 W,QCDL,X1X2MIN
      INTEGER NJET,N_PRIMARY_PART,NF,I_STR,IBEAM,NPART,NORDER,
     &  NEVENT,N_SECOND_PART
      REAL*8 HIGGS_MASS,PLAB,PCEM,WHAT,
     &X1,X2,
     &XMW,XMZ,XMW2,XMZ2,AEM,GAMW,GAMZ,GAMW2,GAMZ2,
     &SWS,CWS,XMBOT,XMTOP,XMB2,XMT2,GWEAK2,
     &XMCHARM,XMCHARM2,XMHIGGS,GAMH,GAMH2,
     &XMASS,X_DECAYMASS,
     &Y_M1,Y_M2,Y_EMIN,
     &GAMT

      COMMON/AA/W,QCDL,X1X2MIN,NJET,N_PRIMARY_PART,NF,I_STR,IBEAM,NPART,
     &  NEVENT,NORDER,N_SECOND_PART
      COMMON/HIGGS/HIGGS_MASS
      COMMON/MOM/PLAB(4,12),PCEM(4,12),WHAT
      COMMON/PMF/X1,X2
      COMMON/STAND0/XMW,XMZ,XMW2,XMZ2,AEM,GAMW,GAMZ,GAMW2,GAMZ2
      COMMON/STAND1/SWS,CWS,XMBOT,XMTOP,XMB2,XMT2,GWEAK2
      COMMON/STAND5/XMCHARM,XMCHARM2,XMHIGGS,GAMH,GAMH2
      COMMON/XPART_MASS/XMASS(100),X_DECAYMASS(100)
      COMMON/MYX12/Y_M1,Y_M2,Y_EMIN
      COMMON/HTOP/GAMT
CCAO!
      REAL*8 GAMW_NLO,GAMZ_NLO,ALFA_EM_MZ,ACOUPL,WCOUPL,
     &ZCOUPL,HCOUPL,SW2_EFF,EW_ALFA
      COMMON/STAND7/ GAMW_NLO,GAMZ_NLO,ALFA_EM_MZ,ACOUPL,WCOUPL,
     >       ZCOUPL,HCOUPL,SW2_EFF,EW_ALFA
CCAO! The readin W-boson mass and width
      Real*8 XMW_IN,XMW_GRID,GAMW_IN
      COMMON/WMASS_IN/XMW_IN,XMW_GRID,GAMW_IN
CCAO! 
      REAL*8 SEC_PLAB
      COMMON/DECAY_MOMENTA/SEC_PLAB(4,50)

      REAL*8 YM12MIN,YM12MAX,QTMin,QTMax, TMmin, Tmmax, ETmin
      COMMON/MYYM12/YM12MIN,YM12MAX,QTMin,QTMax, TMmin, Tmmax, ETmin

      REAL*8 Y,WT
      DIMENSION Y(50)
      INTEGER I
      REAL*8 RN
      EXTERNAL RN

      INTEGER*4 IMONTE
      REAL*8 B_QT2,C_QT2,Q2_V,WT_Q2,QT2_V,WT_QT2,WT_YB
      REAL*8 WT_THE,WT_PHI,CTHE,STHE,CPHI,SPHI
      REAL*8 PI,ZTEMP,YM1,YM2,EPS
      REAL*8 ONE,YMW,YGW,YRSH,YMREST,YMW_MIN,YMW_MAX,YMW_VAR,
     >YMW_NORM,YMRES2,YMRES,WT_RES,YBMAX,TAU,YBOOST
      REAL*8 TM2_Q,TM_Q,PHI_Q,QPLUS,QMINUS,R2INV
      REAL*8 P_L1,P_L2,P_Q,P_X,P_Y,P_Z,P_NBAR,P_N,R, P_L1CS
      DIMENSION P_L1(5),P_L2(5),P_Q(5),P_X(5),P_Y(5),
     >P_Z(5),P_NBAR(5),P_N(5),R(2), P_L1CS(5)
      REAL*8 Q_V,QT_V,YB_V,THE_STA,PHI_STA
      COMMON/CSFRAME/Q_V,QT_V,YB_V,THE_STA,PHI_STA

CCPY      REAL*8 QT_V_MAX,YB_V_MAX,Q_V_MIN,Q_V_MAX,QT_V_MXY
CCPY      COMMON/HARDCUT/ QT_V_MAX,YB_V_MAX,Q_V_MIN,Q_V_MAX,QT_V_MXY

      REAL*8 QT_V_MIN, QT_V_MAX,YB_V_MAX,Q_V_MIN,Q_V_MAX,QT_V_MXY
      COMMON/HARDCUT/ QT_V_MIN, QT_V_MAX,YB_V_MAX,
     >  Q_V_MIN,Q_V_MAX,QT_V_MXY

CCPY  USED FOR SEPERATING NLO AND NLO_SIG CONTRIBUTIONS
      REAL*8 QT_V_SEP
      COMMON/NLO/ QT_V_SEP

      REAL*8 T_NBAR,T_N
      DIMENSION T_NBAR(4),T_N(4)
      LOGICAL TESTING
CCPY      REAL*8 QT2_V_MIN,QT2_V_MAX,TEMP0,TEMP1,TEMP,QT_V_MIN
      REAL*8 QT2_V_MIN,QT2_V_MAX,TEMP0,TEMP1,TEMP
      REAL*8 Y_EMAX
      Character*40 Boson
      Common / Boson / Boson
      REAL*8 PCM_L1,ECM_L1,ECM_L2,BETA
      Integer ICalls
      Common / Calls / ICalls
      Integer iYPiece, iReadIn
      Common / YPiece / iYPiece, iReadIn
      Data TESTING /.False./
CCPY      Common / pTMin / QT_V_MIN
      INTEGER*4 I_DECAY,I_HADRONIZE
      COMMON/DEC_AY/I_DECAY,I_HADRONIZE

CCAO!cpn March 1, 2001 
CCAO!      common/cymw/ymw, ygw
      
      Character*2 hDecayMode
      Common / hDM / hDecayMode

CCAO! Theorectical Soft-Cutoff
      REAL*8 DELTA_SOFT
      Common /TheoCutoff/ DELTA_SOFT
      Real*8 EPHMIN
      
      
      Double Precision ScaPro
      External ScaPro
CCAO! FOR H --> V V
      REAL*8 YMW2_MAX,YMW2_MIN,YMW2_NORM,YMW2_VAR,YM34SQ,YM34
      REAL*8 YMW1_MAX,YMW1_MIN,YMW1_NORM,YMW1_VAR,YM56SQ,YM56
      REAL*8 CP1,CP2,CP3,CP4,CP5,CP6
      DIMENSION CP1(5),CP2(5),CP3(5),CP4(5),CP5(5),CP6(5)
      REAL*8 YP1,YP2,YP3,YP4,YP5,YP6
      DIMENSION YP1(5),YP2(5),YP3(5),YP4(5),YP5(5),YP6(5)
      REAL*8 RT(2),P1(5),P2(5),P3(5),P4(5),P5(5),P6(5)
      REAL*8 WT0, WT1,WT2,WT3,WT4      

C      TESTING=.TRUE.
C      If (Testing) Print*, ' in WRESPH '

      PI=4.0*ATAN(1.0D0)
      EPS=1.D-4
CsB * At LHC:
      If (W.GT.10.d3) EPS=1.D-3
      WT=-1.0

      Q_V_MIN=YM12MIN
      Q_V_MAX=YM12MAX

      QT_V_MIN=0.D0

C FIND Q2_V
      YMW=XMHIGGS
      YGW=GAMH
CsB___Choose a scheme to generate PS
      IMONTE=2

CCPY__IN THIS CODE, X1 AND X2 SHOULD BE INTERPRETED AS X_A AND X_B
      IF(IMONTE.EQ.0) THEN ! Use linear shape ------------------------------
        Q2_V = Q_V_Min**2 + (Q_V_MAX**2 - Q_V_MIN**2)*Y(1)
        Q_V=DSQRT(Q2_V)
        WT_Q2=(Q_V_MAX**2 - Q_V_MIN**2)
      ELSEIF(IMONTE.EQ.1) THEN ! USE sigma ~ (1/Q2_V) ----------------------
C_______HARD-WIRED CUT
        If (Q_V_MIN.Eq.0.d0) Q_V_MIN = 10.d-8
        Y_EMIN=Q_V_MIN
        Y_EMAX=Q_V_MAX
        Q2_V=Y_EMIN**2*(Y_EMAX**2/Y_EMIN**2)**Y(1)
        Q_V=DSQRT(Q2_V)
        WT_Q2=Q2_V*DLOG(Y_EMAX**2/Y_EMIN**2)
      ELSE IF(IMONTE.EQ.2) THEN ! BREIT-WIGNER SHAPE ------------------------
        ONE=1.D0
        YRSH=W
        YMREST=0.D0
CCPY
C        YMW_MIN=((YM1+YM2)**2-YMW*YMW)/(YMW*YGW)
C        YMW_MAX=((YRSH-YMREST)**2-YMW*YMW)/(YMW*YGW)
C_______HARD-WIRED CUT
        YMW_MIN=((Q_V_MIN)**2-YMW*YMW)/(YMW*YGW)
        YMW_MAX=((Q_V_MAX)**2-YMW*YMW)/(YMW*YGW)
C
        YMW_VAR=DATAN(YMW_MAX)*Y(1)+(ONE-Y(1))*DATAN(YMW_MIN)
        YMW_NORM=DATAN(YMW_MAX)-DATAN(YMW_MIN)
        YMRES2=YMW*YMW+YMW*YGW*DTAN(YMW_VAR)
        YMRES=DSQRT(YMRES2)
        WT_RES=YMW*YGW*YMW_NORM/(DCOS(YMW_VAR)*DCOS(YMW_VAR))
        Q2_V=YMRES2
        Q_V=YMRES
        WT_Q2=WT_RES
      ENDIF

C FIND QT2_V
      QT_V=0.D0
      QT2_V=0.D0

C GET YBOOST
C      Print*, ' CALL YMAXIMUM '
      CALL YMAXIMUM(W,Q_V,QT_V,YBMAX)
C      Print*, ' CALLED YMAXIMUM '
C HARD-WIRED CUT
      IF(YBMAX.GT.YB_V_MAX) YBMAX=YB_V_MAX

      WHAT=Q_V
      TAU=(WHAT/W)**2
      YBOOST=YBMAX*(2.0D0*Y(2)-1.0D0)
      YB_V=YBOOST
      WT_YB=2.0D0*YBMAX
      X1=DSQRT(TAU)*DEXP(YBOOST)
      X2=TAU/X1

C THIS IS FOR d(Q^2)d(Y)d(cos)d(phi)
      WT=WT_Q2*WT_YB

CCAO GENERATE 4-BODY PHASE SPACE FOR 4-fermion decay
      IF(hDecayMode.Eq.'ZZ')then
        YMW=XMZ
        YGW=GAMZ
      ElseIf(hDecayMode.Eq.'WW')then
        YMW=XMW
        YGW=GAMW
      Endif

      YMW1_MIN=(-YMW*YMW)/(YMW*YGW)
      YMW1_MAX=((WHAT)**2-YMW*YMW)/(YMW*YGW)
      YMW1_VAR=DATAN(YMW1_MAX)*Y(3)+(ONE-Y(3))*DATAN(YMW1_MIN)
      YMW1_NORM=DATAN(YMW1_MAX)-DATAN(YMW1_MIN)
      YM34SQ=YMW*YMW+YMW*YGW*DTAN(YMW1_VAR)
      YM34=DSQRT(YM34SQ)
      WT0=YMW*YGW*YMW1_NORM/(DCOS(YMW1_VAR)*DCOS(YMW1_VAR))

      YMW2_MIN=(-YMW*YMW)/(YMW*YGW)
      YMW2_MAX=((WHAT-YM34)**2-YMW*YMW)/(YMW*YGW)
      YMW2_VAR=DATAN(YMW2_MAX)*Y(4)+(ONE-Y(4))*DATAN(YMW2_MIN)
      YMW2_NORM=DATAN(YMW2_MAX)-DATAN(YMW2_MIN)
      YM56SQ=YMW*YMW+YMW*YGW*DTAN(YMW2_VAR)
      YM56=DSQRT(YM56SQ)
      WT1=YMW*YGW*YMW2_NORM/(DCOS(YMW2_VAR)*DCOS(YMW2_VAR))

      IF(WHAT.LT.(YM34+YM56))THEN
        PRINT*,' WHAT < (YM34+YM56) in HVVRESPH'
        PRINT*,' WHAT,YM34,YM56'
        PRINT*, WHAT,YM34,YM56
        WT=0.d0
        RETURN
      ENDIF

C GENERATE HIGGS DECAY INTO W-BOSON PAIR:  H --> W1 + W2
      DO I=1,2
        RT(I)=Y(I+4)
      ENDDO
      CP1(5)=YM34
      CP2(5)=YM56
      CALL X2PHASE(RT,WHAT,CP1(5),CP2(5),YP1,YP2,WT2)
      DO I=1,4
        CP1(I)=YP1(I)
        CP2(I)=YP2(I)
      ENDDO

C GENERATE THE DECAY OF W1 --> L1 + N1
      DO I=1,2
        RT(I)=Y(I+6)
      ENDDO
      CP3(5)=0.D0
      CP4(5)=0.D0
      CALL X2PHASE(RT,YM34,CP3(5),CP4(5),YP3,YP4,WT3)
      CALL BOOST_5(YP3,CP1,CP3)
      DO I=1,4
        CP4(I)=CP1(I)-CP3(I)
      ENDDO

C GENERATE THE DECAY OF W2 --> L2 + N2
      DO I=1,2
        RT(I)=Y(I+8)
      ENDDO
      CP5(5)=0.D0
      CP6(5)=0.D0
      CALL X2PHASE(RT,YM56,CP5(5),CP6(5),YP5,YP6,WT4)
      CALL BOOST_5(YP5,CP2,CP5)
      DO I=1,4
        CP6(I)=CP2(I)-CP5(I)
      ENDDO
      
CCAO FILL IN MOMENTUM IN CS FRAME
      DO 410 I=1,2
        PCEM(I,1)=0.D0
        PCEM(I,2)=0.D0
410   CONTINUE
      PCEM(3,1)=Q_V/2.D0
      PCEM(4,1)=PCEM(3,1)
      PCEM(3,2)=-PCEM(3,1)
      PCEM(4,2)=PCEM(3,1)
      DO 414 I=1,4
        PCEM(I,3)=CP3(I)  ! L1
        PCEM(I,4)=CP4(I)  ! L2
        PCEM(I,5)=CP5(I)  ! L3
        PCEM(I,6)=CP6(I)  ! L4
414   CONTINUE

          
      WT=WT*            ! PS of g g --> H
     &   WT0*WT1*WT2*WT3*WT4/(2.0*PI)**2* ! PS of H --> 4f
     &   1.D0/(2*Q_V**2)*
     &   1.D0/(2*PI)       ! Factor compensate the legacy grid 
C In the Collins-Soper frame
C GIVE A RANDOM CHOICE OF AZIMUTHAL ANGLE FOR W-BOSON, IN LAB FRAME
      PHI_Q=2.0*PI*RN(2)
      TM2_Q=Q2_V+QT2_V
      TM_Q=SQRT(TM2_Q)

CsB 1&2 (vector boson) 4 momentum in lab frame
      P_Q(1)=QT_V*COS(PHI_Q)
      P_Q(2)=QT_V*SIN(PHI_Q)
      P_Q(3)=TM_Q*(EXP(YB_V)-EXP(-YB_V))/2.0D0
      P_Q(4)=TM_Q*(EXP(YB_V)+EXP(-YB_V))/2.0D0

      P_Q(5)=SQRT(-P_Q(1)**2-P_Q(2)**2-P_Q(3)**2+
     >        P_Q(4)**2)
CCAO RETURN TO LAB FRAME
      CALL CS2LAB(CP3,P3,P_Q)
      CALL CS2LAB(CP4,P4,P_Q)     
      CALL CS2LAB(CP5,P5,P_Q)
      CALL CS2LAB(CP6,P6,P_Q)
     
      DO I=1,4
        PLAB(I,3)=P3(I)  ! E+
        PLAB(I,4)=P4(I)  ! Nu
        PLAB(I,5)=P5(I)  ! E-
        PLAB(I,6)=P6(I)  ! Nu
        
        SEC_PLAB(I,1)=P3(I)  ! E+
        SEC_PLAB(I,2)=P4(I)  ! Nu
        SEC_PLAB(I,3)=P5(I)  ! E-
        SEC_PLAB(I,4)=P6(I)  ! Nu
      ENDDO
      
      
      RETURN
      END ! WglorESPH


      FUNCTION SC(P1,P2)
      IMPLICIT NONE
CCAO SCALAR PRODUCTION OF P1 AND P2      
      REAL*8 SC, P1,P2
      DIMENSION P1(4),P2(4)
      REAL*8 TEMP_SC
      
      TEMP_SC=P1(1)*P2(1)+ P1(2)*P2(2)+P1(3)*P2(3)      
      SC=P1(4)*P2(4)-TEMP_SC
      
      RETURN
      END
