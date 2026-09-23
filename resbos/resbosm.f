!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%                                  %%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%             ResBosM              %%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%                                  %%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%         Version 06.19.20cpy      %%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%         Version 11.22.16cpy      %%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%         Version 08.11.16cpy      %%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%         Version 06.16.10cpy      %%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%         Version 06.09.08pn       %%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%                                  %%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
C
C    ResBosM.f  :  contains all the modules needed for the
C                    calculations. (See ResBos.f for details.)
C Authors                       ,
C -------                 C. Balazs & C.-P. Yuan
C                   Department of Physics and Astronomy
C                        Michigan State University
C                   balazs@hep.anl.gov & yuan@pa.msu.edu
C               with substantial modifications by Pavel Nadolsky
C                         nadolsky@pa.msu.edu
C
C     THE EVENT GENERATOR WAS CREATED BY MODIFYING THE POPULAR
C     MONTE CARLO PROGRAM PAPAGENO (VERSION 3.07) WRITTEN BY
C     IAN HINCHLIFFE.
C
C Warning:
C           Watch out the hard-wired cuts build in the subroutine WRESPH,
C           they are optimized for W+ production at the 1.8 TeV Tevatron.
C________________________________________________________________________
CCPY MARCH 2008: ADD H-->ZZ->4 FERMION CHANNEL
CCPY June 2009:  Improve gg -> H calculation with HWIDTH().
CCPY June 2010: refine SW2_EFF in SUBROUTINE DY_COUPL

C------------------------------------------------------------------------
      SUBROUTINE GENERATOR
C------------------------------------------------------------------------
C
C THIS FILE CONTAINS THE MAIN PROGRAMS FOR THE EVENT GENERATOR.
C IT IS MODIFIED FROM PAPAGENO 3.07
C
      IMPLICIT NONE

      REAL*8 W,QCDL,X1X2MIN
      INTEGER*4 NJET,N_PRIMARY_PART,NF,I_STR,IBEAM,NPART,NEVENT,NORDER,
     &IW,N_SECOND_PART
      REAL*8 ETOTMN,PTMIN,ETAMAX,ETASEP
      INTEGER*4 IAPPROX
      INTEGER*4 I_RESONANCE
      REAL*8 FINAL_P
      INTEGER*4 NP_TYPE,N_PARTS
      REAL*8 XL,XU,ACC
      INTEGER*4 NDIM,NCALL,ITMX,NPRN
      REAL*8 XI,SI,SI2,SWGT,SCHI,ALPH
      INTEGER*4 NDO,IT
      LOGICAL IBINNER
      INTEGER*4 IVEGOUT
      LOGICAL I_READIN,I_WRITEOUT
      CHARACTER*50 SAVE_FILE,RESTORE_FILE
      REAL*8 HIGGS_MASS
      INTEGER*4 ITSCALE
      REAL*8 XMW,XMZ,XMW2,XMZ2,AEM,GAMW,GAMZ,GAMW2,GAMZ2,
     &SWS,CWS,XMBOT,XMTOP,XMB2,XMT2,GWEAK2,
     &XMCHARM,XMCHARM2,XMHIGGS,GAMH,GAMH2,
     &XMASS,X_DECAYMASS,
     &Y_M1,Y_M2,Y_EMIN
      INTEGER*4 NCALL1,ITMX1,NCALL2,ITMX2
      REAL*8 A_SUSY,B_SUSY,XMHM,GAMHM,XMH0,GAMH0
      INTEGER*4 NIN,NOUT,I_SEED

      COMMON/AA/W,QCDL,X1X2MIN,NJET,N_PRIMARY_PART,NF,I_STR,IBEAM,NPART,
     &  NEVENT,NORDER,N_SECOND_PART
      COMMON/CHARGE/IW
      COMMON/CUTVAL/ETOTMN,PTMIN,ETAMAX,ETASEP,IAPPROX
      COMMON/RESONANCE/I_RESONANCE
      COMMON/OUT_MOMENTA/FINAL_P(4,50),NP_TYPE(50),N_PARTS
      COMMON/BVEG1/XL(50),XU(50),ACC,NDIM,NCALL,ITMX,NPRN
      COMMON/BVEG2/XI(50,50),SI,SI2,SWGT,SCHI,NDO,IT,ALPH
      COMMON/NININF/IBINNER
      COMMON/DETAILS/IVEGOUT,I_READIN,I_WRITEOUT,SAVE_FILE,RESTORE_FILE
      COMMON/HIGGS/HIGGS_MASS
      COMMON/PAT_SCALE/ITSCALE
      COMMON/STAND0/XMW,XMZ,XMW2,XMZ2,AEM,GAMW,GAMZ,GAMW2,GAMZ2
      COMMON/STAND1/SWS,CWS,XMBOT,XMTOP,XMB2,XMT2,GWEAK2
      COMMON/STAND5/XMCHARM,XMCHARM2,XMHIGGS,GAMH,GAMH2
      COMMON/XPART_MASS/XMASS(100),X_DECAYMASS(100)
      COMMON/MYX12/Y_M1,Y_M2,Y_EMIN
      COMMON/VEGASI/NCALL1,ITMX1,NCALL2,ITMX2
      COMMON/MY_SUSY/A_SUSY,B_SUSY,XMHM,GAMHM,XMH0,GAMH0
      COMMON/EFFIC/NIN,NOUT
      COMMON/SEED/I_SEED
      Character*145 Filename, YFilename
      Common/GridName/ Filename, YFilename
      Character*40 Boson
      Common / Boson / Boson
!      Integer N_SUD_AB,N_SUD_A,N_SUD_B,N_WIL_C,I_SUD
!      Common / Sudakov / N_SUD_AB,N_SUD_A,N_SUD_B,N_WIL_C,I_SUD
      Integer LTOpt
      Common / Process / LTOpt
      Integer iYPiece, iReadIn
      Common / YPiece / iYPiece, iReadIn
      Integer iProc
      Common / PartonProcess / iProc
      Character*2 hDecayMode
      Common / hDM / hDecayMode
      INTEGER*4 I_DECAY,I_HADRONIZE
      COMMON/DEC_AY/I_DECAY,I_HADRONIZE

CCAO! Specify the lepton flavor in W decay
      INTEGER*4 ID_LEPTON, ID_ELE,ID_MU,ID_TAU
      COMMON/LEPTONID/ID_LEPTON,ID_ELE,ID_MU,ID_TAU
CCAO! Input the lepton masses
      REAL*8 XM_ELE,XM_MUON,XM_TAU,XM_NU
      COMMON/LEPTONMASS/ XM_ELE,XM_MUON,XM_TAU,XM_NU

      INTEGER*4 NBINIT, iwgt
      REAL*8 MAX_WGT,NEW_MAX_WGT,WGT_FACTOR
      COMMON/UNWEIGHT/MAX_WGT,NEW_MAX_WGT,WGT_FACTOR,iwgt,NBINIT
CCPY
      INTEGER IFLAV
      REAL*8 CAA_E,CAA_O,CZA_E,CZA_O,CZZ_E,CZZ_O
      COMMON/DY_IFLAV/ CAA_E,CAA_O,CZA_E,CZA_O,CZZ_E,CZZ_O,IFLAV

CCPY      CHARACTER*50 TYPE(50),STRUCFNS(20),SCALETYPE(5)
csb       CHARACTER*50 TYPE(50),STRUCFNS(100),SCALETYPE(5)
      CHARACTER*50 TYPE(50),SCALETYPE(8)
      REAL*8 EMIN,AVGI,SD,CHI2A,EFF
      INTEGER*4 I,J
      REAL*8 FXN
      EXTERNAL FXN

cpn2008 to check the correctness of the reweighting
      real*8 avgi_old
      integer nin_old, nout_old
      real*8 weight_new, weight_old
      integer icut_mismatch
      common/reweight/weight_new, weight_old, icut_mismatch
      data icut_mismatch /0/


CCPY JUMP OVER THE CONVERSION FROM LTOPT TO NPART
      IF(NPART.GE.15 .AND. NPART.LE.20) GOTO 123 
      IF(NPART.EQ.50) GOTO 123  

CsB___NPart -> Switches between LO, NLO and resummed processes.
C             11 does Resummed, 12 does NLO, 14 does LO.
      If (LTOpt.Eq.-1) then
        NPart = 14
CCPY      Else If (LTOpt.Eq.0 .or. LTOpt.Eq.4) then
      Else If (LTOpt.Eq.0 .or. LTOpt.Eq.4 .or. LTOpt.Eq.3) then
        NPart = 11
        If (Boson.Eq.'AA') then ! We'll also calc. built in matrix elements
          XMW = 80.36
          XMZ = 91.187
          Print*, ' Hard wired MW, MZ used for built in amplitudes. '
        EndIf
      Else If (LTOpt.Eq.2) then
        NPart = 12
      Else If (LTOpt.Eq.1) then
        NPart = 13
      Else If (LTOpt.Eq.9) then
CsB_____This is the case when we use the built in matrix elements.
CsB_____Since the EW parameters are not read from the grid:
        XMW = 80.36
        XMZ = 91.187
C        Print*, ' mt, mH = ', XMTop, XMHiggs
        Print*, ' Hard wired MW, MZ are used. '
cpn2005 Npart for built-in matrix elements is set in resbos.for
      Else
        Print *, ' Process is not installed. LTOpt = ', LTOpt
        Stop
      End If

123   CONTINUE

      IF(I_SEED.EQ.0) I_SEED=3219873

CsB___Set up SM parameters, masses, etc.
      CALL STANDARD

 
CsB___Additional setup for G G -> H
      If (Boson.Eq.'H0') Call SetKappa

c      DO I=1,50
c        TYPE(I)=' NO SUCH PROCESS '
c      End Do
C THESE FOUR SETS HAVE THE SAME SETUP
C NOTE: AS OF 10/1995, NPART=11 ALSO DOES W-, Z0 AND A0 DECAYS.
C       THE LABELS OF THE LEPTONS ARE DISTINGUISHED BY EITHER
C       FERMION OR ANTI-FERMION LINE.
c      TYPE(11)=' V^+ --> l l` (RESUMMED) '
c      TYPE(12)=' V^+ --> l l` (NLO) '
c      TYPE(13)=' V^+ --> l l` (NLO_SIG) '
c      TYPE(14)=' V^+ --> l l` (LO) '

CCPY Implement tree level processes
      TYPE(21)=' q qbar --> photon photon g  '
      TYPE(23)=' h1 h2 --> Z Z X  '

      SCALETYPE(1)=' SQRT(S-HAT) FOR PARTON SCALE'
      SCALETYPE(2)=' AVERAGE VALUE OF TRANSVERSE MASS FOR PARTON SCALE'
      SCALETYPE(3)=' TOTAL TRANSVERSE MASS FOR PARTON SCALE'
      SCALETYPE(4)=' VECTOR BOSON MASS XMW OR XMZ '
      SCALETYPE(5)=' INVARIANT MASS OF PARTICLES 3 AND 4'
CPN      SCALETYPE(5)=' Sqrt(sHat) '
      SCALETYPE(6)=' HALF INVARIANT MASS OF PARTICLES 3 AND 4 '
      SCALETYPE(6)=' TWO INVARIANT MASSES OF PARTICLES 3 AND 4 '
      SCALETYPE(8)=' (pT(3)+pT(4))/2 '

C================================
C TYpe_V: W+, W+_RA, W-, W-_RA, Z0, Z0_RA
C
C     NPART = 11   FOR RESUMMED; use WRESPH and WPRES
C     NPart = 12   for "NLO"     (Asymptotic piece AND PERTURBATIVE PIECE);
C                               use WRESPH and WPRES 
C
C     NPart = 13   for "NLO_Sig" (Delta sigma piece); use WLOPH and WPRES
C     NPART = 14   FOR LEADING ORDER INITIAL STATE; 
C                       WITH FINAL STATE NLO QED (S+V) CORRECTIONS; 
C                       use WLOPH and WPRES_RA
C
C     NPART = 15   RESUMMED INITIAL STATE (CSS PIECE ONLY) 
C                  WITH FINAL STATE NLO QED (REAL) EMISSION; 
C                  use WRES3PH and WGRES
C
C     NPART = 16   LEADING ORDER INITIAL STATE 
C                  WITH FINAL STATE NLO QED (REAL) EMISSION;
C                  use WLO3PH and WGRES

C================================
C     NPART = 20    G  G -> photon photon X NLO, WPRESPH-like phase space
C     NPART = 21    q qB -> photon photon X NEXT-TO-LEADING ORDER
C     NPART = 22    G  G -> photon photon X NLO, standard phase space routine
C     NPART = 23    q qB -> Z Z NLO (O(alpha_s), 1 loop)
C     NPART = 24    G  G -> H0 -> Z Z NLO (O(alpha_s^3), 2 loops)

C     NPART = 50    MANIPULATING LEGACY GRIDS
C================================

CsB___Number of primary (and secondary) produced particles
      IF ( NPART.EQ. 1 .OR. NPART.EQ. 2 .OR. NPART.EQ.3 .OR.
     >     NPART.EQ.11 .OR. NPART.EQ.12 .OR.
     >     NPART.EQ.13 .OR. NPART.EQ.14 ) THEN
        N_PRIMARY_PART=2
CCPY!        If (Boson.Eq.'H0' .and. hDecayMode.Eq.'ZZ' .and. i_Decay.Eq.1)
CCPY!     .    N_SECOND_PART = 4
      If (Boson.Eq.'H0' .and. hDecayMode.Eq.'ZZ')
     .    N_SECOND_PART = 4
      If (Boson.Eq.'H0' .and. hDecayMode.Eq.'WW')
     .    N_SECOND_PART = 4
      ELSEIF(NPART.EQ. 4 .OR. NPART.EQ. 5 .OR. NPART.EQ. 6 .OR.
     >       (NPART.GE.21 .AND. nPart.LE.24) ) THEN
        N_PRIMARY_PART=3
CCPY: Nov 2015
      ELSEIF(NPART.EQ.22 .OR. NPART.EQ.23) THEN
        N_PRIMARY_PART=3
      ELSEIF(NPART.EQ.15.OR.NPART.EQ.16) THEN
        N_PRIMARY_PART=3
      ELSEIF(NPART.EQ.50) THEN
        N_PRIMARY_PART=1
      ENDIF
!      Print*,'#PRIMARY,#SECOND=',N_PRIMARY_PART,N_SECOND_PART

C_____Set the number of dimensions for phase-space integration
      IF(NPART.EQ.11 .OR. NPART.EQ.12) THEN ! RESUMMED OR NLO PROCESSES ==== 
        NDIM=4
ctemp        NDIM=3
        If(iYPiece.Eq.1) then
          nDim = 5 ! The Y piece depends on phi...
          If( Boson.EQ.'AA' .or. Boson.EQ.'AG'
     .   .or. Boson.EQ.'ZZ' .or. Boson.EQ.'ZG'
     .   .or. Boson.EQ.'H0')
     .    NDim = 4 ! ... except for AA (AG) or ZZ (ZG) or H0.
        Endif
CCPY
        If((Boson.Eq.'H0').and.
     >     (hDecayMode.Eq.'WW'.Or.hDecayMode.Eq.'ZZ'))then
CCAO_hww
C     3 random #s for Q_V, QT_V and YB_V
C     8 random #s for Vector boson 4-body decay
           NDIM=11
        EndIf

      ELSEIF(NPART.EQ.13 .OR. NPART.EQ.14) THEN ! LO or NLO_Sig processes --
CCPY (NPART.EQ.14) IS FOR LEADING ORDER
C    (NPART.EQ.13) IS FOR _NLO_SIG
        NDIM=3
        If((Boson.Eq.'H0').and.
     >     (hDecayMode.Eq.'WW'.Or.hDecayMode.Eq.'ZZ'))then
CCAO_hww
C     2 random #s for Q_V and YB_V
C     8 random #s for Vector boson 4-body decay
           NDIM=10
        EndIf
CCAO
      ELSEIF(NPART.EQ.15) THEN ! RESUMMED FOR 'W+A(R)' OR 'W-A(R)' PROCESS ====
        NDIM=8
      ELSEIF(NPART.EQ.16) THEN ! LO FOR 'W+A(R)' OR 'W-A(R)' PROCESS --
        NDIM=7
      elseif(npart.eq.21) then !AG->AA+JET with AARESPH
        ndim = 6
      ELSEIF(NPART.EQ.50) THEN ! MANIPULATING LEGACY GRIDS
CCPY THIS ONLY WORKS FOR ONE Q_V VALUE 
        NDIM=2

      ELSE ! Processes with built in amplitudes ----------------------------
        NDIM=3*N_PRIMARY_PART-2
CCPY WE SET PHI_ANGLE IN C.M. FRAME TO BE DUMMY VARIABLE
CsB Not for Phase_Gen
C        NDIM=NDIM-1
      ENDIF ! nPart ========================================================

      NJET=0

      I_RESONANCE=0 ! I_RESONANCE for double EW bosons

      EMIN=0.D0

C_____Masses of produced particles for processes with built in amplitudes
C      Print*, ' nPart = ', nPart
      If (NPART.eq.21.or.NPART.eq.22) THEN 
                                     ! photon photon X production
        XMASS(1)=0.0D0 ! photon 1
        XMASS(2)=0.0D0 ! photon 2
        XMASS(3)=0.0D0 ! parton
      Else If (nPart.Eq.23) then ! q Q -> Z Z X production
        XMASS(1) = xMz ! Z 1
        XMASS(2) = xMz ! Z 2
        XMASS(3) =0.D0 ! parton
      Else If (nPart.Eq.24) then ! h1 h2 -> H X -> Z Z X production
        XMASS(1) = xMz ! Z 1
        XMASS(2) = xMz ! Z 2
        XMASS(3) =0.D0 ! parton
CCAO  
      ELSE IF(NPART.EQ.15.OR.NPART.EQ.16) THEN
        IF(BOSON.EQ.'W+A(R)'.OR.BOSON.EQ.'W-A(R)')THEN
           IF(ID_LEPTON.EQ.ID_ELE)THEN
             XMASS(1) = XM_NU           ! NU
             XMASS(2) = 0.D0            ! PHOTON
             XMASS(3) = XM_ELE          ! POSITRON
           ELSEIF(ID_LEPTON.EQ.ID_MU)THEN
             XMASS(1) = XM_NU           ! NU
             XMASS(2) = 0.D0            ! PHOTON
             XMASS(3) = XM_MUON         ! MUON
           ELSEIF(ID_LEPTON.EQ.ID_TAU)THEN
             XMASS(1) = XM_NU           ! NU
             XMASS(2) = 0.D0            ! PHOTON
             XMASS(3) = XM_TAU          ! TAU
           ENDIF
        ELSEIF(BOSON.EQ.'Z0A(R)')THEN
           IF(ID_LEPTON.EQ.ID_ELE)THEN
             XMASS(1) = XM_ELE          ! ELECTRON
             XMASS(2) = 0.D0            ! PHOTON
             XMASS(3) = XM_ELE          ! POSITION
           ELSEIF(ID_LEPTON.EQ.ID_MU)THEN
             XMASS(1) = XM_MUON          ! MUON^-
             XMASS(2) = 0.D0            ! PHOTON
             XMASS(3) = XM_MUON          ! MUON^+
           ELSEIF(ID_LEPTON.EQ.ID_TAU)THEN
             XMASS(1) = XM_TAU           ! TAU^-
             XMASS(2) = 0.D0            ! PHOTON
             XMASS(3) = XM_TAU          ! TAU^+
           ENDIF 
        ENDIF
      ENDIF

CsB___Tests
!      Call TestDiLog
!      Print*, ' iProc = ', iProc
!      Call TestPYALEM
!      Call TestALPHAS
!      Call TestBrH2AA
!      Call TestPDF

      DO 175 I=1,N_PRIMARY_PART
        EMIN= EMIN+XMASS(I)
175   CONTINUE
      EMIN=DMAX1(EMIN,ETOTMN)

CCPY
      IF(NPART.EQ.1 .OR. NPART.EQ.2 .OR. NPART.EQ.3) THEN
        I_RESONANCE=1
        EMIN=DMAX1(EMIN,ETOTMN)
      ELSEIF(NPART.EQ.4 .OR. NPART.EQ.5 .OR. NPART.EQ.6) THEN
        I_RESONANCE=1
        EMIN=DMAX1(EMIN,ETOTMN)
      ELSEIF(NPART.EQ.11 .OR. NPART.EQ.12 .OR. NPART.EQ.13
     >.OR. NPART.EQ.14) THEN
        I_RESONANCE=1
        EMIN=DMAX1(EMIN,ETOTMN)
CCAO 
      ELSEIF(NPART.EQ.15.OR.NPART.EQ.16) THEN
        I_RESONANCE=1
        EMIN=DMAX1(EMIN,ETOTMN)
CCPY Nov. 2015
C      ELSEIF(NPART.eq.21.or.NPART.eq.22) THEN !gamma gamma production
      ELSEIF(NPART.eq.21.or.NPART.eq.22.or.NPART.eq.23) THEN !gamma gamma production
        I_RESONANCE=1
        EMIN=DMAX1(EMIN,ETOTMN)
      ELSEIF(NPART.EQ.50) THEN
        I_RESONANCE=1
      END IF

      Y_EMIN=EMIN
      N_PARTS=N_PRIMARY_PART

C_____VEGAS PARAMETERS
CBOX
      ACC=-1.D0
      DO 1111 I=1,50
       DO 1113 J=1,50
        XI(I,J)=1.0D0
1113   CONTINUE
       XL(I)=0.0D0
       XU(I)=1.0D0
1111  CONTINUE
C
cpn2008 
      if (iwgt.eq.-1.or.iwgt.eq.-2) 
     >  open(unit=35,file='weights.dat')

      NIN=0
      NOUT=0
      IF(IVEGOUT.EQ.1) THEN
        NPRN=0
      ELSE
        NPRN=-1
      END IF
      IBINNER=.FALSE.
      IF(ITMX1.EQ.0) THEN
       ITMX=1
      ELSE
       ITMX=ITMX1
      END IF
      IF(NCALL1.EQ.0) THEN
       NCALL=100
      ELSE
       NCALL=NCALL1
      END IF
CBOX
      IF(I_READIN) THEN
        CALL RESTR(NDIM,RESTORE_FILE)
      ELSE
        CALL VEGAS(FXN,AVGI,SD,CHI2A)
      END IF
      IBINNER=.TRUE.
      IF(ITMX2.EQ.0) THEN
       ITMX=ITMX
      ELSE
       ITMX=ITMX2
      END IF
      IF(NCALL2.EQ.0)THEN
       NCALL=NEVENT/ITMX
      ELSE
       NCALL=NCALL2
       NEVENT=NCALL*ITMX
      END IF
      NIN=0
      NOUT=0
      CALL VEGAS1(FXN,AVGI,SD,CHI2A)
CBOX
      IF(I_WRITEOUT) THEN
        CALL SAVE(NDIM,SAVE_FILE)
      END IF

      EFF=1.0D0*NIN/(NIN+NOUT)
      PRINT 51
   51 FORMAT('   ')
      PRINT 81
   81 FORMAT(
     &' ------------------------------------------------------------- ')
C      PRINT 81
C      PRINT 51
      PRINT 52
   52 FORMAT ('    EVENT SUMMARY ')
      PRINT 81
      PRINT 51
      IF(IBEAM.EQ.1)PRINT 71, W
      IF(IBEAM.EQ.0)PRINT 170, W
      IF(IBEAM.EQ.-1)PRINT 72, W
      IF(IBEAM.EQ.-2)PRINT 171
   71 FORMAT('    PROTON - PROTON COLLISIONS, SQRT(S)=',F12.3,' GEV')
   72 FORMAT('    PROTON-ANTIPROTON COLLISIONS, SQRT(S)=',F12.3,' GEV')
  170 FORMAT('    Proton-Nucleon collisions, SQRT(S)=',F12.3,' GEV')
  171 FORMAT('    Pion-Nucleon collisions, SQRT(S)=',F12.3,' GEV')
C      PRINT 30,W
C   30 FORMAT(/,'    [SQRT(S)]  :  ',F12.3,'  GEV')
csb      PRINT 73, STRUCFNS(I_STR)
      PRINT 70, Filename
      PRINT 75, YFilename

      IF(ITSCALE.EQ.4)THEN
        Print*, '  '
        PRINT 76
C        PRINT 77,HIGGS_MASS
        PRINT 77, xMHiggs
        Print*, '  '
      ELSE
        PRINT 73, SCALETYPE(ITSCALE)
      END IF
   76 FORMAT('    Renormalization and factorization scales = Q ')
   77 FORMAT('    HIGGS MASS = ', F10.3, ' GEV')
      PRINT 81
      PRINT 51
   70 FORMAT(/,2X,'  Low  Q_T grid: ', A40)
   75 FORMAT(  2X,'  High Q_T grid: ', A40)
   73 FORMAT(/,2X,'  Scale type: ', A60)
      Print*, '    TYPE OF PRODUCED PARTICLES:'
CCPY

CsB  Print process type.
        IF (NJET.EQ.0) then
CCPY
CsB         PRINT 18, TYPE(NPART)
          If (NPART.Eq.11) then ! ------------------------------------------
            If (Boson.Eq.'W+') Print*,
     >                       '       W^+ --> nu_e e^+ (RESUMMED)'
            If (Boson.Eq.'W-') Print*,
     >                       '       W^- --> nu^Bar_e e^- (RESUMMED)'
            If (Boson.Eq.'Z0') Print*,
     >                       '       Z^0 --> e^+ e^- (RESUMMED)'
            If (Boson.Eq.'W+_RA') Print*,
     >                       '       W^+ --> nu_e e^+ (RESUMMED)'
            If (Boson.Eq.'W-_RA') Print*,
     >                       '       W^- --> nu^Bar_e e^- (RESUMMED)'
            If (Boson.Eq.'Z0_RA') Print*,
     >                       '       Z^0 --> e^+ e^- (RESUMMED)'
            If (Boson.Eq.'HP') Print*,
     >                       '       W+* --> W+ h (RESUMMED)'

            If (Boson.Eq.'HM') Print*,
     >                       '       W-* --> W- h (RESUMMED)'

            If (Boson.Eq.'HZ') Print*,
     >                       '       Z0* --> Z0 h (RESUMMED)'

            If (Boson.Eq.'DY'.AND. IFLAV.EQ.1) Print*,
     >           '  U_UB --> Z^0,A^0 --> e^+ e^- (RESUMMED DY Pair)'

            If (Boson.Eq.'DY'.AND. IFLAV.EQ.2) Print*,
     >           '  D_DB --> Z^0,A^0 --> e^+ e^- (RESUMMED DY Pair)'

            If (Boson.Eq.'A0') Print*,
     >                       '       A^0 --> e^+ e^- (RESUMMED)'
            If (Boson.Eq.'AA') Print*,
     >                       '       q qB(G) --> A^0 A^0 (RESUMMED)'
            If (Boson.Eq.'AG') Print*,
     >                       '       G G --> A^0 A^0 (RESUMMED)'
            If (Boson.Eq.'ZZ') Print*,
     >                       '       q qB(G) --> Z^0 Z^0 (RESUMMED)'
            If (Boson.Eq.'ZG') Print*,
     >                       '       G G --> Z^0 Z^0 (RESUMMED)'
CCPY!            If (Boson.Eq.'H0') Print*,
CCPY!     >                       '       G G --> H^0 (RESUMMED) '
            If (Boson.Eq.'H0') then
                Print*,      '       G G --> H^0 (RESUMMED) '
                If (hDecayMode.Eq.'ZZ') Print*,
     >                        '       H --> Z^0 Z^0 --> 4 fermion'
                If (hDecayMode.Eq.'WW') Print*,
     >                        '       H --> W^+ W^- --> 4 fermion'
                If (hDecayMode.Eq.'AA') Print*,
     >                        '       H --> A^0 A^0              '
            End If

            If (Boson.Eq.'GL') Print*,
     >   '  q qbar -> GL --> t tbar (partial RESUMMED, 
     > not exact calculation)'

            If (Boson.Eq.'GG') Print*,
     >   '  g g --> t tbar (partial RESUMMED, not exact calculation)'

            If (Boson.Eq.'W+A(SV)') Print*,
     >                 '    W^+ A --> nu_e e^+  (s+v) (RESUMMED)'
            If (Boson.Eq.'W-A(SV)') Print*,
     >                 '   A W^- --> nu^Bar_e e^- (s+v) (RESUMMED)'
            If (Boson.Eq.'Z0A(SV)') Print*,
     >                 '   A Z0 --> e^+ e^- (s+v) (RESUMMED)'
          Else If (NPART.Eq.12) then ! -------------------------------------
            If (Boson.Eq.'W+') Print*,
     >                       '       W^+ --> nu_e e^+ (NLO)'
            If (Boson.Eq.'W-') Print*,
     >                       '       W^- --> nu^Bar_e e^- (NLO)'
            If (Boson.Eq.'Z0') Print*,
     >                       '       Z^0 --> e^+ e^- (NLO)'
            If (Boson.Eq.'W+_RA') Print*,
     >                       '       W^+ --> nu_e e^+ (NLO)'
            If (Boson.Eq.'W-_RA') Print*,
     >                       '       W^- --> nu^Bar_e e^- (NLO)'
            If (Boson.Eq.'Z0_RA') Print*,
     >                       '       Z^0 --> e^+ e^- (NLO)'
            If (Boson.Eq.'A0') Print*,
     >                       '       A^0 --> e^+ e^- (NLO)'
            If (Boson.Eq.'AA') Print*,
     >                       '       q qB(G) --> photon photon (NLO)'

C            If (Boson.Eq.'GL') Print*,
C     >     '  q qbar -> GL --> t tbar (partial NLO, not exact calculation)'

          Else If (NPART.Eq.13) then ! -------------------------------------
            If (Boson.Eq.'W+') Print*,
     >                       '       W^+ --> nu_e e^+ (Delta Sigma)'
            If (Boson.Eq.'W-') Print*,
     >                       '       W^- --> nu^Bar_e e^- (Delta Sigma)'
            If (Boson.Eq.'Z0') Print*,
     >                       '       Z^0 --> e^+ e^- (Delta Sigma)'
            If (Boson.Eq.'W+_RA') Print*,
     >                       '       W^+ --> nu_e e^+ (Delta Sigma)'
            If (Boson.Eq.'W-_RA') Print*,
     >                       '       W^- --> nu^Bar_e e^- (Delta Sigma)'
            If (Boson.Eq.'Z0_RA') Print*,
     >                       '       Z^0 --> e^+ e^- (Delta Sigma)'
            If (Boson.Eq.'A0') Print*,
     >                       '       A^0 --> e^+ e^- (Delta Sigma)'
C            If (Boson.Eq.'GL') Print*,
C     > 'q qbar -> GL --> t tbar (partial, Delta Sigma, not exact calculation)'

CCAO 
          Else If (NPART.EQ.15) then
            If (Boson.Eq.'W+A(R)') Print*,
     >                       '       W^+ A --> nu_e e^+  A (RESUMMED)'
            If (Boson.Eq.'W-A(R)') Print*,
     >                       '   A W^- --> A  nu^Bar_e e^- (RESUMMED)'
            If (Boson.Eq.'Z0A(R)') Print*,
     >                       '       Z0 A --> e^- e^+  A (RESUMMED)'
CCAO
          Else If (NPART.EQ.16) then
            If (Boson.Eq.'W+A(R)') Print*,
     >                       '       W^+ A --> nu_e e^+  A (LO)'
            If (Boson.Eq.'W-A(R)') Print*,
     >                       '   A W^- --> A  nu^Bar_e e^- (LO)'
            If (Boson.Eq.'Z0A(R)') Print*,
     >                       '       Z0 A --> e^- e^+  A (LO)'
          Else If (NPART.eq.21.or.npart.eq.22) then ! --------------------

            if (Boson.eq.'AA') then
              If(iProc.Eq.0) then
                Print*, '       q qB(G) --> photon photon X (NLO)'
              ElseIf(iProc.Eq.1) then
                Print*, '       q qB --> photon photon X (NLO)'
              ElseIf(iProc.Eq.2) then
                Print*, '       q G --> photon photon X (NLO)'
              End If            !iProc
            elseif (Boson.eq.'AG') then
              Print*, '       G G --> photon photon G (NLO)'
            endif               !Boson
            
          Else If (NPART.Eq.23) then ! -------------------------------------
            If(iProc.Eq.0) then
                               Print*,
     >                       '       q qB(G) -> Z Z X  (NLO)'
            ElseIf(iProc.Eq.1) then
                               Print*,
     >                       '       q qB -> Z Z X  (NLO)'
            ElseIf(iProc.Eq.2) then
                               Print*,
     >                       '       q G --> Z Z X  (NLO)'
            End If
          Else If (NPART.Eq.24) then ! -------------------------------------
            Print*,
     >      '   h1 h2 -> H0 -> Z0 Z0 X   O(alpha alpha_s^3))'
          Else If (NPART.Eq.50) then ! -------------------------------------
            If(iProc.Eq.0) then
                               Print*,
     >                 ' Legacy grid files including q Q and q G '
            ElseIf(iProc.Eq.1) then
                               Print*,
     >                 ' Legacy grid files including only q Q '
            ElseIf(iProc.Eq.2) then
                               Print*,
     >                 ' Legacy grid files including only q G '
            End If
          End If
        End If

   18 FORMAT ('      ',A50)

CCPY
C        IF(NPART.NE.3 .AND. NPART.NE.6) THEN
C          IF(IW.EQ.1) THEN
C            PRINT 65
C   65       FORMAT('   POSITIVE CHARGE: ')
C          ELSEIF(IW.EQ.-1) THEN
C            PRINT 66
C   66       FORMAT('   NEGATIVE CHARGE: ')
C          ELSE
C            PRINT 67
C   67       FORMAT('   CHARGES SUMMED ')
C          END IF
C        END IF

cpn2008 during the reweighting, save the total cross section and numbers
c      of events passing and failing the cuts for an additional cross
c      check during the iwgt=-2 runs
      if (iwgt.eq.-1) then
        write (35,7) avgi, nin, nout
      else if (iwgt.eq.-2) then
        read (35,7) avgi_old, nin_old, nout_old
        if (abs(avgi/avgi_old - 1d0).gt.1e-1.or.nin.ne.nin_old
     >    .or.nout.ne.nout_old) then
          print *,'Error: the Monte-Carlo integration events for ',
     >      'iwgt=-1 and iwgt = -2 do not match'
          print *,
     >      'iwgt, total cross section, accepted events, failed events:'
          print *,-1, avgi_old, nin_old, nout_old
          print *,-2, avgi, nin, nout
          call quit
       endif                   !abs(...
 7     format(g20.12,1x,i10,1x,i10)
      endif                     !iwgt

      Print 51
      PRINT 60, (XMASS(I),I=1,N_PRIMARY_PART)
   60 FORMAT('    FINAL STATE PRIMARY PARTICLE MASSES:',6F10.2)
      PRINT 51
      PRINT 81
      PRINT 10, AVGI,2.*SD
   10 FORMAT(/,'    TOTAL CROSS SECTION =',G16.7,'  PB +/-',G16.7)
      PRINT 20,NIN,NOUT
   20 FORMAT(/,' ',I10,' EVENTS PASSING CUTS',I10,' EVENTS FAILING')
      PRINT 23,EFF
   23 FORMAT('    EFFICIENCY :',G16.8)
      PRINT 40, EMIN
   40 FORMAT(/,'    MINIMUM Sqrt(sHat) = ',1G16.8,'  GeV',/)

CCPY
      WRITE(16,*) 'C*****************************   '
      WRITE(16,10) AVGI,2.*SD
      WRITE(16,73) SCALETYPE(ITSCALE)
      WRITE(16,70)  Filename
      WRITE(16,75) YFilename
      WRITE(16,20) NIN,NOUT
      WRITE(16,23) EFF
      WRITE(16,40) EMIN
      WRITE(16,*) 'C*****************************   '

cpn2008 
      if (iwgt.eq.-1.or.iwgt.eq.-2) then
        print *,'NOTE: reweighting is turned on (un/reweight = ',
     >    iwgt, ' in resbos.in).'
        print *,
     >  'The above total cross sections correspond to the saved event',
     >  'weights, not the weights from the current run'
        write(16,*) 'NOTE: reweighting is turned on (un/reweight = ',
     >    iwgt, ' in resbos.in).'
        write(16,*)
     >  'The above total cross sections correspond to the saved event',
     >  'weights, not the weights from the current run'
        if (icut_mismatch.ne.0) 
     >    print *,'Fixed ',icut_mismatch,' events with mismatched ',
     >    'cuts in the iwgt=-1 and -2 runs (nothing to worry)'
        close(35) 
      endif                     !iwgt

        RETURN
      END ! Generator

C --------------------------------------------------------------------------
      SUBROUTINE MAKE_EVENT(Y,WEIGHT,IPASS)
C --------------------------------------------------------------------------
      IMPLICIT NONE

      REAL*8 W,QCDL,X1X2MIN
      INTEGER*4 NJET,N_PRIMARY_PART,NF,I_STR,IBEAM,NPART,NEVENT,NORDER,
     &ITSCALE,N_SECOND_PART
      REAL*8 XMW,XMZ,XMW2,XMZ2,AEM,GAMW,GAMZ,GAMW2,GAMZ2,
     &SWS,CWS,XMBOT,XMTOP,XMB2,XMT2,GWEAK2,
     &XMCHARM,XMCHARM2,XMHIGGS,GAMH,GAMH2,
     &XMASS,X_DECAYMASS,
     &Y_M1,Y_M2,Y_EMIN,
     &ZUP,ZDOWN,XMGW2,XMGZ2,XMGH2,XMH2,VSQ,
     &GAMT,
     &PLAB,PCEM,WHAT
      INTEGER*4 NP_TYPE,N_PARTS
      REAL*8 FINAL_P,
     &Q1,Q2,
     &PT,ETA,DR,ETOT,
     &Q1Q,Q2Q,Q1W,Q2W,Q_HOLD,
     &X1,X2,
     &SEC_PT,SEC_ETA
      INTEGER*4 I_RESONANCE
      REAL*8 Q
      LOGICAL IBINNER
      REAL*8 XL,XU,ACC
      INTEGER*4 NDIM,NCALL,ITMX,NPRN,NDO,IT
      REAL*8 XI,SI,SI2,SWGT,SCHI,ALPH

      COMMON/AA/W,QCDL,X1X2MIN,NJET,N_PRIMARY_PART,NF,I_STR,IBEAM,NPART,
     &  NEVENT,NORDER,N_SECOND_PART
      COMMON/PAT_SCALE/ITSCALE
      COMMON/STAND0/XMW,XMZ,XMW2,XMZ2,AEM,GAMW,GAMZ,GAMW2,GAMZ2
      COMMON/STAND1/SWS,CWS,XMBOT,XMTOP,XMB2,XMT2,GWEAK2
      COMMON/STAND5/XMCHARM,XMCHARM2,XMHIGGS,GAMH,GAMH2
      COMMON/XPART_MASS/XMASS(100),X_DECAYMASS(100)
      COMMON/MYX12/Y_M1,Y_M2,Y_EMIN
      COMMON/HIGGS_STUFF/ZUP,ZDOWN,XMGW2,XMGZ2,XMGH2,XMH2,VSQ
      COMMON/HTOP/GAMT
      COMMON/MOM/PLAB(4,12),PCEM(4,12),WHAT
      COMMON/OUT_MOMENTA/FINAL_P(4,50),NP_TYPE(50),N_PARTS
      integer nflmax
      parameter (nflmax=5)
      COMMON/STRUCTURE/q1(-nflmax:nflmax),q2(-nflmax:nflmax)
      COMMON/VARS/PT(10),ETA(10),DR(10,10),ETOT
      COMMON/SEC_VARS/SEC_PT(10),SEC_ETA(10)
      COMMON/QGWTB/Q1Q(8),Q2Q(8),Q1W(8),Q2W(8),Q_HOLD
      COMMON/PMF/X1,X2
      COMMON/RESONANCE/I_RESONANCE
      COMMON/SCALE/Q
      COMMON/NININF/IBINNER
      COMMON/BVEG1/XL(50),XU(50),ACC,NDIM,NCALL,ITMX,NPRN
      COMMON/BVEG2/XI(50,50),SI,SI2,SWGT,SCHI,NDO,IT,ALPH
      LOGICAL SMEARING
      REAL*8 JET_SM_RANGE
c      COMMON/JET_SMEAR/SMEARING,JET_SM_RANGE
      COMMON/JET_SMEAR/JET_SM_RANGE,SMEARING

      REAL*8 PTECUT,YECUT,DELR_CUT,PTNCUT,yNcut, ywucut, ywlcut
      COMMON/MYCUT/PTECUT,YECUT,DELR_CUT,PTNCUT,yNcut, ywucut, ywlcut
      REAL*8 YM12,YPT12,YP12
      COMMON/MY1/YM12,YPT12,YP12(4)
      REAL*8 YXP,YXETA,YXPT,YXDELR,YXPP
      COMMON/SECCONV/ YXP(4,30),YXETA(30),YXPT(30),YXDELR(30,30),
     >YXPP(30)
      REAL*8 YP,YETA,YPT,YDELR
      COMMON/PASSCONV/ YP(4,30),YETA(30),YPT(30),YDELR(30,30)
      REAL*8 ETOP,THE_STA,PHI_STA,RFWP,RFWM,RCWP,RCWM,
     >RFEP,RFEM,RCEP,RCEM
      COMMON/CMS/ETOP,THE_STA,PHI_STA,RFWP,RFWM,RCWP,RCWM,
     >RFEP,RFEM,RCEP,RCEM

      REAL*8 PI
      REAL*8 CONVERSION,PJ,PBOO,PL,PCM,PREVENT_UNDERFLOW
      REAL*8 WEIGHT,Y,ZTEMP,JACOBIAN,PHASE_WT,WTME,FLUX
c      REAL*8 E_BEFORE,E_AFTER,E_SIGMA,WTA,WTAA
      REAL*8 WTAA
      INTEGER*4 I,J,IPART
      LOGICAL IPASS,PASS_CUTS
      REAL*8 Y_ECM,Y_P1PL,Y_P1MN,Y_P2PL,Y_P2MN,Y_FLUX,Y_Q1CM
      DIMENSION PJ(4,100),PBOO(4),PL(4),PCM(4)
      DIMENSION Y(50)

      Real*8 tmp
      Dimension tmp(4)
      Real*8 QT_V_SEP
      Common / NLO / QT_V_SEP
      Character*40 Boson
      Common / Boson / Boson
      REAL*8 SEC_PLAB
      COMMON/DECAY_MOMENTA/SEC_PLAB(4,50)
      Character*2 hDecayMode
      Common / hDM / hDecayMode

cpn2008 for reweighting
      real*8 phase_wt_old, weight_old, weight_new
      logical ipass_old
      integer icut_mismatch
      common/reweight/weight_new, weight_old, icut_mismatch

      real*8 MAX_WGT,NEW_MAX_WGT,WGT_FACTOR
      INTEGER*4 NBINIT, iwgt
      COMMON/UNWEIGHT/MAX_WGT,NEW_MAX_WGT,WGT_FACTOR,IWGT,NBINIT


      DATA CONVERSION/389385.73D3/,PREVENT_UNDERFLOW/1.D8/
      Data Pi / 3.141592653589793 /

      IF(IBINNER)THEN
        DO 92 I=1,4
          DO 91 J=1,50
            FINAL_P(I,J)=0.D0
91        CONTINUE
92      CONTINUE
      END IF
      JACOBIAN=1.
      IPASS=.FALSE.
      WEIGHT=0.D0

cpn2008 Initialize the matrix element weight for PDF reweighting
      wtme=0d0

!     Print*, I_RESONANCE, NPART
CsB___Generate the phase space
      IF(I_RESONANCE.EQ.1) THEN ! processes with grids ----------------------
        IF(NPART.EQ.11 .OR. NPART.EQ.12) THEN
CCPY!!          Print *,' Calling WRESPH '
CCPY!          CALL WRESPH(Y,PHASE_WT) ! PS for resummed and NLO real emissions
          IF(Boson.Eq.'H0'.and.
     >     (hDecayMode.Eq.'WW'.or.hDecayMode.Eq.'ZZ'))then
            CALL HVVRESPH(Y,PHASE_WT)
          Else
            CALL WRESPH(Y,PHASE_WT) ! PS for W->e nu with initial state resum
          EndIf          

        ELSEIF(NPART.EQ.13 .OR. NPART.EQ.14) THEN
          CALL WLOPH(Y,PHASE_WT)  ! PS for LO and NLO virtual corrections
        elseif (npart.eq.21) then      
          CALL AARESPH(Y,PHASE_WT) ! PS for AA+JET, WRESPH-like routine
CCPY NOv. 2015:
C        ELSEIF(NPART.EQ.22) THEN
        ELSEIF(NPART.EQ.22 .OR. NPART.EQ.23) THEN
          CALL AAJPH(Y,PHASE_WT) ! PS for AA+JET, standard phase space routine
CCAO
        ELSEIF(NPART.EQ.15) THEN
          CALL WRES3PH(Y,PHASE_WT)
        ELSEIF(NPART.EQ.16) THEN
          CALL WLO3PH(Y,PHASE_WT)
        ELSEIF(NPART.EQ.50) THEN
          CALL WGRIDPH(Y,PHASE_WT)
        ELSE
          WRITE (UNIT=*,FMT='(A)')
     >  ' YOU HAVE ASKED FOR A NON-EXISTENT PROCESS'
          STOP
        END IF
        IF(PHASE_WT.LT.0.D0) goto 320
        Y_FLUX=.5/(X1*X2*W*W)
      ELSE ! processes with built-in amplitudes -----------------------------
        IF(Y_M1.LT.0.1) Y_M1=0.0
        IF(Y_M2.LT.0.1) Y_M2=0.0
        X1X2MIN=(Y_EMIN**2-Y_M1**2-Y_M2**2)/W**2
        ZTEMP=1.D0/((1./X1X2MIN-1.)*Y(1)+1.)
        X1=ZTEMP**Y(2)
        X2=ZTEMP/X1
        JACOBIAN=-ZTEMP**2*LOG(ZTEMP)*(1./X1X2MIN-1.)
        WHAT=SQRT(X1*X2*W*W+Y_M1**2+Y_M2**2)
        Y_ECM=0.5*W
        Y_P1PL=SQRT(2.0)*Y_ECM*X1 ! p1+ = x1 Sqrt(S/2)
        Y_P2MN=SQRT(2.0)*Y_ECM*X2 ! p2- = x2 Sqrt(S/2)
        Y_P1MN=0.5*Y_M1**2/Y_P1PL ! p1-
        Y_P2PL=0.5*Y_M2**2/Y_P2MN ! p2+
        Y_FLUX=1.0/4.0/SQRT((Y_P1PL*Y_P2MN+Y_P1MN*Y_P2PL)**2-
     >  Y_M1**2*Y_M2**2)
        Y_Q1CM=SQRT(((WHAT**2+Y_M1**2-Y_M2**2)/2.D0/WHAT)**2-Y_M1**2)
        PLAB(1,1)=0.D0
        PLAB(2,1)=0.D0
        PLAB(3,1)=(Y_P1PL-Y_P1MN)/SQRT(2.0)
        PLAB(4,1)=(Y_P1PL+Y_P1MN)/SQRT(2.0)
        PLAB(1,2)=0.D0
        PLAB(2,2)=0.D0
        PLAB(3,2)=(Y_P2PL-Y_P2MN)/SQRT(2.0)
        PLAB(4,2)=(Y_P2PL+Y_P2MN)/SQRT(2.0)
        PCEM(1,1)=0.D0
        PCEM(2,1)=0.D0
        PCEM(1,2)=0.D0
        PCEM(2,2)=0.D0
        PCEM(3,1)=Y_Q1CM
        PCEM(4,1)=SQRT(Y_Q1CM**2+Y_M1**2)
        PCEM(3,2)=-Y_Q1CM
        PCEM(4,2)=SQRT(Y_Q1CM**2+Y_M2**2)
        CALL PHASE_GEN(Y,N_PRIMARY_PART,WHAT,XMASS,PJ,PHASE_WT)
        PHASE_WT=PHASE_WT*(2.*PI)**(4.-3*N_PRIMARY_PART)
        DO 53 I=1,4
          PBOO(I)=PLAB(I,1)+PLAB(I,2)
53      CONTINUE
        DO 57 J=1,N_PRIMARY_PART
          DO 61 I=1,4
            PCM(I)=PJ(I,J)
            PCEM(I,J+2)=PJ(I,J)
61        CONTINUE
          CALL LORENTZ(PBOO,PCM,PL)
          Do I=1,4
            PLAB(I,J+2)=PL(I)
          End Do
57      CONTINUE
      END IF ! --------------------------------------------------------------

CsB___Construct pT and y for primary particles
      Do IPART=1,N_PRIMARY_PART
        PT(IPART)=SQRT(PLAB(1,IPART+2)**2+PLAB(2,IPART+2)**2)
        IF((PLAB(4,IPART+2)-PLAB(3,IPART+2)).GT.0.0.AND.
     .  (PLAB(4,IPART+2)+PLAB(3,IPART+2)).GT.0.0)THEN
          ETA(IPART)=.5*LOG((PLAB(4,IPART+2)+PLAB(3,IPART+2))
     .    /(PLAB(4,IPART+2)-PLAB(3,IPART+2)))
        ELSEIF((PLAB(4,IPART+2)-PLAB(3,IPART+2)).LE.0.)THEN
          ETA(IPART)=20.
        ELSEIF((PLAB(4,IPART+2)+PLAB(3,IPART+2)).LE.0.)THEN
          ETA(IPART)=-20.
        END IF
      End Do

CsB___Construct pT and y for secondary decay particles
      Do IPART=1,N_SECOND_PART
        SEC_PT(IPART)=SQRT(SEC_PLAB(1,IPART)**2+SEC_PLAB(2,IPART)**2)
        IF((SEC_PLAB(4,IPART)-SEC_PLAB(3,IPART)).GT.0.0.AND.
     .     (SEC_PLAB(4,IPART)+SEC_PLAB(3,IPART)).GT.0.0)THEN
          SEC_ETA(IPART)=.5*LOG((SEC_PLAB(4,IPART)+SEC_PLAB(3,IPART))
     .                         /(SEC_PLAB(4,IPART)-SEC_PLAB(3,IPART)))
        ELSEIF((SEC_PLAB(4,IPART)-SEC_PLAB(3,IPART)).LE.0.)THEN
          SEC_ETA(IPART)=20.
        ELSEIF((SEC_PLAB(4,IPART)+SEC_PLAB(3,IPART)).LE.0.)THEN
          SEC_ETA(IPART)=-20.
        END IF
      End Do

CsB___Fill lab momenta for 3&4 when built in amplitudes used
      If (NPart.ge.21.AND.nPart.LE.24) then 
                                 ! built-in double EW gauge-boson production
        Do i=1,4
          PLab(i,9) = PLab(i,3) + PLab(i,4)
        End Do
cdump
c        Print*, ' PLab(3+4) = '
c        Print*, (PLab(j,9), j=1,4)
c        Print*,'Sqrt(PLab(1,9)**2+PLab(2,9)**2) =',
c     > Sqrt(PLab(1,9)**2+PLab(2,9)**2)
c        Print*, ' QT_V_SEP = ', QT_V_SEP
C_______Cut on QT3&4 (QT) minimum
        If (Sqrt(PLab(1,9)**2+PLab(2,9)**2).LT.QT_V_SEP) goto 320
      End If

CCPY APRIL 14, 2004: MAKE SURE THAT PHASE_WT IS ALWAYS POSITIVE FOR 
C NORMAL PROCESSES.
C WE RETURN PHASE_WT < 0 TO IMPLEMENT A SELECTION OF THE MONTE CARLO EVENT
C SIMILAR TO SUBROUTINE PRIMARY_CUTS. (FOR EXAMPLE, WHEN USING PHASE SPACE
C SLICING METHOD IN DOING NLO CALCUALTION INVOLVING REAL EMISSION PROCESS.
      IF(PHASE_WT.LT.0.D0) goto 320

!      Print*, ' Calling PRIMARY_CUTS '
      CALL PRIMARY_CUTS(PASS_CUTS)
      IF(.NOT.PASS_CUTS) goto 320 !to the end of the subroutine
C      Print*, ' Passed primary cuts in Make_Event '

CsB___Renormalization and factorization scale choice
C     Print*, ' ITScale =', ITSCALE
      IF(ITSCALE.EQ.1) THEN
        Q=WHAT
      ELSEIF(ITSCALE.EQ.2) THEN
        ETOT=0.D0
        Do J=1,N_PRIMARY_PART
          ETOT=ETOT+SQRT(PT(J)**2+XMASS(J)**2)
        End Do
        Q=ETOT/FLOAT(N_PRIMARY_PART)
      ELSEIF(ITSCALE.EQ.3)  THEN
        ETOT=0.D0
        DO 217 J=1,N_PRIMARY_PART
          ETOT=ETOT+SQRT(PT(J)**2+XMASS(J)**2)
217     CONTINUE
        Q=ETOT

      ELSEIF(ITSCALE.EQ.4) THEN
        If (Boson(1:2).EQ.'W+' .Or. Boson(1:2).EQ.'W-') Q=XMW
        If (Boson(1:2).EQ.'Z0') Q = xMz
        If (Boson.EQ.'HP' .Or. Boson.EQ.'HM') Q=XMW+xMHiggs
        If (Boson.EQ.'HZ') Q=XMZ+xMHiggs
        If (Boson.EQ.'H0') Q = xMHiggs

      ELSEIF(ITSCALE.EQ.5) THEN
CPN____Factorization scale = renormalization scale = 
C                           = invariant mass of 3&4
        Do i=1,4
          tmp(i) = PLab(i,3) + PLab(i,4)
        End Do
        Q = tmp(4)**2 - tmp(3)**2 - tmp(2)**2 - tmp(1)**2
        Q = DSQRT(Q)

      ELSEIF(ITSCALE.EQ.6) THEN
CsB     Factorization scale = renormalization scale = 1/2*invariant mass of 3&4
        Do i=1,4
          tmp(i) = PLab(i,3) + PLab(i,4)
        End Do
        Q = tmp(4)**2 - tmp(3)**2 - tmp(2)**2 - tmp(1)**2
        Q = DSQRT(Q)/2.d0

      ELSEIF(ITSCALE.EQ.7) THEN
CsB     Factorization scale = renormalization scale = 2*invariant mass of 3&4
        Do i=1,4
          tmp(i) = PLab(i,3) + PLab(i,4)
        End Do
        Q = tmp(4)**2 - tmp(3)**2 - tmp(2)**2 - tmp(1)**2
        Q = DSQRT(Q)*2.d0

      ELSEIF(ITSCALE.EQ.8) THEN
CsB     Factorization scale = renormalization scale = (pT(3)+pT(4))/2
        Q = DSqrt(PLAB(1,3)**2 + PLAB(2,3)**2) +
     +      DSQRT(PLAB(1,4)**2 + PLAB(2,4)**2)
        Q = Q/2.
      ELSE
        PRINT*,' NO SUCH ITSCALE'
        STOP
      END IF

CCPY DUMMY VARIABLE
C      NF=4
C      IF(Q.GT.2.*XMBOT) NF=5

      FLUX=Y_FLUX
CCPY
      IF ( NPART.EQ.11 .OR. NPART.EQ.12 .OR. NPART.EQ.13
     >.OR. NPART.EQ.14 .OR. NPART.EQ.15 .OR. NPART.EQ.16) GOTO 224 ! No calls to built in Amp's/PDF's

      IF(NPART.EQ.50) GOTO 224

CsB___Get the parton luminosities
      CALL YPARTONS
      DO 223 I=-nflmax, nflmax
        Q1(I)=Q1(I)*PREVENT_UNDERFLOW
        Q2(I)=Q2(I)*PREVENT_UNDERFLOW
223   CONTINUE

224   CONTINUE

C     NPART = 11   FOR RESUMMED; use WRESPH and WPRES
C     NPart = 12   for "NLO"     (Asymptotic piece AND PERTURBATIVE PIECE);
C                               use WRESPH and WPRES 
C
C     NPart = 13   for "NLO_Sig" (Delta sigma piece); use WLOPH and WPRES
C     NPART = 14   FOR LEADING ORDER INITIAL STATE; 
C                       WITH FINAL STATE NLO QED (S+V) CORRECTIONS; 
C                       use WLOPH and WPRES_RA
C
C     NPART = 15   RESUMMED INITIAL STATE (CSS PIECE ONLY) 
C                  WITH FINAL STATE NLO QED (REAL) EMISSION; 
C                  use WRES3PH and WGRES
C
C     NPART = 16   LEADING ORDER INITIAL STATE 
C                  WITH FINAL STATE NLO QED (REAL) EMISSION;
C                  use WLO3PH and WGRES
C
CsB___Call the part of the code that calculates the amplitude squared
C      Print*, ' nPart =', nPart
      IF(NPART.EQ.11 .OR. NPART.EQ.12 .OR.
     >   NPART.EQ.13 .OR. NPART.EQ.14) THEN
CCAO
        IF(BOSON.EQ.'W+_RA' .OR. BOSON.EQ.'W-_RA' .OR.
     >     BOSON.EQ.'Z0_RA' .OR. 
     >     BOSON.EQ.'W-A(SV)' .OR. BOSON.EQ.'W+A(SV)'
     >     .OR. BOSON.EQ.'Z0A(SV)') THEN
          CALL WPRES_RA(WTME)
CCPY
        ELSEIF(Boson.Eq.'H0'.and.
     >         (hDecayMode.Eq.'ZZ'.or.hDecayMode.Eq.'WW'))then
          CALL WPRES_HVV(WTME)
        ELSE  
          CALL WPRES(WTME)
c          wtme=1
        ENDIF
CCAO 
      ELSEIF(NPART.EQ.15.OR.NPART.EQ.16) THEN
        IF(BOSON.EQ.'W+A(R)' .OR. BOSON.EQ.'W-A(R)') THEN
          CALL WGRES(WTME)
        ELSEIF(BOSON.EQ.'Z0A(R)') THEN
          CALL ZGRES(WTME)
        ENDIF        
      ELSEIF(NPART.EQ.21.or.NPART.eq.22) THEN ! G G -> A A G  NLO
        
        if (Boson.Eq.'AA') then                               
c          CALL QQBZZG(WTME)     ! Built-in amplitude**2 for q Q(G) -> A A X  NLO 
          CALL QQAAG(WTME)   ! Built-in amplitude**2 for q Q(G) -> A A X  NLO 
c          wtme=1d12
        elseif (Boson.Eq.'AG') then
          CALL GGAAG(WTME)      ! Built-in amplitude**2 for G G -> A A X  NLO 
c          wtme=1d12
        endif                   ! Boson

      ELSEIF(NPART.EQ.23) THEN  ! q Q -> Z Z X  NLO
        CALL QQBZZG(WTME)      ! Call built in amplitude**2
      ELSEIF(NPART.EQ.24) THEN ! parton parton' -> H0 X -> Z Z X (NLO)
        CALL pp_H0X_ZZX(WTME)  ! Call built in amplitude**2
      ELSEIF(NPART.EQ.50) THEN ! Manipulating Legacy Grid files
        CALL WGRID(WTME)      ! Call built in amplitude**2
      ELSE
        WRITE (UNIT=*,FMT='(A)')
     .  ' YOU HAVE ASKED FOR A NON-EXISTENT PROCESS'
        STOP
      END IF


C No smearing of the jet momenta to take into account detector
C resolution
      WTAA=1.

      CALL SECONDARY_CUTS(PASS_CUTS)
      IF(.NOT.PASS_CUTS) goto 320     !FAILED
      IPASS=.TRUE.                   ! PASSED CUTS

CsB___Calculate the weight: phase weight * amplitude squared
CCAO  add in npart=15 and npart=16
      IF(NPART.EQ.11 .OR. NPART.EQ.12 .OR. NPART.EQ.13 .OR.
     >   NPART.EQ.14 .OR. NPART.EQ.15 .OR. NPART.EQ.16) THEN ! Processes with grids
         WEIGHT = PHASE_WT * WTME
C yfu test
c         print*,phase_WT,WTME
C         If (Abs(WEIGHT).LT.1.d-9)
C     >     Print*, ' PHASE_WT,WTME :', PHASE_WT,WTME
      ELSEIF(NPART.EQ.50) THEN      
         WEIGHT = PHASE_WT * WTME
      ELSE                     ! Processes with built-in phase space routines
C        Jacobian   = J
C        Flux       = 1/4/(s/2)=1/2/s
C        Conversion = 389385.73d3, from 1/GeV^2 to picobarns
C        Weight     = 1 J WTp WTa C 1/2/s/x1/x2 1/UF^2
C        Note: In ResBos the PDF returns x*f (momentum distribution)
C              and the 1/x1/x2 factor will be cancelled by the x1 x2 factor
C              in the WTa which contains the x1*f(x1)*x2*f(x2).
         WEIGHT=WTAA*JACOBIAN*PHASE_WT*FLUX*CONVERSION
     1          /X1/X2*WTME
     1          /PREVENT_UNDERFLOW/PREVENT_UNDERFLOW
      ENDIF

C      Print*, ' WTAA,PHASE_WT,WTME,WEIGHT '
C      Print*,   WTAA,PHASE_WT,WTME,WEIGHT

cpn2008 If the PDF reweighting is turned on
 320  if (iwgt.eq.-1) then !store PDF weights 
        write (35, 321) weight, phase_wt, ipass
      else if (iwgt.eq.-2) then 
cpn2008 Read the phase space weight, weight, and the flag ipass (specifying
c       if the cuts were passed) from the run for iwgt=-1. Stop if the
c       old (iwgt=-1) and new (iwgt=-2) phase_weights do not agree.
c       Save the new weight into a variable weight_new
cp      and then overwrite the variable weight with the old weight value. In this
c       way, an old value of the weight is passed to VEGAS, thus
c       ensuring that the sequence of the MC events in the new run is
c       exactly the same as in the old run. The variable weight_new is
c       then passed into the subroutine PASS_OUTPUT to write the correct
c       new output into the ntuple.
c  
c       In a very small number of events, the kinematical cuts may be
c       passed differently in the old and new runs because of roundoff
c       errors in PRIMARY_CUTS. To identify and fix such mismatch, compare 
c       the old and new ipass values, and then always overwrite the new
c       ipass with the old ipass. If the new and old ipass differ, write
c       an old weight into weight_new rather than the new weight.

        read (35, 321,err=322,end=322) weight_old, phase_wt_old,
     >    ipass_old

C_yfu for AFB factorization //1e-6 -> 1e-1
       if (abs((phase_wt-phase_wt_old)/max(1e-16,phase_wt_old))
     >    .gt.1e-1) then
          print *,
     >      'STOP: the new and old phase space weights differ:'
          print 321, phase_wt,phase_wt_old
          print *,
     >      'Did the input parameters change compared to the old run?'
          call kinemprint(weight_old,phase_wt_old, weight, phase_wt)
          call quit
        endif                   !abs(...


      if (ipass.neqv.ipass_old) then
c        print *,'make_event: ipass != ipass_old;','
c     >    ipass, ipass_old, weight, weight_old = ',
c     >    ipass, ipass_old, weight, weight_old
        weight_new = weight_old 
        icut_mismatch = icut_mismatch + 1
      else
        weight_new = weight
      endif                     !
      
      ipass = ipass_old; weight = weight_old
      endif                     !iwgt
      

 321  format (g22.16,2x,g22.16,2x,l1) 

      RETURN
 322  print *,'I/O error: check the file weights.dat'
      call Quit
      return 
      
      END ! MAKE_EVENT

CsB Initialize certain variables
C --------------------------------------------------------------------------
      SUBROUTINE STARTHERE
C --------------------------------------------------------------------------
      IMPLICIT NONE

      REAL*8 W,QCDL,X1X2MIN
      INTEGER*4 NJET,N_PRIMARY_PART,NF,I_STR,IBEAM,NPART,NEVENT,NORDER,
     ,          N_SECOND_PART
      REAL*8 XL,XU,ACC
      REAL*8 XI,SI,SI2,SWGT,SCHI,ALPH
      INTEGER*4 NDIM,NCALL,ITMX,NPRN,NDO,IT,
     &IW,IAPPROX
      REAL*8 ETOTMN,PTMIN,ETAMAX,ETASEP
      INTEGER*4 I_DECAY,IVEGOUT
      LOGICAL I_READIN,I_WRITEOUT,I_HADRONIZE
      CHARACTER*50 SAVE_FILE,RESTORE_FILE
      REAL*8 HIGGS_MASS
      INTEGER*4 ITSCALE
      REAL*8 SWS,CWS,XMBOT,XMTOP,XMB2,XMT2,GWEAK2,
     &XMASS,X_DECAYMASS,
     &MAX_WGT,NEW_MAX_WGT,WGT_FACTOR
      LOGICAL UNWGT_FLAG
      INTEGER*4 NBINIT, iwgt
      LOGICAL CHARM_ONLY,BOTTOM_ONLY,
     &SMEARING
      REAL*8 JET_SM_RANGE,
     &Y_M1,Y_M2,Y_EMIN
      INTEGER*4 I_RESONANCE

      COMMON/AA/W,QCDL,X1X2MIN,NJET,N_PRIMARY_PART,NF,I_STR,IBEAM,NPART,
     &  NEVENT,NORDER,N_SECOND_PART
      COMMON/BVEG1/XL(50),XU(50),ACC,NDIM,NCALL,ITMX,NPRN
      COMMON/BVEG2/XI(50,50),SI,SI2,SWGT,SCHI,NDO,IT,ALPH
      COMMON/CHARGE/IW
      COMMON/CUTVAL/ETOTMN,PTMIN,ETAMAX,ETASEP,IAPPROX
      COMMON/DEC_AY/I_DECAY,I_HADRONIZE
      COMMON/DETAILS/IVEGOUT,I_READIN,I_WRITEOUT,SAVE_FILE,RESTORE_FILE
      COMMON/HIGGS/HIGGS_MASS
      COMMON/PAT_SCALE/ITSCALE
      COMMON/STAND1/SWS,CWS,XMBOT,XMTOP,XMB2,XMT2,GWEAK2
      COMMON/XPART_MASS/XMASS(100),X_DECAYMASS(100)
      COMMON/UNWEIGHT/MAX_WGT,NEW_MAX_WGT,WGT_FACTOR,IWGT,NBINIT
      COMMON/HEAVY_FLAGS/CHARM_ONLY,BOTTOM_ONLY
c      COMMON/JET_SMEAR/SMEARING,JET_SM_RANGE
      COMMON/JET_SMEAR/JET_SM_RANGE,SMEARING
      COMMON/MYX12/Y_M1,Y_M2,Y_EMIN
      COMMON/RESONANCE/I_RESONANCE

      INTEGER I

C INITILIZE ALL THE VARIABLES
      NJET=0
      IW=1
      I_STR=15
      ITSCALE=1
      DO 53 I=1,100
        XMASS(I)=0.0
53    CONTINUE
CsB___Higgs and top masses are read from grid or set in SUBROUTINE GENERATOR
C      HIGGS_MASS=100.0
C      XMTOP=175.0
      PTMIN=0.0
      ETOTMN=0.0
      ETAMAX=99999.0
      ETASEP=0.0
      IAPPROX=0
      I_RESONANCE=0
      I_HADRONIZE=.FALSE.
      NF=4
      IVEGOUT=1
      ITMX=1
      NEVENT=200
      I_READIN=.FALSE.
      I_WRITEOUT=.FALSE.
      UNWGT_FLAG=.FALSE.
      WGT_FACTOR=1.
      CHARM_ONLY=.FALSE.
      BOTTOM_ONLY=.FALSE.
      SMEARING=.FALSE.
      Y_M1=0.0
      Y_M2=0.0
      IF(JET_SM_RANGE.LT.0.5) JET_SM_RANGE=20.0

      RETURN
      END ! STARTHERE

C --------------------------------------------------------------------------
      SUBROUTINE STANDARD_old
C --------------------------------------------------------------------------
CsB   Updated: 2004 Apr. 20
      IMPLICIT NONE
      INTEGER*4 I,J
      REAL*8 TEMP,A_W,XMQX

      REAL*8 HIGGS_MASS,
     &ZUP,ZDOWN,XMGW2,XMGZ2,XMGH2,XMH2,VSQ,
     &XMW,XMZ,XMW2,XMZ2,AEM,GAMW,GAMZ,GAMW2,GAMZ2,
     &SWS,CWS,XMBOT,XMTOP,XMB2,XMT2,GWEAK2,
     &XLU,XLD,XRU,XRD,QU,QD,GFERMI,XLE,XRE,
     &XLU2,XLD2,XRU2,XRD2,QU2,QD2,XLE2,XRE2,
     &XLU4,XLD4,XRU4,XRD4,QU4,QD4,PI,
     &XMCHARM,XMCHARM2,XMHIGGS,GAMH,GAMH2,
     &CURRENTW_MASS,CURRENTZ_MASS,
     &GAMT,
     &g_fVf_P_u,g_fVf_M_u,g_fVf_P_d,g_fVf_M_d

      COMMON/HIGGS/HIGGS_MASS
      COMMON/HIGGS_STUFF/ZUP,ZDOWN,XMGW2,XMGZ2,XMGH2,XMH2,VSQ
      COMMON/STAND0/XMW,XMZ,XMW2,XMZ2,AEM,GAMW,GAMZ,GAMW2,GAMZ2
      COMMON/STAND1/SWS,CWS,XMBOT,XMTOP,XMB2,XMT2,GWEAK2
      COMMON/STAND2/XLU,XLD,XRU,XRD,QU,QD,GFERMI,XLE,XRE
      COMMON/STAND3/XLU2,XLD2,XRU2,XRD2,QU2,QD2,XLE2,XRE2
      COMMON/STAND4/XLU4,XLD4,XRU4,XRD4,QU4,QD4,PI
      COMMON/STAND5/XMCHARM,XMCHARM2,XMHIGGS,GAMH,GAMH2
      Common/Stand6/g_fVf_P_u,g_fVf_M_u,g_fVf_P_d,g_fVf_M_d
      COMMON/RUNNING_BOSON/CURRENTW_MASS,CURRENTZ_MASS
      COMMON/HTOP/GAMT
      REAL*8 KM_L,KM_A,KM_RHO,KM_ETA,VKM
      COMMON/CKM/VKM(6,6)

      REAL*8 YM12MIN,YM12MAX,QTMin,QTMax, TMmin, Tmmax, ETmin
      COMMON/MYYM12/YM12MIN,YM12MAX,QTMin,QTMax, TMmin, Tmmax, ETmin
CCAO!
      REAL*8 GAMW_NLO,GAMZ_NLO,ALFA_EM_MZ,ACOUPL,WCOUPL,
     &ZCOUPL,HCOUPL,SW2_EFF,EW_ALFA
      COMMON/STAND7/ GAMW_NLO,GAMZ_NLO,ALFA_EM_MZ,ACOUPL,WCOUPL,
     >       ZCOUPL,HCOUPL,SW2_EFF,EW_ALFA

CCAO! THE COUPLING OF Z-F-F
      REAL*8
     >   EL_Z,ER_Z,QL_Z,QR_Z,EL_A,ER_A,QL_A,QR_A
      COMMON/ZFF_INIT/
     >   EL_Z,ER_Z,QL_Z(2),QR_Z(2),EL_A,ER_A,QL_A(2),QR_A(2)
CCAO! The readin W-boson mass and width
      Real*8 XMW_IN,XMW_GRID,GAMW_IN
      COMMON/WMASS_IN/XMW_IN,XMW_GRID,GAMW_IN
      Real*8 XMZ_IN,XMZ_GRID,GAMZ_IN
      COMMON/ZMASS_IN/XMZ_IN,XMZ_GRID,GAMZ_IN
     
      Integer BIPDF, ISet
      Common / BuildInPDF / BIPDF, ISet
      
      REAL*8 Y_EW,Y_PW
      REAL*8 TOP_MASS
      Real*8 GMu
      Real*8 Pi2,R2,
     &  EW_DELA,EW_DELRHO,ALPI_MT,EW_RHO,EW_DELKA,EW_DELR,
     &  Sw_Eff,Cw_Eff,Tw_Eff,T3u,T3d

      REAL*8 GAMW_MW,GAMZ_MZ
      REAL*8 GAMW_NLO_R,GAMZ_NLO_R
      REAL*8 BR_Z_EE_NLO_R,BR_W_EN_NLO_R
      REAL*8 BR_Z_EE,BR_W_EN
      LOGICAL FIRST
      DATA FIRST/.TRUE./
CCAO! 
      REAL*8 UNIT,HBARC2      
      REAL*8 GWEAK,EW_AMZ,EW_A,AEM_MZ
      REAL*8 ZCOUPL_OLD,ZCOUPL_MZ
      REAL*8 ALPHAS_MW,ALPHAS_MZ
      REAL*8 QE,VUB
      REAL*8 SW,CW,SW2_EFF_OLD,SW2_EFF_MZ
      REAL*8 A02

      REAL*8 ALPHAS,ALFAS
      EXTERNAL ALPHAS

      REAL*8 H2BB,H2TT,RUN_MASS_TOT
      EXTERNAL RUN_MASS_TOT
      
      REAL*8 GAMH_TMP, HWIDTH
      
CCHEN!********************** 
      REAL*8 CWW,H2WW,CZZ,H2ZZ,HV,V
      REAL*8 EWC_HWW,COEFF,AA1,AA2,FI1,XX,
     .       TERM1,TERM2,TERM3,TERM4,TERM5,DEL12SQ
      EXTERNAL FI1
      HV(V)=3.D0*(1.D0-8.D0*V+20.D0*V**2)/DSQRT((4.D0*V-1.D0))
     .      *DACOS((3.D0*V-1.D0)/2.D0/DSQRT(V**3))
     .      -(1.D0-V)*(47.D0/2.D0*V-13.D0/2.D0+1.D0/V)
     .      -3.D0/2.D0*(1.D0-6.D0*V+4.D0*V**2)*DLOG(V)
c********************************
C      print*,' inside standard: XMW, XMZ =',XMW,XMZ

CCAO!      Data PI/3.14159 26535 89793 24D0/
      PI=4.0D0*ATAN(1.0D0)      

C_____Parameters
      Pi2 = Pi**2
      R2 = Sqrt(2.d0)
CCAO!
C UNIT: to convert energy from  GeV  to TeV, then unit=1.0d-3
c WARNING!!! You must keep UNIT=1.0D0 in this version of the code!!!
      UNIT=1.0D0
C Convert from GeV to pb
      HBARC2=UNIT**2*0.38937966D0*1.0D9
      
CsB   EW couplings
      GMU = 1.16637D-5/UNIT**2
      EW_ALFA = 1.0D0/137.0359895D0

CsB___Quark masses
C      XMCHARM=1.5D0
      XMCHARM=1.4D0
      XMCHARM2=XMCHARM**2
CCPY      XMBOT=5.0D0
C      XMBOT=4.5D0
      XMBOT=4.75D0
      XMB2=XMBOT**2
      TOP_MASS = xMTop
C      Print*, ' TOP_MASS =', TOP_MASS

CsB___W and Z mass
C     The W and Z masses are read from the grid.
      XMZ2=XMZ**2
      XMW2=XMW**2
      CURRENTZ_MASS=XMZ
      CURRENTW_MASS=XMW
C_____Test W, Z and top masses
!     Print*, ' <STANDARD> W, Z masses:', XMW,XMZ

C_____EW_DELA=0.0595 +/- 0.0009
CCAO!      EW_DELA=0.0595+40.0D0/9.0D0*EW_ALFA/PI*DLOG(xMZ/91.187D0)
C          EW_DELA=0.02767+0.03150-0.00007=0.0591
      EW_DELA=0.0591+40.0D0/9.0D0*EW_ALFA/PI*DLOG(XMZ/91.187D0/UNIT)

C_____Rho and Delta rho
C EW_ALFA AT  MZ
      EW_AMZ=EW_ALFA/(1.0D0-EW_DELA)
      EW_A=PI*EW_ALFA/R2/GMU*UNIT**2
      ALPI_MT=ALPHAS(XMTOP)/PI      
C     PRINT*,' XMTOP, ALPHAS(XMTOP)=',XMTOP, ALPI_MT*PI
      EW_DELRHO=3.0D0*GMU*TOP_MASS**2/8.0D0/R2/PI2
      EW_DELRHO=EW_DELRHO*(1.D0-2.0/3.0*(PI2/3.0+1.0)*ALPI_MT)
      EW_RHO=1.0D0/(1.0D0-EW_DELRHO)

CsB___Sin^2 and Cos^2 of the weak angle - in the on shell scheme
      CWS=XMW2/XMZ2
      SWS=1.D0-CWS
      SW=DSQRT(SWS)
      CW=DSQRT(CWS)
CCPY  EW_DELKA=CWS/SWS*EW_DELRHO
C     EW_DELR=EW_DELA-EW_DELKA
C     print*,' old ew_delr,EW_DELA,EW_DELKA =',ew_delr,EW_DELA,EW_DELKA
CPN___Calculate delta_r using mW
      A02 = PI*EW_AMZ/R2/GMU
      EW_DELR = 1 - A02/XMW2/SWS 
      EW_DELKA =  EW_DELA - EW_DELR
C     print*,' new ew_delr,EW_DELKA  =',ew_delr,EW_DELKA 
CCAO! SW2_EFF=SWS*(1.0d0 + EW_DELKA)
      SW2_EFF_OLD=SWS*(1.0d0 + EW_DELKA)

CCPY FIXED A SCHEME HERE
CCPY!      SW2_EFF_MZ=0.23143
      SW2_EFF_MZ=0.231676     ! Zfitter (mt=172, mH=115)
cBIN to test the cross section
c      SW2_EFF_MZ=0.221588

      SW2_EFF=SW2_EFF_MZ
      Sw_Eff = Sqrt(SW2_EFF)
      Cw_Eff = Sqrt(1.d0 - SW2_EFF)
      Tw_Eff = Sw_Eff/Cw_Eff

C_____alpha = e^2/(4 Pi)
C      AEM=1.D0/128.0
      AEM = Sqrt(2.d0)*XMW2*GMU*SWS/Pi
      AEM_MZ=1.0/128.937
      ALFA_EM_MZ=EW_AMZ
C_____Test alpha_em
      WRITE(16,*)'1/AEM, 1/EW_AMZ, 1/AEM_MZ =',1/AEM,1/EW_AMZ,1/AEM_MZ

C_____Weak couplings: g_w and Fermi constant
      GWEAK2=AEM*4.d0*PI/SWS
      GWEAK = Sqrt(GWEAK2)
      GFERMI=PI*AEM/SQRT(2.d0)/SWS/XMW2

C_____Quark charges
      Qu  =  2.d0/3.d0
      Qd  = -1.d0/3.d0
      QU2=QU**2
      QU4=QU2**2
      QD2=QD**2
      QD4=QD2**2

C_____Quark weak-isospins
      T3u =  1.d0/2.d0
      T3d = -1.d0/2.d0

C_____Couplings for fermion-Z-fermion vertex in PRD43('91)3626 Eq.(3)
      g_fVf_P_u = - Qu*Tw_Eff                   ! for f=u type quark & V = Z
      g_fVf_M_u = T3u/Sw_Eff/Cw_Eff - Qu*Tw_Eff ! for f=u type quark & V = Z
      g_fVf_P_d = - Qd*Tw_Eff                   ! for f=d type quark & V = Z
      g_fVf_M_d = T3d/Sw_Eff/Cw_Eff - Qd*Tw_Eff ! for f=d type quark & V = Z


C_____Lepton and quark couplings
CCPY: ERROR      XLE=1.D0+2.*SWS
      QE=-1.0
      XLE=-1.D0-2.*SW2_EFF*QE
      XRE=-2.D0*QE*SW2_EFF
      XLU=1.D0-2.*SW2_EFF*QU
      XLD=-1.D0-2.*SW2_EFF*QD
      XRU=-2.*QU*SW2_EFF
      XRD=-2.*QD*SW2_EFF
      XLE2=XLE**2
      XRE2=XRE**2
      XLU2=XLU**2
      XLD2=XLD**2
      XRU2=XRU**2
      XRD2=XRD**2
      XLU4=XLU**4
      XLD4=XLD**4
      XRU4=XRU**4
      XRD4=XRD**4

C_____Test top mass - it is read from the grid or from the input file
C      Print*, ' <STANDARD> top mass:', XMTOP,TOP_MASS
      IF(XMTOP.LT.XMBOT) XMTOP=TOP_MASS
      XMT2=XMTOP**2
C_____Calculate top width
      IF(XMTOP.GT.XMW+XMBOT) THEN
        Y_EW=(XMTOP**2+XMW**2-XMBOT**2)/2.0/XMTOP
        Y_PW=SQRT(Y_EW**2-XMW**2)
        GAMT=GFERMI*Y_PW/2./PI/SQRT(2.)*
     >      (2.0*Y_PW**2+3.0*XMW**2*(1.0-Y_EW/XMTOP))
      ELSE
        GAMT=3.0*9.0*XMTOP**5*GFERMI**2/(192.0*PI**3)
      ENDIF

C     ACOUPL IS (e**2)/2
      ACOUPL=R2*XMW2*GMU*2.0D0*SWS
C     WCOUPL IS (g**2)/8
      WCOUPL=R2*XMW2*GMU/2.0D0
C     WE IGNORE THE DIFFERENCE BETWEEN b-quark AND OTHER FERMIONS
C     ZCOUPL IS (g/Cos_w)**2/4
      ZCOUPL_OLD=R2*XMZ2*GMU*(1.0+EW_DELRHO)
      ZCOUPL_MZ=PI*ALFA_EM_MZ/SWS/CWS
      ZCOUPL=ZCOUPL_MZ      
c     WRITE(*,*) 'ZCOUPL_OLD, ZCOUPL_MZ, ZCOUPL =',
c    > ZCOUPL_OLD, ZCOUPL_MZ, ZCOUPL
C     HCOUPL IS FOR HIGGS
      HCOUPL=R2*GMU*PI/576.0D0

C_____W width
      GAMW=GFERMI*XMW**3/6./PI/SQRT(2.)*(9.)
      GAMW2=GAMW**2
C     NLO QCD CORRECTION
CCPY  ALPHAS_MW=ALPI(XMW)*PI      
      ALPHAS_MW=ALPHAS(XMW)
C     PRINT*,' XMW, ALPHAS_MW =',XMW, ALPHAS_MW 
      GAMW_NLO_R=GAMW*(1.0+2.0*ALPHAS_MW/PI/3.0)

C_____Z width & Br
CCPY Aug 2016      GAMW_MW=2.093*(XMW/80.385)**3      
CCPY PDG 2105      GAMW_MW=2.085*(XMW/80.385)**3      
CCPY      GAMW_MW=2.118*(XMW/80.426)**3      
      GAMW_MW=GAMW_IN

      WRITE(*,*) ' ALPHAS_MW, GAMW, GAMW_NLO_R, GAMW_MW =',
     >            ALPHAS_MW, GAMW, GAMW_NLO_R, GAMW_MW  
CCPY      GAMW_NLO=GAMW_MW
      GAMW_NLO=GAMW_NLO_R

      TEMP=XMZ**3/DSQRT(2.d0)/12.d0/PI*GFERMI
      GAMZ=TEMP*(3.d0+3.d0*(XLE**2+XRE**2)+
     +     6.d0*(XLU**2+XRU**2)+9.*(XLD**2+XRD**2))
      GAMZ2=GAMZ**2     
C     NLO QCD CORRECTION
CCPY  ALPHAS_MZ=ALPI(XMZ)*PI
      ALPHAS_MZ=ALPHAS(XMZ)
C     PRINT*,' XMZ, ALPHAS_MZ =',XMZ, ALPHAS_MZ 
      GAMZ_NLO_R=TEMP*(3.d0+3.d0*(XLE**2+XRE**2)+
     >     (6.d0*(XLU**2+XRU**2)+9.*(XLD**2+XRD**2))*
     >     (1.0+ALPHAS_MZ/PI))
CCPY Aug 2016:      GAMZ_MZ=2.4960*(XMZ/91.1875)**3     
CCPY PDG 2015:      GAMZ_MZ=2.4952*(XMZ/91.1876)**3     
CCPY      GAMZ_MZ=2.4952*(XMZ/91.1876)**3     
      GAMZ_MZ=GAMZ_IN

      WRITE(*,*) ' ALPHAS_MZ, GAMZ, GAMZ_NLO_R, GAMZ_MZ =',
     >            ALPHAS_MZ, GAMZ, GAMZ_NLO_R, GAMZ_MZ 
CCPY      GAMZ_NLO=GAMZ_MZ
      GAMZ_NLO=GAMZ_NLO_R
  
CCPY Aug 2016:
C      BR_W_EN=1.0/9.0/(1.0+2.0*ALPHAS_MW/PI/3.0)
C      BR_Z_EE=TEMP*(XLE**2+XRE**2)/GAMZ_NLO

      BR_W_EN_NLO_R=1.0/9.0/(1.0+2.0*ALPHAS_MW/PI/3.0)
      BR_Z_EE_NLO_R=TEMP*(XLE**2+XRE**2)/GAMZ_NLO_R

      BR_W_EN=1.0/9.0/(1.0+2.0*ALPHAS_MW/PI/3.0)*(GAMW_NLO_R/GAMW_MW)
      BR_Z_EE=TEMP*(XLE**2+XRE**2)/GAMZ_MZ

C==========================
C_____Higgs width
C      Higgs mass is read from the grid or from the input file
C      IF(HIGGS_MASS.EQ.0)HIGGS_MASS=100.
C      XMHIGGS=HIGGS_MASS
      HIGGS_MASS = XMHIGGS
C      Print*, ' Higgs mass =', XMHIGGS
      XMH2=XMHIGGS**2
      GAMH=0.D0
        
        IF(XMHIGGS.GT.(2.*XMBOT))THEN
          XMQX=RUN_MASS_TOT(XMHIGGS,5)
          H2BB=GFERMI*3.*XMHIGGS*XMQX**2/4./SQRT(2.)/PI*
     1 (1.-4.*XMQX**2/XMHIGGS**2)**1.5
          GAMH=GAMH+H2BB

CCPY
C      Print*, ' Gamma(Higgs)_bb =', GAMH

        END IF
        IF(XMHIGGS.GT.(2.*XMTOP))THEN
          XMQX=RUN_MASS_TOT(XMHIGGS,6)
          H2TT=GFERMI*3.*XMHIGGS*XMQX**2/4./SQRT(2.)/PI*
     1 (1.-4.*XMQX**2/XMHIGGS**2)**1.5
          GAMH=GAMH+H2TT
        END IF

57    CONTINUE


CCPY!****************************************************
CCPY!      IF(XMHIGGS.GT.(2.*XMW))THEN
CCPY!        A_W=4.*XMW**2/XMHIGGS**2
CCPY!        GAMH=GAMH+GFERMI*XMHIGGS**3/32./SQRT(2.)/PI*
CCPY!     1 SQRT(1.-A_W)*(4.-4.*A_W+3.*A_W**2)
CCPY!      END IF
CCPY!      IF(XMHIGGS.GT.(2.*XMZ))THEN
CCPY!        A_W=4.*XMZ**2/XMHIGGS**2
CCPY!        GAMH=GAMH+GFERMI*XMHIGGS**3/64./SQRT(2.)/PI*
CCPY!     1 SQRT(1.-A_W)*(4.-4.*A_W+3.*A_W**2)
CCPY!      END IF
CCPY!      GAMH2=GAMH**2
CCPY!****************************************************

c********************************
CCHEN!CCHEN! ADD IN H-->W* W AND H--> Z* Z
CCHEN!CCHEN! ADD IN EW CORRECTION FROM TOP AND BOTTOM QUARK
CCHEN!C      REF. PLB211,200(1988)
CCHEN!      COEFF=8.0*GFERMI/R2*3.0/(16.0*PI2)
CCHEN!      AA1=XMHIGGS**2/XMTOP**2
CCHEN!      AA2=XMHIGGS**2/XMBOT**2
CCHEN!      DEL12SQ=XMTOP**2-XMBOT**2
CCHEN!      
CCHEN!      
CCHEN!CCHEN!      PRINT*,'XX=3:', FI1(3)
CCHEN!CCHEN!      PRINT*,'XX=5:',FI1(5)
CCHEN!CCHEN!      PRINT*
CCHEN!      
CCHEN!      EWC_HWW=1.0+TERM1+TERM2+TERM3+TERM4+TERM5


      IF(XMHIGGS.GT.(2.*XMW))THEN
        A_W=4.*XMW**2/XMHIGGS**2
        H2WW=GFERMI*XMHIGGS**3/32./SQRT(2.)/PI*
     1 SQRT(1.-A_W)*(4.-4.*A_W+3.*A_W**2)
      ELSEIF(XMHIGGS.GT.XMW .AND. XMHIGGS.LT.(2.*XMW))THEN
       CWW=3.D0*GFERMI**2*XMW**4/16.D0/PI**3
       H2WW=HV(XMW**2/XMHIGGS**2)*CWW*XMHIGGS
      END IF
      GAMH=GAMH+H2WW

      IF(XMHIGGS.GT.(2.*XMZ))THEN
        A_W=4.*XMZ**2/XMHIGGS**2
       H2ZZ = GFERMI*XMHIGGS**3/64./SQRT(2.)/PI*
     1         SQRT(1.-A_W)*(4.-4.*A_W+3.*A_W**2)
      ELSEIF(XMHIGGS.GT.XMZ .AND. XMHIGGS.LT.(2.*XMZ))THEN
       CZZ=GFERMI**2*XMZ**4/64.D0/PI**3
     >       *(7-40/3.D0*SWS+160/9.D0*SWS**2)
       H2ZZ=HV(XMZ**2/XMHIGGS**2)*CZZ*XMHIGGS
      END IF
      
CCPY JUNE 2009:      GAMH=GAMH+H2ZZ
      GAMH_TMP=GAMH+H2ZZ
      
C==========================
CCHEN!: PRINT OUT BR(H -> ZZ)
c      PRINT*,'BR(H -> ZZ)', H2ZZ/GAMH_TMP
c      PRINT*,'BR(H -> WW)', H2WW/GAMH_TMP
C*********************************************

CCPY USE THE WIDTH FROM HDECAY CODE
      GAMH=HWIDTH(XMHIGGS)
      
CCPY      Print*,' GAMH=GAMH_HDECAY, GAMH_TMP =',GAMH,GAMH_TMP

      GAMH2=GAMH**2

C==========================
C (KM MATRIX ELEMENT)**2
C FLAVOR LABEL:
C UP=1, DN=2, ST=3, CH=4 BT=5 TP=6, GL=0
C UB=-1, DB=-2, SB=-3, CB=-4 BB=-5 TB=-6
C TOP QUARK CHANNEL IS SWITCHED OFF
C THIS ONLY WORKS IN    QCD   CALCULATIONS
CCAO!      KM_L=0.22D0
CCAO!      KM_A=1.0D0
CCAO!      KM_RHO=-0.4D0
CCAO!      KM_ETA=0.22D0

      KM_L=0.224D0
      KM_A=0.83
      KM_ETA=0.35D0
      VUB=35.7D-4      
C     KM_RHO=0.155D0
      KM_RHO=DSQRT((VUB/KM_A/KM_L**3)**2-KM_ETA**2)



c     WRITE(*,*) ' CKM PARAMETERS:    KM_L, KM_A, KM_ETA, KM_RHO = ',
c    >         KM_L, KM_A, KM_ETA, KM_RHO     

      DO I=1,6
       DO J=1,6
        IF(I.EQ.J) THEN
         VKM(I,J)=1.0D0
        ELSE
         VKM(I,J)=0.0D0
        ENDIF
       ENDDO
      ENDDO
CBIN see paper:hep-ph/0406184
c      VKM(1,2)=(1.0D0-KM_L**2/2.0D0-KM_L**4/8.0D0)**2
c      VKM(1,3)=(KM_L)**2
c      VKM(1,5)=(KM_A*KM_L**3)**2*(KM_RHO**2+KM_ETA**2)
c      VKM(4,2)=(-KM_L+0.5*KM_A**2*KM_L**5*
c     >          DSQRT((1.0-2.0*KM_RHO)**2+KM_ETA**2))**2
c      VKM(4,3)=(1.0D0-KM_L**2/2.0D0-KM_L**4/8.0D0*
c     >         (1.0+4.0*KM_A**2))**2
c      VKM(4,5)=(KM_A*KM_L**2)**2

      VKM(1,2)=0.97427**2
      VKM(1,3)=0.22536**2
      VKM(1,5)=0.00355**2
      VKM(4,2)=0.22522**2
      VKM(4,3)=0.97343**2
      VKM(4,5)=0.0414**2
      VKM(2,1)=VKM(1,2)
      VKM(3,1)=VKM(1,3)
      VKM(5,1)=VKM(1,5)
      VKM(2,4)=VKM(4,2)
      VKM(3,4)=VKM(4,3)
      VKM(5,4)=VKM(4,5)

      IF(FIRST) THEN
        WRITE(16,*) 'HIGGS_MASS, GAMH, XMTOP,GAMT,YM12MIN,YM12MAX'
        WRITE(16,*) HIGGS_MASS, GAMH, XMTOP,GAMT,YM12MIN,YM12MAX
        WRITE(16,*) 'XMW,XMZ,GAMW,GAMZ,XMBOT,XMTOP,SWS,CWS,AEM,GWEAK2'
        WRITE(16,*) XMW,XMZ,GAMW,GAMZ,XMBOT,XMTOP,SWS,CWS,AEM,GWEAK2

        WRITE(16,*) 'VKM(1,2),VKM(1,3),VKM(1,5),SUM_1'
        WRITE(16,*) VKM(1,2),VKM(1,3),VKM(1,5), 
     >              VKM(1,2)+VKM(1,3)+VKM(1,5)
        WRITE(16,*) 'VKM(4,2),VKM(4,3),VKM(4,5),SUM_4'
        WRITE(16,*) VKM(4,2),VKM(4,3),VKM(4,5), 
     >              VKM(4,2)+VKM(4,3)+VKM(4,5)
CCAO!
        WRITE(16,*) 'GFermi,Gmu =',GFERMI,GMU
        WRITE(16,*) 'WCOUPL,ZCOUPL,HCOUPL,ACOUPL'
        WRITE(16,*)  WCOUPL,ZCOUPL,HCOUPL,ACOUPL
        WRITE(16,*)'1/ALFA_EM_MZ  (ALPHA_{EM} AT MZ) =',1.0/ALFA_EM_MZ
        WRITE(16,*) ' SW2_EFF =',SW2_EFF
        WRITE(16,*) ' ALPHAS_MW=', ALPHAS_MW
        WRITE(16,*) ' ALPHAS_MZ=', ALPHAS_MZ
        WRITE(16,*) 'XMW, XMZ =',XMW, XMZ
        WRITE(16,*) 'GAMW_NLO_R, GAMZ_NLO_R =',GAMW_NLO_R, GAMZ_NLO_R 
        WRITE(16,*) 'GAMW_MW, GAMZ_MZ =',GAMW_MW, GAMZ_MZ 
        WRITE(16,*) 'GAMW_NLO, GAMZ_NLO =',GAMW_NLO, GAMZ_NLO 
        WRITE(16,*) 'XMW_IN,GAMW_IN =',XMW_IN,GAMW_IN
        WRITE(16,*) 'XMZ_IN,GAMZ_IN =',XMZ_IN,GAMZ_IN
        WRITE(16,*) 'BR_Z_EE_NLO_R =', BR_Z_EE_NLO_R
        WRITE(16,*) 'BR_W_EN_NLO_R =', BR_W_EN_NLO_R
        WRITE(16,*) 'BR_Z_EE =', BR_Z_EE
        WRITE(16,*) 'BR_W_EN =', BR_W_EN
        WRITE(16,*) 'XMHIGGS, GAMH =',XMHIGGS, GAMH

CBIN
        WRITE(16,*) ' g_fvf_p_u, g_fvf_m_u =',g_fvf_p_u, g_fvf_m_u
        WRITE(16,*) ' g_fvf_p_d, g_fvf_m_d =',g_fvf_p_d, g_fvf_m_d
        WRITE(16,*) 'XLE, XRE, XLU, XLD, XRU, XRD=',XLE, XRE, XLU,
     &               XLD, XRU, XRD
CCPY!  Watch out the difference between GAMH and GAMH_TMP inside 
CCPY!  SUBROUTINE STANDARD 
CCPY!        WRITE(16,*) 'BR(H -> BB) =', H2BB/GAMH
CCPY!        WRITE(16,*) 'BR(H -> TT) =', H2TT/GAMH
CCPY!        WRITE(16,*) 'BR(H -> ZZ) =', H2ZZ/GAMH
CCPY!        WRITE(16,*) 'BR(H -> WW) =', H2WW/GAMH

        WRITE(16,*) ' GAMH=GAMH_HDECAY, GAMH_TMP =',GAMH,GAMH_TMP
        FIRST=.FALSE.
      ENDIF

CCPY FOR CLARITY, I CALL AN EXTERNAL SUBROUTINE FOR 'DY' PROCESS
      CALL DY_COUPL_old

      RETURN
      END ! Standard

C===========================
CCPY SET UP COUPLINGS FOR 'DY' 
C     IF(BOSON.EQ.'ZU') THEN
C THIS IS FOR UP-TYPE QUARKS IN THE INITIAL STATE
C       IFLAV=1
C       BOSON='DY'
C     ELSEIF(BOSON.EQ.'ZD') THEN
C THIS IS FOR DOWN-TYPE QUARKS IN THE INITIAL STATE
C       IFLAV=2
C       BOSON='DY'
C     ENDIF
C============================
      SUBROUTINE DY_COUPL_old
C-----------------------------------------------------------
      IMPLICIT NONE

      INTEGER IFLAV
      REAL*8 CAA_E,CAA_O,CZA_E,CZA_O,CZZ_E,CZZ_O
      COMMON/DY_IFLAV/ CAA_E,CAA_O,CZA_E,CZA_O,CZZ_E,CZZ_O,IFLAV
              
      REAL*8 XMW,XMZ,XMW2,XMZ2,AEM,GAMW,GAMZ,GAMW2,GAMZ2,
     &SWS,CWS,XMBOT,XMTOP,XMB2,XMT2,GWEAK2
      COMMON/STAND0/XMW,XMZ,XMW2,XMZ2,AEM,GAMW,GAMZ,GAMW2,GAMZ2
      COMMON/STAND1/SWS,CWS,XMBOT,XMTOP,XMB2,XMT2,GWEAK2
CCAO!
      REAL*8 GAMW_NLO,GAMZ_NLO,ALFA_EM_MZ,ACOUPL,WCOUPL,
     &ZCOUPL,HCOUPL,SW2_EFF,EW_ALFA
      COMMON/STAND7/ GAMW_NLO,GAMZ_NLO,ALFA_EM_MZ,ACOUPL,WCOUPL,
     >       ZCOUPL,HCOUPL,SW2_EFF,EW_ALFA

      REAL*8 EL_Z,ER_Z,QL_Z,QR_Z,EL_A,ER_A,QL_A,QR_A
      COMMON/ZFF_INIT/
     >   EL_Z,ER_Z,QL_Z(2),QR_Z(2),EL_A,ER_A,QL_A(2),QR_A(2)

      REAL*8 SW2_EFF_u,SW2_EFF_d

      REAL*8 ALPHAS
      EXTERNAL ALPHAS

      REAL*8 FAC_A,FAC
      
      LOGICAL TESTING

C ACOUPL IS (e**2)/2
C      ACOUPL=R2*MW2*GMU*2.0D0*SW2
C ZCOUPL IS (g/Cos_w)**2/4
C      ZCOUPL=R2*MZ2*GMU*(1.0+EW_DELRHO)

C COUPLINGS OF A0-F-F.  (g*Sin_w)*(\gamma_\mu)*(electric charge)
      FAC_A=DSQRT(2.0*ACOUPL)/2.d0
      QL_A(1) = (2.0/3.0)*FAC_A
      QR_A(1) = (2.0/3.0)*FAC_A
      QL_A(2) = (-1.0/3.0)*FAC_A
      QR_A(2) = (-1.0/3.0)*FAC_A
      EL_A = -1.0*FAC_A
      ER_A = -1.0*FAC_A

C COUPLINGS OF Z0-F-F.  -(g/2/Cos_w)*(\gamma_\mu)*(V^I-A^I*\gamma_5)
C  WHICH IS DEFINED AS  -(g/2/Cos_w)*(\gamma_\mu)*(FL*(1-\GAMMA_5) + FR*(1+\GAMMA_5))
      FAC = DSQRT(ZCOUPL)
CCPY!      QL_Z(1) = (0.5-2.0/3.0*SW2_EFF)*FAC
CCPY!      QR_Z(1) = -2.0/3.0*SW2_EFF*FAC
CCPY!      QL_Z(2) = (-0.5+SW2_EFF/3.0)*FAC
CCPY!      QR_Z(2) = (SW2_EFF/3.0)*FAC
      EL_Z = (-0.5+SW2_EFF)*FAC
      ER_Z = SW2_EFF*FAC
CCPY June 2010: Add the difference in SW2_EFF_u and SW2_EFF_d
C For up-type quark
      SW2_EFF_u=SW2_EFF-0.0001
      SW2_EFF_d=SW2_EFF-0.0002
      
      IF(IFLAV.EQ.1) THEN
C THIS IS FOR UP-TYPE QUARKS IN THE INITIAL STATE
        QL_Z(1) = (0.5-2.0/3.0*SW2_EFF_u)*FAC
        QR_Z(1) = -2.0/3.0*SW2_EFF_u*FAC
        QL_Z(2) = (-0.5+SW2_EFF_u/3.0)*FAC
        QR_Z(2) = (SW2_EFF_u/3.0)*FAC

        CAA_E = ( QL_A(1)*QL_A(1) + QR_A(1)*QR_A(1) )*
     >          ( EL_A*EL_A + ER_A*ER_A )
        CAA_O = ( QL_A(1)*QL_A(1) - QR_A(1)*QR_A(1) )*
     >          ( EL_A*EL_A - ER_A*ER_A )

        CZZ_E = ( QL_Z(1)*QL_Z(1) + QR_Z(1)*QR_Z(1) )*
     >          ( EL_Z*EL_Z + ER_Z*ER_Z )
        CZZ_O = ( QL_Z(1)*QL_Z(1) - QR_Z(1)*QR_Z(1) )*
     >          ( EL_Z*EL_Z - ER_Z*ER_Z )

        CZA_E = ( QL_Z(1)*QL_A(1) + QR_Z(1)*QR_A(1) )*
     >          ( EL_Z*EL_A + ER_Z*ER_A )
        CZA_O = ( QL_Z(1)*QL_A(1) - QR_Z(1)*QR_A(1) )*
     >          ( EL_Z*EL_A - ER_Z*ER_A )

      ELSE
C THIS IS FOR DOWN-TYPE QUARKS IN THE INITIAL STATE
        QL_Z(1) = (0.5-2.0/3.0*SW2_EFF_d)*FAC
        QR_Z(1) = -2.0/3.0*SW2_EFF_d*FAC
        QL_Z(2) = (-0.5+SW2_EFF_d/3.0)*FAC
        QR_Z(2) = (SW2_EFF_d/3.0)*FAC

        CAA_E = ( QL_A(2)*QL_A(2) + QR_A(2)*QR_A(2) )*
     >          ( EL_A*EL_A + ER_A*ER_A )
        CAA_O = ( QL_A(2)*QL_A(2) - QR_A(2)*QR_A(2) )*
     >          ( EL_A*EL_A - ER_A*ER_A )

        CZZ_E = ( QL_Z(2)*QL_Z(2) + QR_Z(2)*QR_Z(2) )*
     >          ( EL_Z*EL_Z + ER_Z*ER_Z )
        CZZ_O = ( QL_Z(2)*QL_Z(2) - QR_Z(2)*QR_Z(2) )*
     >          ( EL_Z*EL_Z - ER_Z*ER_Z )

        CZA_E = ( QL_Z(2)*QL_A(2) + QR_Z(2)*QR_A(2) )*
     >          ( EL_Z*EL_A + ER_Z*ER_A )
        CZA_O = ( QL_Z(2)*QL_A(2) - QR_Z(2)*QR_A(2) )*
     >          ( EL_Z*EL_A - ER_Z*ER_A )


      ENDIF

      WRITE(16,*) 'EL_Z, ER_Z, QL_Z1, QR_Z1, QL_Z2, QR_Z2'
      WRITE(16,*)EL_Z,ER_Z, QL_Z(1),QR_Z(1),QL_Z(2),QR_Z(2)

      WRITE(16,*) ' CZZ_E,CZZ_O,CZA_E,CZA_O,CAA_E,CAA_O,IFLAV'
      WRITE(16,*) CZZ_E,CZZ_O,CZA_E,CZA_O,CAA_E,CAA_O,IFLAV

      RETURN
      END

CBIN Sep, 2016
C --------------------------------------------------------------------------
      SUBROUTINE STANDARD
C --------------------------------------------------------------------------
CsB   Updated: 2004 Apr. 20
      IMPLICIT NONE
      INTEGER*4 I,J
      REAL*8 TEMP,A_W,XMQX

      REAL*8 HIGGS_MASS,
     &ZUP,ZDOWN,XMGW2,XMGZ2,XMGH2,XMH2,VSQ,
     &XMW,XMZ,XMW2,XMZ2,AEM,GAMW,GAMZ,GAMW2,GAMZ2,
     &SWS,CWS,XMBOT,XMTOP,XMB2,XMT2,GWEAK2,
     &XLU,XLD,XRU,XRD,QU,QD,GFERMI,XLE,XRE,
     &XLU2,XLD2,XRU2,XRD2,QU2,QD2,XLE2,XRE2,
     &XLU4,XLD4,XRU4,XRD4,QU4,QD4,PI,
     &XMCHARM,XMCHARM2,XMHIGGS,GAMH,GAMH2,
     &CURRENTW_MASS,CURRENTZ_MASS,
     &GAMT,
     &g_fVf_P_u,g_fVf_M_u,g_fVf_P_d,g_fVf_M_d

      COMMON/HIGGS/HIGGS_MASS
      COMMON/HIGGS_STUFF/ZUP,ZDOWN,XMGW2,XMGZ2,XMGH2,XMH2,VSQ
      COMMON/STAND0/XMW,XMZ,XMW2,XMZ2,AEM,GAMW,GAMZ,GAMW2,GAMZ2
      COMMON/STAND1/SWS,CWS,XMBOT,XMTOP,XMB2,XMT2,GWEAK2
      COMMON/STAND2/XLU,XLD,XRU,XRD,QU,QD,GFERMI,XLE,XRE
      COMMON/STAND3/XLU2,XLD2,XRU2,XRD2,QU2,QD2,XLE2,XRE2
      COMMON/STAND4/XLU4,XLD4,XRU4,XRD4,QU4,QD4,PI
      COMMON/STAND5/XMCHARM,XMCHARM2,XMHIGGS,GAMH,GAMH2
      Common/Stand6/g_fVf_P_u,g_fVf_M_u,g_fVf_P_d,g_fVf_M_d
      COMMON/RUNNING_BOSON/CURRENTW_MASS,CURRENTZ_MASS
      COMMON/HTOP/GAMT
      REAL*8 KM_L,KM_A,KM_RHO,KM_ETA,VKM
      COMMON/CKM/VKM(6,6)

      REAL*8 YM12MIN,YM12MAX,QTMin,QTMax, TMmin, Tmmax, ETmin
      COMMON/MYYM12/YM12MIN,YM12MAX,QTMin,QTMax, TMmin, Tmmax, ETmin

CBIN
      REAL*8 YM12,YPT12,YP12
      COMMON/MY1/YM12,YPT12,YP12(4)
CCAO!
      REAL*8 GAMW_NLO,GAMZ_NLO,ALFA_EM_MZ,ACOUPL,WCOUPL,
     &ZCOUPL,HCOUPL,SW2_EFF,EW_ALFA
      COMMON/STAND7/ GAMW_NLO,GAMZ_NLO,ALFA_EM_MZ,ACOUPL,WCOUPL,
     >       ZCOUPL,HCOUPL,SW2_EFF,EW_ALFA

CCAO! THE COUPLING OF Z-F-F
      REAL*8
     >   EL_Z,ER_Z,QL_Z,QR_Z,EL_A,ER_A,QL_A,QR_A
      COMMON/ZFF_INIT/
     >   EL_Z,ER_Z,QL_Z(2),QR_Z(2),EL_A,ER_A,QL_A(2),QR_A(2)
CCAO! The readin W-boson mass and width
      Real*8 XMW_IN,XMW_GRID,GAMW_IN
      COMMON/WMASS_IN/XMW_IN,XMW_GRID,GAMW_IN
      Real*8 XMZ_IN,XMZ_GRID,GAMZ_IN
      COMMON/ZMASS_IN/XMZ_IN,XMZ_GRID,GAMZ_IN
     
      Integer BIPDF, ISet
      Common / BuildInPDF / BIPDF, ISet
      
      REAL*8 Y_EW,Y_PW
      REAL*8 TOP_MASS
      Real*8 GMu
      Real*8 Pi2,R2,
     &  EW_DELA,EW_DELRHO,ALPI_MT,EW_RHO,EW_DELKA,EW_DELR,
     &  Sw_Eff,Cw_Eff,Tw_Eff,T3u,T3d

      REAL*8 GAMW_MW,GAMZ_MZ
      REAL*8 GAMW_NLO_R,GAMZ_NLO_R
      REAL*8 BR_Z_EE_NLO_R,BR_W_EN_NLO_R
      REAL*8 BR_Z_EE,BR_W_EN
      LOGICAL FIRST
      DATA FIRST/.TRUE./
CCAO! 
      REAL*8 UNIT,HBARC2      
      REAL*8 GWEAK,EW_AMZ,EW_A,AEM_MZ
      REAL*8 ZCOUPL_OLD,ZCOUPL_MZ
      REAL*8 ALPHAS_MW,ALPHAS_MZ
      REAL*8 QE,VUB
      REAL*8 SW,CW,SW2_EFF_OLD,SW2_EFF_MZ
      REAL*8 A02

      REAL*8 ALPHAS,ALFAS
      EXTERNAL ALPHAS

      REAL*8 H2BB,H2TT,RUN_MASS_TOT
      EXTERNAL RUN_MASS_TOT
      
      REAL*8 GAMH_TMP, HWIDTH

      REAL*8 PYALEM
      EXTERNAL PYALEM

CBIN
      REAL*8 rho_e,rho_u,rho_d
      REAL*8 ZCOUPL_e,ZCOUPL_u,ZCOUPL_d
      REAL*8 sw2_eff_e,sw2_eff_u,sw2_eff_d,
     >      cw2_eff_e,cw2_eff_u,cw2_eff_d
      COMMON/ZFITTER/rho_e,rho_u,rho_d,
     >      ZCOUPL_e,ZCOUPL_u,ZCOUPL_d,
     >     sw2_eff_e,sw2_eff_u,sw2_eff_d,
     >     cw2_eff_e,cw2_eff_u,cw2_eff_d   
CBIN 
C contral different electroweak parameters
      InTEGER*4 EW_FLAG,qed_FLAG,CKM_FLAG  
      COMMON/FLAG/EW_FLAG,qed_FLAG,CKM_FLAG   

CFEWZ0
      REAL*8 gvu,gau,gvd,gad,gvl,gal
      common/FEWZ0/ gvu,gau,gvd,gad,gvl,gal
      REAL*8 ALFAFE_W,SW2FE_W  
CFEWZ1
CBin =====================================
      complex*8 XMZ2FE,XMW2FE,CWSFE,SWSFE,SWFE,CWFE,TWFE,
     &  g_fVf_P_uFE,g_fVf_M_uFE,g_fVf_P_dFE,g_fVf_M_dFE,
     &  XLEFE,XREFE,XLUFE,XLDFE,XRUFE,XRDFE,ZCOUPLFE
      real*8 XLE2FE,XRE2FE,XLU2FE,XLD2FE,XRU2FE,XRD2FE,XLU4FE,
     &  XLD4FE,XRU4FE,XRD4FE,ALFA_EM_MZFE,ACOUPLFE,WCOUPLFE
      common/FEWZ/SWSFE
      
CCHEN!********************** 
      REAL*8 CWW,H2WW,CZZ,H2ZZ,HV,V
      REAL*8 EWC_HWW,COEFF,AA1,AA2,FI1,XX,
     .       TERM1,TERM2,TERM3,TERM4,TERM5,DEL12SQ
      EXTERNAL FI1
C_yfu add for changing STW
      CHARACTER*40 envsin2w

      HV(V)=3.D0*(1.D0-8.D0*V+20.D0*V**2)/DSQRT((4.D0*V-1.D0))
     .      *DACOS((3.D0*V-1.D0)/2.D0/DSQRT(V**3))
     .      -(1.D0-V)*(47.D0/2.D0*V-13.D0/2.D0+1.D0/V)
     .      -3.D0/2.D0*(1.D0-6.D0*V+4.D0*V**2)*DLOG(V)
c********************************

C EW_FLAG=3, all the parameters are read from FEWZ3.1 input file
        IF(EW_FLAG.eq.3) THEN
C          PRINT*,' XMW, XMZ, GAMW_IN, GAMZ_IN ARE ALL FIXED 
C     >TO BE THE SAME AS FEWZ3.1 INPUT'
           XMW=80.403d0
           XMZ=91.1876d0
           GAMW_IN=2.141d0
      	    GAMZ_IN=2.4952d0
         ENDIF

C      print*,' inside standard: XMW, XMZ =',XMW,XMZ

CCAO!      Data PI/3.14159 26535 89793 24D0/
      PI=4.0D0*ATAN(1.0D0)      

C_____Parameters
      Pi2 = Pi**2
      R2 = Sqrt(2.d0)
CCAO!
C UNIT: to convert energy from  GeV  to TeV, then unit=1.0d-3
c WARNING!!! You must keep UNIT=1.0D0 in this version of the code!!!
      UNIT=1.0D0
C Convert from GeV to pb
      HBARC2=UNIT**2*0.38937966D0*1.0D9
      
CsB   EW couplings
      GMU = 1.16637D-5/UNIT**2
      EW_ALFA = 1.0D0/137.0359895D0

CsB___Quark masses
C      XMCHARM=1.5D0
      XMCHARM=1.4D0
      XMCHARM2=XMCHARM**2
CCPY      XMBOT=5.0D0
C      XMBOT=4.5D0
      XMBOT=4.75D0
      XMB2=XMBOT**2
      TOP_MASS = xMTop
C      Print*, ' TOP_MASS =', TOP_MASS

CsB___W and Z mass
C     The W and Z masses are read from the grid.
      XMZ2=XMZ**2
      XMW2=XMW**2
      CURRENTZ_MASS=XMZ
      CURRENTW_MASS=XMW
C_____Test W, Z and top masses
!     Print*, ' <STANDARD> W, Z masses:', XMW,XMZ

C_____EW_DELA=0.0595 +/- 0.0009
CCAO!      EW_DELA=0.0595+40.0D0/9.0D0*EW_ALFA/PI*DLOG(xMZ/91.187D0)
C          EW_DELA=0.02767+0.03150-0.00007=0.0591
      EW_DELA=0.0591+40.0D0/9.0D0*EW_ALFA/PI*DLOG(XMZ/91.187D0/UNIT)

C_____Rho and Delta rho
C EW_ALFA AT  MZ
      EW_AMZ=EW_ALFA/(1.0D0-EW_DELA)
      EW_A=PI*EW_ALFA/R2/GMU*UNIT**2
      ALPI_MT=ALPHAS(XMTOP)/PI      
C     PRINT*,' XMTOP, ALPHAS(XMTOP)=',XMTOP, ALPI_MT*PI
      EW_DELRHO=3.0D0*GMU*TOP_MASS**2/8.0D0/R2/PI2
      EW_DELRHO=EW_DELRHO*(1.D0-2.0/3.0*(PI2/3.0+1.0)*ALPI_MT)
      EW_RHO=1.0D0/(1.0D0-EW_DELRHO)

CsB___Sin^2 and Cos^2 of the weak angle - in the on shell scheme
      CWS=XMW2/XMZ2
      SWS=1.D0-CWS
      SW=DSQRT(SWS)
      CW=DSQRT(CWS)
CCPY  EW_DELKA=CWS/SWS*EW_DELRHO
C     EW_DELR=EW_DELA-EW_DELKA
C     print*,' old ew_delr,EW_DELA,EW_DELKA =',ew_delr,EW_DELA,EW_DELKA
CPN___Calculate delta_r using mW
      A02 = PI*EW_AMZ/R2/GMU
      EW_DELR = 1 - A02/XMW2/SWS 
      EW_DELKA =  EW_DELA - EW_DELR
C     print*,' new ew_delr,EW_DELKA  =',ew_delr,EW_DELKA 
CCAO! SW2_EFF=SWS*(1.0d0 + EW_DELKA)
      SW2_EFF_OLD=SWS*(1.0d0 + EW_DELKA)

CCPY FIXED A SCHEME HERE
CCPY!      SW2_EFF_MZ=0.23143
c      SW2_EFF_MZ=0.231676     ! Zfitter (mt=172, mH=115)

c      SW2_EFF=SW2_EFF_MZ
c      Sw_Eff = Sqrt(SW2_EFF)
c      Cw_Eff = Sqrt(1.d0 - SW2_EFF)
c      Tw_Eff = Sw_Eff/Cw_Eff

CBIN Sept,2016 Read the Flag for EW and QED couplings
C EW_FLAG=0, original electroweak parameters
C EW_FLAG=1, electroweak parameters are read from ZFITTER
C EW_FLAG=2, electroweak parameters are used the effective born approximation formula
c            based on the original Resbos code
C EW_FLAG=3, electroweak parameters from FEWZ0 scheme
c EW_FLAG=4, electroweak parameters from FEWZ1 scheme (Gmu-scheme, complex mass)
c Note EW_FLAG=3,4 qed coupling is also fixed.

C qed_FLAG=0, original electroweak parametes
C qed_FLAG=1, qed coupling alpha uses running coupling
c qed_FLAG=2, qed coupling alpha uses G_mu scheme formula

C CKM_FLAG=0, original CKM parameters
C CKM_FLAG=1, 2016 PDG number
C Note the grids used CKM_FLAG=0 parameters


      IF(EW_FLAG.eq.0) then
         rho_e=1.005245
         rho_u=1.005940
         rho_d=1.006907
         sw2_eff_e=0.231676  ! Zfitter (mt=172,mH=115)
         sw2_eff_u=0.231676
         sw2_eff_d=0.231676
         SW2_EFF_u=SW2_EFF_u-0.0001
         SW2_EFF_d=SW2_EFF_d-0.0002
      ELSEIF(EW_FLAG.eq.1) then
C For pole mass m_t = 174 GeV (PDG, 2016)
         rho_e=1.005245
         rho_u=1.005940
         rho_d=1.006907
         sw2_eff_e=0.231514
         sw2_eff_u=0.231407
         sw2_eff_d=0.231280
      ELSEIF(EW_FLAG.eq.2) then
         rho_e=(1.0+EW_DELRHO)**2
         rho_u=(1.0+EW_DELRHO)**2
         rho_d=(1.0+EW_DELRHO)**2
C_yfu add for changing STW
C         sw2_eff_e=0.231676
C         sw2_eff_u=0.231676
C         sw2_eff_d=0.231676
         call getenv('SINT2W',envsin2w)
         if ( envsin2w(1:1).eq.' ') then
           sw2_eff_e=0.2315
           sw2_eff_u=0.2315
           sw2_eff_d=0.2315
         else
           read (envsin2w,*) sw2_eff_e
           read (envsin2w,*) sw2_eff_u
           read (envsin2w,*) sw2_eff_d
         endif
         SW2_EFF_u=SW2_EFF_u-0.0001
         SW2_EFF_d=SW2_EFF_d-0.0002
       ENDIF

       IF(EW_FLAG.eq.3) then
C        PRINT*,' ALL THE COUPLINGS ARE FIXED TO BE THE SAME AS FEWZ3.1 INPUT'
       ENDIF

       IF(EW_FLAG.eq.4) then
        XMZ2FE=XMZ**2-dcmplx(0d0,1d0)*GAMZ_MZ*XMZ
        XMW2FE=XMW**2-dcmplx(0d0,1d0)*GAMW_MW*XMW
        CWSFE=XMW2FE/XMZ2FE
        SWSFE=1.D0-CWSFE
        SWFE=SQRT(Real(SWSFE))
        CWFE=SQRT(Real(CWSFE))
        TWFE=SWFE/CWFE
   
        CWS=REAL(CWSFE)
        SWS=REAL(SWSFE)
        SW=REAL(SWFE)
        CW=REAL(CWFE)
       ENDIF

      cw2_eff_e=1-sw2_eff_e
      cw2_eff_u=1-sw2_eff_u
      cw2_eff_d=1-sw2_eff_d
C_____alpha = e^2/(4 Pi)
C      AEM=1.D0/128.0
      AEM = Sqrt(2.d0)*XMW2*GMU*SWS/Pi
      AEM_MZ=1.0/128.937
      ALFA_EM_MZ=EW_AMZ
C_____Test alpha_em
      WRITE(16,*)'1/AEM, 1/EW_AMZ, 1/AEM_MZ =',1/AEM,1/EW_AMZ,1/AEM_MZ

      WRITE(16,*)'1/AEM_MZp=',1/PYALEM(XMZ2)

    

C_____Weak couplings: g_w and Fermi constant
      GWEAK2=AEM*4.d0*PI/SWS
      GWEAK = Sqrt(GWEAK2)
      GFERMI=PI*AEM/SQRT(2.d0)/SWS/XMW2

C_____Quark charges
      Qu  =  2.d0/3.d0
      Qd  = -1.d0/3.d0
      QU2=QU**2
      QU4=QU2**2
      QD2=QD**2
      QD4=QD2**2

C_____Quark weak-isospins
      T3u =  1.d0/2.d0
      T3d = -1.d0/2.d0

C_____Couplings for fermion-Z-fermion vertex in PRD43('91)3626 Eq.(3)
c      g_fVf_P_u = - Qu*Tw_Eff                   ! for f=u type quark & V = Z
c      g_fVf_M_u = T3u/Sw_Eff/Cw_Eff - Qu*Tw_Eff ! for f=u type quark & V = Z
c      g_fVf_P_d = - Qd*Tw_Eff                   ! for f=d type quark & V = Z
c      g_fVf_M_d = T3d/Sw_Eff/Cw_Eff - Qd*Tw_Eff ! for f=d type quark & V = Z

      IF(EW_FLAG.lE.2) THEN
         g_fVf_P_u = - Qu*Sqrt(sw2_eff_e/cw2_eff_e)                   ! for f=u type quark & V = Z
         g_fVf_M_u = T3u/Sqrt(sw2_eff_u)/Sqrt(cw2_eff_u) -
     >            Qu*Sqrt(sw2_eff_u/cw2_eff_u)                     ! for f=u type quark & V = Z
         g_fVf_P_d = - Qd*Sqrt(sw2_eff_d/cw2_eff_d)                      ! for f=d type quark & V = Z
         g_fVf_M_d = T3d/Sqrt(sw2_eff_d)/Sqrt(cw2_eff_d) -
     >            Qd*Sqrt(sw2_eff_d/cw2_eff_d) 


C_____Lepton and quark couplings
CCPY: ERROR      XLE=1.D0+2.*SWS
      QE=-1.0
c      XLE=-1.D0-2.*SW2_EFF*QE
c      XRE=-2.D0*QE*SW2_EFF
c      XLU=1.D0-2.*SW2_EFF*QU
c      XLD=-1.D0-2.*SW2_EFF*QD
c      XRU=-2.*QU*SW2_EFF
c      XRD=-2.*QD*SW2_EFF



         XLE=-1.D0-2.*SW2_EFF_e*QE
         XRE=-2.D0*QE*SW2_EFF_e
         XLU=1.D0-2.*SW2_EFF_u*QU
         XLD=-1.D0-2.*SW2_EFF_d*QD
         XRU=-2.*QU*SW2_EFF_u
         XRD=-2.*QD*SW2_EFF_d
         XLE2=XLE**2
         XRE2=XRE**2
         XLU2=XLU**2
         XLD2=XLD**2
         XRU2=XRU**2
         XRD2=XRD**2
         XLU4=XLU**4
         XLD4=XLD**4
         XRU4=XRU**4
         XRD4=XRD**4
C FEWZ-0
      ELSEIF(EW_FLAG.EQ.3) THEN 
c==============================================
C EW_FLAG=3, all the parameters are read from FEWZ3.1 input file
cBin the default input value of FEWZ3.1 (z boson)
c g/4cw[...]
         gvu=0.4091d0
         gvd=-0.7045d0
         gvl=-0.11360d0
         gau=-1d0
         gad=1.0d0
         gal=1.0d0
c Resbos renormalization factor is
c g/2cw[...]
         gvu=gvu/2.
         gvd=gvd/2.
         gvl=gvl/2.
         gau=gau/2.
         gad=gad/2.
         gal=gal/2.
        SW2_eff=0.22255d0
        SW_EFF=Sqrt(SW2_eff)
        CW_EFF=Sqrt(1-SW2_eff)
c fL=gv-ga,  fr=gv+ga 
        g_fVf_P_u = (gvu+gau)/(2.0*SW_EFF*CW_EFF)                 ! for f=u type quark & V = Z
        g_fVf_M_u = (gvu-gau)/(2.0*SW_EFF*CW_EFF)! for f=u type quark & V = Z
        g_fVf_P_d = (gvd+gad)/(2.0*SW_EFF*CW_EFF)                  ! for f=d type quark & V = Z
        g_fVf_M_d = (gvd-gad)/(2.0*SW_EFF*CW_EFF) ! for f=d type quark & V = Z

        XLE=gvl-gal
        XRE=gvl+gal
        XLU=gvu-gau
        XLD=gvd-gad
        XRU=gvu+gau
        XRD=gvd+gad

        XLE2=XLE**2
        XRE2=XRE**2
        XLU2=XLU**2
        XLD2=XLD**2
        XRU2=XRU**2
        XRD2=XRD**2
        XLU4=XLU**4
        XLD4=XLD**4
        XRU4=XRU**4
        XRD4=XRD**4
      ENDIF

      IF(EW_FLAG.eq.4) then
   
      g_fVf_P_uFE = - Qu*TWFE                   ! for f=u type quark & V = Z
      g_fVf_M_uFE = T3u/SWFE/CWFE - Qu*TWFE! for f=u type quark & V = Z
      g_fVf_P_dFE = - Qd*TWFE                   ! for f=d type quark & V = Z
      g_fVf_M_dFE = T3d/SWFE/CWFE - Qd*TWFE ! for f=d type quark & V = Z

      g_fVf_P_u = REAL(g_fVf_P_uFE)                   ! for f=u type quark & V = Z
      g_fVf_M_u = REAL(g_fVf_M_uFE)                  !for f=u type quark & V = Z
      g_fVf_P_d = REAL(g_fVf_P_dFE)                  ! for f=d type quark & V = Z
      g_fVf_M_d = REAL(g_fVf_M_dFE)                 ! for f=d type quark & V = Z

      XLEFE=-1.D0-2.*SWSFE*QE
      XREFE=-2.D0*QE*SWSFE
      XLUFE=1.D0-2.*SWSFE*QU
      XLDFE=-1.D0-2.*SWSFE*QD
      XRUFE=-2.*QU*SWSFE
      XRDFE=-2.*QD*SWSFE
      XLE2FE=dble(XLEFE*conjg(XLEFE))
      XRE2FE=dble(XREFE*conjg(XREFE))
      XLU2FE=dble(XLUFE*conjg(XLUFE))
      XLD2FE=dble(XLDFE*conjg(XLDFE))
      XRU2FE=dble(XRUFE*conjg(XRUFE))
      XRD2FE=dble(XRDFE*conjg(XRDFE))

      XLU4FE=dble(XLUFE*conjg(XLUFE))**2
      XLD4FE=dble(XLDFE*conjg(XLUFE))**2
      XRU4FE=dble(XRUFE*conjg(XRUFE))**2
      XRD4FE=dble(XRDFE*conjg(XRDFE))**2


      XLE=REAL(XLEFE)
      XRE=REAL(XREFE)
      XLU=REAL(XLUFE)
      XLD=REAL(XLDFE)
      XRU=REAL(XRUFE)
      XRD=REAL(XRDFE)
      XLE2=XLE2FE
      XRE2=XRE2FE
      XLU2= XLU2FE
      XLD2=XLD2FE
      XRU2=XRU2FE
      XRD2=XRD2FE
      XLU4=XLU4FE
      XLD4=XLD4FE
      XRU4=XRU4FE
      XRD4=XRD4FE
      ENDIF

C_____Test top mass - it is read from the grid or from the input file
C      Print*, ' <STANDARD> top mass:', XMTOP,TOP_MASS
      IF(XMTOP.LT.XMBOT) XMTOP=TOP_MASS
      XMT2=XMTOP**2
C_____Calculate top width
      IF(XMTOP.GT.XMW+XMBOT) THEN
        Y_EW=(XMTOP**2+XMW**2-XMBOT**2)/2.0/XMTOP
        Y_PW=SQRT(Y_EW**2-XMW**2)
        GAMT=GFERMI*Y_PW/2./PI/SQRT(2.)*
     >      (2.0*Y_PW**2+3.0*XMW**2*(1.0-Y_EW/XMTOP))
      ELSE
        GAMT=3.0*9.0*XMTOP**5*GFERMI**2/(192.0*PI**3)
      ENDIF

C     ACOUPL IS (e**2)/2
c      ACOUPL=R2*XMW2*GMU*2.0D0*SWS
C     WCOUPL IS (g**2)/8
c      WCOUPL=R2*XMW2*GMU/2.0D0
C     WE IGNORE THE DIFFERENCE BETWEEN b-quark AND OTHER FERMIONS
C     ZCOUPL IS (g/Cos_w)**2/4
      ZCOUPL_OLD=R2*XMZ2*GMU*(1.0+EW_DELRHO)
      ZCOUPL_MZ=PI*ALFA_EM_MZ/SWS/CWS
C_yfu open up the ZCOUPL
      ZCOUPL=ZCOUPL_MZ 

CBIN Zcoupling
C EW_FLAG=0, original electroweak parameters
C EW_FLAG=1, electroweak parameters are read from ZFITTER
C EW_FLAG=2, electroweak parameters are used the effective born approximation formula
c            based on the original Resbos code
C qed_FLAG=0, original electroweak parametes
C qed_FLAG=1, qed coupling alpha uses running coupling 
c qed_FLAG=2, qed coupling alpha uses G_mu scheme formula
      IF(EW_FLAG.eq.0) then 
         ZCOUPL_e=ZCOUPL_MZ 
         ZCOUPL_u=ZCOUPL_MZ 
         ZCOUPL_d=ZCOUPL_MZ 
      ELSEIF(EW_FLAG.eq.1) then
         ZCOUPL_e=R2*XMZ2*GMU*Sqrt(rho_e)
         ZCOUPL_u=R2*XMZ2*GMU*Sqrt(rho_u)
         ZCOUPL_d=R2*XMZ2*GMU*Sqrt(rho_d)      
      ELSEIF(EW_FLAG.eq.2) then
         ZCOUPL_e=R2*XMZ2*GMU*Sqrt(rho_e)
         ZCOUPL_u=R2*XMZ2*GMU*Sqrt(rho_u)
         ZCOUPL_d=R2*XMZ2*GMU*Sqrt(rho_d)
      ENDIF
CBIN QED coupling
      IF(qed_FLAG.eq.0) then
          ACOUPL=R2*XMW2*GMU*2.0D0*SWS
      ELSEIF(qed_FLAG.eq.1) then
          ACOUPL=2.0*PI*PYALEM(YM12)
      ELSEIF(qed_FLAG.eq.2) then
          ACOUPL=2.0*PI*AEM
      ENDIF

c      ACOUPL=2.0*PI*ALFA_EM_MZ
CBIN running coupling at scale Q^2
c      ACOUPL=2.0*PI*PYALEM(YM12)

CBIN Gmu scheme

c      ACOUPL=2.0*PI*AEM
      WCOUPL=R2*XMW2*GMU/2.0D0
c      ZCOUPL=ZCOUPL_OLD    
c     WRITE(*,*) 'ZCOUPL_OLD, ZCOUPL_MZ, ZCOUPL =',
c    > ZCOUPL_OLD, ZCOUPL_MZ, ZCOUPL
C     HCOUPL IS FOR HIGGS
      HCOUPL=R2*GMU*PI/576.0D0

CNOTE: FE_FLAG=3, all the couplings are fixed, so QED flag is invalid
      IF(EW_FLAG.EQ.3) THEN
         ACOUPL=4.0*pi*AEM_MZ/2.0
         ZCOUPL=ZCOUPL_MZ  
         ALFAFE_W=0.0078125
         SW2FE_W=SW2_EFF
         WCOUPL=ALFAFE_W*PI/2.0/SW2FE_W

  	 ZCOUPL_e=ZCOUPL
	 ZCOUPL_u=ZCOUPL
	 ZCOUPL_d=ZCOUPL
      ENDIF

      IF(EW_FLAG.eq.4) THEN
CBIN
C NOTE: For the W coupling in FEWZ3.1, the sw definition is different with z boson
c production process. All of the parameters are read from the input file
      ALFAFE_W=0.0078125
      SW2FE_W=0.22255
      ALFA_EM_MZFE=DSQRT(2.d0)/PI*GMU*XMW2*(1-XMW2/XMZ2)
      ACOUPLFE=(4.0*PI*ALFA_EM_MZFE)/2.0
      WCOUPLFE=ALFAFE_W*PI/2/SW2FE_W
      ZCOUPLFE=PI*ALFA_EM_MZFE/SWSFE/CWSFE

      ALFA_EM_MZ=ALFA_EM_MZFE
      ACOUPL=ACOUPLFE
      WCOUPL=WCOUPLFE
      ZCOUPL=REAL(ZCOUPLFE) 

      ZCOUPL_e=ZCOUPL
      ZCOUPL_u=ZCOUPL
      ZCOUPL_d=ZCOUPL
      ENDIF


      write(16,*)'Zcouple, Zcouplu, Zcoupld=',
     >           ZCOUPL_e,ZCOUPL_u,ZCOUPL_d

C_____W width
      GAMW=GFERMI*XMW**3/6./PI/SQRT(2.)*(9.)
      GAMW2=GAMW**2
C     NLO QCD CORRECTION
CCPY  ALPHAS_MW=ALPI(XMW)*PI      
      ALPHAS_MW=ALPHAS(XMW)
C     PRINT*,' XMW, ALPHAS_MW =',XMW, ALPHAS_MW 
      GAMW_NLO_R=GAMW*(1.0+2.0*ALPHAS_MW/PI/3.0)

C_____Z width & Br
CCPY Aug 2016      GAMW_MW=2.093*(XMW/80.385)**3      
CCPY PDG 2105      GAMW_MW=2.085*(XMW/80.385)**3      
CCPY      GAMW_MW=2.118*(XMW/80.426)**3      
      GAMW_MW=GAMW_IN

      WRITE(*,*) ' ALPHAS_MW, GAMW, GAMW_NLO_R, GAMW_MW =',
     >            ALPHAS_MW, GAMW, GAMW_NLO_R, GAMW_MW  
CCPY      GAMW_NLO=GAMW_MW
      GAMW_NLO=GAMW_NLO_R

      TEMP=XMZ**3/DSQRT(2.d0)/12.d0/PI*GFERMI
CBIN
      IF(EW_FLAG.eq.0) then 
         TEMP=TEMP
      ELSE
         TEMP=TEMP*(1.0+EW_DELRHO)
      ENDIF

      GAMZ=TEMP*(3.d0+3.d0*(XLE**2+XRE**2)+
     +     6.d0*(XLU**2+XRU**2)+9.*(XLD**2+XRD**2))
      GAMZ2=GAMZ**2     
C     NLO QCD CORRECTION
CCPY  ALPHAS_MZ=ALPI(XMZ)*PI
      ALPHAS_MZ=ALPHAS(XMZ)
C     PRINT*,' XMZ, ALPHAS_MZ =',XMZ, ALPHAS_MZ 
      GAMZ_NLO_R=TEMP*(3.d0+3.d0*(XLE**2+XRE**2)+
     >     (6.d0*(XLU**2+XRU**2)+9.*(XLD**2+XRD**2))*
     >     (1.0+ALPHAS_MZ/PI))
CCPY Aug 2016:      GAMZ_MZ=2.4960*(XMZ/91.1875)**3     
CCPY PDG 2015:      GAMZ_MZ=2.4952*(XMZ/91.1876)**3     
CCPY      GAMZ_MZ=2.4952*(XMZ/91.1876)**3     
      GAMZ_MZ=GAMZ_IN

      WRITE(*,*) ' ALPHAS_MZ, GAMZ, GAMZ_NLO_R, GAMZ_MZ =',
     >            ALPHAS_MZ, GAMZ, GAMZ_NLO_R, GAMZ_MZ 
CCPY      GAMZ_NLO=GAMZ_MZ
      GAMZ_NLO=GAMZ_NLO_R
  
CCPY Aug 2016:
C      BR_W_EN=1.0/9.0/(1.0+2.0*ALPHAS_MW/PI/3.0)
C      BR_Z_EE=TEMP*(XLE**2+XRE**2)/GAMZ_NLO

      BR_W_EN_NLO_R=1.0/9.0/(1.0+2.0*ALPHAS_MW/PI/3.0)
      BR_Z_EE_NLO_R=TEMP*(XLE**2+XRE**2)/GAMZ_NLO_R

      BR_W_EN=1.0/9.0/(1.0+2.0*ALPHAS_MW/PI/3.0)*(GAMW_NLO_R/GAMW_MW)
      BR_Z_EE=TEMP*(XLE**2+XRE**2)/GAMZ_MZ

C==========================
C_____Higgs width
C      Higgs mass is read from the grid or from the input file
C      IF(HIGGS_MASS.EQ.0)HIGGS_MASS=100.
C      XMHIGGS=HIGGS_MASS
      HIGGS_MASS = XMHIGGS
C      Print*, ' Higgs mass =', XMHIGGS
      XMH2=XMHIGGS**2
      GAMH=0.D0
        
        IF(XMHIGGS.GT.(2.*XMBOT))THEN
          XMQX=RUN_MASS_TOT(XMHIGGS,5)
          H2BB=GFERMI*3.*XMHIGGS*XMQX**2/4./SQRT(2.)/PI*
     1 (1.-4.*XMQX**2/XMHIGGS**2)**1.5
          GAMH=GAMH+H2BB

CCPY
C      Print*, ' Gamma(Higgs)_bb =', GAMH

        END IF
        IF(XMHIGGS.GT.(2.*XMTOP))THEN
          XMQX=RUN_MASS_TOT(XMHIGGS,6)
          H2TT=GFERMI*3.*XMHIGGS*XMQX**2/4./SQRT(2.)/PI*
     1 (1.-4.*XMQX**2/XMHIGGS**2)**1.5
          GAMH=GAMH+H2TT
        END IF

57    CONTINUE


CCPY!****************************************************
CCPY!      IF(XMHIGGS.GT.(2.*XMW))THEN
CCPY!        A_W=4.*XMW**2/XMHIGGS**2
CCPY!        GAMH=GAMH+GFERMI*XMHIGGS**3/32./SQRT(2.)/PI*
CCPY!     1 SQRT(1.-A_W)*(4.-4.*A_W+3.*A_W**2)
CCPY!      END IF
CCPY!      IF(XMHIGGS.GT.(2.*XMZ))THEN
CCPY!        A_W=4.*XMZ**2/XMHIGGS**2
CCPY!        GAMH=GAMH+GFERMI*XMHIGGS**3/64./SQRT(2.)/PI*
CCPY!     1 SQRT(1.-A_W)*(4.-4.*A_W+3.*A_W**2)
CCPY!      END IF
CCPY!      GAMH2=GAMH**2
CCPY!****************************************************

c********************************
CCHEN!CCHEN! ADD IN H-->W* W AND H--> Z* Z
CCHEN!CCHEN! ADD IN EW CORRECTION FROM TOP AND BOTTOM QUARK
CCHEN!C      REF. PLB211,200(1988)
CCHEN!      COEFF=8.0*GFERMI/R2*3.0/(16.0*PI2)
CCHEN!      AA1=XMHIGGS**2/XMTOP**2
CCHEN!      AA2=XMHIGGS**2/XMBOT**2
CCHEN!      DEL12SQ=XMTOP**2-XMBOT**2
CCHEN!      
CCHEN!      
CCHEN!CCHEN!      PRINT*,'XX=3:', FI1(3)
CCHEN!CCHEN!      PRINT*,'XX=5:',FI1(5)
CCHEN!CCHEN!      PRINT*
CCHEN!      
CCHEN!      EWC_HWW=1.0+TERM1+TERM2+TERM3+TERM4+TERM5


      IF(XMHIGGS.GT.(2.*XMW))THEN
        A_W=4.*XMW**2/XMHIGGS**2
        H2WW=GFERMI*XMHIGGS**3/32./SQRT(2.)/PI*
     1 SQRT(1.-A_W)*(4.-4.*A_W+3.*A_W**2)
      ELSEIF(XMHIGGS.GT.XMW .AND. XMHIGGS.LT.(2.*XMW))THEN
       CWW=3.D0*GFERMI**2*XMW**4/16.D0/PI**3
       H2WW=HV(XMW**2/XMHIGGS**2)*CWW*XMHIGGS
      END IF
      GAMH=GAMH+H2WW

      IF(XMHIGGS.GT.(2.*XMZ))THEN
        A_W=4.*XMZ**2/XMHIGGS**2
       H2ZZ = GFERMI*XMHIGGS**3/64./SQRT(2.)/PI*
     1         SQRT(1.-A_W)*(4.-4.*A_W+3.*A_W**2)
      ELSEIF(XMHIGGS.GT.XMZ .AND. XMHIGGS.LT.(2.*XMZ))THEN
       CZZ=GFERMI**2*XMZ**4/64.D0/PI**3
     >       *(7-40/3.D0*SWS+160/9.D0*SWS**2)
       H2ZZ=HV(XMZ**2/XMHIGGS**2)*CZZ*XMHIGGS
      END IF
      
CCPY JUNE 2009:      GAMH=GAMH+H2ZZ
      GAMH_TMP=GAMH+H2ZZ
      
C==========================
CCHEN!: PRINT OUT BR(H -> ZZ)
c      PRINT*,'BR(H -> ZZ)', H2ZZ/GAMH_TMP
c      PRINT*,'BR(H -> WW)', H2WW/GAMH_TMP
C*********************************************

CCPY USE THE WIDTH FROM HDECAY CODE
      GAMH=HWIDTH(XMHIGGS)
      
CCPY      Print*,' GAMH=GAMH_HDECAY, GAMH_TMP =',GAMH,GAMH_TMP

      GAMH2=GAMH**2

C==========================
C (KM MATRIX ELEMENT)**2
C FLAVOR LABEL:
C UP=1, DN=2, ST=3, CH=4 BT=5 TP=6, GL=0
C UB=-1, DB=-2, SB=-3, CB=-4 BB=-5 TB=-6
C TOP QUARK CHANNEL IS SWITCHED OFF
C THIS ONLY WORKS IN    QCD   CALCULATIONS
CCAO!      KM_L=0.22D0
CCAO!      KM_A=1.0D0
CCAO!      KM_RHO=-0.4D0
CCAO!      KM_ETA=0.22D0
      KM_L=0.224D0
      KM_A=0.83
      KM_ETA=0.35D0
      VUB=35.7D-4      
C     KM_RHO=0.155D0
      KM_RHO=DSQRT((VUB/KM_A/KM_L**3)**2-KM_ETA**2)
c     WRITE(*,*) ' CKM PARAMETERS:    KM_L, KM_A, KM_ETA, KM_RHO = ',
c    >         KM_L, KM_A, KM_ETA, KM_RHO     

      DO I=1,6
       DO J=1,6
        IF(I.EQ.J) THEN
         VKM(I,J)=1.0D0
        ELSE
         VKM(I,J)=0.0D0
        ENDIF
       ENDDO
      ENDDO
     
c CKM_FLAG=0, default CKM
c CKM_FLAG=1, 2016 PDG
      IF(CKM_FLAG==0) then  
          VKM(1,2)=(1.0D0-KM_L**2/2.0D0-KM_L**4/8.0D0)**2
          VKM(1,3)=(KM_L)**2
          VKM(1,5)=(KM_A*KM_L**3)**2*(KM_RHO**2+KM_ETA**2)
          VKM(4,2)=(-KM_L+0.5*KM_A**2*KM_L**5*
     >          DSQRT((1.0-2.0*KM_RHO)**2+KM_ETA**2))**2
          VKM(4,3)=(1.0D0-KM_L**2/2.0D0-KM_L**4/8.0D0*
     >         (1.0+4.0*KM_A**2))**2
          VKM(4,5)=(KM_A*KM_L**2)**2
       ELSEIF(CKM_FLAG==1) then
          VKM(1,2)=0.97434**2
          VKM(1,3)=0.22506**2
          VKM(1,5)=0.00357**2
          VKM(4,2)=0.22492**2
          VKM(4,3)=0.97351**2
          VKM(4,5)=0.0411**2
	  VKM(6,2)=0.00875**2
	  VKM(6,3)=0.0403**2
          VKM(6,5)=0.99915**2
       ENDIF

      VKM(2,1)=VKM(1,2)
      VKM(3,1)=VKM(1,3)
      VKM(5,1)=VKM(1,5)
      VKM(2,4)=VKM(4,2)
      VKM(3,4)=VKM(4,3)
      VKM(5,4)=VKM(4,5)
      VKM(2,6)=VKM(6,2)
      VKM(3,6)=VKM(6,3)
      VKM(5,6)=VKM(6,5)

      IF(FIRST) THEN
        WRITE(16,*) 'HIGGS_MASS, GAMH, XMTOP,GAMT,YM12MIN,YM12MAX'
        WRITE(16,*) HIGGS_MASS, GAMH, XMTOP,GAMT,YM12MIN,YM12MAX
        WRITE(16,*) 'XMW,XMZ,GAMW,GAMZ,XMBOT,XMTOP,SWS,CWS,AEM,GWEAK2'
        WRITE(16,*) XMW,XMZ,GAMW,GAMZ,XMBOT,XMTOP,SWS,CWS,AEM,GWEAK2

        WRITE(16,*) 'VKM(1,2),VKM(1,3),VKM(1,5),SUM_1'
        WRITE(16,*) VKM(1,2),VKM(1,3),VKM(1,5), 
     >              VKM(1,2)+VKM(1,3)+VKM(1,5)
        WRITE(16,*) 'VKM(4,2),VKM(4,3),VKM(4,5),SUM_4'
        WRITE(16,*) VKM(4,2),VKM(4,3),VKM(4,5), 
     >              VKM(4,2)+VKM(4,3)+VKM(4,5)
CCAO!
        WRITE(16,*) 'GFermi,Gmu =',GFERMI,GMU
        WRITE(16,*) 'WCOUPL,ZCOUPL_e, ZCOUPL_u,ZCOUPL_d,HCOUPL,ACOUPL'
        WRITE(16,*)  WCOUPL,ZCOUPL_e,ZCOUPL_u,ZCOUPL_d,HCOUPL,ACOUPL
        WRITE(16,*)'1/ALFA_EM_MZ  (ALPHA_{EM} AT MZ) =',1.0/ALFA_EM_MZ
        WRITE(16,*) ' SW2_EFF_e, SW2_EFF_u, SW2_EFF_d'
        WRITE(16,*)  SW2_EFF_e, SW2_EFF_u, SW2_EFF_d
        WRITE(16,*) ' ALPHAS_MW=', ALPHAS_MW
        WRITE(16,*) ' ALPHAS_MZ=', ALPHAS_MZ
        WRITE(16,*) 'XMW, XMZ =',XMW, XMZ
        WRITE(16,*) 'GAMW_NLO_R, GAMZ_NLO_R =',GAMW_NLO_R, GAMZ_NLO_R 
        WRITE(16,*) 'GAMW_MW, GAMZ_MZ =',GAMW_MW, GAMZ_MZ 
        WRITE(16,*) 'GAMW_NLO, GAMZ_NLO =',GAMW_NLO, GAMZ_NLO 
        WRITE(16,*) 'XMW_IN,GAMW_IN =',XMW_IN,GAMW_IN
        WRITE(16,*) 'XMZ_IN,GAMZ_IN =',XMZ_IN,GAMZ_IN
        WRITE(16,*) 'BR_Z_EE_NLO_R =', BR_Z_EE_NLO_R
        WRITE(16,*) 'BR_W_EN_NLO_R =', BR_W_EN_NLO_R
        WRITE(16,*) 'BR_Z_EE =', BR_Z_EE
        WRITE(16,*) 'BR_W_EN =', BR_W_EN
        WRITE(16,*) 'XMHIGGS, GAMH =',XMHIGGS, GAMH

CBIN
        WRITE(16,*) ' g_fvf_p_u, g_fvf_m_u =',g_fvf_p_u, g_fvf_m_u
        WRITE(16,*) ' g_fvf_p_d, g_fvf_m_d =',g_fvf_p_d, g_fvf_m_d
        WRITE(16,*) 'XLE, XRE, XLU, XLD, XRU, XRD=',XLE, XRE, XLU,
     &               XLD, XRU, XRD
CCPY!  Watch out the difference between GAMH and GAMH_TMP inside 
CCPY!  SUBROUTINE STANDARD 
CCPY!        WRITE(16,*) 'BR(H -> BB) =', H2BB/GAMH
CCPY!        WRITE(16,*) 'BR(H -> TT) =', H2TT/GAMH
CCPY!        WRITE(16,*) 'BR(H -> ZZ) =', H2ZZ/GAMH
CCPY!        WRITE(16,*) 'BR(H -> WW) =', H2WW/GAMH

        WRITE(16,*) ' GAMH=GAMH_HDECAY, GAMH_TMP =',GAMH,GAMH_TMP
        FIRST=.FALSE.
      ENDIF

CCPY FOR CLARITY, I CALL AN EXTERNAL SUBROUTINE FOR 'DY' PROCESS
      CALL DY_COUPL

      RETURN
      END ! Standard

C===========================
CCPY SET UP COUPLINGS FOR 'DY' 
C     IF(BOSON.EQ.'ZU') THEN
C THIS IS FOR UP-TYPE QUARKS IN THE INITIAL STATE
C       IFLAV=1
C       BOSON='DY'
C     ELSEIF(BOSON.EQ.'ZD') THEN
C THIS IS FOR DOWN-TYPE QUARKS IN THE INITIAL STATE
C       IFLAV=2
C       BOSON='DY'
C     ENDIF
C============================
      SUBROUTINE DY_COUPL
C-----------------------------------------------------------
      IMPLICIT NONE

      INTEGER IFLAV
      REAL*8 CAA_E,CAA_O,CZA_E,CZA_O,CZZ_E,CZZ_O
      COMMON/DY_IFLAV/ CAA_E,CAA_O,CZA_E,CZA_O,CZZ_E,CZZ_O,IFLAV
              
      REAL*8 XMW,XMZ,XMW2,XMZ2,AEM,GAMW,GAMZ,GAMW2,GAMZ2,
     &SWS,CWS,XMBOT,XMTOP,XMB2,XMT2,GWEAK2
      COMMON/STAND0/XMW,XMZ,XMW2,XMZ2,AEM,GAMW,GAMZ,GAMW2,GAMZ2
      COMMON/STAND1/SWS,CWS,XMBOT,XMTOP,XMB2,XMT2,GWEAK2
CCAO!
      REAL*8 GAMW_NLO,GAMZ_NLO,ALFA_EM_MZ,ACOUPL,WCOUPL,
     &ZCOUPL,HCOUPL,SW2_EFF,EW_ALFA
      COMMON/STAND7/ GAMW_NLO,GAMZ_NLO,ALFA_EM_MZ,ACOUPL,WCOUPL,
     >       ZCOUPL,HCOUPL,SW2_EFF,EW_ALFA

      REAL*8 EL_Z,ER_Z,QL_Z,QR_Z,EL_A,ER_A,QL_A,QR_A
      COMMON/ZFF_INIT/
     >   EL_Z,ER_Z,QL_Z(2),QR_Z(2),EL_A,ER_A,QL_A(2),QR_A(2)

c      REAL*8 SW2_EFF_u,SW2_EFF_d

      REAL*8 ALPHAS
      EXTERNAL ALPHAS

      REAL*8 FAC_A,FAC
      
      LOGICAL TESTING

CBIN
      REAL*8 rho_e,rho_u,rho_d
      REAL*8 ZCOUPL_e,ZCOUPL_u,ZCOUPL_d
      REAL*8 sw2_eff_e,sw2_eff_u,sw2_eff_d,
     >      cw2_eff_e,cw2_eff_u,cw2_eff_d
      COMMON/ZFITTER/rho_e,rho_u,rho_d,
     >      ZCOUPL_e,ZCOUPL_u,ZCOUPL_d,
     >     sw2_eff_e,sw2_eff_u,sw2_eff_d,
     >     cw2_eff_e,cw2_eff_u,cw2_eff_d   

C contral different electroweak parameters
      InTEGER*4 EW_FLAG,qed_FLAG,CKM_FLAG  
      COMMON/FLAG/EW_FLAG,qed_FLAG,CKM_FLAG 

CFEWZ0
      REAL*8 gvu,gau,gvd,gad,gvl,gal
      common/FEWZ0/ gvu,gau,gvd,gad,gvl,gal

CFEWZ1
      REAL*8 FAC_AFE
      REAL*8 QL_AFE(2),QR_AFE(2),EL_AFE,ER_AFE
      COMPLEX*8 QL_ZFE(2),QR_ZFE(2),EL_ZFE,ER_ZFE,FACFE,SWSFE
      REAL*8 CAA_EFE,CAA_OFE,CZZ_EFE,CZZ_OFE,CZA_EFE,CZA_OFE
       common/FEWZ/SWSFE

C ACOUPL IS (e**2)/2
C      ACOUPL=R2*MW2*GMU*2.0D0*SW2
C ZCOUPL IS (g/Cos_w)**2/4
C_yfu open up ZCOUPL
C      ZCOUPL=R2*MZ2*GMU*(1.0+EW_DELRHO)

C COUPLINGS OF A0-F-F.  (g*Sin_w)*(\gamma_\mu)*(electric charge)
      FAC_A=DSQRT(2.0*ACOUPL)/2.d0
      QL_A(1) = (2.0/3.0)*FAC_A
      QR_A(1) = (2.0/3.0)*FAC_A
      QL_A(2) = (-1.0/3.0)*FAC_A
      QR_A(2) = (-1.0/3.0)*FAC_A
      EL_A = -1.0*FAC_A
      ER_A = -1.0*FAC_A

C COUPLINGS OF Z0-F-F.  -(g/2/Cos_w)*(\gamma_\mu)*(V^I-A^I*\gamma_5)
C  WHICH IS DEFINED AS  -(g/2/Cos_w)*(\gamma_\mu)*(FL*(1-\GAMMA_5) + FR*(1+\GAMMA_5))
      FAC = DSQRT(ZCOUPL)
CCPY!      QL_Z(1) = (0.5-2.0/3.0*SW2_EFF)*FAC
CCPY!      QR_Z(1) = -2.0/3.0*SW2_EFF*FAC
CCPY!      QL_Z(2) = (-0.5+SW2_EFF/3.0)*FAC
CCPY!      QR_Z(2) = (SW2_EFF/3.0)*FAC
c      EL_Z = (-0.5+SW2_EFF)*FAC
c      ER_Z = SW2_EFF*FAC
     
      IF(EW_FLAG.le.2) THEN
         EL_Z = (-0.5+SW2_EFF_e)*DSQRT(ZCOUPL_e)
         ER_Z = SW2_EFF_e*DSQRT(ZCOUPL_e)
      ELSEIF(EW_FLAG.eq.3) THEN
         EL_Z = (gvl-gal)*FAC/2.
         ER_Z = (gvl+gal)*FAC/2.
      ENDIF

      IF(EW_FLAG.eq.4) THEN
      FAC_AFE=DSQRT(2.0*ACOUPL)/2.d0
      QL_AFE(1) = (2.0/3.0)*FAC_AFE
      QR_AFE(1) = (2.0/3.0)*FAC_AFE
      QL_AFE(2) = (-1.0/3.0)*FAC_AFE
      QR_AFE(2) = (-1.0/3.0)*FAC_AFE
      EL_AFE = -1.0*FAC_AFE
      ER_AFE = -1.0*FAC_AFE

      FAC_A=FAC_AFE
      QL_A(1) = QL_AFE(1)
      QR_A(1) = QR_AFE(1)
      QL_A(2) = QL_AFE(2)
      QR_A(2) = QR_AFE(2)
      EL_A = EL_AFE
      ER_A = ER_AFE

      FACFE = DSQRT(ZCOUPL)  
      FAC   =REAL(FACFE)  

      EL_ZFE = (-0.5+SWSFE)*FACFE
      ER_ZFE = SWSFE*FACFE
      EL_Z = REAL(EL_ZFE) 
      ER_Z = REAL(ER_ZFE) 
      ENDIF

CCPY June 2010: Add the difference in SW2_EFF_u and SW2_EFF_d
C For up-type quark
c      SW2_EFF_u=SW2_EFF-0.0001
c      SW2_EFF_d=SW2_EFF-0.0002
      
      IF(IFLAV.EQ.1) THEN
C THIS IS FOR UP-TYPE QUARKS IN THE INITIAL STATE
c        QL_Z(1) = (0.5-2.0/3.0*SW2_EFF_u)*FAC
c        QR_Z(1) = -2.0/3.0*SW2_EFF_u*FAC
c        QL_Z(2) = (-0.5+SW2_EFF_u/3.0)*FAC
c        QR_Z(2) = (SW2_EFF_u/3.0)*FAC
 	IF(EW_FLAG.le.2) THEN
           QL_Z(1) = (0.5-2.0/3.0*SW2_EFF_u)*DSQRT(ZCOUPL_u)
           QR_Z(1) = -2.0/3.0*SW2_EFF_u*DSQRT(ZCOUPL_u)
           QL_Z(2) = (-0.5+SW2_EFF_u/3.0)*DSQRT(ZCOUPL_u)
           QR_Z(2) = (SW2_EFF_u/3.0)*DSQRT(ZCOUPL_u)
        ELSEIF(EW_FLAG.eq.3) THEN
           QL_Z(1) = (gvu-gau)*FAC/2.
           QR_Z(1) = (gvu+gau)*FAC/2.
           QL_Z(2) = (gvd-gad)*FAC/2.
           QR_Z(2) = (gvd+gad)*FAC/2.
        ENDIF

        CAA_E = ( QL_A(1)*QL_A(1) + QR_A(1)*QR_A(1) )*
     >          ( EL_A*EL_A + ER_A*ER_A )
        CAA_O = ( QL_A(1)*QL_A(1) - QR_A(1)*QR_A(1) )*
     >          ( EL_A*EL_A - ER_A*ER_A )

        CZZ_E = ( QL_Z(1)*QL_Z(1) + QR_Z(1)*QR_Z(1) )*
     >          ( EL_Z*EL_Z + ER_Z*ER_Z )
        CZZ_O = ( QL_Z(1)*QL_Z(1) - QR_Z(1)*QR_Z(1) )*
     >          ( EL_Z*EL_Z - ER_Z*ER_Z )

        CZA_E = ( QL_Z(1)*QL_A(1) + QR_Z(1)*QR_A(1) )*
     >          ( EL_Z*EL_A + ER_Z*ER_A )
        CZA_O = ( QL_Z(1)*QL_A(1) - QR_Z(1)*QR_A(1) )*
     >          ( EL_Z*EL_A - ER_Z*ER_A )

        IF(EW_FLAG.eq.4) THEN
        QL_ZFE(1) = (0.5-2.0/3.0*SWSFE)*FACFE
        QR_ZFE(1) = -2.0/3.0*SWSFE*FACFE
        QL_ZFE(2) = (-0.5+SWSFE/3.0)*FACFE
        QR_ZFE(2) = (SWSFE/3.0)*FACFE

        CAA_EFE = ( QL_AFE(1)*QL_AFE(1) + QR_AFE(1)*QR_AFE(1) )*
     >          ( EL_AFE*EL_AFE + ER_AFE*ER_AFE )
        CAA_OFE = ( QL_AFE(1)*QL_AFE(1) - QR_AFE(1)*QR_AFE(1) )*
     >          ( EL_AFE*EL_AFE - ER_AFE*ER_AFE )

        CZZ_EFE = ( QL_ZFE(1)*CONJG(QL_ZFE(1)) + 
     >              QR_ZFE(1)*CONJG(QR_ZFE(1) ))*
     >          ( EL_ZFE*CONJG(EL_ZFE) + ER_ZFE*CONJG(ER_ZFE) )
        CZZ_OFE = ( QL_ZFE(1)*CONJG(QL_ZFE(1)) -
     >              QR_ZFE(1)*CONJG(QR_ZFE(1) ))*
     >          ( EL_ZFE*CONJG(EL_ZFE) - ER_ZFE*CONJG(ER_ZFE) )

        CZA_EFE =REAL( ( QL_ZFE(1)*QL_AFE(1) + QR_ZFE(1)*QR_AFE(1) )*
     >          ( EL_ZFE*EL_AFE + ER_ZFE*ER_AFE ))
        CZA_OFE =REAL( ( QL_ZFE(1)*QL_AFE(1) - QR_ZFE(1)*QR_AFE(1) )*
     >          ( EL_ZFE*EL_AFE - ER_ZFE*ER_AFE ))



        QL_Z(1) =REAL( QL_ZFE(1))
        QR_Z(1) =REAL( QR_ZFE(1))
        QL_Z(2) =REAL( QL_ZFE(2))
        QR_Z(2) =REAL( QR_ZFE(2))

        CAA_E = CAA_EFE
        CAA_O = CAA_OFE

        CZZ_E = CZZ_EFE
        CZZ_O = CZZ_OFE

        CZA_E = CZA_EFE
        CZA_O = CZA_OFE 
        ENDIF

      ELSE
C THIS IS FOR DOWN-TYPE QUARKS IN THE INITIAL STATE
c        QL_Z(1) = (0.5-2.0/3.0*SW2_EFF_d)*FAC
c        QR_Z(1) = -2.0/3.0*SW2_EFF_d*FAC
c        QL_Z(2) = (-0.5+SW2_EFF_d/3.0)*FAC
c        QR_Z(2) = (SW2_EFF_d/3.0)*FAC
 
        IF(EW_FLAG.LE.2) THEN
           QL_Z(1) = (0.5-2.0/3.0*SW2_EFF_d)*DSQRT(ZCOUPL_d)
           QR_Z(1) = -2.0/3.0*SW2_EFF_d*DSQRT(ZCOUPL_d)
           QL_Z(2) = (-0.5+SW2_EFF_d/3.0)*DSQRT(ZCOUPL_d)
           QR_Z(2) = (SW2_EFF_d/3.0)*DSQRT(ZCOUPL_d)
        ELSEIF(EW_FLAG.eq.3) THEN
           QL_Z(1) = (gvu-gau)*FAC/2.
           QR_Z(1) = (gvu+gau)*FAC/2.
           QL_Z(2) = (gvd-gad)*FAC/2.
           QR_Z(2) = (gvd+gad)*FAC/2.
        ENDIF

        CAA_E = ( QL_A(2)*QL_A(2) + QR_A(2)*QR_A(2) )*
     >          ( EL_A*EL_A + ER_A*ER_A )
        CAA_O = ( QL_A(2)*QL_A(2) - QR_A(2)*QR_A(2) )*
     >          ( EL_A*EL_A - ER_A*ER_A )

        CZZ_E = ( QL_Z(2)*QL_Z(2) + QR_Z(2)*QR_Z(2) )*
     >          ( EL_Z*EL_Z + ER_Z*ER_Z )
        CZZ_O = ( QL_Z(2)*QL_Z(2) - QR_Z(2)*QR_Z(2) )*
     >          ( EL_Z*EL_Z - ER_Z*ER_Z )

        CZA_E = ( QL_Z(2)*QL_A(2) + QR_Z(2)*QR_A(2) )*
     >          ( EL_Z*EL_A + ER_Z*ER_A )
        CZA_O = ( QL_Z(2)*QL_A(2) - QR_Z(2)*QR_A(2) )*
     >          ( EL_Z*EL_A - ER_Z*ER_A )

       IF(EW_FLAG.eq.4) THEN

        QL_ZFE(1) = (0.5-2.0/3.0*SWSFE)*FACFE
        QR_ZFE(1) = -2.0/3.0*SWSFE*FACFE
        QL_ZFE(2) = (-0.5+SWSFE/3.0)*FACFE
        QR_ZFE(2) = (SWSFE/3.0)*FACFE

        CAA_EFE = ( QL_AFE(2)*QL_AFE(2) + QR_AFE(2)*QR_AFE(2) )*
     >          ( EL_AFE*EL_AFE + ER_AFE*ER_AFE )
        CAA_OFE = ( QL_AFE(2)*QL_AFE(2) - QR_AFE(2)*QR_AFE(2) )*
     >          ( EL_AFE*EL_AFE - ER_AFE*ER_AFE )

        CZZ_EFE = ( QL_ZFE(2)*CONJG(QL_ZFE(2)) + 
     >              QR_ZFE(2)*CONJG(QR_ZFE(2) ))*
     >          ( EL_ZFE*CONJG(EL_ZFE) + ER_ZFE*CONJG(ER_ZFE) )
        CZZ_OFE = ( QL_ZFE(2)*CONJG(QL_ZFE(2)) - 
     >              QR_ZFE(2)*CONJG(QR_ZFE(2) ))*
     >          ( EL_ZFE*CONJG(EL_ZFE) - ER_ZFE*CONJG(ER_ZFE) )

        CZA_EFE =REAL( ( QL_ZFE(2)*QL_AFE(2) + QR_ZFE(2)*QR_AFE(2) )*
     >          ( EL_ZFE*EL_AFE + ER_ZFE*ER_AFE ))
        CZA_OFE =REAL( ( QL_ZFE(2)*QL_AFE(2) - QR_ZFE(2)*QR_AFE(2) )*
     >          ( EL_ZFE*EL_AFE - ER_ZFE*ER_AFE ))



        QL_Z(1) =REAL( QL_ZFE(1))
        QR_Z(1) =REAL( QR_ZFE(1))
        QL_Z(2) =REAL( QL_ZFE(2))
        QR_Z(2) =REAL( QR_ZFE(2))

        CAA_E = CAA_EFE
        CAA_O = CAA_OFE

        CZZ_E = CZZ_EFE
        CZZ_O = CZZ_OFE

        CZA_E = CZA_EFE
        CZA_O = CZA_OFE 
	ENDIF


      ENDIF

      WRITE(16,*) 'EL_Z, ER_Z, QL_Z1, QR_Z1, QL_Z2, QR_Z2'
      WRITE(16,*)EL_Z,ER_Z, QL_Z(1),QR_Z(1),QL_Z(2),QR_Z(2)

      WRITE(16,*) ' CZZ_E,CZZ_O,CZA_E,CZA_O,CAA_E,CAA_O,IFLAV'
      WRITE(16,*) CZZ_E,CZZ_O,CZA_E,CZA_O,CAA_E,CAA_O,IFLAV

C_yfu print out the electroweak parameter
c      print*,'gvl = ',gvl,' gal = ',gal,' FAC = ',DSQRT(ZCOUPL_e)
c      print*,'QL_Z = ',QL_Z(1),' QL_A = ',QL_A(1),' QR_Z = ',QR_Z(1),
c     >' QR_A = ',QR_A(1)
c      print*,'EL_Z = ',EL_Z,' EL_A = ',EL_A,' ER_Z = ',ER_Z,
c     >' ER_A = ',ER_A

c      print*,'QL_Z(1)*QL_Z(1) + QR_Z(1)*QR_Z(1): ',
c     >QL_Z(1)*QL_Z(1) + QR_Z(1)*QR_Z(1)
c      print*,'QL_Z(2)*QL_Z(2) + QR_Z(2)*QR_Z(2): ',
c     >QL_Z(2)*QL_Z(2) + QR_Z(2)*QR_Z(2)
c      print*,'EL_Z*EL_Z + ER_Z*ER_Z: ',EL_Z*EL_Z + ER_Z*ER_Z


      RETURN
      END


C---------------------------------------------


      FUNCTION FXN(X,WGT)
C--------------------------------------------------
CsB___Function integrated by Vegas
      IMPLICIT NONE

      REAL*8 XL,XU,ACC
      INTEGER*4 NDIM,NCALL,ITMX,NPRN,NIN,NOUT,NBINIT, iwgt
      REAL*8 XI,SI,SI2,SWGT,SCHI,ALPH
      INTEGER*4 NDO,IT
      LOGICAL IBINNER
      REAL*8 MAX_WGT,NEW_MAX_WGT,WGT_FACTOR
      LOGICAL UNWGT_FLAG
      REAL*8 SEC_PLAB,PLAB,PCEM,WHAT

      COMMON/BVEG1/XL(50),XU(50),ACC,NDIM,NCALL,ITMX,NPRN
      COMMON/BVEG2/XI(50,50),SI,SI2,SWGT,SCHI,NDO,IT,ALPH
      COMMON/EFFIC/NIN,NOUT
      COMMON/NININF/IBINNER
      COMMON/UNWEIGHT/MAX_WGT,NEW_MAX_WGT,WGT_FACTOR,IWGT,NBINIT
      COMMON/DECAY_MOMENTA/SEC_PLAB(4,50)
      COMMON/MOM/PLAB(4,12),PCEM(4,12),WHAT

      REAL*8 YXP,YXETA,YXPT,YXDELR,YXPP
      COMMON/SECCONV/ YXP(4,30),YXETA(30),YXPT(30),YXDELR(30,30),
     >YXPP(30)
      REAL*8 YP,YETA,YPT,YDELR
      COMMON/PASSCONV/ YP(4,30),YETA(30),YPT(30),YDELR(30,30)

      REAL*8 ZERO
      REAL*8 RN
      EXTERNAL RN
      REAL*8 X,WT,WGT,FXN,PASSER
      LOGICAL IDUMP,IN_IT
      DIMENSION X(50)
      DATA IN_IT/.TRUE./

      save

      IF(IN_IT)THEN
        IN_IT=.FALSE.
        NBINIT=0
        IF(WGT_FACTOR.EQ.0.)WGT_FACTOR=1.D0
      END IF
      FXN=0.D0

CsB___Generate weight
      CALL MAKE_EVENT(X,WT,IDUMP)
      
      IF(.NOT.IDUMP)THEN
        IF(.NOT.IBINNER) THEN
          RETURN
        ELSE
          NOUT=NOUT+1
          RETURN
        END IF
      END IF

CsB___Pass the weight to Vegas
      FXN=WT

      IF(UNWGT_FLAG) iwgt=1
      IF((.NOT.IBINNER).AND. (iwgt.eq.1) .AND. IT .EQ. ITMX)THEN
        MAX_WGT = MAX(MAX_WGT,WT*WGT*FLOAT(NCALL))
      ELSE IF(IBINNER)THEN
        NIN=NIN+1
        IF(iwgt.eq.1)THEN ! unweighted events
          NEW_MAX_WGT = MAX(NEW_MAX_WGT,WT*WGT*FLOAT(NCALL))
          IF(MAX_WGT*WGT_FACTOR .LT. WT*WGT*FLOAT(NCALL))THEN
            NBINIT = NBINIT + 1
            CALL PASS_OUTPUT(FXN*WGT/FLOAT(ITMX),ZERO) ! record the event
          ELSE IF(RN(NBINIT)*MAX_WGT*WGT_FACTOR .LT.
     >            WT*WGT*FLOAT(NCALL))THEN
            CALL PASS_OUTPUT(MAX_WGT*WGT_FACTOR/FLOAT(NCALL*ITMX),
     >        FXN*WGT/FLOAT(ITMX)) ! pass the usual weight in WT1
            NBINIT = NBINIT + 1
          END IF
        ELSE !iwgt. ne. 1
          PASSER=FXN*WGT/FLOAT(ITMX)
          IF(PASSER.NE.0.) CALL PASS_OUTPUT(PASSER, ZERO) ! record the event
        END IF                  !iwgt
      END IF

      RETURN
      END ! FXN

C ---------------------------------------------------------------------------
      SUBROUTINE PHASE_GEN(Y,N,ET,XM,P,WT)
C ---------------------------------------------------------------------------
CCPY
C PHASE_GEN HAS NO PI'S IN IT, TO INCLUDE THIS FACTOR, WE SHOULD ADD
C        WT=WT*(2.*PI)**(4.-3*N)
C IN THE END OF THIS SUBROUTINE
C
CsB___The (2.*PI)**(4.-3*N) is included after PHASE_GEN is called.
C
      IMPLICIT NONE
      REAL*8 Z,XM,P,Q,R,XM2,E,V,WT,ACC,PI2LOG,XMT,ET,XMAX
      REAL*8 C,S,F,FAC,RMAS,X,ACCU,WT2,WT3,WTM,P2
      REAL*8 F0,G0,X2,EXTRA_WEIGHT,Y,PI
      INTEGER*4 IN_IT,ITMAX,K,N,NM,ITER,I
      DIMENSION XM(100),P(4,100),Q(4,100),R(4),Z(200),
     .   P2(100),XM2(100),E(100),V(100),Y(50)
      DATA ACC/1.D-14/,ITMAX/6/,IN_IT/0/
      DATA PI2LOG/1.837877D0/, Pi / 3.141592653589793 /
C PI2LOG IS LOG(2*PI)
C
      IF(IN_IT.EQ.0) THEN
        IN_IT=1
        Z(1)=0.D0
        DO 83 K=2,200
          Z(K)=Z(K-1)+LOG(FLOAT(K))
83      CONTINUE
      END IF
C COUNT NONZERO MASSES
      XMT=0.D0
      NM=0
      DO 87 I=1,N
        IF(XM(I).GT.1.E-6) NM=NM+1
        XMT=XMT+ABS(XM(I))
87    CONTINUE

C GENERATE N MASSLESS MOMENTA IN INFINITE PHASE SPACE
      EXTRA_WEIGHT=1.D0
      DO 89 I=1,N-1
        IF(I.EQ.1)THEN
          C=2.D0*Y(3)-1.
          S=SQRT(1.-C*C)
          F=0.D0
C THERE IS ONE REDUNDANT OVERALL PHI
          FAC=Y(4)
        ELSE
          C=2.D0*Y(2+3*(I-1))-1.
          S=SQRT(1.-C*C)
          F=2.D0*PI*Y(3+3*(I-1))
          FAC=Y(4+3*(I-1))
        END IF
        Q(4,I)=-LOG(FAC)
        EXTRA_WEIGHT=EXTRA_WEIGHT*Q(4,I)
        Q(3,I)=Q(4,I)*C
        Q(2,I)=Q(4,I)*S*COS(F)
        Q(1,I)=Q(4,I)*S*SIN(F)
89    CONTINUE
C CALCULATE THE PARAMETERS OF THE SCALE TRANSFORMATION
      DO 91 I=1,4
        R(I)=0.D0
91    CONTINUE
      DO 93 I=1,N-1
        DO 93 K=1,4
          R(K)=R(K)+Q(K,I)
93    CONTINUE
C CALCULATE THE N TH VECTOR
      DO 97 K=1,3
        Q(K,N)=-R(K)
97    CONTINUE
      Q(4,N)=SQRT(Q(1,N)**2+Q(2,N)**2+Q(3,N)**2)
      R(4)=R(4)+Q(4,N)
      RMAS=R(4)
      IF(RMAS.EQ.0.)THEN
        WT=0.D0
        RETURN
      END IF
      X=ET/RMAS
C SCALE THE Q'S TO THE P'S
      DO 111 I=1,N
        DO 111 K=1,4
          P(K,I)=X*Q(K,I)
111   CONTINUE

C RETURN FOR WEIGHTED MASSLESS MOMENTA
      WT=PI2LOG*(N-1)
     1 -LOG(ET-P(4,N))*2.*(1-N)-LOG(ET)-LOG(2.)-LOG(P(4,N))-Z(2*N-3)
      IF(NM.EQ.0) THEN
        WT=EXTRA_WEIGHT*EXP(WT)
        RETURN
      END IF

C MASSIVE PARTICLES: RESCALE THE MOMENTA BY A FACTOR X
      XMAX=SQRT(1.-(XMT/ET)**2)
      DO 113 I=1,N
        XM2(I)=XM(I)**2
        P2(I)=P(4,I)**2
113   CONTINUE
      ITER=0
      X=XMAX
      ACCU=ET*ACC
      DO WHILE (ITER.LE.ITMAX)
        F0=-ET
        G0=0.D0
        X2=X*X
        DO 117 I=1,N
          E(I)=SQRT(XM2(I)+X2*P2(I))
          F0=F0+E(I)
          IF(E(I).GT.0.)THEN
            G0=G0+P2(I)/E(I)
          END IF
117     CONTINUE
        IF(ABS(F0).GT.ACCU.AND.ITER.LE.ITMAX) THEN
          ITER=ITER+1
          X=X-F0/(X*G0)
        ELSE IF(ABS(F0).LE.ACCU) THEN
          ITER=ITMAX+1
        END IF
      END DO
      DO 121 I=1,N
        V(I)=X*P(4,I)
        DO 127 K=1,3
          P(K,I)=X*P(K,I)
127     CONTINUE
        P(4,I)=E(I)
121   CONTINUE

C CALCULATE THE MASS-EFFECT WEIGHT FACTOR
      WT2=1.D0
      WT3=0.D0
      DO 131 I=1,N
        IF(E(I).GT.0.)THEN
          WT2=WT2*V(I)/E(I)
          WT3=WT3+V(I)**2/E(I)
        END IF
131   CONTINUE
      WTM=(2.*N-3.)*LOG(X)+LOG(WT2/WT3*ET)

C RETURN FOR  WEIGHTED MASSIVE MOMENTA
      WT=WT+WTM
      WT=DEXP(WT)*EXTRA_WEIGHT
      IF(XMT.GT.ET)THEN
        WT=0.
        RETURN
      END IF
      RETURN
      END ! PHASE_GEN

C --------------------------------------------------------------------------
      SUBROUTINE LORENTZ(PB,PC,PL)
C --------------------------------------------------------------------------
      IMPLICIT NONE
C THIS PERFORMS A LORENTZ TRANSFORMATION
      REAL*8 PB(4),PC(4),PL(4),Q,F
      INTEGER*4 J
      Q=SQRT(PB(4)**2-PB(1)**2-PB(2)**2-PB(3)**2)
      PL(4)=(PB(4)*PC(4)+PB(3)*PC(3)
     1       +PB(2)*PC(2)+PB(1)*PC(1))/Q
      F=(PL(4)+PC(4))/(Q+PB(4))
      DO  137 J=1,3
        PL(J)=PC(J)+F*PB(J)
137   CONTINUE
      RETURN
      END ! LORENTZ

C -----------------------------------------------------------------------

      subroutine kinemprint(wt_old,ps_old, wt, ps)
      implicit none
      real*8 wt_old,ps_old, wt, ps
      REAL*8 Q_V,QT_V,YB_V,THE_STA,PHI_STA
      COMMON/CSFRAME/Q_V,QT_V,YB_V,THE_STA,PHI_STA
      save

      print *,'wt_old,ps_old, wt, ps'
      print *, wt_old, ps_old, wt, ps
      print *
      print *, 'Q_V,QT_V,YB_V,THE_STA,PHI_STA ='
      print *,Q_V,QT_V,YB_V,THE_STA,PHI_STA
      
      return 
      end !kinemprint

C -----------------------------------------------------------------------

C Originally: THIS IS FOR THE RESUMMED W^+
C AS OF 10/1995, THIS WORKS FOR W+, W-, Z0 AND A0.
C This routine calculates the weight for processes with grids.
C This weight contains everything except the phase-weight.
C This routine is correct only in the Collins-Soper frame.
C --------------------------------------------------------------------------
      SUBROUTINE WPRES(WT)
C --------------------------------------------------------------------------
      IMPLICIT NONE

      REAL*8 WT
      REAL*8 W,QCDL,X1X2MIN
      INTEGER*4 NJET,N_PRIMARY_PART,NF,I_STR,IBEAM,NPART,NEVENT,NORDER,
     ,          N_SECOND_PART
      INTEGER*4 IW
      REAL*8 PLAB,PCEM,WHAT,X1,X2,Q,Q1,Q2
      INTEGER*4 IERR_WT
      REAL*8 XMASS,X_DECAYMASS

      COMMON/AA/W,QCDL,X1X2MIN,NJET,N_PRIMARY_PART,NF,I_STR,IBEAM,NPART,
     &  NEVENT,NORDER,N_SECOND_PART
      COMMON/CHARGE/IW
      COMMON/MOM/PLAB(4,12),PCEM(4,12),WHAT
      COMMON/PMF/X1,X2
      COMMON/SCALE/Q

      Real*8 qQBUL,qQBDL,QBqUL,QBqDL,qGUL,qGDL,GqUL,GqDL,GGL,LParton
      Real*8 ECM
      Common /ECMC/ ECM
      Integer KinCorr
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

CCAO! The read W-boson mass and width
      Real*8 XMW_IN,XMW_GRID,GAMW_IN
      COMMON/WMASS_IN/XMW_IN,XMW_GRID,GAMW_IN
      Real*8 XMZ_IN,XMZ_GRID,GAMZ_IN
      COMMON/ZMASS_IN/XMZ_IN,XMZ_GRID,GAMZ_IN

      REAL*8 Q_V,QT_V,YB_V,THE_STA,PHI_STA
      COMMON/CSFRAME/Q_V,QT_V,YB_V,THE_STA,PHI_STA
      REAL*8 YMW_AMP,YGW_AMP
      COMMON/VMASS/YMW_AMP,YGW_AMP

      REAL*8 CONST,CONST1,CONST2,CONST3,PI,PROP,WTS,WTA

CsB   3DInterp:
      Integer MD1,MD2,MD3, ID1,ID2,ID3,  IY1,IY2,IY3
      Real*8 X1A,X2A,X3A, YAS,YAA, YDAS,YDAA, X1Y,X2Y,X3Y, YAY
      Parameter (MD1 = 200, MD2 = 300, MD3 = 100)
      Dimension X1A(0:MD1), X2A(0:MD2), X3A(0:MD3)
      Dimension YAS(0:MD1,0:MD2,0:MD3), YAA(0:MD1,0:MD2,0:MD3)
      Dimension YDAS(0:MD1,0:MD2,0:MD3), YDAA(0:MD1,0:MD2,0:MD3)
      Dimension X1Y(0:MD1), X2Y(0:MD2), X3Y(0:MD3)
      Dimension YAY(0:MD1,0:MD2,0:MD3,0:7)
      Common / Arrays / X1A,X2A,X3A, YAS,YAA,  X1Y,X2Y,X3Y, YAY
      Common / Dims / iD1,iD2,iD3, iY1,iY2,iY3
      INTEGER N_QT_INTP,N_Y_INTP,N_Q_INTP

      Character*40 Boson
      Common / Boson / Boson
      Integer N_SUD_AB,N_SUD_A,N_SUD_B,N_WIL_C,I_SUD
      Common / Sudakov / N_SUD_AB,N_SUD_A,N_SUD_B,N_WIL_C,I_SUD
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
      Dimension SigY (0:5)
      Real*8 fL,fR, fLA,fRA,fLZ,fRZ, gWeak, sWeak, cWeak
      Real*8 S_Q
      Real*8 Wt1, KFac
      Real*8 sH,tH,uH, CJJ1XFJ_1_AA,SIGS_AA, SIGS_AG,SIGA_AG, 
     &  CJJ1XFJ_1_ZZ,SIGS_ZZ,SigR_ZZ,SigP_ZZ, x,Mz2,Mz4
      External CJJ1XFJ_1_AA, CJJ1XFJ_1_ZZ
      Real*8 XMAT_PHOT
      Real*8 GGAA
      Real*8 ggAA_A2,ggAA_D2
      REAL*8 ETOP,BTOP,CTHE,STHE,BCTH

      Real*8 QT_V_Mat
      Common / Matching / QT_V_Mat
      Integer iOneQ
      Common / OneQ / iOneQ
      Integer iYPSw
      Common / YPertSwitch / iYPSw
      Real*8 PP0,PP3,YP0,YP3
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
      
      Real*8 Wwidth
      Common / InputWwidth / Wwidth

      Real*8 Propagator, DeNom, H2v, Sigma0

      Integer iResScheme, iLV
      Real*8 h1v, h1L
      Common/ ReSumScheme / iResScheme,iLV,h1v,h1L
      Real*8 C1,C2,C3,C4
      Common / LegacyCi / C1,C2,C3,C4

      Double Precision vin,h1ggAA,fmlo2, fmlo2prime
CCAO  
      INTEGER*4 IFLAG
      COMMON/RET/ IFLAG

      Data Pi / 3.141592653589793 /
      Data (SigY(ii), ii=0,5) / 0.d0,0.d0,0.d0,0.d0,0.d0,0.d0/
      Data FacN / 1.d0 /
      Real*8 CF
      Data CF / 1.333333333333333 /

CsB   
      REAL*8 MAX_WGT,NEW_MAX_WGT,WGT_FACTOR
      LOGICAL UNWGT_FLAG
      INTEGER*4 NBINIT, iwgt
      COMMON/UNWEIGHT/MAX_WGT,NEW_MAX_WGT,WGT_FACTOR,IWGT,NBINIT
      Real*8 SgnWt
      Common /SignWeight/ SgnWt

CCPY
      LOGICAL Y_ONE_SIDE
      COMMON/GRID_FILE_FORM/ Y_ONE_SIDE

CCPY
      INTEGER IFLAV
      REAL*8 CAA_E,CAA_O,CZA_E,CZA_O,CZZ_E,CZZ_O
      COMMON/DY_IFLAV/ CAA_E,CAA_O,CZA_E,CZA_O,CZZ_E,CZZ_O,IFLAV

CCPY jAN 16, 2011
      REAL*8 ASY_TC
      COMMON/GL_ASY/ ASY_TC
      
CPY
      REAL*8 ZZ_PROP,AA_PROP,ZA_PROP,ALFA_Q2,FACTOR
      REAL*8 ANGFUNC_EVEN,ANGFUNC_ODD,COEFF_EVEN,COEFF_ODD

      real*8 test_1,test_2
      
      LOGICAL TESTING

CJI 2019
      INTEGER kloc
      REAL*8 a0k, a1k, a2k, a4k
      Dimension a0k(1:13), a1k(1:13), a2k(1:13), a4k(1:13)
      Data a0k / 1.0, 1.04942, 0.99937, 0.98386, 
     &       0.98072, 0.98752, 0.99214, 0.99911,
     &       0.99931, 1.00636, 1.00651, 1.00556, 1.00584 /
      Data a1k / 1.0, 1.38766, 1.32044, 1.27543,
     &       1.23318, 1.19356, 1.15258, 1.10204,
     &       1.05422, 1.01613, 0.97663, 0.93760, 0.890199 /
      Data a2k / 1.0, 0.67204, 0.72555, 0.76808,
     &       0.79788, 0.82768, 0.85473, 0.88185,
     &       0.90417, 0.91681, 0.92883, 0.94021, 0.95305 /
      Data a4k / 1.0, 1.02092, 1.02134, 1.02168,
     &       1.03259, 1.02571, 1.01893, 1.00042,
     &       0.98253, 0.97951, 0.94404, 0.92094, 0.90040 /

CJI NOV 2020
CCPY OCT 2019
      INTEGER*4 KFAC_ANGFUNC
      COMMON/ANGFUNC/ KFAC_ANGFUNC

      IFLAG=0

CsB___Initialization
      SigS = 0.d0
      SigA = 0.d0
      SigY(0) = 0.d0
      SigY(1) = 0.d0
      SigY(2) = 0.d0
      SigY(3) = 0.d0
      SigY(4) = 0.d0
      SigY(5) = 0.d0

CJI NOV 2020
CCPY OCT 2019
C      KFAC_ANGFUNC=0

CCPY THIRD ORDER OF POLYNOMIAL USED IN INTERPOLATION
      N_QT_INTP=4
      N_Y_INTP=4
      N_Q_INTP=4

CCPY Sept 2009: Assuming X2A(1)=0.0d0, i.e. y=0.0 is the first y-bin in 
C one-side y-grid file.

      IF(Y_ONE_SIDE) THEN
CCPY Use 2nd polynomial to do the interpolation.
        IF(DABS(YB_V).LE.X2A(2)) N_Y_INTP=4
c        print*,' YB_V,X2A(1),X2A(2) =', YB_V,X2A(1),X2A(2)
      EndIf

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

      If (nPart.NE.11) then ! interpolate the Asymptotic/Delta-sigma/LO piece
        Call ThreeDInterp ( N_QT_INTP,N_Y_INTP,N_Q_INTP,
     >      YDAS,YDAA, QT_V,yB_V,Q_V,SigS,SigA, ErrSgS,ErrSgA )
        If (iYpiece.Eq.1) then ! Interpolate the Y piece
C         Print*, ' Calling ThreeYInterp '
cC        If (QT_V.LT.QT_V_Sep) then
cCCsB_______If iYPSw=1 then send the Pert piece back
cC          iYPSw = 1
cC        Else
cCCsB_______If iYPSw=0 then send the Y piece back
            iYPSw = 0
cC        End If
          Call ThreeYInterp ( N_QT_INTP,N_Y_INTP,N_Q_INTP,
     >        YDAS,YDAA, QT_V,yB_V,Q_V, SigY, ErrSgY )
          YP0 = SigY(0)
          YP3 = SigY(3)
CsB       Checking interpolation
!         If (jj.LT.10 .Or. (jj.LT.10010 .and. jj.GT.10000)) then
!           Print*, Q_V,QT_V,yB_V,SigY(0)
!         End If
!         jj = jj + 1
!         If (Q_V.GT.36.8 .AND.  Q_V.LT.38.4) then
!           Print*, 'Q_V,QT_V,yB_V,SigY(0),SigS:'
!           Print*, Q_V,QT_V,yB_V,SigY(0),SigS
!         End If
        End If
      Else ! nPart.EQ.11 = interpolate the W (resummed) piece
        If (iYpiece.Eq.1) then
CsB_______For the matched resummed case the CSS piece is included in the Y grid
C         This way the CSS piece interpolation does not have to be called.
          SigS = 0.d0
          SigA = 0.d0

CsB_______The part below is based on the following observations:
C         1. The Y piece is small for the qQ contribution in the low QT
C         2. The fake Y approximates well the Y piece in the low QT
C         3. After the matching point the fake and real Pert curves have
C            the same slope.
C         Using the above we can use resbos to match the qQ contribution.
C
CsB_______For the ZZ case we need specially treat the "qqB_delta(1-z)" piece 
C         which is stored in SigA. This piece is needed in the low QT region.
C         If it is included after the matching point at high QT then it ruins
C         the Pert piece.

          If (Boson.Eq.'ZZ') then
            FacN     = 1.d0
            FacMatch = 1.d0
            
            If (W.GT.2.5d3) then ! LHC
              If (iProc.Eq.0) then ! q Q + q G
                If (QT_V.LT.150.) then
CsB_______________Interpolate the SigA piece
                  Call ThreeDInterp ( N_QT_INTP,N_Y_INTP,N_Q_INTP,
     >              YDAS,YDAA, QT_V,yB_V,Q_V,SigS,SigA, ErrSgS,ErrSgA )
C                 This is a K-factor applied to the fake Pert piece
                  FacMatch = 1.d0
                Else
CCPY              FacMatch = 1.d0-QT_V/800.
                  FacMatch = 1.d0
                End If
              Else If (iProc.Eq.1) then ! q Q 
                If (QT_V.LT.150) then
CsB_______________Interpolate the SigA piece
                  Call ThreeDInterp ( N_QT_INTP,N_Y_INTP,N_Q_INTP,
     >              YDAS,YDAA, QT_V,yB_V,Q_V,SigS,SigA, ErrSgS,ErrSgA )
C                 This is a K-factor applied to the fake Pert piece
                  FacMatch = 1.d0
CCPY              FacN     = 1./(1.+QT_V/400.)**5
                  FacN     = 1.d0
                Else
CCPY              FacMatch = 1.-QT_V/800.
                  FacMatch = 1.d0
CCPY              FacN     = 1./(1.+QT_V/400.)**5
                  FacN     = 1.d0
                End If
              End If
            Else If (W.LT.2.5d3 .and. QT_V.LT.100.) then ! Tevatron 
CsB___________Interpolate the SigA piece
              Call ThreeDInterp ( N_QT_INTP,N_Y_INTP,N_Q_INTP,
     >            YDAS,YDAA, QT_V,yB_V,Q_V,SigS,SigA, ErrSgS,ErrSgA )
C             The "qqB_delta(1-z)" piece is turned off gradually
CCPY          SigA = SigA*(1.d0-QT_V/100.)**2
              SigA = SigA*1.d0
C             This is a K-factor applied on the fake Pert piece
CCPY          FacMatch = 1.d0-QT_V/4000.
              FacMatch = 1.d0
            End If              !if (W.gt.1d4) ! LHC


          Else If (Boson.Eq.'AA'.or.Boson.Eq.'AG') then
cpn2006     Interpolate SigA (the piece proportional to the 
cpn2006     approximate NLO hard-vertex correction); dismiss SigS

            Call ThreeDInterp ( N_QT_INTP,N_Y_INTP,N_Q_INTP,
     >        YDAS,YDAA, QT_V,yB_V,Q_V,SigS,SigA, ErrSgS,ErrSgA )
cC                 This is a K-factor applied to the fake Pert piece
c                  FacMatch = 1.d0
c                  FacN     = 1./(1.+QT_V/400.)**5
c                Else
c                  FacMatch = 1.-QT_V/800.
c                  FacN     = 1./(1.+QT_V/400.)**5c
c                End If
c              End If
c            End If

          Else  ! for all the other processes

          End If

CCPY Sept 2009: 
CsB_______For the matched resummed case the CSS piece is included in the Y grid
C         This way the CSS piece interpolation does not have to be called.
          SigS = 0.d0

CsB_______Interpolate the Y piece
          iYPSw = 0

          Call ThreeYInterp ( N_QT_INTP,N_Y_INTP,N_Q_INTP,
     >        YDAS,YDAA, QT_V,yB_V,Q_V, SigY, ErrSgY )
          YP0 = SigY(0)
          YP3 = SigY(3)
cdump
C        print*,'QT_V,yB_V,Q_V, SigY ',QT_V,yB_V,Q_V, SigY 

        Else                    !iYPiece.eq.0
CsB_______Interpolate the CSS piece
          Call ThreeDInterp ( N_QT_INTP,N_Y_INTP,N_Q_INTP,
     >      YDAS,YDAA, QT_V,yB_V,Q_V,SigS,SigA, ErrSgS,ErrSgA )
        End If !iYpiece.Eq.1
      End If ! nPart.=!=.11

CCAO  
!     IF(IFLAG.EQ.-997)THEN
!       WT=0.D0
!       RETURN
!     ENDIF

CCPY SIGA CAN BE NEGATIVE for large enough pT when the Y piece is not
C    included or for the delta sigma piece.
c        If ( SigS.LT.0d0 .and. NPart.NE.13) then
c          If ( SIGS.LT.-10.d0*Abs(ErrSgS) .and. SigS.LT.-10.d-2 ) then
c            Write(14,*)' Q_V, QT_V, YB_V, SIGS, ErrSgS '
c            Write(14,*)  Q_V, QT_V, YB_V, SIGS, ErrSgS
c            SIGS=0.D0
c          End If
c        End If

CsB___As of 7/1/97 the pert piece is matched with the CSS at the grid level
CCPY

C ============================
CCPY; Jan 16, 2011
      If (Boson.Eq.'GL' .OR. Boson.Eq.'GG') THEN
        SigS=SigY(0)
        GOTO 900 
      ENDIF
C =============================

CCAO! USE THE INPUT W-BOSON MASS TO REWEIGHT THE READIN GRID
C     WHEN WE CALCULATE THE GRID IN LEGACY, THE COUPLING IS PROPOTIONAL TO 
C                GFERMI*XMW_GRID**2
C     SO WE MULTIPLY THE READIN GRID BY
C                GRID_VALUE*(XMW_IN/XMW_GRID)**2
C     NOTE: 
C          WE ALREADY RESET THE XMW IN THE SUBROUTINE READIN(**), SO NOW
C                  XMW=XMW_IN
C
C      print*, ' XMW, XMW_GRID: ', XMW, XMW_GRID
      IF(Boson(1:2).Eq.'W+'.or.Boson(1:2).Eq.'W-')THEN
        SIGS   = SIGS*(XMW/XMW_GRID)**2
        SIGA   = SIGA*(XMW/XMW_GRID)**2
        SIGY(0)= SIGY(0)*(XMW/XMW_GRID)**2
        SIGY(1)= SIGY(1)*(XMW/XMW_GRID)**2
        SIGY(2)= SIGY(2)*(XMW/XMW_GRID)**2
        SIGY(3)= SIGY(3)*(XMW/XMW_GRID)**2
        SIGY(4)= SIGY(4)*(XMW/XMW_GRID)**2
      ElseIf (Boson(1:2).Eq.'Z0') then
        SIGS   = SIGS*(XMZ/XMZ_GRID)**2
        SIGA   = SIGA*(XMZ/XMZ_GRID)**2
        SIGY(0)= SIGY(0)*(XMZ/XMZ_GRID)**2
        SIGY(1)= SIGY(1)*(XMZ/XMZ_GRID)**2
        SIGY(2)= SIGY(2)*(XMZ/XMZ_GRID)**2
        SIGY(3)= SIGY(3)*(XMZ/XMZ_GRID)**2
        SIGY(4)= SIGY(4)*(XMZ/XMZ_GRID)**2
!       Print*, ' SigS,SigA,SIGY: ', SigS,SigA,SIGY
!       Print*, ' XMZ,XMZ_GRID: ', XMZ,XMZ_GRID
      ENDIF

CsB___Weak coupling
      gWeak = DSqrt(GWEAK2)
      sWeak = DSqrt(SWS)
      cWeak = DSqrt(CWS)

CsB___Running alpha_ew implemented on 11/26/97
      eEM   = Sqrt(4.d0*Pi*PYALEM(Q_V**2))

C_____Vector boson - lepton - lepton couplings: fL and fR
C     and mass and width of the vector boson: YMW_AMP, YGW_AMP
      If (Boson.Eq.'W+' .or. Boson.Eq.'W-') then !==========================
CCAO!   fL = gWeak/2.d0/Sqrt(2.d0)
        fL = DSQRT(WCOUPL)
        fR = 0.d0
        YMW_AMP = XMW
        YGW_AMP = GAMW_IN
        If (GAMW_IN.LT.0.) YGW_AMP = GAMW_NLO
      Else If (Boson.Eq.'HP' .or. Boson.Eq.'HM') then ! --------------------
C_____From hep-ph/9703224 eq. after Eq.(2.3)
CCAO!   fL = gWeak/2.d0
        fL = DSQRT(WCOUPL*2.d0)   !     g/2
        fR = 0.d0
        YMW_AMP = XMW
CCAO!   YGW_AMP = GAMW
        YGW_AMP = GAMW_NLO
      Else If (Boson.Eq.'HZ') then ! --------------------
C_____From hep-ph/9703224 eq. after Eq.(2.3)
        fL = DSQRT(WCOUPL*2.0/(1.d0-SW2_EFF))  !   g/2/cw
        fR = 0.d0
        YMW_AMP = XMZ
CCAO!   YGW_AMP = GAMW
        YGW_AMP = GAMZ_NLO
      Else If (Boson.Eq.'Z0') then ! --------------------
CCAO!   fL = gWeak/2.d0/cWeak * (-1.d0/2.d0 + SWS)
CCAO!   fR = gWeak/2.d0/cWeak *               SWS
CCAO!   YMW_AMP = XMZ
CCAO!   YGW_AMP = GAMZ
        FL = DSQRT(ZCOUPL) * (-1.D0/2.D0 + SW2_EFF)
        FR = DSQRT(ZCOUPL) *               SW2_EFF
        YMW_AMP = XMZ
        YGW_AMP = GAMZ_IN
        If (GAMZ_IN.LT.0.) YGW_AMP = GAMZ_NLO
      Else If (Boson.Eq.'DY') then ! ---------------------------------------
!       fLA = -eEM/2.d0
!       fRA = -eEM/2.d0
!       fLZ = gWeak/2.d0/cWeak * (-1.d0/2.d0 + SWS)
!       fRZ = gWeak/2.d0/cWeak *               SWS
C        fLA = -DSQRT(2.0*ACOUPL)/2.d0
C        fRA =  fLA
C        FLZ =  DSQRT(ZCOUPL) * (-1.d0/2.d0 + SW2_EFF)
C        FRZ =  DSQRT(ZCOUPL) *               SW2_EFF
CCPY  CALL DY_COUPL
        YMW_AMP = XMZ
        YGW_AMP = GAMZ_IN
        If (GAMZ_IN.LT.0.) YGW_AMP = GAMZ_NLO
      Else If (Boson.Eq.'A0') then ! ---------------------------------------
C       fL = fR = -e/2 = -g sW/2
C       fL = -gWeak*sWeak/2.d0
C        fL = -eEM/2.d0
        fL = -DSQRT(2.0*ACOUPL)/2.d0
        fR =  fL
        YMW_AMP = 0.d0
        YGW_AMP = 0.d0
      Else If (Boson.Eq.'AA') then ! ---------------------------------------
C       The charge of the quarks (in units of the positron charge) is already
C       included in Legacy in ZFF2.
C       fL = fR = e/2 = g sW/2
C       fL = gWeak*sWeak/2.d0
CCAO!   fL = -eEM/2.d0
CCAO!   fR = fL
        fL = -DSQRT(2.0*ACOUPL)/2.d0
        fR = fL
        YMW_AMP = 0.d0
        YGW_AMP = 0.d0
      Else If (Boson.Eq.'ZZ') then ! ---------------------------------------
C       The following parts of the coupling are included in Legacy in ZFF2:
CsB_____See: Ohnemus, Owens PRD43 (91) 3627 (3)
C       We write Qq^4 = g+^4/2 + g-^4/2 with g+ = g- = Qq for qQ->AA,
C       and replace g+ and g- for ZZ.
C ....
C        ZFF2(1) = (T3u/Sw_Eff/Cw_Eff-Qu*Tw_Eff)**4/2.d0 +
C     +            (-Qu*Tw_Eff)**4/2.d0
C        ZFF2(2) = (T3d/Sw_Eff/Cw_Eff-Qd*Tw_Eff)**4/2.d0 +
C     +            (-Qd*Tw_Eff)**4/2.d0
C ....
C       With ACOUPL = (e**2)/2
C           CONST = 2.d0/3.d0*PI * ACOUPL
C       Then only an overall factor ~ e^2 is necessary
CCAO!   fL = -eEM/2.d0
CCAO!   fR = fL
        fL = -DSQRT(2.0*ACOUPL)/2.d0
        fR = fL
        YMW_AMP = 0.d0
        YGW_AMP = 0.d0
      Else If (Boson.Eq.'AG') then ! ----------------------------------------
C       All factors are included in WT
C_____
      Else If (Boson.Eq.'ZG') then ! ----------------------------------------
C       All factors are included in WT
C_____
      Else If (Boson.Eq.'H0') then ! ----------------------------------------
C       All couplings and other factors are included in WT
        YMW_AMP = xMHiggs
        YGW_AMP = GAMH
C        Print*, 'xMHiggs, GAMH =',xMHiggs, GAMH
C_____
      ELSE
        Print *, ' Wrong vector boson in WPRES: ', Boson
        Stop
      End If ! ==============================================================
!     Print*, 'YMW_AMP,YGW_AMP =', YMW_AMP,YGW_AMP 
CsB___End setting couplings, etc.

CsB___Set the weight
      If (Boson.Eq.'AA') then ! q qBar or q(Bar) G -> photon photon X =======
C       Parton level Mandelstam variables
        sH =  Q_V**2
        tH = -sH/2.d0*(1.d0-Cos(the_sta))
        uH = -sH/2.d0*(1.d0+Cos(the_sta))

C         Rescaling the C function O(alpha_s) qqB delta(1-z) contribution
C         - returned in SIGA by Legacy.
C         In SIGS the full contribution was returned by Legacy but with a
C         kinematics independent delta(1-z) piece.  This latter is in SIGA.
C         The rescaling corrects the kinematics dependence.
C         Depending on the version of Legacy, the rescaling is different.
cpn2006 Compute the correction to the (matched) resummed cross section
cpn2006 due  to angular dependence of the NLO hard-vertex coefficient
        If (iLV.Eq.1) then      ! v2.2 and older
CCPY          DeNom = Pi**2 - 8.d0
          PRINT*,' GENERATE THE GRID FILES USING C++ CODE'
          CALL EXITGracefully
        Else                    ! v2.3 and later
          DeNom = 4.d0*h1L/CF
        End If
cpn2006 Compute the correction to W-A from the angle-dependent 
c      hard vertex at Q_T below the switching point from W+Y to PERT

cdump
c      print*,' qt_v,sigy(0),CJJ1XFJ_1_AA(sH,tH,uH),DeNom,h1L =',
c     > qt_v,sigy(0),CJJ1XFJ_1_AA(sH,tH,uH),DeNom,h1L
      
      
CCPY Sept 2009:        if (abs(sigs).gt.1e-5)
CCPY        if (abs(SIGY(0)).gt.1e-5)
CCPY     >    SIGS_AA = SIGS-SIGA + CJJ1XFJ_1_AA(sH,tH,uH)/DeNom*SIGA

        SIGS_AA = SIGS-SIGA + CJJ1XFJ_1_AA(sH,tH,uH)/DeNom*SIGA

        If (nPart.Eq.13) then ! Delta-sigma piece
          Sigma0 = Pi**2*PyAlEM(Q_V**2)*4.d0/3.d0
          h1v = CF*CJJ1XFJ_1_AA(sH,tH,uH)/4.d0
          KinCorr = 0
          If (KinCorr.Eq.0) then
            X1 = Q_V/ECM*DEXP(yB_V)
            X2 = Q_V/ECM/DEXP(yB_V)
          Else If (KinCorr.Eq.1) then
            X1 = Sqrt(Q_V**2+QT_V**2)/ECM*DEXP(yB_V)
            X2 = Sqrt(Q_V**2+QT_V**2)/ECM/DEXP(yB_V)
          End If
          Q = Q_V ! factorization scale = Q !!! hard wired !!!
          Call YPartons
          
          Call PartonMomentumLumi(qQBUL,qQBDL,QBqUL,QBqDL,
     ,                          qGUL,qGDL,GqUL,GqDL,GGL)
          LParton = (qQBUL + qQBDL + QBqUL + QBqDL)/x1/x2
!         Print*, iBeam,qQBUL/X1/X2,qQBDL/X1/X2,QBqUL/X1/X2,QBqDL/X1/X2,
!    ,                  qGUL/X1/X2,qGDL/X1/X2,GqUL/X1/X2,GqDL/X1/X2
!         Print*, SigS,Sigma0/sH * 2.d0*alphaS(Q_V)/Pi*h1v,
!    ,            Sigma0,sH,alphaS(Q_V),Pi,h1v
!         Stop
          SigS_AA = SigS + 
     >      Sigma0/sH*2.d0*alphaS(Q_V)/Pi*(h1v-h1L)*LParton
        End If

        WT =
C           Identical particle final state and spin average factors
     &      1.d0/2.d0/2.d0/2.d0*
C           Amplitude squared, flux and phase space factors
     &      32.d0*(fL**2+fR**2)*(tH/uH+uH/tH)*(SigS_AA+SigY(0))

CsB___An attempt to calculate the PRT weight in the same run.
C     Cannot be done. Vegas calls the different phase spaces. 
C     Before they are called, # of dimensions, etc. has to be set.
C     These parameters are different for resummed and fixed order. 
C      Call QQBZZG(WTPrt)
C      Print*, ' WT, WtPrt = (Stopping here)', WT, WtPrt
C      Stop

CsB_____To obtain the real Y piece from the fake one we use the K-factor
C       (only to PLOT the matched QT dist'n at the LHC)
        KFac = 1.d0

CCPY Sept 2009
CCPY!
CCPY!        If (iYPiece.Eq.1 .And. W.GT.1.d4 .And. QT_V.GT.50. .And. 
CCPY!     .      QT_V.LT.400.) then
CCPY!          KFac = Exp(-(QT_V-50.)*Log(.21/.003)/(150.-50.))
CCPY!          WT = KFac*WT
CCPY!        EndIf
CCPY!
      Else If (Boson.Eq.'ZZ') then ! q qBar or q(Bar) G --> Z0 Z0 G ----------
        Mz2 = XMZ**2
        Mz4 = XMZ**4

C_______Beta of Z pair
        x  = Sqrt(1.d0-4.d0*Mz2/Q_V**2)

C________Testing the Mz -> 0 (di-photon) case.
C        Mz2 = 1.d0
C        Mz4 = 1.d0
C        x   = Sqrt(1.d0 - 4.d0*Mz2/Q_V**2)
C        Print*, ' Mz2 =', Mz2

C_______Parton level Mandelstam variables
        sH =  Q_V**2
        tH = -sH/2.d0*(1.d0 - x*Cos(the_sta)) + Mz2
        uH = -sH/2.d0*(1.d0 + x*Cos(the_sta)) + Mz2

        If (iLV.Eq.1) then      ! v2.2 and older
CCPY          DeNom = Pi**2 - 8.d0
          PRINT*,' GENERATE THE GRID FILES USING C++ CODE'
          CALL EXITGracefully
        Else                    ! v2.3 and later
          DeNom = 4.d0*h1L/CF
        End If

C Rescaling the C function O(alpha_s) qqB delta(1-z) contribution
C - returned in SIGA by Legacy.
C In SIGS the full contribution was returned by Legacy.
C
C Comparing PRD43('91)3632 (12) with PRD46('92)2018 (11) we see that
C 1/(2 Pi s) * F^{virt}(v) of the diphoton case corresponds to
C 8 (F(t,u)+F(u,t)) of the diZ case.
C In other words:
C F^{virt}(v) of the diphoton case corresponds to
C 16 Pi s (F(t,u)+F(u,t)) of the diZ case.
C

cdump
c      test_1=(CJJ1XFJ_1_ZZ(sH,tH,uH,x)+CJJ1XFJ_1_ZZ(sH,uH,tH,x))*CF/4
c      print*,'qt_v,sigy(0),test_1,DeNom,h1L'
c      print*,qt_v,sigy(0),test_1,DeNom,h1L
      
C Then the CSS piece corrected with the qqB_delta(1-z) piece is:

CCPY Sept 2009:
CCPY!        SigR_ZZ = SigS-SIGA + 
CCPY!     +            (CJJ1XFJ_1_ZZ(sH,tH,uH,x) + CJJ1XFJ_1_ZZ(sH,uH,tH,x))/
CCPY!     /             (Pi**2 - 8.d0)*SIGA +
CCPY!     +            SigY(0)*FacMatch
CCPY!

        SigR_ZZ = SigS-SIGA + 
     +            (CJJ1XFJ_1_ZZ(sH,tH,uH,x) + CJJ1XFJ_1_ZZ(sH,uH,tH,x))/
     /             Denom*SIGA +
     +            SigY(0)*FacMatch

C_______The last term plays no role when there's no Y piece.
C       On the other hand it gives a possibility for just DRAWING a smooth,
C       matched curve. Same for FacN.

        SigS_ZZ = FacN*SigR_ZZ
        
C Testing the Mz -> 0 case
C        Print*, CJJ1XFJ_1_AA(sH,tH,uH), CJJ1XFJ_1_ZZ(sH,tH,uH,x),
C     &    CJJ1XFJ_1_ZZ(sH,uH,tH,x)
C        Print*, SIGS-SIGA,
C     &    CJJ1XFJ_1_ZZ(sH,tH,uH,x) + CJJ1XFJ_1_ZZ(sH,uH,tH,x)
C        Stop
C This test showes that CJJ1XFJ_1_ZZ(sH,tH,uH,x) + CJJ1XFJ_1_ZZ(sH,uH,tH,x)
C for Mz -> 0 reproduces CJJ1XFJ_1_AA(sH,tH,uH) up to the fifth digit.

        WT =
C_______Identical particle final state and spin average factors
     &      1.d0/2.d0/2.d0/2.d0*
C_______Amplitude squared (see PRD43('91)3626 Eq.(4)):
     &      32.d0*(fL**2+fR**2)*
     *            (tH/uH + uH/tH - Mz4*(1.d0/tH**2 + 1.d0/uH**2) +
     +             4.d0*Mz2*sH/tH/uH )*SigS_ZZ


      ElseIf (Boson.Eq.'AG') then ! G G --> photon photon X ------------------
        sH =  Q_V**2
        tH = -sH/2.d0*(1.d0-Cos(the_sta))
        uH = -sH/2.d0*(1.d0+Cos(the_sta))
        vin = -tH/sH

        Call ggAA_A(sH,tH,uH) ! only necessary if ggAA_A2 is used below

        If (iLV.Eq.1) then      ! < v2.2 or Asym,DSig piece
           SIGA = 0.d0          !no second dipole term
           h1L = 1.
           h1v = h1L
        Else
           If (nPart.Eq.12) then
              h1v = h1L         ! approximate h1 - no need for by Asym
           Else
              h1v = h1ggAA(vin,Nf)
           End If
        End If                  !iLV
        
        If (nPart.Eq.12) then   ! Asymptotic piece
           SigS_AG = SIGS
           SigA_AG = SIGA
        Else                    ! resummed piece
cpn2006 Compute the correction to W-A from the angle-dependent 
c      hard vertex at Q_T below the switching point from W+Y to PERT
CCPY!*******************************************************
CCPY!cdump
CCPY!      print*,' npart =',npart
CCPY!      print*,' SIGA,SIGS,sigy(0) =', SIGA,SIGS,sigy(0) 
CCPY!      print*,' qt_v,h1v,h1L =', qt_v,h1v,h1L 
CCPY!
CCPY!*******************************************************
CCPY Sept 2009:        if (abs(sigs).gt.1e-5)
CCPY!          if (abs(SIGY(0)).gt.1e-5) 
CCPY!     >      SigS_AG = SIGS-SIGA + h1v/h1L*SIGA ! see comments at Boson.Eq.'AA'

          SigS_AG = SIGS-SIGA + h1v/h1L*SIGA ! see comments at Boson.Eq.'AA'
          SigA_AG = 0.d0         

        End If

        If (nPart.Eq.13) then ! Delta-sigma piece
          Sigma0 = 2.d0*Pi/8.d0/8.d0*Pi**2*(alphaS(Q_V)/Pi)**2
          KinCorr = 0
          If (KinCorr.Eq.0) then
            X1 = Q_V/ECM*DEXP(yB_V)
            X2 = Q_V/ECM/DEXP(yB_V)
          Else If (KinCorr.Eq.1) then
            X1 = Sqrt(Q_V**2+QT_V**2)/ECM*DEXP(yB_V)
            X2 = Sqrt(Q_V**2+QT_V**2)/ECM/DEXP(yB_V)
          End If
          Q = Q_V ! factorization scale = Q !!! hard wired !!!
          Call YPartons
          Call PartonMomentumLumi(qQBUL,qQBDL,QBqUL,QBqDL,
     ,                          qGUL,qGDL,GqUL,GqDL,GGL)
          LParton = GGL/x1/x2
!         Print*, iBeam,GGL/X1/X2
!         Print*, SigS,Sigma0/sH * 2.d0*alphaS(Q_V)/Pi*h1v,
!    ,            Sigma0,sH,alphaS(Q_V),Pi,h1v
!         Stop
          SigS_AG = SigS 
     >      + Sigma0/sH*2.d0*alphaS(Q_V)/Pi*(h1v-h1L)*LParton
        End If

        WT =
CsB         Legacy includes in SigS: Pi*alpha_s^2/32
C___________Identical particle final state and spin average factors
     &      1.d0/2.d0/2.d0/2.d0*
C___________Amplitude squared
!    &      2.d0*(gWeak*sWeak)**4/(16.0*Pi**2)*
     &      2.d0*eEM**4/(16.0*Pi**2)*
!!    &      2.d0*DSqrt(2.0*ACOUPL)**4/(16.0*Pi**2)*
CsB_________Quark masses are zero in loop:
!    &      (11.0/9.0)**2*XMAT_PHOT(sH,tH,uH)*(SigS_AG+SigY(0))
CsB_________t and b quark masses are non-zero in loop:
c     &      GGAA(sH,tH,uH)*(SigS_AG+SigY(0))
CsB_________Implements Eq(3) from hep-ph/0206194 
     &      (ggAA_A2(sH,tH,uH)*(SigS_AG+SigY(0)) +
     +       ggAA_D2(sH,tH,uH)* SigA_AG )
cpn2005!     Alternative computations using equivalents of subroutines
cpn2005!      ggAA_A2 and ggAA_D2 (fMLO2 and fMLO2prime) that utilize
cpn2005!     matrix elements used in the h1ggAA subroutine      
cpn2005     &      (fMLO2(vin,Q_V)*(SigS_AG+SigY(0)) +
cpn2005     +       fMLO2prime(vin,Q_V,phi_sta)* SigA_AG )
cpn2005     >       0)
CCPY!******************************************************************************
CCPY!
CCPY!       Print*, 'sH,tH,uH=',sH,tH,uH
CCPY!       Print*, 'SIGA_AG =',SIGA_AG
CCPY!       Print*, 'GGAA(sH,tH,uH)=',GGAA(sH,tH,uH)
CCPY!       Print*, 'XMAT_PHOT=',(11.0/9.0)**2*XMAT_PHOT(sH,tH,uH)
CCPY!       Print*, 'ggAA_A2(sH,tH,uH)=',ggAA_A2(sH,tH,uH)
CCPY!       Write(*,'(A,2f12.6)')'ratios=',GGAA(sH,tH,uH)/ggAA_A2(sH,tH,uH),
CCPY!     ,               (11.0/9.0)**2*XMAT_PHOT(sH,tH,uH)/ggAA_A2(sH,tH,uH)
CCPY!       Print*, 'ggAA_D2(sH,tH,uH)=',ggAA_D2(sH,tH,uH)
CCPY!       Print*, '---'
CCPY!       Stop
CCPY!******************************************************************************
CCPY!
      ElseIf (Boson.Eq.'ZG') then ! G G --> photon photon X ------------------
        Mz2 = XMZ**2 ! Z mass squared
        x2 = 1.d0 - 4.d0*Mz2/Q_V**2
        If (x2.Lt.0.d0) then
          If (x2.LT.-1.d-2) Print*, ' x2 = ', x2
          x = 1.d-6
        Else
          x = Sqrt(x2)
        End If

C_______Parton level Mandelstam variables
        sH =  Q_V**2
        tH = -sH/2.d0*(1.d0 - x*Cos(the_sta)) + Mz2
        uH = -sH/2.d0*(1.d0 + x*Cos(the_sta)) + Mz2

CCPY Sept 2009
CCPY! First, the function ggAA only holds for gg -> AA, 
CCPY! but not for gg -> ZZ because of different external leg masses. 
CCPY! Furthermore, the couplings of (internal) quarks to Z are different 
CCPY! from those to photon. I did not see either of 
CCPY! these being implemented in the code. Hence,
CCPY! I concluded that the gg -> ZZ calculation 
CCPY! has not been correctly completed inside this code. 

        WT =
C_______Identical particle final state and spin average factors
     &      1.d0/2.d0/2.d0/2.d0*
C_______Amplitude squared
     &      2.d0*eEM**4/(16.0*Pi**2)*
C_______Compared to the di-photon production the f-VB-f couplings have to be
C       changed for Z Z production.
     &      GGAA(sH,tH,uH)*(SigS+SigY(0))

      ElseIf (Boson.Eq.'H0') then ! G G --> H0 -------------------------------
CsB_____The general formula for fixed order scalar production and decay:
C       sigma = 4 Pi Gamma(production)*Scalar Propagator*Gamma(decay)

CsB_____The following is implemented in res.for:
C        ELSEIF(TYPE_V.EQ.'H0') THEN
CC         THIS IS THE (\sigma_0)*(Q_V)^2 IN MY NOTES
CC         In SUBROUTINE STANDARD: HCOUPL=R2*GMU*PI/576.0D0 =>
CC         Const = Sqrt(2) Gmu Pi/576 * alpha_s^2/pi^2 Q^2 =
CC                 Sqrt(2) Gmu alpha_s^2 Q^2 / (576 Pi)
C          CONST=HCOUPL*(ALPI(S_Q)*Q_V)**2
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
CCPY!        If (QT_V.GT.130. .and. QT_V.LT.400.) 
CCPY!     >    KFac = 1.d0 - .5*(QT_V-130.)/(400.-130.)

CsB     Electroweak K-factor could be added from hep-ph/0404071


ccpy===================================================
CCPY THIS WAS IMPLEMENTED FOR THE RESBOS VERSION WITHOUT INCLUDING 
C MODULES IN HVVRESPH AND WPRES_HVV (PRIOR TO Jan 2006)
CCPY        If (hDecayMode.Eq.'ZZ') then
CCPY!********************************************************************************
CCPY!          WT =
CCPY!C_________Effective coupling of the Higgs to gluons in the heavy top quark limit
CCPY!C         see Eq.(13) hep-ph/9611272 by Kramer-Laenen-Spira
CCPY!!     *         (1.0d0+KappaSoft-KappaNLO) *
CCPY!C_________Factor which takes into account the quark masses in the g g -> H 
CCPY!C         fermion loop
CCPY!     *         ENFS(Q_V) *
CCPY!C_________d(sigma[G G -> H])/d(Q^2)d(y)d(Q_T) from Legacy
CCPY!     *         (SigS + KFac*SigY(0)) *
CCPY!C_________Higgs propagator
CCPY!     *         1.d0/((Q_V**2-YMW_AMP**2)**2 + (YMW_AMP*YGW_AMP)**2) *
CCPY!C_________H -> Z Z amplitude squared and averaged
CCPY!     *         Sqrt(2.d0)*GFermi *
CCPY!     *         (12.d0*xMz**4 - 4.d0*(Q_V*xMz)**2 + Q_V**4) *
CCPY!C_________Extra 2*Q^2 factor to balance the one which is in the phase space
CCPY!     *         2.0*Q_V**2 *
CCPY!C_________symmetry factor for identical particles in final state
CCPY!     *         1.0/2.0
C
CCPY!C_________If they decay, then include the Z -> l l' branching ratio for both Z's
CCPY!          If (I_DECAY.EQ.1) WT = (.0336)**2 * WT
CCPY        Else If (hDecayMode.Eq.'AA') then
ccpy===================================================


        If (hDecayMode.Eq.'AA') then
CCPY!********************************************************************************
CCPY June 2009
C   Definition of the VB propagator
c       Print*, ' YMW_AMP,YGW_AMP = ', YMW_AMP,YGW_AMP
        If (YMW_AMP.Eq.0.e0) then
          Propagator = (Q_V**2-YMW_AMP**2)**2+(YMW_AMP*YGW_AMP)**2
        Else
          Propagator = 
     &       (Q_V**2-YMW_AMP**2)**2+Q_V**4*(YGW_AMP/YMW_AMP)**2
        End If

            WT = 
C_________Factor which takes into account the quark masses in the g g -> H 
C         fermion loop
     *         ENFS(Q_V) *
C_________d(sigma[G G -> H])/d(Q^2)d(y)d(Q_T) from Legacy
     *         (SigS + KFac*SigY(0)) *
C_________Higgs propagator
     *         1.d0/Propagator  *
C_________Factor to balance |M|^2
     *         32.d0*Pi*Q_V**2/Q_V * YGW_AMP *
C_________Extra Q^2 factor to balance the one which is in the phase space
     *         2.0*Q_V**2 * 
C_________symmetry factor for identical particles in final state
     *         1.0/2.0
     
C_________Br(H -> photon photon)
            If (I_DECAY.EQ.1) WT=WT* BrH2AA(Q_V)    
                      
        End If ! hDecayMode.Eq.'AA'

      Else If (Boson.Eq.'HP' .or. Boson.Eq.'HM'
     >.or. Boson.Eq.'HZ') then ! p p(B) -> W+- h or Z h -----

C      If (SigS.lt.0.d0 .or. SigA.lt.0.d0) then
C        Print*, 'SigS,SigA, QT_V,yB_V in WT: ', SigS,SigA, QT_V,yB_V
C      End If

C_______Comparing PRD56('97)5558 Eq.(1) and hep-ph/9703224 Eq.(2.1)
C       (plus see my notes on hep-ph/9703224)

        pStar = 
     =   Sqrt((Q_V**2-(YMW_AMP+xMHiggs)**2)*
     >   (Q_V**2-(YMW_AMP-xMHiggs)**2))/
     /   (2.d0*Q_V)

CCPY Nov. 2009; Corrrect an overall factor of 4
C_______Note: beta is included in the phase-space WRESPH
        WT =
C________propagator:
     *   1.D0/((Q_V**2-YMW_AMP**2)**2+(YMW_AMP*YGW_AMP)**2) *
C________factors for W+- production
     *   4.d0*Q_V**4 *
C________couplings
     *   ( (fL**2.d0 + fR**2.d0) * (
C________angular dependence:
     *      YMW_AMP**2/Q_V**2*
     *      (2.d0 + pStar**2/YMW_AMP**2*Sin(the_sta)**2) * ! H(s)
     *      (SigS + SigY(0)) ))

      Else If (Boson.Eq.'DY') then 
CCPY Oct 2008, INCLUDE 'DRELL-YAN PAIR 'DY' PROCESS
C        TESTING=.TRUE.
        TESTING=.FALSE.

        Sign1 = 1.d0

        Propagator = 
     &       (Q_V**2-YMW_AMP**2)**2+Q_V**4*(YGW_AMP/YMW_AMP)**2

CBIN test the running width effects
c        Propagator = 
c     &       (Q_V**2-YMW_AMP**2)**2+YMW_AMP**4*(YGW_AMP/YMW_AMP)**2

c test the interference
        ZZ_PROP = Q_V**4/Propagator 
        AA_PROP = 1.D0

c        ZZ_PROP = 0.d0
c        AA_PROP = 0.D0
        ZA_PROP = Q_V**2*(Q_V**2-YMW_AMP**2)/Propagator

        IF(TESTING) THEN
          ALFA_Q2=EW_ALFA
          ZZ_PROP = ZZ_PROP
          ZA_PROP = 0.d0
          AA_PROP = 0.d0
        ELSE
C        EW_ALFA = 1.0D0/137.0359895D0
          ALFA_Q2=PYALEM(Q_V**2)
          ZA_PROP = ZA_PROP*(ALFA_Q2/EW_ALFA)
          AA_PROP = AA_PROP*(ALFA_Q2/EW_ALFA)**2
        ENDIF

        IF(KFAC_ANGFUNC.EQ.1) THEN
CJI Add in k-factors for pert piece from 1708.00008
                if(QT_V .lt. 11.4) then
                    kloc = 1
                else if(QT_V .lt. 18.5) then
                    kloc = 2
                else if(QT_V .lt. 25.5) then
                    kloc = 3
                else if(QT_V .lt. 32.6) then
                    kloc = 4
                else if(QT_V .lt. 40.4) then
                    kloc = 5
                else if(QT_V .lt. 50.2) then
                    kloc = 6
                else if(QT_V .lt. 63.9) then
                    kloc = 7
                else if(QT_V .lt. 85.4) then
                    kloc = 8
                else if(QT_V .lt. 105.0) then
                    kloc = 9
                else if(QT_V .lt. 132.0) then
                    kloc = 10 
                else if(QT_V .lt. 173.0) then
                    kloc = 11
                else if(QT_V .lt. 253.0) then
                    kloc = 12
                else
                    kloc = 13
                end if
        End If   ! KFAC_ANGFUNC.EQ.1

CCPY Oct 2019
       IF(KFAC_ANGFUNC.EQ.0) THEN
         ANGFUNC_EVEN = (
     &         (1.d0 + Cos(the_sta)**2) *      (SigS + SigY(0)) +       ! L0
     &         (1.d0-3.d0*Cos(the_sta)**2)/2.d0*SigY(2) +      ! A0
     &         Sin(2.d0*the_sta)*Cos(phi_sta) *SigY(1) +       ! A1
     &         Sin(the_sta)**2*Cos(2.d0*phi_sta)/2.d0*SigY(2)) ! A2
        ANGFUNC_ODD = (
     &         2.d0*Cos(the_sta) *             (SigA + SigY(3)) +    ! A3
     &         Sin(the_sta)*Cos(phi_sta) *      SigY(4))   ! A4

       ELSEIF(KFAC_ANGFUNC.EQ.1) THEN
         ANGFUNC_EVEN = (
     &         (1.d0 + Cos(the_sta)**2) *      (SigS + SigY(0)) +       ! L0
     &         (1.d0-3.d0*Cos(the_sta)**2)/2.d0*SigY(2)*a0k(kloc) +     ! A0
     &         Sin(2.d0*the_sta)*Cos(phi_sta) *SigY(1)*a1k(kloc) +      ! A1
     &         Sin(the_sta)**2*Cos(2.d0*phi_sta)/2.d0*SigY(2)*a2k(kloc))! A2
        ANGFUNC_ODD = (
     &         2.d0*Cos(the_sta) *             (SigA + SigY(3)) +    ! A3
     &         Sin(the_sta)*Cos(phi_sta) *      SigY(4)*a4k(kloc))   ! A4
       ELSEIF(KFAC_ANGFUNC.EQ.2) THEN
        ANGFUNC_EVEN = (SigS + SigY(0)) * (
     &         (1.d0 + Cos(the_sta)**2) +                               ! L0
     &         (1.d0-3.d0*Cos(the_sta)**2)/2.d0*SigY(5) +               ! A0
     &         Sin(2.d0*the_sta)*Cos(phi_sta) *SigY(1) +                ! A1
     &         Sin(the_sta)**2*Cos(2.d0*phi_sta)/2.d0*SigY(2))          ! A2
        ANGFUNC_ODD = (
     &         2.d0*Cos(the_sta) *             (SigA + SigY(3)) +    ! A3
     &         Sin(the_sta)*Cos(phi_sta) * (SigS+SigY(0))*SigY(4))   ! A4
       ELSE
         WRITE(*,*) 'THIS KFAC_ANGFUNC NOT IMPLEMENTED: ',KFAC_ANGFUNC
         WRITE(16,*) 'THIS KFAC_ANGFUNC NOT IMPLEMENTED: ',KFAC_ANGFUNC
         CALL EXITGracefully
       ENDIF

       COEFF_EVEN = AA_PROP*CAA_E + 2.D0*ZA_PROP*CZA_E + ZZ_PROP*CZZ_E
       COEFF_ODD = AA_PROP*CAA_O + 2.D0*ZA_PROP*CZA_O + ZZ_PROP*CZZ_O

C_yfu print out the propagator
c       print*,'MZ = ',YMW_AMP,' gamZ = ',YGW_AMP
c       print*,'Q = ', Q_V,' AA_PROP = ',AA_PROP, ' ZA_PROP = ',ZA_PROP,
c     >' ZZ_PROP = ',ZZ_PROP
c       print*,'CAA_E = ',CAA_E,' CZA_E = ',CZA_E,' CZZ_E = ',CZZ_E
c       print*,'CAA_O = ',CAA_O,' CZA_O = ',CZA_O,' CZZ_O = ',CZZ_O
c       print*,'COEFF_EVEN = ',COEFF_EVEN
c       print*,'COEFF_ODD = ',COEFF_ODD
C_______ADDITIONAL FACTORS:
C Here, 4 for 4*Q^2, and 2 for ZFF
        FACTOR=2.D0*(PI/3.D0)*(4.D0*2.D0)
c       print*,'COEFF_EVEN(times factor) = ',COEFF_EVEN*FACTOR
c       print*,'COEFF_ODD(times factor) = ',COEFF_ODD*FACTOR

        WT = FACTOR*( COEFF_EVEN*ANGFUNC_EVEN + 
     >Sign1*COEFF_ODD*ANGFUNC_ODD)
      
 
! cdump
!       print*,' q_v,qt_v,yb_v,sigs,siga,sigy(0)'
!       print*,q_v,qt_v,yb_v,sigs,siga,sigy(0)
 !      print*,' YGW_AMP, YMW_AMP =', YGW_AMP, YMW_AMP
!       print*,'ALFA_Q2,EW_ALFA, Q_V'
!       print*,ALFA_Q2,EW_ALFA, Q_V
!       print*,' CAA_E,CAA_O,CZA_E,CZA_O,CZZ_E,CZZ_O,IFLAV '
!       print*,CAA_E,CAA_O,CZA_E,CZA_O,CZZ_E,CZZ_O,IFLAV 
!       print*,'ZZ_PROP,AA_PROP,ZA_PROP,ALFA_Q2,FACTOR'
!       print*,ZZ_PROP,AA_PROP,ZA_PROP,ALFA_Q2,FACTOR
!       print*,'ANGFUNC_EVEN,ANGFUNC_ODD,COEFF_EVEN,COEFF_ODD'
!       print*,ANGFUNC_EVEN,ANGFUNC_ODD,COEFF_EVEN,COEFF_ODD


      Else ! p p(B) -> (single) 'W+','W-','Z0' or 'A0'
CsB_____Matching now is determined from the resummed and the Y piece grids.
C
C       From P. Arnold's code for W and Z in the range 20 GeV < QT < 60 GeV
C       the O(al_s^2)/O(al_s) = K factors (for p pBar, E_cm = 1.8 TeV) are:
C       1.45 > K > 1.49.

CsB This fixes instabilities in the A0 Y piece temporarily.
        If (Boson.Eq.'A0' .And. QT_V.LT.QT_V_Mat) then
c          SigY(0) = 0.d0
cc          SigY(1) = 0.d0
c          SigY(2) = 0.d0
cc          SigY(3) = 0.d0
cc          SigY(4) = 0.d0
C When SigY(2) is turned on there are spikes in the E_l distribution.
C SigY(0) causes spikes in most distributions.
C Small spikes in the cut Q_T distribution when SigY(1,3,4) turned on.
          Do ii = 0,4
            SigY(ii) = 0.d0
          End Do
        End If

C      If (SigS.lt.0.d0 .or. SigA.lt.0.d0) then
C        Print*, 'SigS,SigA, QT_V,yB_V in WT: ', SigS,SigA, QT_V,yB_V
C      End If

        Sign1 = 1.d0
        If (Boson.Eq.'W-') Sign1 = -1.d0
CsB General fromula from notes of 10/95 p7.
C   Extended with the Y piece. See PLB355('95)553 (16).

C   Definition of the VB propagator
C       Print*, ' YMW_AMP,YGW_AMP = ', YMW_AMP,YGW_AMP
C       Print*, ' fL fR = ', fL,fR
        If (YMW_AMP.Eq.0.e0) then
          Propagator = (Q_V**2-YMW_AMP**2)**2+(YMW_AMP*YGW_AMP)**2
        Else
          Propagator = 
     &       (Q_V**2-YMW_AMP**2)**2+Q_V**4*(YGW_AMP/YMW_AMP)**2
        End If

       if(KFAC_ANGFUNC.ne.2) then
        WT =
     &   4.d0*Q_V**4 / Propagator *
C Couplings and angular dependence:
     &   (   (fL**2.d0 + fR**2.d0) * (
     &         (1.d0 + Cos(the_sta)**2) *      (SigS + SigY(0)) +    ! L0
     &         (1.d0 - 3.d0*Cos(the_sta)**2)/2.d0 *     SigY(2) +    ! A0
     &         Sin(2.d0*the_sta)*Cos(phi_sta) *         SigY(1) +    ! A1
     &         Sin(the_sta)**2*Cos(2.d0*phi_sta)/2.d0 * SigY(2) ) +  ! A2
     & Sign1*(fL**2.d0 - fR**2.d0) * (
     &         2.d0*Cos(the_sta) *             (SigA + SigY(3)) +    ! A3
     &         Sin(the_sta)*Cos(phi_sta) *              SigY(4) ) )  ! A4
        else
            WT =
     &   4.d0*Q_V**4 / Propagator * ((SigS + SigY(0)) *
C Couplings and angular dependence:
     &   (   (fL**2.d0 + fR**2.d0) * (
     &         (1.d0 + Cos(the_sta)**2)                         +      ! L0
     &         (1.d0 - 3.d0*Cos(the_sta)**2)/2.d0 *     SigY(5) +      ! A0
     &         Sin(2.d0*the_sta)*Cos(phi_sta) *         SigY(1) +      ! A1
     &         Sin(the_sta)**2*Cos(2.d0*phi_sta)/2.d0 * SigY(2) ) +    ! A2
     & Sign1*(fL**2.d0 - fR**2.d0) * (
     &         Sin(the_sta)*Cos(phi_sta) *              SigY(4) ) )    ! A4
     &   +
     & Sign1*(fL**2.d0 - fR**2.d0) * (
     &         2.d0*Cos(the_sta) * (SigA + SigY(3)) )) ! A3
        endif

C yfu test
c       print*,WT,Q_V,Propagator,fL,fR,SigS,SigY(0)
CCPY Add the running effect from QED interaction
        IF (Boson.Eq.'A0') THEN
          ALFA_Q2=PYALEM(Q_V**2)
          WT=WT*(ALFA_Q2/EW_ALFA)**2
        ENDIF

C yfu test
c       if(Boson.eq.'A0') then
c         WT = 3.d0 / 8.d0 * 
c     & ( (1.d0 + Cos(the_sta)**2)*(SigS +SigY(0)) +
c     &         Sin(2.d0*the_sta)*Cos(phi_sta) *         SigY(1) +  !A1
c     &         Sin(the_sta)**2*Cos(2.d0*phi_sta)/2.d0 * SigY(2) )  !A2
c       endif

CsB The wrong way to do it, to test (1 + Cos[th])^2 distribution:
C        WT =
C     &   4.d0*Q_V**4 *
C VB Propagator:
C     &   1.D0/((Q_V**2-YMW_AMP**2)**2+(YMW_AMP*YGW_AMP)**2) *
C Couplings and angular dependence:
C     &   ( (fL**2.d0 + fR**2.d0) * (
C     &       (1.d0 + Cos(the_sta)**2) *      (SigS + SigY(0)) ) +    ! L0
C     &     (fL**2.d0 - fR**2.d0) * (
C     &       2.d0*Cos(the_sta) *             (SigS + SigY(0)) ) )    ! A3

      End If ! p p(B) -> (single) 'W+','W-','Z0' or 'A0'

CsB   Test the weight:
!     Print*, ' WT,Propagator,Q_V: ', WT,Propagator,Q_V
!     Print*, ' fL,fR:', fL,fR
!     Print*, ' Cos(the_sta),Cos(phi_sta): ',Cos(the_sta),Cos(phi_sta)
!     Stop

CsB___Weight is set +++======================================================

CsB___No need to prevent negative weigths passed to Vegas.
C      If (WT.LE.0.d0) WT = 1.d-9

CsB___For number of Q points in the grid = 1.
      If (iOneQ.Eq.1) then
CCPY
        If (Boson.Eq.'W+' .or. Boson.Eq.'W-' 
     >.or. Boson.Eq.'Z0'.or. Boson.Eq.'DY') then
CsB The numerical factor below is only for V = W,Z and Q = M_V GeV.
C   It also assumes Delta Q (set in ReadIn) 0.003.
C   The empirical factor 0.980925 was put in by hand.
CCPY June 2009:          WT = 0.980925*Pi/2.d0*YGW_AMP/0.003d0*WT
          WT = Pi/2.d0*YGW_AMP/0.003d0*WT
        End If
      End If

      RETURN

  900 CONTINUE

C FOR BOSON = GL or GG -----------------------------------------------------
CCPY 'GL' PROCESS: (This is to estimate the effect from initial state 
C radiation using Drell-Yan Sudakov factor in the calculaiton. 
C  q qbar -> GL --> t tbar (partial RESUMMED, initial state)
CCPY 'GG' PROCESS: (This is to estimate the effect from initial state 
C radiation using gg -> H0 Sudakov factor in the calculaiton. 
C  gluon gluon --> t tbar (partial RESUMMED, initial state)
C

      IF(ABS(XMASS(1)-XMASS(2)).GT.1.D-8) THEN
        WRITE(*,*) ' XMASS(1) .NE. XMASS(2) IN WPRES '
        CALL EXITGracefully
      ENDIF

      S_Q=Q_V
CSM
      ALFAS=ALPHAS(S_Q)
      IF(Boson.Eq.'GL') THEN ! ---------------------------------------------
CCPY Jan 16, 2011
C In leading order QCD, the asymmetry of top: ASY_TC = 0.d0.
       
       WT =
     &   4.0*PI*ALFAS*2.0*SIGS*
     &   ( 1.d0 + Cos(the_sta)**2+
     &     4.0*XMASS(1)**2/Q_V**2*(1.d0 - Cos(the_sta)**2)+ 
     &     ASY_TC*Cos(the_sta) )
     
      ELSEIF(Boson.Eq.'GG') THEN ! -----------------------------------------
       ETOP=SQRT(Q_V**2)/2.0
       BTOP=SQRT(1.-XMASS(1)**2/ETOP**2)
       CTHE=COS(THE_STA)
       STHE=SQRT(1.-CTHE**2)
       BCTH=BTOP*CTHE
       WT=2.*(XMASS(1)/ETOP)**2*((1.+BTOP)**2+(1.-BTOP)**2)
     $  +2.*(BTOP*STHE)**2*((1.+CTHE)**2+(1.-CTHE)**2)
     $  +4.*(XMASS(1)/ETOP*BTOP*STHE**2)**2
       WT=WT*( 16./3.*(1./(1.-BCTH)**2+1./(1.+BCTH)**2)
     $ -4./3./(1.-BCTH**2) )
       WT = WT*
     &   4.0*PI*ALFAS/2.0*SIGS

      ENDIF
C ============================
CCPY; Jan 16, 2011
C The calculations for 'GL' and 'GG" do not contain exact 
C matrix elements, hence, WT could become negative. Here, we 
C force WT to be always positive. 
      IF(WT.LT.0.D0 .AND.(Boson.Eq.'GL' .OR. Boson.Eq.'GG')) THEN
         WT=0.D0
      ENDIF

C =============================

      RETURN
      END ! WPRES

C ---------------------------------------------------------------------------
       Function CJJ1XFJ_1_AA (sH,tH,uH)
C ---------------------------------------------------------------------------
CsB Returns the finite terms of the O(alpha_s^1) virtual corrections
C   for double photon production.
C   See PRD46('92)2021 (12) (Ohnemus-Owens).
      Implicit Double Precision (a-h,o-z)
      Data PI/3.14159 26535 89793 24D0/

C Definition of v and v1 = v - 1 (s = Q^2)
      v = -uH/sH
      v1 = 1.d0 - v ! = -tH/sH
C Introduce the combinations v/(1-v) and (1-v)/v
      vv1 = v/v1
      v1v = v1/v
C From PRD46('92)2021 (12)
      F_Virt = (2.d0 + vv1)*Log(v)**2 + (2.d0+v1v)*Log(v1)**2 +
     +      (vv1+v1v)*(Log(v)**2+Log(v1)**2-3.d0+2.d0*Pi**2/3.d0) +
     +      2.d0*(Log(v)+Log(v1)+1.d0) + 3.d0*(vv1*Log(v1)+v1v*Log(v))
CsB     +      2.d0*(Log(v)+Log(v1)+1.d0) + 3.d0*(v1v*Log(v1)+v1v*Log(v))

C C function for DY: \delta(1-Z)*CF*(PI^2-8)/4
c      CJJ1XFJ_1 = 0.25*CF*(PI**2-8.0)*FPDF
c
C C function for gamma gamma:
C \delta(1-Z)*CF*(Pi^2 - 8  + 4 - 2 pi^2/3 + 1/(u/t+t/u)(F_Virt-2))/4
C in addition to the DY:
c      CJJ1XFJ_1_AA = 0.25*CF*FPDF*
c     *    ( Pi**2/3.d0 - 4.d0 + 1.d0/(uH/tH+tH/uH)*(F_Virt - 2.d0) )
C Here we return CJJ1XFJ_1_AA without the 0.25*CF*FPDF factor, so in order to
C rescale the DY contribution we have to devide that by (PI**2-8.0) and
C multiply by CJJ1XFJ_1_AA.
      CJJ1XFJ_1_AA = ! 0.25*CF*FPDF*
     *    ( Pi**2/3.d0 - 4.d0 + 1.d0/(uH/tH+tH/uH)*(F_Virt - 2.d0) )

      Return
      End ! CJJ1XFJ_1_AA

C ---------------------------------------------------------------------------
       Function CJJ1XFJ_1_ZZ (s,t,u,x)
C ---------------------------------------------------------------------------
CsB Returns the finite terms of the O(alpha_s^1) virtual corrections
C   for double Z0 production.
C   See Ohnemus-Owens PRD43('91)3626 (12) and (C1).
C   (Compare to photon case: Ohnemus-Owens PRD46('92)2021 (11).
C    Only the coupling constant and definition of F_Virt differ.)

      Implicit Double Precision (a-h,o-z)

      REAL*8 Mz, Mz2,Mz4,Mz6,Mz8,M2B2,M0B2,NN,m0,m2,m41,m42
      COMMON/STAND0/XMW,Mz,XMW2,XMZ2,AEM,GAMW,GAMZ,GAMW2,GAMZ2
      Real*8 Li2
      PARAMETER (Pi = 3.14159 26535 89793 24D0)
      Logical Testing
      Data Testing /.False./

C Testing the Mz -> 0 (di-photon) case.
      If (Testing) then
        Mz = Sqrt(1.d0)
        Print*, ' Mz =', Mz
      End If

C Miscellaneous definitions
      Mz2 = Mz**2
      Mz4 = Mz**4
      Mz6 = Mz**6
      Mz8 = Mz**8
      t2  = t**2
      t3  = t**3
      u2  = u**2
      u3  = u**3
      st  =  s * t
      su  =  s * u
      tu  =  t * u
      t2u = t2 * u
      t3u = t3 * u
      If (Testing) then
        Print*, ' Mz8 =', Mz8
        Print*, ' t3,u3 =', t3,u3
        Print*, ' t2u,t3u =', t2u,t3u
      End If

C Definition of Ai from PRD43('91)3626 (C2)
      A1 = 4.d0 + 4.d0*t/u + 2.d0*u/t + 16.d0*Mz4/tu -
     +     8.d0*Mz2/t - 4.d0*Mz4/t2 - 16.d0*Mz2/u
      A2 = t/(Mz2-t)*
     *     ( -2.d0*t/u + 10.d0*Mz**2/u - 2.d0 + u*Mz2/t2 +
     +        8.d0*Mz6/(t2u) - 2.d0*Mz4/t2 - 2.d0*Mz6/t3 -
     +        u/t - 16.d0*Mz4/tu + 6.d0*Mz2/t )
      A3 = 4.d0 + 4.d0*t/u + 4.d0*u/t + 12.d0*Mz4/tu - 8.d0*Mz2/t -
     +     8.d0*Mz4/t2 - 8.d0*Mz2/u +
     +     t/(Mz2-t)*
     *     ( 3.d0*t/u - 5.d0*Mz2/u + 8.d0 - 4.d0*u*Mz2/t2 -
     +       5.d0*Mz6/t2u - 4.d0*Mz4/t2 + 8.d0*Mz6/t3 + 10.d0*u/t +
     +       7.d0*Mz4/tu - 18.d0*Mz2/t )
      A4 = 6.d0 + 3.d0*t/u + 3.d0*u/t + 5.d0*Mz4/tu -
     +     8.d0*Mz2/t - 3.d0*Mz4/t2 - 10.d0*Mz2/u +
     +     t2/(Mz2-t)**2 *
     *     ( -4.d0*t/u + 24.d0*Mz2/u - 4.d0 - 46.d0*Mz4/tu +
     +       18.d0*Mz2/t + 2.d0*u*Mz2/t2 + 28.d0*Mz6/t2u -18.d0*Mz4/t2 )
      A5 = t2/(Mz2-t)**2 *
     *     ( 4.d0 + 4.d0*u/t + 6.d0*Mz4/tu - 14.d0*Mz2/t +
     +       4.d0*u*Mz4/t3 + 16.d0*Mz8/t3u - 16.d0*Mz6/t3 -
     +       6.d0*u*Mz2/t2 - 20.d0*Mz6/t2u + 22.d0*Mz4/t2 ) +
     +     3.d0*t/u + Mz4/tu - 2.d0*Mz2/t - 6.d0*Mz4/t2 +
     +     t/(s-4.d0*Mz2)*
     *     ( -4.d0*Mz2/u - 1.d0 + u2/t2 - 2.d0*u*Mz2/t2 +
     +       4.d0*u/t + 2.d0*Mz2/t) +
     +     t2/(s-4.d0*Mz2)**2*
     *     ( -t/u - 8.d0*Mz2/u + u/t + 8.d0*Mz2/t )
      A6 = 1.d0/3.d0*
     *     (-6.d0 - 7.d0*t/u - 16.d0*Mz4/tu + 28.d0*Mz2/t + 4.d0*Mz4/t2)
      A7 = t2/su - 4.d0*t*Mz2/su - 3.d0*t/s + 2.d0*Mz4/st +
     +     4.d0*Mz2/s + 4.d0*t/u - 8.d0*Mz2/t +
     +     2.d0*Mz2/(s-4.d0*Mz2)*
     *     ( t2/su - 4.d0*t*Mz2/su - t/s + 4.d0*Mz2/s +
     +       t/u + Mz2/t - 1.d0 ) +
     +     6.d0*Mz4/(s-4.d0*Mz2)**2 * t/s * ( 1.d0 - t/u )
      A8 = t2/(Mz2-t)**2 *
     *     ( -4.d0*t/u + 28.d0*Mz2/u - 4.d0 - 70.d0*Mz4/tu +
     +       22.d0*Mz2/t - 2.d0*u*Mz4/t3 - 28.d0*Mz8/t3u +
     +       18.d0*Mz6/t3 + 2.d0*u*Mz2/t2 + 74.d0*Mz6/t2u -
     +       36.d0*Mz4/t2) +
     +     t/(s-4.d0*Mz2)*(1.d0-t/u) + 6.d0 + 3.d0*Mz4/tu -
     +     6.d0*Mz2/t + 5.d0*Mz4/t2 - 8.d0*Mz2/u
      If (Testing) then
        Print*, ' Ai ='
        Print*, A1, A2, A3
        Print*, A4, A5, A6
        Print*, A7, A8
      End If

C Definition of Bi from PRD43('91)3626 (B3)
      B1 = Li2((Mz2-u)/s) - Li2((s-Mz2)/s) - Li2(-st/(s-Mz2)/(Mz2-t)) +
     +     Li2(-t/(s-Mz2)) + Li2(Mz2/(Mz2-t)) -
     -     .5d0*Log(-st/(s-Mz2)/(Mz2-t))**2 + .5d0*Log(-t/(s-Mz2))**2 +
     +     .25d0*Log(Mz2/(Mz2-t))**2 + .25d0*Log((Mz2-t)/s)**2 -
     +     .25d0*Log(s*Mz2/t2)**2 + Log((Mz2-t)/s)*Log(-(Mz2-u)/t) -
     +     .5d0*Log((Mz2-t)/s)*Log(s*Mz2/t2)
      B7 = 1.d0/x * ( -4.d0*Li2((1.d0-x)/2) + 2.d0*Log((1.d0-x)/2)**2 -
     +                 Log(s/Mz2)**2 + Pi**2/3.d0 )
      If (Testing) then
        Print*, ' Bi ='
        Print*, B1, B7
      End If

C From PRD43('91)3626 (C1)
      F = A1*B1 + A2*( Pi**2 + Log(s/(-t))**2 - Log(Mz**2/s)**2 ) +
     +    A3*Log(Mz**2/(-t)) + A4*Log(s/(-t)) + A5*Log(Mz**2/s) +
     +    A6*Pi**2 + A7*B7 + A8
      If (Testing) then
        Print*, ' F =', F
      End If
      
C Definitions of |M_0^Born|^2/N and |M_2^Born|^2/N, 
C (cf. Eqs.(6) and (7) of PRD43 3626 and my notes on it), 
C where N = N_C e^4 ((g-)^44 + (g+)^4).
C Since we define M0B2 and M2B2 normalized to N, thus in the definition of
C CJJ1XFJ_1_ZZ we can substitute NN = 1. (Putting it differently: NN always
C cancels, in both terms of CJJ1XFJ_1_ZZ.)

      m0  = t/u + u/t
      m2  = Mz2*s/tu
      m41 = Mz4*(1.d0/t2 + 1.d0/u2)
      m42 = 1.d0 - Mz4/tu
      M0B2 = 4.d0*(m0 - m41 + 4.d0*m2)
      M2B2 = 4.d0*(m0 - m41 + 2.d0*m42)
      NN   = 1.d0

C C function for DY: \delta(1-Z)*CF*(PI^2-8)/4
C      CJJ1XFJ_1 = 0.25*CF*(PI**2-8.0)*FPDF
C
C C function for gamma gamma:
C \delta(1-Z)*CF*(Pi^2 - 8  + 4 - 2 pi^2/3 + 1/(u/t+t/u)(F_Virt-2))/4
C in addition to the DY:
c      CJJ1XFJ_1_AA = 0.25*CF*FPDF*
c     *    ( Pi**2/3.d0 - 4.d0 + 1.d0/(uH/tH+tH/uH)*(F_Virt - 2.d0) )
C Here we return CJJ1XFJ_1_AA without the 0.25*CF*FPDF factor, so in order to
C rescale the DY contribution we have to devide that by (PI**2-8.0) and
C multiply by CJJ1XFJ_1_AA.
C      CJJ1XFJ_1_AA = ! 0.25*CF*FPDF*
C     *    ( Pi**2/3.d0 - 4.d0 + 1.d0/(uH/tH+tH/uH)*(F_Virt - 2.d0) )
C
C Comparing PRD43('91)3632 (12) with PRD46('92)2018 (11) we see that
C 1/(2 Pi s) * F^{virt}(v) of the diphoton case corresponds to
C 8 (F(t,u)+F(u,t)) of the diZ case.
C On the other hand Eq.(12) of PRD43('91)3626 has a typo in it.
C I show in my notes that in the Mz -> 0 limit F(t,u)+F(u,t) -> F^{virt}.
C This meas that F^{virt} of the diphoton case corresponds to
C F(t,u)+F(u,t) of the diZ case.

CsB Bug corrected: 09/26/2000
C   For details see the Erratum of PRDD59 114007 and my notes on PRD43 3632.
C      CJJ1XFJ_1_ZZ = ! 0.25*CF*FPDF*
C     *.5d0*( Pi**2/3.d0 - 4.d0 + 1.d0/(u/t+t/u)*(2.d0*F - 2.d0) )
      CJJ1XFJ_1_ZZ = ! 0.25*CF*FPDF*
     *.5d0*(Pi**2/3.d0-M2B2/M0B2+4.d0*NN*2.d0*F/M0B2-3.d0)

C The extra .5 and the 2 are there because we should replace F^{virt} in the
C C funcction by F(t,u)+F(u,t). In order to preserve the coding structure we do
C the following: CJJ1XFJ_1_ZZ(t,u) + CJJ1XFJ_1_ZZ(u,t), which results in
C Pi**2/3.d0 - 4.d0 + 1.d0/(u/t+t/u)*(F(t,u)+F(u,t) - 2.d0)

      Return
      End ! CJJ1XFJ_1_ZZ

CsB -------------------------------------------------------------------------

      Subroutine PartonMomentumLumi(qQBUL,qQBDL,QBqUL,QBqDL,
     >  qGUL,qGDL,GqUL,GqDL,GGL)
cpn2006 Returns parton _momentum_ luminosities for hadrons 1 and 2, 
c       i.e., weighted products of PDFs x1*pdf(ipdf1,x1,Q)*x2*
c       pdf(ipdf2,x2,Q). The types of 1 and 2 are defined in YPARTONS

      Real*8 Frac_n,Frac_p,
     .       qQBUL,qQBDL,QBqUL,QBqDL,qGUL,qGDL,GqUL,GqDL,GGL
      REAL*8 SF(33),SFp(33),SFn(33)
      REAL*8 Q1,Q2
      integer nflmax
      parameter (nflmax=5)
      COMMON/STRUCTURE/q1(-nflmax:nflmax),q2(-nflmax:nflmax)

cpn2006 Parton luminosities for hadrons 1 and 2, normalized by 
c       The type of the hadron 2 for given ibeam is set in YPARTONS
c       Flavor index conventions in Q1, Q2:
c       tbar bbar cbar sbar dbar ubar g  u  d  s  c  b  t
c        -6  -5   -4   -3    -2  -1   0  1  2  3  4  5  6
C uU, cC:
        SF(11)=Q1(1)*Q2(-1)+Q1(4)*Q2(-4)
C dD, sS, bB:
        SF(12)=Q1(2)*Q2(-2)+Q1(3)*Q2(-3)+Q1(5)*Q2(-5)
C Uu, Cc:
        SF(21)=Q1(-1)*Q2(1)+Q1(-4)*Q2(4)
C Dd, Ss, Bb:
        SF(22)=Q1(-2)*Q2(2)+Q1(-3)*Q2(3)+Q1(-5)*Q2(5)

C (u + c) G
        SF(13)=(Q1(1)+Q1(4))*Q2(0)
C (d + s + b) G
        SF(14)=(Q1(2)+Q1(3)+Q1(5))*Q2(0)
C G (u + c)
        SF(23)=Q1(0)*(Q2(1)+Q2(4))
C G (d + s + b)
        SF(24)=Q1(0)*(Q2(2)+Q2(3)+Q2(5))

C G (U + C)
        SF(15)=Q1(0)*(Q2(-1)+Q2(-4))
C G (D + S + B)
        SF(16)=Q1(0)*(Q2(-2)+Q2(-3)+Q2(-5))
C (U + C) G
        SF(25)=(Q1(-1)+Q1(-4))*Q2(0)
C (D + S + B) G
        SF(26)=(Q1(-2)+Q1(-3)+Q1(-5))*Q2(0)

        qQBUL = SF(11)
        qQBDL = SF(12)
        QBqUL = SF(21)
        QBqDL = SF(22)
        qGUL  = SF(13) + SF(25)
        qGDL  = SF(14) + SF(26)
        GqUL = SF(15) + SF(23)
        GqDL = SF(16) + SF(24)

        GGL = Q1(0)*Q2(0)

      Return
      End

     
CsB -------------------------------------------------------------------------

C THIS IS WRESPH.FOR
C "PHASE SPACE" FOR WPRES, WMRES
C THIS IS NOT A REAL PHASE SPACE ROUTINE.
C AS OF 10/1995, THIS WORKS FOR W+, W-, Z0 AND A0.
C
CsB This phase space routine does dQ^2 dy dQT dOmega and includes a 2 Q^2
C   specific to the VBP processes.
C
C WATCH OUT THE HARD-WIRED CUTS BUILD IN THIS SUBROUTINE, THEY
C ARE OPTIMISED FOR W+ PRODUCTION AT 1.8 TEV TEVATRON.

C --------------------------------------------------------------------------
      SUBROUTINE WRESPH(Y,WT)
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
CCAO! The read W-boson mass and width
      Real*8 XMW_IN,XMW_GRID,GAMW_IN
      COMMON/WMASS_IN/XMW_IN,XMW_GRID,GAMW_IN
      Real*8 XMZ_IN,XMZ_GRID,GAMZ_IN
      COMMON/ZMASS_IN/XMZ_IN,XMZ_GRID,GAMZ_IN

      REAL*8 YM12MIN,YM12MAX,QTMin,QTMax, TMmin, Tmmax, ETmin
      COMMON/MYYM12/YM12MIN,YM12MAX,QTMin,QTMax, TMmin, Tmmax, ETmin

      REAL*8 Y,WT
      DIMENSION Y(50)
      INTEGER I
      REAL*8 RN
      EXTERNAL RN

      INTEGER*4 IMONTE
      REAL*8 B_QT2,C_QT2,Q2_V,WT_Q2,QT2_V,WT_QT2,WT_YB,alpha1
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

      REAL*8 QT_V_MIN, QT_V_MAX,YB_V_MAX,Q_V_MIN,Q_V_MAX,QT_V_MXY
      COMMON/HARDCUT/ QT_V_MIN, QT_V_MAX,YB_V_MAX,
     >  Q_V_MIN,Q_V_MAX,QT_V_MXY

CCPY  USED FOR SEPERATING NLO AND NLO_SIG CONTRIBUTIONS
      REAL*8 QT_V_SEP
      COMMON/NLO/ QT_V_SEP

      REAL*8 T_NBAR,T_N
      DIMENSION T_NBAR(4),T_N(4)
      LOGICAL TESTING
      REAL*8 QT2_V_MIN,QT2_V_MAX,TEMP0,TEMP1,TEMP
      REAL*8 Y_EMAX
      Character*40 Boson
      Common / Boson / Boson
!      Integer N_SUD_AB,N_SUD_A,N_SUD_B,N_WIL_C,I_SUD
!      Common / Sudakov / N_SUD_AB,N_SUD_A,N_SUD_B,N_WIL_C,I_SUD
      REAL*8 PCM_L1,ECM_L1,ECM_L2,BETA
      Integer ICalls
      Common / Calls / ICalls
      Integer iYPiece, iReadIn
      Common / YPiece / iYPiece, iReadIn
      Data TESTING /.False./
      INTEGER*4 I_DECAY,I_HADRONIZE
      COMMON/DEC_AY/I_DECAY,I_HADRONIZE

!      Real*8 QT,MT
      
      Character*2 hDecayMode
      Common / hDM / hDecayMode

      Real*8 Reson,xNum
      External Reson
      
      Double Precision ScaPro
      External ScaPro

C      TESTING=.TRUE.
      If (Testing) Print*, ' in WRESPH '

      PI=4.0*ATAN(1.0D0)
      EPS=1.D-4
CsB * At LHC:
      If (W.GT.10.d3) EPS=1.D-3
      WT=-1.0

      Q_V_MIN=YM12MIN
      Q_V_MAX=YM12MAX

CCPY Feb 2010:
      IF(NPART.EQ.11) THEN ! FOR RESUMMED
        if (Boson.eq.'AA'.or.Boson.eq.'AG') then
          QT_V_MIN=QT_V_SEP
        else
          QT_V_MIN=0.D0          
        endif
      ELSEIF(NPART.EQ.12 .or.(NPart.ge.21 .and. NPart.le.24)) THEN ! FOR NLO
          QT_V_MIN=QT_V_SEP
      ENDIF

CsB For processes requiring grids:
C   QT_V_MAX and YB_V_MAX are determined in ReadIn and passed through COMMON
C   QT_V_MAX is the maximal qT for the resummed grid.

C FIND Q2_V
      If (Boson.Eq.'H0' .And. hDecayMode.Eq.'ZZ' .And. 
     .    I_DECAY.EQ.1) then ! H -> Z Z -> leptons
        xNum = 9.d0
        YM1 = Reson(XMZ,GamZ,xNum)
        YM2 = Reson(XMZ,GamZ,xNum)
      Else
        YM1=XMASS(1)
        YM2=XMASS(2)
      EndIf
      
CsB * Set VB (3&4) mass and width
      YMW=0.d0
      YGW=0.d0
      If (Boson.Eq.'W+' .or. Boson.Eq.'W-' .or.
     >    Boson.Eq.'W+_RA' .or. Boson.Eq.'W-_RA' .or.
     >    Boson.Eq.'W+A(SV)' .or. Boson.Eq.'W-A(SV)' .or.
     >    Boson.Eq.'HP' .or. Boson.Eq.'HM') then
        YMW=XMW
        YGW = GAMW_IN
        If (GAMW_IN.LT.0.) YGW = GAMW_NLO
      ELSE IF (Boson.Eq.'Z0' .or. Boson.Eq.'Z0_RA' .or.
     .         Boson.Eq.'DY' .or. Boson.Eq.'HZ') then
        YMW=XMZ
        YGW = GAMZ_IN
CCPY Aug 2016:        If (GAMW_IN.LT.0.) YGW = GAMZ_NLO
        If (GAMZ_IN.LT.0.) YGW = GAMZ_NLO
      ELSE IF (Boson.Eq.'H0') then
        YMW=XMHIGGS
        YGW=GAMH
      ELSE IF (Boson.Eq.'A0' .or. Boson.Eq.'GL'
     .    .or. Boson.Eq.'AA' .or. Boson.Eq.'AG'
     .    .or. Boson.Eq.'ZZ' .or. Boson.Eq.'ZG') then
        YMW=0.d0
        YGW=0.d0
      ENDIF
!     Print*, 'YMW,YGW =', YMW,YGW 

CsB___Choose a scheme to generate PS
CCPY Oct 2008, INCLUDE 'DRELL-YAN PAIR 'DY' PROCESS
      If (Boson.Eq.'W+' .or. Boson.Eq.'W-' .or.
     .    Boson.Eq.'W+_RA' .or. Boson.Eq.'W-_RA' .or.
     .    Boson.Eq.'W+A(SV)' .or. Boson.Eq.'W-A(SV)' .or.
     .    Boson.Eq.'Z0' .or. Boson.Eq.'Z0_RA' .or. 
     .    Boson.Eq.'H0' .or. Boson.Eq.'DY') then
        IMONTE=2
      else if (Boson.Eq.'AA'.or.Boson.Eq.'AG') then
        imonte=3
      else if (Boson.Eq.'A0') then
        imonte=1
      ELSE
c        imonte = 1
        IMONTE=0
      END IF

CCPY__IN THIS CODE, X1 AND X2 SHOULD BE INTERPRETED AS X_A AND X_B

      IF(IMONTE.EQ.0.or.IMONTE.EQ.3) THEN ! Use linear shape ---------------
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
CsB The high qT end of the Y piece is read from the Y grid.
cC      If (QT_V_MXY.LE.1.d-6) QT_V_MXY = 100.d0
      IF(NPART.EQ.11) THEN ! RESUMMED CASE
C ASSUME THE SHAPE AS d(sigma)/d(QT2_V) = c* EXP(-b*QT2_V)
C WHERE b IS FITTED WITH RESUMMED RESULT OF UNPOLARIZED W-BOSON CALCULATIONS
C QT2_V_MAX IS THE MAXIMUM QT2 OF V
CCPY INPUT FROM FIT TO W^+ AT 1.8 TEV, P-PBAR
        If (iYPiece.Eq.0) then
          B_QT2=(5.9-1.2)/(QT_V_MAX**2-0.1)
          QT2_V_MAX=QT_V_MAX**2
        Else If (iYPiece.Eq.1) then
          B_QT2=(5.9-1.2)/(QT_V_MXY**2-0.1)
C_yfu          QT2_V_MAX=QT_V_MXY**2
          QT2_V_MAX=MIN(QT_V_MXY**2,QTMAX**2)
        End If

C HARD-WIRED CUT
c        print *, qt_v_min, qt_v_sep
        QT2_V_MIN=QT_V_MIN**2
        TEMP0=EXP(-B_QT2*QT2_V_MIN)
        TEMP=TEMP0-EXP(-B_QT2*QT2_V_MAX)
        C_QT2=B_QT2/TEMP
        QT2_V=-DLOG(TEMP0-TEMP*Y(2))/B_QT2
        QT_V=DSQRT(QT2_V)
        WT_QT2=1.0/C_QT2/(TEMP0-TEMP*Y(2))

cpn2005 Alternative qT sampling function
cpn2005 the qT^2 range is stretched as z~ (qT^2)^alpha1,
cpn2005 with alpha1 << 1, and qT^2_min <= qT^2 <= qT^2_max corresponding
cpn2005 to 0 <= z <= 1
c        QT2_V_MIN=QT_V_MIN**2
c        alpha1=0.2
c        temp0 = (qt2_v_min/qt2_v_max)**alpha1
c        temp1 = temp0 + (1d0-temp0)*y(2)
c        temp =  temp1**(1d0/alpha1)
c        qt2_v = qt2_v_max*temp
c        qt_v =  dsqrt(qt2_v)
c        wt_qt2=qt2_v_max*(1d0-temp0)*temp/alpha1/temp1

      Elseif( NPART.EQ.12 .or. (NPart.ge.21 .and. nPart.le.24)) THEN ! NLO CASE
        If (iYPiece.Eq.0) then
          QT2_V_MAX=QT_V_MAX**2
        Else If (iYPiece.Eq.1) then
C_yfu          QT2_V_MAX=QT_V_MXY**2
          QT2_V_MAX=MIN(QT_V_MXY**2,QT_V_MAX**2)
        End If
        QT2_V_MIN=QT_V_MIN**2

        If (iMonte.Eq.0) then   ! assumes d(sigma)/d(QT^2) = 1/X * LOG(X) 
                                !WITH X = QT^2/Q^2            
          TEMP0=LOG(QT2_V_MIN/Q2_V)**2
          TEMP=LOG(QT2_V_MAX/Q2_V)**2-TEMP0
          C_QT2=2.0/TEMP
          TEMP1=SQRT(TEMP0+TEMP*Y(2))
          QT2_V=Q2_V*EXP(-TEMP1)
          QT_V=DSQRT(QT2_V)
          WT_QT2=-1.0/C_QT2/TEMP1*QT2_V
        elseif (iMonte.eq.2 .or. iMonte.eq.1) then ! assumes d(sigma)/d(QT^2) ~ (1/QT2_V) ! yfu add iMonte .eq. 1
          QT2_V=QT2_V_MIN*(QT2_V_MAX/QT2_V_MIN)**Y(2)
          QT_V=DSQRT(QT2_V)
          WT_QT2=QT2_V*DLOG(QT2_V_MAX/QT2_V_MIN)
        Else if (iMonte.Eq.3) then !optimal qT sampling for AA production
cpn2005 the qT^2 range is stretched as z~ (qT^2)^alpha1,
cpn2005 with alpha1 << 1, and qT^2_min <= qT^2 <= qT^2_max corresponding
cpn2005 to 0 <= z <= 1
          alpha1=0.2
          temp0 = (qt2_v_min/qt2_v_max)**alpha1
          temp1 = temp0 + (1d0-temp0)*y(2)
          temp =  temp1**(1d0/alpha1)
          qt2_v = qt2_v_max*temp
          qt_v =  dsqrt(qt2_v)
          wt_qt2=qt2_v_max*(1d0-temp0)*temp/alpha1/temp1
        End if                  !iMonte
      ENDIF

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

C DO THE LEPTON KINEMATICS IN THE REST FRAME OF W
C       CHANGE THE VARIABLES ( COS(THE_STAR), PHI_STAR ) IN C.M.
C       FRAME TO ( R(1), R(2) ).
      R(1)=Y(4)
ctemp      R(1)=RN(2)

CsB___
      If (iYPiece.Eq.1) then
        R(2)=Y(5)  ! The Y piece depends on phi...
        If (Boson.Eq.'AA' .or. Boson.Eq.'AG'
     . .or. Boson.Eq.'ZZ' .or. Boson.Eq.'ZG'
     . .or. Boson.Eq.'H0')
     .  R(2)=RN(1) ! ... except for AA (AG), ZZ (ZG) or H0.
      Else
        R(2)=RN(1)
!        If (Boson.Eq.'AG') R(2)=Y(5)
      End If

      THE_STA=DACOS(2.D0*R(1)-1.D0)
      PHI_STA=2.D0*PI*R(2)
c      PHI_STA=0.D0

C FOR THE_STA
      WT_THE=2.D0
C FOR PHI_STA
      WT_PHI=2.D0*PI

C ADD THE FACTOR FROM d^4(Q): PI/2
      WT=PI/2.0D0*WT_QT2*WT_Q2*WT_YB*WT_THE*WT_PHI
ctemp      WT=PI/2.0D0*WT_QT2*WT_Q2*WT_YB

CsB___Check that the decay is kinematically allowed: Q > m1 + m2
      If (Q_V.LT.YM1+YM2) then
!        Print*, ' Q < m1 + m2 '
        WT = 0.
        Return
      EndIf

C FOR MASSIVE FINAL STATE PARTICLES (ARISING FROM 2-BODY PHASE SPACE FACTOR).
      If (TESTING) Print*,'THE_STA,PHI_STA',THE_STA,PHI_STA
      CTHE=COS(THE_STA)
      STHE=SIN(THE_STA)
      CPHI=COS(PHI_STA)
      SPHI=SIN(PHI_STA)
      IF(TESTING) Print*,'Q_V,YM1,YM2=',Q_V,YM1,YM2
      PCM_L1=DSQRT((Q_V**2-(YM1+YM2)**2)*(Q_V**2-(YM1-YM2)**2))/
     >2.0/Q_V
      ECM_L1=(Q_V**2+YM1**2-YM2**2)/2.0/Q_V
      ECM_L2=(Q_V**2-YM1**2+YM2**2)/2.0/Q_V
      IF(TESTING) THEN
        PRINT*,'YM1,YM2,PCM_L1,ECM_L1,ECM_L2'
        PRINT*,YM1,YM2,PCM_L1,ECM_L1,ECM_L2
      ENDIF
C      BETA=PCM_L1/ECM_L2
      BETA=2.0*PCM_L1/Q_V

C ADD THE OTHER FACTORS
      WT=WT*
     &   1.d0/(2.d0*Pi)**3.d0/8.d0 * BETA * ! <= 1 -> 2 phase space dOmega
     &   1.d0/(Pi*qT_V) *                   ! <= d4Q = pi/2 dqT^2 dQ dy
     &   1.d0/(2*Q_V**2)

ctemp      WT=WT*
c     &   1.d0/(Pi*qT_V) *                   ! <= d4Q = pi/2 dqT^2 dQ dy
c     &   1.d0/(2*Q_V**2)

      IF(TESTING) THEN
        PRINT*,'X1,X2,TAU,YBOOST'
        PRINT*,X1,X2,TAU,YBOOST
        PRINT*,'Q_V,QT_V,YB_V,THE_STA,PHI_STA'
        PRINT*,Q_V,QT_V,YB_V,THE_STA,PHI_STA
        PRINT*,'WT_QT2,WT_Q2,WT_YB,WT_THE,WT_PHI,WT'
        PRINT*,WT_QT2,WT_Q2,WT_YB,WT_THE,WT_PHI,WT
      ENDIF

C**************************************
C In the Collins-Soper frame

      R2INV=1.0D0/SQRT(2.0D0)

      P_NBAR(1)=0.D0
      P_NBAR(2)=0.D0
      P_NBAR(3)=R2INV
      P_NBAR(4)=R2INV

      P_N(1)=0.D0
      P_N(2)=0.D0
      P_N(3)=-R2INV
      P_N(4)=R2INV

C GIVE A RANDOM CHOICE OF AZIMUTHAL ANGLE FOR W-BOSON, IN LAB FRAME
      PHI_Q=2.0*PI*RN(2)
      TM2_Q=Q2_V+QT2_V
      TM_Q=SQRT(TM2_Q)

CsB   Massive lepton(1) 4 momentum in CS frame
CsB   To test differences between W rest frame being the CS frame and the
C     one with the z direction in the W motion direction.
      P_L1CS(1) = PCM_L1 * SThe * CPhi
      P_L1CS(2) = PCM_L1 * SThe * SPhi
      P_L1CS(3) = PCM_L1 * CThe
      P_L1CS(4) = ECM_L1
      P_L1CS(5) = YM1

c      P_L1CS(5) = Sqrt(-P_L1CS(1)**2 -P_L1CS(2)**2 -P_L1CS(3)**2 +
c     >                  P_L1CS(4)**2)

CsB 1&2 (vector boson) 4 momentum in lab frame
      P_Q(1)=QT_V*COS(PHI_Q)
      P_Q(2)=QT_V*SIN(PHI_Q)
      P_Q(3)=TM_Q*(EXP(YB_V)-EXP(-YB_V))/2.0D0
      P_Q(4)=TM_Q*(EXP(YB_V)+EXP(-YB_V))/2.0D0

      P_Q(5)=SQRT(-P_Q(1)**2-P_Q(2)**2-P_Q(3)**2+
     >        P_Q(4)**2)

      QPLUS=R2INV*(P_Q(4)+P_Q(3))
      QMINUS=R2INV*(P_Q(4)-P_Q(3))

cpn2006 X,Y,Z unity vectors in the CS frame, given in the lab frame
      P_X(1)=TM_Q*COS(PHI_Q)/Q_V
      P_X(2)=TM_Q*SIN(PHI_Q)/Q_V
      P_X(3)=P_Q(3)*QT_V/TM_Q/Q_V
      P_X(4)=P_Q(4)*QT_V/TM_Q/Q_V

      P_Z(1)=0d0; P_Z(2)=0d0; P_Z(3)=P_Q(4)/TM_Q; P_Z(4)=P_Q(3)/TM_Q;

cpn This is original calculation; the equation for P_X loses 
cpn precision at very small Q_T
c      DO 70 I=1,4
c         P_X(I)=-Q_V/QT_V/TM_Q*(QPLUS*P_NBAR(I)+QMINUS*P_N(I)-
c     >   (TM2_Q/Q2_V)*P_Q(I))
c         P_Z(I)=1.0/TM_Q*(QPLUS*P_NBAR(I)-QMINUS*P_N(I))
c70    CONTINUE

      CALL FINDY(Q_V,P_Q,P_X,P_Z,P_Y)

      P_X(5)=-P_X(1)**2-P_X(2)**2-P_X(3)**2+
     >        P_X(4)**2
      P_Y(5)=-P_Y(1)**2-P_Y(2)**2-P_Y(3)**2+
     >        P_Y(4)**2
      P_Z(5)=-P_Z(1)**2-P_Z(2)**2-P_Z(3)**2+
     >        P_Z(4)**2

      IF(TESTING) THEN
        PRINT*,' P_Q, Q_V = ', P_Q,Q_V
        PRINT*,' P_X =',P_X
        PRINT*,' P_Y =',P_Y
        PRINT*,' P_Z =',P_Z
      ENDIF

      IF(ABS(P_Q(5)-Q_V).GT.EPS) THEN
       PRINT*,' ERROR IN P_Q'
C       CALL EXITGracefully
       RETURN
      ENDIF

      IF(ABS(P_X(5)+1.0).GT.EPS .OR. ABS(P_Y(5)+1.0).GT.EPS
     >.OR. ABS(P_Z(5)+1.0).GT.EPS) THEN
       PRINT*,' POTENTIAL ERROR IN MASS: P_X, P_Y  OR  P_Z'
C       Call ExitGracefully
       RETURN
      ELSE
       P_X(5)=SQRT(-P_X(5))
       P_Y(5)=SQRT(-P_Y(5))
       P_Z(5)=SQRT(-P_Z(5))
      ENDIF

      DO 80 I=1,4
CsB   Boost from CS frame to lab frame.
CCPY  Only works for massless particles (i.e. YM1=YM2=0 ).
C         P_L1(I)=Q_V/2.0D0*(P_Q(I)/Q_V+STHE*CPHI*P_X(I)+
C     >   STHE*SPHI*P_Y(I)+CTHE*P_Z(I))
C         P_L2(I)=Q_V/2.0D0*(P_Q(I)/Q_V-STHE*CPHI*P_X(I)-
C     >   STHE*SPHI*P_Y(I)-CTHE*P_Z(I))
C
         P_L1(I)=ECM_L1*P_Q(I)/Q_V+PCM_L1*(STHE*CPHI*P_X(I)+
     >   STHE*SPHI*P_Y(I)+CTHE*P_Z(I))
         P_L2(I)=ECM_L2*P_Q(I)/Q_V+PCM_L1*(-STHE*CPHI*P_X(I)-
     >   STHE*SPHI*P_Y(I)-CTHE*P_Z(I))
80    CONTINUE
CsB   Note: the above is not the same as the one using BOOST_5
!      CALL BOOST_5 (P_L1CS, P_Q, P_L1)
CsB   Check CS -> Lab boost
!      QT = Sqrt(P_Q(1)**2 + P_Q(2)**2)
!      MT = Sqrt(P_Q(5)**2 + QT**2)
!      Print *,' Q,QT,MT = ',P_Q(5),QT,MT
!      Print *,' P_Q = ', P_Q
!      Print *,' P_L1(4)',P_L1(4),1.d0/P_Q(5)/MT*(P_Q(4)*MT*P_L1CS(4) +
!     +  P_Q(4)*QT*P_L1CS(1) + P_Q(3)*P_Q(5)*P_L1CS(3))
!      Print *,' P_L1(1)',P_L1(1),1.d0/P_Q(5)/MT*(QT*MT*P_L1CS(4) + 
!     +  MT**2*P_L1CS(1) )
!      Print *,' P_L1(2)',P_L1(2),P_L1CS(2)
!      Print *,' P_L1(3)',P_L1(3),1.d0/P_Q(5)/MT*(P_Q(3)*MT*P_L1CS(4) +
!     +  P_Q(3)*QT*P_L1CS(1) + P_Q(4)*P_Q(5)*P_L1CS(3))
!      Stop

C BOOSTED TO LAB FRAME
      P_L1(5)=-P_L1(1)**2-P_L1(2)**2-P_L1(3)**2+
     >        P_L1(4)**2
      P_L2(5)=-P_L2(1)**2-P_L2(2)**2-P_L2(3)**2+
     >        P_L2(4)**2

      IF(TESTING) PRINT*,' P_L1(5),P_L2(5) = ',P_L1(5),P_L2(5)

      P_L1(5)=SQRT(ABS(P_L1(5)))
      P_L2(5)=SQRT(ABS(P_L2(5)))

CsB   To test differences between W rest frame being the CS frame and the
C     one with the z direction in the W motion direction.
c      If (.False.) then
c        If (ICalls.Eq.0) then
c          Print*, ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! '
c          Print*, ' In this run the W rest frame is NOT the CS frame '
c          Print*, ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! '
c          ICalls = 1
c        End If
c        Call Boost_5 (P_L1CS, P_Q, P_L1)
c         DO I=1,4
c          P_L2(I)=P_Q(I)-P_L1(I)
c         ENDDO
c         P_L2(5)=YM2
c      End If

CsB   Check lepton masses in Lab frame.
      IF(YM1.GT.EPS .OR. YM2.GT.EPS) THEN
        IF(ABS(P_L1(5)-YM1).GT.P_L1(5)*EPS .OR.
     >  ABS(P_L2(5)-YM2).GT.P_L2(5)*EPS) THEN
          PRINT*,' POTENTIAL ERROR IN MASS: P_L1 OR P_L2 (3)'
CCPY       Call ExitGracefully
          PRINT *,'P_L1 =  ',P_L1
          PRINT *,'P_L2 =  ',P_L2
          PRINT *,'P_Q  =  ',P_Q
          RETURN
        ENDIF
      ELSE
        IF(P_L1(5).GT.4.d0*EPS .OR. P_L2(5).GT.4.d0*EPS) THEN
          PRINT*,' POTENTIAL ERROR IN MASS: P_L1 OR P_L2 (4)'
CCPY       Call ExitGracefully
          PRINT *,'P_L1 =  ',P_L1
          PRINT *,'P_L2 =  ',P_L2
          PRINT *,'P_Q  =  ',P_Q
          RETURN
        ENDIF
      ENDIF

CCPY  TESTING PROGRAM
      IF(TESTING) THEN
        DO 73 I=1,4
           T_NBAR(I)=TM_Q/QPLUS/2.0*(P_Z(I)-QT_V/Q_V*P_X(I)
     >               +TM_Q/Q_V/Q_V*P_Q(I))
           T_N(I)=TM_Q/QMINUS/2.0*(-P_Z(I)-QT_V/Q_V*P_X(I)
     >               +TM_Q/Q_V/Q_V*P_Q(I))
73      CONTINUE
        PRINT*,' T_NBAR =',T_NBAR
        PRINT*,' T_N =',T_N
      ENDIF

CsB Fill the lab momenta PLAB for final particles 1,2 and 9=1&2
      DO 90 I=1,4
C NO INFORMATION ON TWO INCOMING PARTONS: SET THEM TO BE ZERO
        PLAB(I,1)=0.D0
        PLAB(I,2)=0.D0
C FERMION
        PLAB(I,3)=P_L1(I)
C ANTI-FERMION
        PLAB(I,4)=P_L2(I)
C VECTOR BOSON
        PLAB(I,9)=P_Q(I)
90    CONTINUE

      IF(TESTING) THEN
C       PRINT *,'((PLAB(J1,J2),J1=1,4),J2=3,4)'
C       PRINT *,((PLAB(J1,J2),J1=1,4),J2=3,4)
        PRINT *,'P_L1 =  ',P_L1
        PRINT *,'P_L2 =  ',P_L2
        PRINT *,'P_Q  =  ',P_Q
      ENDIF
CsB___Check |PLAB(_,3)|,|PLAB(_,4)| 
!        Print *,' |PLAB(_,3)|,|PLAB(_,4)| = ',
!     >          Sqrt(ScaPro(3,3)),Sqrt(ScaPro(4,4))

CsB___For the case when each of the two final state particles furter decay 
C     into two massles particles, like for H -> Z0 Z0 -> e+ e- mu+ mu-
      IF(Boson.Eq.'H0' .and. hDecayMode.Eq.'ZZ' .and. i_Decay.Eq.1) THEN 
        CALL HZZDECAY (P_Q,P_L1CS,P_L1,P_L2,
     ,P_X,P_Y,P_Z,STHE,CTHE,SPHI,CPHI,Q_V)
      ENDIF

      RETURN
      END ! WRESPH

CsB * *******************************************************************

      SUBROUTINE YMAXIMUM(ECM,Q_V,QT_V,YMAX)
      IMPLICIT NONE

      REAL*8 ECM,Q_V,QT_V,YMAX
      REAL*8 X,ACOSH,T1,T2,T3

      ACOSH(X)=DLOG(X+DSQRT(X**2-1.0))

      T1=ECM**2+Q_V**2
C      Print*, ' QT_V,Q_V ', QT_V,Q_V
      T2=DSQRT(QT_V**2+Q_V**2)
      T3=2.0*ECM*T2
C      Print*, ' ECM, X**2-1.0 ', ECM, (T1/T3)**2-1.0
      YMAX=ACOSH(T1/T3)

      RETURN
      END ! YMAXIMUM

CsB * *******************************************************************

      SUBROUTINE FINDY(Q_V,Q,X,Z,Y)
      IMPLICIT NONE
      REAL*8 Q,X,Z,Y,Q_V
      DIMENSION Q(5),X(5),Z(5),Y(5)
      INTEGER*4 I

      Y(1)=Q(4)*X(3)*Z(2) - Q(3)*X(4)*Z(2) - Q(4)*X(2)*Z(3) +
     >     Q(2)*X(4)*Z(3) + Q(3)*X(2)*Z(4) - Q(2)*X(3)*Z(4)

      Y(2)=-Q(4)*X(3)*Z(1) + Q(3)*X(4)*Z(1) + Q(4)*X(1)*Z(3) -
     >      Q(1)*X(4)*Z(3) - Q(3)*X(1)*Z(4) + Q(1)*X(3)*Z(4)

      Y(3)=Q(4)*X(2)*Z(1) - Q(2)*X(4)*Z(1) - Q(4)*X(1)*Z(2) +
     >     Q(1)*X(4)*Z(2) + Q(2)*X(1)*Z(4) - Q(1)*X(2)*Z(4)

      Y(4)=-Q(3)*X(2)*Z(1) + Q(2)*X(3)*Z(1) + Q(3)*X(1)*Z(2) -
     >      Q(1)*X(3)*Z(2) - Q(2)*X(1)*Z(3) + Q(1)*X(2)*Z(3)

C CORRECTED TO BJORKEN-DRELL MATRIC
      Y(4)=-Y(4)
      DO I=1,4
        Y(I)=Y(I)/Q_V
      ENDDO

      RETURN
      END

CsB * *******************************************************************

      Subroutine ThreeDInterp ( iDD1,iDD2,iDD3, YDAS,YDAA, qT,yIn,Qv,
     >                          sigS,sigA, ErrSgS,ErrSgA )
C
C     The following convention apply:
C     X1A (X2A, X3A) contains the qT (y, Q) vector
C     iD1+1 (iD2+1, iD3+1) is the number of qT (y, Q) values in the grid
C     iDD1 (iDD2, iDD3) is the number of qT (y, Q)
C          values used at the interpolation
C     The 1st indices of YAS and YAA belong to qT, the 2nd to y, the 3rd to Q.
      Implicit Double Precision (A-H,O-Z)
C     The dimensions below give the maximal possible gridsize.
      Parameter (MD1 = 200, MD2 = 300, MD3 = 100)
      Dimension X1A(0:MD1), X2A(0:MD2), X3A(0:MD3)
      Dimension YAS(0:MD1,0:MD2,0:MD3), YAA(0:MD1,0:MD2,0:MD3)
      Dimension YDAS(iDD3,iDD1,iDD2), YDAA(iDD3,iDD1,iDD2)
      Dimension X1Y(0:MD1), X2Y(0:MD2), X3Y(0:MD3)
      Dimension YAY(0:MD1,0:MD2,0:MD3,0:7)
      Logical Testing
      Common /ECMC/ ECM
      Common / Arrays / X1A,X2A,X3A, YAS,YAA,  X1Y,X2Y,X3Y, YAY
      Common / Dims / iD1,iD2,iD3, iY1,iY2,iY3
      Character*40 Boson
      Common / Boson / Boson
!      Integer N_SUD_AB,N_SUD_A,N_SUD_B,N_WIL_C,I_SUD
!      Common / Sudakov / N_SUD_AB,N_SUD_A,N_SUD_B,N_WIL_C,I_SUD
      Common / nyEdge / nyEdge
      COMMON/HARDCUT/ QT_V_MIN, QT_V_MAX,YB_V_MAX,
     >  Q_V_MIN,Q_V_MAX,QT_V_MXY
      Integer LTOpt
      Common / Process / LTOpt
      Common / Matching / QT_V_Mat
      COMMON/AA/W,QCDL,X1X2MIN,NJET,N_PRIMARY_PART,NF,I_STR,IBEAM,NPART,
     &  NEVENT,NORDER,N_SECOND_PART

CCPY
      LOGICAL Y_ONE_SIDE
      COMMON/GRID_FILE_FORM/ Y_ONE_SIDE
      INTEGER iDD2_FIT

      Common / TI / nTI
      Parameter (yMx = 5.d0)
      Data eps / 1.d-9 /
      Data Testing / .False. /
      
cdump
C      Testing = .True.

CsB Set flag that tells Polin3 that it was called from ThreeDInterp
      nTI = 0

CsB * Testing
      If (Testing) then
        Print*, ' iD1, iD2, iD3 = ', iD1, iD2, iD3
        Print*, ' iDD1/2, iDD2/2, iDD3/2 = ', iDD1/2, iDD2/2, iDD3/2
      EndIf

cCsB The examination of the matching condition is moved to before the
cC   subroutine call of ThreeDInterp
cCCsB If e.g. qT_W > QT_V_Mat then the CSS piece is set to be zero
cCC   (simple matching).
cC      If (qT.GT.QT_V_Mat) then
cC        If (Testing) Print*, 'QT_V_Mat = ', QT_V_Mat
cC        SigS   = 0.d0
cC        ErrSgS = 0.d0
cC        SigA   = 0.d0
cC        ErrSgA = 0.d0
cC        Return
cC      End If

      y = yIn
CsB * Take care of the photon case when the rapidity spectrum is symmetric.
C     (p N collision is an exception)
CCPY Sept 2009: Add iBeam.Eq.1 option
      If (Boson.Eq.'AA' .and. iBeam.Eq.0) Goto 160

      IF(Y_ONE_SIDE) THEN
        y = Abs(yIn)
      EndIf
  160 Continue

CsB * Take special care of the p p case when the L0 rapidity spectrum is
C     symmetric and A3 is anti-symmetric wrt y -> -y.
C     At the LHC
CsB   Commented out: 2005 02 15.  Not needed when grid spans both signs of y.
C     If (iBeam.Eq.1 .and. Boson.Eq.'W+') y = Abs(yIn)

CsB * If Q, qT or y right on the edge of the grid
c      If (qT.Eq.X1A(0)) qT = qT + eps
c      If ( y.Eq.X2A(0))  y =  y + eps
c      If (Qv.Eq.X3A(0)) Qv = Qv + eps
c      If (qT.Eq.X1A(iD1)) qT = qT - eps
c      If ( y.Eq.X2A(iD2))  y =  y - eps
c      If (Qv.Eq.X3A(iD3)) Qv = Qv - eps

CsB * Find out the location of Q, qT, y
      call LOCATE(X1A, iD1, qT, LqT)
      call LOCATE(X2A, iD2,  y, Ly)
      call LOCATE(X3A, iD3, Qv, LQ)
      If (Testing) Print *, ' LqT, Ly, LQ = ', LqT, Ly, LQ

CsB * If not on the grid
      If ( Abs((X1A(iD1)-qT)*(qT-X1A(0))).LT.0.e-9 .or.
     .     Abs((X2A(iD2)- y)*( y-X2A(0))).LT.0.e-9 .or.
     .     Abs((X3A(iD3)-Qv)*(Qv-X3A(0))).LT.0.e-9 ) then
        Print*, ' Warning from 3DInterp: Extrapolation used'
        Print*, ' Q, qT, y = ', Qv,qT,y
        Print*, ' LQ,LqT,Ly= ', LQ,LqT,Ly
        Print*, ' iD3,iD1,iD2= ', iD3,iD1,iD2
      EndIf

CsB * If close to the ends
CsB   ( This works for even iDD for sure )
      If (LqT.LT.iDD1/2)     LqT = iDD1/2
      If (LqT.GT.iD1-iDD1/2) LqT = iD1-iDD1/2+1
      If ( Ly.LT.iDD2/2)      Ly = iDD2/2
      If ( Ly.GT.iD2-iDD2/2)  Ly = iD2-iDD2/2+1
      If ( LQ.LT.iDD3/2)      LQ = iDD3/2
      If ( LQ.GT.iD3-iDD3/2)  LQ = iD3-iDD3/2+1
      If (Testing) Print *, ' LqT, Ly, LQ = ', LqT, Ly, LQ
CsB   For the asymptotic piece we can not use the lowest qT point
C     (which is a dummy sigma(qT=0)=0).
      If (nPart.Eq.12 .and. LqT.Eq.iDD1/2) LqT = iDD1/2 + 1

CsB * Take care of the end points of the rapidity - if necessary
C     ( The right way to do it would be:
C       If the interpolation interval includes yMax
C       than (within Polin3) move the closest point of X2A
C       where the xsect = 0 into yMax.
C       The problem with the above approach is that if the interpolation
C       point will fall very close to the edge, then Polint gives back
C       a weird negative value. So may be, we should set a limit how close
C       to the edge we interpolate using Polint. )
      If (nyEdge.EQ.1) then
        If (Ly.LT.5 .or. Ly.GT.iD2-5) then
          Call YMAXIMUM (ECM,Qv,qT,yMax)
          If (Testing) then
            Call Locate (X2A, iD2, yMax, LyMax)
            Print*, ' yMax, LyMax = ', yMax, LyMax
          EndIf
          Print*, ' Fix yMx. Should be sg like X2A(iD2) ... '
          Call ExitGracefully
          If (Abs(yMax).LT.yMx) then
            X2A(0)  = -yMx
            X2A(1)  = -yMax
            X2A(13) =  yMax
            X2A(14) =  yMx
          Else
            X2A(0)  = -yMax
            X2A(1)  = -yMx
            X2A(13) =  yMx
            X2A(14) =  yMax
          End If
        End If
      End If

CsB * Fill the dummy array(s) you want to interpolate
C     ( This works for even iDD for sure )
C     Since the LO and DSi pieces have only onr qT (no qT integral in the
C     phase space the order of interpolation [1->y 2->qT 3->Q] is preserved
C     so that the interpolation does not have to be changed.)

cCsB   The YD arrays filled up this way so that y gets interpolated
cCsB   first, then qT, and finally Q
c        Do i = 1,iDD1
c          Do j = 1,iDD2
c            Do k = 1,iDD3
c              YDAS(k,i,j) = YAS(LqT-iDD1/2+i-1,Ly-iDD2/2+j-1,
c     >          LQ-iDD3/2+k-1)
c              YDAA(k,i,j) = YAA(LqT-iDD1/2+i-1,Ly-iDD2/2+j-1,
c     >          LQ-iDD3/2+k-1)
c              If (Testing) then
ccc                Print*, ' i,j,k, YDAS(i,j,k),YDAA(i,j,k) = ',
ccc     >                    i,j,k, YDAS(i,j,k), YDAA(i,j,k)
cc                    Print*, X3A(LQ-iDD3/2+k-1), X1A(LqT-iDD1/2+i-1),
cc     >                  X2A(Ly-iDD2/2+j-1), YDAS(k,i,j),YDAA(k,i,j)
c              End If
c            End Do
c          End Do
c        End Do

CCPY Nov 2009:
      IF(Y_ONE_SIDE) THEN
        iDD2_FIT=iDD2-1
      ELSE
        iDD2_FIT=iDD2
      EndIf


        Do i = 1,iDD1
CCPY          Do j = 1,iDD2
          Do j = 1,iDD2_FIT
            Do k = 1,iDD3
CsB   Note: the YD arrays filled up this way so that y gets interpolated
CsB   first, then qT, and finally Q.
CCAO          If (nPart.Eq.13 .or. nPart.Eq.14) then ! DSi or LO
              If (nPart.Eq.13 .or. nPart.Eq.14 .or. 
     .            nPart.eq.16) then ! DSi or LO
                YDAS(k,i,j) = YAS(LqT-iDD1/2+i-1,Ly-iDD2/2+j-1,
     >            LQ-iDD3/2+k-1)
                YDAA(k,i,j) = YAA(LqT-iDD1/2+i-1,Ly-iDD2/2+j-1,
     >            LQ-iDD3/2+k-1)
CCAO          Else If (nPart.Eq.11 .or. nPart.Eq.12) then ! CSS or Asy
              Else If (nPart.Eq.11 .or. nPart.Eq.12 .or. 
     .                 nPart.eq.15) then ! CSS or Asy
CsB   Note: the YD arrays filled up this way so that y gets interpolated
CsB   first, then Q, and finally qT.
                YDAS(i,k,j) = YAS(LqT-iDD1/2+i-1,Ly-iDD2/2+j-1,
     >            LQ-iDD3/2+k-1)
                YDAA(i,k,j) = YAA(LqT-iDD1/2+i-1,Ly-iDD2/2+j-1,
     >            LQ-iDD3/2+k-1)
              End If
              If (Testing) then
cc                Print*, ' i,j,k, YDAS(i,j,k),YDAA(i,j,k) = ',
cc     >                    i,j,k, YDAS(i,j,k), YDAA(i,j,k)
                    Print*, i,j,k,
     >                  X3A(LQ-iDD3/2+k-1), X1A(LqT-iDD1/2+i-1),
     >                  X2A(Ly-iDD2/2+j-1), YDAS(k,i,j),YDAA(k,i,j)
              End If
            End Do
          End Do
        End Do

CsB * Interpolate
      If (Testing) then
        Print*, ' X1A(LqT-iDD1/2), X2A(Ly-iDD2/2), X3A(LQ-iDD3/2) ='
        Print*,   X1A(LqT-iDD1/2), X2A(Ly-iDD2/2), X3A(LQ-iDD3/2)
        Print*, ' Q, qT, y = ', Qv, qT, y
        Print*, ' iDD1, iDD2, iDD3 = ', iDD1, iDD2, iDD3
      EndIf

      If (nyEdge.Eq.0) then
CCAO    If (nPart.Eq.11 .or. nPart.Eq.12) then ! CSS or Asy
        If (nPart.Eq.11 .or. nPart.Eq.12 .or.
     .      nPart.eq.15) then ! CSS or Asy
CsB   Note: Polin is called this way so that y gets interpolated
C     first, then Q, and finally qT.
          Call Polin3 ( X1A(LqT-iDD1/2), X3A(LQ-iDD3/2), X2A(Ly-iDD2/2),
     >           YDAS, iDD1,iDD3,iDD2_FIT, qT,Qv,y, DsigS,DErrSgS )
CCAO      
!         IF(IFLAG.EQ.-997) GOTO 33
          Call Polin3 ( X1A(LqT-iDD1/2), X3A(LQ-iDD3/2), X2A(Ly-iDD2/2),
     >           YDAA, iDD1,iDD3,iDD2_FIT, qT,Qv,y, DsigA,DErrSgA )
CCAO      
!         IF(IFLAG.EQ.-997) GOTO 33       
        Else If (nPart.Eq.13 .or. nPart.Eq.14 .or.
     .           nPart.Eq.16) then ! DSi or LO
CsB   Note: Polin is called this way so that y gets interpolated
C     first, then qT, and finally Q.
          Call Polin3 ( X3A(LQ-iDD3/2), X1A(LqT-iDD1/2), X2A(Ly-iDD2/2),
     >           YDAS, iDD3,iDD1,iDD2_FIT, Qv,qT,y, DsigS,DErrSgS )
CCAO      
!         IF(IFLAG.EQ.-997) GOTO 33       
          Call Polin3 ( X3A(LQ-iDD3/2), X1A(LqT-iDD1/2), X2A(Ly-iDD2/2),
     >           YDAA, iDD3,iDD1,iDD2_FIT, Qv,qT,y, DsigA,DErrSgA )
CCAO      
!         IF(IFLAG.EQ.-997) GOTO 33       
        End If
      EndIf
      If (nyEdge.Eq.2) then
        Call Polin3m( X3A(LQ-iDD3/2), X1A(LqT-iDD1/2), X2A(Ly-iDD2/2),
     >                YDAS, iDD3, iDD1, iDD2, Qv, qT, y, DsigS,DErrSgS )
        Call Polin3m( X3A(LQ-iDD3/2), X1A(LqT-iDD1/2), X2A(Ly-iDD2/2),
     >                YDAA, iDD3, iDD1, iDD2, Qv, qT, y, DsigA,DErrSgA )
      EndIf
      SigS   = DSigS
      ErrSgS = DErrSgS
      SigA   = DSigA
      ErrSgA = DErrSgA

CsB * Take special care of the p p case when the L0 rapidity spectrum is
C     symmetric and A3 is anti-symmetric wrt y -> -y.
CsB   Commented out: 2005 02 15.  Not needed when grid spans both signs of y.
CCPY SEPT 2009:

      iSignY = 1.d0
      If (yIn.NE.0.d0) iSignY = y/yIn

      IF(BOSON.NE.'AA'.AND.BOSON.NE.'AG'.AND.
     >BOSON.NE.'ZZ'.AND.BOSON.NE.'ZG') THEN
CCPY Sept 2009: 
C The meaning of the column under "SigA" is to store 
C     the CSS delta(1-z) piece for 'AA", 'AG', 'ZZ" And 'ZG'.
C In other cases, it stores the coefficient for the 
C angular function 2*cos(theta).
     
        If (Y_ONE_SIDE) SigA = iSignY*DSigA

      ENDIF
      
      If (Testing) then
        Print*, ' SigS = ', SigS, ' +-', ErrSgS
        Print*, ' SigA = ', SigA, ' +-', ErrSgA
      EndIf

33    CONTINUE

      Return
      End

CsB * *******************************************************************

      Subroutine ThreeYInterp ( iDD1,iDD2,iDD3, YDAS,YDAA, qT,yIn,Qv,
     >                          SigY, ErrSgS )
C
C     The following notation apply:
C     X1A (X2A, X3A) contains the qT (y, Q) values
C     iD1+1 (iD2+1, iD3+1) is the number of qT (y, Q) values in the grid
C     iDD1 (iDD2, iDD3) is the number of qT (y, Q)
C          values used at the interpolation
C     The 1st indices of YAS and YAA belong to qT, the 2nd to y, the 3rd to Q.
      Implicit Double Precision (A-H,O-Z)
C     The dimensions below give the maximal possible gridsize.
      Parameter (MD1 = 200, MD2 = 300, MD3 = 100)
      Dimension X1A(0:MD1), X2A(0:MD2), X3A(0:MD3)
      Dimension YAS(0:MD1,0:MD2,0:MD3), YAA(0:MD1,0:MD2,0:MD3)
      Dimension YDAS(iDD3,iDD1,iDD2), YDAA(iDD3,iDD1,iDD2)
      Dimension X1Y(0:MD1), X2Y(0:MD2), X3Y(0:MD3)
      Dimension YAY(0:MD1,0:MD2,0:MD3,0:7)
      Dimension SigY (0:5)
      Logical Testing
      Common /ECMC/ ECM
      Common / Arrays / X1A,X2A,X3A, YAS,YAA,  X1Y,X2Y,X3Y, YAY
      Common / Dims / iD1,iD2,iD3, iY1,iY2,iY3
      Character*40 Boson
      Common / Boson / Boson
!      Integer N_SUD_AB,N_SUD_A,N_SUD_B,N_WIL_C,I_SUD
!      Common / Sudakov / N_SUD_AB,N_SUD_A,N_SUD_B,N_WIL_C,I_SUD
      Common / nyEdge / nyEdge
      COMMON/HARDCUT/ QT_V_MIN, QT_V_MAX,YB_V_MAX,
     >  Q_V_MIN,Q_V_MAX,QT_V_MXY
      Common / Matching / QT_V_MAT
      COMMON/AA/W,QCDL,X1X2MIN,NJET,N_PRIMARY_PART,NF,I_STR,IBEAM,NPART,
     &  NEVENT,NORDER,N_SECOND_PART

CCPY
      LOGICAL Y_ONE_SIDE
      COMMON/GRID_FILE_FORM/ Y_ONE_SIDE

      Common / TI / nTI
      Integer iYPSw
      Common / YPertSwitch / iYPSw
      Parameter (yMx = 5.d0)
      Data eps / 1.d-9 /
      Data Testing /.False./

cdump
C      Testing = .True.

      DO I=0,5
        SigY(I)=0.D0
      ENDDO

CsB Set flag that tells Polin3 that it was called from ThreeYInterp
      nTI = 1

CsB * Testing
      If (Testing) then
        Print*, ' iY1, iY2, iY3 = ', iY1, iY2, iY3
        Print*, ' iDD1/2, iDD2/2, iDD3/2 = ', iDD1/2, iDD2/2, iDD3/2
C        Do K = 0, iY3+1
C        Do J = 1, iY2-1
C        Do I = 1, iY1
C            Print*, X3Y(K), X1Y(I), X2Y(J)
C            Print*, YAY(I,J,K,5), YAY(I,J,K,6)
C            Print*, YAY(I,J,K,1), YAY(I,J,K,4)
C            Print*, YAY(I,J,K,2), YAY(I,J,K,7)
C        End Do
C        End Do
C        End Do
      EndIf

      y = yIn
CsB * Take care of the photon case when the rapidity spectrum is symmetric.
C     (p N collision is an exception)
CCPY Sept 2009: Add iBeam.Eq.1 option
      If (Boson.Eq.'AA' .and. iBeam.Eq.0) Goto 160

      IF(Y_ONE_SIDE) THEN
        y = Abs(yIn)
      EndIf
  160 Continue

CsB * Take special care of the p p case when the L0 rapidity spectrum is
C     symmetric and A3 is anti-symmetric wrt y -> -y.
C     At the LHC
CsB   Commented out: 2005 02 15.  Not needed when grid spans both signs of y.
C     If (iBeam.Eq.1 .and. Boson.Eq.'W+') y = Abs(yIn)

cCsB * If Q, qT or y right on the edge of the grid
c      If (qT.Eq.X1Y(0)) qT = qT + eps
c      If ( y.Eq.X2Y(0))  y =  y + eps
c      If (Qv.Eq.X3Y(0)) Qv = Qv + eps
c      If (qT.Eq.X1Y(iY1)) qT = qT - eps
c      If ( y.Eq.X2Y(iY2))  y =  y - eps
c      If (Qv.Eq.X3Y(iY3)) Qv = Qv - eps

CsB * Find out the location of Q, qT, y
      call LOCATE(X1Y, iY1, qT, LqT)
      call LOCATE(X2Y, iY2,  y, Ly)
      call LOCATE(X3Y, iY3, Qv, LQ)
      If (Testing) Print *, ' LqT, Ly, LQ = ', LqT, Ly, LQ

CsB * If not on the grid
      If ( Abs((X1Y(iY1)-qT)*(qT-X1Y(0))).LT.0.e-9 .or.
     .     Abs((X2Y(iY2)- y)*( y-X2Y(0))).LT.0.e-9 .or.
     .     Abs((X3Y(iY3)-Qv)*(Qv-X3Y(0))).LT.0.e-9 ) then
        Print*, ' Warning from 3YInterp: Extrapolation used'
C        Print*, ' LqT,Ly,LQ = ', LqT,Ly,LQ
        Print*, ' Q, qT, y = ', Qv,qT,y
C        Print*, ' (qT-X1Y(0))   =',(qT-X1Y(0))
C        Print*, ' (X1Y(iY1)-qT) =',(X1Y(iY1)-qT)
C        Print*, ' (X2Y(iY2)- y) =',(X2Y(iY2)- y)
C        Print*, ' ( y-X2Y(0))   =',( y-X2Y(0))
C        Print*, ' (X3Y(iY3)-Qv)  =',(X3Y(iY3)-Qv)
C        Print*, ' (Qv-X3Y(0))    =',(Qv-X3Y(0))
C        Print*, ' >>> ', X3Y(0), X3A(0)
C        Print*, X3Y(iY3), X1Y(iY1)
C        Print*, X2Y
c        Print*, YAY(I,J,K,5), YAY(I,J,K,6)
c        Print*, YAY(I,J,K,1), YAY(I,J,K,2), YAY(I,J,K,4)
c      Stop
      EndIf

CsB * If close to the ends
CsB   ( This works for even iDD for sure )
      If (LqT.LT.iDD1/2)     LqT = iDD1/2
      If (LqT.GT.iY1-iDD1/2) LqT = iY1-iDD1/2+1
      If ( Ly.LT.iDD2/2)      Ly = iDD2/2
      If ( Ly.GT.iY2-iDD2/2)  Ly = iY2-iDD2/2+1
      If ( LQ.LT.iDD3/2)      LQ = iDD3/2
      If ( LQ.GT.iY3-iDD3/2)  LQ = iY3-iDD3/2+1
      If (Testing) Print *, ' LqT, Ly, LQ = ', LqT, Ly, LQ

CsB L0 and A1 pieces -------------------------------------------------------
CsB If iYPSw=0 then send the Y piece back
C   If iYPSw=1 then send the Pert piece back
      nPiece = 0                        ! regular L0 piece (Y piece)
CsB      If (qT.GT.QT_V_MAT) nPiece = 5 ! perturbative L0 piece
      If (iYPSw.Eq.1) nPiece = 5 ! perturbative L0 piece
CsB   The YD arrays filled up this way so that y gets interpolated
CsB   first, then qT, and finally Q
      Do i = 1,iDD1
        Do j = 1,iDD2
          Do k = 1,iDD3
            YDAS(i,k,j) = YAY(LqT-iDD1/2+i-1,Ly-iDD2/2+j-1,
     >        LQ-iDD3/2+k-1,nPiece)
            YDAA(i,k,j) = YAY(LqT-iDD1/2+i-1,Ly-iDD2/2+j-1,
     >        LQ-iDD3/2+k-1,1)
            If (Testing) then
c              Print*, ' i,j,k, YDAS(i,k,j),YDAA(i,k,j) = ',
c     >                  i,j,k, YDAS(i,k,j), YDAA(i,k,j)
C                  Print*, X3Y(LQ-iDD3/2+k-1), X1Y(LqT-iDD1/2+i-1),
C     >                X2Y(Ly-iDD2/2+j-1), YDAS(i,k,j),YDAA(i,k,j)
            EndIf
          EndDo
        EndDo
      EndDo

CsB * Interpolate
      If (Testing) then
        Print*, ' X1Y(LqT-iDD1/2), X2Y(Ly-iDD2/2), X3Y(LQ-iDD3/2) ='
        Print*,   X1Y(LqT-iDD1/2), X2Y(Ly-iDD2/2), X3Y(LQ-iDD3/2)
        Print*, ' Q, qT, y =', Qv, qT, y
c          Print*, ' iDD1, iDD2, iDD3 = ', iDD1, iDD2, iDD3
      EndIf
      If (nyEdge.Eq.0) then
CsB * Note: Polin is called this way so that y gets interpolated
C     first, then Q, and finally qT.
        Call Polin3 ( X1Y(LqT-iDD1/2), X3Y(LQ-iDD3/2), X2Y(Ly-iDD2/2),
     >                YDAS, iDD1, iDD3, iDD2, qT, Qv, y, DsigS,DErrSgS )
        Call Polin3 ( X1Y(LqT-iDD1/2), X3Y(LQ-iDD3/2), X2Y(Ly-iDD2/2),
     >                YDAA, iDD1, iDD3, iDD2, qT, Qv, y, DsigA,DErrSgA )
CsB * Note: Polin is called this way so that y gets interpolated
C     first, then qT, and finally Q.
C        Call Polin3 ( X3A(LQ-iDD3/2), X1A(LqT-iDD1/2), X2A(Ly-iDD2/2),
C     >                YDAS, iDD3, iDD1, iDD2, Qv, qT, y, DsigS,DErrSgS )
C        Call Polin3 ( X3A(LQ-iDD3/2), X1A(LqT-iDD1/2), X2A(Ly-iDD2/2),
C     >                YDAA, iDD3, iDD1, iDD2, Qv, qT, y, DsigA,DErrSgA )
      EndIf
      SigY(0) = DSigS
      ErrSgY  = DErrSgS
      SigY(1) = DSigA
      If (DErrSgA.GT.ErrSgY) ErrSgY  = DErrSgA

CCPY Sept 2009: A1 does not change sign.
C===========
CCPY!CsB * Take special care of the p p case when the L0 rapidity spectrum is
CCPY!C     symmetric and A1 is anti-symmetric wrt y -> -y.
CCPY!C     For W+ at LHC
CCPY!      iSignY = 1.d0
CCPY!      If (yIn.NE.0.d0) iSignY = y/yIn
CCPY!CsB   Commented out: 2005 02 15.  Not needed when grid spans both signs of y.
CCPY!C     If (iBeam.Eq.1 .and. Boson.Eq.'W+') SigY(1) = iSignY*DSigA
CCPY!
C===========

      If (Testing) then
        Print*, ' SigY(0) = ', SigY(0), ' +-', ErrSgY
        Print*, ' SigY(1) = ', SigY(1), ' +-', ErrSgY
      EndIf

      If (Testing) then
        Print*, 'Boson,iYPSw,nyEdge = ',Boson,iYPSw,nyEdge 
      EndIf


ccpy sept 2009      If (Boson.Eq.'ZZ') Return
      If (Boson.Eq.'ZZ' .OR. Boson.Eq.'ZG') Return

CsB A0 = A2 and A3 pieces --------------------------------------------------
CsB If iYPSw=0 then send the Y piece back
C   If iYPSw=1 then send the Pert piece back
C   (simple matching).
      nPiece = 3                     ! regular A3 piece (A3 of Y piece)
CsB      If (qT.GT.QT_V_Mat) nPiece = 6 ! perturbative A3 piece
      If (iYPSw.Eq.1) nPiece = 6 ! perturbative A3 piece
      Do i = 1,iDD1
        Do j = 1,iDD2
          Do k = 1,iDD3
            YDAS(i,k,j) = YAY(LqT-iDD1/2+i-1,Ly-iDD2/2+j-1,
     >        LQ-iDD3/2+k-1,2)
            YDAA(i,k,j) = YAY(LqT-iDD1/2+i-1,Ly-iDD2/2+j-1,
     >        LQ-iDD3/2+k-1,nPiece)
            If (Testing) then
c              Print*, ' i,j,k, YDAS(i,k,j),YDAA(i,k,j) = ',
c     >                  i,j,k, YDAS(i,k,j), YDAA(i,k,j)
c                  Print*, X3A(LQ-iDD3/2+k-1), X1A(LqT-iDD1/2+i-1),
c     >                X2A(Ly-iDD2/2+j-1), YDAS(i,k,j),YDAA(i,k,j)
            EndIf
          EndDo
        EndDo
      EndDo

CsB * Interpolate
      If (Testing) then
c        Print*, ' X1Y(LqT-iDD1/2), X2Y(Ly-iDD2/2), X3Y(LQ-iDD3/2) ='
c        Print*,   X1Y(LqT-iDD1/2), X2Y(Ly-iDD2/2), X3Y(LQ-iDD3/2)
c          Print*, ' iDD1, iDD2, iDD3 = ', iDD1, iDD2, iDD3
      EndIf

      If (nyEdge.Eq.0) then
CsB   Note: Polin is called this way so that y gets interpolated
CsB   first, then Q, and finally qT.
        Call Polin3 ( X1Y(LqT-iDD1/2), X3Y(LQ-iDD3/2), X2Y(Ly-iDD2/2),
     >                YDAS, iDD1, iDD3, iDD2, qT, Qv, y, DsigS,DErrSgS )
CsB_____This is a general problem (should be fixed before this stage).
C       Reason for this fix here: 
C       Polint/Ratint cannot interpolate on an all zero Y vector.
C       If (Boson.Eq.'AA' .or. Boson.Eq.'AG') Goto 199
C       If (Boson.Eq.'AA') Goto 199
C       As of 01/2004 Polin3 was fixed to take care of this.
        Call Polin3 ( X1Y(LqT-iDD1/2), X3Y(LQ-iDD3/2), X2Y(Ly-iDD2/2),
     >                YDAA, iDD1, iDD3, iDD2, qT, Qv, y, DsigA,DErrSgA )
 199  Continue
CsB   Note: Polin is called this way so that y gets interpolated
CsB   first, then qT, and finally Q.
C        Call Polin3 ( X3Y(LQ-iDD3/2), X1Y(LqT-iDD1/2), X2Y(Ly-iDD2/2),
C     >                YDAS, iDD3, iDD1, iDD2, Qv, qT, y, DsigS,DErrSgS )
C        Call Polin3 ( X3Y(LQ-iDD3/2), X1Y(LqT-iDD1/2), X2Y(Ly-iDD2/2),
C     >                YDAA, iDD3, iDD1, iDD2, Qv, qT, y, DsigA,DErrSgA )
      EndIf
      SigY(2) = DSigS
      If (DErrSgS.GT.ErrSgY) ErrSgY  = DErrSgS
      SigY(3) = DSigA
      If (DErrSgA.GT.ErrSgY) ErrSgY  = DErrSgA

CsB * Take special care of the p p case when the L0 rapidity spectrum is
C     symmetric and A3 is anti-symmetric wrt y -> -y.
CsB   Commented out: 2005 02 15.  Not needed when grid spans both signs of y.
CCPY SEPT 2009:

      iSignY = 1.d0
      If (yIn.NE.0.d0) iSignY = y/yIn
      If (Y_ONE_SIDE) SigY(3) = iSignY*DSigA

      If (Testing) then
        Print*, ' SigY(2) = ', SigY(2), ' +-', ErrSgY
        Print*, ' SigY(3) = ', SigY(3), ' +-', ErrSgY
      EndIf

CsB      If (iYPSw.Eq.1) Return ! Only the Pert L0 and A3 pieces are needed

CsB A4 piece  --------------------------------------------------------------
      Do i = 1,iDD1
        Do j = 1,iDD2
          Do k = 1,iDD3
            YDAS(i,k,j) = YAY(LqT-iDD1/2+i-1,Ly-iDD2/2+j-1,
     >        LQ-iDD3/2+k-1,4)
            If (Testing) then
              Print*, ' i,j,k, YDAS(i,k,j) = ',
     >                  i,j,k, YDAS(i,k,j)
                  Print*, X3A(LQ-iDD3/2+k-1), X1A(LqT-iDD1/2+i-1),
     >                X2A(Ly-iDD2/2+j-1), YDAS(i,k,j)
            EndIf
          EndDo
        EndDo
      EndDo

CsB * Interpolate
      If (Testing) then
c        Print*, ' X1Y(LqT-iDD1/2), X2Y(Ly-iDD2/2), X3Y(LQ-iDD3/2) ='
c        Print*,   X1Y(LqT-iDD1/2), X2Y(Ly-iDD2/2), X3Y(LQ-iDD3/2)
c          Print*, ' iDD1, iDD2, iDD3 = ', iDD1, iDD2, iDD3
      EndIf
      If (nyEdge.Eq.0) then
CsB * Note: Polin is called this way so that y gets interpolated
C     first, then Q, and finally qT.
        Call Polin3 ( X1Y(LqT-iDD1/2), X3Y(LQ-iDD3/2), X2Y(Ly-iDD2/2),
     >                YDAS, iDD1, iDD3, iDD2, qT, Qv, y, DsigS,DErrSgS )
CCPY Sept 2009: Remove a few wrong codes below.

      EndIf

      SigY(4) = DSigS
      If (DErrSgA.GT.ErrSgY) ErrSgY = DErrSgS

CsB * Take special care of the p p case when the L0 rapidity spectrum is
C     symmetric and A3 is anti-symmetric wrt y -> -y.
CsB   Commented out: 2005 02 15.  Not needed when grid spans both signs of y.
CCPY SEPT 2009:

      iSignY = 1.d0
      If (yIn.NE.0.d0) iSignY = y/yIn
      If (Y_ONE_SIDE) SigY(4) = iSignY*DSigS


      If (Testing) then
        Print*, ' SigY(4) = ', SigY(4), ' +-', ErrSgY
      EndIf

CsB A0 piece  --------------------------------------------------------------
      Do i = 1,iDD1
        Do j = 1,iDD2
          Do k = 1,iDD3
            YDAS(i,k,j) = YAY(LqT-iDD1/2+i-1,Ly-iDD2/2+j-1,
     >        LQ-iDD3/2+k-1,7)
            If (Testing) then
              Print*, ' i,j,k, YDAS(i,k,j) = ',
     >                  i,j,k, YDAS(i,k,j)
                  Print*, X3A(LQ-iDD3/2+k-1), X1A(LqT-iDD1/2+i-1),
     >                X2A(Ly-iDD2/2+j-1), YDAS(i,k,j)
            EndIf
          EndDo
        EndDo
      EndDo

CsB * Interpolate
      If (Testing) then
c        Print*, ' X1Y(LqT-iDD1/2), X2Y(Ly-iDD2/2), X3Y(LQ-iDD3/2) ='
c        Print*,   X1Y(LqT-iDD1/2), X2Y(Ly-iDD2/2), X3Y(LQ-iDD3/2)
c          Print*, ' iDD1, iDD2, iDD3 = ', iDD1, iDD2, iDD3
      EndIf
      If (nyEdge.Eq.0) then
CsB * Note: Polin is called this way so that y gets interpolated
C     first, then Q, and finally qT.
        Call Polin3 ( X1Y(LqT-iDD1/2), X3Y(LQ-iDD3/2), X2Y(Ly-iDD2/2),
     >                YDAS, iDD1, iDD3, iDD2, qT, Qv, y, DsigS,DErrSgS )
CCPY Sept 2009: Remove a few wrong codes below.

      EndIf

      SigY(5) = DSigS
      If (DErrSgA.GT.ErrSgY) ErrSgY = DErrSgS

CsB * Take special care of the p p case when the L0 rapidity spectrum is
C     symmetric and A3 is anti-symmetric wrt y -> -y.
CsB   Commented out: 2005 02 15.  Not needed when grid spans both signs of y.
CCPY SEPT 2009:

      iSignY = 1.d0
      If (yIn.NE.0.d0) iSignY = y/yIn
      If (Y_ONE_SIDE) SigY(5) = iSignY*DSigS


      If (Testing) then
        Print*, ' SigY(5) = ', SigY(5), ' +-', ErrSgY
      EndIf

      Return
      End

CsB * *******************************************************************

      Subroutine ThreeKInterp ( iDD1,iDD2,iDD3, YDAS,YDAA, qT,yIn,Q,
     >                          sigS,sigA, ErrSgS,ErrSgA )
C
C     The following convention apply:
C     X1A (X2A, X3A) contains the qT (y, Q) vector
C     iD1+1 (iD2+1, iD3+1) is the number of qT (y, Q) values in the grid
C     iDD1 (iDD2, iDD3) is the number of qT (y, Q)
C          values used at the interpolation
C     The 1st indices of YAS and YAA belong to qT, the 2nd to y, the 3rd to Q.

CsB   Note: Arguments passed through common block:
C     X1A, X2A and X3A: the arrrays that hold the coordinate points
C     YAS and YAA: the arrrays on which we interpolate

      Implicit Double Precision (A-H,O-Z)
C     The dimensions below give the maximal possible gridsize.
      Parameter (MD1 = 200, MD2 = 300, MD3 = 100)
      Dimension X1A(0:MD1), X2A(0:MD2), X3A(0:MD3)
      Dimension YAS(0:MD1,0:MD2,0:MD3), YAA(0:MD1,0:MD2,0:MD3)
      Dimension YDAS(iDD3,iDD1,iDD2), YDAA(iDD3,iDD1,iDD2)
      Dimension X1Y(0:MD1), X2Y(0:MD2), X3Y(0:MD3)
      Dimension YAY(0:MD1,0:MD2,0:MD3,0:7)
      Logical Testing
      Common /ECMC/ ECM
      Common / ArraysK / X1A,X2A,X3A, YAS,  X1Y,X2Y,X3Y, YAA
      Common / DimsK / iD1,iD2,iD3, iY1,iY2,iY3
      Character*40 Boson
      Common / Boson / Boson
!      Integer N_SUD_AB,N_SUD_A,N_SUD_B,N_WIL_C,I_SUD
!      Common / Sudakov / N_SUD_AB,N_SUD_A,N_SUD_B,N_WIL_C,I_SUD
      Common / nyEdge / nyEdge
      COMMON/HARDCUT/ QT_V_MIN, QT_V_MAX,YB_V_MAX,
     >  Q_V_MIN,Q_V_MAX,QT_V_MXY
      Integer LTOpt
      Common / Process / LTOpt
      Common / Matching / QT_V_Mat
      COMMON/AA/W,QCDL,X1X2MIN,NJET,N_PRIMARY_PART,NF,I_STR,IBEAM,NPART,
     &  NEVENT,NORDER,N_SECOND_PART
 
CCPY
      LOGICAL Y_ONE_SIDE
      COMMON/GRID_FILE_FORM/ Y_ONE_SIDE

      Common / TI / nTI
      Parameter (yMx = 5.d0)
      Data eps / 1.d-9 /
      Data Testing / .False. /

C      Testing = .True.
C      Print*, X3A(2),X2A(2),X1A(2),YAS(2,2,2),YAA(2,2,2)
C      Print*, X3A(2),X2A(3),X1A(4),YAS(4,3,2),YAA(4,3,2)
C      Print*, X3A(4),X2A(3),X1A(2),YAS(2,3,4),YAA(2,3,4)
C      Stop


CsB Set flag that tells Polin3 that it was called from ThreeKInterp
      nTI = 0

CsB * Testing
C      If (.True.) then
      If (Testing) then
        Print*, ' iD1, iD2, iD3 = ', iD1, iD2, iD3
        Print*, ' iDD1/2, iDD2/2, iDD3/2 = ', iDD1/2, iDD2/2, iDD3/2
      EndIf

      y = yIn
CsB * Take care of the photon case when the rapidity spectrum is symmetric.
C     (p N collision is an exception)
CCPY Sept 2009: Add iBeam.Eq.1 option
      If (Boson.Eq.'AA' .and. iBeam.Eq.0) Goto 160

      IF(Y_ONE_SIDE) THEN
        y = Abs(yIn)
      End If
  160 Continue

CsB * Take special care of the p p case when the L0 rapidity spectrum is
C     symmetric and A3 is anti-symmetric wrt y -> -y.
C     For W+ at LHC
CsB   Commented out: 2005 02 15.  Not needed when grid spans both signs of y.
C     If (iBeam.Eq.1 .and. Boson.Eq.'W+') y = Abs(yIn)

CsB * If Q, qT or y right on the edge of the grid
      If (qT.Eq.X1A(0)) qT = qT + eps
      If ( y.Eq.X2A(0))  y =  y + eps
      If ( Q.Eq.X3A(0))  Q =  Q + eps
      If (qT.Eq.X1A(iD1)) qT = qT - eps
      If ( y.Eq.X2A(iD2))  y =  y - eps
      If ( Q.Eq.X3A(iD3))  Q =  Q - eps

CsB * Find out the location of Q, qT, y
      call LOCATE(X1A, iD1, qT, LqT)
      call LOCATE(X2A, iD2,  y, Ly)
      call LOCATE(X3A, iD3,  Q, LQ)
      If (Testing) Print *, ' LqT, Ly, LQ = ', LqT, Ly, LQ
C      If (.True.) Print *, ' LqT, Ly, LQ = ', LqT, Ly, LQ

CsB * If not on the grid
      If ( Abs((X1A(iD1)-qT)*(qT-X1A(0))).LT.0.e-9 .or.
     .     Abs((X2A(iD2)- y)*( y-X2A(0))).LT.0.e-9 .or.
     .     Abs((X3A(iD3)- Q)*( Q-X3A(0))).LT.0.e-9 ) then
        Print*, ' Warning from 3KInterp: Extrapolation used'
        Print*, ' Q, qT, y = ', Q,qT,y
        Print*, ' LQ,LqT,Ly= ', LQ,LqT,Ly
        Print*, ' iD3,iD1,iD2= ', iD3,iD1,iD2
      EndIf

CsB * If close to the ends
CsB   ( This works for even iDD for sure )
      If (LqT.LT.iDD1/2)     LqT = iDD1/2
      If (LqT.GT.iD1-iDD1/2) LqT = iD1-iDD1/2+1
      If ( Ly.LT.iDD2/2)      Ly = iDD2/2
      If ( Ly.GT.iD2-iDD2/2)  Ly = iD2-iDD2/2+1
      If ( LQ.LT.iDD3/2)      LQ = iDD3/2
      If ( LQ.GT.iD3-iDD3/2)  LQ = iD3-iDD3/2+1
      If (Testing) Print *, ' LqT, Ly, LQ = ', LqT, Ly, LQ

CsB * Take care of the end points of the rapidity - if necessary
C     ( The right way to do it would be:
C       If the interpolation interval includes yMax
C       than (within Polin3) move the closest point of X2A
C       where the xsect = 0 into yMax.
C       The problem with the above approach is that if the interpolation
C       point will fall very close to the edge, then Polint gives back
C       a weird negative value. So may be, we should set a limit how close
C       to the edge we interpolate using Polint. )
      If (nyEdge.EQ.1) then
        If (Ly.LT.5 .or. Ly.GT.iD2-5) then
          Call YMAXIMUM (ECM,Q,qT,yMax)
          If (Testing) then
            Call Locate (X2A, iD2, yMax, LyMax)
            Print*, ' yMax, LyMax = ', yMax, LyMax
          EndIf
          Print*, ' Fix yMx. Should be sg like X2A(iD2) ... '
          Call ExitGracefully
          If (Abs(yMax).LT.yMx) then
            X2A(0)  = -yMx
            X2A(1)  = -yMax
            X2A(13) =  yMax
            X2A(14) =  yMx
          Else
            X2A(0)  = -yMax
            X2A(1)  = -yMx
            X2A(13) =  yMx
            X2A(14) =  yMax
          End If
        End If
      End If

CsB * Fill the dummy array(s) you want to interpolate
C     ( This works for even iDD for sure )
C     Since the LO and DSi pieces have only onr qT (no qT integral in the
C     phase space the order of interpolation [1->y 2->qT 3->Q] is preserved
C     so that the interpolation does not have to be changed.)

cCsB   The YD arrays filled up this way so that y gets interpolated
cCsB   first, then qT, and finally Q
c        Do i = 1,iDD1
c          Do j = 1,iDD2
c            Do k = 1,iDD3
c              YDAS(k,i,j) = YAS(LqT-iDD1/2+i-1,Ly-iDD2/2+j-1,
c     >          LQ-iDD3/2+k-1)
c              YDAA(k,i,j) = YAA(LqT-iDD1/2+i-1,Ly-iDD2/2+j-1,
c     >          LQ-iDD3/2+k-1)
c              If (Testing) then
ccc                Print*, ' i,j,k, YDAS(i,j,k),YDAA(i,j,k) = ',
ccc     >                    i,j,k, YDAS(i,j,k), YDAA(i,j,k)
cc                    Print*, X3A(LQ-iDD3/2+k-1), X1A(LqT-iDD1/2+i-1),
cc     >                  X2A(Ly-iDD2/2+j-1), YDAS(k,i,j),YDAA(k,i,j)
c              End If
c            End Do
c          End Do
c        End Do

        Do i = 1,iDD1
          Do j = 1,iDD2
            Do k = 1,iDD3
CsB   Note: the YD arrays filled up this way so that y gets interpolated
CsB   first, then qT, and finally Q.
c              If (nPart.Eq.13 .or. nPart.Eq.14) then ! DSi or LO
c                YDAS(k,i,j) = YAS(LqT-iDD1/2+i-1,Ly-iDD2/2+j-1,
c     >            LQ-iDD3/2+k-1)
c                YDAA(k,i,j) = YAA(LqT-iDD1/2+i-1,Ly-iDD2/2+j-1,
c     >            LQ-iDD3/2+k-1)
c                Print*, X3A(LQ-iDD3/2+k-1), X2A(Ly-iDD2/2+j-1),
c     >                  X1A(LqT-iDD1/2+i-1), YDAS(k,i,j),YDAA(k,i,j)
c              Else If (nPart.Eq.11 .or. nPart.Eq.12) then ! CSS or Asy
CsB   Note: the YD arrays filled up this way so that y gets interpolated
CsB   first, then Q, and finally qT.
                YDAS(i,k,j) = YAS(LqT-iDD1/2+i-1,Ly-iDD2/2+j-1,
     >            LQ-iDD3/2+k-1)
                YDAA(i,k,j) = YAA(LqT-iDD1/2+i-1,Ly-iDD2/2+j-1,
     >            LQ-iDD3/2+k-1)
c               End If
              If (.True.) then
C              If (Testing) then
C                Print*, ' i,j,k, YDAS(i,j,k),YDAA(i,j,k) = ',
C     >                    i,j,k, YDAS(i,j,k), YDAA(i,j,k)
C                    Print*, X3A(LQ-iDD3/2+k-1), X1A(LqT-iDD1/2+i-1),
C     >                  X2A(Ly-iDD2/2+j-1), YDAS(i,k,j),YDAA(i,k,j)
              End If
            End Do
          End Do
        End Do
ccc                If (X2A(Ly-iDD2/2).Eq.-1.d0) Stop

CsB * Interpolate
      If (Testing) then
        Print*, ' X1A(LqT-iDD1/2), X2A(Ly-iDD2/2), X3A(LQ-iDD3/2) ='
        Print*,   X1A(LqT-iDD1/2), X2A(Ly-iDD2/2), X3A(LQ-iDD3/2)
        Print*, ' Q, qT, y = ', Q, qT, y
c          Print*, ' iDD1, iDD2, iDD3 = ', iDD1, iDD2, iDD3
      EndIf
c      If (nyEdge.Eq.0) then
c        If (nPart.Eq.11 .or. nPart.Eq.12) then ! CSS or Asy
CsB   Note: Polin is called this way so that y gets interpolated
C     first, then Q, and finally qT.
          Call Polin3 ( X1A(LqT-iDD1/2), X3A(LQ-iDD3/2), X2A(Ly-iDD2/2),
     >                  YDAS, iDD1,iDD3,iDD2, qT,Q,y, DsigS,DErrSgS )
          Call Polin3 ( X1A(LqT-iDD1/2), X3A(LQ-iDD3/2), X2A(Ly-iDD2/2),
     >                  YDAA, iDD1,iDD3,iDD2, qT,Q,y, DsigA,DErrSgA )
c        Else If (nPart.Eq.13 .or. nPart.Eq.14) then ! DSi or LO
cCsB   Note: Polin is called this way so that y gets interpolated
cC     first, then qT, and finally Q.
c          Call Polin3 ( X3A(LQ-iDD3/2), X1A(LqT-iDD1/2), X2A(Ly-iDD2/2),
c     >                  YDAS, iDD3,iDD1,iDD2, Q,qT,y, DsigS,DErrSgS )
c          Call Polin3 ( X3A(LQ-iDD3/2), X1A(LqT-iDD1/2), X2A(Ly-iDD2/2),
c     >                  YDAA, iDD3,iDD1,iDD2, Q,qT,y, DsigA,DErrSgA )
c        End If
c      EndIf
c      If (nyEdge.Eq.2) then
c        Call Polin3m( X3A(LQ-iDD3/2), X1A(LqT-iDD1/2), X2A(Ly-iDD2/2),
c     >                YDAS, iDD3, iDD1, iDD2, Q, qT, y, DsigS, DErrSgS )
c        Call Polin3m( X3A(LQ-iDD3/2), X1A(LqT-iDD1/2), X2A(Ly-iDD2/2),
c     >                YDAA, iDD3, iDD1, iDD2, Q, qT, y, DsigA, DErrSgA )
c      EndIf
      SigS   = DSigS
      ErrSgS = DErrSgS
      SigA   = DSigA
      ErrSgA = DErrSgA

CsB * Take special care of the p p case when the L0 rapidity spectrum is
C     symmetric and A3 is anti-symmetric wrt y -> -y.
CsB   Commented out: 2005 02 15.  Not needed when grid spans both signs of y.
CCPY SEPT 2009:

      iSignY = 1.d0
      If (yIn.NE.0.d0) iSignY = y/yIn
      If (Y_ONE_SIDE) SigA = iSignY*DSigA

      If (Testing) then
        Print*, ' SigS = ', SigS, ' +-', ErrSgS
        Print*, ' SigA = ', SigA, ' +-', ErrSgA
      EndIf

      Return
      End

CsB * Vegas **************************************************************

      SUBROUTINE VEGAS(FXN,AVGI,SD,CHI2A)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      REAL*8 XL,XU,ACC,XI,SI,SI2,SWGT,SCHI,ALPH
      INTEGER*4 NDIM,NCALL,ITMX,NPRN,NDO,IT
      COMMON/BVEG1/XL(50),XU(50),ACC,NDIM,NCALL,ITMX,NPRN
      COMMON/BVEG2/XI(50,50),SI,SI2,SWGT,SCHI,NDO,IT,ALPH
      COMMON/SEED/I_SEED
C
      DIMENSION D(50,50),DI(50,50),XIN(50),R(50),DX(50),DT(50),
     1   X(50),KG(50),IA(50)
      DIMENSION RAND(50)
CBOX
      REAL*8 Y_ALPH
      DATA NDMX/50/,Y_ALPH/.2D0/,ONE/1./,MDS/0/

      NDO=1
      DO 1 J=1,NDIM
    1 XI(1,J)=ONE

CBOX
      IF(ALPH.LT.1.D-4) ALPH=Y_ALPH

      ENTRY VEGAS1(FXN,AVGI,SD,CHI2A)
      IT=0
      SI=0.
      SI2=SI
      SWGT=SI
      SCHI=SI
      ENTRY VEGAS2(FXN,AVGI,SD,CHI2A)
      ND=NDMX
      NG=1
      IF(MDS.EQ.0)GO TO 2
      NG=(NCALL/2.)**(1./NDIM)
      MDS=1
      IF((2*NG-NDMX).LT.0)GO TO 2
      MDS=-1
      NPG=NG/NDMX+1
      ND=NG/NPG
      NG=NPG*ND
    2 K=NG**NDIM
      NPG=NCALL/K
      IF(NPG.LT.2)NPG=2
      CALLS=NPG*K
      DXG=ONE/NG
      DV2G=(CALLS*DXG**NDIM)**2/NPG/NPG/(NPG-ONE)
      XND=ND
      NDM=ND-1
      DXG=DXG*XND
      XJAC=ONE/CALLS
      DO 3 J=1,NDIM
        DX(J)=XU(J)-XL(J)
    3 XJAC=XJAC*DX(J)
      IF(ND.EQ.NDO)GO TO 8
      RC=NDO/XND
      DO 7 J=1,NDIM
        K=0
        XN=0.
        DR=XN
        I=K
    4   K=K+1
        DR=DR+ONE
        XO=XN
        XN=XI(K,J)
    5   IF(RC.GT.DR)GO TO 4
        I=I+1
        DR=DR-RC
        XIN(I)=XN-(XN-XO)*DR
        IF(I.LT.NDM)GO TO 5
        DO 6 I=1,NDM
    6   XI(I,J)=XIN(I)
    7 XI(ND,J)=ONE
      NDO=ND
CsB    8 IF(NPRN.GE.0)WRITE(6,200)NDIM,CALLS,IT,ITMX,ACC
C       1   ,MDS,ND,(XL(J),XU(J),J=1,NDIM)
    8 IF(NPRN.GE.0) Then
        WRITE(6,200)NDIM,CALLS,IT,ITMX
        Write(6,  *)
     > ' Iteration, Integral, Accumulated integral, Chi**2/it''n '
      End If
      ENTRY VEGAS3(FXN,AVGI,SD,CHI2A)
    9 IT=IT+1
      TI=0.
      TSI=TI
      DO 10 J=1,NDIM
        KG(J)=1
        DO 10 I=1,ND
          D(I,J)=TI
   10 DI(I,J)=TI
   11 FB=0.
      F2B=FB
      K=0
   12 K=K+1
      CALL YRANDA(NDIM,RAND)
      WGT=XJAC
      DO 15 J=1,NDIM
        XN=(KG(J)-RAND(J))*DXG+ONE
        IA(J)=XN
        IF(IA(J).GT.1)GO TO 13
        XO=XI(IA(J),J)
        RC=(XN-IA(J))*XO
        GO TO 14
   13   XO=XI(IA(J),J)-XI(IA(J)-1,J)
        RC=XI(IA(J)-1,J)+(XN-IA(J))*XO
   14   X(J)=XL(J)+RC*DX(J)
   15 WGT=WGT*XO*XND
      F=WGT
      F=F*FXN(X,WGT)
CsB
C      If (Abs(WGT).LT.1e-9) Print*, ' WGT =', WGT
C      If (Abs(F).LT.1e-9) Print*, ' F =', F
CsB
      F2=F*F
      FB=FB+F
      F2B=F2B+F2
      DO 16 J=1,NDIM
        DI(IA(J),J)=DI(IA(J),J)+F
   16 IF(MDS.GE.0)D(IA(J),J)=D(IA(J),J)+F2
      IF(K.LT.NPG) GO TO 12
  888 FORMAT(1X,'F',G14.6,'F2',G14.6,'FB',G14.6,'F2B',G14.6)
      F2B= SQRT(F2B*NPG)
      F2B=(F2B-FB)*(F2B+FB)
 1661 FORMAT(1X,'F2B',G14.6,'NPG',  I10)
      TI=TI+FB
      TSI=TSI+F2B
   33 FORMAT(1X,'TSI',G14.6,'F2B',G14.6)
      IF(MDS.GE.0)GO TO 18
      DO 17 J=1,NDIM
   17 D(IA(J),J)=D(IA(J),J)+F2B
   18 K=NDIM
   19 KG(K)=MOD(KG(K),NG)+1
      IF(KG(K).NE.1)GO TO 11
      K=K-1
      IF(K.GT.0)GO TO 19
      TSI=TSI*DV2G
      TI2=TI*TI
      WGT=TI2/TSI
      SI=SI+TI*WGT
      SI2=SI2+TI2
      SWGT=SWGT+WGT
      SCHI=SCHI+TI2*WGT
  995 FORMAT(1X,'SWGT',G14.6,'SI2',G14.6)
      AVGI=SI/SWGT
      SD=SWGT*IT/SI2
      CHI2A=SD*(SCHI/SWGT-AVGI*AVGI)/(IT-.999)
      SD=SQRT(ONE/SD)
      IF(NPRN.LT.0)GO TO 21
      TSI=SQRT(TSI)
      WRITE(6,201)IT,TI,TSI,AVGI,SD,CHI2A
      IF(NPRN.EQ.0)GO TO 21
      DO 20 J=1,NDIM
   20 WRITE(6,202) J,(XI(I,J),DI(I,J),D(I,J),I=1,ND)
   21 DO 23 J=1,NDIM
        XO=D(1,J)
        XN=D(2,J)
        D(1,J)=(XO+XN)/2.
        DT(J)=D(1,J)
        DO 22 I=2,NDM
          D(I,J)=XO+XN
          XO=XN
          XN=D(I+1,J)
          D(I,J)=(D(I,J)+XN)/3.
   22   DT(J)=DT(J)+D(I,J)
        D(ND,J)=(XN+XO)/2.
   23 DT(J)=DT(J)+D(ND,J)
      DO 28 J=1,NDIM
        RC=0.
        DO 24 I=1,ND
          R(I)=0.
          IF(D(I,J).LE.0.)GO TO 24
          XO=DT(J)/D(I,J)
          R(I)=((XO-ONE)/XO/LOG(XO))**ALPH
   24   RC=RC+R(I)
        RC=RC/XND
        K=0
        XN=0.
        DR=XN
        I=K
   25   K=K+1
        DR=DR+R(K)
        XO=XN
        XN=XI(K,J)
   26   IF(RC.GT.DR)GO TO 25
        I=I+1
        DR=DR-RC
        XIN(I)=XN-(XN-XO)*DR/R(K)
        IF(I.LT.NDM)GO TO 26
        DO 27 I=1,NDM
   27   XI(I,J)=XIN(I)
   28 XI(ND,J)=ONE
      IF(IT.LT.ITMX.AND.ACC*ABS(AVGI).LT.SD)GO TO 9
CsB  200 FORMAT(1X,'0INPUT PARAMETERS FOR VEGAS:  NDIM=',I3,
C     1  '   NCALL=',F8.0/28X,'  IT=',I5,'    ITMX=',I5/28X,
C     2  '  ACC=',G9.3/28X,'  MDS=',I3,'     ND=',I4/28X,
C     3  '  (XL,XU)=',(T40,'( ',G12.6,' , ',G12.6,' )'))
  200 FORMAT(1X,'0INPUT PARAMETERS FOR VEGAS:  NDIM=',I3,
     2  '   NCALL=',F8.0,', IT=',I5,', ITMX=',I5)
CsB  201 FORMAT(///' INTEGRATION BY VEGAS' / '0ITERATION NO.',I3,
C     1  ':  INTEGRAL=',G14.8/21X,'STD DEV =',G10.4 /
C     2  ' ACCUMULATED RESULTS:   INTEGRAL=',G14.8/
C     3  24X,'STD DEV =',G10.4 / 24X,'CHI**2 PER IT''N =',G10.4)
  201 FORMAT('  ', I4, '  ', G14.8,'+/-', G10.4, '  ',
     >                       G14.8,'+/-', G10.4, '  ', G10.4)
  202 FORMAT(1X,'0DATA FOR AXIS',I2,/,' ',6X,'X',7X,'  DELT I ',
     1  2X,'CONV','CE   ',11X,'X',7X,'  DELT I ',2X,'CONV','CE  '
     2  ,11X,'X',7X,'   DELT I ',2X,'CONV','CE  ',/,
     3  (1X,' ',3G12.4,5X,3G12.4,5X,3G12.4))
      RETURN
      END

CBOX
      SUBROUTINE YRANDA(N,RAND)
      REAL*8 RN
      INTEGER*4 I,N
      REAL*8 RAND
      DIMENSION RAND(50)
      EXTERNAL RN

      DO 10 I=1,N
        RAND(I)=RN(I)
10    CONTINUE

      RETURN
      END

C       FUNCTION RN(IIN)
      FUNCTION RNDM(DUMMY)
C       IBM FORTRAN VERSION OF RANDOM NUMBER GENERATOR   CERN V104
      REAL*8 RN
      COMMON /ZZRNDM/ IY
      IY = IY * 69069
      IF (IY)  5,4,6
    4 IY = 65539
      IX = IY * 69069
      GO TO 6
    5 IY = IY + 2147483647 + 1
    6 YFL = IY
      RNDM = YFL * .4656613E-9
      RN=RNDM
      RETURN
      END


c      SUBROUTINE RDMIN(IXRAN)
c      COMMON /ZZRNDM/ IX
c      IX = IXRAN
c      RETURN
c      ENTRY RDMOUT(IXRAN)
c      IXRAN = IX
c      RETURN
c      END


      FUNCTION RN(IIN)
C      FUNCTION RN_VAX(IIN)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C THIS DOES NOT WORK ON BOX.   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT NONE

      REAL*8 RN_VAX
      INTEGER*4 IIN
      INTEGER*4 I_SEED
      LOGICAL INIT
      DOUBLE PRECISION    RN
      DOUBLE PRECISION    DS(2),    DM(2),    DSEED
      DOUBLE PRECISION    DX24,     DX48
      DOUBLE PRECISION    DL,       DC,       DU,       DR
      COMMON/SEED/I_SEED
      DATA      DS     /  1665 1885., 286 8876.  /
      DATA      DM     /  1518 4245., 265 1554.  /
      DATA      DX24   /  1677 7216.  /
      DATA      DX48   /  281 4749 7671 0656.  /
      DATA INIT/.TRUE./
CCPY MAKE SURE THAT RN() CAN NEVER BE ZERO
1     CONTINUE
      IF(INIT)THEN
        INIT=.FALSE.
        DSEED=FLOAT(I_SEED)
        DS(2)  =  DINT(DSEED/DX24)
        DS(1)  =  DSEED - DS(2)*DX24
      END IF
      DL  =  DS(1) * DM(1)
      DC  =  DINT(DL/DX24)
      DL  =  DL - DC*DX24
      DU  =  DS(1)*DM(2) + DS(2)*DM(1) + DC
      DS(2)  =  DU - DINT(DU/DX24)*DX24
      DS(1)  =  DL
      RN=  (DS(2)*DX24 + DS(1)) / DX48
      IF(RN.EQ.0.D0) GOTO 1

      RETURN
      END

      SUBROUTINE SAVE(NDIM,SAVE_FILE)
      IMPLICIT NONE
      INTEGER*4 I,J,NBINIT,NDO,IT,NDIM,iwgt
      REAL *8 SI,SI2,SWGT
      REAL*8 MAX_WGT,NEW_MAX_WGT,WGT_FACTOR,XI,SCHI,ALPH
      CHARACTER*(*) SAVE_FILE,RESTORE_FILE
      LOGICAL I_READIN,I_WRITEOUT
      COMMON/BVEG2/XI(50,50),SI,SI2,SWGT,SCHI,NDO,IT,ALPH
      COMMON/UNWEIGHT/MAX_WGT,NEW_MAX_WGT,WGT_FACTOR,iwgt,NBINIT
      OPEN(UNIT=7,FILE=SAVE_FILE,FORM='FORMATTED',STATUS='NEW')
      PRINT *, NDO,IT,SI,SI2,SWGT,SCHI
CCPY
      WRITE(7,*) NDO,IT,SI,SI2,SWGT,SCHI,
     1       ((XI(J,I),I=1,NDIM),J=1,NDO)
      IF(iwgt.eq.1)WRITE(7,*)MAX_WGT
      CLOSE(UNIT=7)
      RETURN
      ENTRY RESTR(NDIM,RESTORE_FILE)
      OPEN(UNIT=8,FILE=RESTORE_FILE,FORM='FORMATTED',STATUS='OLD')
CCPY
      READ(8,*) NDO,IT,SI,SI2,SWGT,SCHI,
     1    ((XI(J,I),I= 1,NDIM),J=1,NDO)
      IF(iwgt.eq.1) READ(8,*)MAX_WGT
      CLOSE(UNIT=8)
      RETURN
      END ! SAVE


CsB *
CsB * LOCATE routine, from Numerical Recipes  **************************
CsB *

      SUBROUTINE LOCATE(XX, N, X, J)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION XX(N)

      JL = 0
      JU = N + 1
  10      IF (JU - JL.GT.1) THEN
        JM = (JU + JL)/2
        IF ((XX(N).GT.XX(1)).EQV.(X.GT.XX(JM))) THEN
          JL = JM
        ELSE
          JU = JM
        ENDIF
      GO TO 10
      ENDIF
      J = JL
      RETURN
      END

CsB * Polin3's **************************************************************
C ---------------------------------------------------------------------------
      Subroutine Polin3 (X1A,X2A,X3A,YA,L,M,N,X1,X2,X3,Y,DY)
C ---------------------------------------------------------------------------
C * Modified from Polin2 of Numerical recepies by Cs.B. 1/95.
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (LMAX=20, MMAX=20, NMAX=20)
      DIMENSION X1A(L), X2A(M), X3A(N), YA(L,M,N),
     >          Y1TMP(LMAX), Y2TMP(MMAX), Y3TMP(NMAX)
      Common / ECMC / ECM
      COMMON/AA/W,QCDL,X1X2MIN,NJET,N_PRIMARY_PART,NF,I_STR,IBEAM,NPART,
     &  NEVENT,NORDER,N_SECOND_PART
      Common / TI / nTI
      Character*40 Boson
      Common / Boson / Boson

      Logical Testing
      Data Testing /.False./

cdump
C      Testing = .True.

      If(Testing) then
        Print*,' L,X1A =', L,X1A
        Print*,' M,X2A =', M,X2A
        Print*,' N,X3A =', N,X3A
      endif

      If (.False.) then
        Print*, ' Input y grid to RatInt: '
        Print*, X3A
        Print*, ' Input Q grid to middle interpolation: '
        Print*, X2A
        Print*, ' Input QT grid to last interpolation: '
        Print*, X1A
        Print*, ' Input Xsect values to RatInt: '
        Print*, YA
        Print*, ' Rest of input to Polin3 (...,QT,Q,y): '
        Print*, L,M,N,X1,X2,X3
        Stop
      End If

      DO 12 i1=1,L
        DO 11 i2=1,M
          Do 10 i3 = 1,N
            Y3TMP(i3)=YA(i1,i2,i3)
C            Print*, i1,i2,i3, YA(i1,i2,i3)
10        Continue
C         Print*, ' Input "y" values to RatInt: '
C         Print*, Y3TMP(1), Y3TMP(2), Y3TMP(3), Y3TMP(4)
C         Interpolate in y direction using RatInt
C         Note: Ratint seem to perform better than PolInt,
C         at the edges. It 'shoots over' less.
C          Call RatInt(X3A,Y3TMP,N,X3,Y2TMP(i2),DY)
C         Ratint fails for the pion grids.
C         This "If" fix is to deal with the flat interpolation case
C         I'm not sure it improves PolInt.  It could be moved into Pol/RatInt
!          If (Y3TMP(1).Eq.Y3TMP(2) .and. Y3TMP(3).Eq.Y3TMP(4) .and.
!     .        Y3TMP(1).Eq.Y3TMP(3)) then
!            Y2TMP(i2) = Y3TMP(2)
!            GoTo 11
!          Else
            Call PolInt(X3A,Y3TMP,N,X3,Y2TMP(i2),DY)
!          End If
CsB       If PolInt fails, try RatInt.
          If (Y2TMP(i2).lt.0. .and. Y3TMP(2).GT.0. .and. 
     .        Y3TMP(3).GT.0.) then
            If(iBeam.NE.-2 .and. iBeam.NE.0) then
!           Diagnosis:
!           Print*, ' Q_T = ', X1A(i1),' Q = ', X2A(i2)
!           Print*, '-------------'
!           Print*, ' y = ', X3A
!           Print*, '-------------'
!           Print*, Y3TMP
!           Print*, '-------------'
!           Print*, X3,Y2TMP(i2)
!           Stop
              Call RatInt(X3A,Y3TMP,N,X3,Y2TMP(i2),DY)
            Else
c             If ( Y3TMP(1).GT.0. .and. Y3TMP(4).LE.0. ) Y2TMP(i2)=0.d0
            End If
          End If
CsB       Testing the y interpolation
          If (Testing .and. Y2TMP(i2).lt.0.0) then
            Print*, ' Input  of y interpolation: ', X3A, Y3TMP
            Print*, ' Result of y interpolation: ', X3,X2,X1, Y2TMP(i2)
            Print*, ' Stopped in y interpolation '
            Stop
          End If ! End Testing y interpolation
11      CONTINUE
cc        Print*, L, M, N, X2A(1), X2A(M), Y2TMP(1), Y2TMP(M)
C       This "If" fix is to deal with the flat interpolation case
C       I'm not sure it improves PolInt.  It could be moved into Pol/RatInt
!        If (Y2TMP(1).Eq.Y2TMP(2) .and. Y2TMP(3).Eq.Y2TMP(4) .and.
!     .      Y2TMP(1).Eq.Y2TMP(3)) then
!          Y1TMP(i1) = Y2TMP(2)
!          GoTo 12
!        Else
          CALL POLINT(X2A,Y2TMP,M,X2,Y1TMP(i1),DY)
!        End If
CsB     Testing the middle interpolation
        If (Testing .and. Y2TMP(i1).lt.0.0) then
          Print*, ' Input  of mid interp.: ', X2A, Y2TMP
          Print*, ' Result of mid interp.: ', X3, X2, X1, Y1TMP(i1)
          Print*, ' Stopped in mid interpolation '
          Stop
        End If ! End Testing middle interpolation
12    CONTINUE
c     Print*, X1A(1), X1A(L), Y1TMP(1), Y1TMP(L)
C
CsB We concluded that the qT interpolation for the low qT portion of the
C   Asymptotic grid fails in lack of many gridpoints.
C   On the other hand we know the simple qT shape for the Asymptotic case
C   to be 1/qT (a[Q,y] Log[Q^2/qT^2] + b[Q,y]).
C   Here we use this to determine a and b and calculate the cross section
C   analytically instead of interpolating.
C
C      If (.False.) then ! asymptotic case, low pT part, do not call it yet:
c         have to figure out what to do for high (~20) qT's.
      If (W.GE.900 .and. nPart.eq.12 .and. nTI.eq.0 .and.
     .    X1.LT.5.d0 .and. Y1Tmp(3).GT.1.d0) then
CsB Asymptotic case, low qT part, small qT, large enough cross section
C   If the warning ' Low precision in the analytical qT shape.'
C   persists just turn the analytic formula off.
C   (AnalAsym needs to know Q (X2) to give the cross section at qT.
C    The y value (X3) just passed to check things whithin.)
        Call AnalAsym(X1A,Y1TMP,L,X1,X2,X3,Y,DY)
      Else ! resum, delta sigma, LO cases
C       This "If" fix is to deal with the flat interpolation case
C       I'm not sure it improves PolInt.  It could be moved into Pol/RatInt
!        If (Y1TMP(1).Eq.Y1TMP(2) .and. Y1TMP(3).Eq.Y1TMP(4) .and.
!     .      Y1TMP(1).Eq.Y1TMP(3)) then
!          Y = Y1TMP(2)
!        Else
          CALL POLINT(X1A,Y1TMP,L,X1,Y,DY)
!        End If
      End If

      RETURN
      END

CsB -----------------------------------------------------------------------
      Subroutine Polin3m(X1A,X2A,X3A,YA,L,M,N,X1,X2,X3,Y,DY)
CsB -----------------------------------------------------------------------
CsB * Modified from Polin2 of Numerical recepies by C. Balazs at 1/95.
CsB *
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (LMAX=20, MMAX=20, NMAX=20)
      DIMENSION X1A(L), X2A(M), X3A(N), YA(L,M,N),
     >          Y1TMP(LMAX), Y2TMP(MMAX), Y3TMP(NMAX)
      Logical Testing
      Common / ECMC / ECM

      Testing = .True.
      Testing = .False.
            
CsB * Whether we are at the left or the right part of the y spectrum?
c      iSgn = Sign(1,X3)
      iSgn = Int(X3/Abs(X3))
c      If (Testing) Print*, X3, iSgn

      DO 12 i1=1,L
        DO 11 i2=1,M
          Do 10 i3 = 1,N
            Y3TMP(i3)=YA(i1,i2,i3)
CsB            If (Testing) Print*, i1,i2,i3, YA(i1,i2,i3)
10        Continue
c          Print*, X3A(1), X3A(2), X3A(3), X3A(4),
c     Print*, Y3TMP(1), Y3TMP(2), Y3TMP(3), Y3TMP(4)
CsB * Take care of yMax
          Call yMaximum (ECM, X1A(i1), X2A(i2), yMax)
          Call Locate (Y3TMP, N, iSgn*yMax, LyM)
          If (LyM.GT.0 .And. LyM.LT.N) Y3TMP(LyM+(iSgn+1)/2) = yMax
c
c          If ( (-yMax.GT.Y3TMP(1) .And. -yMax.LT.Y3TMP(N)) .Or.
c     >         (-yMax.GT.Y3TMP(1) .And.  yMax.LT.Y3TMP(N))  )   then
c
c          EndIf
CsB * Interpolate in y direction
C         This "If" fix is to deal with the flat interpolation case
C         I'm not sure it improves PolInt.  It could be moved into Pol/RatInt
!         If (Y3TMP(1).Eq.Y3TMP(2) .and. Y3TMP(3).Eq.Y3TMP(4) .and.
!    .        Y3TMP(1).Eq.Y3TMP(3)) then
!           Y2TMP(i2) = Y3TMP(2)
!           GoTo 11
!         Else
            Call Polint(X3A,Y3TMP,N,X3,Y2TMP(i2),DY)
!         End If
c          Print*, Y2TMP(i2)
11      CONTINUE
c       Print*, L, M, N, X2A(1), X2A(M), Y2TMP(1), Y2TMP(M)
C       This "If" fix is to deal with the flat interpolation case
C       I'm not sure it improves PolInt.  It could be moved into Pol/RatInt
!       If (Y2TMP(1).Eq.Y2TMP(2) .and. Y2TMP(3).Eq.Y2TMP(4) .and.
!    .      Y2TMP(1).Eq.Y2TMP(3)) then
!         Y1TMP(i1) = Y2TMP(2)
!         GoTo 12
!       Else
          CALL POLINT(X2A,Y2TMP,M,X2,Y1TMP(i1),DY)
!       End If
12    CONTINUE
c     Print*, X1A(1), X1A(L), Y1TMP(1), Y1TMP(L)
C     This "If" fix is to deal with the flat interpolation case
C     I'm not sure it improves PolInt.  It could be moved into Pol/RatInt
!     If (Y1TMP(1).Eq.Y1TMP(2) .and. Y1TMP(3).Eq.Y1TMP(4) .and.
!    .    Y1TMP(1).Eq.Y1TMP(3)) then
!       Y = Y1TMP(2)
!     Else
        CALL POLINT(X1A,Y1TMP,L,X1,Y,DY)
!     End If
      RETURN
      END

C -------------------------------------------------------------------------
      SUBROUTINE POLINT(XA,YA,N,X,Y,DY)
C -------------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      PARAMETER (NMAX=10)
      DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)

      If (.True.) then
!       Print*, XA
!       Print*, YA
!       Print*, N,X,Y,DY
!       Stop
      EndIf

      NS=1
      DIF=ABS(X-XA(1))
      DO 11 I=1,N
        DIFT=ABS(X-XA(I))
        IF (DIFT.LT.DIF) THEN
          NS=I
          DIF=DIFT
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
          IF(DEN.EQ.0.) Print*, ' DEN = 0 in Polint. '
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE
      RETURN
      END

CsB *** RatInt from Numerical Recepies ******************************
C -------------------------------------------------------------------------
      SUBROUTINE RATINT(XA,YA,N,X,Y,DY)
C -------------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (NMAX=10,TINY=1.E-25)
      Common /RatIntErr/ iDD0
      DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)
CCAO  
!     INTEGER*4 IFLAG
!     COMMON/RET/ IFLAG
!     IFLAG=0
      iDD0 = 0
      NS=1
      HH=ABS(X-XA(1))
      DO 11 I=1,N
        H=ABS(X-XA(I))
        IF (H.EQ.0.)THEN
          Y=YA(I)
          DY=0.0
          RETURN
        ELSE IF (H.LT.HH) THEN
          NS=I
          HH=H
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)+TINY
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          W=C(I+1)-D(I)
          H=XA(I+M)-X
          T=(XA(I)-X)*D(I)/H
          DD=T-C(I+1)
CCAO      IF DUE TO NUMERICAL ACCUARACY, DD IS INDEED EQUAL TO 0
C         THEN SET IFLAG=-997 AND REJECT THIS EVENT
C         IFLAG IS PASSED THROUGH COMMON BLOCK /RET/
!         IF(DD.EQ.0) THEN 
!           IFLAG=-997
!           GOTO 15
!         ENDIF
          IF(DD.EQ.0.) then
            Print*, ' >>> DD = 0 in RatInt. <<< '
CsB         I think this happens when an identically 0 input vector is interpolated.
            iDD0 = 1
          ENDIF
          DD=W/DD
          D(I)=C(I+1)*DD
          C(I)=T*DD
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE
CCAO  ----------------
!15   CONTINUE
      RETURN
      END ! RatInt

C---------------------------------------------------------------------------
      Subroutine AnalAsym(X1A,Y1TMP,L,X1,X2,X3,Y,DY)
C---------------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
CsB X1A() contains the qT values
C   Y1TMP() contains the cross section values
C   L is the dim of X1A (& Y1TMP ?)
C   X1 is the qT point where we are looking for the cross section
C   X2 is Q
C   The interpolated value is returned in Y
C   Te error of interpolation is returned in DY
      DIMENSION X1A(L), Y1TMP(L)
      Common / PrintLevel / lPr
      Logical Testing
      Data Testing /.False./

C      Testing = .True.
      If (Testing) then
CsB Check the arguments
        Print*, 'X1A,Y1TMP,L,X1,X2,Y,DY'
        Print*,  X1A,Y1TMP,L,X1,X2,Y,DY
      End If

CsB As reference points we choose the first and the third one.
C   This is because we want to avoid the last, which may be small.
C   If the last point is < ~1, then the a and b is calculated with
C   very little precision.
      qT1 = X1A(1)
      qT2 = X1A(3)
      DsDqT1 = Y1TMP(1)
      DsDqT2 = Y1TMP(3)

CsB Solution of the generic pT shape
      Q = X2
      a = (-DsDqT1*qT1 + DsDqT2*qT2)/Log(qT1**2/qT2**2)
      b = (DsDqT2*qT2*Log(Q**2/qT1**2) - DsDqT1*qT1*Log(Q**2/qT2**2))/
     /     Log(qT2**2/qT1**2)

CsB First test the formula on a known grid point, on point #2.
      qT = X1A(2)
      Y  = 1/qT*(a*Log(Q**2/qT**2) + b)
      DY = (Y-Y1TMP(2))/Y1TMP(2)
CsB If the analytic results differs from the grid point more than 0.01%
C   return zero.
        ii = ii + 1
      If (Abs(DY).GT.1e-2 .and. lPr.GE.2) then
        Write(14,*) ' Low precision in the analytical qT shape. '
        Write(14,*) ' Inputs: X1A,Y1TMP,L,ii,X1,X2,X3,Y,DY '
        Write(14,*)   X1A,Y1TMP,L,ii,X1,X2,X3,Y,DY
        Write(14,*) ' qT, Anal. value, error to grid :'
        Write(14,*)   qT, Y, DY
        Return
C        Stop
      End If

CsB If the formula works then calculate it for the qT point needed.
      qT = X1
      Y  = 1/qT*(a*Log(Q**2/qT**2) + b)
C      Print*, ' "Interpolated" value: ', Y
C      Stop

      Return
      End

CsB *** MISC.FOR ******************************************************

C---------------------------------------------------------------------------
      LOGICAL FUNCTION ASK(QUERY)
C---------------------------------------------------------------------------
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
C
      CHARACTER QUEND*11,  QUERY*(*), CH*1
C
      PARAMETER (QUEND= ' (Y OR N)? ' )
C
      CALL RTB (QUERY, LEN)
1     WRITE(NOUT, 90) QUERY(1:LEN), QUEND
      READ(NIN, 91) CH
      CALL UPCASE (CH)
      ASK = .FALSE.
      IF ( (CH.EQ.'Y') .OR. (CH.EQ.' ') ) THEN
         ASK = .TRUE.
         RETURN
      ELSE IF (CH.EQ.'N') THEN
         RETURN
      ELSE
         WRITE(NOUT, 92)
      ENDIF
      GOTO 1
 90     FORMAT ('$', 3A)
 91     FORMAT (A1)
 92     FORMAT (' BAD ANSWER--TRY AGAIN ')
C               *************************
      END

      SUBROUTINE RTB (CH, LENCH)
C                        Set LENCH = length of CH, not counting trailing
C                        blanks and nulls
      CHARACTER*(*) CH
C
      CHARACTER*1 NEXT, NULL
C
      NULL = CHAR(0)
      MAX = LEN(CH)
      DO 1  I = 1, MAX
            LENCH = MAX + 1 - I
            NEXT = CH(LENCH:LENCH)
            IF ( (NEXT .NE. ' ') .AND. (NEXT .NE. NULL) )  RETURN
 1    CONTINUE
      LENCH = 0
      RETURN
C               *************************
      END

        SUBROUTINE UC(A)
C            Converts A to all upper case.
C            System dependent: assumes ICHAR(uc letter) - ICHAR(lc letter)
C            is constant.
        CHARACTER A*(*), C*(1)
        INTEGER I
        ENTRY UPCASE(A)
        DO 1 I=1, LEN(A)
                C = A(I:I)
C
C                       Use ASCII ordering for detecting lc:
                IF ( LGE(C, 'a') .AND. LLE(C, 'z') )THEN
                        A(I:I) = CHAR(ICHAR(C)+ICHAR('A')-ICHAR('a'))
                        ENDIF
 1              CONTINUE
        RETURN
C               *************************
        END
C

      SUBROUTINE WARNI (IWRN, NWRT, MSG, NMVAR, IVAB,
     >                  IMIN, IMAX, IACT)
C                                                  Integer version of above
      CHARACTER*(*) MSG, NMVAR

      IW = IWRN
      IV = IVAB

      WRITE (NWRT,'(1X,A/1X,2A, I10)') MSG, NMVAR, ' = ', IV

      IF  (IW .EQ. 0) THEN
         PRINT '(1X,A/1X, 2A,I10 /A,I4)', MSG, NMVAR, ' = ', IV,
     >         ' For all warning messages, check file unit #', NWRT
         IF (IACT .EQ. 1) THEN
         PRINT       '(A/2I10)', ' The limits are: ', IMIN, IMAX
         WRITE (NWRT,'(A/2I10)') ' The limits are: ', IMIN, IMAX
         ENDIF
      ENDIF

      IWRN = IW + 1

      RETURN
C               *************************
      END

      SUBROUTINE WARNR (IWRN, NWRT, MSG, NMVAR, VARIAB,
     >                  VMIN, VMAX, IACT)

C       Subroutine to handle warning messages.  Writes the (warning) message
C       and prints out the name and value of an offending variable to SYS$OUT
C       the first time, and to output file unit # NWRT in subsequent times.
C
C       The switch IACT decides whether the limits (VMIN, VMAX) are active or
C       not.

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)

      CHARACTER*(*) MSG, NMVAR

      IW = IWRN
      VR = VARIAB

      WRITE (NWRT,'(I5, 2A/1X,2A,1PD16.7)') IW, '   ', MSG,
     >                  NMVAR, ' = ', VR

      IF  (IW .EQ. 0) THEN
         PRINT '(1X, A/1X,2A,1PD16.7/A,I4)', MSG, NMVAR, ' = ', VR,
     >         ' For all warning messages, check file unit #', NWRT
         IF (IACT .EQ. 1) THEN
         PRINT       '(A/2(1PE15.4))', ' The limits are: ', VMIN, VMAX
         WRITE (NWRT,'(A/2(1PE15.4))') ' The limits are: ', VMIN, VMAX
         ENDIF
      ENDIF

      IWRN = IW + 1

      RETURN
      END
C                       ****************************

CsB * GBook *************************************************************

C YUAN : I MODIFY THIS SO THAT THE MAXIMUM BIN NUMBER IS 60 INSTEAD OF
C        100. THE CHANGE IS IN THE WRITE_FORMAT ::
C        WHERE A56 IS SUPPOSED TO BE A108 FOR WIDTH=132 SCREEN,
C        IN OTHER WORDS, FOR 100 BINS.
C       REMEMBER TO SET ALL THE 'CONTENT VARIABLES' OF GBOOK
C       TO BE   REAL*4   IN THE DOUBLE PRECISION PROGRAM.
C*********************************************************************
C***  HEADER  ********************************************************
C
C  GBOOK - A SMALL HISTOGRAM PACKAGE
C  WRITTEN DECEMBER 1979, LAST CHANGED MARCH 1983
C  AUTHOR: TORBJORN SJOSTRAND, DEPT. OF THEORETICAL PHYSICS,
C  UNIVERSITY OF LUND, SOLVEGATAN 14A, S-223 62 LUND, SWEDEN
C  PRESENT ADDRESS: FNAL, THEORY, TEL. 312-840-3753
C  PLEASE REPORT ANY ERRORS TO THE AUTHOR
C
C***  INTRODUCTION  **************************************************
C
C  GBOOK IS A SMALL SUBROUTINE PACKAGE FOR GETTING ONE- OR TWO-
C  DIMENSIONAL HISTOGRAMS ON AN ORDINARY LINE PRINTER OR TERMINAL.
C  IT CAN BE USED IN A MANNER VERY SIMILAR TO HBOOK, BUT HAS BEEN
C  WRITTEN COMPLETELY INDEPENDENTLY. THE CODE IS FULLY IN FORTRAN77,
C  WITH A MINIMUM OF MACHINE DEPENDENCE. THE USAGE CAN BE DIVIDED INTO
C  FOUR STEPS: BOOKING, FILLING, EDITING AND PRINTING. ALL SUBROUTINES
C  HAVE NAMES SIX CHARACTERS LONG AND BEGINNING WITH G.
C
C***  BOOKING  *******************************************************
C
C  THERE ARE 100 HISTOGRAMS AT THE DISPOSAL OF THE USER, EACH GIVEN BY
C  A NUMBER BETWEEN 1 AND 100. BEFORE A HISTOGRAM CAN BE USED, SPACE
C  MUST BE RESERVED FOR IT.
C
C  CALL GBOOK1(ID,TITLE,NX,XL,XU)
C  PURPOSE: BOOK A ONE-DIMENSIONAL HISTOGRAM.
C  ID: HISTOGRAM NUMBER, INTEGER BETWEEN 1 AND 100.
C  TITLE: HISTOGRAM TITLE, CAN BE GIVEN EITHER AS A CHARACTER STRING
C     OF AT MOST 60 CHARACTERS OR AS A CHARACTER*60 VARIABLE.
C  NX: NUMBER OF BINS (IN X DIRECTION) IN HISTOGRAM, INTEGER BETWEEN
C     1 AND 100.
C  XL, XU: LOWER AND UPPER BOUND, RESPECTIVELY, ON THE (X) RANGE
C     COVERED BY THE HISTOGRAM.
C
C  CALL GBOOK2(ID,TITLE,NX,XL,XU,NY,YL,YU)
C  PURPOSE: BOOK A TWO-DIMENSIONAL HISTOGRAM.
C  ID: HISTOGRAM NUMBER, INTEGER BETWEEN 1 AND 100.
C  TITLE: HISTOGRAM TITLE, SEE GBOOK1.
C  NX: NUMBER OF BINS IN X DIRECTION, INTEGER BETWEEN 1 AND 50.
C  XL, XU: LOWER AND UPPER BOUND ON X RANGE OF HISTOGRAM.
C  NY: NUMBER OF BINS IN Y DIRECTION, ARBITRARY POSITIVE INTEGER.
C  YL, YU: LOWER AND UPPER BOUND OF Y RANGE OF HISTOGRAM.
C
C***  FILLING  *******************************************************
C
C  FOR BOOKED HISTOGRAMS WEIGHTS CAN BE FILLED AT GIVEN COORDINATES.
C
C  CALL GFILL1(ID,X,W)
C  PURPOSE: FILL IN A ONE-DIMENSIONAL HISTOGRAM.
C  ID: HISTOGRAM NUMBER.
C  X: X COORDINATE OF POINT.
C  W: WEIGHT TO BE ADDED IN THIS POINT.
C
C  CALL GFILL2(ID,X,Y,W)
C  PURPOSE: FILL IN A TWO-DIMENSIONAL HISTOGRAM.
C  ID: HISTOGRAM NUMBER.
C  X: X COORDINATE OF POINT.
C  Y: Y COORDINATE OF POINT.
C  W: WEIGHT TO BE ADDED IN THIS POINT.
C
C***  EDITING  *******************************************************
C
C  FOR EDITING OF HISTOGRAMS BEFORE PRINTOUT TWO ROUTINES ARE AVAILABLE.
C
C  CALL GSCALE(ID,F)
C  PURPOSE: RESCALE THE CONTENTS OF A HISTOGRAM.
C  ID: HISTOGRAM NUMBER.
C  F: RESCALING FACTOR, I.E. FACTOR THAT ALL BIN CONTENTS (INCLUDING
C     OVERFLOW ETC.) ARE MULTIPLIED BY.
C  REMARK: A TYPICAL RESCALING FACTOR FOR A ONE-DIMENSIONAL HISTOGRAM
C     COULD BE F = 1/(BIN SIZE * NUMBER OF EVENTS) =
C     = NX/(XU-XL) * 1/(NUMBER OF EVENTS).
C
C  CALL GOPERA(ID1,OPER,ID2,ID3,F1,F2)
C  PURPOSE: THIS IS A GENERAL PURPOSE ROUTINE FOR EDITING ONE OR SEVERAL
C     HISTOGRAMS, WHICH ALL ARE ASSUMED TO HAVE THE SAME NUMBER OF
C     BINS. OPERATIONS ARE CARRIED OUT BIN BY BIN, INCLUDING OVERFLOW
C     BINS ETC.
C  OPER: GIVES THE TYPE OF OPERATION TO BE CARRIED OUT, A ONE-CHARACTER
C     STRING OR A CHARACTER*1 VARIABLE.
C  OPER= '+', '-', '*', '/': ADD, SUBRACT, MULTIPLY OR DIVIDE THE
C     CONTENTS IN ID1 AND ID2 AND PUT THE RESULT IN ID3. F1 AND F2, IF
C     NOT 1., GIVE FACTORS BY WHICH THE ID1 AND ID2 BIN CONTENTS ARE
C     MULTIPLIED BEFORE THE INDICATED OPERATION. (DIVISION WITH A
C     VANISHING BIN CONTENT WILL GIVE 0.)
C  OPER= 'A', 'S', 'L': FOR 'S' THE SQUARE ROOT OF THE CONTENT IN ID1
C     IS TAKEN (RESULT 0 FOR NEGATIVE BIN CONTENTS) AND FOR 'L' THE
C     10-LOGARITHM IS TAKEN (A NONPOSITIVE BIN CONTENT IS BEFORE THAT
C     REPLACED BY 0.8 TIMES THE SMALLEST POSITIVE BIN CONTENT).
C     THEREAFTER, IN ALL THREE CASES, THE CONTENT IS MULTIPLIED BY F1
C     AND ADDED WITH F2, AND THE RESULT IS PLACED IN ID3. THUS ID2
C     IS DUMMY IN THESE CASES.
C  OPER= 'M': INTENDED FOR STATISTICAL ANALYSIS, BIN-BY-BIN MEAN AND
C     STANDARD DEVIATION OF A VARIABLE, ASSUMING THAT ID1 CONTAINS
C     ACCUMULATED WEIGHTS, ID2 ACCUMULATED WEIGHT*VARIABLE AND
C     ID3 ACCUMULATED WEIGHT*VARIABLE-SQUARED. AFTERWARDS ID2 WILL
C     CONTAIN THE MEAN VALUES (=ID2/ID1) AND ID3 THE STANDARD
C     DEVIATIONS (=SQRT(ID3/ID1-(ID2/ID1)**2)). IN THE END, F1
C     MULTIPLIES ID1 (FOR NORMALIZATION PURPOSES), WHILE F2 IS DUMMY.
C
C***  PRINTING  ******************************************************
C
C  AT PRINTING AXES ARE CHOSEN SUCH THAT HISTOGRAMS ARE SUPPOSED
C  TO FIT INTO ONE PAGE. FOR ONE-DIMENSIONAL HISTOGRAMS NUMBERS
C  SMALLER IN MAGNITUDE THAN 10**(-10) ARE CONSIDERED TO BE 0,
C  OTHERWISE SCALES ARE CHOSEN AUTOMATICALLY FOR MAXIMUM RESOLUTION.
C  TWO-DIMENSIONAL HISTOGRAMS ARE STRONGLY ORIENTED TOWARDS HAVING
C  WEIGHTS IN THE ORDER OF OR BIGGER THAN UNITY : SIGNS 1 - 9 ARE
C  CHOSEN IN A LINEAR SCALE UP TO WEIGHTS 9.5, WHEREAS A - Z WILL
C  EITHER BE CHOSEN IN A LINEAR SCALE WITH STEP 1 OR, IF MAXIMUM
C  WEIGHT IS LARGER THAN 36.5, IN A LOGARITHMICALLY EVEN SCALE.
C  NEGATIVE NUMBERS ARE ALLOWED BOTH IN ONE- AND TWO-DIMENSIONAL
C  HISTOGRAMS, AND ARE PRECEDED BY A - SIGN.
C
C  CALL GCLEAR
C  PURPOSE: PRINT OUT ALL HISTOGRAMS THAT HAVE BEEN FILLED, AND
C     RESET THEM THEREAFTER TO 0.
C
C  CALL GPRINT(ID)
C  PURPOSE: PRINT OUT A SINGLE HISTOGRAM.
C  ID: HISTOGRAM TO BE PRINTED.
C
C  CALL GRESET(ID)
C  PURPOSE: RESET ALL BIN CONTENTS, INCLUDING OVERFLOW ETC., TO 0.
C  ID: HISTOGRAM TO BE RESET.
C
C***  COMMON BLOCK AND SPACE REQUIREMENTS  ***************************
C
C  A COMMONBLOCK
C  COMMON /GBOOK/ A(10000)
C  IS USED TO STORE HISTOGRAM INFORMATION. THE HISTOGRAM INDEX TAKES
C  102 POSITIONS. EACH BOOKED ONE- (TWO-) DIMENSIONAL HISTOGRAM TAKES
C  AN ADDITIONAL 38+NX (38+NX*NY)  POSITIONS. THE PROGRAM HAS TO BE
C  RECOMPILED WITH CHANGED COMMONBLOCK IF MORE IS REQUIRED.
C
C***  END  DESCRIPTION  **********************************************
C*********************************************************************

      SUBROUTINE GBOOK1(ID,TITLE,NX,XL,XU)
      COMMON /GBOOK/ A(10000)
      CHARACTER TITLE*(*),TITFX*60
      EQUIVALENCE (REQ,IEQ)
      A(ID+2)=A(1)+A(2)
      A(2)=A(2)+38+NX
      IS=A(ID+2)+0.5
      A(IS+1)=NX
      A(IS+2)=XL
      A(IS+3)=XU
      A(IS+4)=(XU-XL)/NX
      A(IS+5)=1
      CALL GRESET(ID)
      TITFX=TITLE//' '
      DO 100 IT=1,20
      IEQ=256**2*ICHAR(TITFX(3*IT-2:3*IT-2))+256*ICHAR(TITFX(3*IT-1:
     &3*IT-1))+ICHAR(TITFX(3*IT:3*IT))
      A(IS+18+NX+IT)=REQ
  100  continue
      RETURN
      END

C*********************************************************************

      SUBROUTINE GBOOK2(ID,TITLE,NX,XL,XU,NY,YL,YU)
      COMMON /GBOOK/ A(10000)
      CHARACTER TITLE*(*),TITFX*60
      EQUIVALENCE (REQ,IEQ)
      A(ID+2)=A(1)+A(2)
      A(2)=A(2)+38+NX*NY
      IS=A(ID+2)+0.5
      A(IS+1)=NX
      A(IS+2)=XL
      A(IS+3)=XU
      A(IS+4)=(XU-XL)/NX
      A(IS+5)=NY
      A(IS+6)=YL
      A(IS+7)=YU
      A(IS+8)=(YU-YL)/NY
      CALL GRESET(ID)
      TITFX=TITLE//' '
      DO 100 IT=1,20
      IEQ=256**2*ICHAR(TITFX(3*IT-2:3*IT-2))+256*ICHAR(TITFX(3*IT-1:
     &3*IT-1))+ICHAR(TITFX(3*IT:3*IT))
  100 A(IS+18+NX*NY+IT)=REQ
      RETURN
      END

C*********************************************************************

      SUBROUTINE GFILL1(ID,X,W)
      COMMON /GBOOK/ A(10000)
      IS=A(ID+2)+0.5
      A(IS+9)=A(IS+9)+1.
      IOX=2
      IF(X.LT.A(IS+2)) IOX=1
      IF(X.GE.A(IS+3)) IOX=3
      A(IS+12+IOX)=A(IS+12+IOX)+W
      IF(IOX.NE.2) RETURN
cpn2006 Prevent the case when ix=Nbins because of the roundoff error,
c       corrupting the title of the next histogram
c      IX=(X-A(IS+2))/A(IS+4)
      IX=min( (X-A(IS+2))/A(IS+4), a(is+1)-1)
      A(IS+19+IX)=A(IS+19+IX)+W
      RETURN
      END

C*********************************************************************

      SUBROUTINE GFILL2(ID,X,Y,W)
      COMMON /GBOOK/ A(10000)
      IS=A(ID+2)+0.5
      A(IS+9)=A(IS+9)+1.
      IOX=2
      IF(X.LT.A(IS+2)) IOX=1
      IF(X.GE.A(IS+3)) IOX=3
      IOY=2
      IF(Y.LT.A(IS+6)) IOY=1
      IF(Y.GE.A(IS+7)) IOY=3
      A(IS+6+3*IOY+IOX)=A(IS+6+3*IOY+IOX)+W
      IF(IOX.NE.2.OR.IOY.NE.2) RETURN
      IX=(X-A(IS+2))/A(IS+4)
      IY=(Y-A(IS+6))/A(IS+8)
      IC=INT(A(IS+1)+0.5)*IY+IX
      A(IS+19+IC)=A(IS+19+IC)+W
      RETURN
      END

C*********************************************************************

      SUBROUTINE GSCALE(ID,F)
      COMMON /GBOOK/ A(10000)
      IS=A(ID+2)+0.5
      DO 100 IC=IS+10,IS+18+INT(A(IS+1)+0.5)*INT(A(IS+5)+0.5)
  100 A(IC)=F*A(IC)
      RETURN
      END

C*********************************************************************

      SUBROUTINE GOPERA(ID1,OPER,ID2,ID3,F1,F2)
      COMMON /GBOOK/ A(10000)
      CHARACTER OPER*1
      IS1=A(ID1+2)+0.5
      IS2=A(ID2+2)+0.5
      IS3=A(ID3+2)+0.5
      NC=INT(A(IS3+1)+0.5)*INT(A(IS3+5)+0.5)
      IF(OPER.EQ.'+'.OR.OPER.EQ.'-'.OR.OPER.EQ.'*'.OR.OPER.EQ.'/')
     &A(IS3+9)=A(IS1+9)+A(IS2+9)
      IF(OPER.EQ.'A'.OR.OPER.EQ.'S'.OR.OPER.EQ.'L') A(IS3+9)=A(IS1+9)
      IF(OPER.EQ.'+') THEN
      DO 100 IC=10,18+NC
  100 A(IS3+IC)=F1*A(IS1+IC)+F2*A(IS2+IC)
      ELSEIF(OPER.EQ.'-') THEN
      DO 110 IC=10,18+NC
  110 A(IS3+IC)=F1*A(IS1+IC)-F2*A(IS2+IC)
      ELSEIF(OPER.EQ.'*') THEN
      DO 120 IC=10,18+NC
  120 A(IS3+IC)=F1*A(IS1+IC)*F2*A(IS2+IC)
      ELSEIF(OPER.EQ.'/') THEN
      DO 130 IC=10,18+NC
      FA2=F2*A(IS2+IC)
      IF(ABS(FA2).LE.1E-10) A(IS3+IC)=0.
  130 IF(ABS(FA2).GT.1E-10) A(IS3+IC)=F1*A(IS1+IC)/FA2
      ELSEIF(OPER.EQ.'A') THEN
      DO 140 IC=10,18+NC
  140 A(IS3+IC)=F1*A(IS1+IC)+F2
      ELSEIF(OPER.EQ.'S') THEN
      DO 150 IC=10,18+NC
  150 A(IS3+IC)=F1*SQRT(MAX(0.,A(IS1+IC)))+F2
      ELSEIF(OPER.EQ.'L') THEN
      ZMIN=1E30
      DO 160 IC=19,18+NC
  160 IF(A(IS1+IC).LT.ZMIN.AND.A(IS1+IC).GT.1E-20) ZMIN=0.8*A(IS1+IC)
      DO 170 IC=10,18+NC
  170 A(IS3+IC)=F1*ALOG10(MAX(A(IS1+IC),ZMIN))+F2
      ELSEIF(OPER.EQ.'M') THEN
      DO 180 IC=10,18+NC
      IF(ABS(A(IS1+IC)).LE.1E-10) A(IS2+IC)=0.
      IF(ABS(A(IS1+IC)).GT.1E-10) A(IS2+IC)=A(IS2+IC)/A(IS1+IC)
      IF(ID3.NE.0.AND.ABS(A(IS1+IC)).LE.1E-10) A(IS3+IC)=0.
      IF(ID3.NE.0.AND.ABS(A(IS1+IC)).GT.1E-10) A(IS3+IC)=
     &SQRT(MAX(A(IS3+IC)/A(IS1+IC)-A(IS2+IC)**2,0.))
  180 A(IS1+IC)=F1*A(IS1+IC)
      ENDIF
      RETURN
      END

C*********************************************************************

      SUBROUTINE GCLEAR
      COMMON /GBOOK/ A(10000)
      DO 100 ID=1,INT(A(1)+0.5)
      IS=A(ID+2)+0.5
      IF(IS.EQ.0.OR.A(IS+9).LT.0.5) GOTO 100
      CALL GPRINT(ID)
      CALL GRESET(ID)
  100 CONTINUE
      RETURN
      END

C*********************************************************************

      SUBROUTINE GPRINT(ID)
      COMMON /GBOOK/ A(10000)
      CHARACTER TITLE*60, CTIME*8, OUT*56, CHA(40)*1
      DIMENSION IROW(100), IFRA(100), DYAC(10), EV(20)
      EQUIVALENCE (REQ,IEQ)
      DATA DYAC/.04,.05,.06,.08,.10,.12,.15,.20,.25,.30/, LIN/41/
      DATA CHA/' ','0','1','2','3','4','5','6','7','8','9','A','B',
     &'C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q',
     &'R','S','T','U','V','W','X','Y','Z','*','-','!'/
CCPY PRINT OUT TOPDRAW FORMAT:
      IO_TD=36
CCPY FOR TOPDRAW OUTPUT
C      I_TOPDRAW=0
      I_TOPDRAW=1

      IS=A(ID+2)+0.5
      IF(A(IS+9).LT.0.5) WRITE(16,400) ID
      IF(A(IS+9).LT.0.5) RETURN
      NX=INT(A(IS+1)+0.5)
      NY=INT(A(IS+5)+0.5)
      DO 100 IT=1,20
      EV(IT)=0.
      REQ=A(IS+18+NX*NY+IT)
  100 TITLE(3*IT-2:3*IT)=CHAR(IEQ/256**2)//CHAR(MOD(IEQ,256**2)/256)
     &//CHAR(MOD(IEQ,256))
C YUAN
C      CALL IDATE(IMON,IDAY,IYEAR)
C      CALL TIME(CTIME)
C      WRITE(16,410) ID, TITLE, IYEAR, IMON, IDAY, CTIME(1:5)
      WRITE(16,410) ID, TITLE
CCPY
      IF(I_TOPDRAW.EQ.1) THEN
        WRITE(IO_TD,410) ID, TITLE
      ENDIF
      IF(NY.EQ.1) THEN
      YMIN=A(IS+19)
      YMAX=A(IS+19)
      DO 120 IX=IS+20,IS+18+NX
      IF(A(IX).LT.YMIN) YMIN=A(IX)
  120 IF(A(IX).GT.YMAX) YMAX=A(IX)
      IF(YMAX-YMIN.GT.LIN*DYAC(1)*1E-9) THEN
      IF(YMIN.GT.0..AND.YMIN.LT.YMAX/10.) YMIN=0.
      IF(YMAX.LT.0..AND.YMAX.GT.YMIN/10.) YMAX=0.
      IPOT=INT(ALOG10(YMAX-YMIN)+10.)-10
      IF(YMAX-YMIN.LT.LIN*DYAC(1)*10.**IPOT) IPOT=IPOT-1
      IF(YMAX-YMIN.GT.LIN*DYAC(10)*10.**IPOT) IPOT=IPOT+1
      DELY=DYAC(1)
      DO 130 IDEL=1,9
  130 IF(YMAX-YMIN.GE.LIN*DYAC(IDEL)*10.**IPOT) DELY=DYAC(IDEL+1)
      DY=DELY*10.**IPOT
      DO 140 IX=1,NX
      CTA=ABS(A(IS+18+IX))/DY
      IROW(IX)=SIGN(CTA+0.95,A(IS+18+IX))
  140 IFRA(IX)=10.*(CTA+1.25-AINT(CTA+0.95))
      IRMI=SIGN(ABS(YMIN)/DY+0.95,YMIN)
      IRMA=SIGN(ABS(YMAX)/DY+0.95,YMAX)
      DO 160 IR=IRMA,IRMI,-1
      IF(IR.EQ.0) GOTO 160
      OUT=' '
      DO 150 IX=1,NX
      IF(IR.EQ.IROW(IX)) OUT(IX:IX)=CHA(IFRA(IX)+23*(IFRA(IX)/12))
  150 IF(IR*(IROW(IX)-IR).GT.0) OUT(IX:IX)=CHA(35)
      WRITE(16,420) IR*DELY, IPOT, OUT
  160 CONTINUE
      IPOT=INT(ALOG10(MAX(YMAX,-YMIN))+10.0001)-10
      OUT=' '
      DO 170 IX=1,NX
      IF(A(IS+18+IX).LT.-10.**(IPOT-4)) OUT(IX:IX)=CHA(39)
  170 IROW(IX)=10.**(3-IPOT)*ABS(A(IS+18+IX))+0.5
      WRITE(16,430) OUT
      DO 190 IR=4,1,-1
      DO 180 IX=1,NX
  180 OUT(IX:IX)=CHA(2+MOD(IROW(IX),10**IR)/10**(IR-1))
  190 WRITE(16,440) IPOT+IR-4, OUT
      IPOT=INT(ALOG10(MAX(-A(IS+2),A(IS+3)-A(IS+4)))+10.0001)-10
      OUT=' '
      DO 200 IX=1,NX
      IF(A(IS+2)+(IX-1)*A(IS+4).LT.-10**(IPOT-3)) OUT(IX:IX)=CHA(39)
  200 IROW(IX)=10.**(2-IPOT)*ABS(A(IS+2)+(IX-1)*A(IS+4))+0.5
      WRITE(16,450) OUT
      DO 220 IR=3,1,-1
      DO 210 IX=1,NX
  210 OUT(IX:IX)=CHA(2+MOD(IROW(IX),10**IR)/10**(IR-1))
  220 WRITE(16,440) IPOT+IR-3, OUT
      ENDIF
      IF(I_TOPDRAW.EQ.1) THEN
CsB -> 3
        Write(IO_TD,*) ' '
        Write(IO_TD,*) '    x (central value), sigma(x)/BinSize'
CsB      BinSize = A(IS+4)
      ENDIF
      DO 230 IX=1,NX
      CTA=ABS(A(IS+18+IX))
      X=A(IS+2)+(IX-0.5)*A(IS+4)
      IF(I_TOPDRAW.EQ.1) THEN
CsB -> 1
        Write(IO_TD,*) X, CTA/A(IS+4)
      ENDIF
      EV(1)=EV(1)+CTA
      EV(2)=EV(2)+CTA*X
  230 EV(3)=EV(3)+CTA*X**2
      XMEAN=EV(2)/MAX(EV(1),1E-20)
      XRMS=SQRT(MAX(0.,EV(3)/MAX(EV(1),1E-20)-XMEAN**2))
      WRITE(16,460) INT(A(IS+9)+0.5),XMEAN,A(IS+13),A(IS+2),A(IS+14),
     &XRMS,A(IS+15),A(IS+3)
      ELSE
      ZMAX=35.4
      DO 240 IC=IS+19,IS+NX*NY+18
  240 IF(ABS(A(IC)).GT.ZMAX) ZMAX=1.0001*ABS(A(IC))
      ZFAC=26./ALOG(ZMAX/9.5)
      DO 290 IY=1,MAX(NY+4,38)
      OUT=' '
      IF(IY.EQ.1.OR.IY.EQ.NY+2) THEN
      OUT(4:4)=CHA(38)
      DO 250 IX=5,2*NX+5
  250 OUT(IX:IX)=CHA(39)
      OUT(2*NX+6:2*NX+6)=CHA(38)
      ELSEIF(IY.EQ.NY+3) THEN
      DO 260 IX=6,2*NX+4,2
  260 OUT(IX:IX)=CHA(2+(IX-3)/20)
      ELSEIF(IY.EQ.NY+4) THEN
      DO 270 IX=6,2*NX+4,2
  270 OUT(IX:IX)=CHA(IX/2-10*((IX-3)/20))
      ELSEIF(IY.LE.NY+1) THEN
      OUT(1:1)=CHA(2+(NY+2-IY)/10)
      OUT(2:2)=CHA(4+NY-IY-10*((NY+2-IY)/10))
      OUT(4:4)=CHA(40)
      DO 280 IX=6,2*NX+4,2
      CT=A(IS+16+NX*(NY+1-IY)+IX/2)
      IF(CT.LE.-0.05) OUT(IX-1:IX-1)=CHA(39)
      CTA=ABS(CT)
      IF(CTA.GE.0.1) OUT(IX:IX)=CHA(2)
      IF(CTA.GE.0.5.AND.CTA.LT.36.) OUT(IX:IX)=CHA(2+INT(CTA+0.5))
  280 IF(ZMAX.GT.35.5.AND.CTA.GE.9.5) OUT(IX:IX)=
     &CHA(12+INT(ZFAC*ALOG(CTA/9.5)))
      OUT(2*NX+6:2*NX+6)=CHA(40)
      ENDIF
      IF(IY.GE.2.AND.IY.LE.38) THEN
      WMIN=0.1
      IF(IY.GE.3) WMIN=IY-2.5
      IF(ZMAX.GT.35.5.AND.IY.GE.12) WMIN=9.5*EXP((IY-12)/ZFAC)
      ENDIF
      IF(IY.EQ.1) WRITE(16,470) OUT
CCPY
C      IF(IY.GT.1.AND.IY.LE.38) WRITE(16,480)CHA(IY),WMIN,OUT
C  480 FORMAT(4X,A1,1X,F7.1,'-',3X,A56)

      IF(IY.GT.1.AND.IY.LE.38) WRITE(16,*)CHA(IY),WMIN,OUT

  290 IF(IY.GE.39) WRITE(16,490) OUT
      DO 300 IY=1,NY
      Y=A(IS+6)+(IY-0.5)*A(IS+8)
      DO 300 IX=1,NX
      X=A(IS+2)+(IX-0.5)*A(IS+4)
      CTA=ABS(A(IS+18+NX*(IY-1)+IX))
      EV(1)=EV(1)+CTA
      EV(2)=EV(2)+CTA*X
      EV(3)=EV(3)+CTA*X**2
      EV(4)=EV(4)+CTA*Y
      EV(5)=EV(5)+CTA*Y**2
  300 EV(6)=EV(6)+CTA*X*Y
      XMEAN=EV(2)/MAX(EV(1),1E-20)
      XRMS=SQRT(MAX(0.,EV(3)/MAX(EV(1),1E-20)-XMEAN**2))
      YMEAN=EV(4)/MAX(EV(1),1E-20)
      YRMS=SQRT(MAX(0.,EV(5)/MAX(EV(1),1E-20)-YMEAN**2))
      XYCOR=(EV(6)/MAX(EV(1),1E-20)-XMEAN*YMEAN)/MAX(1E-20,XRMS*YRMS)
      WRITE(16,500) INT(A(IS+9)+0.5),(A(IS+J), J=16,18),XMEAN,A(IS+2),
     &XRMS, A(IS+3), (A(IS+J), J=13,15), YMEAN
CBOX
        WRITE(16,502) A(IS+6), YRMS, A(IS+7),
     &(A(IS+J), J=10,12), XYCOR
C
      ENDIF
      RETURN
  400 FORMAT(/1X,'HISTOGRAM NO',I4,' : NO ENTRIES')
C  410 FORMAT('1'/5X,'HISTOGRAM NO',I4,4X,A60,'19',I2,'-',I2,'-',I2,
C     &1X,A5/)
  410 FORMAT(/1X,'HISTOGRAM NO',I4,4X,A56/)
  420 FORMAT(1X,F7.2,'*10**',I2,4X,A56)
  430 FORMAT(/7X,'CONTENTS',4X,A56)
  440 FORMAT(8X,'*10**',I2,4X,A56)
  450 FORMAT(/7X,'LOW EDGE',4X,A56)
  460 FORMAT(/1X,'ENTRIES=',I10,1P,1X,'MEAN=',E10.3,1X,'UNDERFL='
     &,E10.3,1X,'LOW EDGE =',E10.3/1X,'ALL CHAN=',E10.3,1X,
     &'RMS=',E10.3,1X,'OVERFL =',E10.3,1X,'HIGH EDGE=',E10.3)
CCPY
  470 FORMAT(1X,'SCALE',1X,A56)
  480 FORMAT(4X,A1,1X,F7.1,'-',3X,A56)
  490 FORMAT(17X,A56)
CBOX
  500 FORMAT(/1X,'ENTRIES=',I10,1P,1X,E10.3,' ? ',E10.3,' ? ',E10.3,
     &1X,'XMEAN =',E10.3/1X,'XMIN   =',E10.3,1X,36('-'),1X,'XRMS  =',
     &E10.3/1X,'XMAX   =',E10.3,1X,'OVFL',1X,E10.3,' ? ',E10.3,
     &' ? ',E10.3,1X,'YMEAN =',E10.3)
  502 FORMAT(1X,'YMIN   =',E10.3,1X,36('-'),
     &1X,'YRMS  =',E10.3/1X,'YMAX   =',E10.3,1X,E10.3,' ? ',E10.3,
     &' ? ',E10.3,1X,'XYCOR =',E10.3)
      END

C*********************************************************************

      SUBROUTINE GRESET(ID)
      COMMON /GBOOK/ A(10000)
      IS=A(ID+2)+0.5
      DO 110 IC=IS+9,IS+18+INT(A(IS+1)+0.5)*INT(A(IS+5)+0.5)
  110 A(IC)=0.
      RETURN
      END

C*********************************************************************

      BLOCK DATA GDATA
      COMMON /GBOOK/ A(10000)
      DATA (A(J), J=1,102)/100.,2.,100*0./
      END

C***  END GBOOK  *****************************************************

      SUBROUTINE ERROR(IPL)
      INTEGER LOOP(100)
      DATA MODE/0/
      IF(IPL.EQ.0) THEN
      MODE=0
      ELSEIF(IPL.LT.0) THEN
      MODE=1
      DO 100 J=1,100
  100 LOOP(J)=0
      NST=-IPL
      NEV=1
      ELSEIF(MODE.EQ.1.AND.NST*(NEV/NST).NE.NEV) THEN
      LOOP(IPL)=LOOP(IPL)+1
      NEV=NEV+1
      ELSEIF(MODE.EQ.1) THEN
      LOOP(IPL)=LOOP(IPL)+1
      WRITE(16,120) (J,J=0,10),(J1,(LOOP(J1+J2),J2=1,10),J1=0,90,10)
      DO 110 J=1,100
  110 LOOP(J)=0
      NEV=1
      ENDIF
      RETURN
  120 FORMAT(//20X,'ERROR DUMP OF LOOP COUNTS',11(/5X,11I6)/)
      END

CsB ******************************************************************

      SUBROUTINE EGHIST
      INTEGER*4 MXHIST,NHSHFT,ID
      LOGICAL HEXIST
      COMMON/MYHIST/MXHIST,NHSHFT,HEXIST(200)
      REAL*4 A
      COMMON /GBOOK/ A(10000)
      DO 100 ID=1,MXHIST
       IF(.NOT.HEXIST(ID)) GOTO 100
       CALL GPRINT(ID)
       CALL GRESET(ID)
100   CONTINUE
      RETURN
      END

CsB *******************************************************************

      DOUBLE PRECISION FUNCTION ALPHAS(AMU)
C     DEFINED IN THE MS_BAR SCHEME, TWO-LOOP ORDER
C     IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT NONE

      REAL*8 AMU
      REAL*8 W,QCDL,X1X2MIN
      INTEGER*4 NJET,N_PRIMARY_PART,NF,I_STR,IBEAM,NPART,NEVENT,NORDER,
     ,          N_SECOND_PART
      REAL*8 XMW,XMZ,XMW2,XMZ2,AEM,GAMW,GAMZ,GAMW2,GAMZ2,
     &SWS,CWS,XMBOT,XMTOP,XMB2,XMT2,GWEAK2,
     &XMCHARM,XMCHARM2,XMHIGGS,GAMH,GAMH2
      COMMON/AA/W,QCDL,X1X2MIN,NJET,N_PRIMARY_PART,NF,I_STR,IBEAM,NPART,
     &  NEVENT,NORDER,N_SECOND_PART
      COMMON/STAND0/XMW,XMZ,XMW2,XMZ2,AEM,GAMW,GAMZ,GAMW2,GAMZ2
      COMMON/STAND1/SWS,CWS,XMBOT,XMTOP,XMB2,XMT2,GWEAK2
      COMMON/STAND5/XMCHARM,XMCHARM2,XMHIGGS,GAMH,GAMH2
      CHARACTER*10 PDF_EVL_NAME
      COMMON/PDFFILE/ PDF_EVL_NAME
      Integer BIPDF,iSet
      Common / BuildInPDF / BIPDF,iSet

      INTEGER ISTART
      DATA ISTART/0/
      REAL*8 B0,B1,XRAT,RAT
      REAL*8 AL4,AL5,AL6,XRATMIN
      SAVE AL4,AL5,AL6,XRATMIN

C     Print*, ' Inside Function alphas: AMU =', AMU

      IF(ISTART.EQ.0)THEN
CCPY!*****************************************************************
CCPY!CsB     QCD Lambda for 4 flavors
CCPY!        If (BIPDF.EQ.31) then ! CTEQ3M
CCPY!          NORDER=2
CCPY!          QCDL=0.239
CCPY!        Else If (BIPDF.EQ.32) then ! CTEQ3L
CCPY!          NORDER=1
CCPY!          QCDL=0.177
CCPY!        Else If (BIPDF.EQ.41) then ! CTEQ4M
CCPY!          NORDER=2
CCPY!          QCDL=0.298
CCPY!        Else If (BIPDF.EQ.43) then ! CTEQ4L
CCPY!          NORDER=1
CCPY!          QCDL=0.236
CCPY!        Else If (BIPDF.EQ.51) then ! CTEQ5M
CCPY!          NORDER=2
CCPY!          QCDL=0.326
CCPY!CCPY        Else If (BIPDF.EQ.61 .or. iSet.Eq.1601) then ! CTEQ6M
CCPY!        Else If (BIPDF.GE.200 .AND. BIPDF.LT.500) then ! CTEQ6
CCPY!          NORDER=2
CCPY!          QCDL=0.326
CCPY!CCPY
CCPY        Else If (BIPDF.EQ.903) then ! USING pdf.pds file
CCPY!*****************************************************************
CCPY Aug 2016:
      If (BIPDF.EQ.903) then ! USING pdf.pds file
          NORDER=2
          QCDL=0.326
      ELSE
        Print*, ' Lambda_QCD is set for BIPDF: ', BIPDF
        Print*, ' Try to set PDF to be 903 in the resbos.in file '
        Stop
      ENDIF

!     Print*, ' XMBOT, XMTOP =', XMBOT, XMTOP
!     PRINT*,' NORDER, QCDL =',NORDER, QCDL

        AL4=QCDL
        IF(NORDER.EQ.1) THEN
          AL5=AL4*(AL4/XMBOT)**(2./23.)
          AL6=AL5*(AL5/XMTOP)**(2./21.)
        ELSEIF(NORDER.EQ.2) THEN
          AL5=AL4*(AL4/XMBOT)**(2./23.)*
     1      (LOG(XMBOT**2/AL4**2))**(-963./13225.)
          AL6=AL5*(AL5/XMTOP)**(2./21.)*
     1      (LOG(XMTOP**2/AL5**2))**(-321./3381.)
        ENDIF
        XRATMIN=5.
        B0=12.*3.14159/(33.-8.)
        B1=6.*(153.-19.*4.)/(33.-8.)**2
        ISTART=1
      END IF
      XRAT=AMU**2
      IF(XRAT.LT.XRATMIN)XRAT=XRATMIN
      IF(DSQRT(XRAT).GT.XMTOP)THEN
        RAT=DLOG(XRAT/AL6**2)
        B0=12.*3.14159/21.
        B1=6.*(153.-19.*6.)/21.**2
      ELSEIF(DSQRT(XRAT).GT.XMBOT)THEN
        RAT=DLOG(XRAT/AL5**2)
        B0=12.*3.14159/23.
        B1=6.*(153.-19.*5.)/23.**2
      ELSE
        RAT=DLOG(XRAT/AL4**2)
        B0=12.*3.14159/25.
        B1=6.*(153.-19.*4.)/25.**2
      END IF
      IF(NORDER.EQ.1) THEN
        ALPHAS=B0/RAT
      ELSEIF(NORDER.EQ.2) THEN
        ALPHAS=B0/RAT*(1.-B1*DLOG(RAT)/RAT)
      ENDIF

!     Print*, AL4,AL5,Al6
!     Print*, B0,B1
!     Print*, ALPHAS
            
      RETURN
      END

C**********************************
C THIS IS WLOPH.FOR
C "PHASE SPACE" FOR WPLO, WMLO
C THIS IS NOT A REAL PHASE SPACE ROUTINE.
C AS OF 10/1995, THIS WORKS FOR W+, W-, Z0 AND A0.
C WATCH OUT THE HARD-WIRED CUTS BUILD IN THIS SUBROUTINE, THEY
C ARE OPTIMISED FOR W+ PRODUCTION AT 1.8 TEV TEVATRON.
C --------------------------------------------------------------------------
      SUBROUTINE WLOPH(Y,WT)
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
      REAL*8 P_L1,P_L2,P_Q,P_X,P_Y,P_Z,P_NBAR,P_N,R
      DIMENSION P_L1(5),P_L2(5),P_Q(5),P_X(5),P_Y(5),
     >P_Z(5),P_NBAR(5),P_N(5),R(2)
      REAL*8 Q_V,QT_V,YB_V,THE_STA,PHI_STA
      COMMON/CSFRAME/Q_V,QT_V,YB_V,THE_STA,PHI_STA

      REAL*8 QT_V_MIN, QT_V_MAX,YB_V_MAX,Q_V_MIN,Q_V_MAX,QT_V_MXY
      COMMON/HARDCUT/ QT_V_MIN, QT_V_MAX,YB_V_MAX,
     >  Q_V_MIN,Q_V_MAX,QT_V_MXY

CCPY USED FOR SEPERATING NLO AND NLO_SIG CONTRIBUTIONS
      REAL*8 QT_V_SEP
      COMMON/NLO/ QT_V_SEP
CCAO!
      REAL*8 GAMW_NLO,GAMZ_NLO,ALFA_EM_MZ,ACOUPL,WCOUPL,
     &ZCOUPL,HCOUPL,SW2_EFF,EW_ALFA
      COMMON/STAND7/ GAMW_NLO,GAMZ_NLO,ALFA_EM_MZ,ACOUPL,WCOUPL,
     >       ZCOUPL,HCOUPL,SW2_EFF,EW_ALFA

CCAO! The read W-boson mass and width
      Real*8 XMW_IN,XMW_GRID,GAMW_IN
      COMMON/WMASS_IN/XMW_IN,XMW_GRID,GAMW_IN
      Real*8 XMZ_IN,XMZ_GRID,GAMZ_IN
      COMMON/ZMASS_IN/XMZ_IN,XMZ_GRID,GAMZ_IN

cb      INTEGER J1,J2
      REAL*8 T_NBAR,T_N
      DIMENSION T_NBAR(4),T_N(4)
      LOGICAL TESTING
      REAL*8 QT2_V_MIN,QT2_V_MAX
      REAL*8 Y_EMAX
      Character*40 Boson
      Common / Boson / Boson
!      Integer N_SUD_AB,N_SUD_A,N_SUD_B,N_WIL_C,I_SUD
!      Common / Sudakov / N_SUD_AB,N_SUD_A,N_SUD_B,N_WIL_C,I_SUD
      REAL*8 PCM_L1,ECM_L1,ECM_L2,BETA
      Data Pi/ 3.1415922654 /

C      TESTING=.TRUE.
      TESTING=.FALSE.

C      PI=4.0*ATAN(1.0D0)
      EPS=1.3D-4
      WT=-1.0

CsB * Set 3&4 invariant mass range
      Q_V_MIN=YM12MIN
      Q_V_MAX=YM12MAX

CsB For processes requiring grids:
C   QT_V_MAX and YB_V_MAX are determined in ReadIn and passed through COMMON
C   QT_V_MAX is the maximal qT for the resummed grid.

C FIND Q2_V
      YM1=XMASS(1)
      YM2=XMASS(2)

      YMW=0.d0
      YGW=0.d0
      If (Boson.Eq.'W+' .or. Boson.Eq.'W-' .or.
     .    Boson.Eq.'W+_RA' .or. Boson.Eq.'W-_RA' .or. 
     .    Boson.Eq.'W+A(SV)' .or. Boson.Eq.'W-A(SV)' .or. 
     .    Boson.Eq.'HP' .or. Boson.Eq.'HM') then
        YMW=XMW
        YGW = GAMW_IN
        If (GAMW_IN.LT.0.) YGW = GAMW_NLO
      ELSE IF (Boson.Eq.'Z0' .or. Boson.Eq.'Z0_RA'
     .    .or. Boson.Eq.'DY' .or. Boson.Eq.'HZ') then
        YMW=XMZ
        YGW = GAMZ_IN
        If (GAMZ_IN.LT.0.) YGW = GAMZ_NLO
      ELSE IF (Boson.Eq.'H0') then
        YMW=XMHIGGS
        YGW=GAMH
      ELSE IF (Boson.Eq.'A0' .or. Boson.Eq.'GL'
     .    .or. Boson.Eq.'AA' .or. Boson.Eq.'AG'
     .    .or. Boson.Eq.'ZZ' .or. Boson.Eq.'ZG') then
        YMW=0.d0
        YGW=0.d0
      ENDIF
!     Print*, 'YMW,YGW =', YMW,YGW
CCPY Oct 2008, INCLUDE 'DRELL-YAN PAIR 'DY' PROCESS
      If  (Boson.Eq.'W+' .or. Boson.Eq.'W-' .or.
     .     Boson.Eq.'W+_RA' .or. Boson.Eq.'W-_RA' .or.
     .     Boson.Eq.'W+A(SV)' .or. Boson.Eq.'W-A(SV)' .or.
     .     Boson.Eq.'Z0' .or. Boson.Eq.'Z0_RA' .or. 
     .     Boson.Eq.'H0' .or. Boson.Eq.'DY') then
        IMONTE=2
      else if (Boson.Eq.'AA'.or.Boson.Eq.'AG') then
        imonte=3
      else if (Boson.Eq.'A0') then
        imonte=1
      ELSE
C       IMONTE=1
        IMONTE=0
      END IF

CCPY IN THIS CODE, X1 AND X2 SHOULD BE INTERPRETED AS X_A AND X_B
      IF(IMONTE.EQ.0.or.imonte.eq.3) THEN
C Use linear shape
        Q2_V = Q_V_Min**2 + (Q_V_MAX**2 - Q_V_MIN**2)*Y(1)
        Q_V=DSQRT(Q2_V)
        WT_Q2=(Q_V_MAX**2 - Q_V_MIN**2)
      ELSEIF(IMONTE.EQ.1) THEN
C USE sigma ~ (1/Q2_V)
C HARD-WIRED CUT
        If (Q_V_MIN.Eq.0.d0) Q_V_MIN = 10.d-8
        Y_EMIN=Q_V_MIN
        Y_EMAX=Q_V_MAX
        Q2_V=Y_EMIN**2*(Y_EMAX**2/Y_EMIN**2)**Y(1)
        Q_V=DSQRT(Q2_V)
        WT_Q2=Q2_V*DLOG(Y_EMAX**2/Y_EMIN**2)
      ELSE IF(IMONTE.EQ.2) THEN
C BREIT-WIGNER SHAPE
        ONE=1.D0
        YRSH=W
        YMREST=0.D0
CCPY
C        YMW_MIN=((YM1+YM2)**2-YMW*YMW)/(YMW*YGW)
C        YMW_MAX=((YRSH-YMREST)**2-YMW*YMW)/(YMW*YGW)
C HARD-WIRED CUT
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
      CALL YMAXIMUM(W,Q_V,QT_V,YBMAX)
C HARD-WIRED CUT
      IF(YBMAX.GT.YB_V_MAX) YBMAX=YB_V_MAX

      WHAT=Q_V
      TAU=(WHAT/W)**2
      YBOOST=YBMAX*(2.0D0*Y(2)-1.0D0)
      YB_V=YBOOST
      WT_YB=2.0D0*YBMAX
      X1=DSQRT(TAU)*DEXP(YBOOST)
      X2=TAU/X1

C DO THE LEPTON KINEMATICS IN THE REST FRAME OF W
C       CHANGE THE VARIABLES ( COS(THETA_STAR), PHI_STAR ) IN C.M.
C       FRAME TO ( R(1), R(2) ).
      R(1)=Y(3)
CCPY R(2): FOR PHI_ANGLE IS NOT KNOWN TO VEGAS
      R(2)=RN(1)
      THE_STA=DACOS(2.D0*R(1)-1.D0)
      PHI_STA=2.D0*PI*R(2)
C FOR THE_STA
      WT_THE=2.D0
C FOR PHI_STA
      WT_PHI=2.D0*PI

C THIS IS FOR d(Q^2)d(Y)d(cos)d(phi)
      WT=WT_Q2*WT_YB*WT_THE*WT_PHI

C FOR MASSIVE FINAL STATE PARTICLES (ARISING FROM 2-BODY PHASE SPACE FACTOR).
      CTHE=COS(THE_STA)
      STHE=SIN(THE_STA)
      CPHI=COS(PHI_STA)
      SPHI=SIN(PHI_STA)
         PCM_L1=DSQRT((Q_V**2-(YM1+YM2)**2)*(Q_V**2-(YM1-YM2)**2))/
     >   2.0/Q_V
         ECM_L1=(Q_V**2+YM1**2-YM2**2)/2.0/Q_V
         ECM_L2=(Q_V**2-YM1**2+YM2**2)/2.0/Q_V
      IF(TESTING) THEN
        PRINT*,'YM1,YM2,PCM_L1,ECM_L1,ECM_L2'
        PRINT*,YM1,YM2,PCM_L1,ECM_L1,ECM_L2
      ENDIF
C      BETA=PCM_L1/ECM_L2
      BETA=2.0*PCM_L1/Q_V

C ADD THE OTHER FACTORS
      WT=WT*
     &   1.d0/(2.d0*Pi)**3.d0/8.d0 * BETA *
     &   1.d0/(2*Q_V**2)

      IF(TESTING) THEN
        PRINT*,'X1,X2,TAU,YBOOST'
        PRINT*,X1,X2,TAU,YBOOST
        PRINT*,'Q_V,QT_V,YB_V,THE_STA,PHI_STA'
        PRINT*,Q_V,QT_V,YB_V,THE_STA,PHI_STA
!        PRINT*,'WT_QT2,WT_Q2,WT_YB,WT_THE,WT_PHI,WT'
!        PRINT*,WT_QT2,WT_Q2,WT_YB,WT_THE,WT_PHI,WT
      ENDIF

C**************************************
C In the Collins-Soper frame

      R2INV=1.0D0/SQRT(2.0D0)

      P_NBAR(1)=0.D0
      P_NBAR(2)=0.D0
      P_NBAR(3)=R2INV
      P_NBAR(4)=R2INV

      P_N(1)=0.D0
      P_N(2)=0.D0
      P_N(3)=-R2INV
      P_N(4)=R2INV

C GIVE A RANDOM CHOICE OF AZIMUTHAL ANGLE FOR W-BOSON, IN LAB FRAME
      PHI_Q=2.0*PI*RN(2)
      TM2_Q=Q2_V+QT2_V
      TM_Q=SQRT(TM2_Q)

      P_Q(1)=QT_V*COS(PHI_Q)
      P_Q(2)=QT_V*SIN(PHI_Q)
      P_Q(3)=TM_Q*(EXP(YB_V)-EXP(-YB_V))/2.0D0
      P_Q(4)=TM_Q*(EXP(YB_V)+EXP(-YB_V))/2.0D0

      P_Q(5)=SQRT(-P_Q(1)**2-P_Q(2)**2-P_Q(3)**2+
     >        P_Q(4)**2)

      QPLUS=R2INV*(P_Q(4)+P_Q(3))
      QMINUS=R2INV*(P_Q(4)-P_Q(3))

CsB * Find X, Y, Z vectors
      DO 70 I=1,4
C         P_X(I)=-Q_V/QT_V/TM_Q*(QPLUS*P_NBAR(I)+QMINUS*P_N(I)-
C     >   (TM2_Q/Q2_V)*P_Q(I))
CCPY FOR QT_V=0
         P_X(I)=0.D0
         P_Z(I)=1.0/TM_Q*(QPLUS*P_NBAR(I)-QMINUS*P_N(I))
70    CONTINUE
      P_X(1)=1.0D0

CsB For qT_V = 0: Y = (0,1,0,0)
      CALL FINDY(Q_V,P_Q,P_X,P_Z,P_Y)

      P_X(5)=-P_X(1)**2-P_X(2)**2-P_X(3)**2+
     >        P_X(4)**2
      P_Y(5)=-P_Y(1)**2-P_Y(2)**2-P_Y(3)**2+
     >        P_Y(4)**2
      P_Z(5)=-P_Z(1)**2-P_Z(2)**2-P_Z(3)**2+
     >        P_Z(4)**2

      IF(TESTING) THEN
        PRINT*,'--- P_Q, Q_V = ',P_Q,Q_V
        PRINT*,' P_X =',P_X
        PRINT*,' P_Y =',P_Y
        PRINT*,' P_Z =',P_Z
      ENDIF

      IF(ABS(P_Q(5)-Q_V).GT.EPS) THEN
       PRINT*,' ERROR IN P_Q'
C       Call ExitGracefully
       RETURN
      ENDIF

      IF(ABS(P_X(5)+1.0).GT.EPS .OR. ABS(P_Y(5)+1.0).GT.EPS
     >.OR. ABS(P_Z(5)+1.0).GT.EPS) THEN
       PRINT*,' POTENTIAL ERROR IN MASS: P_X, P_Y  OR  P_Z'
C       Call ExitGracefully
       RETURN
      ELSE
       P_X(5)=SQRT(-P_X(5))
       P_Y(5)=SQRT(-P_Y(5))
       P_Z(5)=SQRT(-P_Z(5))
      ENDIF

CsB * Construct lepton momneta
      DO 80 I=1,4
CCPY Only works for massless particles (i.e. YM1=YM2=0 ).
C         P_L1(I)=Q_V/2.0D0*(P_Q(I)/Q_V+STHE*CPHI*P_X(I)+
C     >   STHE*SPHI*P_Y(I)+CTHE*P_Z(I))
C         P_L2(I)=Q_V/2.0D0*(P_Q(I)/Q_V-STHE*CPHI*P_X(I)-
C     >   STHE*SPHI*P_Y(I)-CTHE*P_Z(I))
C
         P_L1(I)=ECM_L1*P_Q(I)/Q_V+PCM_L1*(STHE*CPHI*P_X(I)+
     >   STHE*SPHI*P_Y(I)+CTHE*P_Z(I))
         P_L2(I)=ECM_L2*P_Q(I)/Q_V+PCM_L1*(-STHE*CPHI*P_X(I)-
     >   STHE*SPHI*P_Y(I)-CTHE*P_Z(I))
80    CONTINUE

CsB * Lepton masses
C-
      P_L1(5)=-P_L1(1)**2-P_L1(2)**2-P_L1(3)**2+
     >        P_L1(4)**2
      P_L2(5)=-P_L2(1)**2-P_L2(2)**2-P_L2(3)**2+
     >        P_L2(4)**2

      IF(TESTING) THEN
        PRINT*,'P_L1(5),P_L2(5)'
        PRINT*,P_L1(5),P_L2(5)
      ENDIF

      P_L1(5)=SQRT(ABS(P_L1(5)))
      P_L2(5)=SQRT(ABS(P_L2(5)))

      IF(YM1.GT.EPS .OR. YM2.GT.EPS) THEN
        IF(ABS(P_L1(5)-YM1).GT.P_L1(5)*EPS .OR.
     >  ABS(P_L2(5)-YM2).GT.P_L2(5)*EPS) THEN
          PRINT*,' POTENTIAL ERROR IN MASS: P_L1 OR P_L2 (1)'
CCPY       Call ExitGracefully
          PRINT *,'P_L1 =  ',P_L1
          PRINT *,'P_L2 =  ',P_L2
          PRINT *,'P_Q =  ',P_Q
          RETURN
        ENDIF
      ELSE
        IF(P_L1(5).GT.EPS .OR. P_L2(5).GT.EPS) THEN
          PRINT*,' POTENTIAL ERROR IN MASS: P_L1 OR P_L2 (2)'
CCPY       Call ExitGracefully
          PRINT *,'P_L1 =  ',P_L1
          PRINT *,'P_L2 =  ',P_L2
          PRINT *,'P_Q =  ',P_Q
          RETURN
        ENDIF
      ENDIF
C-

CCPY TESTING PROGRAM
      IF(TESTING) THEN
        DO 73 I=1,4
           T_NBAR(I)=TM_Q/QPLUS/2.0*(P_Z(I)-QT_V/Q_V*P_X(I)
     >               +TM_Q/Q_V/Q_V*P_Q(I))
           T_N(I)=TM_Q/QMINUS/2.0*(-P_Z(I)-QT_V/Q_V*P_X(I)
     >               +TM_Q/Q_V/Q_V*P_Q(I))
73      CONTINUE
        PRINT*,' T_NBAR =',T_NBAR
        PRINT*,' T_N =',T_N
      ENDIF

CsB * Fill in lab momenta
C-
      DO 90 I=1,4
C NO INFORMATION ON TWO INCOMING PARTONS: SET THEM TO BE ZERO
        PLAB(I,1)=0.D0
        PLAB(I,2)=0.D0
C FERMION
        PLAB(I,3)=P_L1(I)
C ANTI-FERMION
        PLAB(I,4)=P_L2(I)
C VECTOR BOSON
        PLAB(I,9)=P_Q(I)
90    CONTINUE
C-

      IF(TESTING) THEN
C       PRINT *,'((PLAB(J1,J2),J1=1,4),J2=3,4)'
C       PRINT *,((PLAB(J1,J2),J1=1,4),J2=3,4)
        PRINT *,'P_L1 =  ',P_L1
        PRINT *,'P_L2 =  ',P_L2
        PRINT *,'P_Q =  ',P_Q
      ENDIF

      RETURN
      END ! WLOPH

CsB ------------------------------------------------------------------------
      SUBROUTINE BOOST_5(P,R,Q)
CsB ------------------------------------------------------------------------
C
CCPY  THE AXES IN  K  AND  KP   FRAMES ARE PARALLEL, BUT THE VELOCITY
C     OF THE FRAME  KP   IN FRAM   K   IS IN AN ARBITARY DIRECTION.
C
C     THE FOUR VECTOR P IS ASSUMED TO BE GIVEN IN THE REST FRAME OF R,
C     WHICH MUST BE A TIMELIKE VECTOR.
C     THE FIFTH COMPONENT IS MASS.
C     OUTPUT Q IS THE VECTOR P BOOSTED TO THE FRAME IN WHICH R IS GIVEN.

      IMPLICIT REAL*8 (A-H,J-Z)
      DIMENSION P(5),R(5),Q(5),BETA(3)
      X = 0D0
      Y = 0D0
      DO 10 I = 1,3
      BETA(I) = R(I)/R(4)
      X = X + BETA(I)**2
   10 Y = Y + BETA(I)*P(I)
      IF (X.LT.1D-16.OR.X.GE.(1D0-1D-12)) GOTO 30
      GAMMA = 1D0/DSQRT(1D0-X)
      DO 20 I = 1,3
   20 Q(I) = P(I)+BETA(I)*(Y*(GAMMA-1D0)/X + GAMMA*P(4))
      Q(4) = GAMMA*(P(4) + Y)
      RETURN
   30 CONTINUE
      DO 40 I = 1,4
   40 Q(I) = P(I)
      IF(X.GE.(1D0-1D-12)) WRITE(*,1000) R
      RETURN
 1000 FORMAT (' THE REFERENCE VECTOR ',4F10.3,'   IS NOT TIMELIKE.')
      END ! BOOST_5

CsB ------------------------------------------------------------------------
      Subroutine ExitGracefully
CsB ------------------------------------------------------------------------
      Print*, ' Stopped in Exit'
      Stop
      End

C --------------------------------------------------------------------------
      FUNCTION XMAT_PHOT(S,T,U)
C --------------------------------------------------------------------------
C This routine is borrowed from Papageno.
      IMPLICIT NONE
      REAL*8 XMAT_PHOT,S,T,U,XMRP_PPP,XMRP_MPM,XMIP_MPM,
     &  XMRM_PMP,XMIM_PMP,PI
      PI=4.0*ATAN(1.0D0)
      XMRP_PPP=-4.*(2.+2.*(T-U)/S*LOG(T/U)+(T**2+U**2)/S**2
     1 *((LOG(T/U))**2+PI**2))
      XMRP_MPM=-4.*(2.+2.*(S-T)/U*LOG(-S/T)
     1+(S**2+T**2)/U**2*(LOG(-S/T))**2)
      XMIP_MPM=8.*PI*((S-T)/U+(S**2+T**2)/U**2*LOG(-S/T))
      XMRM_PMP=-4.*(2.+2.*(S-U)/T*LOG(-S/U)
     1+(S**2+U**2)/T**2*(LOG(-S/U))**2)
      XMIM_PMP=8.*PI*((S-U)/T+(S**2+U**2)/T**2*LOG(-S/U))
      XMAT_PHOT=2.D0*( XMRP_PPP**2+64.*5.+XMRP_MPM**2+
     1 XMRM_PMP**2+XMIP_MPM**2+XMIM_PMP**2)
      RETURN
      END

CsB --- Begin PYTHIA insert ------------------------------------------------
      FUNCTION GGAA(SH,TH,UH)
C --------------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
C...Stolen shamelessly from Pythia
C...g + g -> gamma + gamma or g + g -> g + gamma
C......
      INTEGER KCHG(6)
      REAL*8 PMAS(6),PARU(1)
      DATA KCHG/-1,2,-1,2,-1,2/
      DATA PMAS/.001D0,.001D0,.2D0,1.6D0,5D0,175d0/

      PARU(1)=4D0*ATAN(1D0)
      A0STUR=0D0
      A0STUI=0D0
      A0TSUR=0D0
      A0TSUI=0D0
      A0UTSR=0D0
      A0UTSI=0D0
      A1STUR=0D0
      A1STUI=0D0
      A2STUR=0D0
      A2STUI=0D0
      ALST=LOG(-SH/TH)
      ALSU=LOG(-SH/UH)
      ALTU=LOG(TH/UH)
      IMAX=5
      SH2=SH*SH
      TH2=TH*TH
      UH2=UH*UH
      DO 1180 I=1,IMAX
        EI=KCHG(I)/3D0
        EIWT=EI**2
        SQMQ=PMAS(I)**2
        EPSS=4D0*SQMQ/SH
        EPST=4D0*SQMQ/TH
        EPSU=4D0*SQMQ/UH
        IF(EPSS.LT.1.D-4) THEN
          B0STUR=1D0+(TH-UH)/SH*ALTU+0.5D0*(TH2+UH2)/SH2*(ALTU**2+
     &    PARU(1)**2)
          B0STUI=0D0
          B0TSUR=1D0+(SH-UH)/TH*ALSU+0.5D0*(SH2+UH2)/TH2*ALSU**2
          B0TSUI=-PARU(1)*((SH-UH)/TH+(SH2+UH2)/TH2*ALSU)
          B0UTSR=1D0+(SH-TH)/UH*ALST+0.5D0*(SH2+TH2)/UH2*ALST**2
          B0UTSI=-PARU(1)*((SH-TH)/UH+(SH2+TH2)/UH2*ALST)
          B1STUR=-1D0
          B1STUI=0D0
          B2STUR=-1D0
          B2STUI=0D0
        ELSE
          CALL PYWAUX(1,EPSS,W1SR,W1SI)
          CALL PYWAUX(1,EPST,W1TR,W1TI)
          CALL PYWAUX(1,EPSU,W1UR,W1UI)
          CALL PYWAUX(2,EPSS,W2SR,W2SI)
          CALL PYWAUX(2,EPST,W2TR,W2TI)
          CALL PYWAUX(2,EPSU,W2UR,W2UI)
          CALL PYI3AU(EPSS,TH/UH,Y3STUR,Y3STUI)
          CALL PYI3AU(EPSS,UH/TH,Y3SUTR,Y3SUTI)
          CALL PYI3AU(EPST,SH/UH,Y3TSUR,Y3TSUI)
          CALL PYI3AU(EPST,UH/SH,Y3TUSR,Y3TUSI)
          CALL PYI3AU(EPSU,SH/TH,Y3USTR,Y3USTI)
          CALL PYI3AU(EPSU,TH/SH,Y3UTSR,Y3UTSI)
          B0STUR=1D0+(1D0+2D0*TH/SH)*W1TR+(1D0+2D0*UH/SH)*W1UR+
     &    0.5D0*((TH2+UH2)/SH2-EPSS)*(W2TR+W2UR)-
     &    0.25D0*EPST*(1D0-0.5D0*EPSS)*(Y3SUTR+Y3TUSR)-
     &    0.25D0*EPSU*(1D0-0.5D0*EPSS)*(Y3STUR+Y3UTSR)+
     &    0.25D0*(-2D0*(TH2+UH2)/SH2+4D0*EPSS+EPST+EPSU+
     &    0.5D0*EPST*EPSU)*(Y3TSUR+Y3USTR)
          B0STUI=(1D0+2D0*TH/SH)*W1TI+(1D0+2D0*UH/SH)*W1UI+
     &    0.5D0*((TH2+UH2)/SH2-EPSS)*(W2TI+W2UI)-
     &    0.25D0*EPST*(1D0-0.5D0*EPSS)*(Y3SUTI+Y3TUSI)-
     &    0.25D0*EPSU*(1D0-0.5D0*EPSS)*(Y3STUI+Y3UTSI)+
     &    0.25D0*(-2D0*(TH2+UH2)/SH2+4D0*EPSS+EPST+EPSU+
     &    0.5D0*EPST*EPSU)*(Y3TSUI+Y3USTI)
          B0TSUR=1D0+(1D0+2D0*SH/TH)*W1SR+(1D0+2D0*UH/TH)*W1UR+
     &    0.5D0*((SH2+UH2)/TH2-EPST)*(W2SR+W2UR)-
     &    0.25D0*EPSS*(1D0-0.5D0*EPST)*(Y3TUSR+Y3SUTR)-
     &    0.25D0*EPSU*(1D0-0.5D0*EPST)*(Y3TSUR+Y3USTR)+
     &    0.25D0*(-2D0*(SH2+UH2)/TH2+4D0*EPST+EPSS+EPSU+
     &    0.5D0*EPSS*EPSU)*(Y3STUR+Y3UTSR)
          B0TSUI=(1D0+2D0*SH/TH)*W1SI+(1D0+2D0*UH/TH)*W1UI+
     &    0.5D0*((SH2+UH2)/TH2-EPST)*(W2SI+W2UI)-
     &    0.25D0*EPSS*(1D0-0.5D0*EPST)*(Y3TUSI+Y3SUTI)-
     &    0.25D0*EPSU*(1D0-0.5D0*EPST)*(Y3TSUI+Y3USTI)+
     &    0.25D0*(-2D0*(SH2+UH2)/TH2+4D0*EPST+EPSS+EPSU+
     &    0.5D0*EPSS*EPSU)*(Y3STUI+Y3UTSI)
          B0UTSR=1D0+(1D0+2D0*TH/UH)*W1TR+(1D0+2D0*SH/UH)*W1SR+
     &    0.5D0*((TH2+SH2)/UH2-EPSU)*(W2TR+W2SR)-
     &    0.25D0*EPST*(1D0-0.5D0*EPSU)*(Y3USTR+Y3TSUR)-
     &    0.25D0*EPSS*(1D0-0.5D0*EPSU)*(Y3UTSR+Y3STUR)+
     &    0.25D0*(-2D0*(TH2+SH2)/UH2+4D0*EPSU+EPST+EPSS+
     &    0.5D0*EPST*EPSS)*(Y3TUSR+Y3SUTR)
          B0UTSI=(1D0+2D0*TH/UH)*W1TI+(1D0+2D0*SH/UH)*W1SI+
     &    0.5D0*((TH2+SH2)/UH2-EPSU)*(W2TI+W2SI)-
     &    0.25D0*EPST*(1D0-0.5D0*EPSU)*(Y3USTI+Y3TSUI)-
     &    0.25D0*EPSS*(1D0-0.5D0*EPSU)*(Y3UTSI+Y3STUI)+
     &    0.25D0*(-2D0*(TH2+SH2)/UH2+4D0*EPSU+EPST+EPSS+
     &    0.5D0*EPST*EPSS)*(Y3TUSI+Y3SUTI)
          B1STUR=-1D0-0.25D0*(EPSS+EPST+EPSU)*(W2SR+W2TR+W2UR)+
     &    0.25D0*(EPSU+0.5D0*EPSS*EPST)*(Y3SUTR+Y3TUSR)+
     &    0.25D0*(EPST+0.5D0*EPSS*EPSU)*(Y3STUR+Y3UTSR)+
     &    0.25D0*(EPSS+0.5D0*EPST*EPSU)*(Y3TSUR+Y3USTR)
          B1STUI=-0.25D0*(EPSS+EPST+EPSU)*(W2SI+W2TI+W2UI)+
     &    0.25D0*(EPSU+0.5D0*EPSS*EPST)*(Y3SUTI+Y3TUSI)+
     &    0.25D0*(EPST+0.5D0*EPSS*EPSU)*(Y3STUI+Y3UTSI)+
     &    0.25D0*(EPSS+0.5D0*EPST*EPSU)*(Y3TSUI+Y3USTI)
          B2STUR=-1D0+0.125D0*EPSS*EPST*(Y3SUTR+Y3TUSR)+
     &    0.125D0*EPSS*EPSU*(Y3STUR+Y3UTSR)+
     &    0.125D0*EPST*EPSU*(Y3TSUR+Y3USTR)
          B2STUI=0.125D0*EPSS*EPST*(Y3SUTI+Y3TUSI)+
     &    0.125D0*EPSS*EPSU*(Y3STUI+Y3UTSI)+
     &    0.125D0*EPST*EPSU*(Y3TSUI+Y3USTI)
        ENDIF
        A0STUR=A0STUR+EIWT*B0STUR
        A0STUI=A0STUI+EIWT*B0STUI
        A0TSUR=A0TSUR+EIWT*B0TSUR
        A0TSUI=A0TSUI+EIWT*B0TSUI
        A0UTSR=A0UTSR+EIWT*B0UTSR
        A0UTSI=A0UTSI+EIWT*B0UTSI
        A1STUR=A1STUR+EIWT*B1STUR
        A1STUI=A1STUI+EIWT*B1STUI
        A2STUR=A2STUR+EIWT*B2STUR
        A2STUI=A2STUI+EIWT*B2STUI
 1180 CONTINUE
      ASQSUM=A0STUR**2+A0STUI**2+A0TSUR**2+A0TSUI**2+A0UTSR**2+
     &    A0UTSI**2+4D0*A1STUR**2+4D0*A1STUI**2+A2STUR**2+A2STUI**2

      GGAA=2D0*ASQSUM*64D0
      RETURN
      END
CCC
C
CCC
      SUBROUTINE PYWAUX(IAUX,EPS,WRE,WIM)

C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      REAL*8 PARU(1)

      ASINH(X)=LOG(X+SQRT(X**2+1D0))
      ACOSH(X)=LOG(X+SQRT(X**2-1D0))
      PARU(1)=4D0*ATAN(1D0)
      IF(EPS.LT.0D0) THEN
        IF(IAUX.EQ.1) WRE=2D0*SQRT(1D0-EPS)*ASINH(SQRT(-1D0/EPS))
        IF(IAUX.EQ.2) WRE=4D0*(ASINH(SQRT(-1D0/EPS)))**2
        WIM=0D0
      ELSEIF(EPS.LT.1D0) THEN
        IF(IAUX.EQ.1) WRE=2D0*SQRT(1D0-EPS)*ACOSH(SQRT(1D0/EPS))
        IF(IAUX.EQ.2) WRE=4D0*(ACOSH(SQRT(1D0/EPS)))**2-PARU(1)**2
        IF(IAUX.EQ.1) WIM=-PARU(1)*SQRT(1D0-EPS)
        IF(IAUX.EQ.2) WIM=-4D0*PARU(1)*ACOSH(SQRT(1D0/EPS))
      ELSE
        IF(IAUX.EQ.1) WRE=2D0*SQRT(EPS-1D0)*ASIN(SQRT(1D0/EPS))
        IF(IAUX.EQ.2) WRE=-4D0*(ASIN(SQRT(1D0/EPS)))**2
        WIM=0D0
      ENDIF

      RETURN
      END
CCC
C
CCC
      SUBROUTINE PYI3AU(EPS,RAT,Y3RE,Y3IM)

C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)

      REAL*8 PARU(1)

      PARU(1)=4D0*ATAN(1D0)
      BE=0.5D0*(1D0+SQRT(1D0+RAT*EPS))
      IF(EPS.LT.1D0) GA=0.5D0*(1D0+SQRT(1D0-EPS))

      IF(EPS.LT.0D0) THEN
        IF(ABS(EPS).LT.1.D-4.AND.ABS(RAT*EPS).LT.1.D-4) THEN
          F3RE=PYSPEN(-0.25D0*EPS/(1D0+0.25D0*(RAT-1D0)*EPS),0D0,1)-
     &    PYSPEN((1D0-0.25D0*EPS)/(1D0+0.25D0*(RAT-1D0)*EPS),0D0,1)+
     &    PYSPEN(0.25D0*(RAT+1D0)*EPS/(1D0+0.25D0*RAT*EPS),0D0,1)-
     &    PYSPEN((RAT+1D0)/RAT,0D0,1)+0.5D0*(LOG(1D0+0.25D0*RAT*EPS)**2-
     &    LOG(0.25D0*RAT*EPS)**2)+LOG(1D0-0.25D0*EPS)*
     &    LOG((1D0+0.25D0*(RAT-1D0)*EPS)/(1D0+0.25D0*RAT*EPS))+
     &    LOG(-0.25D0*EPS)*LOG(0.25D0*RAT*EPS/(1D0+0.25D0*(RAT-1D0)*
     &    EPS))
        ELSEIF(ABS(EPS).LT.1.D-4.AND.ABS(RAT*EPS).GE.1.D-4) THEN
          F3RE=PYSPEN(-0.25D0*EPS/(BE-0.25D0*EPS),0D0,1)-
     &    PYSPEN((1D0-0.25D0*EPS)/(BE-0.25D0*EPS),0D0,1)+
     &    PYSPEN((BE-1D0+0.25D0*EPS)/BE,0D0,1)-
     &    PYSPEN((BE-1D0+0.25D0*EPS)/(BE-1D0),0D0,1)+
     &    0.5D0*(LOG(BE)**2-LOG(BE-1D0)**2)+
     &    LOG(1D0-0.25D0*EPS)*LOG((BE-0.25D0*EPS)/BE)+
     &    LOG(-0.25D0*EPS)*LOG((BE-1D0)/(BE-0.25D0*EPS))
        ELSEIF(ABS(EPS).GE.1.D-4.AND.ABS(RAT*EPS).LT.1.D-4) THEN
          F3RE=PYSPEN((GA-1D0)/(GA+0.25D0*RAT*EPS),0D0,1)-
     &    PYSPEN(GA/(GA+0.25D0*RAT*EPS),0D0,1)+
     &    PYSPEN((1D0+0.25D0*RAT*EPS-GA)/(1D0+0.25D0*RAT*EPS),0D0,1)-
     &    PYSPEN((1D0+0.25D0*RAT*EPS-GA)/(0.25D0*RAT*EPS),0D0,1)+
     &    0.5D0*(LOG(1D0+0.25D0*RAT*EPS)**2-LOG(0.25D0*RAT*EPS)**2)+
     &    LOG(GA)*LOG((GA+0.25D0*RAT*EPS)/(1D0+0.25D0*RAT*EPS))+
     &    LOG(GA-1D0)*LOG(0.25D0*RAT*EPS/(GA+0.25D0*RAT*EPS))
        ELSE
          F3RE=PYSPEN((GA-1D0)/(GA+BE-1D0),0D0,1)-
     &    PYSPEN(GA/(GA+BE-1D0),0D0,1)+PYSPEN((BE-GA)/BE,0D0,1)-
     &    PYSPEN((BE-GA)/(BE-1D0),0D0,1)+0.5D0*(LOG(BE)**2-
     &    LOG(BE-1D0)**2)+LOG(GA)*LOG((GA+BE-1D0)/BE)+
     &    LOG(GA-1D0)*LOG((BE-1D0)/(GA+BE-1D0))
        ENDIF
        F3IM=0D0
      ELSEIF(EPS.LT.1D0) THEN
        IF(ABS(EPS).LT.1.D-4.AND.ABS(RAT*EPS).LT.1.D-4) THEN
          F3RE=PYSPEN(-0.25D0*EPS/(1D0+0.25D0*(RAT-1D0)*EPS),0D0,1)-
     &    PYSPEN((1D0-0.25D0*EPS)/(1D0+0.25D0*(RAT-1D0)*EPS),0D0,1)+
     &    PYSPEN((1D0-0.25D0*EPS)/(-0.25D0*(RAT+1D0)*EPS),0D0,1)-
     &    PYSPEN(1D0/(RAT+1D0),0D0,1)+LOG((1D0-0.25D0*EPS)/
     &    (0.25D0*EPS))*LOG((1D0+0.25D0*(RAT-1D0)*EPS)/
     &    (0.25D0*(RAT+1D0)*EPS))
          F3IM=-PARU(1)*LOG((1D0+0.25D0*(RAT-1D0)*EPS)/
     &    (0.25D0*(RAT+1D0)*EPS))
        ELSEIF(ABS(EPS).LT.1.D-4.AND.ABS(RAT*EPS).GE.1.D-4) THEN
          F3RE=PYSPEN(-0.25D0*EPS/(BE-0.25D0*EPS),0D0,1)-
     &    PYSPEN((1D0-0.25D0*EPS)/(BE-0.25D0*EPS),0D0,1)+
     &    PYSPEN((1D0-0.25D0*EPS)/(1D0-0.25D0*EPS-BE),0D0,1)-
     &    PYSPEN(-0.25D0*EPS/(1D0-0.25D0*EPS-BE),0D0,1)+
     &    LOG((1D0-0.25D0*EPS)/(0.25D0*EPS))*
     &    LOG((BE-0.25D0*EPS)/(BE-1D0+0.25D0*EPS))
          F3IM=-PARU(1)*LOG((BE-0.25D0*EPS)/(BE-1D0+0.25D0*EPS))
        ELSEIF(ABS(EPS).GE.1.D-4.AND.ABS(RAT*EPS).LT.1.D-4) THEN
          F3RE=PYSPEN((GA-1D0)/(GA+0.25D0*RAT*EPS),0D0,1)-
     &    PYSPEN(GA/(GA+0.25D0*RAT*EPS),0D0,1)+
     &    PYSPEN(GA/(GA-1D0-0.25D0*RAT*EPS),0D0,1)-
     &    PYSPEN((GA-1D0)/(GA-1D0-0.25D0*RAT*EPS),0D0,1)+
     &    LOG(GA/(1D0-GA))*LOG((GA+0.25D0*RAT*EPS)/
     &    (1D0+0.25D0*RAT*EPS-GA))
          F3IM=-PARU(1)*LOG((GA+0.25D0*RAT*EPS)/
     &    (1D0+0.25D0*RAT*EPS-GA))
        ELSE
          F3RE=PYSPEN((GA-1D0)/(GA+BE-1D0),0D0,1)-
     &    PYSPEN(GA/(GA+BE-1D0),0D0,1)+PYSPEN(GA/(GA-BE),0D0,1)-
     &    PYSPEN((GA-1D0)/(GA-BE),0D0,1)+LOG(GA/(1D0-GA))*
     &    LOG((GA+BE-1D0)/(BE-GA))
          F3IM=-PARU(1)*LOG((GA+BE-1D0)/(BE-GA))
        ENDIF
      ELSE
        RSQ=EPS/(EPS-1D0+(2D0*BE-1D0)**2)
        RCTHE=RSQ*(1D0-2D0*BE/EPS)
        RSTHE=SQRT(MAX(0D0,RSQ-RCTHE**2))
        RCPHI=RSQ*(1D0+2D0*(BE-1D0)/EPS)
        RSPHI=SQRT(MAX(0D0,RSQ-RCPHI**2))
        R=SQRT(RSQ)
        THE=ACOS(MAX(-0.999999D0,MIN(0.999999D0,RCTHE/R)))
        PHI=ACOS(MAX(-0.999999D0,MIN(0.999999D0,RCPHI/R)))
        F3RE=PYSPEN(RCTHE,RSTHE,1)+PYSPEN(RCTHE,-RSTHE,1)-
     &  PYSPEN(RCPHI,RSPHI,1)-PYSPEN(RCPHI,-RSPHI,1)+
     &  (PHI-THE)*(PHI+THE-PARU(1))
        F3IM=PYSPEN(RCTHE,RSTHE,2)+PYSPEN(RCTHE,-RSTHE,2)-
     &  PYSPEN(RCPHI,RSPHI,2)-PYSPEN(RCPHI,-RSPHI,2)
      ENDIF

      Y3RE=2D0/(2D0*BE-1D0)*F3RE
      Y3IM=2D0/(2D0*BE-1D0)*F3IM

      RETURN
      END
CCC
C
CCC
      FUNCTION PYSPEN(XREIN,XIMIN,IREIM)

C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      REAL*8 PARU(1)

C...Local array and data.
      DIMENSION B(0:14)
      DATA B/
     &1.000000D+00,        -5.000000D-01,         1.666667D-01,
     &0.000000D+00,        -3.333333D-02,         0.000000D+00,
     &2.380952D-02,         0.000000D+00,        -3.333333D-02,
     &0.000000D+00,         7.575757D-02,         0.000000D+00,
     &-2.531135D-01,         0.000000D+00,         1.166667D+00/

      PARU(1)=4D0*ATAN(1D0)
      XRE=XREIN
      XIM=XIMIN
      IF(ABS(1D0-XRE).LT.1.D-6.AND.ABS(XIM).LT.1.D-6) THEN
        IF(IREIM.EQ.1) PYSPEN=PARU(1)**2/6D0
        IF(IREIM.EQ.2) PYSPEN=0D0
        RETURN
      ENDIF

      XMOD=SQRT(XRE**2+XIM**2)
      IF(XMOD.LT.1.D-6) THEN
        IF(IREIM.EQ.1) PYSPEN=0D0
        IF(IREIM.EQ.2) PYSPEN=0D0
        RETURN
      ENDIF

      XARG=SIGN(ACOS(XRE/XMOD),XIM)
      SP0RE=0D0
      SP0IM=0D0
      SGN=1D0
      IF(XMOD.GT.1D0) THEN
        ALGXRE=LOG(XMOD)
        ALGXIM=XARG-SIGN(PARU(1),XARG)
        SP0RE=-PARU(1)**2/6D0-(ALGXRE**2-ALGXIM**2)/2D0
        SP0IM=-ALGXRE*ALGXIM
        SGN=-1D0
        XMOD=1D0/XMOD
        XARG=-XARG
        XRE=XMOD*COS(XARG)
        XIM=XMOD*SIN(XARG)
      ENDIF
      IF(XRE.GT.0.5D0) THEN
        ALGXRE=LOG(XMOD)
        ALGXIM=XARG
        XRE=1D0-XRE
        XIM=-XIM
        XMOD=SQRT(XRE**2+XIM**2)
        XARG=SIGN(ACOS(XRE/XMOD),XIM)
        ALGYRE=LOG(XMOD)
        ALGYIM=XARG
        SP0RE=SP0RE+SGN*(PARU(1)**2/6D0-(ALGXRE*ALGYRE-ALGXIM*ALGYIM))
        SP0IM=SP0IM-SGN*(ALGXRE*ALGYIM+ALGXIM*ALGYRE)
        SGN=-SGN
      ENDIF

      XRE=1D0-XRE
      XIM=-XIM
      XMOD=SQRT(XRE**2+XIM**2)
      XARG=SIGN(ACOS(XRE/XMOD),XIM)
      ZRE=-LOG(XMOD)
      ZIM=-XARG

      SPRE=0D0
      SPIM=0D0
      SAVERE=1D0
      SAVEIM=0D0
      DO 100 I=0,14
        IF(MAX(ABS(SAVERE),ABS(SAVEIM)).LT.1D-30) GOTO 110
        TERMRE=(SAVERE*ZRE-SAVEIM*ZIM)/DBLE(I+1)
        TERMIM=(SAVERE*ZIM+SAVEIM*ZRE)/DBLE(I+1)
        SAVERE=TERMRE
        SAVEIM=TERMIM
        SPRE=SPRE+B(I)*TERMRE
        SPIM=SPIM+B(I)*TERMIM
  100 CONTINUE

  110 IF(IREIM.EQ.1) PYSPEN=SP0RE+SGN*SPRE
      IF(IREIM.EQ.2) PYSPEN=SP0IM+SGN*SPIM

      RETURN
      END
CsB --- End PYTHIA insert -------------------------------------------------

CsB For G G -> A A G
C -------------------------------------------------------------------------
      SUBROUTINE GGAAG(WT)
C -------------------------------------------------------------------------
      Implicit none

      DOUBLE PRECISION W,QCDL,X1X2MIN
      INTEGER*4 NJET,N_PRIMARY_PART,NF,I_STR,IBEAM,NPART,NEVENT,NORDER,
     &IW,N_SECOND_PART
      COMMON/AA/W,QCDL,X1X2MIN,NJET,N_PRIMARY_PART,NF,I_STR,IBEAM,NPART,
     &  NEVENT,NORDER,N_SECOND_PART
      INTEGER*4 IERR_WT
      COMMON/MYAMP/IERR_WT
      DOUBLE PRECISION AmpSq, tmp, PentagonGGAAG, PentagonQGAAQ
      integer iproc
      Common / PartonProcess / iProc

      double precision qQBUL,qQBDL,QBqUL,QBqDL,qGUL,qGDL,GqUL,GqDL,GGL
      double precision WT, WT1, WT2, WT3

      DOUBLE PRECISION Q_V,QT_V,YB_V,THE_STA,PHI_STA
      COMMON/CSFRAME/Q_V,QT_V,YB_V,THE_STA,PHI_STA

      Logical Debug
      Data Debug /.False./
C      Debug = .True.

C THIS IS FOR THE PROCESS
C
C     G(1) + G(2) --> A(3) A(4) G(5)
C
C HERE IS THE AVERAGED (AMPLITUDE)**2, THE AMPLITUDE INCLUDES BOTH THE
C COLOR AND SPIN FACTORS. A factor 1/2 is also included for
C identical particles in the final state.
C
C    U=A(1), D=A(2), UBAR=A(3),DBAR=A(9),GLUE=A(4),S=SBAR=A(5),
C    C=CBAR=A(6), B=BBAR=A(7), T=TBAR=A(8).
C    ALL DISTRIBUTIONS IN THE PROTON.
C

cpn2006 Parton luminosities for hadrons 1 and 2 
      call PartonMomentumLumi(qQBUL,qQBDL,QBqUL,QBqDL,
     >  qGUL,qGDL,GqUL,GqDL,GGL)


      WT1=0d0; WT2=0d0; WT3=0d0
CsB___Initialize flag for error of the weight
      IERR_WT=-1

      if (iproc.ne.2) then 
C_____Calculating the squared matrix element of the real emission process:
C     For G G -> A A G:
        AmpSq = 0d0
cpn2006 PentagonGGAAG(i1,i2,i3,i4,i5) is the amplitude for 
c G(i1)+G(i2) -> A(i3)+A(i4)+G(i5)
        AmpSq=PentagonGGAAG(1,2,3,4,5)
!     Print*, ' AmpSq1 = ', AmpSq
!     tmp = AmpSq
!      Call gggAAA2(AmpSq)
!     Print*, ' AmpSq2/AmpSq1 = ',AmpSq/tmp
!     Call TestggaagRoutines
!     Stop

CsB____G G -> A A G cross section
        WT1 = GGL*AmpSq
      endif !iproc.ne.2

      if (iproc.ne.1) then 
cpn2006 PentagonQGAAQ(i1,i2,i3,i4,i5) is the amplitude for 
c Q(i4)+G(i5) -> A(i1)+A(i2)+G(i3)

CPN____Q G -> A A Q cross section
        AmpSq=PentagonQGAAQ(1,2,3,4,5)
        WT2 = (qGUL+qGDL)*AmpSq

CPN____G Q -> A A Q cross section
        AmpSq=PentagonQGAAQ(2,1,3,4,5)
        WT3 = (GqUL+GqDL)*AmpSq
      endif !iproc.ne.2


      If (iProc.Eq.0) then
        WT = WT1+WT2+WT3        ! All Contributions
      ElseIf (iProc.Eq.1) then  ! GG only
        WT = WT1
      ElseIf (iProc.Eq.2) then ! q G + Gq  only
        WT = WT2 + WT3
      ElseIf (iProc.Eq.3) then ! q G  only
        WT = WT2
      ElseIf (iProc.Eq.4) then ! G q only
        WT = WT3
      Else
        Print*, ' iProc = ', iProc, ' value in not defined. '
        Stop
      EndIf
      If (Debug) Print*, ' WT1 = ', WT1
!     If (Debug) Stop

      RETURN
      END

CsB -----------------------------------------------------------------------

CsB For q Q or q G -> di-EW gauge-boson + parton
C -------------------------------------------------------------------------
      SUBROUTINE QQBZZG(WT)
C -------------------------------------------------------------------------
      IMPLICIT NONE

      REAL*8 WT
      REAL*8 W,QCDL,X1X2MIN
      INTEGER*4 NJET,N_PRIMARY_PART,NF,I_STR,IBEAM,NPART,NEVENT,NORDER,
     ,          N_SECOND_PART
      INTEGER*4 IW
      REAL*8 PLAB,PStore(4,5),PCEM,WHAT,X1,X2,Q,Q1,Q2
      INTEGER*4 IERR_WT
      REAL*8 XMASS,X_DECAYMASS

      COMMON/AA/W,QCDL,X1X2MIN,NJET,N_PRIMARY_PART,NF,I_STR,IBEAM,NPART,
     &  NEVENT,NORDER,N_SECOND_PART
      COMMON/CHARGE/IW
      COMMON/MOM/PLAB(4,12),PCEM(4,12),WHAT
      COMMON/PMF/X1,X2
      COMMON/SCALE/Q
      COMMON/MYAMP/IERR_WT
      COMMON/XPART_MASS/XMASS(100),X_DECAYMASS(100)

      REAL*8 XMW,XMZ,XMW2,XMZ2,AEM,GAMW,GAMZ,GAMW2,GAMZ2,
     >SWS,CWS,XMBOT,XMTOP,XMB2,XMT2,GWEAK2,GAMT
      COMMON/STAND0/XMW,XMZ,XMW2,XMZ2,AEM,GAMW,GAMZ,GAMW2,GAMZ2
      COMMON/STAND1/SWS,CWS,XMBOT,XMTOP,XMB2,XMT2,GWEAK2
      COMMON/HTOP/GAMT

      REAL*8 SF,SFp,SFn
      DIMENSION SF(50),SFp(50),SFn(50)
      REAL*8 ALPHAS,ALFAS,PKEEP(4)
      REAL*8 PI,WT1,WT2,WT3,WT4
      INTEGER J1
      EXTERNAL ALPHAS

      REAL*8 P1,P2,P3,P4,P5,AMP2_up,AMP2_dn,spin,color,fact_iden
      Real*8 Amp2_Sub_up,Amp2_Sub_dn
      Real*8 CPAMP2_UP,CPAMP2_DN
      INTEGER I,J
      Real*8 g_fVf_P, g_fVf_M
      Common / Couplings / g_fVf_P, g_fVf_M
      Character*40 Boson
      Common / Boson / Boson
      Real*8 VVMass
      Common / VVMass / VVMass
      INTEGER IJUMP, iProc
      Common / PartonProcess / iProc
      Integer i1,i2,i3,i4,i5
      Common / Indices / i1,i2,i3,i4,i5
      Integer iSub
      Real*8 cutM45, SubtractionShape
      logical SubQGFlag
      Common / Subtraction / cutM45, SubtractionShape(3:4),iSub,
     >  SubQGFlag(3:4)
CsB___q(B) G color factor: sum and average
      Real*8 COLORqG
C      Data COLORqG / 0.16666666666666666 / ! COLORqG = 4.0/3.0/8.0
      Real*8 Frac_n,Frac_p,
     .       qQBUL,qQBDL,QBqUL,QBqDL,qGUL,qGDL,GqUL,GqDL,GGL
      Real*8 g_fVf_P_u,g_fVf_M_u,g_fVf_P_d,g_fVf_M_d
      Common / Stand6 / g_fVf_P_u,g_fVf_M_u,g_fVf_P_d,g_fVf_M_d
      Logical Debug
      Data i1,i2,i3,i4,i5 / 1,2,3,4,5 /
      Data Pi /3.141592653589793/
      Data Debug /.False./
C      Debug = .True.

CsB___Store PLAB
cpn      If (Debug) then
c        Print*, ' Storing momenta... '
        Do i = 1,5
          Do j = 1,4
            PStore(j,i) = PLab(j,i)
          End Do
        End Do
c      End If

C THIS IS FOR THE PROCESS
C
C     q(1) + qbar(2) --> Z(3) Z(4) g(5)
C
C HERE IS THE AVERAGED (AMPLITUDE)**2, THE AMPLITUDE INCLUDES BOTH THE
C COLOR AND SPIN FACTORS. A factor 1/2 is also included for
C identical particles in the final state.
C

      call PartonMomentumLumi(qQBUL,qQBDL,QBqUL,QBqDL,
     >  qGUL,qGDL,GqUL,GqDL,GGL)

CsB___Initialize flag for error of the weight
      IERR_WT=-1

CsB___The strong coulping
      ALFAS=ALPHAS(Q)

CsB___Fermion-VB-fermion couplings and individual mass for the produced
C     EW-gauge-boson
      If (Boson.Eq.'AA') then
        g_fVf_P = 2.d0/3.d0 ! for f = u type quark and V = photon
        g_fVf_M = 2.d0/3.d0 ! for f = u type quark and V = photon
        VVMass = 0.d0
      Else If (Boson.Eq.'ZZ') then
        g_fVf_P = g_fVf_P_u ! for f = u type quark and V = Z
        g_fVf_M = g_fVf_M_u ! for f = u type quark and V = Z
        VVMass = xMz
      End If

C Calculating the squared matrix element of the real emission process:
C    q(p(1)) qB(p(2)) -> Z(p(3)) Z(p(4)) G(p(5)).
C Based on PRD43('91)3637 (Ohnemus-Owens).

CsB Interchange p3 and p4 to test symmetry of amplitude
c      Do J1=1,4
c        PKEEP(J1)=PLAB(J1,3)
c        PLAB(J1,3)=PLAB(J1,4)
c        PLAB(J1,4)=PKEEP(J1)
c      End Do

CsB Check indices
      If (Debug) Print*, ' i1,i2,i3,i4,i5: ', i1,i2,i3,i4,i5

CsB For q qB -> V V G the formula (D3) apply without change:
      CALL Amp2qqBZZG (alfas, Amp2_up,Amp2_dn)
      If (Debug) then
        Print*, ' Amp2_up, Amp2_dn    = ', Amp2_up,Amp2_dn
C        Call AMP2_QQVVG(CPAMP2_UP,CPAMP2_DN,PI,ALFAS,XMZ,GWEAK2,SWS,
C     >                  PLAB)
CsB C.-P.'s amplitude square contains the spin and color factors
C        CPAMP2_UP = CPAMP2_UP/(1.0/2.0/2.0*4.0/3.0/3.0)
C        CPAMP2_DN = CPAMP2_DN/(1.0/2.0/2.0*4.0/3.0/3.0)
C        If (Abs(Amp2_up/CPAMP2_UP - 1.0).GT.0.001) then
C          Print*, ' CPAMP2_UP,CPAMP2_DN = ',
C     >              CPAMP2_UP,CPAMP2_DN,AMP2_UP/CPAMP2_UP
C        End If
      End If

C_____AVERAGE OVER SPIN AND COLOR
      SPIN=1.0/2.0/2.0
C_____Color factor: q qB color sum and average
      COLOR=4.0/3.0/3.0
C_____FOR IDENTITICAL PARTICLES IN FINAL STATE
      FACT_IDEN=1.0/2.0

CsB____q/p qB/pB cross section
       WT1 = qQBUL*Amp2_up + qQBDL*Amp2_dn
       WT1=WT1*SPIN*COLOR*FACT_IDEN

CCPY FOR DI-PHOTON PRODUCTION
       IJUMP=1
       IF(IJUMP.EQ.1) THEN
CsB qB/p q/pB matrix element
         WT2 = QBqUL*Amp2_up + QBqDL*Amp2_dn
         WT2=WT2*SPIN*COLOR*FACT_IDEN
         GO TO 888
       ENDIF

C For PPBAR, switch the incoming parton legs
C Has to be done for non-identical particle final states
      DO 92 J1=1,4
       PKEEP(J1)=PLAB(J1,1)
       PLAB(J1,1)=PLAB(J1,2)
       PLAB(J1,2)=PKEEP(J1)
92    CONTINUE
CCsB Csossing could also be done switching indices: 1<->2
C      i1 = 2
C      i2 = 1

      CALL Amp2qqBZZG (alfas, Amp2_up,Amp2_dn)
C      If (Debug) then
C        Print*, ' Amp2_up,Amp2_dn     : ', Amp2_up,Amp2_dn
C        Call AMP2_QQVVG(CPAMP2_UP,CPAMP2_DN,PI,ALFAS,XMZ,GWEAK2,SWS,
C     >                  PLAB)
C        CPAMP2_UP = CPAMP2_UP/(1.0/2.0/2.0*4.0/3.0/3.0)
C        CPAMP2_DN = CPAMP2_DN/(1.0/2.0/2.0*4.0/3.0/3.0)
C        Print*, ' CPAMP2_UP,CPAMP2_DN = ',
C     >            CPAMP2_UP,CPAMP2_DN,AMP2_UP/CPAMP2_UP
C      End If

C AVERAGE OVER SPIN AND COLOR
      SPIN=1.0/2.0/2.0
C check this color factor
      COLOR=4.0/3.0/3.0
C FOR IDENTITICAL PARTICLES IN FINAL STATE
      FACT_IDEN=1.0/2.0

       WT2=SF(21)*Amp2_up+SF(22)*Amp2_dn
       WT2=WT2*SPIN*COLOR*FACT_IDEN

c      IERR_WT=-1
c      IF(WT1.LT.0.0 .OR. WT2.LT.0.0) THEN
c       IERR_WT=1
c       WT=0.0D0
c      ENDIF

C PUT BACK THE ORIGINAL PLAB'S
      DO 94 J1=1,4
        PKEEP(J1)=PLAB(J1,1)
        PLAB(J1,1)=PLAB(J1,2)
        PLAB(J1,2)=PKEEP(J1)
94    CONTINUE

      IF(WT2.LT.0.0) THEN
        IERR_WT=1
        WT=0.0D0
        RETURN
      ENDIF

  888 CONTINUE
      IF(WT1.LT.0.0) THEN
        IERR_WT=1
        WT=0.0D0
        RETURN
      ENDIF

CsB Before doing crossing, do a cross check:
C   Switching momentum values should be the same as swithing the indices
C   in formulae (D3).
      If (Debug) then
        CALL Amp2qqBZZG (alfas, Amp2_up,Amp2_dn)
        Print*, ' Amp2_up,Amp2_dn     : ', Amp2_up,Amp2_dn
C Now cross momenta: p2 <-> p5 and ...
        Do J1=1,4
          PKEEP(J1)=PLAB(J1,2)
          PLAB(J1,2)=PLAB(J1,5)
          PLAB(J1,5)=PKEEP(J1)
        End Do
C Show momenta:
C        Do i=1,5
C        Do j=1,4
C          Print*, ' PLAB(',j,',',i,') = ', PLAB(j,i)
C        End Do
C        End Do
C ... cross indices: 2 <-> 5:
        Print*, ' i1,i2,i3,i4,i5: ', i1,i2,i3,i4,i5
        i2 = 5
        i5 = 2
        CALL Amp2qqBZZG (alfas, Amp2_up,Amp2_dn)
        Print*, ' Amp2_up,Amp2_dn     : ', Amp2_up,Amp2_dn
Csb The amplitudes are the same as the original.
C   This suggests that the crossing is implemented correctly.
C
ccC If indices are not crossed the amplitude blows up, because the appropriate
ccC invariant tij (this case t15) = 0.
cc        i2 = 2
cc        i5 = 5
cc        CALL Amp2qqBZZG (alfas, Amp2_up,Amp2_dn)
cc        Print*, ' Amp2_up,Amp2_dn     : ', Amp2_up,Amp2_dn
cdump        Stop
      End If

CsB Get q G amplitude from q qB by crossing (D5)
C   Switch p2 <-> -p5
      Do J1=1,4
        PKEEP(J1)=PLAB(J1,2)
        PLAB(J1,2)=-PLAB(J1,5)
        PLAB(J1,5)=-PKEEP(J1)
      End Do
C      Do J1=1,4
C        PLAB(J1,2)=-pStore(J1,5)
C        PLAB(J1,5)=-pStore(J1,2)
C      End Do
C Crossing could be done crossing indices 2 <-> 5...
C      i2 = 5
C      i5 = 2
C ... and the signs of p(2) and p(5):
C      Do J1=1,4
C        PLAB(J1,2)=-PLAB(J1,2)
C        PLAB(J1,5)=-PLAB(J1,5)
C      End Do
C This change should also be restored after the crossed amplitude calculated.
C From now on we only do crossing switching the appropriate momentum values.

      CALL Amp2qqBZZG (alfas, Amp2_up,Amp2_dn)
C      If (Debug) then
C        Print*, ' Amp2_up,Amp2_dn     : ', Amp2_up,Amp2_dn
C      End If

CsB 9/9/97
C   For the q G -> photon photon q process we define a subtraction term which
C   subtracts the collinear singular term when one of the final state photon is
C   emitted from the final state quark leg.
C   This contribution is added back but the q -> q photon splitting function
C   replaced with the LL q -> q photon fragmentation function and an
C   additional cut on the invariant mass of the collinear pair (M > 1-2 GeV).
C   We use the fragmentation function by Bailey, Ohnemus and Owens.
CsB___iSub=0 no subtraction is done
C     iSub=1 a collinear piece is subtracted from the qG->AAq amplitude
C            unconditionally
C     iSub=2 a collinear piece is subtracted from the qG->AAq amplitude
C            if the invariant mass of the collinear pair < cutM45 GeV
C     iSub=3 a collinear piece is subtracted from the qG->AAq amplitude
C            and a fragmenation piece is added
C     iSub=4 a collinear piece is subtracted from the qG->AAq amplitude
C            if the invariant mass of the collinear pair < cutM45 GeV
C            and a fragmenation piece is added
cpn2005 Negative values of iSub correspond to smooth cone isolation and
cpn2005 are  handled by PRIMARY_CUTS
C     iSub note: iSub=2 achieves Y = P - A ~ 0 at QT ~ 0 both at the 
C            Tevatron and the LHC.
cpn2005      If (iSub.NE.0) then

CCPY SEPT 2009
      IF(BOSON.EQ.'ZZ') ISUB=0
      
      If (iSub.gt.0) then
        Call SubtractColl (2, alfas, Amp2_Sub_up,Amp2_Sub_dn, j1)
        If (j1.Eq.1) Print*, Amp2_up, Amp2_up - Amp2_Sub_up,
     ,                           100.*(Amp2_up - Amp2_Sub_up)/Amp2_up
        Amp2_up = Amp2_up - Amp2_Sub_up
        Amp2_dn = Amp2_dn - Amp2_Sub_dn
      End If

CsB q/p G/pB + G/p q/pB cross section
      WT3 = qGUL*Amp2_up + qGDL*Amp2_dn
CsB !!! Temporary settings !!!
CsB !!! Only to test cancellation between Pert and Asympt !!!
C   Keep up-type quark only
C      WT3 = (SF(13))*Amp2_up
      COLORqG = 4.0/3.0/8.0
CsB The overall minus comes with the crossing
      WT3 = -WT3*SPIN*COLORqG*FACT_IDEN

      Do j = 1,4
        Plab(j,1) = -PStore(j,5)
        Plab(j,2) =  PStore(j,2)
        Plab(j,5) = -PStore(j,1)
      End Do

C Switching the indices 1 <-> 5 would be...
C      i1 = 5
C      i5 = 1
C      i2 = 2

      CALL Amp2qqBZZG (alfas, Amp2_up,Amp2_dn)
      If (Debug) then
        Print*, ' Amp2_up, Amp2_dn    = ', Amp2_up,Amp2_dn
      End If

CsB 9/9/97
C   For the q G -> photon photon q process we define a subtraction term which
C   subtracts the collinear singular term when one of the final state photon is
C   emitted from the final state quark leg.
C   This contribution will be included in the q G cross section which we
C   calculate using PYTHIA.
cpn2005      If (iSub.NE.0) then
      If (iSub.gt.0) then
        Call SubtractColl (1, alfas, Amp2_Sub_up,Amp2_Sub_dn,j1)
        Amp2_up = Amp2_up - Amp2_Sub_up
        Amp2_dn = Amp2_dn - Amp2_Sub_dn
      End If

CsB G/p qB/pB + qB/p G/pB cross section
      WT4 = GqUL*Amp2_up + GqDL*Amp2_dn
CsB The overall minus comes with the crossing
      WT4 = -WT4*SPIN*COLORqG*FACT_IDEN

CsB Restore pLab
C   Switch p1 <-> -p5
      Do J1=1,4
        PKEEP(J1)=PLAB(J1,1)
        PLAB(J1,1)=-PLAB(J1,5)
        PLAB(J1,5)=-PKEEP(J1)
      End Do

CsB Check pLab
      If (Debug) then
C        Print*, ' Checking momenta againts stored... '
        Do J1=1,4
          If (PLAB(J1,1).NE.PStore(j1,1)) Print*, ' Problem: PLAB 1 '
          If (PLAB(J1,2).NE.PStore(j1,2)) Print*, ' Problem: PLAB 2 '
          If (PLAB(J1,3).NE.PStore(j1,3)) Print*, ' Problem: PLAB 3 '
          If (PLAB(J1,4).NE.PStore(j1,4)) Print*, ' Problem: PLAB 4 '
          If (PLAB(J1,5).NE.PStore(j1,5)) Print*, ' Problem: PLAB 5 '
        End Do
C        Stop
      End If

CsB q qB: WT1
C   qB q: WT2
C   q  G: WT3
C   qB G: WT4
CsB iProc=0/1/2<->all/qqB/qG contribution
      If (iProc.Eq.0) then
        WT = WT1 + WT2 + WT3 + WT4   ! Both Contributions
      ElseIf (iProc.Eq.1) then
        WT = WT1 + WT2              ! q qB + qB q only
      ElseIf (iProc.Eq.2) then
        WT = WT3 + WT4              ! q G + qB G  only
      ElseIf (iProc.Eq.3) then
        WT = WT3                    ! q G  only
      elseif (iproc.eq.4) then      !qB G only
        WT = WT4
      Else
        Print*, ' iProc = ', iProc, ' value in not defined. '
        Stop
      EndIf
      If (Debug) Print*, ' WTi = ', WT1, WT2, WT3, WT4
C      If (Debug) Stop

      RETURN
      END

C --------------------------------------------------------------------------
      Subroutine SubtractColl (j5, alfas, Amp2Subup, Amp2Subdn, iDummy)
C --------------------------------------------------------------------------

CsB p(j5) is the momentum which p(5) is swithed with.
C   In order to be consistent with Amp2qqBZZG (to use the same constant
C   numerical factors, etc.) we use the scaleton of Amp2qqBZZG to construct the
C   subtraction piece.

      Implicit Double Precision (a-h,o-z)
      Dimension AmpSq(6,6), SP(5,5), t(5,5)
      COMMON/STAND0/XMW,XMZ,XMW2,XMZ2,AEM,GAMW,GAMZ,GAMW2,GAMZ2
      COMMON/STAND1/SWS,CWS,XMBOT,XMTOP,XMB2,XMT2,GWEAK2
      Common / ScalarProducts / SP
      Common / tij / t
      Common / Couplings / g_fVf_P, g_fVf_M
C     Common / Indices / i1,i2,i3,i4,i5
      COMMON/MOM/PLAB(4,12),PCEM(4,12),WHAT
CsB   In order to decide whether 3 or 4 is closer to 5 compare Delta R's.
C     At this stage of the code Primary Cuts was already called.
C     So we can get Delta R's from there.
      REAL*8 YXP,YXETA,YXPT,YXDELR,YXPP
      COMMON/SECCONV/ YXP(4,30),YXETA(30),YXPT(30),YXDELR(30,30),
     >YXPP(30)
CsB   We also need the 3&4 invariant mass (YM12)
      REAL*8 YM12,YPT12,YP12
      COMMON/MY1/YM12,YPT12,YP12(4)
      Integer iSub
      Real*8 cutM45, SubtractionShape
      logical SubQGFlag
      Common / Subtraction / cutM45, SubtractionShape(3:4),iSub,
     >  SubQGFlag(3:4)

      External AlphaS, ScaPro
      Logical Debug
      Data Debug / .False. /
      Data Pi / 3.141592653589793 /
C     Debug = .True.

CsB   p(j5) is the momentum which p(5) is swithed with.
C     If (Debug) Print*, j5
CsB   Obtain the momenta for the q G -> A A q process undoing the j5 <-> -5 cross:
      Do ii=1,4
        If (j5.Eq.1) then
          YXP(ii,1)=-PLAB(ii,5)
          YXP(ii,2)= PLAB(ii,2)
        Else If (j5.Eq.2.) then
          YXP(ii,1)= PLAB(ii,1)
          YXP(ii,2)=-PLAB(ii,5)
        Else
          Print*, ' Stopped in SubtractColl. '
          Stop
        End If
        YXP(ii,3)= PLAB(ii,3)
        YXP(ii,4)= PLAB(ii,4)
        YXP(ii,5)=-PLAB(ii,j5)
      End Do

CsB   Construct the necessary scalar products
C     Momenta (pLab) are passed through the Common Block in function ScaPro
C     pLab(i,j) is the i-th coordinate of the momentum p(j)
C     Evaluate scalar products once and for all
C     If (Debug) Print*, ' Scalar products: '
      Do ii = 1, 5
      Do jj = 1, ii
        SP(ii,jj) = ScaPro(ii,jj)
        SP(jj,ii) = SP(ii,jj)
C       If (Debug) Print*, ii, jj, SP(ii,jj)
      End Do
      End Do
CsB Undo the 1 <-> 5 or 2 <-> 5 cross:
      If (j5.Eq.1) then
        P12 = -SP(2,5)
        P13 = -SP(3,5)
        P14 = -SP(4,5)
      Else If (j5.Eq.2) then
        P12 = -SP(1,5)
        P13 =  SP(1,3)
        P14 =  SP(1,4)
      End If
      P35 = -SP(j5,3)
      P45 = -SP(j5,4)

CsB___Running alpha_ew implemented on 11/26/97
      eEM = Sqrt( 4.d0*Pi*PYALEM(2.d0*SP(3,4)) )

C     Print*, 2.d0*SP(3,4), eEM, gWeak2*SwS
CsB   We use the same prefactor here as for the amplitude in (D2)
C     PreFac = 16.d0*Pi*Alfas*(gWeak2*SwS)**2*(g_fVf_P**4+g_fVf_M**4)
      PreFac = 16.d0*Pi*Alfas*eEM**4*(g_fVf_P**4+g_fVf_M**4)

CsB The subtraction term is defined as the limit of the M^2(q G -> A A q) as
C   the photon which emitted from the final state quark and the quark becomes
C   collinear:
C   M^2_Sub = lim M^2(q G -> A A q) = 4 Pi alpha P(x,q <- q) M^2(q G -> A q)
C   (see Catani,Seymour hep-ph/9605323 and my notes including qG-AAq.nb)
C   with x = E(q)/(E(A)+E(q)) where p(q) and p(A) becomes parallel.
C   Note: either q (when calculating W3) or qBar (for W4) is labelled as 5.

CsB   In order to decide whether 3 or 4 is closer to 5 we compare the 3&5 and 4&5
C     invariant masses.
C     We construct the invariant masses below:
      YM35  = DSQRT( (YXP(4,3)+YXP(4,5))**2 - (YXP(3,3)+YXP(3,5))**2 -
     -               (YXP(2,3)+YXP(2,5))**2 - (YXP(1,3)+YXP(1,5))**2 )
      YM45  = DSQRT( (YXP(4,4)+YXP(4,5))**2 - (YXP(3,4)+YXP(3,5))**2 -
     -               (YXP(2,4)+YXP(2,5))**2 - (YXP(1,4)+YXP(1,5))**2 )
C     At this stage of the code Primary Cuts was already called.
C     We can get other variables from there if needed.
C     E.g. the 3&4 inavriant mass:
      YM34 = YM12

cCsB Also, construct the light-cone momentum components for 4 and 5 to define x:
cC   (We do not include Sqrt[2]'s because we use the ratios only.)
c      pp3  = YXP(4,3) + YXP(3,3)
c      pm3  = YXP(4,3) - YXP(3,3)
c      pp4  = YXP(4,4) + YXP(3,4)
c      pm4  = YXP(4,4) - YXP(3,4)
c      pp5  = YXP(4,5) + YXP(3,5)
c      pm5  = YXP(4,5) - YXP(3,5)
c      pp35 = pp3 + pp5
c      pm35 = pm3 + pm5
c      pp45 = pp4 + pp5
c      pm45 = pm4 + pm5

      If (YM45.LE.YM35) then ! 4 and 5 becomes collinear
        If ((iSub.Eq.2 .or. iSub.Eq.4) .and. YM45.GT.cutM45) then
          Amp2Subup = 0.d0
          Amp2Subdn = 0.d0
          Return
        End If

CsB     Covariant definition of x (after Catani-Seymour):
        x = YM35/(YM34 + YM35 + YM45)
        IF(x.LE.0.d0) x=1.d-6
        IF(x.GE.1.d0) x=.999999d0

CsB     Note that compared to my notes there's an overall "-" sign difference with
C       the expression below. This is because in my notes I already include the "-"
C       which comes from the crossing. On the other hand in this code the "-" sign
C       is applied after this routine is called. So the "-" in my notes has to be
C       taken out.
        AmpS = -4.d0*(1.d0+x**2)/(1.d0-x)*(P12**2+P13**2)/P12/P13/P45

CsB     The fragmentation piece is defined as the subtraction piece, except
C       P[q->q photon] replaced by D^{LL}[q->q photon].
C       (See Bailey, Owens and Ohnemus, PRD46,2026(1992).
C       Note: their z is our 1 - x. Also, we use only the non-singlet term here.
C       The difference is small for small z. See my Fragmentation.nb for further
C       details.)
C       The non-singlet piece:
        z = 1.d0 - x
        AmpFNS = -4.d0*
     *                 (2.21d0 - 1.28d0*z + 1.29d0*z**2) * z**0.049d0/
     /                 (z * (1.d0-1.63d0*Log(x)))*
     *           (P12**2+P13**2)/P12/P13/P45
C       The singlet piece:
        AmpFS  = -4.d0*
     *                 0.0020* x**2 *z**(-1.54d0)*
     *           (P12**2+P13**2)/P12/P13/P45

      Else If (YM45.GT.YM35) then ! 3 and 5 becomes collinear
        If ((iSub.Eq.2 .or. iSub.Eq.4) .and. YM35.GT.5.d0) then
          Amp2Subup = 0.d0
          Amp2Subdn = 0.d0
          Return
        End If

CsB     Covariant definition of x (after Catani-Seymour):
        x = YM45/(YM34 + YM35 + YM45)
        IF(x.LE.0.d0) x=1.d-6
        IF(x.GE.1.d0) x=.999999d0
        AmpS   = -4.d0*(1.d0+x**2)/(1.d0-x)*(P12**2+P14**2)/P12/P14/P35
        z = 1.d0 - x
        AmpFNS = -4.d0*
     *                 (2.21d0 - 1.28d0*z + 1.29d0*z**2) * z**0.049d0/
     /                 (z * (1.d0-1.63d0*Log(x)))*
     *           (P12**2+P14**2)/P12/P14/P35
        AmpFS  = -4.d0*
     *                 0.0020* x**2 *z**(-1.54d0)*
     *           (P12**2+P14**2)/P12/P14/P35
      End If

C For up type quark:
C ((D2) contains the color sum factor 4.0. We will put it in later)

      If (iSub.Eq.3 .or. iSub.Eq.4) then ! Subtraction is done AND frag piece is added
CsB___Subtract and add back fragmentation:
C     (Remember:
C      1. the term below will be SUBTRACTED from the full 2 -> 3 amp^2
C      2. PreFac contains Q_u^4 )
        Amp2Subup = PreFac*(AmpS - AmpFNS - AmpFS/(2.d0/3.d0)**2)/
     /              4.d0
        Amp2Subdn = PreFac*(AmpS - AmpFNS - AmpFS/(-1.d0/3.d0)**2)/
     /              4.d0/16.d0
        Return
      End If

CsB___Only subtraction is done:
      Amp2Subup = PreFac*AmpS/4.d0

C g_fVf_P and g_fVf_M are defined for up type quark: 2/3 and 0
C (for double photon production).
C For down type quark: Amp2d = (Q_d/Q_u)^4 * Amp2up
      Amp2Subdn = Amp2Subup/16.d0

CsB -----
CsB This part below is just tests the collinear pieces. It selects kinematics
C   when M35 or M45 is small and prints M2(2->3) and |M(Subtraction)|^2.
C   According to my tests the two terms agree within 1%, which numerically
C   proves that the subtraction terms are correct.
      If (Debug) then
        iDummy = 0
        If (YM35.LT..1) then
          Print*, ' YM35,Amp2Subup:', YM35, Amp2Subup
          iDummy = 1
        End If
        If (YM45.LT..1) then
          Print*, ' YM45,Amp2Subup:', YM45, Amp2Subup
          iDummy = 1
        End If
      End If

      Return
      End ! SubtractColl

cpn2005-----------------------------------------------------------------
      function ConeETmax(Ephys,delr)
cpn2005 The function returns the maximal energy of the unobserved jet
cpn2005 allowed by the smooth cone isolation criterion a la Frixione
cpn2005 Input: isub -- defines the type of the isolation function 
cpn2005        Ephys -- physical energy that controls the isolation
cpn2005 condition, for instance the transverse energy of the photon 
cpn2005 or the diphoton's transverse momentum
cpn2005        delr -- angular separation of the photon and the jet
      implicit none
      real*8 ConeETmax, Ephys, delr, eps_c, ChiDelta
      Integer iSub
      Real*8 cutM45, SubtractionShape
      logical SubQGFlag
      Common / Subtraction / cutM45, SubtractionShape(3:4),iSub,
     >  SubQGFlag(3:4)

cpn2005 Isolation cone parameters
      integer iqqt
      real*8 DelRIso, ETmaxIso
      common/isolation/DelRIso, ETMaxIso,iqqt

      data eps_c /1.0/      

cpn2005 This is Frixione's Chi(delta) function
c      ChiDelta(delr)= Ephys*eps_c*
c     >  (1d0-cos(delr))/(1d0-cos(DelRIso))
cpn   Modified Frixione's Chi(delta) v.2
      ChiDelta(delr)= ETmaxIso*
     >  ((1d0-cos(delr))/(1d0-cos(DelRIso)))**2

      
      
      if (isub.eq.-1) then      !smooth cone isolation
        ConeETmax = ChiDelta(delr)
      else if (isub.eq.-2) then !smooth + conventional isolation
        ConeETmax = min(ChiDelta(delr),ETmaxIso)
      else                      !conventional isolation cone
        ConeETmax = ETmaxIso
      endif                     !isub

      return 
      end                       !ConeETmax

 
      
CsB ------------------------------------------------------------------------

C Calculating the squared matrix element of the real emission process:
C    q(p(1)) qB(p(2)) -> Z(p(3)) Z(p(4)) G(p(5)).
C Based on PRD43('91)3637 (Ohnemus-Owens).

C --------------------------------------------------------------------------
      Subroutine Amp2qqBZZG (alfas, Amp2up, Amp2dn)
C --------------------------------------------------------------------------
      Implicit Double Precision (a-h,o-z)

      Dimension AmpSq(6,6), SP(5,5), t(5,5)
      COMMON/STAND0/XMW,XMZ,XMW2,XMZ2,AEM,GAMW,GAMZ,GAMW2,GAMZ2
      COMMON/STAND1/SWS,CWS,XMBOT,XMTOP,XMB2,XMT2,GWEAK2
      Common / ScalarProducts / SP
      Common / tij / t
      Common / Couplings / g_fVf_P, g_fVf_M
c      Common / VVMass / xMV
      Common / Indices / i1,i2,i3,i4,i5
      External AlphaS, ScaPro
      Character*40 Boson
      Common / Boson / Boson
c      Real*8 SW2_EFF
c      Common/Stand6/SW2_EFF
c      Real*8 g_fVf_P_u,g_fVf_M_u,g_fVf_P_d,g_fVf_M_d
      Common / Stand6 / g_fVf_P_u,g_fVf_M_u,g_fVf_P_d,g_fVf_M_d
      Data Pi / 3.141592653589793 /
      Logical Debug
      Data Debug / .False. /

C      Debug = .True.

Csb Check indices
C      Print*, ' i1,i2,i3,i4,i5: ', i1,i2,i3,i4,i5

C Momenta (pLab) are passed through the Common Block in function ScaPro
C pLab(i,j) is the i-th coordinate of the momentum p(j)
C Evaluate scalar products once and for all
C      If (Debug) Print*, ' Scalar products: '
      Do 100 ii = 1, 5
      Do 100 jj = 1, ii
        SP(ii,jj) = ScaPro(ii,jj)
        SP(jj,ii) = SP(ii,jj)
C        If (Debug) Print*, ii, jj, SP(ii,jj)
 100  Continue

C Mandelstam variables: t_ij = (p_i - p_j)^2
C Checked against C.-P.'s propagators
C The agreement also suggests that the scalar products are right.
C      If (Debug) Print*, ' Mandelstam variables: '
      Do 130 ii = 1, 5
      Do 130 jj = 1, ii
        t(ii,jj) = SP(ii,ii) - 2.d0*SP(ii,jj) + SP(jj,jj)
        t(jj,ii) = t(ii,jj)
C        If (Debug) Print*, ii, jj, t(ii,jj)
 130  Continue

C Relations (D4)
      AmpSq(1,1) = A2_11(i1,i2,i3,i4,i5)
      AmpSq(1,2) = A2_12(i1,i2,i3,i4,i5)
      AmpSq(1,3) = A2_13(i1,i2,i3,i4,i5)
      AmpSq(1,4) = A2_14(i1,i2,i3,i4,i5)
      AmpSq(1,5) = A2_15(i1,i2,i3,i4,i5)
      AmpSq(1,6) = A2_16(i1,i2,i3,i4,i5)

      AmpSq(2,2) = A2_11(i1,i2,i4,i3,i5)
      AmpSq(2,3) = A2_14(i1,i2,i4,i3,i5)
      AmpSq(2,4) = A2_13(i1,i2,i4,i3,i5)
      AmpSq(2,5) = A2_16(i1,i2,i4,i3,i5)
      AmpSq(2,6) = A2_15(i1,i2,i4,i3,i5)

      AmpSq(3,3) = A2_33(i1,i2,i3,i4,i5)
      AmpSq(3,4) = A2_12(i2,i1,i3,i4,i5)
      AmpSq(3,5) = A2_15(i2,i1,i4,i3,i5)
      AmpSq(3,6) = A2_16(i2,i1,i4,i3,i5)

      AmpSq(4,4) = A2_33(i1,i2,i4,i3,i5)
      AmpSq(4,5) = A2_16(i2,i1,i3,i4,i5)
      AmpSq(4,6) = A2_15(i2,i1,i3,i4,i5)

      AmpSq(5,5) = A2_55(i1,i2,i3,i4,i5)
      AmpSq(5,6) = A2_14(i1,i2,i5,i4,i3)

      AmpSq(6,6) = A2_55(i1,i2,i4,i3,i5)

CsB___Running alpha_ew implemented on 11/26/97
      eEM = Sqrt( 4.d0*Pi*PYALEM(2.d0*SP(3,4)) )

C      PreFac = 16.d0*Pi*Alfas*(gWeak2*SwS)**2*(g_fVf_P**4+g_fVf_M**4)
      PreFac = 16.d0*Pi*Alfas*eEM**4*(g_fVf_P**4+g_fVf_M**4)

      Amp2 = 0.d0
C_____Sum of the squared amplitudes
      If (Debug) Print*, ' Matrix elements (up type): '
      Do 120 ii = 1,6
      Do 120 jj = ii,6
C          If (Debug) Print*, ii,jj, AmpSq(ii,jj)*PreFac
          If (Debug .and. ii.Eq.jj) Print*, ii,jj, AmpSq(ii,jj)*PreFac
          Amp2 = Amp2 + AmpSq(ii,jj)
 120  Continue

C Fromula (D2)
C For up type quark:
C ((D2) contains the color sum factor 4.0. We will put it in later)
      Amp2up = PreFac*Amp2/4.d0

      If (Boson.Eq.'AA') then
C_______g_fVf_P and g_fVf_M are defined for up type quark: 2/3 and 0
C       (for double photon production).
C       For down type quark: Amp2d = (Q_d/Q_u)^4 * Amp2up
        Amp2dn = Amp2up/16.d0
      Else If (Boson.Eq.'ZZ') then
C_______For d type quark: Amp2d = ((gM_d^4+gP_d^4)/(gM_u^4+gP_u^2)) * Amp2up
        Amp2dn = (g_fVf_M_d**4 + g_fVf_P_d**4)/
     /           (g_fVf_M_u**4 + g_fVf_P_u**4) * Amp2up
      End If

      Return
      End

C --------------------------------------------------------------------------
      Function A2R (i1,i2,i3,i4,i5)
C --------------------------------------------------------------------------
C Expression (D3)
      Implicit Double Precision (a-h,o-z)
      Dimension SP(5,5), t(5,5)
      Real*8 MV
      Common / ScalarProducts / SP
      Common / tij / t
      Common / VVMass / MV
      Logical Debug
      Data Debug / .False. /

      Entry A2_11 (i1,i2,i3,i4,i5)
C      Debug = .True.
C      If (Debug) Print*, ' MV = ', MV
C      If (Debug) Print*, ' tij: ', t(i1,i4),t(i1,i5),t(i2,i3),t(i2,i4)
C Based on the scalar products checked this numerical result by hand.
C Seems to be fine.
      A2_11 = 32.d0*SP(i1,i5)*
     *        (2.d0*SP(i2,i4)*SP(i4,i5)-SP(i2,i5)*MV**2)/
     /        ( t(i1,i5)**2*t(i2,i4)**2 )
C      If (Debug) Print*, ' A2_11 = ', A2_11
C      Debug = .False.
      Return

      Entry A2_33 (i1,i2,i3,i4,i5)
c      Debug = .True.
      A2_33 = 32.d0*SP(i2,i5)*
     *        (2.d0*SP(i1,i3)*SP(i3,i5)-SP(i1,i5)*MV**2)/
     /        ( t(i1,i3)**2*t(i2,i5)**2 )
c      If (Debug) Print*, ' A2_33 = ', A2_33
c      Debug = .False.
      Return

      Entry A2_55 (i1,i2,i3,i4,i5)
c      Debug = .True.
      A2_55 = 16.d0*( SP(i1,i2)*MV**4 -
     -                2.d0*SP(i1,i3)*SP(i2,i3)*MV**2 +
     +                4.d0*SP(i1,i3)*SP(i2,i4)*SP(i3,i4) -
     -                2.d0*SP(i1,i4)*SP(i2,i4)*MV**2 )/
     /        ( t(i1,i3)**2*t(i2,i4)**2 )
c      If (Debug) Print*, ' A2_55 = ', A2_55
c      Debug = .False.
      Return

      Entry A2_12 (i1,i2,i3,i4,i5)
c      Debug = .True.
      A2_12 = 128.d0*SP(i1,i5)*SP(i2,i5)*
     *        ( SP(i2,i3) + SP(i2,i4) - SP(i3,i4) )/
     /        ( t(i1,i5)**2*t(i2,i3)*t(i2,i4) )
c      If (Debug) Print*, ' A2_12 = ', A2_12
c      Debug = .False.
      Return

      Entry A2_13 (i1,i2,i3,i4,i5)
c      Debug = .True.
      A2_13 = 64.d0*
     *        ( SP(i1,i2)**2*(SP(i3,i4) - SP(i3,i5) - SP(i4,i5)) +
     +          SP(i1,i2)*(SP(i1,i4)*SP(i3,i5) + SP(i1,i3)*SP(i2,i5) -
     -                   SP(i1,i4)*SP(i2,i3) + SP(i1,i4)*SP(i2,i5) -
     -                   SP(i1,i3)*SP(i2,i4)) +
     +          SP(i1,i2)*(SP(i1,i5)*SP(i2,i3) + SP(i1,i5)*SP(i2,i4) -
     -                   SP(i1,i5)*SP(i3,i4) + SP(i2,i3)*SP(i4,i5) -
     -                   SP(i2,i5)*SP(i3,i4)) +
     +          SP(i1,i3)*SP(i2,i5)*(SP(i2,i4) - SP(i1,i4)) +
     +          SP(i1,i5)*SP(i2,i4)*(SP(i1,i3) - SP(i2,i3)) )/
     /        ( t(i1,i3)*t(i1,i5)*t(i2,i4)*t(i2,i5) )
c      If (Debug) Print*, ' A2_13 = ', A2_13
c      Debug = .False.
      Return

      Entry A2_14 (i1,i2,i3,i4,i5)
c      Debug = .True.
      A2_14 = 128.d0*SP(i1,i2)*
     *        (SP(i1,i2) - SP(i1,i4) - SP(i2,i4) + SP(i4,i4))*
     *        (SP(i1,i2) - SP(i1,i5) - SP(i2,i5) + SP(i5,i5))/
     /        ( t(i1,i4)*t(i1,i5)*t(i2,i4)*t(i2,i5) )
c      If (Debug) Print*, ' A2_14 = ', A2_14
c      Debug = .False.
      Return

      Entry A2_15 (i1,i2,i3,i4,i5)
c      Debug = .True.
      A2_15 = 64.d0*
     *        (SP(i1,i3) + SP(i1,i5) - SP(i3,i5))*
     *        (2.d0*SP(i1,i4)*SP(i2,i4) - SP(i1,i2)*MV**2)/
     /        ( t(i1,i3)*t(i1,i5)*t(i2,i4)**2 )
c      If (Debug) Print*, ' A2_15 = ', A2_15
c      Debug = .False.
      Return

      Entry A2_16 (i1,i2,i3,i4,i5)
c      Debug = .True.

      A2_16 = 64.d0*
     *        ( SP(i1,i2)**2*(SP(i3,i5) - SP(i3,i4) - SP(i4,i5)) +
     +          SP(i1,i2)*(SP(i1,i3)*SP(i2,i4) - SP(i1,i3)*SP(i2,i5) +
     +                     SP(i1,i3)*SP(i4,i5) + SP(i1,i4)*SP(i2,i3) +
     +                     SP(i1,i4)*SP(i2,i5)) +
     +          SP(i1,i2)*(SP(i1,i5)*SP(i2,i4) - SP(i1,i5)*SP(i2,i3) -
     -                     SP(i1,i4)*SP(i3,i5) - SP(i2,i4)*SP(i3,i5) +
     +                     SP(i2,i5)*SP(i3,i4)) +
     +          SP(i1,i5)*SP(i2,i4)*(SP(i2,i3) - SP(i1,i3)) +
     +          SP(i1,i4)*SP(i2,i3)*(SP(i1,i5) - SP(i2,i5)) +
     +          SP(i1,i2)*MV**2*(SP(i1,i2) - SP(i1,i3) - SP(i2,i5)) +
     +          MV**2/2.d0*(SP(i1,i2)*SP(i3,i5) + SP(i1,i3)*SP(i2,i5) -
     -                      SP(i1,i5)*SP(i2,i3)) )/
     /        ( t(i1,i4)*t(i1,i5)*t(i2,i3)*t(i2,i4) )
c      If (Debug) Print*, ' A2_16 = ', A2_16
c      Debug = .False.
      Return

      End

C --------------------------------------------------------------------------
      Function ScaPro(i,j)
C --------------------------------------------------------------------------
      Implicit Double Precision (a-h,o-z)
      COMMON/MOM/PLAB(4,12),PCEM(4,12),WHAT

C Using the metric: (-,-,-,+)
      ScaPro = -pLab(1,i)*pLab(1,j) - pLab(2,i)*pLab(2,j) -
     -          pLab(3,i)*pLab(3,j) + pLab(4,i)*pLab(4,j)

      Return
      End

C --------------------------------------------------------------------------
      Function ScaProS(i,j)
C --------------------------------------------------------------------------
      Implicit Double Precision (a-h,o-z)
      COMMON/DECAY_MOMENTA/SEC_PLAB(4,50)

C Using the metric: (-,-,-,+)
      ScaProS = -Sec_pLab(1,i)*Sec_pLab(1,j) - 
     -           Sec_pLab(2,i)*Sec_pLab(2,j) -
     -           Sec_pLab(3,i)*Sec_pLab(3,j) + 
     +           Sec_pLab(4,i)*Sec_pLab(4,j)
     
      Return
      End

CsB =<======================================================================

      SUBROUTINE YPARTONS
cpn   Return parton densities Q1, Q2 for given x1, x2 values in given hadrons 1 and 2 (protons, antiprotons, nuclei) 
      IMPLICIT NONE
      integer nflmax, ipdf, irt
      parameter (nflmax=5)
      Integer BIPDF, ISet
      double precision Q1, Q2, Q2tmp(-nflmax:nflmax), Frac_n, Frac_p, 
     >  X1, X2, Q, Ctq4Pdf, Ctq5Pdf, Ctq6Pdf
      COMMON/STRUCTURE/Q1(-nflmax:nflmax),Q2(-nflmax:nflmax)
      COMMON/PMF/X1,X2
      COMMON/SCALE/Q
     
      Common / BuildInPDF / BIPDF, ISet

      REAL*8 W,QCDL,X1X2MIN
      INTEGER*4 NJET,N_PRIMARY_PART,NF,I_STR,IBEAM,NPART,NEVENT,NORDER,
     &     IW,N_SECOND_PART
      COMMON/AA/W,QCDL,X1X2MIN,NJET,N_PRIMARY_PART,NF,I_STR,IBEAM,NPART,
     &  NEVENT,NORDER,N_SECOND_PART

      data frac_n /0.55555555/


      If (BIPDF.GE.31 .and. BIPDF.LE.39) then
        call ctq3pds(iset,q1,x1,Q,irt)
        call ctq3pds(iset,q2,x2,Q,irt)
      ElseIf (BIPDF.GE.41 .and. BIPDF.LE.49) then
C   The functions Ctq4Pdf, Ctq5Pdf, Ctq6Pdf
C   returns the parton distribution inside the proton for parton [Iparton]
C   at [X] Bjorken_X and scale [Q] (GeV) in PDF set [Iset].
C   Iparton  is the parton label (5, 4, 3, 2, 1, 0, -1, ......, -5)
C                            for (b, c, s, d, u, g, u_bar, ..., b_bar)
        do ipdf = -nflmax, nflmax
          q1(ipdf) = x1*Ctq4Pdf ( ipdf, X1, Q)
          q2(ipdf) = x2*Ctq4Pdf ( ipdf, X2, Q)
        enddo 

      ElseIf (BIPDF.GE.51 .and. BIPDF.LE.59) then
C   The function Ctq5Pdf (Iparton, X, Q)
C   returns the parton distribution inside the proton for parton [Iparton] 
C   at [X] Bjorken_X and scale [Q] (GeV) in PDF set [Iset].
C   Iparton  is the parton label (5, 4, 3, 2, 1, 0, -1, ......, -5)
C                            for (b, c, s, d, u, g, u_bar, ..., b_bar),
        do ipdf = -nflmax, nflmax
          q1(ipdf) = x1*Ctq5Pdf ( ipdf, X1, Q)
          q2(ipdf) = x2*Ctq5Pdf ( ipdf, X2, Q)
        enddo 
CCPY      ElseIf (BIPDF.GE.61 .and. BIPDF.LE.69) then
        Else If (BIPDF.GE.200 .AND. BIPDF.LT.500) then ! CTEQ6
        do ipdf = -nflmax, nflmax
          q1(ipdf) = x1*Ctq6Pdf ( ipdf, X1, Q)
          q2(ipdf) = x2*Ctq6Pdf ( ipdf, X2, Q)
        enddo 
CCPY
        Else If (BIPDF.EQ.903) then ! USING pdf.pds file
        do ipdf = -nflmax, nflmax
          q1(ipdf) = x1*Ctq6Pdf ( ipdf, X1, Q)
          q2(ipdf) = x2*Ctq6Pdf ( ipdf, X2, Q)
        enddo 
      Else
        Print*, ' YPartons: unknown PDF set for BIPDF = ', BIPDF
        Stop
      End If

cpn   By now, Q1 and Q2 contain PDF's for _protons_ at momentum 
cpn   fractions x1 and x2. Assuming the hadron 1 to be a proton, 
cpn   we keep Q1 as it is. We recompute Q2 according to the type 
cpn   of the hadron 2 
      
      if (ibeam.eq.-1) then     !ppbar collision
        do ipdf=-nflmax, nflmax
          q2tmp(ipdf) = q2(ipdf)
        enddo     

        do ipdf=-nflmax, nflmax
          q2(ipdf) = q2tmp(-ipdf) 
        enddo 
      elseif(ibeam.eq.0) then   !p-nucleus collision

        frac_p = 1d0 - frac_n

        q2tmp(1) = frac_p*q2(1) + frac_n*q2(2)
        q2tmp(2) = frac_p*q2(2) + frac_n*q2(1)
        q2(1) = q2tmp(1); q2(2) = q2tmp(2)

cpn2006 It is not at all obvious that u and d sea quark PDFs are
c        interchanged in protons and neutrons, especially when they are a
c       part of a nucleus; here it is assumed that u and d sea quark PDFs
c       are switched in accordance with the original ResBos convention
        q2tmp(-1) = frac_p*q2(-1) + frac_n*q2(-2)
        q2tmp(-2) = frac_p*q2(-2) + frac_n*q2(-1)
        q2(-1) = q2tmp(-1); q2(-2) = q2tmp(-2)
      endif                     !ibeam

      RETURN
      END                       !YPARTONS

C =========================================================================
C =========================================================================
CsB   Code below is the CTEQ3 PDF. Needed for processes with no input grids.

C     Version 3 CTEQ distribution function in a parametrized form.

C   By: H.L. Lai, J. Botts, J. Huston, J.G. Morfin, J.F. Owens, J. Qiu,
C       W.K. Tung & H. Weerts;  Preprint MSU-HEP/41024, CTEQ 404

C   This file contains three versions of the same CTEQ3 parton distributions:
C
C Two "front-end" subprograms:
C     FUNCTION Ctq3Pd (Iset, Iparton, X, Q, Irt)
C         returns the PROBABILITY density for a GIVEN flavor;
C     SUBROUTINE Ctq3Pds(Iset, Pdf, XX, QQ, Irt)
C         returns an array of MOMENTUM densities for ALL flavors;
C One lower-level subprogram:
C     FUNCTION Ctq3df (Iset, Iprtn, XX, QQ, Irt)
C         returns the MOMENTUM density of a GIVEN valence or sea distribution.

C      One supplementary function to return the QCD lambda parameter
C      concerning these distributions is also included (see below).

C     Although DOUBLE PRECISION is used, conversion to SINGLE PRECISION
C     is straightforward by removing the
C     Implicit Double Precision statements.

C     Since this is an initial distribution of version 3, it is
C     useful for the authors to maintain a record of the distribution
C     list in case there are revisions or corrections.
C     In the interest of maintaining the integrity of this package,
C     please do not freely distribute this program package; instead, refer
C     any interested colleagues to direct their request for a copy to:
C     Lai@cteq11.pa.msu.edu or Tung@msupa.pa.msu.edu.

C   If you have detailed questions concerning these CTEQ3 distributions,
C   or if you find problems/bugs using this initial distribution, direct
C   inquires to Hung-Liang Lai or Wu-Ki Tung.

C     -------------------------------------------
C     Detailed instructions follow.

C     Name convention for CTEQ distributions:  CTEQnSx  where
C        n : version number                      (currently n = 3)
C        S : factorization scheme label: = [M L D] for [MS-bar LO DIS]
c               resp.
C        x : special characteristics, if any
C        (e.g. S(F) for singular (flat) small-x, L for "LEP lambda value")
C        (not applicable to CTEQ3 since only three standard sets are given.)

C    Explanation of functional arguments:

C    Iset is the set label; in this version, Iset = 1, 2, 3
C                           correspond to the following CTEQ global fits:

C          cteq3M  : best fit in the MS-bar scheme
C          cteq3L  : best fit in Leading order QCD
C          cteq3D  : best fit in the DIS scheme

C   Iprtn  is the parton label (6, 5, 4, 3, 2, 1, 0, -1, ......, -6)
C                          for (t, b, c, s, d, u, g, u_bar, ..., t_bar)
C  *** WARNING: We use the parton label 2 as D-quark, and 1 as U-quark which
C               might be different with your labels.

C   X, Q are the usual x, Q;
C   Irt is a return error code (see individual modules for explanation).
C
C     ---------------------------------------------

C  Since the QCD Lambda value for the various sets are needed more often than
C  the other parameters in most applications, a special function
C     Wlamd3 (Iset, Iorder, Neff)                    is provided
C  which returns the lambda value for Neff = 4,5,6 effective flavors as well as
C  the order these values pertain to.

C     ----------------------------------------------
C     The range of (x, Q) used in this round of global analysis is, approxi-
C     mately,  0.01 < x < 0.75 ; and 4 GeV^2 < Q^2 < 400 GeV^2 for fixed target
C     experiments and 0.0001 < x < 0.1 from HERA data.

C    The range of (x, Q) used in the reparametrization of the QCD evolved
C    parton distributions is 10E-6 < x < 1 ; 1.6 GeV < Q < 10 TeV.  The
C    functional form of this parametrization is:

C      A0 * x^A1 * (1-x)^A2 * (1 + A3 * x^A4) * [log(1+1/x)]^A5

C   with the A'coefficients being smooth functions of Q.  For heavy quarks,
C   a threshold factor is applied to A0 which simulates the proper Q-dependence
C   of the QCD evolution in that region according to the renormalization
C   scheme defined in Collins-Tung, Nucl. Phys. B278, 934 (1986).

C   Since this function is positive definite and smooth, it provides sensible
C   extrapolations of the parton distributions if they are called beyond
C   the original range in an application. There is no artificial boundaries
C   or sharp cutoff's.
C    ------------------------------------------------

      FUNCTION Ctq3Pd (Iset, Iparton, X, Q, Irt)

C   This function returns the CTEQ parton distributions f^Iset_Iprtn/proton
C   --- the PROBABILITY density

C   (Iset, Iparton, X, Q): explained above;

C    Irt : return error code: see module Ctq3df for explanation.

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      Ifl = Iparton
      JFL = ABS(Ifl)
C                                                             Valence
      IF (Ifl.Eq.1 .or. Ifl.Eq.2) THEN
        VL = Ctq3df(Iset, Ifl, X, Q, Irt)
      ELSE
        VL = 0.
      ENDIF
C                                                             Sea
      SEA = Ctq3df (Iset, -JFL, X, Q, Irt)
C                                              Full (probability) Distribution
      Ctq3pd = (VL + SEA) / X

      Return
C                         *************************
      END
C
C
      SUBROUTINE Ctq3Pds(Iset, Pdf, X, Q, Irt)

C   This function returns the CTEQ parton distributions xf^Iset_Iprtn/proton
C   --- the Momentum density in array form
c
C    (Iset, X, Q): explained in header comment lines;

C     Irt : return error code -- cumulated over flavors:
C           see module Ctq3df for explanation on individual flavors.
C     Pdf (Iparton);
C         Iparton = -6, -5, ...0, 1, 2 ... 6
C               has the same meaning as explained in the header comment lines.

      Implicit Double Precision (A-H, O-Z)
      Dimension Pdf (-5:5)

      Irt=0
      do 10 I=-5,2
         if(I.le.0) then
            Pdf(I) = Ctq3df(Iset,I,X,Q,Irt1)
            Pdf(-I)= Pdf(I)
         else
            Pdf(I) = Ctq3df(Iset,I,X,Q,Irt1) + Pdf(-I)
         endif
         Irt=Irt+Irt1
  10  Continue

      Return
C                         *************************
      End

      FUNCTION Ctq3df (Iset, Iprtn, XX, QQ, Irt)

C            Returns xf(x,Q) -- the momentum fraction distribution !!
C            Returns valence and sea rather than combined flavor distr.

C            Iset : PDF set label

C            Iprtn  : Parton label:   2, 1 = d_ and u_ valence
C                                     0 = gluon
C                            -1, ... -6 = u, d, s, c, b, t sea quarks

C            XX  : Bjorken-x
C            QQ  : scale parameter "Q"
C      Irt : Return code
C      0 : no error
C      1 : parametrization is slightly negative; reset to 0.0.
C          (This condition happens rarely -- only for large x where the
C          absolute value of the parton distribution is extremely small.)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)
      Parameter (Nst = 3)

      DIMENSION
     >   Iord(Nst), Isch(Nst), Nqrk(Nst),Alm(Nst)
     > , Vlm(4:6,Nst), Qms(4:6, Nst)
     > , Xmn(Nst), Qmn(Nst), Qmx(Nst)

c                                          --------- CTEQ3M
c
      DATA
     >  Isch(1), Iord(1), Nqrk(1), Alm(1) /  1,  2,  6, .239  /
     >  (Vlm(I,1), I=4,6) / .239,    .158,     .063   /
     >  (Qms(I,1), I=4,6) / 1.60,   5.00,  180.0 /
     >  Xmn(1), Qmn(1), Qmx(1) /  1.E-6,  1.60,  1.E4  /

c                                          --------- CTEQ3L
c
      DATA
     >  Isch(2), Iord(2), Nqrk(2), Alm(2) /  1,  1,  6, .177  /
     >  (Vlm(I,2), I=4,6) / .177,    .132,     .066   /
     >  (Qms(I,2), I=4,6) / 1.60,   5.00,  180.0 /
     >  Xmn(2), Qmn(2), Qmx(2) /  1.E-6,  1.60,  1.E4  /

c                                          --------- CTEQ3D
c
      DATA
     >  Isch(3), Iord(3), Nqrk(3), Alm(3) /  1,  2,  6, .247  /
     >  (Vlm(I,3), I=4,6) / .247,    .164,     .066   /
     >  (Qms(I,3), I=4,6) / 1.60,   5.00,  180.0 /
     >  Xmn(3), Qmn(3), Qmx(3) /  1.E-6,  1.60,  1.E4  /


      Data Ist, Lp, Qsto / 0, -10, 1.2345 /

      save Ist, Lp, Qsto
      save SB, SB2, SB3

      X  = XX
      Irt = 0
      if(Iset.eq.Ist .and. Qsto.eq.QQ) then
C                                             if only change is in x:
        if (Iprtn.eq.Lp) goto 100
C                         if change in flv is within "light" partons:
        if (Iprtn.ge.-3 .and. Lp.ge.-3) goto 501
      endif

      Ip = abs(Iprtn)
C                                                  Set up Qi for SB
      If (Ip .GE. 4) then
         If (QQ .LE. Qms(Ip, Iset)) Then
           Ctq3df = 0.0
           Return
         Endif
         Qi = Qms(ip, Iset)
      Else
         Qi = Qmn(Iset)
      Endif
C                   Use "standard lambda" of parametrization program
      Alam = Alm (Iset)

      SBL = LOG(QQ/Alam) / LOG(Qi/Alam)
      SB = LOG (SBL)
      SB2 = SB*SB
      SB3 = SB2*SB

 501  Iflv = 3 - Iprtn

      Goto (1,2,3, 311) Iset

 1    Goto(11,12,13,14,15,16,17,18,19)Iflv
c   Ifl =   2
  11  A0=Exp(-0.7266E+00-0.1584E+01*SB +0.1259E+01*SB2-0.4305E-01*SB3)
      A1= 0.5285E+00-0.3721E+00*SB +0.5150E+00*SB2-0.1697E+00*SB3
      A2= 0.4075E+01+0.8282E+00*SB -0.4496E+00*SB2+0.2107E+00*SB3
      A3= 0.3279E+01+0.5066E+01*SB -0.9134E+01*SB2+0.2897E+01*SB3
      A4= 0.4399E+00-0.5888E+00*SB +0.4802E+00*SB2-0.1664E+00*SB3
      A5= 0.3678E+00-0.8929E+00*SB +0.1592E+01*SB2-0.5713E+00*SB3
      goto 100
c   Ifl =   1
  12  A0=Exp( 0.2259E+00+0.1237E+00*SB +0.3035E+00*SB2-0.2935E+00*SB3)
      A1= 0.5085E+00+0.1651E-01*SB -0.3592E-01*SB2+0.2782E-01*SB3
      A2= 0.3732E+01+0.4901E+00*SB +0.2218E+00*SB2-0.1116E+00*SB3
      A3= 0.7011E+01-0.6620E+01*SB +0.2557E+01*SB2-0.1360E+00*SB3
      A4= 0.8969E+00-0.2429E+00*SB +0.1811E+00*SB2-0.6888E-01*SB3
      A5= 0.8636E-01+0.2558E+00*SB -0.3082E+00*SB2+0.2535E+00*SB3
      goto 100
c   Ifl =   0
  13  A0=Exp(-0.2318E+00-0.9779E+00*SB -0.3783E+00*SB2+0.1037E-01*SB3)
      A1=-0.2916E+00+0.1754E+00*SB -0.1884E+00*SB2+0.6116E-01*SB3
      A2= 0.5349E+01+0.7460E+00*SB +0.2319E+00*SB2-0.2622E+00*SB3
      A3= 0.6920E+01-0.3454E+01*SB +0.2027E+01*SB2-0.7626E+00*SB3
      A4= 0.1013E+01+0.1423E+00*SB -0.1798E+00*SB2+0.1872E-01*SB3
      A5=-0.5465E-01+0.2303E+01*SB -0.9584E+00*SB2+0.3098E+00*SB3
      goto 100
c   Ifl =  -1
  14  A0=Exp(-0.2906E+01-0.1069E+00*SB -0.1055E+01*SB2+0.2496E+00*SB3)
      A1=-0.2875E+00+0.6571E-01*SB -0.1987E-01*SB2-0.1800E-02*SB3
      A2= 0.9854E+01-0.2715E+00*SB -0.7407E+00*SB2+0.2888E+00*SB3
      A3= 0.1583E+02-0.7687E+01*SB +0.3428E+01*SB2-0.3327E+00*SB3
      A4= 0.9763E+00+0.7599E-01*SB -0.2128E+00*SB2+0.6852E-01*SB3
      A5=-0.8444E-02+0.9434E+00*SB +0.4152E+00*SB2-0.1481E+00*SB3
      goto 100
c   Ifl =  -2
  15  A0=Exp(-0.2328E+01-0.3061E+01*SB +0.3620E+01*SB2-0.1602E+01*SB3)
      A1=-0.3358E+00+0.3198E+00*SB -0.4210E+00*SB2+0.1571E+00*SB3
      A2= 0.8478E+01-0.3112E+01*SB +0.5243E+01*SB2-0.2255E+01*SB3
      A3= 0.1971E+02+0.3389E+00*SB -0.5268E+01*SB2+0.2099E+01*SB3
      A4= 0.1128E+01-0.4701E+00*SB +0.7779E+00*SB2-0.3506E+00*SB3
      A5=-0.4708E+00+0.3341E+01*SB -0.3375E+01*SB2+0.1353E+01*SB3
      goto 100
c   Ifl =  -3
  16  A0=Exp(-0.3780E+01+0.2499E+01*SB -0.4962E+01*SB2+0.1936E+01*SB3)
      A1=-0.2639E+00-0.1575E+00*SB +0.3584E+00*SB2-0.1646E+00*SB3
      A2= 0.8082E+01+0.2794E+01*SB -0.5438E+01*SB2+0.2321E+01*SB3
      A3= 0.1811E+02-0.2000E+02*SB +0.1951E+02*SB2-0.6904E+01*SB3
      A4= 0.9822E+00+0.4972E+00*SB -0.8690E+00*SB2+0.3415E+00*SB3
      A5= 0.1772E+00-0.6078E+00*SB +0.3341E+01*SB2-0.1473E+01*SB3
      goto 100
c   Ifl =  -4
  17  A0=SB** 0.1122E+01*Exp(-0.4232E+01-0.1808E+01*SB +0.5348E+00*SB2)
      A1=-0.2824E+00+0.5846E+00*SB -0.7230E+00*SB2+0.2419E+00*SB3
      A2= 0.5683E+01-0.2948E+01*SB +0.5916E+01*SB2-0.2560E+01*SB3
      A3= 0.2051E+01+0.4795E+01*SB -0.4271E+01*SB2+0.4174E+00*SB3
      A4= 0.1737E+00+0.1717E+01*SB -0.1978E+01*SB2+0.6643E+00*SB3
      A5= 0.8689E+00+0.3500E+01*SB -0.3283E+01*SB2+0.1026E+01*SB3
      goto 100
c   Ifl =  -5
  18  A0=SB** 0.9906E+00*Exp(-0.1496E+01-0.6576E+01*SB +0.1569E+01*SB2)
      A1=-0.2140E+00-0.6419E-01*SB -0.2741E-02*SB2+0.3185E-02*SB3
      A2= 0.5781E+01+0.1049E+00*SB -0.3930E+00*SB2+0.5174E+00*SB3
      A3=-0.9420E+00+0.5511E+00*SB +0.8817E+00*SB2+0.1903E+01*SB3
      A4= 0.2418E-01+0.4232E-01*SB -0.1244E-01*SB2-0.2365E-01*SB3
      A5= 0.7664E+00+0.1794E+01*SB -0.4917E+00*SB2-0.1284E+00*SB3
      goto 100
c   Ifl =  -6
  19  A0=SB** 0.1000E+01*Exp(-0.8460E+01+0.1154E+01*SB +0.8838E+01*SB2)
      A1=-0.4316E-01-0.2976E+00*SB +0.3174E+00*SB2-0.1429E+01*SB3
      A2= 0.4910E+01+0.2273E+01*SB +0.5631E+01*SB2-0.1994E+02*SB3
      A3= 0.1190E+02-0.2000E+02*SB -0.2000E+02*SB2+0.1292E+02*SB3
      A4= 0.5771E+00-0.2552E+00*SB +0.7510E+00*SB2+0.6923E+00*SB3
      A5= 0.4402E+01-0.1627E+01*SB -0.2085E+01*SB2-0.6737E+01*SB3
      goto 100

 2    Goto(21,22,23,24,25,26,27,28,29)Iflv
c   Ifl =   2
  21  A0=Exp( 0.1141E+00+0.4764E+00*SB -0.1745E+01*SB2+0.7728E+00*SB3)
      A1= 0.4275E+00-0.1290E+00*SB +0.3609E+00*SB2-0.1689E+00*SB3
      A2= 0.3000E+01+0.2946E+01*SB -0.4117E+01*SB2+0.1989E+01*SB3
      A3=-0.1302E+01+0.2322E+01*SB -0.4258E+01*SB2+0.2109E+01*SB3
      A4= 0.2586E+01-0.1920E+00*SB -0.3754E+00*SB2+0.2731E+00*SB3
      A5=-0.2251E+00-0.5374E+00*SB +0.2245E+01*SB2-0.1034E+01*SB3
      goto 100
c   Ifl =   1
  22  A0=Exp( 0.1907E+00+0.4205E-01*SB +0.2752E+00*SB2-0.3171E+00*SB3)
      A1= 0.4611E+00+0.2331E-01*SB -0.3403E-01*SB2+0.3174E-01*SB3
      A2= 0.3504E+01+0.5739E+00*SB +0.2676E+00*SB2-0.1553E+00*SB3
      A3= 0.7452E+01-0.6742E+01*SB +0.2849E+01*SB2-0.1964E+00*SB3
      A4= 0.1116E+01-0.3435E+00*SB +0.2865E+00*SB2-0.1288E+00*SB3
      A5= 0.6659E-01+0.2714E+00*SB -0.2688E+00*SB2+0.2763E+00*SB3
      goto 100
c   Ifl =   0
  23  A0=Exp(-0.7631E+00-0.7241E+00*SB -0.1170E+01*SB2+0.5343E+00*SB3)
      A1=-0.3573E+00+0.3469E+00*SB -0.3396E+00*SB2+0.9188E-01*SB3
      A2= 0.5604E+01+0.7458E+00*SB -0.5082E+00*SB2+0.1844E+00*SB3
      A3= 0.1549E+02-0.1809E+02*SB +0.1162E+02*SB2-0.3483E+01*SB3
      A4= 0.9881E+00+0.1364E+00*SB -0.4421E+00*SB2+0.2051E+00*SB3
      A5=-0.9505E-01+0.3259E+01*SB -0.1547E+01*SB2+0.2918E+00*SB3
      goto 100
c   Ifl =  -1
  24  A0=Exp(-0.2740E+01-0.7987E-01*SB -0.9015E+00*SB2-0.9872E-01*SB3)
      A1=-0.3909E+00+0.1244E+00*SB -0.4487E-01*SB2+0.1277E-01*SB3
      A2= 0.9163E+01+0.2823E+00*SB -0.7720E+00*SB2-0.9360E-02*SB3
      A3= 0.1080E+02-0.3915E+01*SB -0.1153E+01*SB2+0.2649E+01*SB3
      A4= 0.9894E+00-0.1647E+00*SB -0.9426E-02*SB2+0.2945E-02*SB3
      A5=-0.3395E+00+0.6998E+00*SB +0.7000E+00*SB2-0.6730E-01*SB3
      goto 100
c   Ifl =  -2
  25  A0=Exp(-0.2449E+01-0.3513E+01*SB +0.4529E+01*SB2-0.2031E+01*SB3)
      A1=-0.4050E+00+0.3411E+00*SB -0.3669E+00*SB2+0.1109E+00*SB3
      A2= 0.7470E+01-0.2982E+01*SB +0.5503E+01*SB2-0.2419E+01*SB3
      A3= 0.1503E+02+0.1638E+01*SB -0.8772E+01*SB2+0.3852E+01*SB3
      A4= 0.1137E+01-0.1006E+01*SB +0.1485E+01*SB2-0.6389E+00*SB3
      A5=-0.5299E+00+0.3160E+01*SB -0.3104E+01*SB2+0.1219E+01*SB3
      goto 100
c   Ifl =  -3
  26  A0=Exp(-0.3640E+01+0.1250E+01*SB -0.2914E+01*SB2+0.8390E+00*SB3)
      A1=-0.3595E+00-0.5259E-01*SB +0.3122E+00*SB2-0.1642E+00*SB3
      A2= 0.7305E+01+0.9727E+00*SB -0.9788E+00*SB2-0.5193E-01*SB3
      A3= 0.1198E+02-0.1799E+02*SB +0.2614E+02*SB2-0.1091E+02*SB3
      A4= 0.9882E+00-0.6101E+00*SB +0.9737E+00*SB2-0.4935E+00*SB3
      A5=-0.1186E+00-0.3231E+00*SB +0.3074E+01*SB2-0.1274E+01*SB3
      goto 100
c   Ifl =  -4
  27  A0=SB** 0.1122E+01*Exp(-0.3718E+01-0.1335E+01*SB +0.1651E-01*SB2)
      A1=-0.4719E+00+0.7509E+00*SB -0.8420E+00*SB2+0.2901E+00*SB3
      A2= 0.6194E+01-0.1641E+01*SB +0.4907E+01*SB2-0.2523E+01*SB3
      A3= 0.4426E+01-0.4270E+01*SB +0.6581E+01*SB2-0.3474E+01*SB3
      A4= 0.2683E+00+0.9876E+00*SB -0.7612E+00*SB2+0.1780E+00*SB3
      A5=-0.4547E+00+0.4410E+01*SB -0.3712E+01*SB2+0.1245E+01*SB3
      goto 100
c   Ifl =  -5
  28  A0=SB** 0.9838E+00*Exp(-0.2548E+01-0.7660E+01*SB +0.3702E+01*SB2)
      A1=-0.3122E+00-0.2120E+00*SB +0.5716E+00*SB2-0.3773E+00*SB3
      A2= 0.6257E+01-0.8214E-01*SB -0.2537E+01*SB2+0.2981E+01*SB3
      A3=-0.6723E+00+0.2131E+01*SB +0.9599E+01*SB2-0.7910E+01*SB3
      A4= 0.9169E-01+0.4295E-01*SB -0.5017E+00*SB2+0.3811E+00*SB3
      A5= 0.2402E+00+0.2656E+01*SB -0.1586E+01*SB2+0.2880E+00*SB3
      goto 100
c   Ifl =  -6
  29  A0=SB** 0.1001E+01*Exp(-0.6934E+01+0.3050E+01*SB -0.6943E+00*SB2)
      A1=-0.1713E+00-0.5167E+00*SB +0.1241E+01*SB2-0.1703E+01*SB3
      A2= 0.6169E+01+0.3023E+01*SB -0.1972E+02*SB2+0.1069E+02*SB3
      A3= 0.4439E+01-0.1746E+02*SB +0.1225E+02*SB2+0.8350E+00*SB3
      A4= 0.5458E+00-0.4586E+00*SB +0.9089E+00*SB2-0.4049E+00*SB3
      A5= 0.3207E+01-0.3362E+01*SB +0.5877E+01*SB2-0.7659E+01*SB3
      goto 100

 3    Goto(31,32,33,34,35,36,37,38,39)Iflv
c   Ifl =   2
  31  A0=Exp( 0.3961E+00+0.4914E+00*SB -0.1728E+01*SB2+0.7257E+00*SB3)
      A1= 0.4162E+00-0.1419E+00*SB +0.3680E+00*SB2-0.1618E+00*SB3
      A2= 0.3248E+01+0.3028E+01*SB -0.4307E+01*SB2+0.1920E+01*SB3
      A3=-0.1100E+01+0.2184E+01*SB -0.3820E+01*SB2+0.1717E+01*SB3
      A4= 0.2082E+01-0.2756E+00*SB +0.3043E+00*SB2-0.1260E+00*SB3
      A5=-0.4822E+00-0.5706E+00*SB +0.2243E+01*SB2-0.9760E+00*SB3
      goto 100
c   Ifl =   1
  32  A0=Exp( 0.2148E+00+0.5814E-01*SB +0.2734E+00*SB2-0.2902E+00*SB3)
      A1= 0.4810E+00+0.1657E-01*SB -0.3800E-01*SB2+0.3125E-01*SB3
      A2= 0.3509E+01+0.3923E+00*SB +0.4010E+00*SB2-0.1932E+00*SB3
      A3= 0.7055E+01-0.6552E+01*SB +0.3466E+01*SB2-0.5657E+00*SB3
      A4= 0.1061E+01-0.3453E+00*SB +0.4089E+00*SB2-0.1817E+00*SB3
      A5= 0.8687E-01+0.2548E+00*SB -0.2967E+00*SB2+0.2647E+00*SB3
      goto 100
c   Ifl =   0
  33  A0=Exp(-0.4665E+00-0.7554E+00*SB -0.3323E+00*SB2-0.2734E-04*SB3)
      A1=-0.3359E+00+0.2395E+00*SB -0.2377E+00*SB2+0.7059E-01*SB3
      A2= 0.5451E+01+0.6086E+00*SB +0.8606E-01*SB2-0.1425E+00*SB3
      A3= 0.1026E+02-0.9352E+01*SB +0.4879E+01*SB2-0.1150E+01*SB3
      A4= 0.9935E+00-0.5017E-01*SB -0.1707E-01*SB2-0.1464E-02*SB3
      A5=-0.4160E-01+0.2305E+01*SB -0.1063E+01*SB2+0.3211E+00*SB3
      goto 100
c   Ifl =  -1
  34  A0=Exp(-0.3323E+01+0.2296E+00*SB -0.1109E+01*SB2+0.2223E+00*SB3)
      A1=-0.3410E+00+0.8847E-01*SB -0.1111E-01*SB2-0.5927E-02*SB3
      A2= 0.9753E+01-0.5182E+00*SB -0.4670E+00*SB2+0.1921E+00*SB3
      A3= 0.1977E+02-0.1600E+02*SB +0.9481E+01*SB2-0.1864E+01*SB3
      A4= 0.9818E+00+0.2839E-02*SB -0.1188E+00*SB2+0.3584E-01*SB3
      A5=-0.7934E-01+0.1004E+01*SB +0.3704E+00*SB2-0.1220E+00*SB3
      goto 100
c   Ifl =  -2
  35  A0=Exp(-0.2714E+01-0.2868E+01*SB +0.3700E+01*SB2-0.1671E+01*SB3)
      A1=-0.3893E+00+0.3341E+00*SB -0.3897E+00*SB2+0.1420E+00*SB3
      A2= 0.8359E+01-0.3267E+01*SB +0.5327E+01*SB2-0.2245E+01*SB3
      A3= 0.2359E+02-0.5669E+01*SB -0.4602E+01*SB2+0.3153E+01*SB3
      A4= 0.1106E+01-0.4745E+00*SB +0.7739E+00*SB2-0.3417E+00*SB3
      A5=-0.5557E+00+0.3433E+01*SB -0.3390E+01*SB2+0.1354E+01*SB3
      goto 100
c   Ifl =  -3
  36  A0=Exp(-0.3985E+01+0.2855E+01*SB -0.5208E+01*SB2+0.1937E+01*SB3)
      A1=-0.3337E+00-0.1150E+00*SB +0.3691E+00*SB2-0.1709E+00*SB3
      A2= 0.7968E+01+0.3641E+01*SB -0.6599E+01*SB2+0.2642E+01*SB3
      A3= 0.1873E+02-0.1999E+02*SB +0.1734E+02*SB2-0.5813E+01*SB3
      A4= 0.9731E+00+0.5082E+00*SB -0.8780E+00*SB2+0.3231E+00*SB3
      A5=-0.5542E-01-0.4189E+00*SB +0.3309E+01*SB2-0.1439E+01*SB3
      goto 100
c   Ifl =  -4
  37  A0=SB** 0.1105E+01*Exp(-0.3952E+01-0.1901E+01*SB +0.5137E+00*SB2)
      A1=-0.3543E+00+0.6055E+00*SB -0.6941E+00*SB2+0.2278E+00*SB3
      A2= 0.5955E+01-0.2629E+01*SB +0.5337E+01*SB2-0.2300E+01*SB3
      A3= 0.1933E+01+0.4882E+01*SB -0.3810E+01*SB2+0.2290E+00*SB3
      A4= 0.1806E+00+0.1655E+01*SB -0.1893E+01*SB2+0.6395E+00*SB3
      A5= 0.4790E+00+0.3612E+01*SB -0.3152E+01*SB2+0.9684E+00*SB3
      goto 100
c   Ifl =  -5
  38  A0=SB** 0.9818E+00*Exp(-0.1825E+01-0.7464E+01*SB +0.2143E+01*SB2)
      A1=-0.2604E+00-0.1400E+00*SB +0.1702E+00*SB2-0.8476E-01*SB3
      A2= 0.6005E+01+0.6275E+00*SB -0.2535E+01*SB2+0.2219E+01*SB3
      A3=-0.9067E+00+0.1149E+01*SB +0.1974E+01*SB2+0.4716E+01*SB3
      A4= 0.3915E-01+0.5945E-01*SB -0.9844E-01*SB2+0.2783E-01*SB3
      A5= 0.5500E+00+0.1994E+01*SB -0.6727E+00*SB2-0.1510E+00*SB3
      goto 100
c   Ifl =  -6
  39  A0=SB** 0.1002E+01*Exp(-0.8553E+01+0.3793E+00*SB +0.9998E+01*SB2)
      A1=-0.5870E-01-0.2792E+00*SB +0.6526E+00*SB2-0.1984E+01*SB3
      A2= 0.4716E+01+0.4473E+00*SB +0.1128E+02*SB2-0.1937E+02*SB3
      A3= 0.1289E+02-0.1742E+02*SB -0.1983E+02*SB2-0.9274E+00*SB3
      A4= 0.5647E+00-0.2732E+00*SB +0.1074E+01*SB2+0.5981E+00*SB3
      A5= 0.4390E+01-0.1262E+01*SB -0.9026E+00*SB2-0.9394E+01*SB3
      goto 100

 311  stop 'This option is not currently supported.'

 100  Ctq3df = A0 *(x**A1) *((D1-x)**A2) *(D1+A3*(x**A4))
     $            *(log(D1+D1/x))**A5

      if(Ctq3df.lt.D0) then
        Ctq3df = D0
        Irt=1
      endif

      Ist = Iset

      Lp  = Iprtn
      Qsto = QQ

      Return
C                                  -----------------------
      ENTRY Wlamd3 (Iset, Iorder, Neff)

C     Returns the EFFECTIVE QCD lambda values for order=Iorder and
C     effective # of flavors = Neff for each of the PDF sets.

      Iorder = Iord (Iset)
      Wlamd3 = VLM  (Neff, Iset)

      RETURN

C                         *************************
      END


C ==========================================================================
CsB   Code below is the CTEQ4 PDF. Needed for processes with no input grids.
C============================================================================
C                CTEQ Parton Distribution Functions: Version 4.6
C                             June 21, 1996
C                   Modified: 10/17/96, 1/7/97, 1/15/97
C                             2/17/97, 2/21/97
C                   Last Modified on April 2, 1997
C
C   Ref[1]: "IMPROVED PARTON DISTRIBUTIONS FROM GLOBAL ANALYSIS OF RECENT DEEP
C         INELASTIC SCATTERING AND INCLUSIVE JET DATA"
C   By: H.L. Lai, J. Huston, S. Kuhlmann, F. Olness, J. Owens, D. Soper
C       W.K. Tung, H. Weerts
C       Phys. Rev. D55, 1280 (1997)
C
C   Ref[2]: "CHARM PRODUCTION AND PARTON DISTRIBUTIONS"
C   By: H.L. Lai and W.K. Tung
C       MSU-HEP-61222, CTEQ-622, e-Print Archive: hep-ph/9701256
C       to appear in Z. Phys.
C
C   This package contains 13 sets of CTEQ4 PDF's. Details are:
C ---------------------------------------------------------------------------
C  Iset   PDF        Description       Alpha_s(Mz)  Lam4  Lam5   Table_File
C ---------------------------------------------------------------------------
C   1    CTEQ4M   Standard MSbar scheme   0.116     298   202    cteq4m.tbl
C   2    CTEQ4D   Standard DIS scheme     0.116     298   202    cteq4d.tbl
C   3    CTEQ4L   Leading Order           0.132     236   181    cteq4l.tbl
C   4    CTEQ4A1  Alpha_s series          0.110     215   140    cteq4a1.tbl
C   5    CTEQ4A2  Alpha_s series          0.113     254   169    cteq4a2.tbl
C   6    CTEQ4A3            ( same as CTEQ4M )
C   7    CTEQ4A4  Alpha_s series          0.119     346   239    cteq4a4.tbl
C   8    CTEQ4A5  Alpha_s series          0.122     401   282    cteq4a5.tbl
C   9    CTEQ4HJ  High Jet                0.116     303   206    cteq4hj.tbl
C   10   CTEQ4LQ  Low Q0                  0.114     261   174    cteq4lq.tbl
C ---------------------------------------------------------------------------
C   11   CTEQ4HQ  Heavy Quark             0.116     298   202    cteq4hq.tbl
C   12   CTEQ4HQ1 Heavy Quark:Q0=1,Mc=1.3 0.116     298   202    cteq4hq1.tbl
C   13   CTEQ4F3  Nf=3 FixedFlavorNumber  0.106     (Lam3=385)   cteq4f3.tbl
C   14   CTEQ4F4  Nf=4 FixedFlavorNumber  0.111     292   XXX    cteq4f4.tbl
C ---------------------------------------------------------------------------
C
C   The available applied range is 10^-5 < x < 1 and 1.6 < Q < 10,000 (GeV)
C   except CTEQ4LQ(4HQ1) for which Q starts at a lower value of 0.7(1.0) GeV.
C   Lam5 (Lam4, Lam3) represents Lambda value (in MeV) for 5 (4,3) flavors.
C   The matching alpha_s between 4 and 5 flavors takes place at Q=5.0 GeV,
C   which is defined as the bottom quark mass, whenever it can be applied.
C
C   The Table_Files are assumed to be in the working directory.
C
C   Before using the PDF, it is necessary to do the initialization by
C       Call SetCtq4(Iset)
C   where Iset is the desired PDF specified in the above table.
C
C   The function Ctq4Pdf (Iparton, X, Q)
C   returns the parton distribution inside the proton for parton [Iparton]
C   at [X] Bjorken_X and scale [Q] (GeV) in PDF set [Iset].
C   Iparton  is the parton label (5, 4, 3, 2, 1, 0, -1, ......, -5)
C                            for (b, c, s, d, u, g, u_bar, ..., b_bar),
C      whereas CTEQ4F3 has, by definition, only 3 flavors and gluon;
C              CTEQ4F4 has only 4 flavors and gluon.
C
C   For detailed information on the parameters used, e.q. quark masses,
C   QCD Lambda, ... etc.,  see info lines at the beginning of the
C   Table_Files.
C
C   These programs, as provided, are in double precision.  By removing the
C   "Implicit Double Precision" lines, they can also be run in single
C   precision.
C
C   If you have detailed questions concerning these CTEQ4 distributions,
C   or if you find problems/bugs using this package, direct inquires to
C   Hung-Liang Lai(Lai_H@pa.msu.edu) or Wu-Ki Tung(Tung@pa.msu.edu).
C
C===========================================================================

      Function Ctq4Pdf (Iparton, X, Q)
      Implicit Double Precision (A-H,O-Z)
      Logical Warn
      Common
     > / CtqPar2 / Nx, Nt, NfMx, MxVal
     > / QCDtable /  Alambda, Nfl, Iorder

      Data Warn /.true./
      save Warn

      If (X .lt. 0D0 .or. X .gt. 1D0) Then
	Print *, 'X out of range in Ctq4Pdf: ', X
	Stop
      Endif
      If (Q .lt. Alambda) Then
	Print *, 'Q out of range in Ctq4Pdf: ', Q
	Stop
      Endif
      If ((Iparton .lt. -NfMx .or. Iparton .gt. NfMx)) Then
         If (Warn) Then
C        put a warning for calling extra flavor.
	     Warn = .false.
	     Print *, 'Warning: Iparton out of range in Ctq4Pdf: '
     >              , Iparton
         Endif
         Ctq4Pdf = 0D0
         Return
      Endif

      Ctq4Pdf = PartonX4 (Iparton, X, Q)
      if(Ctq4Pdf.lt.0.D0)  Ctq4Pdf = 0.D0

      Return

C                             ********************
      End

      FUNCTION PartonX4 (IPRTN, X, Q)
C
C   Given the parton distribution function in the array Upd in
C   COMMON / CtqPar1 / , this routine fetches u(fl, x, q) at any value of
C   x and q using Mth-order polynomial interpolation for x and Ln(Q/Lambda).
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      PARAMETER (MXX = 201, MXQ = 25, MXF = 6)
      PARAMETER (MXPQX = (MXF *2 +2) * MXQ * MXX)
      PARAMETER (M= 2, M1 = M + 1)
C
      Common
     > / CtqPar1s / Al, XV(0:MXX), QL(0:MXQ), UPD(MXPQX)
     > / CtqPar2 / Nx, Nt, NfMx, MxVal
     > / XQrange / Qini, Qmax, Xmin
C
      Dimension Fq(M1), Df(M1)
C                                                 Work with Log (Q)
      QG  = LOG (Q/AL)

C                           Find lower end of interval containing X
      JL = -1
      JU = Nx+1
 11   If (JU-JL .GT. 1) Then
         JM = (JU+JL) / 2
         If (X .GT. XV(JM)) Then
            JL = JM
         Else
            JU = JM
         Endif
         Goto 11
      Endif

      Jx = JL - (M-1)/2
      If (X .lt. Xmin) Then
         Print '(A, 2(1pE12.4))',
     >     ' WARNING: X < Xmin, extrapolation used; X, Xmin =', X, Xmin
         If (Jx .LT. 0) Jx = 0
      Elseif (Jx .GT. Nx-M) Then
         Jx = Nx - M
      Endif
C                                    Find the interval where Q lies
      JL = -1
      JU = NT+1
 12   If (JU-JL .GT. 1) Then
         JM = (JU+JL) / 2
         If (QG .GT. QL(JM)) Then
            JL = JM
         Else
            JU = JM
         Endif
         Goto 12
      Endif

      Jq = JL - (M-1)/2
      If (Jq .LT. 0) Then
         Jq = 0
         If (Q .lt. Qini)  Print '(A, 2(1pE12.4))',
     >     ' WARNING: Q < Qini, extrapolation used; Q, Qini =', Q, Qini
      Elseif (Jq .GT. Nt-M) Then
         Jq = Nt - M
         If (Q .gt. Qmax)  Print '(A, 2(1pE12.4))',
     >     ' WARNING: Q > Qmax, extrapolation used; Q, Qmax =', Q, Qmax
      Endif

      If (Iprtn .GE. 3) Then
         Ip = - Iprtn
      Else
         Ip = Iprtn
      EndIf
C                             Find the off-set in the linear array Upd
      JFL = Ip + NfMx
      J0  = (JFL * (NT+1) + Jq) * (NX+1) + Jx
C
C                                           Now interpolate in x for M1 Q's
      Do 21 Iq = 1, M1
         J1 = J0 + (Nx+1)*(Iq-1) + 1
         Call Polint (XV(Jx), Upd(J1), M1, X, Fq(Iq), Df(Iq))
 21   Continue
C                                          Finish off by interpolating in Q
      Call Polint (QL(Jq), Fq(1), M1, QG, Ftmp, Ddf)

      PartonX4 = Ftmp
C
      RETURN
C                        ****************************
      END

      Subroutine SetCtq4 (Iset)
      Implicit Double Precision (A-H,O-Z)
      Parameter (Isetmax=14)
      Character Flnm(Isetmax)*18, Tablefile*40
      Data (Flnm(I), I=1,Isetmax)
     > / 'input/cteq4m.tbl', 'input/cteq4d.tbl', 'input/cteq4l.tbl'
     > , 'input/cteq4a1.tbl', 'input/cteq4a2.tbl', 'input/cteq4m.tbl'
     > , 'input/cteq4a4.tbl', 'input/cteq4a5.tbl', 'input/cteq4hj.tbl'
     > , 'input/cteq4lq.tbl', 'input/cteq4hq.tbl', 'input/cteq4hq1.tbl'
     > , 'input/cteq4f3.tbl', 'input/cteq4f4.tbl' /
      Data Tablefile / 'test.tbl' /
      Data Isetold, Isetmin, Isettest / -987, 1, 911 /
      save

 21   Continue
C             If data file not initialized, do so.
      If(Iset.ne.Isetold) then
	 IU= NextUn()
         If (Iset .eq. Isettest) then
            Print* ,'Opening ', Tablefile
c 21         Open(IU, File=Tablefile, Status='OLD', Err=101)
             Open(IU, File=Tablefile, Status='OLD', Err=101)
         ElseIf (Iset.lt.Isetmin .or. Iset.gt.Isetmax) Then
	    Print *, 'Invalid Iset number in SetCtq4 :', Iset
	    Stop
         Else
            Tablefile=Flnm(Iset)
            Open(IU, File=Tablefile, Status='OLD', Err=100)
	 Endif
         Call ReadTbl4 (IU)
         Close (IU)
	 Isetold=Iset
      Endif
      Return

 100  Print *, ' Data file ', Tablefile, ' cannot be opened '
     >//'in SetCtq4!!'
      Stop
 101  Print*, Tablefile, ' cannot be opened '
      Print*, 'Please input the .tbl file:'
      Read (*,'(A)') Tablefile
      Goto 21
C                             ********************
      End

C ----------------------------------------------------------------------------
      Subroutine ReadTbl4 (Nu)
C ----------------------------------------------------------------------------
      Implicit Double Precision (A-H,O-Z)
      Character Line*80
      PARAMETER (MXX = 201, MXQ = 25, MXF = 6)
      PARAMETER (MXPQX = (MXF *2 +2) * MXQ * MXX)
      Common
     > / CtqPar1s / Al, XV(0:MXX), QL(0:MXQ), UPD(MXPQX)
     > / CtqPar2 / Nx, Nt, NfMx, MxVal
     > / XQrange / Qini, Qmax, Xmin
     > / QCDtable /  Alambda, Nfl, Iorder
     > / Masstbl / Amass(6)

      Read  (Nu, '(A)') Line
      Read  (Nu, '(A)') Line
      Read  (Nu, *) Dr, Fl, Al, (Amass(I),I=1,6)
      Iorder = Nint(Dr)
      Nfl = Nint(Fl)
      Alambda = Al

      Read  (Nu, '(A)') Line
      Read  (Nu, *) NX,  NT, NfMx

      Read  (Nu, '(A)') Line
      Read  (Nu, *) QINI, QMAX, (QL(I), I =0, NT)

      Read  (Nu, '(A)') Line
      Read  (Nu, *) XMIN, (XV(I), I =0, NX)

      Do 11 Iq = 0, NT
         QL(Iq) = Log (QL(Iq) /Al)
   11 Continue
C
C                  Since quark = anti-quark for nfl>2 at this stage,
C                  we Read  out only the non-redundent data points
C     No of flavors = NfMx (sea) + 1 (gluon) + 2 (valence)

      Nblk = (NX+1) * (NT+1)
      Npts =  Nblk  * (NfMx+3)
      Read  (Nu, '(A)') Line
      Read  (Nu, *, IOSTAT=IRET) (UPD(I), I=1,Npts)

      Return
C                        ****************************
      End

      Function NextUn()
C                                 Returns an unallocated FORTRAN i/o unit.
      Logical EX
C
      Do 10 N = 10, 300
         INQUIRE (UNIT=N, OPENED=EX)
         If (.NOT. EX) then
            NextUn = N
            Return
         Endif
 10   Continue
      Stop ' There is no available I/O unit. '
C               *************************
      End
C

CsB
      SUBROUTINE POLINT_CT4M (XA,YA,N,X,Y,DY)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C                                        Adapted from "Numerical Recipes"
      PARAMETER (NMAX=10)
      DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)
      NS=1
      DIF=ABS(X-XA(1))
      DO 11 I=1,N
        DIFT=ABS(X-XA(I))
        IF (DIFT.LT.DIF) THEN
          NS=I
          DIF=DIFT
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
CCPY          IF(DEN.EQ.0.)PAUSE
          IF(DEN.EQ.0.) THEN
            PRINT*,' DEN.EQ.0 '
          ENDIF
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE
      RETURN
      END

C =========================================================================
C                CTEQ Parton Distribution Functions: Version 5.0
C                             Nov. 1, 1999
C
C   Ref: "GLOBAL QCD ANALYSIS OF PARTON STRUCTURE OF THE NUCLEON:
C         CTEQ5 PPARTON DISTRIBUTIONS"
C
C  hep-ph/9903282; to be published in Eur. Phys. J. C 1999.
C
C  These PDF's use quadratic interpolation of attached tables. A parametrized 
C  version of the same PDF's without external tables is under construction.  
C  They will become available later.
C
C   This package contains 7 sets of CTEQ5 PDF's; plus two updated ones.
C   The undated CTEQ5M1 and CTEQHQ1 use an improved evolution code.
C   Both the original and the updated ones fit current data with comparable
C   accuracy.  The CTEQHQ1 set also involve a different choice of scale,
C   hence differs from CTEQHQ slightly more.  It is preferred over CTEQ5HQ.

C   Details are:
C ---------------------------------------------------------------------------
C  Iset   PDF        Description       Alpha_s(Mz)  Lam4  Lam5   Table_File
C ---------------------------------------------------------------------------
C   1    CTEQ5M   Standard MSbar scheme   0.118     326   226    cteq5m.tbl
C   2    CTEQ5D   Standard DIS scheme     0.118     326   226    cteq5d.tbl
C   3    CTEQ5L   Leading Order           0.127     192   146    cteq5l.tbl
C   4    CTEQ5HJ  Large-x gluon enhanced  0.118     326   226    cteq5hj.tbl
C   5    CTEQ5HQ  Heavy Quark             0.118     326   226    cteq5hq.tbl
C   6    CTEQ5F3  Nf=3 FixedFlavorNumber  0.106     (Lam3=395)   cteq5f3.tbl
C   7    CTEQ5F4  Nf=4 FixedFlavorNumber  0.112     309   XXX    cteq5f4.tbl
C         --------------------------------------------------------
C   8    CTEQ5M1  Improved CTEQ5M         0.118     326   226    cteq5m1.tbl
C   9    CTEQ5HQ1 Improved CTEQ5HQ        0.118     326   226    ctq5hq1.tbl
C ---------------------------------------------------------------------------
C   
C  The available applied range is 10^-5 << x << 1 and 1.0 << Q << 10,000 (GeV).
C   Lam5 (Lam4, Lam3) represents Lambda value (in MeV) for 5 (4,3) flavors. 
C   The matching alpha_s between 4 and 5 flavors takes place at Q=4.5 GeV,  
C   which is defined as the bottom quark mass, whenever it can be applied.
C
C   The Table_Files are assumed to be in the working directory.
C   
C   Before using the PDF, it is necessary to do the initialization by
C       Call SetCtq5(Iset) 
C   where Iset is the desired PDF specified in the above table.
C   
C   The function Ctq5Pdf (Iparton, X, Q)
C   returns the parton distribution inside the proton for parton [Iparton] 
C   at [X] Bjorken_X and scale [Q] (GeV) in PDF set [Iset].
C   Iparton  is the parton label (5, 4, 3, 2, 1, 0, -1, ......, -5)
C                            for (b, c, s, d, u, g, u_bar, ..., b_bar),
C      whereas CTEQ5F3 has, by definition, only 3 flavors and gluon;
C              CTEQ5F4 has only 4 flavors and gluon.
C   
C   For detailed information on the parameters used, e.q. quark masses, 
C   QCD Lambda, ... etc.,  see info lines at the beginning of the 
C   Table_Files.
C
C   These programs, as provided, are in double precision.  By removing the
C   "Implicit Double Precision" lines, they can also be run in single 
C   precision.
C   
C   If you have detailed questions concerning these CTEQ5 distributions, 
C   or if you find problems/bugs using this package, direct inquires to 
C   Hung-Liang Lai(lai@phys.nthu.edu.tw) or Wu-Ki Tung(Tung@pa.msu.edu).
C   
C===========================================================================

      Function Ctq5Pdf (Iparton, X, Q)
      Implicit Double Precision (A-H,O-Z)
      Logical Warn
      Common
     > / CtqPar2 / Nx, Nt, NfMx, MxVal
     > / QCDtable /  Alambda, Nfl, Iorder

      Data Warn /.true./
      save Warn

      If (X .lt. 0D0 .or. X .gt. 1D0) Then
	Print *, 'X out of range in Ctq5Pdf: ', X
	Stop
      Endif
      If (Q .lt. Alambda) Then
	Print *, 'Q out of range in Ctq5Pdf: ', Q
	Stop
      Endif
      If ((Iparton .lt. -NfMx .or. Iparton .gt. NfMx)) Then
         If (Warn) Then
C        put a warning for calling extra flavor.
	     Warn = .false.
	     Print *, 'Warning: Iparton out of range in Ctq5Pdf: '
     >              , Iparton
         Endif
         Ctq5Pdf = 0D0
         Return
      Endif

      Ctq5Pdf = PartonX5 (Iparton, X, Q)
      if(Ctq5Pdf.lt.0.D0)  Ctq5Pdf = 0.D0

      Return

C                             ********************
      End

      FUNCTION PartonX5 (IPRTN, X, Q)
C
C   Given the parton distribution function in the array Upd in
C   COMMON / CtqPar1 / , this routine fetches u(fl, x, q) at any value of
C   x and q using Mth-order polynomial interpolation for x and Ln(Q/Lambda).
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      PARAMETER (MXX = 201, MXQ = 25, MXF = 6)
      PARAMETER (MXPQX = (MXF *2 +2) * MXQ * MXX)
      PARAMETER (M= 2, M1 = M + 1)
C
      Logical First
      Common 
     > / CtqPar1s / Al, XV(0:MXX), QL(0:MXQ), UPD(MXPQX)
     > / CtqPar2 / Nx, Nt, NfMx, MxVal
     > / XQrange / Qini, Qmax, Xmin
C
      Dimension Fq(M1), Df(M1)

      Data First /.true./
      save First
C                                                 Work with Log (Q)
      QG  = LOG (Q/AL)

C                           Find lower end of interval containing X
      JL = -1
      JU = Nx+1
 11   If (JU-JL .GT. 1) Then
         JM = (JU+JL) / 2
         If (X .GT. XV(JM)) Then
            JL = JM
         Else
            JU = JM
         Endif
         Goto 11
      Endif

      Jx = JL - (M-1)/2
      If (X .lt. Xmin .and. First ) Then
         First = .false.
         Print '(A, 2(1pE12.4))', 
     >     ' WARNING: X << Xmin, extrapolation used; X, Xmin =', X, Xmin
         If (Jx .LT. 0) Jx = 0
      Elseif (Jx .GT. Nx-M) Then
         Jx = Nx - M
      Endif
C                                    Find the interval where Q lies
      JL = -1
      JU = NT+1
 12   If (JU-JL .GT. 1) Then
         JM = (JU+JL) / 2
         If (QG .GT. QL(JM)) Then
            JL = JM
         Else
            JU = JM
         Endif
         Goto 12
      Endif

      Jq = JL - (M-1)/2
      If (Jq .LT. 0) Then
         Jq = 0
         If (Q .lt. Qini)  Print '(A, 2(1pE12.4))', 
     >     ' WARNING: Q << Qini, extrapolation used; Q, Qini =', Q, Qini
      Elseif (Jq .GT. Nt-M) Then
         Jq = Nt - M
         If (Q .gt. Qmax)  Print '(A, 2(1pE12.4))', 
     >     ' WARNING: Q > Qmax, extrapolation used; Q, Qmax =', Q, Qmax
      Endif

      If (Iprtn .GE. 3) Then
         Ip = - Iprtn
      Else
         Ip = Iprtn
      EndIf
C                             Find the off-set in the linear array Upd
      JFL = Ip + NfMx
      J0  = (JFL * (NT+1) + Jq) * (NX+1) + Jx
C
C                                           Now interpolate in x for M1 Q's
      Do 21 Iq = 1, M1
         J1 = J0 + (Nx+1)*(Iq-1) + 1
         Call Polint (XV(Jx), Upd(J1), M1, X, Fq(Iq), Df(Iq))
 21   Continue
C                                          Finish off by interpolating in Q
      Call Polint_CTEQ5 (QL(Jq), Fq(1), M1, QG, Ftmp, Ddf)

      PartonX5 = Ftmp
C
      RETURN
C                        ****************************
      END

      Subroutine SetCtq5 (Iset)
      Implicit Double Precision (A-H,O-Z)
      Parameter (Isetmax=9)
      Character Flnm(Isetmax)*18, Tablefile*40
      Data (Flnm(I), I=1,Isetmax)
     > /'input/cteq5m.tbl','input/cteq5d.tbl','input/cteq5l.tbl', 
     >  'input/cteq5hj.tbl','input/cteq5hq.tbl','input/cteq5f3.tbl',
     >  'input/cteq5f4.tbl','input/cteq5m1.tbl','input/ctq5hq1.tbl'/
      Data Tablefile / 'test.tbl' /
      Data Isetold, Isetmin, Isettest / -987, 1, 911 /
      save

C             If data file not initialized, do so.
      If(Iset.ne.Isetold) then
	 IU= NextUn()
         If (Iset .eq. Isettest) then
            Print* ,'Opening ', Tablefile
 21         Open(IU, File=Tablefile, Status='OLD', Err=101)
            GoTo 22
 101        Print*, Tablefile, ' cannot be opened '
            Print*, 'Please input the .tbl file:'
            Read (*,'(A)') Tablefile
            Goto 21
 22         Continue
         ElseIf (Iset.lt.Isetmin .or. Iset.gt.Isetmax) Then
	    Print *, 'Invalid Iset number in SetCtq5 :', Iset
	    Stop
         Else
            Tablefile=Flnm(Iset)
            Open(IU, File=Tablefile, Status='OLD', Err=100)
	 Endif
         Call ReadTbl5 (IU)
         Close (IU)
	 Isetold=Iset
      Endif
      Return

 100  Print *, ' Data file ', Tablefile, ' cannot be opened '
     >//'in SetCtq5!!'
      Stop
C                             ********************
      End

      Subroutine ReadTbl5 (Nu)
      Implicit Double Precision (A-H,O-Z)
      Character Line*80
      PARAMETER (MXX = 201, MXQ = 25, MXF = 6)
      PARAMETER (MXPQX = (MXF *2 +2) * MXQ * MXX)
      Common 
     > / CtqPar1s / Al, XV(0:MXX), QL(0:MXQ), UPD(MXPQX)
     > / CtqPar2 / Nx, Nt, NfMx, MxVal
     > / XQrange / Qini, Qmax, Xmin
     > / QCDtable /  Alambda, Nfl, Iorder
     > / Masstbl / Amass(6)
      
      Read  (Nu, '(A)') Line     
      Read  (Nu, '(A)') Line
      Read  (Nu, *) Dr, Fl, Al, (Amass(I),I=1,6)
      Iorder = Nint(Dr)
      Nfl = Nint(Fl)
      Alambda = Al

      Read  (Nu, '(A)') Line 
      Read  (Nu, *) NX,  NT, NfMx

      Read  (Nu, '(A)') Line
      Read  (Nu, *) QINI, QMAX, (QL(I), I =0, NT)

      Read  (Nu, '(A)') Line
      Read  (Nu, *) XMIN, (XV(I), I =0, NX)

      Do 11 Iq = 0, NT
         QL(Iq) = Log (QL(Iq) /Al)
   11 Continue
C
C                  Since quark = anti-quark for nfl>2 at this stage, 
C                  we Read  out only the non-redundent data points
C     No of flavors = NfMx (sea) + 1 (gluon) + 2 (valence) 

      Nblk = (NX+1) * (NT+1)
      Npts =  Nblk  * (NfMx+3)
      Read  (Nu, '(A)') Line
      Read  (Nu, *, IOSTAT=IRET) (UPD(I), I=1,Npts)

      Return
C                        ****************************
      End
C

      SUBROUTINE POLINT_CTEQ5 (XA,YA,N,X,Y,DY)
 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C                                        Adapted from "Numerical Recipes" 
      PARAMETER (NMAX=10)
      DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)
      NS=1
      DIF=ABS(X-XA(1))
      DO 11 I=1,N
        DIFT=ABS(X-XA(I))
        IF (DIFT.LT.DIF) THEN
          NS=I
          DIF=DIFT
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
CCPY          IF(DEN.EQ.0.)PAUSE
          IF(DEN.EQ.0.) THEN
            PRINT*,' DEN.EQ.0 '
          ENDIF
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE
      RETURN
      END

C======================================================================


      Function Ctq6Pdf (Iparton, X, Q)
      Implicit Double Precision (A-H,O-Z)
      Logical Warn
      Common
     > / CtqPar2 / Nx, Nt, NfMx, MxVal
     > / QCDtable /  Alambda, Nfl, Iorder

      Data Warn /.true./
      save Warn

      If (X .lt. 0d0 .or. X .gt. 1D0) Then
        Print *, 'X out of range in Ctq6Pdf: ', X
        Ctq6Pdf = 0D0
        Return
      Endif

      If (Q .lt. Alambda) Then
        Print *, 'Q out of range in Ctq6Pdf: ', Q
        Stop
      Endif

      If ((Iparton .lt. -NfMx .or. Iparton .gt. NfMx)) Then
         If (Warn) Then
C        put a warning for calling extra flavor.
             Warn = .false.
             Print *, 'Warning: Iparton out of range in Ctq6Pdf! '
             Print *, 'Iparton, MxFlvN0: ', Iparton, NfMx
         Endif
         Ctq6Pdf = 0D0
         Return
      Endif

      Ctq6Pdf = PartonX6 (Iparton, X, Q)
      if (Ctq6Pdf.lt.0.D0) Ctq6Pdf = 0.D0

      Return

C                             ********************
      End

      Subroutine SetCtq6 (Iset)
      Implicit Double Precision (A-H,O-Z)
      Parameter (Isetmax0=8)
      Character Flnm(Isetmax0)*6, nn*3, Tablefile*40
      Logical fmtpds
      Data (Flnm(I), I=1,Isetmax0)
     > / 'cteq6m', 'cteq6d', 'cteq6l', 'cteq6l','ctq61.','cteq6s'
     >  ,'ctq65.', 'ctq66.' /
      Data Isetold, Isetmin0, Isetmin1, Isetmax1 /-987,1,100,140/
      Data Isetmin2,Isetmax2 /200,240/
      Data Isetmin3,Isetmax3 /300,340/
      Data Isetmin4,Isetmax4 /400,444/
      Data IsetminS,IsetmaxS /11,15/
      Data IsetmnSp07,IsetmxSp07 /30,34/
      Data IsetmnSm07,IsetmxSm07 /35,37/
      Data IsetmnC07,IsetmxC07 /40,46/
      Data IsetmnC08,IsetmxC08 /450,453/
      Data IsetmnAS08,IsetmxAS08 /460,463/

      Data IsetHQ /21/
      Common /Setchange/ Isetch
      save

C             If data file not initialized, do so.
      If(Iset.ne.Isetold) then
        fmtpds=.true.

        If (Iset.ge.Isetmin0 .and. Iset.le.3) Then
C                                                  Iset = 1,2,3 for 6m, 6d, 6l
          fmtpds=.false.
          Tablefile=Flnm(Iset)//'.tbl'
        Elseif (Iset.eq.4) Then
C                                                             4  (2nd LO fit)
          fmtpds=.false.
          Tablefile=Flnm(Iset)//'1.tbl'
        Elseif (Iset.ge.Isetmin1 .and. Iset.le.Isetmax1) Then
C                                                               101 - 140
          fmtpds=.false.
          write(nn,'(I3)') Iset
          Tablefile=Flnm(1)//nn//'.tbl'
        Elseif (Iset.ge.Isetmin2 .and. Iset.le.Isetmax2) Then
C                                                               200 - 240
          fmtpds=.false.
          write(nn,'(I3)') Iset
          Tablefile=Flnm(5)//nn(2:3)//'.tbl'
        Elseif (Iset.ge.IsetminS .and. Iset.le.IsetmaxS) Then
C                                                               11 - 15
          If(Iset.eq.11) then
            Tablefile=Flnm(6)//'a.pds'
          Elseif(Iset.eq.12) then
            Tablefile=Flnm(6)//'b.pds'
          Elseif(Iset.eq.13) then
            Tablefile=Flnm(6)//'c.pds'
          Elseif(Iset.eq.14) then
            Tablefile=Flnm(6)//'b+.pds'
          Elseif(Iset.eq.15) then
            Tablefile=Flnm(6)//'b-.pds'
          Endif
        Elseif (Iset.eq.IsetHQ) Then
C                                                               21
          TableFile='cteq6hq.pds'
        Elseif (Iset.ge.IsetmnSp07 .and. Iset.le.IsetmxSp07) Then
C                                                    (Cteq6.5S)  30 - 34
          write(nn,'(I2)') Iset
          Tablefile=Flnm(7)//'s+'//nn(2:2)//'.pds'
        Elseif (Iset.ge.IsetmnSm07 .and. Iset.le.IsetmxSm07) Then
C                                                    (Cteq6.5S)  35 - 37
          Is = Iset - 5
          write(nn,'(I2)') Is
          Tablefile=Flnm(7)//'s-'//nn(2:2)//'.pds'
        Elseif (Iset.ge.IsetmnC07 .and. Iset.le.IsetmxC07) Then
C                                                    (Cteq6.5C)  40 - 46
          write(nn,'(I2)') Iset
          Tablefile=Flnm(7)//'c'//nn(2:2)//'.pds'
        Elseif (Iset.ge.Isetmin3 .and. Iset.le.Isetmax3) Then
C                                                    (Cteq6.5)  300 - 340
          write(nn,'(I3)') Iset
          Tablefile=Flnm(7)//nn(2:3)//'.pds'
        Elseif (Iset.ge.Isetmin4 .and. Iset.le.Isetmax4) Then
C                                                    (Cteq6.6)  400 - 444   
          write(nn,'(I3)') Iset
          Tablefile=Flnm(8)//nn(2:3)//'.pds'
        Elseif (Iset.ge.IsetmnC08 .and. Iset.le.IsetmxC08) Then
C                                                   (Cteq6.6C)  451 - 454
          write(nn,'(I3)') Iset 
          Tablefile=Flnm(8)//'c'//nn(3:3)//'.pds'
        Elseif (Iset.ge.IsetmnAS08 .and. Iset.le.IsetmxAS08) Then
C                                                   (Cteq6.6AS)  461 - 454
          write(nn,'(I3)') Iset 
          Tablefile=Flnm(8)//'c'//nn(3:3)//'.pds'
        Elseif (Iset.eq.903) then
          Tablefile='pdf.pds'
        Else
          Print *, 'Invalid Iset number in SetCtq6 :', Iset
          Stop
        Endif
        IU= NextUn()
        Open(IU, File=Tablefile, Status='OLD', Err=100)
 21     Call Readpds (IU,fmtpds)
        Close (IU)
      Print *, ' Data file ', Tablefile, ' opened '
        Isetold=Iset
        Isetch=1
      Endif
      Return

 100  Print *, ' Data file ', Tablefile, ' cannot be opened '
     >  //'in SetCtq6!!'
      Stop
C                             ********************
      End

      Subroutine Readpds (Nu,fmtpds)
      Implicit Double Precision (A-H,O-Z)
      Character Line*80
      Logical fmtpds
      PARAMETER (MXX = 201, MXQ = 25, MXF = 6, MaxVal=4)
      PARAMETER (MXPQX = (MXF+1+MaxVal) * MXQ * MXX)
      Common
     > / CtqPar1 / Al, XV(0:MXX), TV(0:MXQ), UPD(MXPQX)
     > / CtqPar2 / Nx, Nt, NfMx, MxVal
     > / XQrange / Qini, Qmax, Xmin
     > / QCDtable /  Alambda, Nfl, Iorder
     > / Masstbl / Amass(6)

      Read  (Nu, '(A)') Line
      Read  (Nu, '(A)') Line
      Read  (Nu, *) Dr, Fl, Al, (Amass(I),I=1,6)
      Iorder = Nint(Dr)
      Nfl = Nint(Fl)
      Alambda = Al

      Read  (Nu, '(A)') Line
      If(fmtpds) then
C                                               This is the .pds (WKT) format
        Read  (Nu, *) N0, N0, N0, NfMx, MxVal, N0
        If(MxVal.gt.MaxVal) MxVal=3 !old .pds format (read in KF, not MxVal)

        Read  (Nu, '(A)') Line
        Read  (Nu, *) NX,  NT, N0, NG, N0

        Read  (Nu, '(A)') (Line,I=1,NG+2)
        Read  (Nu, *) QINI, QMAX, (aa,TV(I), I =0, NT)

        Read  (Nu, '(A)') Line
        Read  (Nu, *) XMIN, aa, (XV(I), I =1, NX)
        XV(0)=0D0
      Else
C                                               This is the old .tbl (HLL) format
         MxVal=2
         Read  (Nu, *) NX,  NT, NfMx

         Read  (Nu, '(A)') Line
         Read  (Nu, *) QINI, QMAX, (TV(I), I =0, NT)

         Read  (Nu, '(A)') Line
         Read  (Nu, *) XMIN, (XV(I), I =0, NX)

         Do 11 Iq = 0, NT
            TV(Iq) = Log(Log (TV(Iq) /Al))
 11      Continue
      Endif

      Nblk = (NX+1) * (NT+1)
      Npts =  Nblk  * (NfMx+1+MxVal)
      Read  (Nu, '(A)') Line
      Read  (Nu, *, IOSTAT=IRET) (UPD(I), I=1,Npts)

      Return
C                        ****************************
      End

      Function PartonX6 (IPRTN, XX, QQ)

c  Given the parton distribution function in the array U in
c  COMMON / PEVLDT / , this routine interpolates to find
c  the parton distribution at an arbitray point in x and q.
c
      Implicit Double Precision (A-H,O-Z)

      PARAMETER (MXX = 201, MXQ = 25, MXF = 6, MaxVal=4)
      PARAMETER (MXPQX = (MXF+1+MaxVal) * MXQ * MXX)

      Common
     > / CtqPar1 / Al, XV(0:MXX), TV(0:MXQ), UPD(MXPQX)
     > / CtqPar2 / Nx, Nt, NfMx, MxVal
     > / XQrange / Qini, Qmax, Xmin
     > /Setchange/ Isetch

      Dimension fvec(4), fij(4)
      Dimension xvpow(0:mxx)
      Data OneP / 1.00001 /
      Data xpow / 0.3d0 /       !**** choice of interpolation variable
      Data nqvec / 4 /
C      Data ientry / 0 /
      Data X, Q, JX, JQ /-1D0, -1D0, 0, 0/
      Save xvpow
      Save X, Q, JX, JQ, JLX, JLQ
      Save ss, const1, const2, const3, const4, const5, const6
      Save sy2, sy3, s23, tt, t12, t13, t23, t24, t34, ty2, ty3
      Save tmp1, tmp2, tdet

      If((XX.eq.X).and.(QQ.eq.Q)) goto 99
c store the powers used for interpolation on first call...
      if(Isetch .eq. 1) then
         Isetch = 0

         xvpow(0) = 0D0
         do i = 1, nx
            xvpow(i) = xv(i)**xpow
         enddo
      endif

      X = XX
      Q = QQ
      tt = log(log(Q/Al))

c      -------------    find lower end of interval containing x, i.e.,
c                       get jx such that xv(jx) .le. x .le. xv(jx+1)...
      JLx = -1
      JU = Nx+1
 11   If (JU-JLx .GT. 1) Then
         JM = (JU+JLx) / 2
         If (X .Ge. XV(JM)) Then
            JLx = JM
         Else
            JU = JM
         Endif
         Goto 11
      Endif
C                     Ix    0   1   2      Jx  JLx         Nx-2     Nx
C                           |---|---|---|...|---|-x-|---|...|---|---|
C                     x     0  Xmin               x                 1
C
      If     (JLx .LE. -1) Then
        Print '(A,1pE12.4)', 'Severe error: x <= 0 in PartonX6! x = ', x
        Stop
      ElseIf (JLx .Eq. 0) Then
         Jx = 0
      Elseif (JLx .LE. Nx-2) Then

C                For interrior points, keep x in the middle, as shown above
         Jx = JLx - 1
      Elseif (JLx.Eq.Nx-1 .or. x.LT.OneP) Then

C                  We tolerate a slight over-shoot of one (OneP=1.00001),
C              perhaps due to roundoff or whatever, but not more than that.
C                                      Keep at least 4 points >= Jx
         Jx = JLx - 2
      Else
        Print '(A,1pE12.4)', 'Severe error: x > 1 in PartonX6! x = ', x
        Stop
      Endif
C          ---------- Note: JLx uniquely identifies the x-bin; Jx does not.

C                       This is the variable to be interpolated in
      ss = x**xpow

      If (JLx.Ge.2 .and. JLx.Le.Nx-2) Then

c     initiation work for "interior bins": store the lattice points in s...
      svec1 = xvpow(jx)
      svec2 = xvpow(jx+1)
      svec3 = xvpow(jx+2)
      svec4 = xvpow(jx+3)

      s12 = svec1 - svec2
      s13 = svec1 - svec3
      s23 = svec2 - svec3
      s24 = svec2 - svec4
      s34 = svec3 - svec4

      sy2 = ss - svec2
      sy3 = ss - svec3

c constants needed for interpolating in s at fixed t lattice points...
      const1 = s13/s23
      const2 = s12/s23
      const3 = s34/s23
      const4 = s24/s23
      s1213 = s12 + s13
      s2434 = s24 + s34
      sdet = s12*s34 - s1213*s2434
      tmp = sy2*sy3/sdet
      const5 = (s34*sy2-s2434*sy3)*tmp/s12
      const6 = (s1213*sy2-s12*sy3)*tmp/s34

      EndIf

c         --------------Now find lower end of interval containing Q, i.e.,
c                          get jq such that qv(jq) .le. q .le. qv(jq+1)...
      JLq = -1
      JU = NT+1
 12   If (JU-JLq .GT. 1) Then
         JM = (JU+JLq) / 2
         If (tt .GE. TV(JM)) Then
            JLq = JM
         Else
            JU = JM
         Endif
         Goto 12
       Endif

      If     (JLq .LE. 0) Then
         Jq = 0
      Elseif (JLq .LE. Nt-2) Then
C                                  keep q in the middle, as shown above
         Jq = JLq - 1
      Else
C                         JLq .GE. Nt-1 case:  Keep at least 4 points >= Jq.
        Jq = Nt - 3

      Endif
C                                   This is the interpolation variable in Q

      If (JLq.GE.1 .and. JLq.LE.Nt-2) Then
c                                        store the lattice points in t...
      tvec1 = Tv(jq)
      tvec2 = Tv(jq+1)
      tvec3 = Tv(jq+2)
      tvec4 = Tv(jq+3)

      t12 = tvec1 - tvec2
      t13 = tvec1 - tvec3
      t23 = tvec2 - tvec3
      t24 = tvec2 - tvec4
      t34 = tvec3 - tvec4

      ty2 = tt - tvec2
      ty3 = tt - tvec3

      tmp1 = t12 + t13
      tmp2 = t24 + t34

      tdet = t12*t34 - tmp1*tmp2

      EndIf


c get the pdf function values at the lattice points...

 99   If (Iprtn .Gt. MxVal) Then
         Ip = - Iprtn
      Else
         Ip = Iprtn
      EndIf
      jtmp = ((Ip + NfMx)*(NT+1)+(jq-1))*(NX+1)+jx+1

      Do it = 1, nqvec

         J1  = jtmp + it*(NX+1)

       If (Jx .Eq. 0) Then
C                          For the first 4 x points, interpolate x^2*f(x,Q)
C                           This applies to the two lowest bins JLx = 0, 1
C            We can not put the JLx.eq.1 bin into the "interrior" section
C                           (as we do for q), since Upd(J1) is undefined.
         fij(1) = 0
         fij(2) = Upd(J1+1) * XV(1)**2
         fij(3) = Upd(J1+2) * XV(2)**2
         fij(4) = Upd(J1+3) * XV(3)**2
C
C                 Use Polint which allows x to be anywhere w.r.t. the grid

         Call Polint4F (XVpow(0), Fij(1), ss, Fx)

         If (x .GT. 0D0)  Fvec(it) =  Fx / x**2
C                                              Pdf is undefined for x.eq.0
       ElseIf  (JLx .Eq. Nx-1) Then
C                                                This is the highest x bin:

        Call Polint4F (XVpow(Nx-3), Upd(J1), ss, Fx)

        Fvec(it) = Fx

       Else
C                       for all interior points, use Jon's in-line function
C                              This applied to (JLx.Ge.2 .and. JLx.Le.Nx-2)
         sf2 = Upd(J1+1)
         sf3 = Upd(J1+2)

         g1 =  sf2*const1 - sf3*const2
         g4 = -sf2*const3 + sf3*const4

         Fvec(it) = (const5*(Upd(J1)-g1)
     &               + const6*(Upd(J1+3)-g4)
     &               + sf2*sy3 - sf3*sy2) / s23

       Endif

      enddo
C                                   We now have the four values Fvec(1:4)
c     interpolate in t...

      If (JLq .LE. 0) Then
C                         1st Q-bin, as well as extrapolation to lower Q
        Call Polint4F (TV(0), Fvec(1), tt, ff)

      ElseIf (JLq .GE. Nt-1) Then
C                         Last Q-bin, as well as extrapolation to higher Q
        Call Polint4F (TV(Nt-3), Fvec(1), tt, ff)
      Else
C                         Interrior bins : (JLq.GE.1 .and. JLq.LE.Nt-2)
C       which include JLq.Eq.1 and JLq.Eq.Nt-2, since Upd is defined for
C                         the full range QV(0:Nt)  (in contrast to XV)
        tf2 = fvec(2)
        tf3 = fvec(3)

        g1 = ( tf2*t13 - tf3*t12) / t23
        g4 = (-tf2*t34 + tf3*t24) / t23

        h00 = ((t34*ty2-tmp2*ty3)*(fvec(1)-g1)/t12
     &    +  (tmp1*ty2-t12*ty3)*(fvec(4)-g4)/t34)

        ff = (h00*ty2*ty3/tdet + tf2*ty3 - tf3*ty2) / t23
      EndIf

      PartonX6 = ff

      Return
C                                       ********************
      End

      SUBROUTINE POLINT4F (XA,YA,X,Y)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C  The POLINT4 routine is based on the POLINT routine from "Numerical Recipes",
C  but assuming N=4, and ignoring the error estimation.
C  suggested by Z. Sullivan.
      DIMENSION XA(4),YA(4)

      H1=XA(1)-X
      H2=XA(2)-X
      H3=XA(3)-X
      H4=XA(4)-X

      W=YA(2)-YA(1)
      DEN=W/(H1-H2)
      D1=H2*DEN
      C1=H1*DEN

      W=YA(3)-YA(2)
      DEN=W/(H2-H3)
      D2=H3*DEN
      C2=H2*DEN

      W=YA(4)-YA(3)
      DEN=W/(H3-H4)
      D3=H4*DEN
      C3=H3*DEN

      W=C2-D1
      DEN=W/(H1-H3)
      CD1=H3*DEN
      CC1=H1*DEN

      W=C3-D2
      DEN=W/(H2-H4)
      CD2=H4*DEN
      CC2=H2*DEN

      W=CC2-CD1
      DEN=W/(H1-H4)
      DD1=H4*DEN
      DC1=H1*DEN

      If((H3+H4).lt.0D0) Then
         Y=YA(4)+D3+CD2+DD1
      Elseif((H2+H3).lt.0D0) Then
         Y=YA(3)+D2+CD1+DC1
      Elseif((H1+H2).lt.0D0) Then
         Y=YA(2)+C2+CD1+DC1
      ELSE
         Y=YA(1)+C1+CC1+DC1
      ENDIF

      RETURN
C               *************************
      END
C                                                                 April 2009
C End of CTEQ6 routines
C =========================================================================

C =========================================================================
CsB  Routines below are to double check the NLO qq -> vvG process.
C -------------------------------------------------------------------------
      SUBROUTINE AMP2_QQVVG(AMP2UP,AMP2DN,PI,ALFAS,XMZ,
     > GWEAK2,SWS,PLAB)
C -------------------------------------------------------------------------
CCPY APRIL 17, 1996
C THIS IS TO CALCULATE AMPLITUDE**2 FOR
C Q QBAR --> A A, OR Z Z, OR A Z
C
      IMPLICIT NONE
      REAL*8 PI,ALFAS,XMZ,GWEAK2,SWS,PLAB
      REAL*8 AMP2UP,AMP2DN
      DIMENSION PLAB(4,12)
      INTEGER*4 M_u1,M_u2,M_v1,M_v2,M_g,M_zo
      COMMON/QQAAG_INIT/ M_u1,M_u2,M_v1,M_v2,M_g,M_zo

CsB      REAL*8 PLAB
CsB      DIMENSION PLAB(4,12)
      REAL*8 YXP,YXETA,YXPT,YXDELR,YXPP
      COMMON/SECCONV/ YXP(4,30),YXETA(30),YXPT(30),YXDELR(30,30),
     >YXPP(30)

      REAL*8 Q_UP,Q_DN,A_L,A_R,B_L,B_R

C THIS IS FOR THE PROCESS
C     q(u1) + qbar(u2) --> V1(v1) V2(v2) g(g)
C
C HERE IS THE AVERAGED (AMPLITUDE)**2, THE AMPLITUDE INCLUDES BOTH THE
C COLOR AND SPIN FACTORS.
C
      CALL INIT_QQVVG(PI,ALFAS,XMZ,GWEAK2,SWS,PLAB)

C FOR DI-PHOTON PAIR PRODUCTION
C ELECTRIC CHARGE
      Q_UP=2.0/3.0
      Q_DN=-1.0/3.0

      A_L=Q_UP
      A_R=Q_UP
      B_L=Q_UP
      B_R=Q_UP

      CALL AMPS_QQVVG(PI,ALFAS,XMZ,GWEAK2,SWS,A_L,A_R,
     >B_L,B_R,AMP2UP)

C FOR DI-PHOTON PAIR PRODUCTION
      AMP2DN=AMP2UP*(Q_DN/Q_UP)**4

      RETURN
      END

      SUBROUTINE INIT_QQVVG(PI,ALFAS,XMZ,GWEAK2,SWS,PLAB)
      IMPLICIT NONE
      REAL*8 PI,ALFAS,XMZ,GWEAK2,SWS

C THE AMPLITUDE FOR
C     q(u1) + qbar(u2) --> V1(v1) V2(v2) g(g)

      INTEGER*4 M_u1,M_u2,M_v1,M_v2,M_g,M_zo
      COMMON/QQAAG_INIT/ M_u1,M_u2,M_v1,M_v2,M_g,M_zo

      REAL*8 PLAB
      DIMENSION PLAB(4,12)
      REAL*8 YXP,YXETA,YXPT,YXDELR,YXPP
      COMMON/SECCONV/ YXP(4,30),YXETA(30),YXPT(30),YXDELR(30,30),
     >YXPP(30)

      REAL*8 R2,PROP6S,PROP7S,PROP8S,PROP9S,PROP10S,PROP11S
      COMMON/QVGINIT2/ R2,PROP6S,PROP7S,PROP8S,PROP9S,PROP10S,PROP11S
      REAL*8 SPIN,COLORA,CONST1,CONST2,CONST
      COMMON/QVGINIT3/ SPIN,COLORA,CONST1,CONST2,CONST

      REAL*8 PV,PX,PY,PLONG,XMV,TEST
CsB      DIMENSION PV(4),PX(4),PY(4),PLONG(4),xMV
      DIMENSION PV(4),PX(4),PY(4),PLONG(4)

      INTEGER*4 J1,J2
      REAL*8 YXDOT
      EXTERNAL YXDOT

      Logical Debug
      Data Debug /.False./

C      Debug = .True.

      real*8 te1,te2,te3,te4,te5,te6
      real*8 te11,te12,te13,te14,te15,te16

C MOMENTUM ASSIGNMENT
C u1(1), u2(2), v1(3), v2(4), g(5)

      M_u1=1
      M_u2=2
      M_v1=3
      M_v2=4
      M_g=5
C Note: 30 is reserved for zo always
      M_zo=30

      DO 101 J1=1,30
       DO 101 J2=1,4
        YXP(J2,J1)=0.0D0
101   CONTINUE
      DO 110 J1=1,5
       DO 110 J2=1,4
        YXP(J2,J1)=PLAB(J2,J1)
110   CONTINUE
      DO 120 J1=1,4
       YXP(J1,6)=YXP(J1,1)-YXP(J1,5)
       YXP(J1,7)=YXP(J1,4)-YXP(J1,2)
       YXP(J1,8)=YXP(J1,1)-YXP(J1,3)
       YXP(J1,9)=YXP(J1,2)-YXP(J1,5)
       YXP(J1,10)=YXP(J1,3)-YXP(J1,2)
       YXP(J1,11)=YXP(J1,1)-YXP(J1,4)
120   CONTINUE
      DO 122 J1=1,30
        YXPT(J1)=DSQRT(YXP(1,J1)**2+YXP(2,J1)**2)
        IF((YXP(4,J1)-YXP(3,J1)).GT.0.0 .AND.
     >  (YXP(4,J1)+YXP(3,J1)).GT.0.0)THEN
          YXETA(J1)=0.5*DLOG((YXP(4,J1)+YXP(3,J1))
     >    /(YXP(4,J1)-YXP(3,J1)))
        ELSEIF((YXP(4,J1)-YXP(3,J1)).LE.0.)THEN
          YXETA(J1)=20.
        ELSEIF((YXP(4,J1)+YXP(3,J1)).LE.0.)THEN
          YXETA(J1)=-20.
        END IF
122   CONTINUE

      R2=SQRT(2.0D0)

C ASSIGN THE POLARIZATION VECTORS FOR V1(3), V2(4), AND GLUON(5)
C GAUGE BOSONS

C FOR V1
      DO J1=1,4
       PV(J1)=YXP(J1,M_v1)
      ENDDO
      CALL FINDPOL(PV,PX,PY,PLONG,XMV)
      DO J1=1,4
C X-POLARIZATION VECTOR OF V1
        YXP(J1,21)=PX(J1)
C Check Ward Identity: Using unit momentum
C        YXP(J1,21)=PV(J1)/PV(4)
C Y-POLARIZATION VECTOR OF V1
        YXP(J1,22)=PY(J1)
C Check Ward Identity: Using unit momentum
C        YXP(J1,22)=PV(J1)/PV(4)
        IF(XMV.GT.1.D-4) THEN
C LONGITUDINAL-POLARIZATION VECTOR OF V1
          YXP(J1,23)=PLONG(J1)
        ENDIF
      ENDDO

C FOR V2
      DO J1=1,4
       PV(J1)=YXP(J1,M_v2)
      ENDDO
      CALL FINDPOL(PV,PX,PY,PLONG,XMV)
      DO J1=1,4
C X-POLARIZATION VECTOR OF V2
        YXP(J1,24)=PX(J1)
C Check Ward Identity: Using unit momentum
C        YXP(J1,24)=PV(J1)/PV(4)
C Y-POLARIZATION VECTOR OF V2
        YXP(J1,25)=PY(J1)
C Check Ward Identity: Using unit momentum
C        YXP(J1,25)=PV(J1)/PV(4)
        IF(XMV.GT.1.D-4) THEN
C LONGITUDINAL-POLARIZATION VECTOR OF V2
          YXP(J1,26)=PLONG(J1)
        ENDIF
      ENDDO

C FOR G
      DO J1=1,4
       PV(J1)=YXP(J1,M_g)
      ENDDO
      CALL FINDPOL(PV,PX,PY,PLONG,XMV)
      DO J1=1,4
C X-POLARIZATION VECTOR OF G
        YXP(J1,27)=PX(J1)
C Y-POLARIZATION VECTOR OF G
        YXP(J1,28)=PY(J1)
        IF(XMV.GT.1.D-4) THEN
C LONGITUDINAL-POLARIZATION VECTOR OF G
          YXP(J1,29)=PLONG(J1)
        ENDIF
      ENDDO

      If(Debug) then
CCPY
C      PRINT*,' (YXP(J1,M_V1),J1=1,4),
C     &     (YXP(J1,21),J1=1,4),(YXP(J1,22),J1=1,4) '
C      PRINT*,(YXP(J1,M_V1),J1=1,4),
C     &     (YXP(J1,21),J1=1,4),(YXP(J1,22),J1=1,4)
C      PRINT*,' '
C      PRINT*,' (YXP(J1,M_V2),J1=1,4),
C     &     (YXP(J1,24),J1=1,4),(YXP(J1,25),J1=1,4) '
C      PRINT*,(YXP(J1,M_V2),J1=1,4),
C     &     (YXP(J1,24),J1=1,4),(YXP(J1,25),J1=1,4)
C      PRINT*,' '
C      PRINT*,' (YXP(J1,M_G),J1=1,4),
C     &     (YXP(J1,27),J1=1,4),(YXP(J1,28),J1=1,4) '
C      PRINT*,(YXP(J1,M_G),J1=1,4),
C     &     (YXP(J1,27),J1=1,4),(YXP(J1,28),J1=1,4)
C      PRINT*,' '

      te1=yxdot(m_v1,21)
      te2=yxdot(m_v1,22)
      te3=yxdot(m_v2,24)
      te4=yxdot(m_v2,25)
      te5=yxdot(m_g,27)
      te6=yxdot(m_g,28)

      print*,' te1,te2,te3,te4,te5,te6 '
      print*, te1,te2,te3,te4,te5,te6

      te11=yxdot(21,22)
      te12=yxdot(24,25)
      te13=yxdot(27,28)
      te14=yxdot(m_v1,m_v1)
      te15=yxdot(m_v2,m_v2)
      te16=yxdot(m_g,m_g)

      print*,' te11,te12,te13,te14,te15,te16 '
      print*, te11,te12,te13,te14,te15,te16

      ENDIF

      PROP6S=YXDOT(6,6)
      PROP7S=YXDOT(7,7)
      PROP8S=YXDOT(8,8)
      PROP9S=YXDOT(9,9)
      PROP10S=YXDOT(10,10)
      PROP11S=YXDOT(11,11)
CsB
C      Print*, ' PROP6S, PROP7S, PROP8S, PROP9S, PROP10S, PROP11S '
C      Print*,   PROP6S, PROP7S, PROP8S, PROP9S, PROP10S, PROP11S

C COLOR FACTOR FOR Q QBAR --> V1 V2 G
      SPIN=1.0/2.0/2.0
      COLORA=4.0/3.0/3.0
      CONST1=(4.0*PI*ALFAS)*(GWEAK2*SWS)**2
      CONST2=2.0**2*YXP(4,M_u1)*YXP(4,M_u2)
      CONST=CONST1*CONST2

      RETURN
      END

      SUBROUTINE FINDPOL(P,PX,PY,PLONG,XM)
      IMPLICIT NONE
      REAL*8 CTH,STH,Phi,CPHI,SPHI,PI,PMAG,XM,TEST
      INTEGER I
      REAL*8 P,PX,PLONG,PY
      DIMENSION P(4),PX(4),PLONG(4),PY(4)
C HERE WE USE CARTESIAN COMPONENTS FOR THE POLARIZATION VECTORS
C
      PI=4.0*ATAN(1.0D0)
      PMAG=DSQRT(P(1)**2+P(2)**2+P(3)**2)
      CTH=P(3)/P(4)
      STH=DSQRT(1.0D0-CTH**2)
      IF(STH.NE.0.D0) THEN
        CPHI=P(1)/PMAG/STH
        SPHI=P(2)/PMAG/STH
      ELSE
        CPHI=0.D0
        SPHI=0.D0
      ENDIF
C FOR PX
      PX(4)=0.0D0
      PX(1)=CPHI**2*CTH+SPHI**2
      PX(2)=CPHI*SPHI*(CTH-1.0D0)
      PX(3)=-STH*CPHI
C FOR PY
      PY(4)=0.0D0
      PY(1)=CPHI*SPHI*(CTH-1.0D0)
      PY(2)=SPHI**2*CTH + CPHI**2
      PY(3)=-STH*SPHI
C FOR PLONG
      XM=0.D0
      DO I=1,4
        PLONG(I)=0.D0
      ENDDO
      TEST=P(4)**2-PMAG**2
C FOR MASSLESS BOSONS
      IF(DABS(TEST).LT.1.D-4) RETURN
      XM=DSQRT(TEST)
      PLONG(4)=PMAG/XM
      PLONG(1)=P(4)*STH*CPHI/XM
      PLONG(2)=P(4)*STH*SPHI/XM
      PLONG(3)=P(4)*CTH/XM

      RETURN
      END


      SUBROUTINE AMPS_QQVVG(PI,ALFAS,XMZ,GWEAK2,SWS,A_L,A_R,
     >B_L,B_R,AMPSQUARE)
      IMPLICIT NONE
      REAL*8 PI,ALFAS,XMZ,GWEAK2,SWS,A_L,A_R,B_L,B_R,AMPSQUARE

C THE AMPLITUDE FOR
C     q(u1) + qbar(u2) --> V1(v1) V2(v2) g(g)

      INTEGER*4 M_u1,M_u2,M_v1,M_v2,M_g,M_zo
      COMMON/QQAAG_INIT/ M_u1,M_u2,M_v1,M_v2,M_g,M_zo
      REAL*8 YXP,YXETA,YXPT,YXDELR,YXPP
      COMMON/SECCONV/ YXP(4,30),YXETA(30),YXPT(30),YXDELR(30,30),
     >YXPP(30)
      REAL*8 R2,PROP6S,PROP7S,PROP8S,PROP9S,PROP10S,PROP11S
      COMMON/QVGINIT2/ R2,PROP6S,PROP7S,PROP8S,PROP9S,PROP10S,PROP11S
      REAL*8 SPIN,COLORA,CONST1,CONST2,CONST
      COMMON/QVGINIT3/ SPIN,COLORA,CONST1,CONST2,CONST

      COMPLEX*16 BSSSSSK
      EXTERNAL BSSSSSK
      COMPLEX*16 TM1,TM2,TM3,TM4,TM5,TM6,AMP
      REAL*8 AMPS,AMPSL,AMPSR
      INTEGER LU1,LV1,LV2,LG,LU1M

C 1 FOR '-', 2 FOR '+' HELICITY
      AMPSL=0.D0
      AMPSR=0.D0
      AMPS=0.D0

C FOR DI-PHOTON PAIR PRODUCTION
C AMPSQUARE FOR LU1=2 SHOULD BE THE SAME AS LU1=1
C ONLY TRANSVERSE POLARIZATIONS

      DO 998 LU1=1,2
        IF(LU1.EQ.1) THEN
          LU1M=2
        ELSEIF(LU1.EQ.2) THEN
          LU1M=1
        ENDIF
CCPY INPOSING SYMMETRIES FOR DI-PHOTON PRODUCTIONS
        IF(LU1.EQ.2) GOTO 998

        DO 999 LV1=21,22
        DO 999 LV2=24,25
        DO 999 LG=27,28

         TM1=BSSSSSK(2,LU1,LV2,LU1M,7,LU1,LV1,LU1M,6,LU1,LG,LU1M,1,LU1)
         TM1=-TM1/PROP6S/PROP7S
         TM2=BSSSSSK(2,LU1,LV2,LU1M,7,LU1,LG,LU1M,8,LU1,LV1,LU1M,1,LU1)
         TM2=-TM2/PROP8S/PROP7S
         TM3=BSSSSSK(2,LU1,LG,LU1M,9,LU1,LV2,LU1M,8,LU1,LV1,LU1M,1,LU1)
         TM3=TM3/PROP8S/PROP9S
         TM4=BSSSSSK(2,LU1,LV1,LU1M,10,LU1,LV2,LU1M,6,LU1,LG,LU1M,1,LU1)
         TM4=-TM4/PROP6S/PROP10S
         TM5=BSSSSSK(2,LU1,LV1,LU1M,10,LU1,LG,LU1M,11,LU1,LV2,LU1M,1,
     >       LU1)
         TM5=-TM5/PROP10S/PROP11S
         TM6=BSSSSSK(2,LU1,LG,LU1M,9,LU1,LV1,LU1M,11,LU1,LV2,LU1M,1,LU1)
         TM6=TM6/PROP9S/PROP11S

         AMP=TM1+TM2+TM3+TM4+TM5+TM6
CsB To test the matrix element term by term
C         Amp = TM3
C         AmpS = Conjg(Amp)*TM2
         AMPS=CONJG(AMP)*AMP

CDUMP (to check Ward identities)
C      PRINT*,'   '
C      PRINT*,' LV1,LV2,LG,LU1,AMPS '
C      PRINT*,LV1,LV2,LG,LU1,AMPS

        IF(LU1.EQ.1) THEN
           AMPSL=AMPSL+AMPS*(A_L*B_L)**2
        ELSEIF(LU1.EQ.2) THEN
           AMPSR=AMPSR+AMPS*(A_R*B_R)**2
        ENDIF

999     CONTINUE
998   CONTINUE

CDUMP
C      PRINT*,' AMPSL,AMPSR = ',AMPSL,AMPSR

C FOR DI-PHOTON PAIR PRODUCTION

      AMPSR=AMPSL
      AMPSQUARE=(AMPSL+AMPSR)*SPIN*COLORA*CONST

      RETURN
      END

cC -------------------------------------------------------------------------
c      Subroutine TestAmp (TM1,TM2,TM3,TM4,TM5,TM6)
cC -------------------------------------------------------------------------
c      Implicit Double Precision (a-h,o-z)
c
cCsB Check the matrix element term by term
c         dummy(1) = TM1
c         dummy(2) = TM2
c         dummy(3) = TM3
c         dummy(4) = TM4
c         dummy(5) = TM5
c         dummy(6) = TM6
c         Do ii = 1, 6
c           Amp = dummy(ii)
c         AMPS=CONJG(AMP)*AMP
c
c      Return
c      End

C -------------------------------------------------------------------------
        FUNCTION YXDOT(J1,J2)
C -------------------------------------------------------------------------
	IMPLICIT NONE
	REAL*8 YXDOT
	INTEGER*4 J1,J2
        REAL*8 YXP,YXETA,YXPT,YXDELR,YXPP
        COMMON/SECCONV/ YXP(4,30),YXETA(30),YXPT(30),YXDELR(30,30),
     >  YXPP(30)

	YXDOT=YXP(4,J1)*YXP(4,J2)-YXP(3,J1)*YXP(3,J2)-
     >        YXP(2,J1)*YXP(2,J2)-YXP(1,J1)*YXP(1,J2)

	RETURN
	END

C -------------------------------------------------------------------------
	FUNCTION BSK(V1,H1,V3,H3,V2,H2)
C -------------------------------------------------------------------------
	IMPLICIT NONE
	COMPLEX*16 BSK
	INTEGER*4 V1,H1,V3,H3,V2,H2
	REAL*8 A,B,PSL
	DIMENSION PSL(4)
	CHARACTER*1 HELI1,HELI2,HELI3
        REAL*8 YXP,YXETA,YXPT,YXDELR,YXPP
        COMMON/SECCONV/ YXP(4,30),YXETA(30),YXPT(30),YXDELR(30,30),
     >  YXPP(30)
	INTEGER*4 J1

	 IF(H1.EQ.1) THEN
	HELI1='-'
	 ELSEIF(H1.EQ.2) THEN
	HELI1='+'
	 ENDIF
	 IF(H2.EQ.1) THEN
	HELI2='-'
	 ELSEIF(H2.EQ.2) THEN
	HELI2='+'
	 ENDIF
	 IF(H3.EQ.1) THEN
	HELI3='-'
	 ELSEIF(H3.EQ.2) THEN
	HELI3='+'
	 ENDIF
	 DO 10 J1=1,4
	PSL(J1)=YXP(J1,V3)
10     CONTINUE
	 CALL YBRA_SL_KET(V1,HELI1,PSL,HELI3,V2,HELI2,A,B)
	 BSK=CMPLX(A,B)

	RETURN
	END

	FUNCTION BSSK(V1,H1,V3,H3,V4,H4,V2,H2)
	IMPLICIT NONE
	COMPLEX*16 BSSK
	INTEGER*4 V1,H1,V3,H3,V2,H2,V4,H4
	COMPLEX*16 BSK
	EXTERNAL BSK
	INTEGER*4 zo
	zo=30
	 BSSK=BSK(V1,H1,V3,H3,zo,1)*BSK(zo,1,V4,H4,V2,H2)
     >      + BSK(V1,H1,V3,H3,zo,2)*BSK(zo,2,V4,H4,V2,H2)
	RETURN
	END

	FUNCTION BSSSK(V1,H1,V3,H3,V4,H4,V5,H5,V2,H2)
	IMPLICIT NONE
	COMPLEX*16 BSSSK
	INTEGER*4 V1,H1,V3,H3,V2,H2,V4,H4,V5,H5
	COMPLEX*16 BSK,BSSK
	EXTERNAL BSK,BSSK
	INTEGER*4 zo
	zo=30
	 BSSSK=BSK(V1,H1,V3,H3,zo,1)*BSSK(zo,1,V4,H4,V5,H5,V2,H2)
     >      +BSK(V1,H1,V3,H3,zo,2)*BSSK(zo,2,V4,H4,V5,H5,V2,H2)

	RETURN
	END

C -------------------------------------------------------------------------
	Function BSSSSK(V1,H1,V3,H3,V4,H4,V5,H5,v6,h6,V2,H2)
C -------------------------------------------------------------------------
	IMPLICIT NONE
	COMPLEX*16 BSSSSK
	INTEGER*4 V1,H1,V3,H3,V2,H2,V4,H4,V5,H5,v6,h6
	COMPLEX*16 BSK,BSSSK
	EXTERNAL BSK,BSSSK
	INTEGER*4 zo
	zo=30
	 BSSSSK=BSK(V1,H1,V3,H3,zo,1)*BSSSK(zo,1,V4,H4,V5,H5,v6,h6,V2,H2)
     >       +BSK(V1,H1,V3,H3,zo,2)*BSSSK(zo,2,V4,H4,V5,H5,v6,h6,V2,H2)

	RETURN
	END

C -------------------------------------------------------------------------
	Function BSSSSSK(V1,H1,V3,H3,V4,H4,V5,H5,v6,h6,v7,h7,V2,H2)
C -------------------------------------------------------------------------
	IMPLICIT NONE
	COMPLEX*16 BSSSSSK
	INTEGER*4 V1,H1,V3,H3,V2,H2,V4,H4,V5,H5,v6,h6,v7,h7
	COMPLEX*16 BSK,BSSSSK
	EXTERNAL BSK,BSSSSK
	INTEGER*4 zo
	zo=30
	 BSSSSSK = BSK(V1,H1,V3,H3,zo,1)*
     *           BSSSSK(zo,1,V4,H4,V5,H5,v6,h6,v7,h7,V2,H2) +
     +           BSK(V1,H1,V3,H3,zo,2)*
     *           BSSSSK(zo,2,V4,H4,V5,H5,v6,h6,v7,h7,V2,H2)

	RETURN
	END

C -------------------------------------------------------------------------
	SUBROUTINE YBRA_SL_KET(IMOM1,HELI1,P,CHI,IMOM2,HELI2,AA,BB)
C -------------------------------------------------------------------------
	IMPLICIT NONE
	INTEGER*4 IMOM1,IMOM2
	CHARACTER*1 HELI1,HELI2,CHI
	REAL*8 P
	DIMENSION P(4)
	REAL*8 AA,BB,A1,B1,C1,D1,A2,B2,C2,D2,A,B,C,D,E,F,G,H

	CALL YBRA(IMOM1,HELI1,A1,B1,C1,D1)
	CALL YKET(IMOM2,HELI2,A2,B2,C2,D2)
	CALL YSLASH(P,CHI,A,B,C,D,E,F,G,H)

	AA=(A1*A*A2-A1*B*B2+A1*C*C2-A1*D*D2)
     >  -(B1*B*A2+B1*A*B2+B1*D*C2+B1*C*D2)
     >  +(C1*E*A2-C1*F*B2+C1*G*C2-C1*H*D2)
     >  -(D1*F*A2+D1*E*B2+D1*H*C2+D1*G*D2)

	BB=(B1*A*A2-B1*B*B2+B1*C*C2-B1*D*D2
     >  +A1*B*A2+A1*A*B2+A1*D*C2+A1*C*D2
     >  +D1*E*A2-D1*F*B2+D1*G*C2-D1*H*D2
     >  +C1*F*A2+C1*E*B2+C1*H*C2+C1*G*D2)

	RETURN
	END

C -------------------------------------------------------------------------
	SUBROUTINE YKET(IMOM,HELI,A,B,C,D)
C -------------------------------------------------------------------------
	IMPLICIT NONE
	INTEGER*4 IMOM,I
	CHARACTER*1 HELI
	REAL*8 A,B,C,D
	REAL*8 P
	DIMENSION P(4)
	REAL*8 PP,DENO,SMALL

        REAL*8 YXP,YXETA,YXPT,YXDELR,YXPP
        COMMON/SECCONV/ YXP(4,30),YXETA(30),YXPT(30),YXDELR(30,30),
     >  YXPP(30)
C
C zo has to be 30
C
	IF(IMOM.EQ.30) THEN
	 IF(HELI.EQ.'+') THEN
	A=1.0
	B=0.0
	C=0.0
	D=0.0
	 ELSEIF(HELI.EQ.'-') THEN
	A=0.0
	B=0.0
	C=1.0
	D=0.0
	 ENDIF
	 RETURN
	ENDIF

	SMALL=1.D-8
	DO 10 I=1,4
	 P(I)=YXP(I,IMOM)
10    CONTINUE
	PP=SQRT(P(1)**2+P(2)**2+P(3)**2)
	IF(HELI.EQ.'+') THEN
	 IF(PP+P(3).LT.SMALL) THEN
	A=0.0
	B=0.0
	C=-1.0
	D=0.0
	 ELSE
	DENO=SQRT(2.0*PP*(PP+P(3)))
	A=(PP+P(3))/DENO
	B=0.0
	C=P(1)/DENO
	D=P(2)/DENO
	 ENDIF
	ELSE IF(HELI.EQ.'-') THEN
	 IF(PP+P(3).LT.SMALL) THEN
	A=1.0
	B=0.0
	C=0.0
	D=0.0
	 ELSE
	DENO=SQRT(2.0*PP*(PP+P(3)))
	C=(PP+P(3))/DENO
	D=0.0
	A=-P(1)/DENO
	B=P(2)/DENO
	 ENDIF
	ELSE
	 PRINT*,' ERROR IN YKET OR YBRA'
	ENDIF

	RETURN
	END

	SUBROUTINE YBRA(IMOM,HELI,A,B,C,D)
	IMPLICIT NONE
	INTEGER*4 IMOM
	CHARACTER*1 HELI
	REAL*8 A,B,C,D

	CALL YKET(IMOM,HELI,A,B,C,D)
	B=-B
	D=-D
	RETURN
	END

	SUBROUTINE YSLASH(P,CHI,A,B,C,D,E,F,G,H)
	IMPLICIT NONE
	CHARACTER*1 CHI
	REAL*8 P,A,B,C,D,E,F,G,H
	DIMENSION P(4)

	IF(CHI.EQ.'+') THEN
	 A=P(4)+P(3)
	 C=P(1)
	 D=-P(2)
	 G=P(4)-P(3)
	ELSE IF(CHI.EQ.'-') THEN
	 A=P(4)-P(3)
	 C=-P(1)
	 D=P(2)
	 G=P(4)+P(3)
	ELSE
	 PRINT*,' ERROR IN YSLASH'
	ENDIF
	B=0.0
	H=0.0
	E=C
	F=-D

	RETURN
	END

CsB   End of qq -> vvG double check.

CsB___Running alpha_em from PYTHIA
C --------------------------------------------------------------------------
      FUNCTION PYALEM(Q2)
C --------------------------------------------------------------------------
C...Double precision and integer declarations.
C      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      Real*8 PYALEM,Q2
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      REAL PARU(200)
      INTEGER MSTU(200)
C      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
C      SAVE /PYDAT1/

C...Calculate real part of photon vacuum polarization.
C...For leptons simplify by using asymptotic (Q^2 >> m^2) expressions.
C...For hadrons use parametrization of H. Burkhardt et al.
C...See R. Kleiss et al, CERN 89-08, vol. 3, pp. 129-131.
CsB___See also Phys.Lett.B356:398-403,1995

      PARU(1)=3.1415927
      PARU(101)=1./137.04
      MSTU(101)=1
      PARU(104)=1.
      PARU(103)=1./128.8
      AEMPI=PARU(101)/(3D0*PARU(1))
      IF(MSTU(101).LE.0.OR.Q2.LT.2D-6) THEN
        RPIGG=0D0
      ELSEIF(MSTU(101).EQ.2.AND.Q2.LT.PARU(104)) THEN
        RPIGG=0D0
      ELSEIF(MSTU(101).EQ.2) THEN
        RPIGG=1D0-PARU(101)/PARU(103)
      ELSEIF(Q2.LT.0.09D0) THEN
        RPIGG=AEMPI*(13.4916D0+LOG(Q2))+0.00835D0*LOG(1D0+Q2)
      ELSEIF(Q2.LT.9D0) THEN
        RPIGG=AEMPI*(16.3200D0+2D0*LOG(Q2))+
     &  0.00238D0*LOG(1D0+3.927D0*Q2)
      ELSEIF(Q2.LT.1D4) THEN
        RPIGG=AEMPI*(13.4955D0+3D0*LOG(Q2))+0.00165D0+
     &  0.00299D0*LOG(1D0+Q2)
      ELSE
        RPIGG=AEMPI*(13.4955D0+3D0*LOG(Q2))+0.00221D0+
     &  0.00293D0*LOG(1D0+Q2)
      ENDIF

C...Calculate running alpha_em.
      PYALEM=PARU(101)/(1D0-RPIGG)
      PARU(108)=PYALEM

      RETURN
      END

C --------------------------------------------------------------------------
      Subroutine TestPYALEM
C --------------------------------------------------------------------------

      Real*8 Q_V, PYalEM, Pi
      PARAMETER (Pi = 3.14159 26535 89793 24D0)

CCPY      Do Q_V=60.,100.,5
      Do I=1,10
        Q_V=60.0+(I-1)*5.0
        eEM   = Sqrt(4.d0*Pi*PYALEM(Q_V**2))
        Print*, Q_V, 1./eEM, 1./PYALEM(Q_V**2)
      EndDo

      Stop
      End
C --------------------------------------------------------------------------
      Subroutine TestALPHAS
C --------------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      Double Precision ALPHAS
      External ALPHAS
      
      Do ii=1,10
        AMU = ii*20.
        TMP = ALPHAS (AMU)
        Print*,  AMU, TMP
      End Do
        
      Stop
      Return
      End

C --------------------------------------------------------------------------
      SUBROUTINE TestPDF
C --------------------------------------------------------------------------
      IMPLICIT NONE
      
      integer nflmax
      parameter(nflmax=5)
      REAL*8 Q,x1x2,q1(-nflmax:nflmax)
      REAL*8 Ctq6Pdf ! PDF

      INTEGER*4 i

C   The function Ctq5Pdf (Iparton, X, Q)
C   returns the parton distribution inside the proton for parton [Iparton] 
C   at [X] Bjorken_X and scale [Q] (GeV) in PDF set [Iset].
C   Iparton  is the parton label (5, 4, 3, 2, 1, 0, -1, ......, -5)
C                            for (b, c, s, d, u, g, u_bar, ..., b_bar),
      Q = 30.d0
      Do i=1,10
        x1x2 = i/10.d0
        Q1(1) = x1x2*Ctq6Pdf ( 0, x1x2, Q)
        Print*, ' x, gluon: ', x1x2, Q1(1)
        End Do
      Stop

      RETURN
      END

C ==========================================================================
C* $Id: dilog64.F,v 1.1.1.1 1996/04/01 15:02:05 mclareni Exp $
C*
C* $Log: dilog64.F,v $
C* Revision 1.1.1.1  1996/04/01 15:02:05  mclareni
C* Mathlib gen
C*
C*
C#include "gen/pilot.h"
C#if defined(CERNLIB_DOUBLE)
C      FUNCTION DDILOG(X)
C#include "gen/imp64.inc"
C#endif
C#if !defined(CERNLIB_DOUBLE)
C      FUNCTION RDILOG(X)
C#endif
CsB The dilog fn is defined as
C   Li2(x) = - Integral_0^1 Log(1 - tx)/t dt = Sum_{k=1}^{Infinity} x^k/k^2
C   Note the negative sign before the integral.
C   This way Li2(x) >= 0 for any (real) x.
C   The Li2 function returns values which are consistent with the above
C   definition. This was cross checked using Mathematica.
CsB As of April 2005, there's another Li_2 function in this file called fLi2.  
C   fLi2 also takes complex arguments.  Li2(x) should be replaced by fLi2.
C --------------------------------------------------------------------------
      Function Li2(X)
C --------------------------------------------------------------------------
      Implicit Double Precision (a-h,o-z)
      Real*8 Li2
      DIMENSION C(0:19)

      PARAMETER (Z1 = 1, HF = Z1/2)
      PARAMETER (PI = 3.14159 26535 89793 24D0)
CsB      PARAMETER (PI3 = PI**2/3, PI6 = PI**2/6, PI12 = PI**2/12)
      PARAMETER (PI3 = PI**2/3.d0, PI6 = PI**2/6.d0, PI12 = PI**2/12.d0)

      DATA C( 0) / 0.42996 69356 08136 97D0/
      DATA C( 1) / 0.40975 98753 30771 05D0/
      DATA C( 2) /-0.01858 84366 50145 92D0/
      DATA C( 3) / 0.00145 75108 40622 68D0/
      DATA C( 4) /-0.00014 30418 44423 40D0/
      DATA C( 5) / 0.00001 58841 55418 80D0/
      DATA C( 6) /-0.00000 19078 49593 87D0/
      DATA C( 7) / 0.00000 02419 51808 54D0/
      DATA C( 8) /-0.00000 00319 33412 74D0/
      DATA C( 9) / 0.00000 00043 45450 63D0/
      DATA C(10) /-0.00000 00006 05784 80D0/
      DATA C(11) / 0.00000 00000 86120 98D0/
      DATA C(12) /-0.00000 00000 12443 32D0/
      DATA C(13) / 0.00000 00000 01822 56D0/
      DATA C(14) /-0.00000 00000 00270 07D0/
      DATA C(15) / 0.00000 00000 00040 42D0/
      DATA C(16) /-0.00000 00000 00006 10D0/
      DATA C(17) / 0.00000 00000 00000 93D0/
      DATA C(18) /-0.00000 00000 00000 14D0/
      DATA C(19) /+0.00000 00000 00000 02D0/

      IF(X .EQ. 1) THEN
        H=PI6
      ELSEIF(X .EQ. -1) THEN
       H=-PI12
      ELSE
       T=-X
       IF(T .LE. -2) THEN
        Y=-1/(1+T)
        S=1
        A=-PI3+HF*(LOG(-T)**2-LOG(1+1/T)**2)
       ELSEIF(T .LT. -1) THEN
        Y=-1-T
        S=-1
        A=LOG(-T)
        A=-PI6+A*(A+LOG(1+1/T))
       ELSE IF(T .LE. -HF) THEN
        Y=-(1+T)/T
        S=1
        A=LOG(-T)
        A=-PI6+A*(-HF*A+LOG(1+T))
       ELSE IF(T .LT. 0) THEN
        Y=-T/(1+T)
        S=-1
        A=HF*LOG(1+T)**2
       ELSE IF(T .LE. 1) THEN
        Y=T
        S=1
        A=0
       ELSE
        Y=1/T
        S=-1
        A=PI6+HF*LOG(T)**2
       ENDIF
       H=Y+Y-1
       ALFA=H+H
       B1=0
       B2=0
       DO 1 I = 19,0,-1
       B0=C(I)+ALFA*B1-B2
       B2=B1
    1  B1=B0
       H=-(S*(B0-H*B2)+A)
      ENDIF
C#if defined(CERNLIB_DOUBLE)
C      DDILOG=H
C#endif
C#if !defined(CERNLIB_DOUBLE)
C       DILOG=H
C#endif
      Li2 = H
      RETURN
      END

C --------------------------------------------------------------------------
      Subroutine TestDilog
C --------------------------------------------------------------------------
      Real*8 Li2, a,b,c,x

      a = Li2(0.d0)
      Print*, 0.d0, '   ', a
      a = Li2(.5d0)
      Print*, .5d0, '   ', a
      a = Li2(.314159d0)
      Print*, .314159d0, '   ', a
      a = Li2(.88888d0)
      Print*, .88888d0, '   ', a
      a = Li2(.99999d0)
      Print*, .99999d0, '   ', a
      a = Li2(1.d0)
      Print*, 1.d0, '   ', a
      b = .2d0
      c = .3d0
      a = Li2((c-b)/c)
      Print*, (c-b)/c, '   ', a

      Stop
      End

CsB The code below was implemented from the code HJet.
C   It calculates h1 h2 -> H X -> Z Z X
C   The formulae were taken from NPB297('88)221
C   The total cross section agrees with HJet - 
C   after the HJet output is multiplied by the H -> Z Z branching ratio
C   (Br(H -> Z Z) ~ 0.3 for mH = 250GeV).
C   Warning: something is not quite right with the m_top = infinity option.
C   I think the problem comes form HJet. As long as it's not fixed
C   I hard-wired the m_top finite option in resbos.f.
C ==========================================================================
      Subroutine pp_H0X_ZZX (WTA2)
C --------------------------------------------------------------------------
      Implicit Double Precision (a-h,o-z)

      COMMON/STAND0/XMW,XMZ,XMW2,XMZ2,AEM,GAMW,GAMZ,GAMW2,GAMZ2
      COMMON/STAND5/XMCHARM,XMCHARM2,XMHIGGS,GAMH,GAMH2
      COMMON/STAND2/XLU,XLD,XRU,XRD,QU,QD,GFERMI,XLE,XRE
      COMMON/STU/SH,TH,UH

      Data PI/3.14159 26535 89793 24D0/

C      Print*, ' ScaPro(1,1),ScaPro(2,2)=', ScaPro(1,1),ScaPro(2,2)
C      Print*, ' ScaPro(3,3),ScaPro(4,4)=', ScaPro(3,3),ScaPro(4,4)
C      Print*, ' ScaPro(5,5)            =', ScaPro(5,5)
C      Stop

CsB___Construct the parton level Mandelstam variables
C     Momenta (pLab) are passed through the Common Block in function ScaPro
C     pLab(i,j) is the i-th coordinate of the momentum p(j)
C      If (Debug) Print*, ' Scalar products: '
      sH =  2.d0*ScaPro(1,2)
      tH = -2.d0*ScaPro(2,5)
      uH = -2.d0*ScaPro(1,5)
C      If (Debug) Print*, ii, jj, SP(ii,jj)

CsB___Higgs invariant mass (same as YM12 - can be taken from COMMON)
      Q = Sqrt(2.d0*(ScaPro(3,4) + xMz**2))
C      Print*, ' Q =', Q

CsB___Weight for the production: parton parton' -> H parton''
      Call PDFMEH(WTProd)
C      Print*, ' WTProd =', WTProd
C      Stop

      WTA2 =
     =       WTProd *
C_____Higgs propagator
     *       1.d0/((Q**2-xMHiggs**2)**2 + (xMHiggs*GamH)**2) *
C_____H -> Z Z amplitude squared and averaged
     *       Sqrt(2.d0)*GFermi *
     *       (12.d0*xMz**4 - 4.d0*(Q*xMz)**2 + Q**4)

      Return
      End

C --------------------------------------------------------------------------
      SUBROUTINE PDFMEH(PPB)
C --------------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)

      INTEGER SHIFT_CUTS,SWITCH_ORDER,SWITCH_APPROX
      REAL*8 MW,MZ,MT,MH

CsB___Common blocks communicating with the HJET part of the code
change AAH everywhere !!!
      COMMON/AAHJet/PI,Q,QCDLHJet,NFHJet
      COMMON/PAR/CME,MW,MZ,MT,SIN2W,ALFA,ALFAWK
      COMMON/SWSH/SHIFT_CUTS,SWITCH_ORDER,SWITCH_APPROX
c check for this in HJet routines!!!      COMMON/MOM/PLAB(4,8)
      COMMON/BINS/PPB1,PPB2,PPB3
      COMMON/HIG/MH,WIDH,BRZ

CsB___Common blocks communicating with both HJET part and the rest
      COMMON/STU/SH,TH,UH

CsB___Common blocks communicating with the rest (non-HJet part) of the code
      COMMON/AA/W,QCDL,X1X2MIN,NJET,N_PRIMARY_PART,NF,I_STR,IBEAM,NPART,
     &  NEVENT,NORDER,N_SECOND_PART
      COMMON/STAND0/XMW,XMZ,XMW2,XMZ2,AEM,GAMW,GAMZ,GAMW2,GAMZ2
      COMMON/STAND1/SWS,CWS,XMBOT,XMTOP,XMB2,XMT2,GWEAK2
      integer nflmax
      parameter (nflmax=5)
      COMMON/STRUCTURE/q1(-nflmax:nflmax),q2(-nflmax:nflmax)
      COMMON/STAND5/XMCHARM,XMCHARM2,XMHIGGS,GAMH,GAMH2
      COMMON/PMF/X1,X2
      Common / HQApprox / iHQApprox

      Data PI/3.14159 26535 89793 24D0/

CsB___Redundant
      SHIFT_CUTS   = 0

CsB___SWITCH_ORDER = 1 and 2 for H and H + j respectively
C     Note: for SWITCH_ORDER = 1 the phase-space must be for two particle!
      SWITCH_ORDER = 2
CsB___SWITCH_APPROX = 1 and 2 for finite m_t and infinite m_t respectively
      SWITCH_APPROX = iHQApprox

CsB___Higgs invariant mass (same as YM12 - can be taken from COMMON)
      Q = Sqrt(2.d0*(ScaPro(3,4) + xMZ**2))
C      Print*, ' Q = ', Q

CsB___QDC parameters
      QCDLHJet = QCDL
      nFHJet   = NF

CsB___Ecm of the collision
      CME = W

CsB___Masses
      MW = xMw
      MZ = xMz
      MT = xMtop
      MH = xMHiggs

CsB___Decay parameters
      WIDH = GAMH
      BRZ  = 1.d0/3.d0

CsB___EW parameters
      SIN2W  = SWS
      GFMU   = 1.166389D-5
      ALFAWK = DSQRT(2.D0)*GFMU*MW**2/Pi
      ALFA   = ALFAWK*SIN2W

C -- CALL MATRIX ELEMENTS AND ADD INITIAL STATE
C    COLOR AND SPIN AVERAGES, ALSO FACTOR OF 1/2
C    FOR FINAL STATE IDENTICAL PARTICLES WHERE
C    APPROPRIATE.
         IF(SWITCH_ORDER.EQ.1) THEN ! LO (GG contribution)
           IF(SWITCH_APPROX.EQ.1) THEN
             CALL SQAMPH1(GGH)
               GGH = GGH/256.
               QGH = 0.D0
               QQH = 0.D0
            ELSE
             CALL SQAMPH1_ASY(GGH)
               GGH = GGH/256.
               QGH = 0.D0
               QQH = 0.D0
           ENDIF
          ELSE                      ! NLO (GG, qG and QQ contributions)
              IF(SWITCH_APPROX.EQ.1) THEN
                CALL SQAMPH2GG(GGH)
                CALL SQAMPH2QQ(QQH)
                CALL SQAMPH2QG(QGH)
                   GGH = GGH/256.
                   QGH = QGH/96.
                   QQH = QQH/36.
               ELSE
                CALL SQAMPH2GG_ASY(GGH)
                CALL SQAMPH2QQ_ASY(QQH)
                CALL SQAMPH2QG_ASY(QGH)
                   GGH = GGH/256.
                   QGH = QGH/96.
                   QQH = QQH/36.
              ENDIF
         ENDIF
C
        PPB1 = 0.D0
        PPB2 = 0.D0
        PPB3 = 0.D0
C
C        PD1(1)=(U1+SEA1)/X1
C        PD1(2)=(D1+SEA1)/X1
C        PD1(3)=STR1/X1
C        PD1(4)=CHM1/X1
C        PD1(5)= 0.D0
C        PD1(6)=SEA1/X1
C        PD1(7)=SEA1/X1
C        PD1(8)=STR1/X1
C        PD1(9)=CHM1/X1
C        PD1(10)= 0.D0
C        PD1(11)=G1/X1
C
C        PPB1=PPB1+(PD1(1)*PD2(6)+PD1(2)*PD2(7)+PD1(3)*PD2(8)+
C     .            PD1(4)*PD2(9)+PD1(5)*PD2(10) + PD2(1)*PD1(6)+
C     .            PD2(2)*PD1(7)+PD2(3)*PD1(8)+PD2(4)*PD1(9)+
C     .            PD2(5)*PD1(10))*QQH
C        PPB2=PPB2+((PD1(1)+PD1(2)+PD1(3)+PD1(4)+PD1(5)+PD1(6)+PD1(7)+
C     .            PD1(8)+PD1(9)+PD1(10))*PD2(11) +
C     .             (PD2(1)+PD2(2)+PD2(3)+PD2(4)+PD2(5)+PD2(6)+PD2(7)+
C     .            PD2(8)+PD2(9)+PD2(10))*PD1(11))*QGH
C        PPB3 = PPB3 + (PD1(11)*PD2(11))*GGH

CsB___PDF's were already called from MAKE_EVENT

C uU + dD + sS +... + Uu + Dd +...
        PPB1 = PPB1 + ( Q1(1)*Q2(-1)+Q1(2)*Q2(-2)+Q1(3)*Q2(-3)+
     .                  Q1(4)*Q2(-4)+Q1(5)*Q2(-5) +
     .                  Q2(1)*Q1(-1)+Q2(2)*Q1(-2)+Q2(3)*Q1(-3)+
     .                  Q2(4)*Q1(-4)+Q2(5)*Q1(-5)) * QQH

C (u+d+s+...+U+D+S+...)G + G(u+d+...+U+D+...)
        PPB2 = PPB2 + QGH*( (Q1(1)+Q1(2)+Q1(3)+Q1(4)+Q1(5) +
     .                   Q1(-1)+Q1(-2)+Q1(-3)+Q1(-4)+Q1(-5) ) * Q2(0) +
     .                  (Q2(1)+Q2(2)+Q2(3)+Q2(4)+Q2(5) +
     .                   Q2(-1)+Q2(-2)+Q2(-3)+Q2(-4)+Q2(-5))*Q1(0) )

C GG
        PPB3 = PPB3 + Q1(0)*Q2(0)*GGH

      PPB = PPB1 + PPB2 + PPB3

CsB___These were found in agreement with HJet
CC      Print*, ' Q1(0)/P_U,Q2(0)/P_U =', Q1(0)/1d8,Q2(0)/1d8
C      Print*, ' x1,x2 =', x1,x2
C      Print*, ' Q1(0)/P_U/x1,Q2(0)/P_U/x2 =', Q1(0)/x1/1d8,Q2(0)/x2/1d8
CC      Print*, ' GGH =', GGH
C      Print*, ' PPB3/P_U/x1/x2 =', PPB3/x1/x2/1d16
CC      Stop

      RETURN
      END

C --------------------------------------------------------------------------
      SUBROUTINE SQAMPH1(GGH)
C --------------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 MH,MW,MZ,MT
      COMPLEX*16 A1
      COMMON/PARTONXs/X1,X2
      COMMON/AAHJet/PI,Q,QCDL,NF
      COMMON/PAR/CME,MW,MZ,MT,SIN2W,ALFA,ALFAWK
      COMMON/HIG/MH,WIDH,BRZ
      EXTERNAL A1

CsB___NPB297('88)221(A.3)
      COUPLFAC = ALFAWK*ALPHAS(Q)**2*(MH**4/MW**2)/PI
      GGH = COUPLFAC*CDABS(A1(MH**2))**2
C      Print*, ' ALFAWK, Q, ALPHAS(Q) ', ALFAWK, Q, ALPHAS(Q)
C      PRINT*, ' GGH, COUPLFAC =', GGH, COUPLFAC
C      Print*, ' MH =', MH
C      PRINT*, ' CDABS(A1(MH**2)) =', CDABS(A1(MH**2))
      RETURN
      END

C --------------------------------------------------------------------------
      SUBROUTINE SQAMPH1_ASY(GGH)
C --------------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 MH,MW,MZ,MT
      COMMON/PARTONXs/X1,X2
      COMMON/AAHJet/PI,Q,QCDL,NF
      COMMON/PAR/CME,MW,MZ,MT,SIN2W,ALFA,ALFAWK
      COMMON/HIG/MH,WIDH,BRZ

CsB___NPB297('88)221(A.4)
       COUPLFAC = ALFAWK*ALPHAS(Q)**2*(MH**4/MW**2)/PI
          GGH = COUPLFAC*(2.D0/3.D0)**2
c           PRINT*, 'GGH,COUPLFAC'
c           PRINT*, GGH,COUPLFAC
      RETURN
      END

cC --------------------------------------------------------------------------
c      SUBROUTINE GPAIR(SIG,GX,GY)
cC --------------------------------------------------------------------------
c      IMPLICIT REAL*8(A-H,O-Z)
c      DIMENSION QRAN(10)
c      DATA TWOPI,ROOT2/6.2831853D0,1.4142136D0/
c
c      CALL ARAN9(QRAN,5)
c      R=SIG*ROOT2*DSQRT(-DLOG(QRAN(1)))
c      PHI=TWOPI*QRAN(2)
c      GX=R*DCOS(PHI)
c      GY=R*DSIN(PHI)
c      RETURN
c      END

C --------------------------------------------------------------------------
      SUBROUTINE SQAMPH2GG(GGHG)
C --------------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 MH,MW,MZ,MT
      COMPLEX*16 A2,A4
      COMMON/AAHJet/PI,Q,QCDL,NF
      COMMON/PAR/CME,MW,MZ,MT,SIN2W,ALFA,ALFAWK
      COMMON/HIG/MH,WIDH,BRZ
      COMMON/STU/SH,TH,UH
      EXTERNAL A2,A4

CsB___NPB297('88)221(2.2)
      COUPLFAC = ALFAWK*ALPHAS(Q)**3*4.*3.*8.*(MH**8/MW**2)
      GGHG = COUPLFAC*
     *       ( CDABS(A2(SH,TH,UH))**2 + CDABS(A2(UH,SH,TH))**2 
     +       + CDABS(A2(TH,UH,SH))**2 + CDABS(A4(SH,TH,UH))**2 )/
     /       SH/TH/UH
C
      RETURN
      END

C --------------------------------------------------------------------------
      SUBROUTINE SQAMPH2GG_ASY(GGHG)
C --------------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 MH,MW,MZ,MT
      COMMON/AAHJet/PI,Q,QCDL,NF
      COMMON/PAR/CME,MW,MZ,MT,SIN2W,ALFA,ALFAWK
      COMMON/HIG/MH,WIDH,BRZ
      COMMON/STU/SH,TH,UH

CsB___NPB297('88)221(2.3)
       COUPLFAC = ALFAWK*ALPHAS(Q)**3*4.*3.*8./MW**2/9.
       GGHG = COUPLFAC*( MH**8 + SH**4 + TH**4 + UH**4)/SH/TH/UH
C         PRINT*, 'GGHG,SH,TH,UH'
C         PRINT*, GGHG,SH,TH,UH
      RETURN
      END

C --------------------------------------------------------------------------
      SUBROUTINE SQAMPH2QQ(QQHG)
C --------------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 MH,MW,MZ,MT
      COMPLEX*16 A5
      COMMON/PARTONXs/X1,X2
      COMMON/AAHJet/PI,Q,QCDL,NF
      COMMON/PAR/CME,MW,MZ,MT,SIN2W,ALFA,ALFAWK
      COMMON/HIG/MH,WIDH,BRZ
      COMMON/STU/SH,TH,UH

CsB___NPB297('88)221(2.6)
      COUPLFAC = ALFAWK*ALPHAS(Q)**3*8.*(MH**4/2./MW**2)
      QQHG = COUPLFAC*(UH**2 + TH**2)/SH/(UH + TH)**2*
     .                                CDABS(A5(SH,TH,UH))**2
      RETURN
      END

C --------------------------------------------------------------------------
      SUBROUTINE SQAMPH2QQ_ASY(QQHG)
C --------------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 MH,MW,MZ,MT
      COMMON/PARTONXs/X1,X2
      COMMON/AAHJet/PI,Q,QCDL,NF
      COMMON/PAR/CME,MW,MZ,MT,SIN2W,ALFA,ALFAWK
      COMMON/HIG/MH,WIDH,BRZ
      COMMON/STU/SH,TH,UH

CsB___NPB297('88)221(2.7)
      COUPLFAC = ALFAWK*ALPHAS(Q)**3*8.*2./(9.*MW**2)
      QQHG = COUPLFAC*(UH**2 + TH**2)/SH

      RETURN
      END

C --------------------------------------------------------------------------
      SUBROUTINE SQAMPH2QG(QGHQ)
C --------------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 MH,MW,MZ,MT
      COMPLEX*16 A5
      COMMON/PARTONXs/X1,X2
      COMMON/AAHJet/PI,Q,QCDL,NF
      COMMON/PAR/CME,MW,MZ,MT,SIN2W,ALFA,ALFAWK
      COMMON/HIG/MH,WIDH,BRZ
      COMMON/STU/SH,TH,UH

CsB___NPB297('88)221(2.9)
      COUPLFAC = ALFAWK*ALPHAS(Q)**3*8.*(MH**4/2./MW**2)
      QGHQ = -COUPLFAC*(UH**2 + SH**2)/TH/(UH + SH)**2*
     .                            CDABS(A5(TH,SH,UH))**2
C        PRINT*, 'A51'
C        PRINT*,A5(TH,SH,UH)
      RETURN
      END

C --------------------------------------------------------------------------
      SUBROUTINE SQAMPH2QG_ASY(QGHQ)
C --------------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 MH,MW,MZ,MT
      COMMON/PARTONXs/X1,X2
      COMMON/AAHJet/PI,Q,QCDL,NF
      COMMON/PAR/CME,MW,MZ,MT,SIN2W,ALFA,ALFAWK
      COMMON/HIG/MH,WIDH,BRZ
      COMMON/STU/SH,TH,UH

CsB___NPB297('88)221 -(2.7) s<->t
      COUPLFAC = ALFAWK*ALPHAS(Q)**3*8.*2./(9.*MW**2)
      QGHQ = -COUPLFAC*(UH**2 + SH**2)/TH

      RETURN
      END

C --------------------------------------------------------------------------
      DOUBLE COMPLEX FUNCTION A1(X)
C --------------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 W2,w2m
      REAL*8 Mf,MW,MZ,MT,MH
      COMMON/PAR/CME,MW,MZ,MT,SIN2W,ALFA,ALFAWK
      COMMON/HIG/MH,WIDH,BRZ
      EXTERNAL W2

CsB___NPB297('88)221(A.2)
      Mf = MT
      A1 = (Mf**2/MH**2)*(4. - W2(X)*(1. - 4.*Mf**2/MH**2))
c              w2m = W2(X)
c        PRINT*, 'A1,w2m'
c        PRINT*, A1,w2m
      RETURN
      END

C --------------------------------------------------------------------------
      DOUBLE COMPLEX FUNCTION A2(X,Y,Z)
C --------------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 b2
      EXTERNAL b2

CsB___NPB297('88)221(A.9)
        A2 = b2(X,Y,Z) + b2(X,Z,Y)
ccc      PRINT*, 'A2'
ccc      PRINT*, A2
      RETURN
      END

C --------------------------------------------------------------------------
      DOUBLE COMPLEX FUNCTION A4(X,Y,Z)
C --------------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 b4
      EXTERNAL b4

CsB___NPB297('88)221(A.8)
      A4 = (b4(X,Y,Z) + b4(Z,X,Y) + b4(Y,Z,X))
ccc      PRINT*, 'A4,X,Y,Z'
ccc      PRINT*, A4,X,Y,Z
      RETURN
      END

C --------------------------------------------------------------------------
      DOUBLE COMPLEX FUNCTION b2(X,Y,Z)
C --------------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 Mf,MW,MZ,MT,MH
      COMPLEX*16 W1,W2,W3
      COMMON/PAR/CME,MW,MZ,MT,SIN2W,ALFA,ALFAWK
      COMMON/HIG/MH,WIDH,BRZ
      EXTERNAL W1,W2,W3

CsB___NPB297('88)221(A.11)
      Mf = MT
      b2 = (Mf**2/MH**4)*( X*(Z-X)/(Z+X) + 2.*Y*Z*(Z+2.*X)*
     .    (W1(Y) - W1(MH**2))/(Z+X)**2 +(Mf**2 - X/4.)*(W2(X)/2. +
     .               W2(MH**2)/2.- W2(Y) + W3(X,Y,Z,MH**2) )  +
     .    X**2*( 2.*Mf**2/(Z+X)**2 - 1./2./(Z+X))*
     .    ( W2(Y) - W2(MH**2)) + Y*Z*( W2(MH**2) - 2.*W2(Y))/2./X +
     .    (X - 12.*Mf**2 - 4.*Y*Z/X)/8.*W3(Y,X,Z,MH**2) )

      RETURN
      END

C --------------------------------------------------------------------------
      DOUBLE COMPLEX FUNCTION b4(X,Y,Z)
C --------------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 Mf,MW,MZ,MT,MH
      COMPLEX*16 W2,W3
      COMMON/PAR/CME,MW,MZ,MT,SIN2W,ALFA,ALFAWK
      COMMON/HIG/MH,WIDH,BRZ
      EXTERNAL W1,W2,W3

CsB___NPB297('88)221(A.10)
      Mf = MT
      b4 = Mf**2/MH**2*(-2.D0/3.D0+(-1.D0/4.D0 + Mf**2/MH**2)*
     .           (W2(Y) - W2(MH**2) + W3(X,Y,Z,MH**2) ) )
cc       PRINT*, 'b4,MT,MH,Mf'
cc       PRINT*, b4,MT,MH,Mf
      RETURN
      END

C --------------------------------------------------------------------------
      DOUBLE COMPLEX FUNCTION A5(X,Y,Z)
C --------------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 Mf,MW,MZ,MT,MH
      COMPLEX*16 W1,W2
      COMMON/PAR/CME,MW,MZ,MT,SIN2W,ALFA,ALFAWK
      COMMON/HIG/MH,WIDH,BRZ

CsB___NPB297('88)221(A.14)
      Mf = MT
      A5 = Mf**2/MH**2*( 4. + 4.*X/(Y+Z)*(W1(X) - W1(MH**2)) +
     .              (1.- 4.*Mf**2/(Y+Z))*(W2(X) - W2(MH**2)) )
ccc        PRINT*, 'A52'
ccc        PRINT*,A5
      RETURN
      END

C --------------------------------------------------------------------------
      DOUBLE COMPLEX FUNCTION W1(X)
C --------------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 Mf,MW,MZ,MT
      COMPLEX*16 I
      COMMON/AAHJet/PI,Q,QCDL,NF
      COMMON/PAR/CME,MW,MZ,MT,SIN2W,ALFA,ALFAWK

CsB___NPB297('88)221(A.19)
      Mf = MT
      I = DCMPLX(0.D0, 1.D0)
         IF(X.LT.0.D0) THEN
            W1 = 2.*DSQRT(1. - 4.*Mf**2/X)*DASINH( DSQRT(-X)/2./Mf)
          ELSEIF(X.LT.4.*Mf**2.AND.X.GT.0.D0) THEN
            W1 = 2.*DSQRT(- 1. + 4.*Mf**2/X)*DASIN( DSQRT(X)/2./Mf)
          ELSEIF(X.GT.4.*Mf**2) THEN
            W1 = DSQRT(1. - 4.*Mf**2/X)*(2.*DACOSH( DSQRT(X)/2./Mf) -
     .                      PI*I)
         ENDIF
C           PRINT*, 'W1'
C           PRINT*, W1
      RETURN
      END

C --------------------------------------------------------------------------
      DOUBLE COMPLEX FUNCTION W2(X)
C --------------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 Mf,MW,MZ,MT
      COMPLEX*16 I
      COMMON/AAHJet/PI,Q,QCDL,NF
      COMMON/PAR/CME,MW,MZ,MT,SIN2W,ALFA,ALFAWK

CsB___NPB297('88)221(A.20)
      Mf = MT
      I = DCMPLX(0.D0, 1.D0)
         IF(X.LT.0.D0) THEN
            W2 = 4.*DASINH( DSQRT(-X)/2./Mf)**2
          ELSEIF(X.LT.4.*Mf**2.AND.X.GT.0.D0) THEN
            W2 = - 4.*DASIN( DSQRT(X)/2./Mf)**2
          ELSEIF(X.GT.4.*Mf**2) THEN
            W2 = 4.*DACOSH( DSQRT(X)/2./Mf)**2 - PI**2 - 4.*PI*I
     .               *DACOSH( DSQRT(X)/2./Mf)
         ENDIF
C           PRINT*, 'W2'
C           PRINT*, W2
      RETURN
      END

C --------------------------------------------------------------------------
      DOUBLE COMPLEX FUNCTION W3(X,Y,Z,W)
C --------------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 I3
      EXTERNAL I3

CsB___NPB297('88)221(A.17)
      W3 = I3(X,Y,Z,W) - I3(X,Y,Z,X) - I3(X,Y,Z,Z)
C           PRINT*, 'W3'
C           PRINT*, W3
      RETURN
      END

C --------------------------------------------------------------------------
      DOUBLE COMPLEX FUNCTION I3(X,Y,Z,W)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 Mf,MW,MZ,MT,phi,theta,beta,gamma,f1,f2
      COMPLEX*16 I,SPENCE,SPENCET
      COMMON/AAHJet/PI,Q,QCDL,NF
      COMMON/PAR/CME,MW,MZ,MT,SIN2W,ALFA,ALFAWK
      DATA SMALL/1.D-8/

CsB___NPB297('88)221(A.21-22)
      Mf = MT
      I = DCMPLX(0.D0, 1.D0)
      beta = (1. + DSQRT(1. + 4.*Y*Mf**2/X/Z))/2.
C
      IF(W.LT.0.D0) THEN
        gamma = (1. + DSQRT(1. - 4.*Mf**2/W))/2.
        f2 = beta + gamma - 1.
        I3 = 2.*(-SPENCE(DCMPLX(gamma/f2)) +
     .              SPENCE( DCMPLX((gamma -1)/f2 )) +
     .              SPENCE(DCMPLX(( beta-gamma)/beta)) -
     .              SPENCE(DCMPLX((beta-gamma)/(beta-1))) +
     .        (DLOG(beta)**2 - DLOG(beta-1)**2)/2. + DLOG(gamma)*
     .       DLOG(f2/beta)+DLOG(gamma-1)*DLOG((beta-1)/f2))/(2.*beta-1)
C
      ELSEIF(W.GT.0.D0.AND.W.LT.4.*Mf**2) THEN
           AA = DSQRT(- 1. + 4.*Mf**2/W)
           RR = DSQRT( (AA**2 + 1.)/(AA**2 + (2.*beta - 1.)**2 ))
           PHIRES = RR*(AA**2 + 2.*beta - 1.)/(1 + AA**2)
           THETARES = RR*(AA**2 - 2.*beta + 1.)/(1 + AA**2)
            IF(DABS(PHIRES).GT.1.D0) THEN
           PHIRES = (DABS(PHIRES) - SMALL)*DABS(PHIRES)/PHIRES
                PRINT*, 'PHIRES'
                PRINT*, PHIRES
              ELSEIF(DABS(THETARES).GT.1.D0) THEN
           THETARES = (DABS(THETARES) - SMALL)*DABS(THETARES)/THETARES
                PRINT*, 'THETARES'
                PRINT*, THETARES
            ENDIF
C
           phi = DACOS( PHIRES)
           theta = DACOS( THETARES )
C
        I3 = 2.*(2.*SPENCET(RR,THETARES) - 2.*SPENCET(RR,PHIRES) +
     .           (phi - theta)*(phi + theta - PI) )/(2.*beta - 1.)
C
              ELSEIF(W.GT.4.*Mf**2) THEN
           gamma = (1. + DSQRT(1. - 4.*Mf**2/W))/2.
            f2 = beta + gamma - 1.
        I3 = 2.*(-SPENCE(DCMPLX(gamma/f2)) +
     .                SPENCE( DCMPLX((gamma -1)/f2 )) +
     .                SPENCE(DCMPLX(gamma/(gamma-beta))) -
     .                SPENCE(DCMPLX((gamma-1)/(gamma-beta))) +
     .            DLOG(gamma/(1.-gamma))*DLOG(f2/(beta-gamma)) -
     .            PI*I*DLOG(f2/(beta-gamma)) )/(2.*beta -1)
              ENDIF
      RETURN
      END

      DOUBLE COMPLEX FUNCTION SPENCET(X,Y)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 YP,YN,SPENCE
      EXTERNAL SPENCE
C
C       This is for the function,
C               -1/2 Integrate[ ln[1 - 2*Y*Z + Z^2]/Z, {Z, 0, X} ]
C
CsB___NPB297('88)221(End of Abstract)
               YP = Y + CDSQRT(DCMPLX(Y)**2 - 1.)
               YN = Y - CDSQRT(DCMPLX(Y)**2 - 1.)
ccc             PRINT*, 'X,Y,YP,YN'
ccc             PRINT*,X,Y,YP,YN
C
          SPENCET = (SPENCE(X/YP) + SPENCE(X/YN))/2.
      RETURN
      END

C     This subroutine calculates the Spence function for any complex
C     value of the argument with double precision. 15 significant
C     figures are guaranteed for any ABS(X) < 1.E6
      DOUBLE COMPLEX FUNCTION SPENCE(X)
      DOUBLE COMPLEX X,XX
      DOUBLE PRECISION BER(10),P,PI,B
      COMMON/BER/BER,PI,P

C ****** These lines should preferably be placed outside the subroutine *****
C     COMMON/BER/BER,PI,P
      PI=4.D0*DATAN(1.D0)
      P=PI*PI/6.D0
      CALL BERNOULLI(BER)
C ***************************************************************************

             NFLAG=0
             XX=X
             LAND=1
             CALL COUNTRY(XX,NFLAG)

         IF (NFLAG.EQ.0) THEN
             XX=1.-X
             LAND=2
             CALL COUNTRY(XX,NFLAG)
         END IF

         IF (NFLAG.EQ.0) THEN
             XX=1./X
             LAND=3
             CALL COUNTRY(XX,NFLAG)
         END IF

         IF (NFLAG.EQ.0) THEN
             XX=1./(1.-X)
             LAND=4
             CALL COUNTRY(XX,NFLAG)
         END IF

      CALL SERIES(XX,SPENCE)

      IF (LAND.EQ.1) RETURN

      IF (LAND.EQ.2) THEN
           SPENCE=-SPENCE+P
           B=ABS(1.-X)
           IF(B.GT.0.D0) SPENCE=SPENCE-LOG(X)*LOG(1.-X)
      END IF

      IF (LAND.EQ.3) THEN
          SPENCE=-SPENCE-P-0.5*LOG(-X)*LOG(-X)
      END IF

      IF (LAND.EQ.4) THEN
          SPENCE=SPENCE+2.*P-LOG(X)*LOG(1.-X)+0.5*(LOG(X-1.))**2
      END IF

      RETURN
      END

      SUBROUTINE COUNTRY(X,NFLAG)
      DOUBLE COMPLEX X,Z
      DOUBLE PRECISION R
      Z=1.-X
      IF(ABS(Z).EQ.0.) RETURN
      Z=-LOG(Z)
      R=ABS(Z)
      IF (R.LE.DLOG(2.D0)) NFLAG=1
      RETURN
      END

      SUBROUTINE SERIES(X,SP)
      DOUBLE COMPLEX SP,TP,X,Z,ZZ
      DOUBLE PRECISION BER(10),P,PI,FAC,A
      COMMON/BER/BER,PI,P
      FAC=0.5
      Z=-LOG(1.-X)
      SP=Z*(1.-0.25*Z)
      ZZ=Z*Z*Z
      DO 1 I=1,10
      A=FAC*BER(I)/DFLOAT(2*I+1)
      TP=SP+A*ZZ
      SP=TP
      ZZ=ZZ*Z*Z
  1   FAC=FAC/DFLOAT(2*(I+1)*(2*I+1))
  2   FORMAT(2E27.18,/,2E27.18)
      RETURN
      END

      SUBROUTINE BERNOULLI(BER)
      DOUBLE PRECISION BER(10)
      DIMENSION N(10),M(10)

      N(1)=1
      N(2)=-1
      N(3)=1
      N(4)=-1
      N(5)=5
      N(6)=-691
      N(7)=7
      N(8)=-3617
      N(9)=43867
      N(10)=-174611

      M(1)=6
      M(2)=30
      M(3)=42
      M(4)=30
      M(5)=66
      M(6)=2730
      M(7)=6
      M(8)=510
      M(9)=798
      M(10)=330

      DO 2 J=1,10
  2   BER(J)=DFLOAT(N(J))/DFLOAT(M(J))
      END
C
C
      DOUBLE PRECISION FUNCTION DASINH(X)
      IMPLICIT REAL*8(A-H,O-Z)
        DASINH=DLOG( X + DSQRT(X**2 + 1) )
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION DACOSH(X)
      IMPLICIT REAL*8(A-H,O-Z)
        DACOSH=DLOG( X + DSQRT(X**2 - 1) )
      RETURN
      END

C
c      SUBROUTINE ARAN9(QRAN,ND)
c      REAL*8 QRAN(10)
c      COMMON/SEEDVAX/NUM1
c      DO 2 I=1,ND
c    1 QRAN(I)=RAN(NUM1)
c      IF(QRAN(I).LE.0.0) GOTO 1
c      IF(QRAN(I).GE.1.0) GOTO 1
c    2 CONTINUE
c      RETURN
c      END

C ==========================================================================
C----------------------------------------------------------------------
      Function ENFS_old(Q_V)
C----------------------------------------------------------------------
C_____Factor which takes into account the quark masses in the loop

      Implicit None
      Real*8 ENFS_Old,ENF1,ENF2,AL,SAL,SQAL,FQL,ANF,BNF,Q_V,Pi
      Integer IQ

      REAL*8 
     &SWS,CWS,XMBOT,XMTOP,XMB2,XMT2,GWEAK2,
     &XMCHARM,XMCHARM2,XMHIGGS,GAMH,GAMH2

      COMMON/STAND1/SWS,CWS,XMBOT,XMTOP,XMB2,XMT2,GWEAK2
      COMMON/STAND5/XMCHARM,XMCHARM2,XMHIGGS,GAMH,GAMH2

      Data PI/3.14159 26535 89793 24D0/

C_______Testing the Mq -> infinity limit
C        XMCHARM = 100000.
C        XMBOT   = 100000.
C        XMTOP   = 100000.

        ENF1 = 0.d0
        ENF2 = 0.d0
        DO IQ=1,3
          IF(IQ.EQ.1)AL=XMCHARM/Q_V
          IF(IQ.EQ.2)AL=XMBOT/Q_V
          IF(IQ.EQ.3)AL=XMTOP/Q_V
          SAL=AL**2
          IF(AL.LE..5) THEN
            SQAL=SQRT(.25-SAL)
            FQL=2.*LOG((.5+SQAL)/AL)
            ANF=-PI**2/2.+.5*FQL**2
            BNF=PI*FQL
          ELSE
            ANF=-2.*(ASIN(.5/AL))**2
            BNF=0.0
          END IF
          ENF1=ENF1+ANF*SAL*(4.*SAL-1.)+2.*SAL
          ENF2=ENF2+BNF*SAL*(4.*SAL-1.)
        END DO
        ENFS_Old=(ENF1**2+ENF2**2) ! 9/4 * A^2 of NPB359('91)283 [see (2.2)]
C        Print*, ' ENFS = ',ENFS ! ENFS goes to 1 when m_q -> infinity
C        ENFS=1.d0

        Return
        End

C----------------------------------------------------------------------
      Function ENFS(Q)
C----------------------------------------------------------------------
C_____Factor which takes into account the quark masses in the loop

      Implicit None

      REAL*8
     &SWS,CWS,XMBOT,XMTOP,XMB2,XMT2,GWEAK2,
     &XMCHARM,XMCHARM2,XMHIGGS,GAMH,GAMH2

      COMMON/STAND1/SWS,CWS,XMBOT,XMTOP,XMB2,XMT2,GWEAK2
      COMMON/STAND5/XMCHARM,XMCHARM2,XMHIGGS,GAMH,GAMH2

      Integer i
      Complex A,f
      Real*8 ENFS,Q,t,t1,t2,t3,Pi

      Data Pi/3.14159 26535 89793 24D0/

C_______Testing the mq -> infinity limit
C        qmas_c  = 1000. ! and Do i=1,1 !

        A = (0.d0, 0.d0)
        Do i=1,3
          If (i.EQ.1) t = Q**2/XMCHARM**2/4.d0
          If (i.EQ.2) t = Q**2/XMB2/4.d0
          If (i.EQ.3) t = Q**2/XMT2/4.d0
          If (t.LE.1.) then
            f = DASin(Sqrt(t))**2
          Else 
            t1 = DSqrt(1.d0-1.d0/t)
            t2 = 1.d0 + t1
            t3 = 1.d0 - t1
            f  = -1.d0/4.d0*(Log(t2/t3) - (0.d0,1.d0)*Pi)**2
          End If
          A = A + 2.d0*(t + (t-1.d0)*f)/t**2
        End Do

        ENFS = CDAbs(3.d0/4.d0*A)**2

C        Print*, ' ENFS = ',ENFS ! ENFS goes to nf^2 when m_q -> infinity
C        Stop

        Return
        End

C----------------------------------------------------------------------
      Subroutine SetKappa
C----------------------------------------------------------------------
C_____Effective coupling of the Higgs to gluons in the heavy top quark limit
C     see Eq.(13) of hep-ph/9611272 by Kramer-Laenen-Spira

      Implicit Real*8 (a-h,o-z)
      Real*8 Scale

      COMMON/STAND1/SWS,CWS,XMBOT,XMTOP,XMB2,XMT2,GWEAK2
      COMMON/STAND5/XMCHARM,XMCHARM2,XMHIGGS,GAMH,GAMH2

      Real*8 KappaNLO,KappaSoft
      Common / KappaHiggs / KappaNLO,KappaSoft

      EXTERNAL ALPHAS


      Data Pi/3.14159 26535 89793 24D0/

CsB___# of light flavors
      nf = 5

CsB___Define alpha_s^(5), Eq.(10) of hep-ph/9611272
      as6mt = ALPHAS(XMTOP + .1d0)
      as5mt = as6mt*(1.d0 + 11.d0/72.d0*(as6mt/Pi)**2)
      as6mh = ALPHAS(XMHiggs) 
      as5mh = as6mh*(1.d0 + 11.d0/72.d0*(as6mh/Pi)**2)
C      Print*, as6mt,as5mt,as6mh,as5mh

CCPY!      Kappa = 1.d0 +
CCPY!     +        11.d0/2.d0 * as5mt/Pi +
CCPY!     +        (3866.d0-201.d0)*nf/144.d0*(as5mt/Pi)**2 +
CCPY!     +        (153.d0-19.d0*nf)/(33.d0-2.d0*nf)*(as5mh-as5mt)/Pi
CCPY!
CCPY! nOTE: THE FACTOR (11.d0/2.d0) INSIDE KappaNLO IS ALREADY 
C included INSIDE FH0Del() which is the 
C delta(1-tau) pieces in the MSBar scheme

      KappaNLO=1.d0 +
     +        11.d0/2.d0 * as5mh/Pi

      KappaSoft=KappaNLO+11.d0/2.d0 * (as5mt/Pi-as5mh/Pi)+
     +        (3866.d0-201.d0)*nf/144.d0*(as5mt/Pi)**2 +
     +        (153.d0-19.d0*nf)/(33.d0-2.d0*nf)*(as5mh-as5mt)/Pi

      Print*, ' XMHiggs,KappaNLO,KappaSoft =',
     >XMHiggs,KappaNLO,KappaSoft

      Return
      End

CCPY ------------------------------------------------------------------------
      Subroutine TestHWIDTH
CCPY------------------------------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Do k=45,1005,5
!      Do k=105,105
        Q  = k
        GAM = HWIDTH(Q)
        Print*, ' Q, HWIDTH = ', Q, GAM
      EndDo

      Stop
      End

CCPY------------------------------------------------------------------------
CCPY ------------------------------------------------------------------------
      Function HWIDTH(Q)
CCPY ------------------------------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Real*8 HWIDTH,Q
      Double Precision MhTbl(96)
      Parameter (nP = 5)
      Dimension GAMTBL(96),XA(nP),YA(nP)

CCPY___Values of SM higgs masses and total decay width of SM Higgs boson 
C     from HDecay VERSION 3.4.
C     This could be updated by EW K-factors from hep-ph/0404071
      Data (MhTbl(ii), ii=1,96) 
     // 
     -  50.,60.,70.,80.,90.,100.,110.,
     -  120.,130.,140.,150.,160.,170.,180.,
     -  190.,200.,210.,220.,230.,240.,250.,
     -  260.,270.,280.,290.,300.,310.,320.,
     -  330.,340.,350.,360.,370.,380.,390.,
     -  400.,410.,420.,430.,
     -  440.,450.,460.,470.,480.,490.,500.,
     -  510.,520.,530.,540.,550.,560.,570.,
     -  580.,590.,600.,610.,620.,630.,640.,
     -  650.,660.,670.,680.,690.,700.,710.,
     -  720.,730.,740.,750.,760.,770.,780.,
     -  790.,800.,810.,820.,830.,840.,850.,
     -  860.,870.,880.,890.,900.,910.,920.,
     -  930.,940.,950.,960.,970.,980.,990.,
     -  1000.
     //
      Data (GAMTBL(ii), ii=1,96)
     // 
     -  0.1478E-02  , 0.1724E-02  , 0.1968E-02  , 
     -  0.2213E-02  , 0.2463E-02  , 0.2738E-02  , 
     -  0.3109E-02  , 0.3762E-02  , 0.5120E-02  , 
     -  0.8246E-02  , 0.1687E-01  , 0.7744E-01  , 
     -  0.3840  , 0.6285  , 1.038  ,
     -  1.426  , 1.838  , 2.297  , 2.813  ,
     -  3.394  , 4.046  , 4.772  ,
     -  5.578  , 6.466  , 7.441  ,
     -  8.506  , 9.665  , 10.92  ,
     -  12.28  , 13.76  , 15.60  ,
     -  17.96  , 20.55  , 23.32  ,
     -  26.24  , 29.30  , 32.50  ,
     -  35.84  , 39.30  , 42.91  ,
     -  46.65  , 50.54  , 54.56  ,
     -  58.74  , 63.07  , 67.56  ,
     -  72.21  , 77.03  , 82.02  ,
     -  87.19  , 92.55  , 98.10  ,
     -  103.8  , 109.8  , 116.0  ,
     -  122.3  , 128.9  , 135.8  ,
     -  142.8  , 150.2  , 157.7  ,
     -  165.6  , 173.7  , 182.0  ,
     -  190.7  , 199.7  , 208.9  ,
     -  218.5  , 228.4  , 238.7  ,
     -  249.2  , 260.2  , 271.5  ,
     -  283.2  , 295.2  , 307.7  ,
     -  320.6  , 333.9  , 347.6  ,
     -  361.8  , 376.5  , 391.7  ,
     -  407.3  , 423.5  , 440.2  ,
     -  457.4  , 475.2  , 493.6  ,
     -  512.6  , 532.2  , 552.5  ,
     -  573.4  , 595.0  , 617.3  ,
     -  640.3  , 664.1
     // 

CCPY___Where is Q in MhTbl?
      Call Locate(MhTbl,96,Q,LQ)
!      Print*, ' Q, LQ = ', Q, LQ

CCPY___If Q is not within [MhTbl(1),MhTbl(96)]
      If ( (MhTbl(1)-Q).GE.0.d0 .Or. (Q-MhTbl(96)).GE.0.d0 ) then
        Print*, ' Extrapolation used in calculating H -> gamma gamma! '
        Print*, ' Q, LQ = ', Q, LQ
      EndIf

CCPY___If close to the ends
      If ( LQ.LT.2)  LQ = 2
      If ( LQ.GT.93) LQ = 93
!      Print *, ' LQ = ', LQ

CCPY___Cut small parts of the tables
      Do i=1,nP
        XA(i) = MhTbl(LQ+i-2)
        YA(i) = GAMTBL(LQ+i-2)
      EndDo

CCPY___Interpolate on the small parts
      Call Polint(XA,YA,nP,Q,GAMAtQ,ErrBr)
      HWIDTH = GAMAtQ

      Return
      End

CCPY ------------------------------------------------------------------------
CsB ------------------------------------------------------------------------
      Subroutine TestBrH2AA
CsB ------------------------------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Do k=45,1005,5
!      Do k=105,105
        Q  = k
        Br = BrH2AA(Q)
        Print*, ' Q, Br(H -> photon photon) = ', Q, Br
      EndDo

      Stop
      End

CsB ------------------------------------------------------------------------
      Function BrH2AA(Q)
CsB ------------------------------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Real*8 BrH2AA,Q
      Double Precision MhTbl(96)
      Parameter (nP = 5)
      Dimension BrTbl(96),XA(nP),YA(nP)

CsB___Values of SM higgs masses and H -> photon photon branching ratios
C     from HDecay VERSION 3.4.
C     This could be updated by EW K-factors from hep-ph/0404071
      Data (MhTbl(ii), ii=1,96) 
     // 
     -  50.,60.,70.,80.,90.,100.,110.,
     -  120.,130.,140.,150.,160.,170.,180.,
     -  190.,200.,210.,220.,230.,240.,250.,
     -  260.,270.,280.,290.,300.,310.,320.,
     -  330.,340.,350.,360.,370.,380.,390.,
     -  400.,410.,420.,430.,
     -  440.,450.,460.,470.,480.,490.,500.,
     -  510.,520.,530.,540.,550.,560.,570.,
     -  580.,590.,600.,610.,620.,630.,640.,
     -  650.,660.,670.,680.,690.,700.,710.,
     -  720.,730.,740.,750.,760.,770.,780.,
     -  790.,800.,810.,820.,830.,840.,850.,
     -  860.,870.,880.,890.,900.,910.,920.,
     -  930.,940.,950.,960.,970.,980.,990.,
     -  1000.
     //
      Data (BrTbl(ii), ii=1,96)
     // 
     -  0.2960E-03  , 0.4427E-03  , 0.6276E-03  , 
     -  0.8568E-03  , 0.1137E-02  , 0.1472E-02  , 
     -  0.1832E-02  , 0.2118E-02  , 0.2176E-02  , 
     -  0.1917E-02  , 0.1388E-02  , 0.5534E-03  , 
     -  0.1514E-03  , 0.1025E-03  , 0.6756E-04  , 
     -  0.5287E-04  , 0.4363E-04  , 0.3683E-04  , 
     -  0.3151E-04  , 0.2721E-04  , 0.2368E-04  , 
     -  0.2074E-04  , 0.1828E-04  , 0.1619E-04  , 
     -  0.1440E-04  , 0.1288E-04  , 0.1156E-04  , 
     -  0.1042E-04  , 0.9424E-05  , 0.8557E-05  , 
     -  0.7352E-05  , 0.5908E-05  , 0.4726E-05  , 
     -  0.3781E-05  , 0.3032E-05  , 0.2436E-05  , 
     -  0.1960E-05  , 0.1579E-05  , 0.1272E-05  , 
     -  0.1023E-05  , 0.8218E-06  , 0.6589E-06 , 
     -  0.5271E-06 , 0.4200E-06 , 0.3337E-06 , 
     -  0.2644E-06 , 0.2092E-06 , 0.1657E-06 , 
     -  0.1318E-06 , 0.1060E-06 , 0.8703E-07 , 
     -  0.7368E-07 , 0.6504E-07 , 0.6032E-07 , 
     -  0.5889E-07 , 0.6025E-07 , 0.6381E-07 , 
     -  0.6924E-07 , 0.7623E-07 , 0.8450E-07 , 
     -  0.9373E-07 , 0.1037E-06 , 0.1144E-06 , 
     -  0.1255E-06 , 0.1370E-06 , 0.1486E-06 , 
     -  0.1604E-06 , 0.1723E-06 , 0.1841E-06 , 
     -  0.1958E-06 , 0.2074E-06 , 0.2189E-06 , 
     -  0.2301E-06 , 0.2411E-06 , 0.2518E-06 , 
     -  0.2622E-06 , 0.2724E-06 , 0.2822E-06 , 
     -  0.2917E-06 , 0.3009E-06 , 0.3098E-06 , 
     -  0.3182E-06 , 0.3263E-06 , 0.3341E-06 , 
     -  0.3416E-06 , 0.3487E-06 , 0.3546E-06 , 
     -  0.3612E-06 , 0.3680E-06 , 0.3737E-06 , 
     -  0.3792E-06 , 0.3843E-06 , 0.3891E-06 , 
     -  0.3936E-06 , 0.3978E-06 , 0.4017E-06 
     // 
CsB___Where is Q in MhTbl?
      Call Locate(MhTbl,96,Q,LQ)
!      Print*, ' Q, LQ = ', Q, LQ

CsB___If Q is not within [MhTbl(1),MhTbl(96)]
      If ( (MhTbl(1)-Q).GE.0.d0 .Or. (Q-MhTbl(96)).GE.0.d0 ) then
        Print*, ' Extrapolation used in calculating H -> gamma gamma! '
        Print*, ' Q, LQ = ', Q, LQ
      EndIf

CsB___If close to the ends
      If ( LQ.LT.2)  LQ = 2
      If ( LQ.GT.93) LQ = 93
!      Print *, ' LQ = ', LQ

CsB___Cut small parts of the tables
      Do i=1,nP
        XA(i) = MhTbl(LQ+i-2)
        YA(i) = BrTbl(LQ+i-2)
      EndDo

CsB___Interpolate on the small parts
      Call Polint(XA,YA,nP,Q,BrAtQ,ErrBr)
      BrH2AA = BrAtQ

      Return
      End

CsB ------------------------------------------------------------------------
CsB ------------------------------------------------------------------------
CsB The routines below generate the Z decays for the process: 
C   H -> Z Z -> 4 leptons.
C
C GENERATE THE DECAY OF Z'S FROM H --> Z Z --> E- E+ MU- MU+
C
C INSIDE THE SUBROUTINE WRESPH, WE ADDED THE FOLLOWING STEP.
C
C      IF(I_DECAY.EQ.1) THEN 
C        CALL HZZDECAY
C      ENDIF
C
C----------------------------------------------------------------------------
          
      SUBROUTINE HZZDECAY (P_Q,P_L1CS,P_L1,P_L2,
CsB___To double check boost from CS to lab frame
     ,P_X,P_Y,P_Z,STHE,CTHE,SPHI,CPHI,Q_V)
CsB___To double check |P_L1|,|P_L2|
!      SUBROUTINE HZZDECAY (P_Q,P_L1CS) !,P_L1,P_L2)
      
      IMPLICIT NONE
      INTEGER I1,I2,I3,I4,J1,J2,J3,J4,ii
      REAL*8 P_Q(5),P_L1CS(5),P_L2CS(5),TEMP,R_00,PI,RN,
     >       THEMU_1,PHIMU_1,EMU_1,THEMU_3,PHIMU_3,EMU_3,P_ZMU_1(5),
     >       P_ZMU_3(5),P_HMU_1(5),P_HMU_3(5),P_MU_1(5),P_MU_2(5),
     >       P_MU_3(5),P_MU_4(5),SEC_PLAB,P_L1(5),P_L2(5),PMU_1,PMU_3,
     >       MMU_1,MMU_2,MMU_3,MMU_4,Zero,One,Dummy,ScaPro,ScaProS
CsB___To double check |P_L1|,|P_L2|
!     REAL*8 PLAB,PCEM,WHAT
CsB___To double check boosts
      REAL*8 PMU_2,EMU_2,P_ZMU_2(5),P_HMU_2(5),
     >       PMU_4,EMU_4,P_ZMU_4(5),P_HMU_4(5)
CsB___To double check boost from CS to lab frame
      REAL*8 P_X(5),P_Y(5),P_Z(5),STHE,CTHE,SPHI,CPHI,Q_V,QT,MT
      REAL*8 SinTh1,CosTh1,SinPh1,CosPh1,SinTh2,CosTh2,SinPh2,CosPh2,tmp
      REAL*8 SinAl,CosAl,PhiCorrZZ
      
      Common /PHIMU3/PHIMU_3
      
      REAL*8 XMW,XMZ,XMW2,XMZ2,AEM,GAMW,GAMZ,GAMW2,GAMZ2,
     &       SWS,CWS,XMBOT,XMTOP,XMB2,XMT2,GWEAK2,
     &       XMCHARM,XMCHARM2,XMHIGGS,GAMH,GAMH2
      COMMON/STAND0/XMW,XMZ,XMW2,XMZ2,AEM,GAMW,GAMZ,GAMW2,GAMZ2
      COMMON/STAND1/SWS,CWS,XMBOT,XMTOP,XMB2,XMT2,GWEAK2
      COMMON/STAND5/XMCHARM,XMCHARM2,XMHIGGS,GAMH,GAMH2

      EXTERNAL RN,ScaProS
      COMMON/DECAY_MOMENTA/SEC_PLAB(4,50)
CsB___To double check |P_L1|,|P_L2|
!      COMMON/MOM/PLAB(4,12),PCEM(4,12),WHAT
      DATA PI/3.14159 26535 89793 24D0/
      DATA I1,I2,I3,I4,J1,J2,J3,J4/1,1,1,1,1,1,1,1/
      DATA Zero,One/0.0,1.0/

CsB___Tests
!      Print *,' Entered SUBROUTINE HZZDECAY '
C     Note: The argument of RN is dummy in this code, i.e. RN doesn't depend
!           on it.
!      Print *,' RN(1),RN(2),RN(3),RN(4) = ', RN(1),RN(2),RN(3),RN(4)
!      Stop

C (1) WE NEED THE 4-MOMENTUM OF Z_1 AND Z_2, WHICH ARE LABELLED AS 
C P_L1CS AND P_L2CS
CsB___The CS frame is an L1+L2 rest frame so
      P_L2CS(1) = -P_L1CS(1)
      P_L2CS(2) = -P_L1CS(2)
      P_L2CS(3) = -P_L1CS(3)
!      P_L2CS(4) =  P_L1CS(4)
!      P_L2CS(5) =  P_L1CS(5)
      P_L2CS(4) =  Sqrt(P_L2(5)**2 + 
     +                  P_L2CS(1)**2 + P_L2CS(2)**2 + P_L2CS(3)**2)
      P_L2CS(5) =  P_L2(5)
CsB___Test L1 and L2 masses
!      Print *,' P_L1CS(5),P_L2CS(5) = ', P_L1CS(5),P_L2CS(5)
      
C (2) DETERMINE THE POLARIZATIONS OF Z'S
      TEMP = 1.0 - 2.0*P_L1CS(4)**2/P_L1CS(5)**2
      R_00 = TEMP**2/(2.0+TEMP**2)

      IF(RN(I1).LT.R_00) THEN
C GENERATE LONGITUDINAL Z'S DECAY
!        Print *, ' HELICTY 0 OF Z`S ' 
C POLAR ANGLE AND ENERGY FOR E- FROM Z_1 --> E-(1) E+(2)
        CALL ZFFB(Zero,RN(I2),P_L1CS(5),1,THEMU_1,EMU_1)
C POLAR ANGLE AND ENERGY FOR MU- FROM Z_2 --> MU-(3) MU+(4)
        CALL ZFFB(Zero,RN(I3),P_L2CS(5),1,THEMU_3,EMU_3)
C GENERATE PHI ANGLES
        PHIMU_1=0.0
!        PHIMU_3=2.0*PI*RN(I4)
        PHIMU_3 = PhiCorrZZ(XMZ,XMHiggs)

      ELSE
C GENERATE TRANSVERSE Z'S DECAY
        IF(RN(J1).LT.0.5) THEN
!          Print *, ' HELICTY +1 OF Z`S ' 
!          Dummy = RN(J2)
!          Print *,' Test 1 -> ',Dummy,0.0,P_L1CS(5),0,THEMU_1,EMU_1
          CALL ZFFB(RN(J2),Zero,P_L1CS(5),0,THEMU_1,EMU_1)
          CALL ZFFB(RN(J3),Zero,P_L2CS(5),0,THEMU_3,EMU_3)
        ELSE
!          Print *, ' HELICTY -1 OF Z`S '
          CALL ZFFB(RN(J2),One,P_L1CS(5),0,THEMU_1,EMU_1)
          CALL ZFFB(RN(J3),One,P_L2CS(5),0,THEMU_3,EMU_3)
        ENDIF
C GENERATE PHI ANGLES
        PHIMU_1=0.0
        PHIMU_3=2.0*PI*RN(J4)

      ENDIF
     
CsB___Check EMU_i
!      Print *,' EMU_1, EMU_3 =',EMU_1, EMU_3
     
C (3) CONSTRUCT THE 4-MOMENTA FOR E- AND MU-, CALLED 
C P_ZMU_1 AND P_ZMU_3, WHICH ARE DEFINED IN THE REST FRAME OF Z.
CsB___E- and MU- are hard wired to be massless here.
      MMU_1 = 0.d0
      PMU_1 = Sqrt(EMU_1**2 - MMU_1**2)
      P_ZMU_1(1) = PMU_1 * Sin(THEMU_1) * Cos(PHIMU_1)
      P_ZMU_1(2) = PMU_1 * Sin(THEMU_1) * Sin(PHIMU_1)
      P_ZMU_1(3) = PMU_1 * Cos(THEMU_1)
      P_ZMU_1(4) = EMU_1
      P_ZMU_1(5) = MMU_1
      MMU_3 = 0.d0
      PMU_3 = Sqrt(EMU_3**2 - MMU_3**2)
      P_ZMU_3(1) = PMU_3 * Sin(THEMU_3) * Cos(PHIMU_3)
      P_ZMU_3(2) = PMU_3 * Sin(THEMU_3) * Sin(PHIMU_3)
      P_ZMU_3(3) = PMU_3 * Cos(THEMU_3)
      P_ZMU_3(4) = EMU_3
      P_ZMU_3(5) = MMU_3

CsB___Introduce P_ZMU_2 and P_ZMU_4
      MMU_2 = 0.d0
      EMU_2 = P_L1CS(5) - EMU_1
      PMU_2 = Sqrt(EMU_2**2 - MMU_2**2)
!!!   PMU_2 = Sqrt(EMU_1**2 - MMU_2**2)
      P_ZMU_2(1) = -PMU_1 * Sin(THEMU_1) * Cos(PHIMU_1)
      P_ZMU_2(2) = -PMU_1 * Sin(THEMU_1) * Sin(PHIMU_1)
      P_ZMU_2(3) = -PMU_1 * Cos(THEMU_1)
      P_ZMU_2(4) =  EMU_1
      P_ZMU_2(5) =  MMU_2
      MMU_4 = 0.d0
      EMU_4 = P_L2CS(5) - EMU_3
      PMU_4 = Sqrt(EMU_4**2 - MMU_4**2)
!!!   PMU_4 = Sqrt(EMU_3**2 - MMU_4**2)
      P_ZMU_4(1) = -PMU_3 * Sin(THEMU_3) * Cos(PHIMU_3)
      P_ZMU_4(2) = -PMU_3 * Sin(THEMU_3) * Sin(PHIMU_3)
      P_ZMU_4(3) = -PMU_3 * Cos(THEMU_3)
      P_ZMU_4(4) =  EMU_3
      P_ZMU_4(5) =  MMU_4

CsB___Check four momentum conservation in the Z1 and Z2 rest frames
!      Call CheckMomenta(P_L1,P_ZMU_1,P_ZMU_2,0,1,0)

CsB___Check mass-shell condition
!      Call CheckMomenta(P_ZMU_1,P_ZMU_3,P_ZMU_3,4,0,0)

C (4) BOOST THEM INTO THE C.M. FRAME OF ZZ PAIR (I.E. HIGGS REST FRAME)
      CALL BOOST_5 (P_ZMU_1, P_L1CS, P_HMU_1)
      CALL BOOST_5 (P_ZMU_2, P_L1CS, P_HMU_2)
      CALL BOOST_5 (P_ZMU_3, P_L2CS, P_HMU_3)
      CALL BOOST_5 (P_ZMU_4, P_L2CS, P_HMU_4)

CsB___Check mass-shell condition
!      Call CheckMomenta(P_HMU_1,P_HMU_3,P_HMU_3,4,0,0)

CsB___Check energy-momentum conservation in the Higgs rest (CS) frame
!      Call CheckMomenta(P_L1CS,P_HMU_1,P_HMU_2,0,1,0)
!      Call CheckMomenta(P_L2CS,P_HMU_3,P_HMU_4,0,1,0)

C (5) BOOST THEM INTO LAB FRAME
CsB___BOOST_5 is not the right CS -> Lab boost!
c      CALL BOOST_5 (P_HMU_1, P_Q, P_MU_1)
c      CALL BOOST_5 (P_HMU_3, P_Q, P_MU_3)
c      CALL BOOST_5 (P_HMU_2, P_Q, P_MU_2)
c      CALL BOOST_5 (P_HMU_4, P_Q, P_MU_4)
c      CALL BOOST_5 (P_L1CS, P_Q, P_L1)
CsB___The right CS -> Lab boost from Balazs,Yuan,PRD'97
      Call CS2Lab (P_HMU_1,P_MU_1,P_Q)
      Call CS2Lab (P_HMU_2,P_MU_2,P_Q)
      Call CS2Lab (P_HMU_3,P_MU_3,P_Q)
      Call CS2Lab (P_HMU_4,P_MU_4,P_Q)

CsB   Check lenght of P_MU_1 + P_MU_2
!      Call CheckMomenta (P_L1,P_MU_1,P_MU_2,0,2,0)

C (6) Rotate in Lab frame around z axis to get phi angles right
      tmp=P_L1(1)*(P_MU_1(1)+P_MU_2(1))+P_L1(2)*(P_MU_1(2)+P_MU_2(2))
      CosAl = tmp/(P_L1(1)**2+P_L1(2)**2)
      tmp=P_L1(1)*(P_MU_1(2)+P_MU_2(2))-P_L1(2)*(P_MU_1(1)+P_MU_2(1))
      SinAl = tmp/(P_L1(1)**2+P_L1(2)**2)
      tmp       =  CosAl*P_MU_1(1) + SinAl*P_MU_1(2)
      P_MU_1(2) = -SinAl*P_MU_1(1) + CosAl*P_MU_1(2)
      P_MU_1(1) =  tmp
      tmp       =  CosAl*P_MU_2(1) + SinAl*P_MU_2(2)
      P_MU_2(2) = -SinAl*P_MU_2(1) + CosAl*P_MU_2(2)
      P_MU_2(1) =  tmp

      tmp=P_L2(1)*(P_MU_3(1)+P_MU_4(1))+P_L2(2)*(P_MU_3(2)+P_MU_4(2))
      CosAl = tmp/(P_L2(1)**2+P_L2(2)**2)
      tmp=P_L2(1)*(P_MU_3(2)+P_MU_4(2))-P_L2(2)*(P_MU_3(1)+P_MU_4(1))
      SinAl = tmp/(P_L2(1)**2+P_L2(2)**2)
!      Print *,' SinAl**2 + CosAl**2 = ', SinAl**2 + CosAl**2
      tmp       =  CosAl*P_MU_3(1) + SinAl*P_MU_3(2)
      P_MU_3(2) = -SinAl*P_MU_3(1) + CosAl*P_MU_3(2)
      P_MU_3(1) =  tmp
      tmp       =  CosAl*P_MU_4(1) + SinAl*P_MU_4(2)
      P_MU_4(2) = -SinAl*P_MU_4(1) + CosAl*P_MU_4(2)
      P_MU_4(1) =  tmp

CsB___Check mass-shell condition
!      Call CheckMomenta(P_MU_1,P_MU_3,P_HMU_3,4,0,0)

CsB___Check energy-momentum conservation in the lab frame
!      Call CheckMomenta(P_L1,P_MU_1,P_MU_2,0,1,0)
!      Call CheckMomenta(P_L2,P_MU_3,P_MU_4,0,1,0)

CsB   Check length of transverse components
!      Call CheckMomenta (P_L1,P_MU_1,P_MU_2,0,2,0)
     
!C P_MU_2 AND P_MU_4 CAN BE OBTAINED USING P_MU_2 = P_L1 - P_MU_1, ETC.
!      Do ii=1,4
!!        P_MU_2(ii) = PLAB(ii,3) - P_MU_1(ii)
!!        P_MU_4(ii) = PLAB(ii,4) - P_MU_3(ii)
!        P_MU_2(ii) = P_L1(ii) - P_MU_1(ii)
!        P_MU_4(ii) = P_L2(ii) - P_MU_3(ii)
!      End Do

CsB___E+ and MU+ are hard wired to be massless here.
!      MMU_2 = 0.d0
!      P_MU_2(4) = Sqrt(P_MU_2(1)**2 + P_MU_2(2)**2 + P_MU_2(3)**2 + 
!     >                 MMU_2**2)
!      MMU_4 = 0.d0
!      P_MU_4(4) = Sqrt(P_MU_4(1)**2 + P_MU_4(2)**2 + P_MU_4(3)**2 + 
!     >                 MMU_4**2)
!      Print *,' P_MU_2(4),P_MU_4(4) =',P_MU_2(4),P_MU_4(4)

CsB___Check |P_L1|,|P_L2|
!      Print *,' |P_L1|,|P_L2| = ', 
!     -          Sqrt(-p_L1(1)*p_L1(1) - p_L1(2)*p_L1(2) -
!     -          p_L1(3)*p_L1(3) + p_L1(4)*p_L1(4)),
!     -          Sqrt(-p_L2(1)*p_L2(1) - p_L2(2)*p_L2(2) -
!     -          p_L2(3)*p_L2(3) + p_L2(4)*p_L2(4))

CsB___Double check |P_L1|,|P_L2|
!        Print *,'|PLAB(_,3)|,|PLAB(_,4)| = ',
!     >          Sqrt(ScaPro(3,3)),Sqrt(ScaPro(4,4))

C USING PAPAGENO'S NOTATION, WE SET SEC_PLAB(1) = P_MU_1, ETC.
      Do ii=1,4
        SEC_PLAB(ii,1) = P_MU_1(ii)
        SEC_PLAB(ii,2) = P_MU_2(ii)
        SEC_PLAB(ii,3) = P_MU_3(ii)
        SEC_PLAB(ii,4) = P_MU_4(ii)
      End Do

CsB___Test SEC_PLAB
!      Print *,' SEC_PLAB(i)^2 = ',
!     >        ScaProS(1,1),ScaProS(2,2),ScaProS(3,3),ScaProS(4,4)
!      Print *,' |SEC_PLAB(1+2)| = ',
!     >        Sqrt(ScaProS(1,1)+ScaProS(2,2)+2.*ScaProS(1,2))
!      Print *,' |SEC_PLAB(3+4)| = ',
!     >        Sqrt(ScaProS(3,3)+ScaProS(4,4)+2.*ScaProS(3,4))

      RETURN
      END     ! OF SUBROUTINE HZZDECAY      

c	***********************************************************
      SUBROUTINE ZFFB(RAND3,RAND4,XMZ0,IZ0POL,THEMU,EMU)
C THE ANGULAR DISTRIBUTION OF FERMION IN Z0 REST FRAME
c INPUT:
c      rand3,rand4=random numbers
c      xmz0=mass of decaying vector boson
c      izpol=select polarization distribution desired
c            izpol=1      longitudinal
c            otherwise, 50% right circularly polarized
c                       50% left circularly polarized
c OUTPUT:
c      themu=polar angle of lepton in decay of boson
c      emu=energy of daughter lepton in CMS of parent boson(=xmz0/2)
c

      REAL*8 AA(4),BB(4),CC(4),RAND3,RAND4,XMZ0,THEMU,EMU
      INTEGER*4 IZ0POL

!      Print *, ' Test 2 -> ',RAND3,RAND4,XMZ0,IZ0POL,THEMU,EMU

      IF( IZ0POL .EQ. 1)THEN
C FOR HELICITY 0 ANGULAR DISTRIBUTION
       CC(1)=  1.D0
       CC(2)=  0.D0
       CC(3)= -3.D0
       CC(4)=  4.D0*RAND4-2.D0
       CALL TRIPROOT(CC,CS0)
       THEMU= ACOS( CS0 )
       GO TO 31
      END IF
      SW2= 0.23D0
      V= -4.D0*SW2 +1.D0
      A= 1.D0+ V**2
      XNOR= 8.D0*A/3.D0
      IF( RAND4 .LE. 0.5)THEN
C FOR HELICITY +1 ANGULAR DISTRIBUTION
       B= 4.D0*V
       AA(1)= A/3.D0
       AA(2)= B/2.D0
       AA(3)= A
       AA(4)= 4.D0*A/3.D0 - B/2.D0 - RAND3*XNOR
       CALL TRIPROOT(AA,CSP1)
       THEMU= ACOS(CSP1)
       GO TO 31
      ELSE
C FOR HELICITY -1 ANGULAR DISTRIBUTION
       B= -4.D0*V
       BB(1)= A/3.D0
       BB(2)= B/2.D0
       BB(3)= A
       BB(4)= 4.D0*A/3.D0 - B/2.D0 - RAND3*XNOR
       CALL TRIPROOT(BB,CSM1)
       THEMU= ACOS(CSM1)
      END IF
31      EMU=XMZ0/2.
      RETURN
      END

      SUBROUTINE TRIPROOT(A,CS)
      REAL*8 A(4),X(3)
C SOLVING X FOR  A(1)*X**3 + A(2)*X**2 + A(3)*X + A(4) = 0.
      CALL DDVETA(A,X,MTYPE)
      IF( MTYPE .EQ. -1 )THEN
       IF(ABS(X(1)) .LE. 1.D0 )THEN
        CS=X(1)
       ELSE IF(ABS(X(2)) .LE. 1.D0)THEN
        CS=X(2)
       ELSE IF(ABS(X(3)) .LE. 1.D0)THEN
        CS=X(3)
       ELSE
        PRINT*,' ERROR-A IN TRIPROOT '
       END IF
      ELSE IF( MTYPE .EQ. 1 )THEN
       IF(ABS(X(1)) .LE. 1.D0 )THEN
        CS=X(1)
       ELSE
        PRINT*,' ERROR-B IN TRIPROOT ',X(1),X(2),X(3)
       END IF
      ELSE
       PRINT*,' ERROR-C IN TRIPROOT '
      END IF
      RETURN
      END
 
      SUBROUTINE DDVETA (D, RT, MTYPE)
      DIMENSION C(4), RT(3), X(2) ,D(4)
      DOUBLE PRECISION  C      , RT     , X      , D      , A
      DOUBLE PRECISION  B      , Q      , DABS   , R      , DSQRT
      DOUBLE PRECISION  BB     , PHI1   , DATAN2 , POD    , EXPO
      DOUBLE PRECISION  DLOG   , DEXP   , DCOS   , DSIN   , CLOD
      DOUBLE PRECISION  PLOD
      DO 10 L=2,4
   10 C(L)=D(L)/D(1)
      A=(3.0D0*C(3)-C(2)**2)/3.0D0
      B=(2.0D0*(C(2))**3-9.0D0*C(2)*C(3)+27.0D0*C(4))/27.0D0
      Q=B**2/4.0D0+A**3/27.0D0
      IF(B) 600,700,600
  700 IF(DABS(A)-0.01D0)12,12,300
  600 IF(DABS(Q)-10.0D0**(-15)*B**2/4.0D0) 12,300,300
  300 IF(Q) 11,12,13
   12 MTYPE=0
      GO TO 14
   13 MTYPE=1
      GO TO 15
   11 MTYPE=-1
      Q =DABS(Q)
      R =DSQRT(Q)
      BB=-B/2.0D0
      PHI1=DATAN2(R,BB)/3.0D0
      POD =DSQRT(B**2/4.0D0+Q)
      IF(POD) 73,70,73
   70 X(1)=0.0D0
      X(2)=0.0D0
      GO TO 17
   73 EXPO=(DLOG(POD))/3.0D0
      POD=DEXP(EXPO)
      X(1)= POD*DCOS(PHI1)
      X(2)= POD*DSIN(PHI1)
      GO TO 17
   14 POD=-B/2.0D0
      IF(POD) 83,80,83
   80 X(1)=0.0D0
      GO TO 84
   83 EXPO=(DLOG(DABS(POD)))/3.0D0
      X(1)=(POD /DABS(POD))*DEXP(EXPO)
   84 RT(1)=2.0D0*X(1)
      RT(2)=-X(1)
      RT(3)=RT(2)
      GO TO 199
   15 Q=DSQRT(Q)
      CLOD=-B/2.0D0+Q
      PLOD=CLOD-2.0D0*Q
      IF (CLOD) 96, 95, 96
   95 X(1)=0.0D0
      GO TO 97
   96 EXPO=(DLOG(DABS(CLOD)))/3.0D0
      X(1)=(CLOD/DABS(CLOD))*DEXP(EXPO)
   97 IF (PLOD) 91, 90, 91
   90 X(2)=0.0D0
      GO TO 16
   91 EXPO=(DLOG(DABS(PLOD)))/3.0D0
      X(2)=(PLOD/DABS(PLOD))*DEXP(EXPO)
   16 RT(1)=X(1)+X(2)-C(2)/3.0D0
      RT(2)=-.50D0*(X(1)+X(2))-C(2)/3.0D0
      RT(3)= .8660254037844385D0*(X(1)-X(2))
      RETURN
   17 RT(1)=2.0D0*X(1)
      RT(2)=-X(1)+X(2)*1.732050807568877D0
      RT(3)=-X(1)-X(2)*1.732050807568877D0
  199 DO 200 L=1,3
  200 RT(L)=RT(L)-C(2)/3.0D0
   18 RETURN
       END
c	***********************************************************
CsB____________________________________________________________________________
C     This routine is to transform momenta from the CS frame to the Lab frame
C     It follows my thesis and/or from the appendix of Balazs,Yuan,PRD'97
C     It implements the formula Lorentz transformation matrix.
C     Note that phi_V = 0 is assumed here.
C     q    is the momentum of which rest frame defines the CS frame
C     pCS  is the monetum to be transfarred from CS to the Lab frame
C     pLab is the result
C
      Subroutine CS2Lab (pCS,pLab,q)
      Implicit Real*8 (a-h,o-z)
      Real*8 pCS(5),pLab(5),q(5),MT

!      Print* ,'pCS = ',pCS

      QT = Sqrt(q(1)**2 + q(2)**2)
      MT = Sqrt(q(5)**2 + QT**2)
      PF = 1.d0/q(5)/MT

      pLab(1) = PF*(QT*MT*pCS(4) + MT**2*pCS(1))
      pLab(2) = pCS(2)
      pLab(3) = PF*(q(3)*MT*pCS(4) + q(3)*QT*pCS(1) + q(4)*q(5)*pCS(3))
      pLab(4) = PF*(q(4)*MT*pCS(4) + q(4)*QT*pCS(1) + q(3)*q(5)*pCS(3))
      pLab(5) = pCS(5)
      
!      Print* ,'pLab = ',pLab
      
      Return
      End
CsB____________________________________________________________________________
C     This routine is to transform momenta from the CS frame to the Lab frame
C     It follows my thesis and/or from the appendix of Balazs,Yuan,PRD'97
C     It implements the formula with the X, Y, and Z vectors. 
C     Note that phi_V = 0 is assumed here.
C     q    is the momentum of which rest frame defines the CS frame
C     pCS  is the monetum to be transfarred from CS to the Lab frame
C     pLab is the result
C
C     ! Only for massless particles !
C     !!! Not tested !!!
C
      Subroutine CS2LabXYZ (pCS,pLab,q,X,Y,Z)
      Implicit Real*8 (a-h,o-z)
      Real*8 pCS(5),pLab(5),q(5),X(5),Y(5),Z(5),
     ,       CosTh,SinTh,CosPh,SinPh,tmp

!      Print *,' P_MU_1(4), P_MU_2(4) = ', P_MU_1(4), P_MU_2(4)

C     Angles of MU_1 in CS frame
      CosTh = pCS(3)/pCS(4)
      SinTh = Sqrt(1.d0 - CosTh**2) ! I'm not sure that this is correct !
      tmp   = pCS(2)/pCS(1)
      CosPh = Cos(ATan(tmp)) ! I'm not sure that this is correct !
      SinPh = Sin(ATan(tmp)) ! I'm not sure that this is correct !

C     Boost
      Do ii = 1,4
        pLab(ii) = pCS(4)*(q(ii)/q(5) +
     *    SinTh*CosPh*X(ii)+SinTh*SinPh*Y(ii)+CosTh*Z(ii))
!       Print *,' P_MU_1(4), P_MU_2(4) = ', P_MU_1(4), P_MU_2(4)
      EndDo
      pLab(5) = pCS(5)
      
      Return
      End
CsB____________________________________________________________________________
      Subroutine CheckMomenta(p1,p2,p3,n1,n2,n3)
      Implicit Real*8 (a-h,o-z)
      Real*8 p1(5),p2(5),p3(5)

CsB___Check mass-shell condition
      If (n1.Eq.1 .Or. n1.Eq.4 .Or. n1.Eq.5) Print 
     **,' m1 = ',-p1(1)**2-p1(2)**2-p1(3)**2+p1(4)**2,p1(5)
      If (n1.Eq.2 .Or. n1.Eq.4 .Or. n1.Eq.5) Print
     **,' m2 = ',-p2(1)**2-p2(2)**2-p2(3)**2+p2(4)**2,p2(5)
      If (n1.Eq.3 .Or. n1.Eq.5) Print
     **,' m3 = ',-p3(1)**2-p3(2)**2-p3(3)**2+p3(4)**2,p3(5)

CsB___Check energy-momentum conservation
      If (n2.Eq.1) then
        Do ii=1,5
          Print *,' p1(',ii,'),p2(',ii,')+p3(',ii,') = ', 
     ,              p1(ii),p2(ii)+p3(ii)
        EndDo
      EndIf

CsB   Check lenght of transverse component
      If (n2.Eq.2) then
        Print *,' (p1(1))**2+(p1(2))**2 = ',
     ,            (p1(1))**2+(p1(2))**2
        Print *,' (p2(1)+p3(1))**2+(p2(2)+p3(2))**2 ',
     ,            (p2(1)+p3(1))**2+(p2(2)+p3(2))**2
      EndIf

      Return
      End
CsB____________________________________________________________________________
	Function RESON(XMW,XGW,XNUM)
	IMPLICIT REAL*8 (A-H,O-Z)
	EXTERNAL RN
C GENERATE BREIT-WIGNER RESONATING MASS FOR M12
C        XMW=80.0D0
C        XGW=2.1D0
C +- XNUM WIDTH
C       XNUM=1.0
        XMW_MIN=-XNUM
        XMW_MAX=XNUM
!        Print *,' XMW_MIN,XMW_MAX =',XMW_MIN,XMW_MAX
        ONE=1.D0
        Y=RN(1)
C        XMW_VAR=DATAN(XMW_MAX)*Y+(ONE-Y)*DATAN(XMW_MIN)
C        XMW_NORM=DATAN(XMW_MAX)-DATAN(XMW_MIN)
C        XM12SQ=XMW*XMW+XMW*XGW*DTAN(XMW_VAR)
C        Reson=DSQRT(XM12SQ)
        XMW_VAR=ATAN(XMW_MAX)*Y+(ONE-Y)*ATAN(XMW_MIN)
        XMW_NORM=ATAN(XMW_MAX)-ATAN(XMW_MIN)
        XM12SQ=XMW*XMW+XMW*XGW*TAN(XMW_VAR)
        Reson=SQRT(XM12SQ)
!        Print *,' Y,Reson,XMW_VAR,XMW_NORM,XM12SQ = '
!     ,           ,Y,Reson,XMW_VAR,XMW_NORM,XM12SQ
        YJACO=1.
	RETURN
	END
CsB____________________________________________________________________________
CsB   This function generates the phi angular correlation between two Z's
C     for H -> Z Z -> 4l. Formulae are based on Duncan,Kane,Repko,WW Physics.
C
      Function PhiCorrZZ (mZ,mH)
      Implicit Real*8 (a-h,o-z)
      Real*8 mZ,mH
      External RN
      DATA PI/3.14159 26535 89793 24D0/
      
 100  Continue

      phi = 2.d0*Pi*RN(1)
      D   = 2.d0*mZ**4/((mH**2-2.d0*mZ**2)**2 + 8.d0*mZ**4)
      F   = 1.d0 + D*Cos(2.d0*phi)
      tmp = RN(1)*(1.d0+D)
      If (tmp.LT.F) then
        PhiCorrZZ = phi
      Else
        Goto 100
      End If
      
!      Print *,' PhiCorrZZ = ', PhiCorrZZ
      
      Return
      End
CsB_________________________________________________________________________
CsB End of generating the Z decays for the process: 
C   H -> Z Z -> 4 leptons.
CsB_________________________________________________________________________
C ==========================================================================
CCPY FEB 2001: ADD IN A SPECIAL PHASE SPACE FOR AA PRODUCTION

      SUBROUTINE AAJPH(Y,WT)
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

      REAL*8 YM12MIN,YM12MAX,QTMin,QTMax, TMmin, Tmmax, ETmin
      COMMON/MYYM12/YM12MIN,YM12MAX,QTMin,QTMax, TMmin, Tmmax, ETmin

      REAL*8 Y,WT
      DIMENSION Y(50)
      INTEGER I
      REAL*8 R,PK1,PK2,YQ1,YQ3
      DIMENSION R(8),PK1(5),PK2(5),YQ1(5),YQ3(5)
      REAL*8 ZTEMP,JACOBIAN
      REAL*8 TAU,YBMAX,WT_WHAT,WT_YB,YBOOST
      INTEGER ICHOICE
      REAL*8 ONE,YMT,YGT,YMT_MIN,YMT_MAX,YMT_VAR,YMT_NORM,YMTOPSQ
      REAL*8 RM21,RMG1,GAMOM1,Q1SQM,Y1MAX,X1MAX,X1MIN,XV1
      REAL*8 YM3,YM4

      IF(YM12MIN.LT.Y_EMIN) YM12MIN=Y_EMIN
CCPY THE MEANING OF YM12MAX IN THIS CODE IS NOT THE USUAL ONE.
      X1X2MIN=Y_EMIN**2/W**2
      ZTEMP=1.D0/((1./X1X2MIN-1.)*Y(6)+1.)
      X1=ZTEMP**Y(7)
      X2=ZTEMP/X1
      WHAT=W*SQRT(X1*X2)
      JACOBIAN=-ZTEMP**2*LOG(ZTEMP)*(1./X1X2MIN-1.)

      PLAB(1,1)=0.D0
      PLAB(2,1)=0.D0
      PLAB(3,1)=.5*W*X1
      PLAB(4,1)=.5*W*X1
      PLAB(1,2)=0.D0
      PLAB(2,2)=0.D0
      PLAB(3,2)=-.5*W*X2
      PLAB(4,2)=.5*W*X2

      DO I=1,5
        R(I)=Y(I)
      ENDDO

      DO I=1,4
        PK1(I)=PLAB(I,1)
        PK2(I)=PLAB(I,2)
      ENDDO

      IF( NPART.eq. 21 .or. NPART.eq.22) THEN
        YM3=0.D0
        YM4=0.D0
      ELSEIF(NPART.EQ.23) THEN
        YM3=XMZ
        YM4=XMZ
      ELSE
        WRITE(*,*) '  AAJPH IS NOT FOR THIS PROCESS '
        CALL ExitGracefully
      ENDIF
            
      CALL AAJXPH(R,PK1,PK2,WT,YM3,YM4,YM12MIN,YM12MAX)

      WT=WT*JACOBIAN

      RETURN
      END
C
      SUBROUTINE AAJXPH(R,PK1,PK2,WT,YM3,YM4,YM34MIN,YM34MAX)
C
C THIS IS TO GENERATE THE PHASE SPACE FOR
C     q(1) + qbar(2) --> V(3) V(4) g(5)
C   QUARK(1) + ANTI-QUARK(2) --> GLUON(5) AA
C                                          \--> A(3) A(4)
C
C R(5)   :  RANDOM NUMBER ARRAY (BETWEEN 0 AND 1)
C WT     :  WEIGHT OF THE PRODUCED PHASE SPACE POINT (GE.0)
C PK1,PK2:  THE MOMENTA OF THE INCOMING PARTICLES
C
C THE FOURTH COMPONENT OF THE MOMENTA IS THE TIME COMPONENT
C THE FIFTH COMPONENT CONTAINS THE MASS
C
C The incoming momenta PK1 and PK2 are assumed along the z-axis

      IMPLICIT NONE
      REAL*8 WT,YM3,YM4,YM34MIN,YM34MAX
      INTEGER I
      REAL*8 R,RT,PK1,PK2,Q1,Q3,P3,P4,P5,P6,P8
      DIMENSION R(8),RT(2),PK1(5),PK2(5),Q1(5),Q3(5),
     >  P3(5),P4(5),P5(5),P6(5),P8(5)
      REAL*8 CM,YP1,YP2,YP3,YP4,YP5,YP6
      DIMENSION CM(5),YP1(5),YP2(5),YP3(5),YP4(5),YP5(5),YP6(5)
      REAL*8 PI,SH,RSH,WT0,WT1,WT2,WT3,WT4,RRSH
      DIMENSION RRSH(2)
      REAL*8 CP1,CP2,CP3,CP4,CP5,CP6
      DIMENSION CP1(5),CP2(5),CP3(5),CP4(5),CP5(5),CP6(5)

      REAL*8 PLAB,PCEM,WHAT
      COMMON/MOM/PLAB(4,12),PCEM(4,12),WHAT

      REAL*8 YM34,WT_YM34
      REAL*8 YXDOT
      EXTERNAL YXDOT
      REAL*8 ZERO,ONE,TWO,THR
      DATA ZERO,ONE,TWO,THR/0.D0,1.D0,2.D0,3.D0/

      
      WT=-1.D0
      PI=4.0*ATAN(1.0D0)
      DO 1 I=1,4
        CM(I)=PK1(I)+PK2(I)
    1 CONTINUE
      SH = CM(4)**2 - CM(3)**2
      RSH = SQRT(SH)

      YM34=YM34MIN+(YM34MAX-YM34MIN)*R(1)
      WT_YM34=2.0D0*YM34*(YM34MAX-YM34MIN)

      IF(RSH.LT.YM34 .OR. YM34.LT.(YM3+YM4))THEN
!        PRINT*,' RSH < YM34  OR  YM34 < YM3+YM4 '
!        PRINT*,' RSH,YM34MIN,YM34MAX,YM34,YM3,YM4'
!        PRINT*, RSH,YM34MIN,YM34MAX,YM34,YM3,YM4 
        RETURN
      ENDIF

C GENERATE AA PRODUCTION:  q Q --> G (AA)
      DO I=1,2
        RRSH(I)=R(I+1)
      ENDDO

      CP1(5)=0.D0
      CP2(5)=YM34

      CALL X2PHASE(RRSH,RSH,CP1(5),CP2(5),YP1,YP2,WT1)

      DO I=1,4
        CP1(I)=YP1(I)
        CP2(I)=YP2(I)
      ENDDO

C GENERATE THE DECAY OF (AA) --> A A
      DO I=1,2
        RT(I)=R(I+3)
      ENDDO

      CALL X2PHASE(RT,YM34,YM3,YM4,YP3,YP4,WT2)

      CALL BOOST_5(YP3,CP2,CP3)

      DO I=1,4
        CP4(I)=CP2(I)-CP3(I)
      ENDDO

C FILL IN VECTORS
C   QUARK(1) + ANTI-QUARK(2) --> GLUON(5) AA
C                                          \--> A(3) A(4)
      DO I=1,2
        PCEM(I,1)=0.D0
        PCEM(I,2)=0.D0
      ENDDO

      PCEM(3,1)=RSH/2.D0
      PCEM(4,1)=PCEM(3,1)
      PCEM(3,2)=-PCEM(3,1)
      PCEM(4,2)=PCEM(3,1)
      DO I=1,4
        PCEM(I,3)=CP3(I)
        PCEM(I,4)=CP4(I)
        PCEM(I,5)=CP1(I)
      ENDDO

C BOOST QB SYSTEM
      CALL BOOST_5(CP1,CM,P5)
      CALL BOOST_5(CP3,CM,P3)
      CALL BOOST_5(CP4,CM,P4)
      DO I=1,4
        PLAB(I,3)=P3(I)
        PLAB(I,4)=P4(I)
        PLAB(I,5)=P5(I)
      ENDDO

      WT=WT_YM34*WT1*WT2/(2.0*PI)

      RETURN
      END

      SUBROUTINE X2PHASE(R,YMW,YM1,YM2,P1,P2,WT)
C
C THIS IS TO GENERATE THE PHASE SPACE FOR
C     M1 --> P1 + P2
C THIS SUBROUTINE ONLY WORKS IN THE REST FRAME OF M1
C
C DEFINE PHASE SPACE VOLUME = (2*PI)**4*D(2_BODY PHASE SPACE)
C IT DOES NOT CONTAIN THE FACTOR OF FLUX = 1/(2*YM1)
C
C R(2): RANDOM NUMBER ARRAY (BETWEEN 0 AND 1)
C WT:   WEIGHT OF THE PRODUCED PHASE SPACE POINT (GE.0)
C
      IMPLICIT NONE
      REAL*8 R,YMW,YM1,YM2,P1,P2,WT
      DIMENSION R(2),P1(5),P2(5)
      REAL*8 PI,YY1,THE_STA,PHI_STA,RSH,SH,TEMP1,
     >Q1CM,SPH1,CPH1,STH1,CTH1,DPHASE2,YJACO1,YJACO2
      INTEGER*4 I
      REAL*8 COSTH
      COMMON/X2PH/COSTH

      PI=4.0*ATAN(1.0D0)
C       CHANGE THE VARIABLES ( COS(THETA_STAR), PHI_STAR ) IN C.M.
C       FRAME TO ( YY1, R(2) ).
C       YY1:   -1 ----> 1
C       CHANGE VARIABLE FROM YY1 TO R(1), AND R(1):  0 ---> 1
      YY1=2.D0*R(1)-1.D0
      COSTH=YY1
      THE_STA=DACOS(YY1)
      PHI_STA=2.D0*PI*R(2)
      RSH=YMW
      SH=RSH**2
      P1(5)=YM1
      P2(5)=YM2
      TEMP1=(SH + YM1**2 - YM2**2)/2.D0/RSH
      Q1CM=DSQRT(TEMP1**2-YM1**2)
      SPH1=DSIN(PHI_STA)
      CPH1=DCOS(PHI_STA)
      CTH1=DCOS(THE_STA)
      STH1=DSIN(THE_STA)
      P1(1)=Q1CM*CPH1*STH1
      P1(2)=Q1CM*SPH1*STH1
      P1(3)=Q1CM*CTH1
      P1(4)=(SH + YM1**2 - YM2**2)/2.D0/RSH
      DO 10 I=1,3
        P2(I)=-P1(I)
   10 CONTINUE
      P2(4)=YMW-P1(4)
C -------------------------------
      DPHASE2=Q1CM/(4.D0*RSH*(2.D0*PI)**6)
C JACOBINS
C FOR THE_STA
      YJACO1=2.D0
C FOR PHI_STA
      YJACO2=2.D0*PI
C ADD IN (2*PI)**4 FROM DELTA-FUNCTION
      WT=YJACO1*YJACO2*(2.D0*PI)**4*DPHASE2
C
      RETURN
      END
CsB_________________________________________________________________________
C
C   END OF SPECIAL PHASE SPACE FOR AA PRODUCTION
CsB_________________________________________________________________________
CsB_________________________________________________________________________
C
C   Begin: h^(1) function for gg->AA
CsB_________________________________________________________________________
cpn   -------------------------------0--------------------------------------
      function h1ggAA(vin,Nf)
cpn 2004 Returns the O(\alpha_S/pi) part h^1 of the hard part function
cpn      H(Q) in the gg -> AA processes. 
cpn   Author: Pavel Nadolsky (nadolsky@hep.anl.gov).
cpn   
cpn   h^(1) for the gg->AA process is defined as h^(1) = V/4,
cpn   where V is the O(\alpha_S/\pi) term in the NLO coefficient
cpn   function C^(1)_{g/g} given in Eq. (9) in the paper P.M. Nadolsky,
cpn   C.R.Schmidt, Phys. Lett. B558, 63 (2008);
cpn
cpn June 2006: The O(alpha_s) hard vertex coefficient is evaluated 
cpn   assuming the scale mu=Q, in accordance with my scale.ps note on
cpn   scale choices in various resummation schemes
cpn
cpn   Inputs: 
cpn      v = -t/s -- a function of parton-level Mandelstam
cpn      parameters t and s in the AA c.m. frame (0<=v<=1)
cpn      N_f -- number of active quark flavors
cpn
      implicit none
      include 'pncommon.f'
      double precision N1ggAA, D1ggAA !numerator and denominator functions
                                      !in the hard part of h1 for gg->AA
      double precision vin, h1ggAA, tmp, Ch
      integer Nf
      external N1ggAA, D1ggAA

cpn   Given v, set up auxiliary arrays with commonly used functions of v
      call Setuph1(vin)

cpn   This is V/4
      tmp = 0.25d0*(Nc*Pi2 + N1ggAA()/D1ggAA())

      h1ggAA = tmp
      return 
      end !h1ggAA ->

cpn   ---------------------------------------------------------------------

      subroutine Setuph1(vin)
cpn   Stores common functions of v in an array vvv1 to speed the
cpn   calculation. Stores the same functions  of (1-v) in another
cpn   array vvv2. The arrays vvv1 and vvv2 can be swapped in subsequent 
cpn   helicity amplitudes when symmetrizing these amplitudes with
cpn   respect to the replacement t <-> u (i.e., with respect to v <-> 1-v).
      implicit none
      include 'pncommon.f'
      double complex fLi2         
      double precision fLi3, fLi4 
      double precision vin

      if (vin.lt.0d0.or.vin.gt.1d0) then
        print *,'STOP in Setuph1: vin out of correct bounds, 0<=vin<=1'
        print *,'vin=',vin
        call Quit
      endif !v

cpn   Relative error of integration in the polylogarithms
      rlierr = 1e-5

      ie=(0d0,1d0)                !imaginary unit
C_____Parameters
      Pi = 3.14159265358979d0
      Pi2 = Pi**2
      Nc = 3 !number of colors
      zeta3=1.20205690315959d0

      v = vin
      onemv = 1-v
      v1mv = v/onemv
      onemvv = onemv/v
                                !Functions of v
      vvv1(1)  = v
      vvv1(2)  = v**2
      vvv1(3)  = dlog(v)
      vvv1(4)  = dlog(v)**2
      vvv1(5)  = 1 + v
      vvv1(6)  = dlog(v1mv)
                                !In the current implementation, 
                                !the polylogarithm is realized as a complex
                                !function, while the polylogarithms of the 
                                !3d and 4th degree are realized as purely real 
                                !functions (valid for v<1) 
      vvv1(7)  = dreal(fLi2(dcmplx(v)))
      vvv1(8)  = fLi3(v)
      vvv1(9)  = fLi4(v)
      vvv1(10) = dreal(fLi2(dcmplx(-v1mv)))
      vvv1(11) = fLi3(-v1mv)
      vvv1(12) = fLi4(-v1mv)

cpn                              Functions of 1-v
      vvv2(1) = onemv
      vvv2(2) = onemv**2
      vvv2(3) = dlog(onemv)
      vvv2(4) = dlog(onemv)**2
      vvv2(5) = 1 + onemv
      vvv2(6) = dlog(onemvv)
      vvv2(7) = dreal(fLi2(dcmplx(onemv)))
      vvv2(8) = fLi3(onemv)
      vvv2(9) = fLi4(onemv)
      vvv2(10) = dreal(fLi2(dcmplx(-onemvv)))
      vvv2(11) = fLi3(-onemvv)
      vvv2(12) = fLi4(-onemvv)

      return 
      end !Setuph1 ->

cpn   ---------------------------------------------------------------------

      subroutine Returnv(vvv1,vvv2)
cpn   Restores values of v, 1-v, log(v), etc. from the arrays vvv1 and
cpn   vvv2. Useful for interchanging t and u (i.e., v and 1-v).

      implicit none
      integer Naux
      parameter (Naux=12)
      double precision vvv1(Naux), vvv2(Naux)

      double precision v, v2, dlogv, dlogv2, onepv, dlogv1mv, v1mv,
     >  fLi2v, fLi3v, fLi4v, fLi2v1mv, fLi3v1mv, fLi4v1mv, 
     >  onemv, onemv2, dlog1mv, dlog1mv2, twomv, dlog1mvv, onemvv, 
     >  fLi21mv, fLi31mv, fLi41mv,fLi21mvv, fLi31mvv, fLi41mvv 
      common/vs/ v, v2, dlogv, dlogv2, onepv, dlogv1mv, v1mv, 
     >  fLi2v, fLi3v, fLi4v, fLi2v1mv, fLi3v1mv, fLi4v1mv, 
     >  onemv, onemv2,dlog1mv, dlog1mv2, twomv, dlog1mvv, onemvv, 
     >  fLi21mv, fLi31mv, fLi41mv,fLi21mvv, fLi31mvv, fLi41mvv 

      v         =  vvv1(1) 
      v2        =  vvv1(2) 
      dlogv     =  vvv1(3) 
      dlogv2    =  vvv1(4) 
      onepv     =  vvv1(5) 
      dlogv1mv =   vvv1(6) 
      fLi2v     =  vvv1(7) 
      fLi3v     =  vvv1(8) 
      fLi4v     =  vvv1(9) 
      fLi2v1mv =  vvv1(10)
      fLi3v1mv =  vvv1(11)
      fLi4v1mv =  vvv1(12)
                                               
      onemv     =  vvv2(1)
      onemv2    =  vvv2(2)
      dlog1mv   =  vvv2(3)
      dlog1mv2  =  vvv2(4)
      twomv     =  vvv2(5)
      dlog1mvv  =  vvv2(6)
      fLi21mv   =  vvv2(7)
      fLi31mv   =  vvv2(8)
      fLi41mv   =  vvv2(9)
      fLi21mvv  =  vvv2(10)
      fLi31mvv  =  vvv2(11)
      fLi41mvv  =  vvv2(12)

      return 
      end ! Returnv

cpn   ------------------------------------------------------------------

      function D1ggAA()
cpn   The denominator of the hard part function h1(v) in gg->AA
cpn   Inputs: v and its simple functions (passed through a common block)
      implicit none
      include 'pncommon.f'
      double precision D1ggAA

      if (v.lt.0d0.or.v.gt.1d0) then
        print *,'STOP in D1ggAA: v out of correct bounds, 0<=v<=1'
        print *,'v=',v
        call Quit
      endif !v

      call Returnv(vvv1,vvv2)

cpn   This is the reprocessed result from my Mathematica notebook
      D1ggAA = (32 + 
     >  (dlog1mv2*(-(dlog1mv*(1 + onemv2)) - 2*twomv*v)**2)/v**4 + 
     >  (dlogv2*(-2*onemv*onepv - dlogv*(1 + v2))**2)/onemv**4 + 
     >  (2*dlogv1mv*(1 - 2*v) + (dlogv1mv**2 + Pi2)*(onemv2 + v2))**2 + 
     >  4*Pi2*((dlog1mv*(1 + onemv2) + twomv*v)**2/v**4 + 
     >  (onemv*onepv + dlogv*(1 + v2))**2/onemv**4) + 
     >  (4*(dlog1mv2*(onemv**4 + onemv2) + 2*dlog1mv*onemv2*twomv*v + 
     >  v2*(2*dlogv*onemv*onepv + dlogv2*(1 + v2) + 
     >  onemv2*(2*dlogv1mv*(1 - 2*v) + dlogv1mv**2*(onemv2 + v2) + 
     >  Pi2*(onemv2 + v2)))))/(onemv2*v2))/2d0

      return 
      end !D1ggAA ->

cpn   ---------------------------------------------------------------------

      function N1ggAA()
cpn   The numerator of the hard part function h1(v) in gg->AA
cpn   Inputs: v and its simple functions (passed through a common block)
      implicit none
      include 'pncommon.f'
      double precision N1ggAA, tmp
      integer i1, i2, i3, i4
      double complex fMLO, FL, FSL

      if (v.lt.0d0.or.v.gt.1d0) then
        print *,'STOP in N1ggAA: v out of correct bounds, 0<=v<=1'
        print *,'v=',v
        call Quit
      endif !v

cpn                                          Sum over the helicities
cpn                  The convention for helicity signs: 1="-", 2="+"
      tmp=0d0
      do i1=1,2
        do i2=1,2
          do i3=1,2
            do i4=1,2
cpn for debugging
c              tmpr = dreal(fSL(i1,i2,i3,i4))
c              tmpi = dimag(fSL(i1,i2,i3,i4))

cpn                   This is the numerator from Eq. (9) in P. Nadolsky,
cpn                   C. R. Schmidt, Phys. Lett. B558, 63 (2003) 
              tmp = tmp + 
     >          2d0*Nc*dreal(conjg(fMLO(i1,i2,i3,i4))*FL(i1,i2,i3,i4)) -
     >          2d0*dreal(conjg(fMLO(i1,i2,i3,i4))*FSL(i1,i2,i3,i4)/Nc)
            enddo 
          enddo 
        enddo 
      enddo 
      
      N1ggAA = tmp
      return 
      end !N1ggAA ->

cpn   ---------------------------------------------------------------------

      function fMLO(i1,i2,i3,i4)
cpn   Returns the leading-order amplitude from Eq. (3.15) in Bern, 
cpn   De Freitas, and Dixon
cpn   Inputs: i1, i2, i3, i4 --  signs of helicities of g's and A's
cpn   The convention for helicity signs: 1="-", 2="+"
      implicit none
      include 'pncommon.f'
      double complex fMLO
      integer i1, i2, i3, i4, ii
      
      ii = i1*1000 + i2*100 + i3*10 + i4
      
      call Returnv(vvv1,vvv2)

      if (ii.eq.1111.or.ii.eq.2222.or.
     >  ii.eq.1222.or.ii.eq.2122.or.ii.eq.2212.or.ii.eq.2221.or.
     >  ii.eq.2111.or.ii.eq.1211.or.ii.eq.1121.or.ii.eq.1112) then
        fMLO = (1d0,0d0)
      else if (ii.eq.1122.or.ii.eq.2211) then
        fMLO=-1d0 - dlogv1mv*(1 - 2*v) - 
     >    ((dlogv1mv**2 + Pi2)*(onemv2 + v2))/2d0
      else if (ii.eq.1212.or.ii.eq.2121) then
        fMLO=-0.5d0*ie*Pi*(2*onemv*onepv + 2*dlogv*(1 + v2))/onemv2 - 
     >    (2*onemv + 2*dlogv*onemv*onepv - 2*onemv*v 
     >    + dlogv2*(1 + v2))/(2d0*onemv2)
      else if (ii.eq.2112.or.ii.eq.1221) then
        fMLO=-0.5d0*ie*Pi*(2*dlog1mv*(1 + onemv2) + 2*twomv*v)/v2 - 
     >    (dlog1mv2*(1 + onemv2) + 2*v - 2*onemv*v + 2*dlog1mv*twomv*v)
     >    /(2d0*v2)
      endif
           
      return 
      end ! fMLO ->

      function fMLO2(v,Q)
cpn2005 Calculate the Born gg->gamma gamma matrix element squared, 
cpn2005 summed over the helicities of the external particles
      implicit none
      double precision fMLO2, v, tmp, Q, sq22
      double complex fMLO
      integer i1, i2, i3, i4
      external fMLO

      call Setuph1(v)
      
      tmp = 0.d0
      do i1=1, 2
        do i2=1,2
          do i3=1,2
            do i4=1,2
              tmp = tmp + 
     >          dreal(fMLO(i1,i2,i3,i4)*conjg(fMLO(i1,i2,i3,i4)))
            enddo 
          enddo 
        enddo 
      enddo 
      
      fMLO2 = 64d0*SQ22(Q)*tmp

      return 
      end !fMLO2

      function fMLO2prime(v,Q,phi)
cpn2005 Calculate the "primed" helicity sum of the Born gg->gamma gamma 
cpn2005 matrix element squared, as introduced in Z. Bern, L. Dixon, and
cpn2005 C. Schmidt, hep-ph/0206194
      implicit none
      double precision fMLO2prime, v, tmp, Q, sq22, phi
      double complex fMLO
      integer i1, i2, i3, i4
      external fMLO

      double precision  onemv, v1mv, onemvv,v2,dlogv,dlogv2,dlogv1mv,
     >  dlog1mv, dlog1mv2, onemv2, pi2
      data pi2 /9.869604401089358618/
 

cpn2005 dumb computation
c$$$      call Setuph1(v)
c$$$
c$$$      tmp = 0.d0
c$$$      do i1=1, 2
c$$$        do i2=1,2
c$$$          do i3=1,2
c$$$            do i4=1,2
c$$$              tmp = tmp + 
c$$$     >          dreal(fMLO(3-i1,i2,i3,i4)*conjg(fMLO(i1,i2,i3,i4)))
c$$$            enddo 
c$$$          enddo 
c$$$        enddo 
c$$$      enddo 


cpn2005 Faster computation
cpn This is 4*dreal(fMLO(1,1,2,2)+fMLO(1,2,1,2) + fMLO(1,2,2,1)+1d0)

      onemv = 1-v
      v1mv = v/onemv
      onemvv = onemv/v

      v2  = v**2
      onemv2 = onemv**2
      dlogv  = dlog(v)
      dlogv2  = dlog(v)**2
      dlog1mv = dlog(onemv)
      dlog1mv2 = dlog1mv**2
      dlogv1mv  = dlog(v1mv)
 
      tmp= 4*(-2 - dlogv1mv*(1 - 2*v) + (dlogv*(-1 - v))/onemv + 
     >  (dlog1mv*(-2 + v))/v - 
     >    (dlog1mv2*(1 + onemv2))/(2.*v2) - 
     >  (dlogv2*(1 + v2))/(2.*onemv2) - 
     >  ((dlogv1mv**2 + Pi2)*(onemv2 + v2))/2.)

       
      fMLO2prime = 64d0*SQ22(Q)*(-dcos(2*phi))*tmp

      return 
      end !fMLO2prime


cpn   ---------------------------------------------------------------------

      function fL(i1,i2,i3,i4)
cpn   Returns the NLO amplitude F^{L} from Eqs. (4.7-4.10) in Bern, 
cpn   De Freitas, and Dixon
cpn   Inputs: i1, i2, i3, i4 --  signs of helicities of g's and A's
cpn   The sign helicity convention: 1="-", 2="+"
      implicit none
      include 'pncommon.f'
      double complex fL, fL1222, fL2212, fL1122, fL1212
      integer i1, i2, i3, i4, ii
      
      ii = i1*1000 + i2*100 + i3*10 + i4
      
      if (ii.eq.1111.or.ii.eq.2222) then
        fL = 0.5d0
      else if (ii.eq.1222.or.ii.eq.2122.or.ii.eq.2111.or.ii.eq.1211) 
     >    then
        fL = fL1222(vvv1,vvv2) + fL1222(vvv2,vvv1)
      else if (ii.eq.2212.or.ii.eq.2221.or.ii.eq.1121.or.ii.eq.1112)
     >    then
        fL = fL2212(vvv1,vvv2) + fL2212(vvv2,vvv1)
      else if (ii.eq.1122.or.ii.eq.2211) then
        fL = fL1122(vvv1,vvv2) + fL1122(vvv2,vvv1)
      else if (ii.eq.1212.or.ii.eq.2121) then
        fL = fL1212(vvv1,vvv2)  !no symmetrization t <-> u
      else if (ii.eq.1221.or.ii.eq.2112) then
        fL = fL1212(vvv2,vvv1)  !M_{-++-}(s,t,u) = M_{-+-+}(s,u,t)
      endif
        
      return 
      end ! fL ->

cpn   ---------------------------------------------------------------------

      function fSL(i1,i2,i3,i4)
cpn   Returns the NLO amplitude F^{SL} from Eqs. (4.11-4.16) in Bern, 
cpn   De Freitas, and Dixon
cpn   Inputs: i1, i2, i3, i4 --  signs of helicities of g's and A's
cpn   The sign helicity convention: 1="-", 2="+"
      implicit none
      include 'pncommon.f'
      double complex fSL, fSL1222, fSL1122, fSL1212
      integer i1, i2, i3, i4, ii
      
      ii = i1*1000 + i2*100 + i3*10 + i4
      
      if (ii.eq.1111.or.ii.eq.2222) then
        fSL = -1.5d0
      else if (ii.eq.1222.or.ii.eq.2122.or.ii.eq.2111.or.ii.eq.1211.or.
     >    ii.eq.2212.or.ii.eq.2221.or.ii.eq.1121.or.ii.eq.1112) then
        fSL = fSL1222(vvv1,vvv2) + fSL1222(vvv2,vvv1)
      else if (ii.eq.1122.or.ii.eq.2211) then
        fSL = fSL1122(vvv1,vvv2) + fSL1122(vvv2,vvv1)
      else if (ii.eq.1212.or.ii.eq.2121) then
        fSL = fSL1212(vvv1,vvv2)  !no symmetrization t <-> u
      else if (ii.eq.1221.or.ii.eq.2112) then
        fSL = fSL1212(vvv2,vvv1)  !M_{-++-}(s,t,u) = M_{-+-+}(s,u,t)
      endif
        
      return 
      end ! fSL ->

cpn   ---------------------------------------------------------------------
      function fL1222(aaa1,aaa2)
      implicit none
      include 'pncommon.f'
      double complex fL1222
      double precision aaa1(Naux), aaa2(Naux)

      call Returnv(aaa1,aaa2)
      
      fL1222 = (-(((dlog1mv - dlogv)**2 + Pi2)*(1 - onemv*v)) + 
     >  (2*(dlogv + ie*Pi)*(-9 + 10*onemv*v))/onemv + 
     >  ((dlogv + ie*Pi)**2 + Pi2)*(2 - (v*(4 + 5*v))/onemv2))/8d0

      return 
      end ! fL1222 ->

cpn   ---------------------------------------------------------------------
      function fL2212(aaa1,aaa2)
      implicit none
      include 'pncommon.f'
      double complex fL2212
      double precision aaa1(Naux), aaa2(Naux)

      call Returnv(aaa1,aaa2)
 
      fL2212 = (-(((dlog1mv - dlogv)**2 + Pi2)*(1 - 2*v)**2) + 
     >  (2*(dlogv + ie*Pi)*(-9 + 8*onemv*v))/onemv + 
     > (((dlogv + ie*Pi)**2 + Pi2)*(2*onemv2 - 3*v*(2 + v)))/onemv2)/8d0
      return
 
      end ! fL2212 ->

cpn   ---------------------------------------------------------------------
      function fL1122(aaa1,aaa2)
      implicit none
      include 'pncommon.f'
      double complex fL1122
      double precision aaa1(Naux), aaa2(Naux)

      call Returnv(aaa1,aaa2)
      
      fL1122 = (12 - 4*Pi2 + 6*((dlogv + ie*Pi)**2 + Pi2)*
     >  (8/onemv + 14*onemv + 9/onemv2 - 14*v) + 
     >  (108*(dlogv + ie*Pi)*(1 + 2*onemv*v))/onemv + 
     >  3*((dlog1mv - dlogv)**2 + Pi2)*(-13 + 38*onemv*v) + 
     >  48*((dlog1mv + dlogv)**4/48d0 + 4*fLi4v + 
     >  ((dlog1mv + dlogv)**3*ie*Pi)/12d0 + (dlogv*ie*Pi**3)/2d0 - 
     >  (109*Pi**4)/720d0 + fLi3v*(dlog1mv - 3*dlogv - 2*ie*Pi) - 
     >  (dlogv2*Pi2)/12d0 + fLi2v*((dlogv + ie*Pi)**2 + Pi2))*(-onemv2 -
     >  v2) + 12*((-dlog1mv + dlogv)**3 + 
     >  3*(dlog1mv + ie*Pi)*((dlog1mv - dlogv)**2 + Pi2))*v2 - 
     >  24*(1 + 3*onemv)*v*((dlog1mv - dlogv)*fLi2v1mv + fLi3v1mv + 
     >  (dlog1mv*((dlog1mv - dlogv)**2 + Pi2))/2d0 - zeta3))/48d0
      return 
      end ! fL1122 ->

cpn   ---------------------------------------------------------------------
      function fL1212(aaa1,aaa2)
      implicit none
      include 'pncommon.f'
      double complex fL1212, tmp
      double precision aaa1(Naux), aaa2(Naux)

      call Returnv(aaa1,aaa2)
      
      tmp = 24*onemv*(dlogv + ie*Pi)*((dlog1mv + ie*Pi)**2 + Pi2) + 
     >  (32*Pi2*v)/onemv - (3*((dlog1mv - dlogv)**2 + Pi2)*(-69 + 95*v +
     >  9*v**3 - 27*v2))/onemv - 
     >  (3*((dlog1mv + ie*Pi)**2 + Pi2)*(-27 + 105*v + 103*v**3 - 189*v2
     >  ))/(onemv*v2) + (6*(dlogv + ie*Pi)*
     >  (dlogv2 + 2*dlogv*ie*Pi + (1 + ie**2)*Pi2)*(5 + 4*v + 3*v2))
     >  /onemv2 + (6*(dlog1mv*onemv + ie*Pi + dlogv*v)*(27 - 108*v + 25
     >  *v2))/(onemv*v) - (3*(dlogv2 + 2*dlogv*ie*Pi + (1 + ie**2)*Pi2)
     >  *(15 + 4*v + 29*v2))/onemv2 + (144*(1 + 6*v + 3*v2)*
     >  (dlog1mv**4/24d0 - (dlog1mv**3*dlogv)/6d0 - fLi41mv + fLi4v + 
     >    fLi4v1mv - (7*Pi**4)/360d0 + ((dlog1mv2 + 2*fLi2v)*Pi2)/12d0 - 
     >    (dlog1mv + ie*Pi)*(fLi3v - zeta3)))/onemv2 - 
     >  (144*(1 + v2)*(-dlogv**4/48d0 + fLi4v - Pi**4/90d0 - 
     >  (Pi2*(3*dlogv2 - 6*fLi2v + Pi2))/36d0 + 
     >  ((dlogv + ie*Pi)**2*((dlogv + ie*Pi)**2 + Pi2))/24d0 - 
     >    ((dlogv + ie*Pi)*(fLi3v - zeta3))/2d0))/onemv2 + 
     >  (24*onepv*(ie*Pi**3 + dlog1mv*Pi2 - 6*zeta3))/onemv + 
     >  (2*(-15 + 4*v**3 + 15*v2)*
     >  (6*fLi3v + 3*dlogv2*ie*Pi + ie*Pi**3 - 6*fLi2v*(dlogv + ie*Pi) -
     >  5*dlogv*Pi2 - 6*zeta3))/onemv2 - 
     >  48*(8 + v + (30*v)/onemv)*(fLi31mv + 
     >  ((dlog1mv + ie*Pi)*(-6*fLi21mv + Pi2))/6d0 - 
     >  (dlogv*((dlog1mv + ie*Pi)**2 + Pi2))/2d0 - zeta3)

       tmp = tmp/72d0
       fL1212 = tmp
      return 
      end ! fL1212 ->

cpn   ---------------------------------------------------------------------
      function fSL1222(aaa1,aaa2)
      implicit none
      include 'pncommon.f'
      double complex fSL1222
      double precision aaa1(Naux), aaa2(Naux)

      call Returnv(aaa1,aaa2)
      
      fSL1222 = (-4*(dlogv + ie*Pi)*(-(1/onemv) + v) + 
     >  (((dlogv + ie*Pi)**2 + Pi2)*(1 + v2))/onemv2 + 
     >  (((dlog1mv - dlogv)**2 + Pi2)*(onemv2 + v2))/2d0)/8d0

      return 
      end ! fSL1222 ->

cpn   ---------------------------------------------------------------------
      function fSL1122(aaa1,aaa2)
      implicit none
      include 'pncommon.f'
      double complex fSL1122
      double precision aaa1(Naux), aaa2(Naux)

      call Returnv(aaa1,aaa2)
      
      fSL1122 = (-6 + 2*Pi2 + 4*(-6*fLi4v1mv + fLi2v*Pi2)*(1 - 2*v) - 
     >  (12*(dlogv + ie*Pi)*(1 + 2*onemv*v))/onemv - 
     >  3*((dlog1mv - dlogv)**2 + Pi2)*(3 + 2*onemv*v) + 
     >  (6*((dlogv + ie*Pi)**2 + Pi2)*(1 - 2*v2))/onemv2 - 
     >  48*(-(dlog1mv*dlogv**3)/3d0 + dlogv**4/12d0 + fLi41mv + fLi4v - 
     >  (2*Pi**4)/45d0 - (fLi31mv + fLi3v)*(dlogv + ie*Pi) + 
     >  (dlog1mv*dlogv*Pi2)/12d0 + 
     >  (dlogv*ie*Pi*(-3*dlog1mv*dlogv + dlogv2 + Pi2))/6d0)*v2 + 
     >  24*v*((-dlog1mv + dlogv)*(dlogv2 + fLi2v1mv) + 2*fLi3v -
     >  fLi3v1mv - 2*fLi2v*(dlogv + ie*Pi) - ie*Pi*(dlog1mv2 + Pi2) + 
     >  ((-5*dlog1mv + 5*dlogv + 18*ie*Pi)*((dlog1mv - dlogv)**2 
     >  + Pi2))/12d0 - (2*dlogv*(dlogv2 + Pi2))/3d0 - 3*zeta3))/24d0

      return 
      end ! fSL1122 ->

cpn   ---------------------------------------------------------------------
      function fSL1212(aaa1,aaa2)
      implicit none
      include 'pncommon.f'
      double complex fSL1212
      double precision aaa1(Naux), aaa2(Naux)

      call Returnv(aaa1,aaa2)
      
      fSL1212 = (-6d0 + 2*Pi2 - (6d0*(onemv2 - 2d0*v)*
     >  (dlog1mv*onemv + ie*Pi 
     >  + dlogv*v))/(onemv*v) - 3d0*((dlogv + ie*Pi)**2 + Pi2)*(3d0 - 
     >  2d0*v/onemv2) - 3d0*(2 - onemv2)*((dlog1mv + ie*Pi)**2 + Pi2)/v2
     >  - 3d0*((dlog1mv - dlogv)**2 + Pi2)*(-onemv2 + 2d0*v2) - 
     >  (24d0*(1d0 + v2)*(-fLi41mv + fLi4v1mv + (7d0*Pi**4)/360d0 + 
     >  (dlog1mv**4 - 4d0*dlog1mv**3*dlogv + 
     >  dlogv**3*(dlogv + 2d0*ie*Pi) + 2d0*dlog1mv2*Pi2)/24d0 + 
     >  ((-2d0*dlog1mv + dlogv - ie*Pi)*(fLi3v - zeta3))/2d0))/onemv2 + 
     >  (24d0*(-1d0 - v)*(-dlogv**4/48d0 + fLi4v - Pi**4/90d0 - 
     >  (Pi2*(3d0*dlogv2 - 6d0*fLi2v + Pi2))/36d0 - 
     >  ((dlogv + ie*Pi)*(fLi3v - zeta3))/2d0))/onemv + 
     >  24d0*(1 + (2*v)/onemv)*(fLi31mv + fLi2v*(dlog1mv + ie*Pi) - 
     >  (dlogv2*(dlogv + 3d0*ie*Pi))/8d0 + 
     >  dlogv*(2d0*dlog1mv2 + Pi2)/4d0- zeta3) +12d0*(-1d0 + 2*v/onemv)*
     >  (-dlogv**3/6d0 + fLi3v - fLi2v*(dlogv + ie*Pi) - 
     >  ((dlog1mv + dlogv)*Pi2)/3d0 + zeta3))/12d0

      return 
      end ! fSL1212 ->

CsB ------------------------------------------------------------------------
      Subroutine Quit
CsB ------------------------------------------------------------------------
      Print*, ' Stopped in Quit'
      Stop
      End
CsB_________________________________________________________________________
C
C   End: h^(1) function for gg->AA
CsB_________________________________________________________________________
C ==========================================================================
C ==========================================================================
C --------------------------------------------------------------------------
      SUBROUTINE AARESPH(Y,WT)
C --------------------------------------------------------------------------
CCPY March 2010 ----------------------------------------------------------------------

cpn2005 ----------------------------------------------------------------------
C THIS IS A WRESPH-like subroutine that returns the phase space for the
C AG process, with integration of the matrix element squared over the
c partonic light-cone momentum fractions x1 and x2. The implementation
c is described in my 2005 notes on diphoton production.
c
c 
c WT here corresponds to dPhi= 1/(64S) * 1/shat * 1/(2*pi)**4 
c                              * dQ^2 dy dQT dOmega 
c                              * dx1 dx2 delta(shat+that +uhat-Q^2),
c with x2 being integrated out in terms of x1 and other parameters by
c using the delta-function
CCPY March 2010 ---------------------------------------------------------------------
C I added the factor beta/Y_FLUX to WT.

      IMPLICIT NONE

      REAL*8 W,QCDL,X1X2MIN
      INTEGER NJET,N_PRIMARY_PART,NF,I_STR,IBEAM,NPART,NORDER,
     &  NEVENT,N_SECOND_PART
      REAL*8 HIGGS_MASS,PLAB,PCEM,WHAT,
     &X1,X2,x1hat,x2hat,
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
      COMMON/x1x2hat/X1hat,X2hat

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
CCAO! The read W-boson mass and width
      Real*8 XMW_IN,XMW_GRID,GAMW_IN
      COMMON/WMASS_IN/XMW_IN,XMW_GRID,GAMW_IN
      Real*8 XMZ_IN,XMZ_GRID,GAMZ_IN
      COMMON/ZMASS_IN/XMZ_IN,XMZ_GRID,GAMZ_IN

      REAL*8 YM12MIN,YM12MAX,QTMin,QTMax, TMmin, Tmmax, ETmin
      COMMON/MYYM12/YM12MIN,YM12MAX,QTMin,QTMax, TMmin, Tmmax, ETmin

      REAL*8 Y,WT
      DIMENSION Y(50)
      INTEGER I
      REAL*8 RN
      EXTERNAL RN

      REAL*8 Y_FLUX
      INTEGER*4 IMONTE
      REAL*8 B_QT2,C_QT2,Q2_V,WT_Q2,QT2_V,WT_QT2,WT_YB
      REAL*8 WT_THE,WT_PHI,CTHE,STHE,CPHI,SPHI
      REAL*8 PI,ZTEMP,YM1,YM2,EPS
      REAL*8 ONE,YMW,YGW,YRSH,YMREST,YMW_MIN,YMW_MAX,YMW_VAR,
     >YMW_NORM,YMRES2,YMRES,WT_RES,YBMAX,TAU,YBOOST
      REAL*8 TM2_Q,TM_Q,PHI_Q,R2INV
      REAL*8 P_L1,P_L2,P_Q,P_X,P_Y,P_Z,P_NBAR,P_N,R, P_L1CS
      DIMENSION P_L1(5),P_L2(5),P_Q(5),P_X(5),P_Y(5),
     >P_Z(5),P_NBAR(5),P_N(5),R(2), P_L1CS(5)
      REAL*8 Q_V,QT_V,YB_V,THE_STA,PHI_STA
      COMMON/CSFRAME/Q_V,QT_V,YB_V,THE_STA,PHI_STA

      REAL*8 QT_V_MIN, QT_V_MAX,YB_V_MAX,Q_V_MIN,Q_V_MAX,QT_V_MXY
      COMMON/HARDCUT/ QT_V_MIN, QT_V_MAX,YB_V_MAX,
     >  Q_V_MIN,Q_V_MAX,QT_V_MXY

CCPY  USED FOR SEPERATING NLO AND NLO_SIG CONTRIBUTIONS
      REAL*8 QT_V_SEP
      COMMON/NLO/ QT_V_SEP

      REAL*8 T_NBAR,T_N
      DIMENSION T_NBAR(4),T_N(4)
      LOGICAL TESTING
      REAL*8 QT2_V_MIN,QT2_V_MAX,TEMP0,TEMP1,TEMP
      REAL*8 Y_EMAX
      Character*40 Boson
      Common / Boson / Boson
!      Integer N_SUD_AB,N_SUD_A,N_SUD_B,N_WIL_C,I_SUD
!      Common / Sudakov / N_SUD_AB,N_SUD_A,N_SUD_B,N_WIL_C,I_SUD
      REAL*8 PCM_L1,ECM_L1,ECM_L2,BETA
      Integer ICalls
      Common / Calls / ICalls
      Integer iYPiece, iReadIn
      Common / YPiece / iYPiece, iReadIn
      Data TESTING /.False./
      INTEGER*4 I_DECAY,I_HADRONIZE
      COMMON/DEC_AY/I_DECAY,I_HADRONIZE

cpn2005 for integration over x1, x2
      double precision rtaup, xplus,xminus,qtWrat, dz1min, dz1max, 
     >  alpha1, z1, wt_z1, zMplus, expy
      double precision W2, W2hat

      Character*2 hDecayMode
      Common / hDM / hDecayMode

      Real*8 Reson,xNum
      External Reson
      
      Double Precision ScaPro
      External ScaPro

C      TESTING=.TRUE.
C      If (Testing) Print*, ' in AARESPH '

      PI=4.0*ATAN(1.0D0)
      EPS=1.D-4
CsB * At LHC:
      If (W.GT.10.d3) EPS=1.D-3
      WT=-1.0

      Q_V_MIN=YM12MIN
      Q_V_MAX=YM12MAX

CCPY This subrouine is for (NPART.EQ.12 .or.(NPart.ge.21 .and. NPart.le.26)) 
      QT_V_MIN=QT_V_SEP

CsB For processes requiring grids:
C   QT_V_MAX and YB_V_MAX are determined in ReadIn and passed through COMMON
C   QT_V_MAX is the maximal qT for the resummed grid.

C FIND Q2_V
      If (Boson.Eq.'H0' .And. hDecayMode.Eq.'ZZ' .And. 
     .    I_DECAY.EQ.1) then ! H -> Z Z -> leptons
        xNum = 9.d0
        YM1 = Reson(XMZ,GamZ,xNum)
        YM2 = Reson(XMZ,GamZ,xNum)
      Else
        YM1=XMASS(1)
        YM2=XMASS(2)
      EndIf

CsB * Set VB (3&4) mass and width
      YMW=0.d0
      YGW=0.d0
      If (Boson.Eq.'W+' .or. Boson.Eq.'W-' .or.
     >    Boson.Eq.'W+_RA' .or. Boson.Eq.'W-_RA' .or.
     >    Boson.Eq.'W+A(SV)' .or. Boson.Eq.'W-A(SV)' .or.
     >    Boson.Eq.'HP' .or. Boson.Eq.'HM') then
        YMW=XMW
        YGW = GAMW_IN
        If (GAMW_IN.LT.0.) YGW = GAMW_NLO
      ELSE IF (Boson.Eq.'Z0' .or. Boson.Eq.'Z0_RA' .or.
     .         Boson.Eq.'DY' .or. Boson.Eq.'HZ') then
        YMW=XMZ
        YGW = GAMZ_IN
        If (GAMW_IN.LT.0.) YGW = GAMZ_NLO
      ELSE IF (Boson.Eq.'H0') then
        YMW=XMHIGGS
        YGW=GAMH
      ELSE IF (Boson.Eq.'A0' .or. Boson.Eq.'GL'
     .    .or. Boson.Eq.'AA' .or. Boson.Eq.'AG'
     .    .or. Boson.Eq.'WW_UUB' .or. Boson.Eq.'WW_DDB'
     .    .or. Boson.Eq.'ZZ' .or. Boson.Eq.'ZG') then
        YMW=0.d0
        YGW=0.d0
      ENDIF
!     Print*, 'YMW,YGW =', YMW,YGW 

CsB___Choose a scheme to generate PS
CCPY Oct 2008, INCLUDE 'DRELL-YAN PAIR 'DY' PROCESS
      If (Boson.Eq.'W+' .or. Boson.Eq.'W-' .or.
     .    Boson.Eq.'W+_RA' .or. Boson.Eq.'W-_RA' .or.
     .    Boson.Eq.'W+A(SV)' .or. Boson.Eq.'W-A(SV)' .or.
     .    Boson.Eq.'Z0' .or. Boson.Eq.'Z0_RA' .or. 
     .    Boson.Eq.'H0' .or. Boson.Eq.'DY') then
        IMONTE=2
      else if (Boson.Eq.'AA'.or.Boson.Eq.'AG') then
        imonte=3
      else if (Boson.Eq.'ZZ') then
        imonte=3
      else if (Boson.Eq.'WW_UUB'.or.Boson.Eq.'WW_DDB') then
        imonte=1
      else if (Boson.Eq.'A0') then
        imonte=1
      ELSE
C        IMONTE=1
        IMONTE=0
      END IF

      IF(IMONTE.EQ.0.or.imonte.eq.3) THEN ! Use linear shape ---------------
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

        QT2_V_MAX=QT_V_MAX**2
        QT2_V_MIN=QT_V_MIN**2

CCPY Hard-wired here:
        IMONTE=3

        If (iMonte.Eq.0) then   ! assumes d(sigma)/d(QT^2) = 1/X * LOG(X) 
                                !WITH X = QT^2/Q^2            
          TEMP0=LOG(QT2_V_MIN/Q2_V)**2
          TEMP=LOG(QT2_V_MAX/Q2_V)**2-TEMP0
          C_QT2=2.0/TEMP
          TEMP1=SQRT(TEMP0+TEMP*Y(2))
          QT2_V=Q2_V*EXP(-TEMP1)
          QT_V=DSQRT(QT2_V)
          WT_QT2=-1.0/C_QT2/TEMP1*QT2_V
        elseif (iMonte.eq.2) then ! assumes d(sigma)/d(QT^2) ~ (1/QT2_V)
          QT2_V=QT2_V_MIN*(QT2_V_MAX/QT2_V_MIN)**Y(2)
          QT_V=DSQRT(QT2_V)
          WT_QT2=QT2_V*DLOG(QT2_V_MAX/QT2_V_MIN)
        Else if (iMonte.Eq.3) then !optimal qT sampling for AA production
cpn2005 the qT^2 range is stretched as z~ (qT^2)^alpha1,
cpn2005 with alpha1 << 1, and qT^2_min <= qT^2 <= qT^2_max corresponding
cpn2005 to 0 <= z <= 1
          alpha1=0.3
          temp0 = (qt2_v_min/qt2_v_max)**alpha1
          temp1 = temp0 + (1d0-temp0)*y(2)
          temp =  temp1**(1d0/alpha1)
          qt2_v = qt2_v_max*temp
          qt_v =  dsqrt(qt2_v)
          wt_qt2=qt2_v_max*(1d0-temp0)*temp/alpha1/temp1
        End if                  !iMonte

C GET YBOOST
C      Print*, ' CALL YMAXIMUM '
      CALL YMAXIMUM(W,Q_V,QT_V,YBMAX)
C      Print*, ' CALLED YMAXIMUM '
C HARD-WIRED CUT
      IF(YBMAX.GT.YB_V_MAX) YBMAX=YB_V_MAX

      W2 = W**2                 !S=c.m. energy squared
      TM2_Q=Q2_V+QT2_V          !boson transverse mass squared
      TM_Q=SQRT(TM2_Q)
      TAU=Q2_V/W2               !canonical (qT-independent) tau =Q^2/S   
      rtaup = tm_q/W        !square root of qt-dependent tau = MT^2/S  
      YBOOST=YBMAX*(2.0D0*Y(3)-1.0D0)
      YB_V=YBOOST
      WT_YB=2.0D0*YBMAX
      expy=dexp(yboost)
      XPLUS=rtaup*expy     !x_{A} = M_T*exp(y)/S^0.5  
      XMINUS=rtaup/expy   !x_{B} = M_T*exp(-y)/S^0.5

C DO THE LEPTON KINEMATICS IN THE REST FRAME OF W
C       CHANGE THE VARIABLES ( COS(THE_STAR), PHI_STAR ) IN C.M.
C       FRAME TO ( R(1), R(2) ).
      R(1)=Y(4)
      R(2)=Y(5)  
 
      THE_STA=DACOS(2.D0*R(1)-1.D0)
      PHI_STA=2.D0*PI*R(2)
c      PHI_STA=0.D0

C FOR THE_STA
      WT_THE=2.D0
C FOR PHI_STA
      WT_PHI=2.D0*PI

cpn2005 Weight factors corresponding to dQ^2 dQT^2 dy dOmega
      WT=WT_QT2*WT_Q2*WT_YB*WT_THE*WT_PHI

CsB___Check that the decay is kinematically allowed: Q > m1 + m2
      If (Q_V.LT.YM1+YM2) then
!        Print*, ' Q < m1 + m2 '
        WT = 0.
        Return
      EndIf

      
C FOR MASSIVE FINAL STATE PARTICLES (ARISING FROM 2-BODY PHASE SPACE FACTOR).
      If (TESTING) Print*,'THE_STA,PHI_STA',THE_STA,PHI_STA
      CTHE=COS(THE_STA)
      STHE=SIN(THE_STA)
      CPHI=COS(PHI_STA)
      SPHI=SIN(PHI_STA)
      IF(TESTING) Print*,'Q_V,YM1,YM2=',Q_V,YM1,YM2
      PCM_L1=DSQRT((Q_V**2-(YM1+YM2)**2)*(Q_V**2-(YM1-YM2)**2))/
     >2.0/Q_V
      ECM_L1=(Q_V**2+YM1**2-YM2**2)/2.0/Q_V
      ECM_L2=(Q_V**2-YM1**2+YM2**2)/2.0/Q_V
      IF(TESTING) THEN
        PRINT*,'YM1,YM2,PCM_L1,ECM_L1,ECM_L2'
        PRINT*,YM1,YM2,PCM_L1,ECM_L1,ECM_L2
      ENDIF
C      BETA=PCM_L1/ECM_L2
      BETA=2.0*PCM_L1/Q_V

cpn2005 Add integration over x1 and x2 (xi_{+} and xi_{-} in my 2005
cpn2005 notes on diphotons). x2 is expressed in terms of x1 and other 
cpn2005 variables  using the delta-function for momentum conservation.
cpn2005 x1 is calculated as x1=xplus + Mplus, where Mplus=z1^{1/alpha1}
cpn2005 is the variation of x1, expressed in terms of the actual 
cpn2005 integration variable z1 in order to stretch the x1 integration region:
cpn2005 z1=Mplus^alpha1, where alpha1 is a small power of order 0.2-0.3. 
cpn2005 The integration limits are 
cpn2005 qT^2/(W(1-xminus)) <= Mplus <= 1 - xplus, and similarly for z1. 

      alpha1 = 0.2         ! the phase space for x is stretched as z=x**alpha1
                           ! the integration results should be independent
                           !of alpha1 in a substantial range
      qtWrat=qt2_v/W2

      dz1min = (qtWrat/(1d0-xminus))**alpha1
      dz1max = (1-xplus)**alpha1      
      z1 = dz1min + (dz1max-dz1min)*y(6) !effective integration variable for x1
      zMplus=exp(log(z1)/alpha1)         !M_{+} in my notes
      x1 = xplus + zMplus                !xi_{+}
      x2 = xminus + qtWrat/zMplus        !xi_{-}
      x1hat=xplus/x1                     !\hat x_{+}
      x2hat=xminus/x2                    !\hat x_{-}
      if (x1.lt.xplus.or.x1.gt.1d0.or.x2.lt.xminus.or.x2.gt.1d0) then
        print *, 'AARESPH: WT=0 for unphysical x1, x2 = ',x1,x2
        wt=0d0
        return
      endif 
      
      wt_z1 = (dz1max-dz1min)/alpha1/z1        !the Jacobian for z1
      wt=wt*wt_z1        
      

cpn2005 Multiply by the remaining phase space factors
      wt = wt/(2d0*Pi)**4/(64d0*x1*x2*W2)/W2

CCPY MARCH 2010: FOR GENERAL MASSIVE PARTICLES
      wt=wt*beta
C THIS SUBROUTINE does not INCLUDE the FLUX FACTOR
C Y_FLUX=.5/(X1*X2*W*W)
      Y_FLUX=0.5/(X1*X2*W*W)
      wt=wt/Y_FLUX

cpn2005 To choose only one collinear branch
c      if (x1.lt.2*xplus) then
cpn      if (x1.gt.3*xplus.or.x2.gt.3*xminus) then
c        wt = 0d0
c        return 
c      endif


      IF(TESTING) THEN
        PRINT*,'XPLUS,XMINUS,TAU,YBOOST'
        PRINT*,XPLUS,XMINUS,TAU,YBOOST
        PRINT*,'Q_V,QT_V,YB_V,THE_STA,PHI_STA'
        PRINT*,Q_V,QT_V,YB_V,THE_STA,PHI_STA
        PRINT*,'WT_QT2,WT_Q2,WT_YB,WT_THE,WT_PHI,WT'
        PRINT*,WT_QT2,WT_Q2,WT_YB,WT_THE,WT_PHI,WT
      ENDIF

C**************************************
C In the Collins-Soper frame

      R2INV=1.0D0/SQRT(2.0D0)

c Light-cone unit vectors p_nbar and p_n
      P_NBAR(1)=0.D0
      P_NBAR(2)=0.D0
      P_NBAR(3)=R2INV
      P_NBAR(4)=R2INV

      P_N(1)=0.D0
      P_N(2)=0.D0
      P_N(3)=-R2INV
      P_N(4)=R2INV

C GIVE A RANDOM CHOICE OF AZIMUTHAL ANGLE FOR W-BOSON, IN LAB FRAME
      PHI_Q=2.0*PI*RN(2)

CsB   Massive lepton(1) 4 momentum in CS frame
CsB   To test differences between W rest frame being the CS frame and the
C     one with the z direction in the W motion direction.
      P_L1CS(1) = PCM_L1 * SThe * CPhi
      P_L1CS(2) = PCM_L1 * SThe * SPhi
      P_L1CS(3) = PCM_L1 * CThe
      P_L1CS(4) = ECM_L1
      P_L1CS(5) = YM1

c      P_L1CS(5) = Sqrt(-P_L1CS(1)**2 -P_L1CS(2)**2 -P_L1CS(3)**2 +
c     >                  P_L1CS(4)**2)

CsB 1&2 (vector boson) 4 momentum in lab frame
      P_Q(1)=QT_V*COS(PHI_Q)
      P_Q(2)=QT_V*SIN(PHI_Q)
      P_Q(3)=TM_Q*(EXP(YB_V)-EXP(-YB_V))/2.0D0
      P_Q(4)=TM_Q*(EXP(YB_V)+EXP(-YB_V))/2.0D0

      P_Q(5)=SQRT(-P_Q(1)**2-P_Q(2)**2-P_Q(3)**2+
     >        P_Q(4)**2)

cpn2006 X,Y,Z unity vectors in the CS frame, given in the lab frame
      P_X(1)=TM_Q*COS(PHI_Q)/Q_V
      P_X(2)=TM_Q*SIN(PHI_Q)/Q_V
      P_X(3)=P_Q(3)*QT_V/TM_Q/Q_V
      P_X(4)=P_Q(4)*QT_V/TM_Q/Q_V

      P_Z(1)=0d0; P_Z(2)=0d0; P_Z(3)=P_Q(4)/TM_Q; P_Z(4)=P_Q(3)/TM_Q;

      CALL FINDY(Q_V,P_Q,P_X,P_Z,P_Y)

      P_X(5)=-P_X(1)**2-P_X(2)**2-P_X(3)**2+
     >        P_X(4)**2
      P_Y(5)=-P_Y(1)**2-P_Y(2)**2-P_Y(3)**2+
     >        P_Y(4)**2
      P_Z(5)=-P_Z(1)**2-P_Z(2)**2-P_Z(3)**2+
     >        P_Z(4)**2

      IF(TESTING) THEN
        PRINT*,' P_Q, Q_V = ', P_Q,Q_V
        PRINT*,' P_X =',P_X
        PRINT*,' P_Y =',P_Y
        PRINT*,' P_Z =',P_Z
      ENDIF

      IF(ABS(P_Q(5)-Q_V).GT.EPS) THEN
       PRINT*,' ERROR IN P_Q'
C       CALL EXITGracefully
       RETURN
      ENDIF

      IF(ABS(P_X(5)+1.0).GT.EPS .OR. ABS(P_Y(5)+1.0).GT.EPS
     >.OR. ABS(P_Z(5)+1.0).GT.EPS) THEN
       PRINT*,' POTENTIAL ERROR IN MASS: P_X, P_Y  OR  P_Z'
C       Call ExitGracefully
       RETURN
      ELSE
       P_X(5)=SQRT(-P_X(5))
       P_Y(5)=SQRT(-P_Y(5))
       P_Z(5)=SQRT(-P_Z(5))
      ENDIF

      DO 80 I=1,4
CsB   Boost from CS frame to lab frame.
CCPY  Only works for massless particles (i.e. YM1=YM2=0 ).
C         P_L1(I)=Q_V/2.0D0*(P_Q(I)/Q_V+STHE*CPHI*P_X(I)+
C     >   STHE*SPHI*P_Y(I)+CTHE*P_Z(I))
C         P_L2(I)=Q_V/2.0D0*(P_Q(I)/Q_V-STHE*CPHI*P_X(I)-
C     >   STHE*SPHI*P_Y(I)-CTHE*P_Z(I))
C
         P_L1(I)=ECM_L1*P_Q(I)/Q_V+PCM_L1*(STHE*CPHI*P_X(I)+
     >   STHE*SPHI*P_Y(I)+CTHE*P_Z(I))
         P_L2(I)=ECM_L2*P_Q(I)/Q_V+PCM_L1*(-STHE*CPHI*P_X(I)-
     >   STHE*SPHI*P_Y(I)-CTHE*P_Z(I))
80    CONTINUE
CsB   Note: the above is not the same as the one using BOOST_5
!      CALL BOOST_5 (P_L1CS, P_Q, P_L1)
CsB   Check CS -> Lab boost
!      QT = Sqrt(P_Q(1)**2 + P_Q(2)**2)
!      MT = Sqrt(P_Q(5)**2 + QT**2)
!      Print *,' Q,QT,MT = ',P_Q(5),QT,MT
!      Print *,' P_Q = ', P_Q
!      Print *,' P_L1(4)',P_L1(4),1.d0/P_Q(5)/MT*(P_Q(4)*MT*P_L1CS(4) +
!     +  P_Q(4)*QT*P_L1CS(1) + P_Q(3)*P_Q(5)*P_L1CS(3))
!      Print *,' P_L1(1)',P_L1(1),1.d0/P_Q(5)/MT*(QT*MT*P_L1CS(4) + 
!     +  MT**2*P_L1CS(1) )
!      Print *,' P_L1(2)',P_L1(2),P_L1CS(2)
!      Print *,' P_L1(3)',P_L1(3),1.d0/P_Q(5)/MT*(P_Q(3)*MT*P_L1CS(4) +
!     +  P_Q(3)*QT*P_L1CS(1) + P_Q(4)*P_Q(5)*P_L1CS(3))
!      Stop

C BOOSTED TO LAB FRAME
      P_L1(5)=-P_L1(1)**2-P_L1(2)**2-P_L1(3)**2+
     >        P_L1(4)**2
      P_L2(5)=-P_L2(1)**2-P_L2(2)**2-P_L2(3)**2+
     >        P_L2(4)**2

      IF(TESTING) PRINT*,' P_L1(5),P_L2(5) = ',P_L1(5),P_L2(5)

      P_L1(5)=SQRT(ABS(P_L1(5)))
      P_L2(5)=SQRT(ABS(P_L2(5)))

CsB   To test differences between W rest frame being the CS frame and the
C     one with the z direction in the W motion direction.
c      If (.False.) then
c        If (ICalls.Eq.0) then
c          Print*, ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! '
c          Print*, ' In this run the W rest frame is NOT the CS frame '
c          Print*, ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! '
c          ICalls = 1
c        End If
c        Call Boost_5 (P_L1CS, P_Q, P_L1)
c         DO I=1,4
c          P_L2(I)=P_Q(I)-P_L1(I)
c         ENDDO
c         P_L2(5)=YM2
c      End If

CsB   Check lepton masses in Lab frame.
      IF(YM1.GT.EPS .OR. YM2.GT.EPS) THEN
        IF(ABS(P_L1(5)-YM1).GT.P_L1(5)*EPS .OR.
     >  ABS(P_L2(5)-YM2).GT.P_L2(5)*EPS) THEN
          PRINT*,' POTENTIAL ERROR IN MASS: P_L1 OR P_L2 (3)'
CCPY       Call ExitGracefully
          PRINT *,'P_L1 =  ',P_L1
          PRINT *,'P_L2 =  ',P_L2
          PRINT *,'P_Q  =  ',P_Q
          RETURN
        ENDIF
      ELSE
        IF(P_L1(5).GT.4.d0*EPS .OR. P_L2(5).GT.4.d0*EPS) THEN
          PRINT*,' POTENTIAL ERROR IN MASS: P_L1 OR P_L2 (4)'
CCPY       Call ExitGracefully
          PRINT *,'P_L1 =  ',P_L1
          PRINT *,'P_L2 =  ',P_L2
          PRINT *,'P_Q  =  ',P_Q
          RETURN
        ENDIF
      ENDIF


cpn2005 Fill the lab momenta PLAB for two incoming partons [plab(i,1) &
cpn2005 plab(i,2)], final-state particles 1,2 [plab(i,3) & plab(i,4)],
cpn2005 unobserved jet (5) and vector boson 9=1&2
      DO 90 I=1,4
cpn2005 Momenta for two incoming partons
        PLAB(I,1)=X1*W*R2INV*P_NBAR(I)
        PLAB(I,2)=X2*W*R2INV*P_N(I)
C FERMION
        PLAB(I,3)=P_L1(I)
C ANTI-FERMION
        PLAB(I,4)=P_L2(I)
cpn2005 Unobserved final-state jet (gluon 5)
        plab(i,5)=plab(i,1)+plab(i,2)-plab(i,3)-plab(i,4)
C VECTOR BOSON
        PLAB(I,9)=P_Q(I)
90    CONTINUE

      IF(TESTING) THEN
C       PRINT *,'((PLAB(J1,J2),J1=1,4),J2=3,4)'
C       PRINT *,((PLAB(J1,J2),J1=1,4),J2=3,4)
        PRINT *,'P_L1 =  ',P_L1
        PRINT *,'P_L2 =  ',P_L2
        PRINT *,'P_Q  =  ',P_Q
      ENDIF
CsB___Check |PLAB(_,3)|,|PLAB(_,4)| 
!        Print *,' |PLAB(_,3)|,|PLAB(_,4)| = ',
!     >          Sqrt(ScaPro(3,3)),Sqrt(ScaPro(4,4))

CsB___For the case when each of the two final state particles furter decay 
C     into two massles particles, like for H -> Z0 Z0 -> e+ e- mu+ mu-
      IF(Boson.Eq.'H0' .and. hDecayMode.Eq.'ZZ' .and. i_Decay.Eq.1) THEN 
        CALL HZZDECAY (P_Q,P_L1CS,P_L1,P_L2,
     ,P_X,P_Y,P_Z,STHE,CTHE,SPHI,CPHI,Q_V)
      ENDIF

      RETURN
      END ! AARESPH

C=============================================================

C --------------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION RUN_MASS_TOT(Q,NF)
C --------------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NN=6)
      PARAMETER (ZETA3 = 1.202056903159594D0)
      DIMENSION AM(NN),YMSB(NN)

      REAL*8 qmas_u,qmas_d,qmas_s,qmas_c,qmas_b,qmas_t

      REAL*8 SWS,CWS,XMBOT,XMTOP,XMB2,XMT2,GWEAK2
      COMMON/STAND1/SWS,CWS,XMBOT,XMTOP,XMB2,XMT2,GWEAK2

c-------------------------------------------------
      B0(NF) = (33.D0-2.D0*NF)/12D0
      B1(NF) = (102D0-38D0/3D0*NF)/16D0
      B2(NF) = (2857D0/2D0-5033D0/18D0*NF+325D0/54D0*NF**2)/64D0
      G0(NF) = 1D0
      G1(NF) = (202D0/3D0-20D0/9D0*NF)/16D0
      G2(NF) = (1249D0-(2216D0/27D0+160D0/3D0*ZETA3)*NF
     .       - 140D0/81D0*NF**2)/64D0
      C1(NF) = G1(NF)/B0(NF) - B1(NF)*G0(NF)/B0(NF)**2
      C2(NF) = ((G1(NF)/B0(NF) - B1(NF)*G0(NF)/B0(NF)**2)**2
     .       + G2(NF)/B0(NF) + B1(NF)**2*G0(NF)/B0(NF)**3
     .       - B1(NF)*G1(NF)/B0(NF)**2 - B2(NF)*G0(NF)/B0(NF)**2)/2D0
      TRAN(X)=1D0+4D0/3D0*ALPI(X)
      CQ(X,NF)=(2D0*B0(NF)*X)**(G0(NF)/B0(NF))
     .            *(1D0+C1(NF)*X+C2(NF)*X**2)
c-------------------------------------------------

      PI=4.0D0*ATAN(1.0D0)
      
CCPY: INPUT:

      qmas_s=0.2
      qmas_c=1.6
      qmas_b=xmbot    
      qmas_t=xmtop    
            
      AMS=qmas_s
      AMC=qmas_c
      AMB=qmas_b
      AMT=qmas_t

      ACC = 1.D-8
      AM(1) = 0
      AM(2) = 0
      AM(3) = AMS
      AM(4) = AMC
      AM(5) = AMB
      AM(6) = AMT

C      Print*, ' AMS(i) =',AM(1),AM(2),AM(3),AM(4),AM(5),AM(6)

      IF(NF.GE.4)THEN
       XMSB = AM(NF)/TRAN(AM(NF))
       XMHAT = XMSB/CQ(ALPI(AM(NF)),NF)
      ELSE
       XMSB = 0
       XMHAT = 0
      ENDIF

C      Print*, ' ALPI(Mi) =',
C     ,          ALPI(AM(4)),ALPI(AM(5)),ALPI(AM(6))

      YMSB(3) = AMS
      IF(NF.EQ.3)THEN
       YMSB(4) = YMSB(3)*CQ(Alpi(AM(4)),3)/
     .                   CQ(Alpi( 1.D0),3)
       YMSB(5) = YMSB(4)*CQ(Alpi(AM(5)),4)/
     .                   CQ(Alpi(AM(4)),4)
       YMSB(6) = YMSB(5)*CQ(Alpi(AM(6)),5)/
     .                   CQ(Alpi(AM(5)),5)
      ELSEIF(NF.EQ.4)THEN
       YMSB(4) = XMSB
       YMSB(5) = YMSB(4)*CQ(Alpi(AM(5)),4)/
     .                   CQ(Alpi(AM(4)),4)
       YMSB(6) = YMSB(5)*CQ(Alpi(AM(6)),5)/
     .                   CQ(Alpi(AM(5)),5)
      ELSEIF(NF.EQ.5)THEN
       YMSB(5) = XMSB
       YMSB(4) = YMSB(5)*CQ(Alpi(AM(4)),4)/
     .                   CQ(Alpi(AM(5)),4)
       YMSB(6) = YMSB(5)*CQ(Alpi(AM(6)),5)/
     .                   CQ(Alpi(AM(5)),5)
      ELSEIF(NF.EQ.6)THEN
       YMSB(6) = XMSB
       YMSB(5) = YMSB(6)*CQ(Alpi(AM(5)),5)/
     .                   CQ(Alpi(AM(6)),5)
       YMSB(4) = YMSB(5)*CQ(Alpi(AM(4)),4)/
     .                   CQ(Alpi(AM(5)),4)
      ENDIF
      IF(Q.LT.AMC)THEN
       N0=3
       Q0 = 1.D0
      ELSEIF(Q.LE.AMB)THEN
       N0=4
       Q0 = AMC
      ELSEIF(Q.LE.AMT)THEN
       N0=5
       Q0 = AMB
      ELSE
       N0=6
       Q0 = AMT
      ENDIF
      RUN_MASS_TOT = YMSB(N0)*CQ(Alpi(Q),N0)/
     .               CQ(Alpi(Q0),N0)
      RETURN
      END

C ---------------------------------------------------------------------------
      FUNCTION ALPI(AMU)
C ---------------------------------------------------------------------------

      IMPLICIT NONE
      REAL*8 ALPI,AMU
      REAL*8 ALPHAS,PI
      EXTERNAL ALPHAS
      
      PI=4.0D0*ATAN(1.0D0)      
      ALPI=ALPHAS(AMU)
      ALPI=ALPI/PI
      
      RETURN
      END
C==================================================
CCPY OCT 2008: move the part  related to gg-->AA in resbos.f to resbosm.for where it belongs.
c
CsB --------------------------------------------------------------------
C     Begin routines for gg -> AA
CsB --------------------------------------------------------------------

      Function ggAA_A2(s,t,u)
CsB   Returns the gg -> AA amplitude square; cf. Eq(2) hep-ph/0206194

      Implicit None
      Real*8 ggAA_A2,s,t,u,SQ22
      Integer i,j,k,l
      Complex M
      Common / ggAA_Amplitude / M(2,2,2,2)

      ggAA_A2 = 0.d0
      Do 101 i = 1, 2
      Do 101 j = 1, 2
      Do 101 k = 1, 2
      Do 101 l = 1, 2
        ggAA_A2 = ggAA_A2 + Conjg(M(i,j,k,l))*M(i,j,k,l)
  101 Continue
      ggAA_A2 = 2.d0*ggAA_A2 ! M---- = M++++ etc.

CsB   Sum over quark flavors in the loop
      ggAA_A2 = SQ22(Sqrt(s)) * ggAA_A2

CsB   This normalization factor should be checked
      ggAA_A2 = 4.d0*16.d0 * ggAA_A2

      Return
      End
CsB --------------------------------------------------------------------

      Function SQ22(Q)
CsB   Returns the gg -> AA amplitude square; cf. Eq(2) hep-ph/0206194

      Implicit None
      Real*8 SQ22,Q
      Real*8        SWS,CWS,XMBOT,XMTOP,XMB2,XMT2,GWEAK2,
     &              XMCHARM,XMCHARM2,XMHIGGS,GAMH,GAMH2
      COMMON/STAND1/SWS,CWS,XMBOT,XMTOP,XMB2,XMT2,GWEAK2
      COMMON/STAND5/XMCHARM,XMCHARM2,XMHIGGS,GAMH,GAMH2

cpn2005 This is not correct
c$$$      If (Q.lt.xmcharm)   then
c$$$        SQ22 = (6.d0/9.d0)**2
c$$$      Else If(Q.lt.xmbot) then
c$$$        SQ22 = (10.d0/9.d0)**2
c$$$      Else If(Q.lt.xmtop) then
c$$$        SQ22 = (11.d0/9.d0)**2
c$$$      Else
c$$$        SQ22 = (15.d0/9.d0)**2
c$$$      End If

cpn2005 The massless quark approximation in the box diagram works well for 
cpn2005 4 m^2_b << Q^2 << 4 m_t^2. This is checked by comparing the
cpn2005 massless approximation vs. the Pythia box diagram with the 
cpn2005 heavy-quark masses included. It is therefore safer simply to use
cpn2005 the fixed number of flavors N_f=5 within the 2m_b < Q < 2m_t
cpn2005 region
      if (Q.lt.2*xmbot.or.Q.gt.2*xmtop) then
        print *,'STOP in ConstInit2:'
        print *,'Q must be in the range 2*m_bottom < Q < 2*m_top'
        print *,'Q,m_bottom,m_top = ',Q,xmbot,xmtop
        stop
      endif

      SQ22 = (11./9.)**2


      Return
      End
CsB --------------------------------------------------------------------

      Function ggAA_D2(s,t,u)
CsB   Returns the gg -> AA second dipole term; cf. Eq(10) hep-ph/0206194

      Implicit None
      Real*8 ggAA_D2,s,t,u,SQ22, tmp1,tmp2
      Integer i,j,k,l,ErrFlag
      Complex M,PreFac15,Prefac25
      Common / ggAA_Amplitude / M(2,2,2,2)
      Complex SpinProd
      Common / HeliSpinProd / SpinProd(2,5,5)
      Complex SM
      Common / AmplPhases / SM(2,2,2,2)
      REAL*8 Q_V,QT_V,YB_V,THE_STA,PHI_STA
      COMMON/CSFRAME/Q_V,QT_V,YB_V,THE_STA,PHI_STA


c$$$      Call ConstructInvariants(ErrFlag)
c$$$      If (ErrFlag.Eq.1) then
c$$$        ggAA_D2 = 0.
c$$$        Print*, ' ErrFlag =', ErrFlag
c$$$        Return
c$$$      End If
c$$$!     Call TestInvariants
c$$$      Call ConstructSpinProducts
c$$$      Call ConstructAmplPhases
c$$$
c$$$CsB   Calculate prefactors
c$$$      PreFac15 = SpinProd(1,1,5)*SpinProd(2,5,2)*SpinProd(1,2,1)/
c$$$     /           SpinProd(2,1,5)/SpinProd(1,5,2)/SpinProd(2,2,1)
c$$$      PreFac25 = SpinProd(1,2,5)*SpinProd(2,5,1)*SpinProd(1,1,2)/
c$$$     /           SpinProd(2,2,5)/SpinProd(1,5,1)/SpinProd(2,1,2)
c$$$!     Print*, ' Prefac15,Prefac25 =',Prefac15,Prefac25
c$$$
c$$$CsB   Check
c$$$      Do 100 i = 1, 2
c$$$      Do 100 j = 1, 2
c$$$      Do 100 k = 1, 2
c$$$      Do 100 l = 1, 2
c$$$        If ( Real(Conjg(M(i,j,k,l))*M(i,j,k,l)) .LT. 1d-10 ) 
c$$$     .      M(i,j,k,l) = M(3-i,3-j,3-k,3-l)
c$$$  100 Continue
c$$$      tmp1 = 0.d0
c$$$      tmp2 = 0.d0
c$$$      Do 101 i = 1, 2
c$$$      Do 101 j = 1, 2
c$$$      Do 101 k = 1, 2
c$$$      Do 101 l = 1, 2
c$$$        tmp1    = tmp1 - Real( PreFac15*SM(3-i,j,k,l)*SM(i,j,k,l) *
c$$$     *                         Conjg(M(3-i,j,k,l))*M(i,j,k,l) )
c$$$        tmp2    = tmp2 - Real( PreFac25*SM(i,3-j,k,l)*SM(i,j,k,l) *
c$$$     *                         Conjg(M(i,3-j,k,l))*M(i,j,k,l) )
c$$$  101 Continue
c$$$      ggAA_D2 = tmp1 + tmp2
c$$$!     Print*, 'tmp1,tmp2,ggAA_D2=',tmp1,tmp2,ggAA_D2
c$$$!     Stop

CsB   See gg-AAg.nb, Dipole terms section 
cpn2005 Added the dependence on cos(2*phi_sta)

c      ggAA_D2 = -4.d0*dcos(2*phi_sta)*
c     >  Real(M(1,1,2,2)+M(1,2,1,2)+M(2,1,1,2)+1.) 
c     >  - (qt_v/q_v)**2*Real(-M(1,1,2,2)+M(1,2,1,2)+M(2,1,1,2)+1.)
      ggAA_D2 = -4*dcos(2*phi_sta)*
     >  Real(M(1,1,2,2)+M(1,2,1,2)+M(2,1,1,2)+1.) 

CsB   Sum over quark flavors in the loop
      ggAA_D2 = SQ22(Sqrt(s)) * ggAA_D2

CsB   This normalization factor should be checked
      ggAA_D2 = 4.d0*16.d0 * ggAA_D2

      Return
      End
CsB --------------------------------------------------------------------

      Subroutine ggAA_A(s,t,u)
CsB   Returns the gg -> AA amplitude square; cf. Eq(3) in hep-ph/0206194

      Real*8 s,t,u,Pi
      Complex M
      Common / ggAA_Amplitude / M(2,2,2,2)
      Complex ii,id
      Data ii,id / (0d0,1d0),(0d0,1d-14) /
      Data Pi / 3.1415926535897932 /

      Do 100 i = 1, 2
      Do 100 j = 1, 2
      Do 100 k = 1, 2
      Do 100 l = 1, 2
        M(i,j,k,l) = 0.d0
  100 Continue

      M(2,2,2,2) = 1.d0
      M(1,2,2,2) = 1.d0
      M(2,1,2,2) = 1.d0
      M(2,2,1,2) = 1.d0
      M(2,2,2,1) = 1.d0
 
      M(1,1,2,2) = -1.d0/2.d0*(t**2+u**2)/s**2
     *              *( CLog(Cmplx((t+id)/(u+id)))**2 + Pi**2 )
     -             -(t-u)/s*CLog(Cmplx((t+id)/(u+id))) - 1.d0
 
      M(1,2,1,2) = -1.d0/2.d0*(t**2+s**2)/u**2
     *              * CLog(Cmplx(-(t+id)/(s+id)))**2
     -             -(t-s)/u*CLog(Cmplx(-(t+id)/(s+id))) - 1.d0
     -             -ii*Pi*( (t**2+s**2)/u**2*CLog(Cmplx(-(t+id)/(s+id)))
     +                     +(t-s)/u )
 
      M(2,1,1,2) = -1.d0/2.d0*(u**2+s**2)/t**2
     *              * CLog(Cmplx(-(u+id)/(s+id)))**2
     -             -(u-s)/t*CLog(Cmplx(-(u+id)/(s+id))) - 1.d0
     -             -ii*Pi*( (u**2+s**2)/t**2*CLog(Cmplx(-(u+id)/(s+id)))
     +                     +(u-s)/t )

      Return
      End
CsB --------------------------------------------------------------------

      Subroutine ConstructSpinProducts

      Implicit Double Precision (A-H,O-Z)
      Real*8 pLab,pCEM,wHat
      Common /Mom/ pLab(4,12),pCEM(4,12),wHat
      Real*8 p,s
      Complex pc
      Common / sij / p(4,5),pc(4,5),s(5,5)
      Real*8 pp(5),SinPhi(5,5),CosPhi(5,5)
      Complex SpinProd
      Common / HeliSpinProd / SpinProd(2,5,5)
      Complex ii
      Data ii / (0d0,1d0) /

CsB   Define + components
      Do i = 1, 5
        pp(i) = p(3,i) + p(4,i)
!       Print*, ' pp(',i,') =', pp(i)
      EndDo

CsB   Define helicity spin products
      Do i = 1, 5
      Do j = 1, 5
        Pre = Sqrt(Abs(s(i,j)))
        Den = Pre*Sqrt(pp(i)*pp(j))
        CosPhi(i,j) = (p(1,i)*pp(j)-p(1,j)*pp(i))/Den
        SinPhi(i,j) = (p(2,i)*pp(j)-p(2,j)*pp(i))/Den                 
        SpinProd(1,i,j) = Pre*( CosPhi(i,j)+ii*SinPhi(i,j)) ! <ij>
        SpinProd(2,i,j) = Pre*(-CosPhi(i,j)+ii*SinPhi(i,j)) ! [ij]
!       Print*, 'Pre, Den =',Pre,Den
!       Print*, 'Cos(phi[',i,',',j,']), Sin(phi[',i,',',j,']) =',
!    ,           CosPhi(i,j),SinPhi(i,j)
!       Print*, 'Cos(phi[',i,',',j,'])^2+Sin(phi[',i,',',j,'])^2 =',
!    ,           CosPhi(i,j)**2+SinPhi(i,j)**2
!       Print*, '<',i,',',j,'>, [',i,',',j,'] =',
!    ,          SpinProd(1,i,j),SpinProd(2,i,j)
      EndDo
      EndDo

      Return
      End

CsB --------------------------------------------------------------------

      Subroutine ConstructAmplPhases

      Implicit Double Precision (A-H,O-Z)
      Complex S
      Common / AmplPhases / S(2,2,2,2)
      Complex SpinProd
      Common / HeliSpinProd / SpinProd(2,5,5)
      Complex ii
      Data ii / (0d0,1d0) /

CsB   Eq.(3.5) of hep-ph/0109078

      S(2,2,2,2) = ii*SpinProd(2,1,2)*SpinProd(2,3,4)/
     /                SpinProd(1,1,2)/SpinProd(1,3,4)

      S(1,2,2,2) = ii*SpinProd(1,1,2)*SpinProd(1,1,4)*SpinProd(2,2,4)/
     /                SpinProd(1,3,4)/SpinProd(1,2,3)/SpinProd(1,2,4)

      S(2,1,2,2) = ii*SpinProd(1,2,3)*SpinProd(1,2,4)*SpinProd(2,3,4)/
     /                SpinProd(1,1,4)/SpinProd(1,3,1)/SpinProd(1,3,4)

      S(2,2,1,2) = ii*SpinProd(1,3,2)*SpinProd(1,3,4)*SpinProd(2,2,4)/
     /                SpinProd(1,1,4)/SpinProd(1,2,1)/SpinProd(1,2,4)

      S(1,1,2,2) = ii*SpinProd(1,1,2)*SpinProd(2,3,4)/
     /                SpinProd(2,1,2)/SpinProd(1,3,4)

      S(1,2,1,2) = ii*SpinProd(1,1,3)*SpinProd(2,2,4)/
     /                SpinProd(2,1,3)/SpinProd(1,2,4)

      S(2,1,1,2) = ii*SpinProd(1,2,3)*SpinProd(2,1,4)/
     /                SpinProd(2,2,3)/SpinProd(1,1,4)

      S(1,2,2,2) = 1.d0
      S(2,1,2,2) = 1.d0
      S(2,2,1,2) = 1.d0

      S(1,1,1,1) = S(2,2,2,2)
      S(2,1,1,1) = S(1,2,2,2)
      S(1,2,1,1) = S(2,1,2,2)
      S(1,1,2,1) = S(2,2,1,2)
      S(1,1,1,2) = 1.d0
      S(2,2,1,1) = S(1,1,2,2)
      S(2,1,2,1) = S(1,2,1,2)
      S(1,2,2,1) = S(2,1,1,2)
      S(2,2,2,1) = 1.d0

CsB   Check
!     Do 100 i = 1, 2
!     Do 100 j = 1, 2
!     Do 100 k = 1, 2
!     Do 100 l = 1, 2
!       Print*, '|S(',i,',',j,',',k,',',l,')|=',Abs(S(i,j,k,l))
! 100 Continue
!     Do 110 i = 1, 2
!     Do 110 j = 1, 2
!     Do 110 k = 1, 2
!     Do 110 l = 1, 2
!       Print*, ' S(',i,',',j,',',k,',',l,') =',S(i,j,k,l)
! 110 Continue

      Return
      End
CsB --------------------------------------------------------------------
C     End routines for gg -> AA
CsB --------------------------------------------------------------------
CsB --------------------------------------------------------------------
C     Begin routines for gg -> AAg
CsB --------------------------------------------------------------------

      Subroutine gggAAA2(amp2)
CsB   Returns the gggAA -> 0 amplitude square
      
      Implicit None
      Complex gggAAA
      Common / gggAAAmplitude / gggAAA(2,2,2,2,2,2)
      Complex Appppp,Ampppp,Appppm,Apppmm,Ammppp,Ampppm,tmp
      Real*8 amp2,Nc,Q,pyAlEM,AlphaS,const
      Integer i1,i2,i3,i4,i5,ErrFlag
      Common / Scale / Q

      Real*8        SWS,CWS,XMBOT,XMTOP,XMB2,XMT2,GWEAK2,
     &              XMCHARM,XMCHARM2,XMHIGGS,GAMH,GAMH2
      COMMON/STAND1/SWS,CWS,XMBOT,XMTOP,XMB2,XMT2,GWEAK2
      COMMON/STAND5/XMCHARM,XMCHARM2,XMHIGGS,GAMH,GAMH2
      REAL*8 W,QCDL,X1X2MIN
      INTEGER*4 NJET,N_PRIMARY_PART,NF,I_STR,NPART,NEVENT,NORDER,iBeam,
     &IW,IAPPROX,I_DECAY,N_SECOND_PART
      COMMON/AA/W,QCDL,X1X2MIN,NJET,N_PRIMARY_PART,NF,I_STR,IBEAM,NPART,
     &  NEVENT,NORDER,N_SECOND_PART
      Real*8 Pi
      Data Pi / 3.1415926535897932 /

CsB   Initialization
      Call ConstInit1(1)     ! This is only needed for fLi2...
      Call ConstructInvariants(ErrFlag)
      If (ErrFlag.Eq.1) then
            amp2 = 0.
!       Print*, ' ErrFlag =', ErrFlag
            Return
      End If
      Call ConstructPolarizationVectors(2,3,4,5,1, 3,4,5,1,2)
C     Easy to check gauge invariance, for example, by:
!     Call ConstructPolarizationVectors(3,4,5,1,2, 4,5,1,2,3)
      Call ConstructEpsProducts

CsB   Calculated directly
      gggAAA(2,2,2,1,1,2) =  Apppmm(1,2,3,4,5)
      gggAAA(2,2,2,2,1,2) =  Appppm(1,2,3,4,5)
      gggAAA(2,2,2,2,2,2) =  Appppp(1,2,3,4,5)
      gggAAA(1,2,2,2,2,2) =  Ampppp(1,2,3,4,5)
      gggAAA(1,1,2,2,2,2) =  Ammppp(1,2,3,4,5)
      gggAAA(1,2,2,2,1,2) =  Ampppm(1,2,3,4,5)

CsB   Obtained by switching the photons ...
      gggAAA(2,1,2,2,2,2) = -Ampppp(2,1,3,4,5)
      gggAAA(2,1,2,2,1,2) = -Ampppm(2,1,3,4,5)

CsB   ... or by cyclic relabeling of the gluons ...
      gggAAA(2,2,1,2,2,2) = -Appppm(1,2,4,5,3)
      gggAAA(2,2,2,1,2,2) = -Appppm(1,2,5,3,4)
      gggAAA(2,2,1,2,1,2) = -Apppmm(1,2,4,5,3)
      gggAAA(2,2,1,1,2,2) = -Apppmm(1,2,5,3,4)
      gggAAA(1,2,1,2,2,2) = -Ampppm(1,2,4,5,3)
      gggAAA(1,2,2,1,2,2) = -Ampppm(1,2,5,3,4)
 
CsB   ... or by doing both
      gggAAA(2,1,2,1,2,2) =  Ampppm(2,1,5,3,4)
      gggAAA(2,1,1,2,2,2) =  Ampppm(2,1,4,5,3)

CsB   From Pentagon
!     gggAAA(2,2,2,1,1,2) = gggAAA(1,1,1,2,2,1)*Sqrt(8.d0*Pi**4)
!     gggAAA(1,1,2,2,2,2) =(gggAAA(1,1,2,2,2,1)+gggAAA(2,2,1,1,1,1))/2d0
!    *                    *Sqrt(8.d0*Pi**4)
!     gggAAA(1,2,2,2,1,2) =(gggAAA(1,2,2,2,1,1)+gggAAA(2,1,1,1,2,1))/2d0
!    *                    *Sqrt(8.d0*Pi**4)
!
!     gggAAA(2,1,2,2,1,2) =(gggAAA(2,1,2,2,1,1)+gggAAA(1,2,1,1,2,1))/2d0
!    *                    *Sqrt(8.d0*Pi**4)
!
!     gggAAA(2,2,1,2,1,2) =(gggAAA(2,2,1,2,1,1)+gggAAA(1,1,2,1,2,1))/2d0
!    *                    *Sqrt(8.d0*Pi**4)
!     gggAAA(2,2,1,1,2,2) =(gggAAA(2,2,1,1,2,1)+gggAAA(1,1,2,2,1,1))/2d0
!    *                    *Sqrt(8.d0*Pi**4)
!     gggAAA(1,2,1,2,2,2) =(gggAAA(1,2,1,2,2,1)+gggAAA(2,1,2,1,1,1))/2d0
!    *                    *Sqrt(8.d0*Pi**4)
!     gggAAA(1,2,2,1,2,2) =(gggAAA(1,2,2,1,2,1)+gggAAA(2,1,1,2,1,1))/2d0
!    *                    *Sqrt(8.d0*Pi**4)
!
!     gggAAA(2,1,2,1,2,2) =(gggAAA(2,1,2,1,2,1)+gggAAA(1,2,1,2,1,1))/2d0
!    *                    *Sqrt(8.d0*Pi**4)
!     gggAAA(2,1,1,2,2,2) =(gggAAA(2,1,1,2,2,1)+gggAAA(1,2,2,1,1,1))/2d0
!    *                    *Sqrt(8.d0*Pi**4)

CsB   Obtained by parity symmetry
      gggAAA(1,1,1,1,1,2) =  gggAAA(2,2,2,2,2,2)
      gggAAA(1,1,1,1,2,2) =  gggAAA(2,2,2,2,1,2)
      gggAAA(2,1,1,1,1,2) =  gggAAA(1,2,2,2,2,2)
      gggAAA(2,2,1,1,1,2) =  gggAAA(1,1,2,2,2,2)
      gggAAA(1,1,1,2,2,2) =  gggAAA(2,2,2,1,1,2)
      gggAAA(2,1,1,1,2,2) =  gggAAA(1,2,2,2,1,2)

      gggAAA(1,2,1,1,1,2) =  gggAAA(2,1,2,2,2,2)
      gggAAA(1,1,2,1,1,2) =  gggAAA(2,2,1,2,2,2)
      gggAAA(1,1,1,2,1,2) =  gggAAA(2,2,2,1,2,2)

      gggAAA(1,1,2,1,2,2) =  gggAAA(2,2,1,2,1,2)
      gggAAA(1,1,2,2,1,2) =  gggAAA(2,2,1,1,2,2)

      gggAAA(1,2,1,1,2,2) =  gggAAA(2,1,2,2,1,2)
      gggAAA(1,2,1,2,1,2) =  gggAAA(2,1,2,1,2,2)
      gggAAA(2,1,2,1,1,2) =  gggAAA(1,2,1,2,2,2)
      gggAAA(1,2,2,1,1,2) =  gggAAA(2,1,1,2,2,2)
      gggAAA(2,1,1,2,1,2) =  gggAAA(1,2,2,1,2,2)

      amp2 = 0.d0
      Do 100 i1 = 1, 2
      Do 100 i2 = 1, 2
      Do 100 i3 = 1, 2
      Do 100 i4 = 1, 2
      Do 100 i5 = 1, 2
        amp2 = amp2
     +       + Conjg(gggAAA(i1,i2,i3,i4,i5,2))*gggAAA(i1,i2,i3,i4,i5,2)
  100 Continue

      Nc = 3.
      amp2 = (pyAlEM(Q**2)**2)*(4.*Pi)**2 
     *      *(AlphaS(Q)**3)*(4.*Pi)**3 
     *      * Nc/(Nc**2-1.d0) / Pi**4 * amp2

CsB   Average over the spin of the incoming gluons and 
C     account for the two identical photons in the final state
      amp2 = amp2/4.d0/2.d0

C     Sum over quark flavors in the loop - this can be done calling SQ22 !
      If (Q.lt.xmcharm)   then
            const = (6./9.)**2
      Else If(Q.lt.xmbot) then
        const = (10./9.)**2
      Else If(Q.lt.xmtop) then
            const = (11./9.)**2
      Else
            const = (15./9.)**2
      End If
      amp2 = const*amp2

!     If (amp2.GT.0.001) then
!     If (W.LT.2000. .and. amp2.GT.0.001   .or. 
!    .    W.GT.2000. .and. amp2.GT.0.00025) then
!       If (amp2.GT.0.1) Call TestAs
!       amp2 = 0.
!       Return
!     End If

!     Print*, 'Q,const=',Q,const
!     Print*, 'gggAAA2:',(pyAlEM(Q**2)**2)
!     Print*, 'gggAAA2:',(AlphaS(Q)**3)

  101 Return
      End

CsB --------------------------------------------------------------------

      Subroutine ConstructInvariants(ErrFlag)

      Implicit None
      Integer ErrFlag,i,j,k,l,m
      Real*8 si(10)
      Real*8 Dot
      Real*8 pLab,pCEM,wHat
      Common / Mom / pLab(4,12),pCEM(4,12),wHat
      Real*8 p,s
      Complex pc
      Common / sij / p(4,5),pc(4,5),s(5,5)
C     Needed for gg->AA Asym
      REAL*8 W,QCDL,X1X2MIN
      INTEGER*4 NJET,N_PRIMARY_PART,NF,I_STR,NPART,NEVENT,NORDER,iBeam,
     &N_SECOND_PART
      COMMON/AA/W,QCDL,X1X2MIN,NJET,N_PRIMARY_PART,NF,I_STR,IBEAM,NPART,
     &  NEVENT,NORDER,N_SECOND_PART
      REAL*8 Q_V,QT_V,YB_V,THE_STA,PHI_STA
      COMMON/CSFRAME/Q_V,QT_V,YB_V,THE_STA,PHI_STA
      Real*8 MT,s123,s124

      ErrFlag = 0

CsB   Converting lab momenta to p_i ...
      If (nPart.Eq.12) then      ! gg->AA Asym
C       ... which are all outgoing
        Do j = 1, 4
cpn2005   The momenta of the incoming gluons are defined below
cpn2005          p(j,1) = -pLab(j,1) ! gluon  1 in
cpn2005          p(j,2) = -pLab(j,2) ! gluon  2 in
          p(j,3) =  pLab(j,3) ! photon 1 out
          p(j,4) =  pLab(j,4) ! photon 2 out
cpn2005          p(j,5) =  pLab(j,5) ! gluon    out
        End Do
C       Define the gluon momenta
        MT = Sqrt(Q_V**2+QT_V**2)
cpn2005 
c$$$        s123 = -Q_V**2*Sin(THE_STA)*Sin(PHI_STA)*QT_V/
c$$$     /          (Q_V*Cos(THE_STA)+2.d0*MT)
c$$$        s124 = -s123
c$$$        p(1,1) = -(Q_V-s123)/2.d0 * QT_V*MT/Q_V
c$$$        p(3,1) =  (Q_V-s123)/2.d0 / MT
c$$$        p(4,1) =  (Q_V-s123)/2.d0 / Q_V
c$$$        p(1,2) = -(Q_V-s124)/2.d0 * QT_V*MT/Q_V
c$$$        p(3,2) = -(Q_V-s124)/2.d0 / MT         
c$$$        p(4,2) =  (Q_V-s124)/2.d0 / Q_V

cpn2005 The momenta of the incoming gluons; the minus sign arises
cpn2005 because the gluons are incoming
        p(1,1) = (-1)*(-0.5d0 * QT_V*MT/Q_V)
        p(2,1) =      0d0
        p(3,1) = (-1)*0.5d0*MT
        p(4,1) = (-1)*0.5d0*MT*MT / Q_V
        p(1,2) = (-1)*(-0.5d0 * QT_V*MT/Q_V)
        p(2,2) =  0d0
        p(3,2) = (-1)*(-0.5d0* MT)     
        p(4,2) = (-1)*0.5d0*MT*MT / Q_V

        p(1,5) = -(p(1,1)+p(1,2)+p(1,3)+p(1,4))
        p(2,5) = -(p(2,1)+p(2,2)+p(2,3)+p(2,4))
        p(3,5) = -(p(3,1)+p(3,2)+p(3,3)+p(3,4))
        p(4,5) = -(p(4,1)+p(4,2)+p(4,3)+p(4,4))


      Else If (npart.eq.21.or.nPart.Eq.22) then ! gg->AA Pert
C       ... which are all incoming
        Do j = 1, 4
          p(j,3) =  pLab(j,1) ! gluon  1 in
          p(j,4) =  pLab(j,2) ! gluon  2 in
          p(j,1) = -pLab(j,3) ! photon 1 out
          p(j,2) = -pLab(j,4) ! photon 2 out
          p(j,5) = -pLab(j,5) ! gluon    out
        End Do
      End If

CsB   Creating a set of "complex momenta": this way eps.p is easy to do
      Do i = 1, 4
      Do j = 1, 5
        pc(i,j) = p(i,j)
      End Do
      End Do

CsB   s_{i,j} = 2 p_i . p_j = s_{j,i}
      Do i = 1, 5
      Do j = i, 5
        s(i,j) = 2.d0*Dot(p(1,i),p(1,j))
        s(j,i) = s(i,j)
      End Do
      End Do

CsB   Check for spurious singularities: |s_{i,j}| - |s_{k,l}| ~ 0
      If (nPart.Eq.12) Return
      k = 1
      Do i = 1, 5
        m = i + 1
        Do j = m, 5
        si(k) = s(i,j)
        k = k + 1
        EndDo
      EndDo
      Do k = 1, 10
        m = k + 1
        Do l = m, 10
          If (Abs((Abs(si(k))-Abs(si(l))))/(Abs(si(k))+Abs(si(l))) .LT.
     .        8d-5) then
!         Print*, 'k,si(k)=',k,Abs(si(k)),Abs(si(l))
          ErrFlag = 1
          Return
          End If
        EndDo
      EndDo
      
  101 Return
      End

CsB --------------------------------------------------------------------

      Subroutine ConstructPolarizationVectors(k1,k2,k3,k4,k5, 
     ,                                        q1,q2,q3,q4,q5)
C     Implements Eq.(6.14) of hep-ph/9507456

      Implicit None
      Integer k1,k2,k3,k4,k5, q1,q2,q3,q4,q5, i,j,j1,j2,j3,j4
      Real*8 re,im, cp,Del,den
      Integer k,q
      Dimension cp(4,5),Del(5),k(5),q(5)
      Real*8 p,s
      Complex pc
      Common / sij / p(4,5),pc(4,5),s(5,5)
      Complex eps
      Common / PolarizationVectors / eps(4,5,2)
      Complex ci
      Data ci / (0.d0,1.d0) /
      Real*8 epsilon
      External epsilon

      k(1) = k1
      k(2) = k2
      k(3) = k3
      k(4) = k4
      k(5) = k5
      q(1) = q1
      q(2) = q2
      q(3) = q3
      q(4) = q4
      q(5) = q5

CsB   Contravariant momenta
      Do i = 1, 5
        Do j = 1, 3
          cp(j,i) = -p(j,i)
        EndDo
        cp(4,i) = p(4,i)
      EndDo

      Do i = 1, 5
        Del(i) = s(i,k(i))*
     *           (2.d0*s(i,q(i))*s(k(i),q(i))-s(q(i),q(i))*s(i,k(i)))
     /           /8.d0
        den = DSqrt(2.d0*Del(i))
        If (Del(i).Eq.0.) Go to 101
        Do j1 = 1, 4
          im = 0.d0
          Do 100 j2 = 1, 4
          Do 100 j3 = 1, 4
          Do 100 j4 = 1, 4
            im = im + ! epsilon^{mu q p k}
     +           epsilon(j1,j2,j3,j4)*cp(j2,q(i))*cp(j3,i)*cp(j4,k(i))
  100     Continue
          re = (s(i,k(i))*p(j1,q(i)) - 
     -         s(q(i),k(i))*p(j1,i) - s(q(i),i)*p(j1,k(i)))
     /         /2.d0
          eps(j1,i,1) = (re - ci*im)/den ! eps-
          eps(j1,i,2) = (re + ci*im)/den ! eps+
        End Do
      End Do
      Goto 102

  101 Continue
      Print*, 'Delta(',i,')=0 in ConstructPolarizationVectors'
      Stop
      
  102 Return
      End

CsB --------------------------------------------------------------------

      Subroutine ConstructEpsProducts

      Implicit None
      Real*8 p,s
      Integer i,j,l
      Complex pc
      Common / sij / p(4,5),pc(4,5),s(5,5)
      Complex eps
      Common / PolarizationVectors / eps(4,5,2)
      Complex peps,epep
      Common / epsProducts / peps(5,5,2),epep(5,2,5,2)
      Complex CDot
      
CsB   p_i.eps_j+/- & eps_i+/-.eps_j+/-
      Do 100 i = 1, 5
      Do 100 j = 1, 5
      Do 100 l = 1, 2
        peps(i,  j,l) = CDot( pc(1,i)  ,eps(1,j,l))
        epep(i,1,j,l) = CDot(eps(1,i,1),eps(1,j,l))
        epep(i,2,j,l) = CDot(eps(1,i,2),eps(1,j,l))
  100 Continue

      Return
      End

CsB --------------------------------------------------------------------

      Function I2n(i,j)
CsB   Implements the finite part of Eq.(12) from hep-ph/0312334

      Implicit None
      Complex I2n
      Integer i,j
      Real*8 p,s
      Complex pc
      Common / sij / p(4,5),pc(4,5),s(5,5)
C     Data EulerGamma / 0.5772156649 / 
      Complex idelta
      Data idelta / (0d0,1d-14) /

CsB   s -> s + i*delta
C     I2n = -CLog(Cmplx(-s(i,j))-idelta) - EulerGamma + 2.

CsB   Checked arguments and values against Mathematica
C     Also checked that the constant term cancels in A+++--, A--+++
C     and A-+++-, with a good numeric precision.  This indicates that the 
C     1/epsilon pole also cancels.
      I2n = -CLog(Cmplx(-s(i,j))-idelta)

      Return
      End 

CsB --------------------------------------------------------------------

      Function F1(i,j,k,l,m,n)
CsB   Implements Eq.(44) of hep-ph/0106243

      Implicit None
      Complex F1
      Integer i,j,k,l,m,n
      Double Complex r1,r2
      Double Complex fLi2
C     Complex t1,t2,t3,t4
C     Real*8 rt
      Real*8 p,s
      Complex pc
      Common / sij / p(4,5),pc(4,5),s(5,5)
      Real*8 Pi
      Data Pi / 3.1415926535897932 /
      Double Complex idelta
      Data idelta / (0.d0,1.d-14) /

CsB   Checked against Mathematica.  Up to the different sign of 
C     the imaginary part of Li2, there's an agreement. 
C     rt = s(i,j)/s(m,n)
C     t1 = fLi2(DCmplx(1.d0 - rt))
C     If (rt.LT.0. .and. s(m,n).GT.0.) t1 = Conjg(t1)
C     t3 = fLi2(DCmplx(1.d0-(s(i,j)+delta)/(s(m,n)+delta)))
C     If (t1.NE.t3) Print*, 't1,t3=',t1,t3
C     rt = s(k,l)/s(m,n)
C     t2 = fLi2(DCmplx(1.d0 - rt))
C     If (rt.LT.0. .and. s(m,n).GT.0.) t2 = Conjg(t2)
C     t4 = fLi2(DCmplx(1.d0-(s(k,l)+delta)/(s(m,n)+delta)))
C     If (t2.NE.t4) Print*, 't1,t3=',t1,t3

      r1 = (s(i,j)+idelta)/(s(m,n)+idelta)
      r2 = (s(k,l)+idelta)/(s(m,n)+idelta)
      F1 = - fLi2(1.d0-r1) - fLi2(1.d0-r2)
     -     - CDLog(r1)*CDLog(r2) + Pi**2/6.d0

CsB   Note: there's a delicate cancellation betwen these logs and other 
C     terms within the A--+++, A+++-- and A-+++- amplitudes.

      Return
      End 

CsB --------------------------------------------------------------------

      Function TrFF(i,j,l)
 
      Implicit None
      Complex TrFF
      Integer i,j,l
      Real*8 p,s
      Complex pc
      Common / sij / p(4,5),pc(4,5),s(5,5)
      Complex peps,epep
      Common / epsProducts / peps(5,5,2),epep(5,2,5,2)
  
      TrFF = 2.d0*peps(i,j,l)*peps(j,i,l) - s(i,j)*epep(i,l,j,l)
      
      Return
      End 

CsB --------------------------------------------------------------------

      Function TrFFF(i,j,k,l)

      Implicit None
      Complex TrFFF,CDot
      Integer i,j,k,l
      Real*8 p,s
      Complex pc
      Common / sij / p(4,5),pc(4,5),s(5,5)
      Complex peps,epep
      Common / epsProducts / peps(5,5,2),epep(5,2,5,2)

      TrFFF = - peps(i,j,l)*peps(j,k,l)*peps(k,i,l) 
     +        + peps(i,k,l)*peps(j,i,l)*peps(k,j,l) 
     +        + epep(j,l,k,l)*s(i,j)*peps(k,i,l)/2.d0 
     -        - epep(i,l,k,l)*s(i,j)*peps(k,j,l)/2.d0 
     +        + epep(i,l,j,l)*s(i,k)*peps(j,k,l)/2.d0 
     -        - epep(j,l,k,l)*s(i,k)*peps(j,i,l)/2.d0 
     -        - epep(i,l,j,l)*s(j,k)*peps(i,k,l)/2.d0 
     +        + epep(i,l,k,l)*s(j,k)*peps(i,j,l)/2.d0 

      Return
      End 

CsB --------------------------------------------------------------------

      Function piFjlpk(i,j,l,k)

      Implicit None
      Complex piFjlpk
      Integer i,j,l,k
      Real*8 p,s
      Complex pc
      Common / sij / p(4,5),pc(4,5),s(5,5)
      Complex peps,epep
      Common / epsProducts / peps(5,5,2),epep(5,2,5,2)

      piFjlpk = (s(i,j)*peps(k,j,l) - s(j,k)*peps(i,j,l))/2.d0
      
      Return
      End 

CsB --------------------------------------------------------------------

      Function Appppp(i,j,k,l,m)

      Implicit None
      Complex Appppp
      Integer i,j,k,l,m
      Complex TrFF,TrFFF
      Real*8 p,s
      Complex pc
      Common / sij / p(4,5),pc(4,5),s(5,5)

      Appppp = TrFF(i,j,2)*TrFFF(k,l,m,2)/
     /        (2.d0*s(k,l)*s(l,m)*s(k,m))

      Return
      End 

CsB --------------------------------------------------------------------

      Function Ampppp(i,j,k,l,m)

      Implicit None
      Complex Ampppp
      Integer i,j,k,l,m
      Complex TrFF,piFjlpk,Cmpppp
      Real*8 p,s
      Complex pc
      Common / sij / p(4,5),pc(4,5),s(5,5)

      Ampppp = TrFF(j,k,2)*TrFF(l,m,2)/(s(j,k)**2*s(l,m)**2)*
     *        (Cmpppp(i,j,k,l,m)*piFjlpk(j,i,1,l) - 
     -         Cmpppp(i,j,k,m,l)*piFjlpk(j,i,1,m))

      Return
      End 

CsB --------------------------------------------------------------------

      Function Appppm(i,j,k,l,m)

      Implicit None
      Complex Appppm
      Integer i,j,k,l,m
      Complex TrFF,piFjlpk,Cppppm
      Real*8 p,s
      Complex pc
      Common / sij / p(4,5),pc(4,5),s(5,5)

Checked
      Appppm = TrFF(i,j,2)*TrFF(k,l,2)/(s(i,j)**2*s(k,l)**2)*
     *        (Cppppm(i,j,k,l,m)*piFjlpk(i,m,1,k) - 
     -         Cppppm(i,j,l,k,m)*piFjlpk(i,m,1,l))

      Return
      End 

CsB --------------------------------------------------------------------

      Function AFmmppp(i,j,k,l,m)

      Implicit None
      Complex AFmmppp
      Integer i,j,k,l,m
      Complex TrFF,piFjlpk,CFmmppp,F1
      Real*8 p,s
      Complex pc
      Common / sij / p(4,5),pc(4,5),s(5,5)

      AFmmppp = TrFF(i,j,1)*TrFF(k,l,2)/(s(i,j)**2*s(k,l)**2)*
     *         (CFmmppp(i,j,k,l,m)*piFjlpk(i,m,2,k) - 
     -          CFmmppp(i,j,l,k,m)*piFjlpk(i,m,2,l))*F1(i,k,i,l,j,m)

      Return
      End 

CsB --------------------------------------------------------------------

      Function ABmmppp(i,j,k,l,m)

      Implicit None
      Complex ABmmppp
      Integer i,j,k,l,m
      Complex TrFF,piFjlpk,CBmmppp,I2n
      Real*8 p,s
      Complex pc
      Common / sij / p(4,5),pc(4,5),s(5,5)

      ABmmppp = TrFF(i,j,1)*TrFF(k,l,2)/(s(i,j)**2*s(k,l)**2)*
     *         (CBmmppp(i,j,k,l,m)*piFjlpk(i,m,2,k) - 
     -          CBmmppp(i,j,l,k,m)*piFjlpk(i,m,2,l))*I2n(i,m)

      Return
      End 

CsB --------------------------------------------------------------------

      Function Ammppp(i,j,k,l,m)

      Implicit None
      Complex Ammppp
      Integer i,j,k,l,m
      Complex TrFF,TrFFF,AFmmppp,ABmmppp,A1mmppp
      Real*8 p,s
      Complex pc
      Common / sij / p(4,5),pc(4,5),s(5,5)

      A1mmppp = TrFF(i,j,1)*TrFFF(k,l,m,2)/(2.d0*s(k,l)*s(l,m)*s(k,m))

      Ammppp =  AFmmppp(i,j,k,l,m) - AFmmppp(i,j,k,m,l)  
     -        - AFmmppp(i,j,m,l,k) + AFmmppp(j,i,k,l,m)
     -        - AFmmppp(j,i,k,m,l) - AFmmppp(j,i,m,l,k)
     +        + ABmmppp(i,j,k,l,m) - ABmmppp(i,j,k,m,l)  
     -        - ABmmppp(i,j,m,l,k) + ABmmppp(j,i,k,l,m)
     -        - ABmmppp(j,i,k,m,l) - ABmmppp(j,i,m,l,k)
     +        + A1mmppp

      Return
      End 

CsB --------------------------------------------------------------------

      Function AF1pppmm(i,j,k,l,m)

      Implicit None
      Complex AF1pppmm
      Integer i,j,k,l,m
      Complex TrFF,piFjlpk,CF1pppmm,F1
      Real*8 p,s
      Complex pc
      Common / sij / p(4,5),pc(4,5),s(5,5)
!     Common / Test / iTest

      AF1pppmm = TrFF(i,j,2)*TrFF(l,m,1)/(s(i,j)**2*s(l,m)**2)*
     *          (CF1pppmm(i,j,k,l,m)*piFjlpk(i,k,2,l) + 
     -           CF1pppmm(j,i,k,l,m)*piFjlpk(j,k,2,l))*F1(i,l,j,l,k,m)

!     If (iTest.Eq.1) then
!       Print*, 'F1pppmm:',i,j,k,l,m,
!    ,           CF1pppmm(i,j,k,l,m),
!    ,           CF1pppmm(j,i,k,l,m)
!     End If

      Return
      End 

CsB --------------------------------------------------------------------

      Function AF2pppmm(i,j,k,l,m)

      Implicit None
      Complex AF2pppmm
      Integer i,j,k,l,m
      Complex TrFF,piFjlpk,CF2pppmm,F1
      Real*8 p,s
      Complex pc
      Common / sij / p(4,5),pc(4,5),s(5,5)

      AF2pppmm = TrFF(l,m,1)*TrFF(i,k,2)/(s(l,m)**2*s(i,k)**2)*
     *          (CF2pppmm(i,j,k,l,m)*piFjlpk(i,j,2,l) - 
     -           CF2pppmm(k,j,i,l,m)*piFjlpk(k,j,2,l))*F1(i,l,k,l,j,m)

      Return
      End 

CsB --------------------------------------------------------------------

      Function AB1pppmm(i,j,k,l,m)

      Implicit None
      Complex AB1pppmm
      Integer i,j,k,l,m
      Complex TrFF,piFjlpk,CB1pppmm,I2n
      Real*8 p,s
      Complex pc
      Common / sij / p(4,5),pc(4,5),s(5,5)

      AB1pppmm = TrFF(i,j,2)*TrFF(l,m,1)/(s(i,j)**2*s(l,m)**2)*
     *          (CB1pppmm(i,j,k,l,m)*piFjlpk(i,k,2,l) + 
     -           CB1pppmm(j,i,k,l,m)*piFjlpk(j,k,2,l))*I2n(k,l)

      Return
      End 

CsB --------------------------------------------------------------------

      Function AB2pppmm(i,j,k,l,m)

      Implicit None
      Complex AB2pppmm
      Integer i,j,k,l,m
      Complex TrFF,piFjlpk,CB2apppmm,CB2bpppmm,I2n
      Real*8 p,s
      Complex pc
      Common / sij / p(4,5),pc(4,5),s(5,5)

      AB2pppmm = TrFF(l,m,1)*TrFF(i,k,2)/(s(l,m)**2*s(i,k)**2)*
     *          (CB2apppmm(i,j,k,l,m)*piFjlpk(i,j,2,l) + 
     -           CB2bpppmm(i,j,k,l,m)*piFjlpk(k,j,2,l))*I2n(i,m)

      Return
      End 

CsB --------------------------------------------------------------------

      Function A1pppmm(i,j,k,l,m)

      Implicit None
      Complex A1pppmm
      Integer i,j,k,l,m
      Complex TrFF,piFjlpk,C1pppmm,F1
      Real*8 p,s
      Complex pc
      Common / sij / p(4,5),pc(4,5),s(5,5)
!     Common / Test / iTest

      A1pppmm = TrFF(i,j,2)*TrFF(l,m,1)/(s(i,j)**2*s(l,m)**2)*
     *         (C1pppmm(i,j,k,l,m)*piFjlpk(j,k,2,l) - 
     -          C1pppmm(i,j,k,m,l)*piFjlpk(j,k,2,m))

!     If (iTest.Eq.1) then
!       Print*, '1pppmm:',i,j,k,l,m,
!    ,          C1pppmm(i,j,k,l,m),
!    ,          C1pppmm(i,j,k,m,l)
!     End If

      Return
      End 

CsB --------------------------------------------------------------------

      Function Apppmm(i,j,k,l,m)

      Implicit None
      Complex Apppmm
      Integer i,j,k,l,m
      Complex AFpppmm,AF1pppmm,AF2pppmm,
     ,        ABpppmm,AB1pppmm,AB2pppmm,A1pppmm
      Real*8 p,s
      Complex pc
      Common / sij / p(4,5),pc(4,5),s(5,5)
!     Common / Test / iTest

      AFpppmm = AF1pppmm(i,j,k,l,m) - AF1pppmm(i,j,k,m,l) 
     +        + AF2pppmm(i,j,k,l,m) + AF2pppmm(j,i,k,l,m) 
     -        - AF2pppmm(i,j,k,m,l) - AF2pppmm(j,i,k,m,l)
      ABpppmm = AB1pppmm(i,j,k,l,m) - AB1pppmm(i,j,k,m,l) 
     +        + AB2pppmm(i,j,k,l,m) + AB2pppmm(j,i,k,l,m) 
     -        - AB2pppmm(i,j,k,m,l) - AB2pppmm(j,i,k,m,l)
      Apppmm = AFpppmm + ABpppmm + A1pppmm(i,j,k,l,m)

!     If (iTest.Eq.1) then
!       Print*, 'pppmm:',AFpppmm,ABpppmm,A1pppmm(i,j,k,l,m)
!     End If

      Return
      End 

CsB --------------------------------------------------------------------

      Function AFmpppm(i,j,k,l,m)

      Implicit None
      Complex AFmpppm
      Integer i,j,k,l,m
      Complex TrFF,piFjlpk,F1,
     ,        CF1mpppm,CF2mpppm,CF3mpppm,CF4ampppm,CF4bmpppm
      Real*8 p,s
      Complex pc
      Common / sij / p(4,5),pc(4,5),s(5,5)

Checked
      AFmpppm = TrFF(i,m,1)*TrFF(k,l,2)/(s(i,m)**2*s(k,l)**2)*
     *          CF1mpppm(i,j,k,l,m)*( piFjlpk(i,j,2,k)*F1(i,k,i,l,j,m) 
     -                               -piFjlpk(i,j,2,k)*F1(k,m,l,m,i,j))
     +     + TrFF(i,m,1)*TrFF(j,k,2)/(s(i,m)**2*s(j,k)**2)*
     *      ( (  CF3mpppm(i,j,k,l,m)*piFjlpk(i,l,2,j)
     +         + CF3mpppm(i,k,j,l,m)*piFjlpk(i,l,2,k))*F1(i,j,i,k,l,m)
     +       +( CF4ampppm(i,j,k,l,m)*piFjlpk(i,l,2,j)
     +         +CF4bmpppm(i,j,k,l,m)*piFjlpk(i,l,2,k))*F1(j,m,k,m,i,l) )

      Return
      End 

CsB --------------------------------------------------------------------

      Function ABmpppm(i,j,k,l,m)

      Implicit None
      Complex ABmpppm
      Integer i,j,k,l,m
      Complex TrFF,piFjlpk,I2n,
     ,        CB1mpppm,CB2mpppm,CB3mpppm,CB4ampppm,CB4bmpppm
      Real*8 p,s
      Complex pc
      Common / sij / p(4,5),pc(4,5),s(5,5)

Checked
      ABmpppm = TrFF(i,m,1)*TrFF(k,l,2)/(s(i,m)**2*s(k,l)**2)*
     *          ( CB1mpppm(i,j,k,l,m)*piFjlpk(i,j,2,k)*I2n(i,j)
     +           +CB2mpppm(i,j,k,l,m)*piFjlpk(i,j,2,k)*I2n(j,m) )
     +        + TrFF(i,m,1)*TrFF(j,k,2)/(s(i,m)**2*s(j,k)**2)*
     *         ( (  CB3mpppm(i,j,k,l,m)*piFjlpk(i,l,2,j)
     +            + CB3mpppm(i,k,j,l,m)*piFjlpk(i,l,2,k) )*I2n(l,m)
     +          +( CB4ampppm(i,j,k,l,m)*piFjlpk(i,l,2,j)
     +            +CB4bmpppm(i,j,k,l,m)*piFjlpk(i,l,2,k) )*I2n(i,l) )

!     If (i.Eq.1 .and. j.Eq.2 .and. k.Eq.4 .and. l.Eq.5 .and. m.Eq.3)
!    , then
!     Print*, TrFF(i,m,1)*TrFF(k,l,2)/(s(i,m)**2*s(k,l)**2)*
!    *          ( CB1mpppm(i,j,k,l,m)*piFjlpk(i,j,2,k)*I2n(i,j)
!    +           +CB2mpppm(i,j,k,l,m)*piFjlpk(i,j,2,k)*I2n(j,m) )
!       Print*, TrFF(i,m,1)*TrFF(j,k,2)/(s(i,m)**2*s(j,k)**2)*
!    *         ( (  CB3mpppm(i,j,k,l,m)*piFjlpk(i,l,2,j)
!    +            + CB3mpppm(i,k,j,l,m)*piFjlpk(i,l,2,k) )*I2n(l,m)
!    +          +( CB4ampppm(i,j,k,l,m)*piFjlpk(i,l,2,j)
!    +            +CB4bmpppm(i,j,k,l,m)*piFjlpk(i,l,2,k) )*I2n(i,l) )
!     Print*, TrFF(i,m,1)*TrFF(k,l,2)/(s(i,m)**2*s(k,l)**2)*
!    *          CB1mpppm(i,j,k,l,m)*piFjlpk(i,j,2,k)*I2n(i,j)
!     Print*, TrFF(i,m,1)*TrFF(k,l,2)/(s(i,m)**2*s(k,l)**2)*
!    *          CB2mpppm(i,j,k,l,m)*piFjlpk(i,j,2,k)*I2n(j,m)
!     Print*, TrFF(i,m,1)*TrFF(k,l,2)/(s(i,m)**2*s(k,l)**2)*
!    *          CB1mpppm(i,j,k,l,m)*piFjlpk(i,j,2,k)*I2n(j,m)
!       Print*, i,j,k,l,m,s(i,m),s(k,l),s(j,k)
!     End If

      Return
      End 

CsB --------------------------------------------------------------------

      Function A1mpppm(i,j,k,l,m)

      Implicit None
      Complex A1mpppm
      Integer i,j,k,l,m
      Complex TrFF,piFjlpk,C1mpppm
      Real*8 p,s
      Complex pc
      Common / sij / p(4,5),pc(4,5),s(5,5)

Checked
      A1mpppm = TrFF(i,m,1)*TrFF(k,l,2)/(s(i,m)**2*s(k,l)**2)*
     *          (  C1mpppm(i,j,k,l,m)*piFjlpk(i,j,2,k)
     -           - C1mpppm(i,j,l,k,m)*piFjlpk(i,j,2,l) )
      
      Return
      End 

CsB --------------------------------------------------------------------

      Function Ampppm(i,j,k,l,m)

      Implicit None
      Complex Ampppm
      Integer i,j,k,l,m
      Complex AFmpppm,ABmpppm,A1mpppm
      Real*8 p,s
      Complex pc
      Common / sij / p(4,5),pc(4,5),s(5,5)

      Ampppm =   AFmpppm(i,j,k,l,m) - AFmpppm(i,j,l,k,m) 
     +         + ABmpppm(i,j,k,l,m) - ABmpppm(i,j,l,k,m)
     +         + A1mpppm(i,j,k,l,m)

!     If (i.Eq.1 .and. j.Eq.2 .and. k.Eq.4 .and. l.Eq.5 .and. m.Eq.3)
!    , then
!       Print*, AFmpppm(i,j,k,l,m)
!       Print*, AFmpppm(i,j,l,k,m) 
!       Print*, ABmpppm(i,j,k,l,m)
!       Print*, ABmpppm(i,j,l,k,m)
!       Print*, AFmpppm(i,j,k,l,m)-AFmpppm(i,j,l,k,m) 
!       Print*, ABmpppm(i,j,k,l,m)-ABmpppm(i,j,l,k,m)
!       Print*, AFmpppm(i,j,k,l,m) - AFmpppm(i,j,l,k,m) 
!    +        + ABmpppm(i,j,k,l,m) - ABmpppm(i,j,l,k,m)
!       Print*, A1mpppm(i,j,k,l,m)
!     End If      

      Return
      End 

CsB --------------------------------------------------------------------
C     Routines generated by Mathematica for gg -> AAg
CsB --------------------------------------------------------------------
      
      Function Cmpppp(i,j,k,l,m)
      
      Implicit None
      Complex Cmpppp
      Integer i,j,k,l,m
      Real*8 p,s
      Complex pc
      Common / sij / p(4,5),pc(4,5),s(5,5)
      
      Cmpppp = s(j,k)/s(j,l) - (s(i,m))/s(k,l) - (s(i,m))/s(k,m) - 
     &  (s(i,j)*s(i,m))/(s(j,l)*s(k,m))
      
      Return
      End 
      
CsB &&&&----------------------------------------------------------------
      
      Function Cppppm(i,j,k,l,m)
      
      Implicit None
      Complex Cppppm
      Integer i,j,k,l,m
      Real*8 p,s
      Complex pc
      Common / sij / p(4,5),pc(4,5),s(5,5)

CsB   Checked by 'hand'
      Cppppm = (s(i,k) - s(k,l))/s(j,k) + (s(i,j) + s(k,l))/s(i,m) - 
     &  ((s(j,k) + s(k,m)))/s(i,k) + 
     &  (s(i,k)*s(i,m))/(s(j,k)*s(l,m)) + 
     &  (s(i,m) - s(j,m))/s(l,m) - 
     &  (s(j,k)*s(j,m))/(s(i,k)*s(l,m)) - 
     &  (s(k,l)*s(l,m))/(s(i,m)*s(j,k)) - 
     &  (s(i,k)*s(l,m))/(s(i,m)*s(k,m)) - 
     &  (s(i,k)*s(i,l)*s(l,m))/(s(i,m)*s(j,l)*s(k,m)) + 
     &  s(l,m)**2/(s(i,m)*s(j,l)) - 
     &  ((s(i,j)**2 - s(i,j)*s(l,m) + s(l,m)**2))/
     &   (s(i,m)*s(k,m))
      
      Return
      End 
                        
CsB &&&&----------------------------------------------------------------
      
      Function CFmmppp(i,j,k,l,m)
      
      Implicit None
      Complex CFmmppp
      Integer i,j,k,l,m
      Real*8 p,s
      Complex pc
      Common / sij / p(4,5),pc(4,5),s(5,5)
      
      CFmmppp = (-s(i,l))/s(k,l) - (s(i,l))/s(k,m) + 
     &  (5.D-1*(s(i,j)**2 - 2D0*s(i,k)*s(i,l)))/(s(i,m)*s(k,m))
      
      Return
      End 
      
CsB &&&&----------------------------------------------------------------
      
      Function CBmmppp(i,j,k,l,m)
      
      Implicit None
      Complex CBmmppp 
      Integer i,j,k,l,m
      Real*8 p,s
      Complex pc
      Common / sij / p(4,5),pc(4,5),s(5,5)
      
      CBmmppp = (s(i,m)*s(j,k)**2)/(s(k,m)**2*(s(i,k) + s(k,m))) - 
     &  ((s(i,k) + s(k,m)))/s(i,m) + (s(i,l) + s(j,l))/s(l,m) + 
     &  (s(i,l)*s(k,m))/(s(i,m)*s(l,m)) - 
     &  (s(i,l)*s(l,m))/(s(i,m)*s(k,m)) + 
     &  (s(i,k)*s(l,m)**2)/(s(i,m)*s(k,m)**2) - 
     &  (s(l,m)**2)/(s(i,m)*s(k,m)) + 
     &  (s(i,j) - s(i,l) - s(k,m))/(s(i,l) + s(l,m)) + 
     &  (2D0*s(i,l)*(s(i,m) + s(l,m)))/s(k,m)**2 - 
     &  (2D0*s(j,k)*(s(i,m) + s(l,m)))/(s(k,m)*(s(i,k) + s(k,m))) + 
     &  (2D0*(s(i,m) + s(l,m))**2)/s(k,m)**2 + 
     &  (s(i,k)*(s(i,m) + 2D0*s(l,m)))/s(k,m)**2 - 
     &  ((s(i,m) + 2D0*s(l,m)))/s(k,m) + 
     &  (s(i,m) + 2D0*s(l,m))/(s(i,k) + s(k,m)) + 
     &  (s(l,m)*((s(i,k) + s(k,m))/(s(i,l) + s(l,m)) + 
     &       (s(i,l) + s(l,m))/(s(i,k) + s(k,m))))/s(i,m)
      
      Return
      End 
                        
CsB &&&&----------------------------------------------------------------
      
      Function CF1pppmm(i,j,k,l,m)
      
      Implicit None
      Complex CF1pppmm
      Integer i,j,k,l,m
      Real*8 p,s
      Complex pc
      Common / sij / p(4,5),pc(4,5),s(5,5)
      
      CF1pppmm = (-3D0*s(j,l))/s(i,j) - 
     &  (2D0*s(j,l)*(s(j,k) + s(j,l)))/s(i,j)**2 - 
     &  (5.D-1*(s(i,l)**2 + s(j,l)**2))/(s(i,k)*s(k,l)) - 
     &  (5.D-1*(2D0*s(i,l) + s(k,l)))/s(i,k) - 
     &  ((2D0*s(i,l) - 2D0*s(j,l)))/s(k,m) + 
     &  (2D0*s(j,l)*(2D0*s(j,l) - s(k,l)))/(s(i,j)*s(k,m)) - 
     &  ((s(i,l)**2 + s(j,l)**2))/(s(k,l)*s(k,m)) - 
     &  (s(k,l))/s(k,m) - (2D0*s(j,l)**2*s(k,l))/(s(i,j)**2*s(k,m))
      
      Return
      End 
      
CsB &&&&----------------------------------------------------------------
      
      Function CF2pppmm(i,j,k,l,m)
      
      Implicit None
      Complex CF2pppmm
      Integer i,j,k,l,m
      Real*8 p,s
      Complex pc
      Common / sij / p(4,5),pc(4,5),s(5,5)
      
      CF2pppmm = (-s(k,l))/s(i,k) + (5.D-1*s(k,l)**2)/(s(i,j)*s(j,l)) + 
     &  (5.D-1*(s(k,l) + s(l,m))**2)/(s(i,j)*s(j,l))
      
      Return
      End 
      
CsB &&&&----------------------------------------------------------------
      
      Function CB1pppmm(i,j,k,l,m)
      
      Implicit None
      Complex CB1pppmm
      Integer i,j,k,l,m
      Real*8 p,s
      Complex pc
      Common / sij / p(4,5),pc(4,5),s(5,5)
      
      CB1pppmm = (-(s(i,j) + s(i,k) - s(j,k)))/s(j,k) - 
     &  (2D0*(s(i,k) - 2D0*s(i,l))*s(i,l)*s(j,k))/(s(i,j)*s(i,k)**2) - 
     &  (2D0*s(i,l)**2*s(j,k))/(s(i,j)*s(i,k)**2) - 
     &  ((2D0*s(i,j) - s(i,l))*s(i,l)**2*s(j,k))/
     &   (s(i,k)**2*(s(i,k) + s(i,l))**2) - 
     &  ((4D0*s(i,j) - s(i,l))*s(i,l)**2*s(j,k))/
     &   (s(i,j)*s(i,k)**2*(s(i,k) + s(i,l))) - 
     &  (3D0*s(i,l)**3*s(j,k))/(s(i,j)*s(i,k)**2*(s(i,k) + s(i,l))) + 
     &  (s(i,l)*(4D0*s(i,j) + 5D0*s(i,l))*s(j,k))/
     &   (s(i,k)**2*(s(i,k) + s(i,l))) - 
     &  ((s(i,j) - s(j,l))*s(j,l))/(s(j,k) + s(j,l))**2 + 
     &  s(i,j)/(s(j,k) + s(j,l)) - 
     &  (2D0*s(j,l)**2)/(s(i,j)*(s(j,k) + s(j,l))) - 
     &  (2D0*s(i,l)*s(j,m))/(s(i,j)*(s(i,k) + s(i,l))) - 
     &  (s(i,j)*s(i,l)*s(j,k))/((s(i,k) + s(i,l))**2*s(k,l)) - 
     &  (s(i,l)*s(j,k)**2)/(s(i,k)**2*s(k,l)) + 
     &  ((s(i,m) - s(j,k))*(s(i,k) + s(j,k)))/(s(j,k)*s(k,l)) + 
     &  (s(j,k)*(s(i,j) + s(i,l) + s(j,k)))/
     &   ((s(i,k) + s(i,l))*s(k,l)) - 
     &  (s(j,k)*(s(i,j) + 3D0*s(i,l) + s(j,k)))/(s(i,k)*s(k,l)) - 
     &  (s(i,j)*s(j,l)**2)/((s(j,k) + s(j,l))**2*s(k,l)) + 
     &  ((s(i,j) - s(i,l))**2*s(i,l)*s(k,l))/
     &   (s(i,k)**2*(s(i,k) + s(i,l))**2) - 
     &  (2D0*s(i,j)*(s(i,j) - 2D0*s(i,l))*s(k,l))/
     &   (s(i,k)**2*(s(i,k) + s(i,l))) + 
     &  (2D0*s(i,j) + 2D0*s(j,k) + 2D0*s(k,l))/s(i,k) - 
     &  ((2D0*s(i,l)*s(j,k) + 2D0*s(i,j)*s(k,l) + s(i,l)*s(k,l)))/
     &   s(i,k)**2 + (s(j,l)*(s(i,j) + 2D0*s(j,l) + s(l,m)))/
     &   ((s(j,k) + s(j,l))*s(k,l))
      
      Return
      End 
      
CsB &&&&----------------------------------------------------------------
      
      Function CB2apppmm(i,j,k,l,m)
      
      Implicit None
      Complex CB2apppmm
      Integer i,j,k,l,m
      Real*8 p,s
      Complex pc
      Common / sij / p(4,5),pc(4,5),s(5,5)
      
      CB2apppmm = (5D0*s(k,l))/s(j,l) - (7D0*s(k,l))/(s(j,k) + s(j,l)) + 
     &  (s(i,j)*(s(j,l) - s(k,l))*s(k,l))/
     &   (s(j,k)*(s(j,k) + s(j,l))**2) - 
     &  (2D0*s(j,l)*s(k,l)**2)/s(i,j)**3 - 
     &  (s(j,l)*s(k,l)**2)/(s(i,j)**2*s(j,k)) + 
     &  (s(j,l)**3*s(k,l)**2)/(s(i,j)**2*s(j,k)*(s(j,k) + s(j,l))**2) + 
     &  (2D0*s(i,j)*s(k,l)**2)/(s(j,k)*s(j,l)*(s(j,k) + s(j,l))) + 
     &  (2D0*s(j,l)**2*s(k,l)**2)/(s(i,j)**3*(s(j,k) + s(j,l))) - 
     &  (s(k,l))/(s(j,k) + s(k,l)) - 
     &  (s(j,l)*s(k,l)*(3D0*s(i,k) - s(j,l) + s(k,l)))/
     &   (s(i,j)**2*(s(j,k) + s(j,l))) - 
     &  (s(i,k)*s(k,l)*(s(j,l) + s(k,l)))/
     &   (s(j,k)*s(j,l)*(s(j,k) + s(k,l))) + 
     &  (s(j,l)*s(k,l)*(2D0*s(i,k) - 2D0*s(j,l) + 3D0*s(k,l)))/
     &   (s(j,k)*(s(j,k) + s(j,l))**2) - 
     &  (s(j,l)**2*s(k,l)*(2D0*s(i,k) - s(j,l) + 3D0*s(k,l)))/
     &   (s(i,j)*s(j,k)*(s(j,k) + s(j,l))**2) - 
     &  (s(k,l)*(4D0*s(i,k) - 3D0*s(j,l) + 6D0*s(k,l)))/
     &   (s(i,j)*(s(j,k) + s(j,l))) + 
     &  (s(k,l)*(3D0*s(i,j) + 5D0*s(i,k) + 7D0*s(k,l)))/
     &   (s(j,l)*(s(j,k) + s(j,l))) - 
     &  (s(k,l)*(2D0*s(j,k) + s(j,l) - 3D0*s(k,l) - 3D0*s(l,m)))/
     &   s(i,j)**2 + (s(k,l)*(2D0*s(i,j) - 3D0*s(k,l) - 2D0*s(l,m)))/
     &   (s(j,k)*s(j,l)) - (s(k,l)*
     &     (6D0*s(j,k) + s(j,l) - 3D0*s(k,l) - 2D0*s(l,m)))/
     &   (s(i,j)*s(j,k)) - (s(l,m)**2)/s(i,j)**2 + 
     &  ((s(i,j) + 2D0*s(j,k) + 2D0*s(j,l) - s(l,m))*s(l,m)**2)/
     &   (s(i,j)**2*s(k,l)) - 
     &  (s(j,k)*s(l,m)**3)/(s(i,j)**2*s(j,l)*s(k,l)) - 
     &  (2D0*s(j,k)**2*s(l,m)**3)/(s(i,j)**3*s(j,l)*s(k,l)) - 
     &  (2D0*(2D0*s(j,k) + s(j,l))*s(l,m)**3)/(s(i,j)**3*s(k,l)) + 
     &  (2D0*s(j,k)*(s(k,l)**2 - 3D0*s(l,m)**2))/s(i,j)**3 - 
     &  (2D0*s(j,k)**2*(s(k,l)**2 + 3D0*s(k,l)*s(l,m) + 3D0*s(l,m)**2))/
     &   (s(i,j)**3*s(j,l)) - 
     &  (s(j,k)*(4D0*s(k,l)**2 + 9D0*s(k,l)*s(l,m) + 
     &       6D0*s(l,m)**2))/(s(i,j)**2*s(j,l))
      
      Return
      End 
      
CsB &&&&----------------------------------------------------------------
      
      Function CB2bpppmm(i,j,k,l,m)
      
      Implicit None
      Complex CB2bpppmm
      Integer i,j,k,l,m
      Real*8 p,s
      Complex pc
      Common / sij / p(4,5),pc(4,5),s(5,5)
      
      CB2bpppmm = (-(2D0*s(i,k) + s(j,k)))/s(j,k) + 
     &  (s(i,l)*(3D0*s(i,l) + s(j,l)))/(s(i,j)*s(j,l)) - 
     &  (2D0*s(i,l))/(s(j,k) + s(j,l)) - 
     &  (s(i,l)*(2D0*s(i,l) - s(k,l)))/
     &   (s(i,j)*(s(j,k) + s(j,l))) - 
     &  (s(i,l)*(s(i,l) + s(j,l))**2)/(s(i,j)*s(j,l)*s(k,l)) + 
     &  (s(i,j)*s(i,l)*s(k,l))/(s(j,k)*s(j,l)**2) - 
     &  (s(i,j)*s(i,l)*s(k,l))/(s(j,k)**2*s(j,l)) - 
     &  (2D0*s(i,l)*s(j,k)*s(k,l))/(s(i,j)**2*s(j,l)) + 
     &  (s(i,l)*s(j,l)*s(k,l))/(s(i,j)*(s(j,k) + s(j,l))**2) + 
     &  (2D0*s(i,l)*s(j,k)*s(k,l))/(s(i,j)**2*(s(j,k) + s(j,l))) - 
     &  (s(i,j)*s(i,l)*s(k,l))/(s(j,l)**2*(s(j,k) + s(j,l))) - 
     &  (s(i,j)**2*s(k,l)**2)/(s(j,k)**2*s(j,l)**2) + 
     &  (s(i,j)**2*s(k,l)**2)/(s(j,l)**2*(s(j,k) + s(j,l))**2) - 
     &  (2D0*s(i,j)*s(k,l)*(s(i,j) + s(k,l)))/(s(j,k)**2*s(j,l)) - 
     &  ((s(i,j) + s(j,l))*s(k,l)**2)/
     &   (s(j,k)**2*(s(j,k) + s(k,l))) - 
     &  ((s(i,k)**2 + 2D0*s(i,j)*s(k,l)))/s(j,k)**2 - 
     &  (s(i,j)*s(k,l)*(s(i,j) - 3D0*s(k,m)))/
     &   (s(j,l)*(s(j,k) + s(j,l))**2) + 
     &  (s(j,l)*(s(j,l) + 2D0*s(k,l) - 2D0*s(k,m)))/
     &   (s(j,k) + s(j,l))**2 - 
     &  (2D0*(s(i,j) + s(j,l))*s(k,l)*s(k,m))/
     &   (s(j,k)**2*(s(j,k) + s(k,l))) - 
     &  (s(i,j)*s(k,l)**2*s(k,m))/
     &   (s(j,k)**2*s(j,l)*(s(j,k) + s(k,l))) + 
     &  s(k,m)**2/(s(j,k) + s(j,l))**2 - 
     &  (s(j,l)*s(k,m)**2)/(s(j,k)**2*(s(j,k) + s(k,l))) - 
     &  (s(k,l)*(s(i,j) + 2D0*s(k,m)))/(s(j,k) + s(j,l))**2 - 
     &  (4D0*s(i,l)*s(j,k)*s(l,m))/(s(i,j)**2*s(j,l)) - 
     &  (2D0*s(i,l)*s(l,m)**2)/(s(i,j)**2*s(k,l)) - 
     &  (2D0*s(i,l)*s(j,k)*s(l,m)**2)/(s(i,j)**2*s(j,l)*s(k,l))
      
      Return
      End 
      
CsB &&&&----------------------------------------------------------------
      
      Function C1pppmm(i,j,k,l,m)
      
      Implicit None
      Complex C1pppmm
      Integer i,j,k,l,m
      Real*8 p,s
      Complex pc
      Common / sij / p(4,5),pc(4,5),s(5,5)
      Integer i1,i2
      Integer iTest
      Common / Test / iTest

Checked
      C1pppmm = (s(i,l) - s(i,m))/s(i,k) + 
     &  (s(i,l)*(s(i,l) + 2D0*s(i,m)))/(s(i,k)*(s(i,k) + s(i,l))) + 
     &  s(i,l)**2/(s(i,k)*s(j,k)) - 
     &  ((s(i,l) + s(i,m))*(s(i,l) - 3D0*s(k,l)))/(s(i,k)*s(j,k)) - 
     &  (s(i,j)*(2D0*s(i,k) + s(i,l) + s(i,m)))/
     &   (s(j,k)*(s(i,m) - s(k,l))) - 
     &  ((s(i,k) + s(i,m) - 3D0*s(j,k) - s(k,l)))/s(i,l) + 
     &  (2D0*s(i,j) + s(i,k))/s(k,l) + 
     &  (s(i,j)**2 + 2D0*s(i,j)*s(i,k) + s(i,k)**2 + s(i,j)*s(j,m) + 
     &     s(i,k)*s(j,m) + s(j,m)**2)/(s(j,k)*s(k,l)) + 
     &  ((s(i,k) - s(i,m))**2*(s(i,k) + s(k,l)))/
     &   (s(i,k)*s(i,l)*s(j,k)) - 
     &  ((s(i,k) - s(i,m))**2*(s(i,k) + s(k,l)))/
     &   (s(i,l)*(s(i,k) + s(i,l))*s(j,k)) + 
     &  ((s(i,k) - s(i,l))**2*(s(i,k) + s(k,l)))/
     &   (s(i,k)*(s(i,k) + s(i,m))*s(j,k)) - 
     &  ((s(i,k) - 2D0*s(i,l) - 3D0*s(i,m) - 3D0*s(k,m)))/s(j,k) + 
     &  (s(i,k)*s(j,k))/(s(i,l)*(s(i,l) - s(k,m))) - 
     &  (s(i,l)*s(i,m)**2)/(s(i,k)*s(j,k)*s(k,m)) - 
     &  ((s(i,k) + s(i,l))*s(j,m))/(s(i,l)*s(k,m)) - 
     &  (s(i,l)*s(i,m)*(s(i,k) - s(i,m) + s(k,l)))/
     &   (s(i,k)*s(j,k)*s(k,m)) + 
     &  (s(i,m)*s(j,k))/((s(i,l) - s(k,m))*s(k,m)) + 
     &  ((2D0*s(j,k) - s(j,l))*s(k,m))/
     &   (s(i,l)*(s(i,l) - s(k,m))) + 
     &  ((s(i,m) + s(j,k))*s(j,m)*(s(i,m) + s(k,m)))/
     &   (s(j,k)*(s(i,m) - s(k,l))*s(k,l)) + 
     &  (s(i,l)*(-s(i,m)**2 + s(i,l)*s(k,m)))/
     &   (s(i,k)*(s(i,k) + s(i,m))*s(k,m))

      If (iTest.Eq.1) then
 
        Print*,'C1pppmm:',(s(i,l) - s(i,m))/s(i,k),
     &  (s(i,l)*(s(i,l) + 2D0*s(i,m)))/(s(i,k)*(s(i,k) + s(i,l))),
     &  s(i,l)**2/(s(i,k)*s(j,k)),
     &  ((s(i,l) + s(i,m))*(s(i,l) - 3D0*s(k,l)))/(s(i,k)*s(j,k)),
     &  (s(i,j)*(2D0*s(i,k) + s(i,l) + s(i,m)))/
     &   (s(j,k)*(s(i,m) - s(k,l))),
     &  ((s(i,k) + s(i,m) - 3D0*s(j,k) - s(k,l)))/s(i,l) + 
     &  (2D0*s(i,j) + s(i,k))/s(k,l),
     &  (s(i,j)**2 + 2D0*s(i,j)*s(i,k) + s(i,k)**2 + s(i,j)*s(j,m) + 
     &     s(i,k)*s(j,m) + s(j,m)**2)/(s(j,k)*s(k,l)),
     &  ((s(i,k) - s(i,m))**2*(s(i,k) + s(k,l)))/
     &   (s(i,k)*s(i,l)*s(j,k)),
     &  ((s(i,k) - s(i,m))**2*(s(i,k) + s(k,l)))/
     &   (s(i,l)*(s(i,k) + s(i,l))*s(j,k)), 
     &  ((s(i,k) - s(i,l))**2*(s(i,k) + s(k,l)))/
     &   (s(i,k)*(s(i,k) + s(i,m))*s(j,k)),
     &  ((s(i,k) - 2D0*s(i,l) - 3D0*s(i,m) - 3D0*s(k,m)))/s(j,k),
     &  (s(i,k)*s(j,k))/(s(i,l)*(s(i,l) - s(k,m))),
     &  (s(i,l)*s(i,m)**2)/(s(i,k)*s(j,k)*s(k,m)),
     &  ((s(i,k) + s(i,l))*s(j,m))/(s(i,l)*s(k,m)),
     &  (s(i,l)*s(i,m)*(s(i,k) - s(i,m) + s(k,l)))/
     &   (s(i,k)*s(j,k)*s(k,m)),
     &  (s(i,m)*s(j,k))/((s(i,l) - s(k,m))*s(k,m)),
     &  ((2D0*s(j,k) - s(j,l))*s(k,m))/
     &   (s(i,l)*(s(i,l) - s(k,m))),
     &  ((s(i,m) + s(j,k))*s(j,m)*(s(i,m) + s(k,m)))/
     &   (s(j,k)*(s(i,m) - s(k,l))*s(k,l)),
     &  (s(i,l)*(-s(i,m)**2 + s(i,l)*s(k,m)))/
     &   (s(i,k)*(s(i,k) + s(i,m))*s(k,m))
       Print*, 'i,j,k,l,m:',i,j,k,l,m
        Do i1=1,5
        Do i2=i1,5
          Print*, 's(',i1,i2,')=',s(i1,i2)
        End Do
        End Do
        Print*, 'p1:',p(1,1),p(2,1),p(3,1),p(4,1)
        Print*, 'p2:',p(1,2),p(2,2),p(3,2),p(4,2)
        Print*, 'p3:',p(1,3),p(2,3),p(3,3),p(4,3)
        Print*, 'p4:',p(1,4),p(2,4),p(3,4),p(4,4)
        Print*, 'p5:',p(1,5),p(2,5),p(3,5),p(4,5)
!       tmp = 2.d0*Dot(p(1,1),p(1,4))
        Print*, 'C1pppmm:',C1pppmm
        Print*, ' --- '
        iTest = 0
      End If 
        
      Return
      End       
                  
CsB &&&&----------------------------------------------------------------
      
      Function CF1mpppm(i,j,k,l,m)
      
      Implicit None
      Complex CF1mpppm
      Integer i,j,k,l,m
      Real*8 p,s
      Complex pc
      Common / sij / p(4,5),pc(4,5),s(5,5)
      
      CF1mpppm = (-s(i,l))/s(j,k) - (s(i,k)*s(i,l))/(s(i,j)*s(j,k)) + 
     &  (5.D-1*s(i,m)**2)/(s(i,j)*s(j,k)) - (s(i,l))/s(k,l)
      
      Return
      End 
      
CsB &&&&----------------------------------------------------------------
      
      Function CF3mpppm(i,j,k,l,m)
      
      Implicit None
      Complex CF3mpppm
      Integer i,j,k,l,m
      Real*8 p,s
      Complex pc
      Common / sij / p(4,5),pc(4,5),s(5,5)
      
      CF3mpppm = s(i,k)/s(j,k) + s(i,k)/s(j,l) + 
     &  (s(i,j)*s(i,k))/(s(i,l)*s(j,l)) - 
     &  (5.D-1*s(i,m)**2)/(s(i,l)*s(j,l)) + 
     &  (2D0*s(i,k)*(1D0 + s(j,l)/s(j,k) - (s(i,j))/s(l,m) - 
     &       (s(i,j)*s(i,m))/(s(j,k)*s(l,m))))/s(j,k) - 
     &  ((-2D0*s(i,j)*s(i,k) + s(i,m)**2))/(s(i,l)*s(l,m))
      
      Return
      End 
      
CsB &&&&----------------------------------------------------------------
      
      Function CF4ampppm(i,j,k,l,m)
      
      Implicit None
      Complex CF4ampppm
      Integer i,j,k,l,m
      Real*8 p,s
      Complex pc
      Common / sij / p(4,5),pc(4,5),s(5,5)
      
      CF4ampppm = (-s(i,k))/s(j,k) - (s(i,k))/s(j,l) + 
     &  (5.D-1*(-2D0*s(i,j)*s(i,k) + s(i,m)**2))/(s(i,l)*s(j,l))
      
      Return
      End 
      
CsB &&&&----------------------------------------------------------------
      
      Function CF4bmpppm(i,j,k,l,m)
      
      Implicit None
      Complex CF4bmpppm
      Integer i,j,k,l,m
      Real*8 p,s
      Complex pc
      Common / sij / p(4,5),pc(4,5),s(5,5)
      
      CF4bmpppm = s(i,j)/s(j,k) + (s(i,j) - s(i,k))/s(k,l) - 
     &  (5.D-1*s(i,k)**2)/(s(i,l)*s(k,l)) - 
     &  (5.D-1*(s(i,k) + s(i,m))**2)/(s(i,l)*s(k,l))
      
      Return
      End 
                        
CsB &&&&----------------------------------------------------------------
      
      Function CB1mpppm(i,j,k,l,m)
      
      Implicit None
      Complex CB1mpppm
      Integer i,j,k,l,m
      Real*8 p,s
      Complex pc
      Common / sij / p(4,5),pc(4,5),s(5,5)
      
      CB1mpppm = (-2D0*(s(i,k) - s(i,m))*(s(i,j) + s(j,k)))/s(j,l)**2 + 
     &  (3D0*s(i,j) + 3D0*s(i,m) + s(j,k))/s(j,l) + 
     &  (s(j,l)**2*(s(i,j) + s(j,l)))/(s(i,j)*s(j,k)**2) + 
     &  (2D0*s(i,j) + 4D0*s(i,m) + 4D0*s(j,l))/s(j,k) - 
     &  ((2D0*s(i,j) + 2D0*s(j,l) - s(k,m)))/
     &   (s(i,k) + s(j,k)) - 
     &  ((s(i,k) - s(j,l))**2*(2D0*s(i,j) + s(j,l)))/
     &   (s(j,l)**2*s(k,m)) - 
     &  (2D0*s(j,k)*(s(i,j) + s(j,k))*s(k,m))/(s(i,j)*s(j,l)**2) - 
     &  (2D0*s(j,l)*s(k,m))/(s(i,j)*(s(i,k) + s(j,k))) - 
     &  (3D0*s(j,l)*s(k,m)**2)/(s(i,j)*s(j,k)**2) - 
     &  (s(j,l)*s(k,m)**2)/(s(i,j)*(s(i,k) + s(j,k))**2) + 
     &  (s(j,l)*(s(i,j) + s(j,l))*(2D0*s(k,m) - 3D0*s(l,m)))/
     &   (s(i,j)*s(j,k)**2) - 
     &  (s(j,k)*(3D0*s(i,l) - s(j,k) - 2D0*s(l,m)))/
     &   (s(i,j)*s(j,l)) - ((2D0*s(j,l) - 2D0*s(l,m)))/
     &   (s(i,l) + s(j,l)) - ((s(k,l) - 2D0*s(l,m)))/s(i,j) - 
     &  ((4D0*s(i,j) - s(i,l))*s(i,l)**2)/(s(j,k)**2*s(l,m)) + 
     &  (s(i,l)*(2D0*s(i,j) + s(j,k) - 2D0*s(j,l)))/(s(j,k)*s(l,m)) - 
     &  (2D0*s(i,l)*(s(i,l) - s(j,l)))/(s(j,k)*s(l,m)) - 
     &  (4D0*s(i,l)*s(l,m))/(s(i,j)*s(j,k)) - 
     &  (s(j,l)*s(k,m)*s(l,m))/(s(i,j)*(s(i,l) + s(j,l))**2) + 
     &  (s(k,m)**2*(3D0*s(j,l) + s(k,m))*s(l,m))/
     &   (s(i,j)*s(j,k)**2*(s(i,k) + s(j,k))) + 
     &  (s(j,l)*s(k,m)*(-s(j,k) + 2D0*s(k,m))*s(l,m))/
     &   (s(i,j)*s(j,k)*(s(i,k) + s(j,k))**2) - 
     &  (s(l,m)**2)/(s(i,j)*s(j,l)) - 
     &  ((s(j,l) - s(k,m))*s(l,m)**2)/
     &   (s(i,j)*(s(i,l) + s(j,l))**2) + 
     &  (s(k,m)*s(l,m)**2)/(s(i,j)*s(j,l)*(s(i,l) + s(j,l))) - 
     &  (s(k,m)*(-s(j,k) + s(k,m))**2*s(l,m)**2)/
     &   (s(i,j)*s(j,k)**2*(s(i,k) + s(j,k))**2) - 
     &  (2D0*s(k,m)*(s(k,m) - 2D0*s(l,m))*(s(j,l) + s(l,m)))/
     &   (s(i,j)*s(j,k)*(s(i,k) + s(j,k))) - 
     &  (s(l,m)*(2D0*s(j,l) + s(l,m)))/(s(i,j)*(s(i,l) + s(j,l))) - 
     &  ((s(i,j) + s(j,l))*(8D0*s(i,l) - s(j,l) + 3D0*s(l,m)))/
     &   s(j,k)**2 - (s(j,l)*(s(i,l) - 2D0*s(j,l) + 4D0*s(l,m)))/
     &   (s(i,j)*s(j,k)) + (2D0*
     &     (s(i,j)*s(j,l) + s(i,j)*s(k,l) - s(j,l)*s(l,m))**3)/
     &   (s(i,j)*s(j,k)**3*(s(i,k) + s(j,k))*s(l,m))
      
      Return
      End 
      
CsB &&&&----------------------------------------------------------------
      
      Function CB2mpppm(i,j,k,l,m)
      
      Implicit None
      Complex CB2mpppm
      Integer i,j,k,l,m
      Real*8 p,s
      Complex pc
      Common / sij / p(4,5),pc(4,5),s(5,5)
      
Checked
      CB2mpppm = s(i,l)/s(j,l) + (s(i,k)**2*s(j,l))/(s(i,j)*s(j,k)**2) - 
     &  (s(j,l)*s(k,m))/s(j,k)**2 - 
     &  (s(i,k)**2*s(i,l)*s(j,l))/
     &   (s(i,j)*s(j,k)**2*(s(j,k) + s(k,m))) - 
     &  (s(j,m)*(s(i,k) + s(k,m))**2)/
     &   (s(j,k)**2*(s(j,k) + s(k,m))) - 
     &  (s(i,l)*(s(i,j)*s(j,l) + s(i,l)*s(j,m)))/
     &   (s(i,j)*s(j,l)*(s(j,l) + s(l,m)))
      
      Return
      End 

CsB &&&&----------------------------------------------------------------
      
      Function CB3mpppm(i,j,k,l,m)
      
      Implicit None
      Complex CB3mpppm
      Integer i,j,k,l,m
      Real*8 p,s
      Complex pc
      Common / sij / p(4,5),pc(4,5),s(5,5)
      
        CB3mpppm = (s(i,l) + s(j,l))/(s(i,j) + s(j,k)) - 
     &  (2D0*s(i,k)*(s(j,k) + s(j,l)))/(s(j,k)*(s(i,j) + s(j,k))) - 
     &  (1D0*s(i,j)*s(j,k)**2)/(s(i,l)*(s(i,k) + s(j,k))*s(k,l)) - 
     &  (1D0*(s(j,k) + s(j,l)))/s(k,l) - 
     &  (1D0*(s(j,k)**2 + s(i,j)*s(j,l) + 2D0*s(j,k)*s(j,l) + 
     &       s(j,l)**2))/(s(i,l)*s(k,l)) + (2D0*s(k,l))/s(j,k) - 
     &  (1D0*s(k,l))/(s(i,k) + s(j,k)) - 
     &  (1D0*s(i,j)**2*s(k,l)**2)/
     &   (s(i,l)*(s(i,j) + s(j,k))*s(j,l)**2) - 
     &  (1D0*s(j,k)*(s(j,k) + s(k,l)))/(s(i,k) + s(j,k))**2 - 
     &  (1D0*s(j,k)**2*(s(j,l) + s(k,l)))/
     &   (s(i,l)*(s(i,k) + s(j,k))**2) + 
     &  (2D0*s(j,k) + 2D0*s(k,l))/s(j,l) - 
     &  (1D0*(s(j,k) - 3D0*s(j,m) + 3D0*s(k,l)))/s(i,l) + 
     &  (s(j,k)*(-1D0*s(i,l)*s(j,k) + 2D0*s(i,j)*s(k,l)))/
     &   ((s(i,j) + s(j,k))*s(j,l)**2) + 
     &  (s(k,l)*(s(j,k) + 3D0*s(j,l) - 1D0*s(k,l) - 3D0*s(k,m)))/
     &   (s(i,l)*s(j,l)) - (1D0*(2D0*s(j,k) + s(j,l))*s(l,m))/
     &   (s(i,j) + s(j,k))**2 + 
     &  (s(j,k)**2*(s(j,l) + s(l,m)))/
     &   (s(i,l)*(s(i,j) + s(j,k))*s(j,l)) + 
     &  (s(j,k)**2*s(l,m)*(s(j,l) + s(l,m)))/
     &   (s(i,l)*(s(i,j) + s(j,k))**2*s(j,l))
      
      Return
      End 
      
CsB &&&&----------------------------------------------------------------

      Function CB4ampppm(i,j,k,l,m)
      
      Implicit None
      Complex CB4ampppm
      Integer i,j,k,l,m
      Real*8 p,s
      Complex pc
      Common / sij / p(4,5),pc(4,5),s(5,5)
      
      CB4ampppm = (-3D0*s(i,m)**2)/s(j,l)**2 - 
     &  (s(i,l)*s(j,k)**3)/(s(j,l)**2*(s(i,j) + s(j,l))**2) + 
     &  (s(j,k)**2*(3D0*s(i,k) + 2D0*s(j,k)))/
     &   (s(j,l)**2*(s(i,j) + s(j,l))) + 
     &  (2D0*s(i,m)**2*(s(i,m) - 3D0*s(k,l)))/s(j,l)**3 + 
     &  (s(j,k)**2*s(k,l))/(s(i,l)*(s(i,j) + s(j,l))**2) - 
     &  ((s(i,m) - s(j,l) + s(k,l)))/s(k,l) - 
     &  (3D0*s(k,l)*(s(i,k) + 2D0*s(k,l)))/(s(i,l)*s(j,l)) - 
     &  (s(k,l)**2*(2D0*s(i,k) + 3D0*s(k,l) - 3D0*s(k,m)))/
     &   (s(i,l)*s(j,l)**2) + 
     &  (s(i,m)**2*(s(i,m) - s(j,l))*(2D0*s(i,l) + s(j,l)))/
     &   (s(j,l)**3*s(k,m)) - 
     &  (s(j,k)*s(k,l))/((s(i,j) + s(j,l))*s(k,m)) + 
     &  (s(k,l)*(s(i,k) + s(k,l)))/(s(i,l)*s(k,m)) - 
     &  (2D0*s(j,k)**3*s(k,m))/(s(j,l)**3*(s(i,j) + s(j,l))) + 
     &  (2D0*s(k,l)**3*s(k,m))/(s(i,l)*s(j,l)**3) - 
     &  ((s(k,l) - s(k,m)))/(s(j,k) + s(k,m)) - 
     &  (s(k,l)*s(k,m))/(s(i,l)*(s(j,k) + s(k,m))) - 
     &  (s(j,l)*(s(j,k) + s(k,m)))/(s(i,l)*s(k,l)) + 
     &  (s(k,l) + s(k,m))/(s(i,j) + s(j,l)) - 
     &  (s(j,k)*(2D0*s(j,k) + s(j,l))*(s(k,l) + s(k,m)))/
     &   (s(j,l)*(s(i,j) + s(j,l))**2) - 
     &  ((s(j,l) + s(k,l) + s(k,m)))/s(i,l)
      
      Return
      End 
                        
CsB &&&&----------------------------------------------------------------
      
      Function CB4bmpppm(i,j,k,l,m)
      
      Implicit None
      Complex CB4bmpppm
      Integer i,j,k,l,m
      Real*8 p,s
      Complex pc
      Common / sij / p(4,5),pc(4,5),s(5,5)
      
      CB4bmpppm = (s(j,k)*s(j,l))/(s(i,l)*(s(i,j) + s(j,l))) + 
     &  (s(j,k)*(s(j,k) + s(j,l)))/(s(i,j) + s(j,l))**2 + 
     &  ((s(j,k) - s(j,l))*(s(j,k) + s(j,l)))/
     &   (s(j,l)*(s(i,j) + s(j,l))) + 
     &  (2D0*s(i,l) + 3D0*s(j,l))/s(k,l) + 
     &  (4D0*s(k,l)*(s(j,k) + s(k,l)))/s(j,l)**2 + 
     &  (3D0*s(j,k) + 6D0*s(j,l) + 9D0*s(k,l))/s(j,l) - 
     &  (s(j,k)**2*s(j,l))/((s(i,j) + s(j,l))**2*s(k,m)) - 
     &  (s(j,k)**2*s(j,l))/(s(i,l)*(s(i,j) + s(j,l))*s(k,m)) - 
     &  (s(i,m)**2*(2D0*s(i,l) + s(j,l)))/(s(j,l)**2*s(k,m)) - 
     &  (2D0*s(k,l)**2*s(k,m))/(s(i,l)*s(j,l)**2) + 
     &  ((s(i,l) + s(j,l))*(s(k,l) - 2D0*s(k,m)))/
     &   (s(k,l)*(s(j,k) + s(k,m))) + 
     &  (s(j,l)*s(k,m))/(s(i,l)*(s(j,k) + s(k,m))) + 
     &  (s(j,l)*(s(j,l) + s(k,m)))/(s(i,l)*s(k,l)) - 
     &  ((3D0*s(j,k) - 2D0*s(j,l) - 3D0*s(k,l) + 3D0*s(k,m)))/
     &   s(i,l) - (s(k,l)*(2D0*s(j,k) - s(k,l) + 5D0*s(k,m)))/
     &   (s(i,l)*s(j,l)) + (s(i,l)*s(j,k) + s(j,k)*s(j,l) + 
     &      s(j,l)*s(k,m))**2/(s(i,l)*s(k,l)**2*(s(j,k) + s(k,m)))
      
      Return
      End 
                        
CsB &&&&----------------------------------------------------------------
      
      Function C1mpppm(i,j,k,l,m)
      
      Implicit None
      Complex C1mpppm
      Integer i,j,k,l,m
      Real*8 p,s
      Complex pc
      Common / sij / p(4,5),pc(4,5),s(5,5)
      
      C1mpppm = (-(2D0*s(j,l) - 2D0*s(j,m)))/s(j,k) + 
     &  (s(j,l)**2*(s(j,l) - 3D0*s(k,m)))/
     &   (s(i,j)*s(j,k)*(s(i,j) + s(j,k))) - 
     &  (s(i,j)*(2D0*s(i,j) + 6D0*s(j,l) - 3D0*s(k,m)))/s(j,k)**2 + 
     &  (s(i,k)*(s(i,k) - 2D0*s(j,l)))/(s(j,l)*s(k,m)) + 
     &  (s(j,k)**2*s(k,m))/(s(i,j)*s(j,l)**2) - 
     &  ((s(j,k) + s(j,l))*s(k,m))/(s(j,k)*(s(i,j) + s(j,l))) - 
     &  (s(j,k)*(s(j,k) + s(j,l))*s(k,m))/
     &   (s(j,l)**2*(s(i,j) + s(j,l))) + 
     &  ((s(j,k) - s(k,m))**2*s(k,m))/
     &   (s(j,k)**2*(s(i,k) + s(j,k))) + 
     &  (2D0*s(i,k) - 2D0*s(j,k) + s(k,m))/s(j,l) + 
     &  (s(j,k) + s(j,l) + s(k,m))/(s(i,j) + s(j,k)) + 
     &  (s(j,l)*(s(j,l) + 2D0*s(k,m)))/(s(j,k)*(s(i,j) + s(j,k))) - 
     &  (s(j,l)**2*(2D0*s(j,l) - 3D0*s(k,m) - 3D0*s(l,m)))/
     &   (s(i,j)*s(j,k)**2) + 
     &  ((s(j,k) + s(j,l))*(s(j,k) - s(j,l) - 2D0*s(l,m)))/
     &   (s(i,j)*(s(i,j) + s(j,l))) - 
     &  (s(j,k)*(2D0*s(j,k) - 3D0*s(k,m) - 2D0*s(l,m)))/
     &   (s(i,j)*s(j,l)) - (3D0*s(j,l)*
     &     (2D0*s(j,l) - 2D0*s(k,m) - s(l,m)))/s(j,k)**2 + 
     &  (s(j,l)*s(k,m)**2)/(s(i,j)*s(l,m)**2) + 
     &  (s(j,l)*s(k,m)**2)/((s(i,k) + s(j,k))*s(l,m)**2) - 
     &  (s(i,l)*(s(i,l) - s(j,k))**2)/(s(j,k)**2*s(l,m)) - 
     &  (s(j,l)*(s(j,k) - s(k,m)))/(s(i,j)*s(l,m)) + 
     &  (s(j,l)*s(k,m))/((s(i,k) + s(j,k))*s(l,m)) - 
     &  (2D0*s(j,l)*s(k,m)**2)/(s(j,k)*(s(i,k) + s(j,k))*s(l,m)) + 
     &  (s(j,l)*(s(j,k) + s(k,m))**2)/
     &   (s(i,j)*(s(i,j) + s(j,k))*s(l,m)) + 
     &  (s(i,j)*s(l,m))/s(j,k)**2 + 
     &  (s(j,l)**3*s(l,m))/(s(i,j)*s(j,k)**2*(s(i,j) + s(j,k))) + 
     &  ((3D0*s(j,l) - s(l,m))*s(l,m))/(s(i,j)*s(k,m)) - 
     &  ((3D0*s(j,k) - s(j,l) + s(l,m)))/s(i,j) - 
     &  ((s(j,l) - s(l,m))*s(l,m))/(s(k,m)*(s(k,l) + s(l,m))) - 
     &  (s(j,l)*s(l,m)**2)/(s(i,j)*s(k,m)*(s(k,l) + s(l,m))) - 
     &  (s(j,l)*(2D0*s(k,m)**2 + 3D0*s(j,l)*s(l,m) + 
     &       3D0*s(k,m)*s(l,m)))/(s(i,j)*s(j,k)*s(l,m))
      
      Return
      End 
      
CsB --------------------------------------------------------------------
C     Supporting routines for gg -> AAg 
CsB --------------------------------------------------------------------

      Function Dot(p,q)

      Implicit None
      Real*8 p,q,Dot
      Dimension p(4),q(4)
!     Common / Test / iTest

      Dot = - p(1)*q(1) - p(2)*q(2) - p(3)*q(3) + p(4)*q(4)

!     If (iTest.Eq.1) then
!       Print*, p
!       Print*, q
!       Print*, -p(3)*q(3),p(4)*q(4)
!     End If

      Return
      End
      
CsB --------------------------------------------------------------------

      Function CDot(p,q)

      Implicit None
      Complex p,q,CDot

      Dimension p(4),q(4)

      CDot = - p(1)*q(1) - p(2)*q(2) - p(3)*q(3) + p(4)*q(4)

!     Print*, 'p =', p(1),p(2),p(3),p(4)
!     Print*, 'q =', q(1),q(2),q(3),q(4)
!     Print*, CDot
!     Print*, '---'

      Return
      End
      
CsB --------------------------------------------------------------------

      Function epsilon(i,j,k,l)
CsB   The fully antisymmetric tensor.  

      Implicit None
      Real*8 epsilon,eps(4,4,4,4)
      Integer i,j,k,l
!     Dimension eps(4,4,4,4)

      eps(1,2,3,4) =  1.d0
      eps(1,2,4,3) = -1.d0
      eps(1,3,2,4) = -1.d0
      eps(1,3,4,2) =  1.d0
      eps(1,4,2,3) =  1.d0
      eps(1,4,3,2) = -1.d0
      eps(2,1,3,4) = -1.d0
      eps(2,1,4,3) =  1.d0
      eps(2,3,1,4) =  1.d0
      eps(2,3,4,1) = -1.d0
      eps(2,4,1,3) = -1.d0
      eps(2,4,3,1) =  1.d0
      eps(3,1,2,4) =  1.d0
      eps(3,1,4,2) = -1.d0
      eps(3,2,1,4) = -1.d0
      eps(3,2,4,1) =  1.d0
      eps(3,4,1,2) =  1.d0
      eps(3,4,2,1) = -1.d0
      eps(4,1,2,3) = -1.d0
      eps(4,1,3,2) =  1.d0
      eps(4,2,1,3) =  1.d0
      eps(4,2,3,1) = -1.d0
      eps(4,3,1,2) = -1.d0
      eps(4,3,2,1) =  1.d0

      epsilon = 0.d0
      If (i.Eq.j .or. i.Eq.k .or. i.Eq.l .or. j.Eq.k .or. j.Eq.l 
     ..or.k.Eq.l) Return 

C     Note that the definition of epsilon also dependens on the fact
C     that momenta are defined as (px,py,pz,E) in ResBos.  Were they 
C     defined as (E,px,py,pz), epsilon would be -epsilon!
      epsilon = eps(i,j,k,l)
      
      Return
      End

CsB --------------------------------------------------------------------
C     Pentagon routines for gg -> AAg
CsB --------------------------------------------------------------------
CsB   These were cut out of pentagon.f

      subroutine ConstInit1(iProc)
C     Set up the initial parameters for GG->AAG via box (iproc=1) or GQ
c     ->AAQ via box (iproc=2)

      implicit NONE
      include 'pncommon.f'
      integer iProc

      zero=1d-14
      Pi= 3.141592653589793238d0
      cg=1./16./pi/pi

      ie=(0.,1.)     !imaginary unit

      nleg=5  !number of outer legs

      MomFlag=1   !read particle momenta
C     MomFlag=2   !read kinematic invariants

      Nc=3       !Number of colors

      if (iproc.eq.1) then      !gg->aag
        inout(1)=-1             !inout(j)=1 if the j-th particle is ingoing
        inout(2)=-1             !and -1 if it is outgoing, according to the  
        inout(3)=-1             !convention assumed in the matrix element
        inout(4)=1
        inout(5)=1

cpn2006 The overall normalization and color prefactor for GG->AAG cross section
c       corresponding to the overall normalization of the 5-leg box
c       amplitude in my Pentagon notes; color-averaging factor for 
c       2 initial gluons is included
        const=1253456.629       !4*(4*pi)**5
        const=const*Nc/(Nc**2.-1)
      elseif (iproc.eq.2) then  !qg->aaq
        inout(1)=-1              !inout(j)=1 if the j-th particle is ingoing
        inout(2)=1             !and -1 if it is outgoing, according to the
        inout(3)=-1              !convention assumed in the matrix element
        inout(4)=-1
        inout(5)=1

cpn2006 The overall normalization and color prefactor for QG->AAQ cross section
c       corresponding to the overall normalization of the box matrix
c       element M_5^3 in Bern-Dixon-Kosower paper, NPB 437, 259 (1995) 
c       the color-averaging factor for the initial-state
c       quark and gluon is included
        const=1253456.629       !4*(4*pi)**5

        const=const/Nc
      else
        print *,'ConstInit1: no such iproc=',iproc
        call ExitGracefully
      endif                     !iproc
      
      const=const/4./2.         !average over the spin of the incoming
                             !gluons and account for 2 identical particles
                                !in the final state
      
      End                       ! ConstInit1


C----------------------------------------------------------------------
      function fLi2(x)
C----------------------------------------------------------------------
C     Dilogarithm function for x+-i0. The precision is about 12 figures. 
C     Translation of Carl Schmidt's C code. The imaginary part for real x>1
C     is generated according to x -> x+i0

      implicit NONE
      include 'pncommon.f'

      double complex fLi2,Out,x,fLi2_2
      double precision czero,xr
      Data epsil / 1.e-10 / 
      external fLi2_2 
  
      czero=1d1*epsil
      pi=3.14159265358979d0

C     Check if the imaginary part of x is (almost) infinitesimally small
      if (dabs(dimag(x)-epsil).gt.czero) then
        print *,'fLi2: too large imaginary part of x=',dimag(x)
        STOP
      endif

      xr=dreal(x)

      if (xr.ge.0.and.xr.le.1.0) then
        Out=fLi2_2(x)
      elseif (xr.lt.0) then
        Out=fLi2_2(1/(1-x))+0.5*(dlog(1-xr))**2-
     &      dlog(-xr)*dlog(1-xr)-pi*pi/6.
      elseif (xr.gt.1.0) then
        Out=-fLi2_2(1/x)+pi*pi/3.-0.5*(dlog(xr))**2

        if (dimag(x).ge.0) then
          Out=Out+ie*pi*dlog(xr)
        else
          Out=Out-ie*pi*dlog(xr)
        endif

      endif

      fLi2=Out

      End ! fLi2

C----------------------------------------------------------------------------

      double complex function fLi2_2(x)
      implicit NONE
      include 'pncommon.f'
     
      double complex x,Out
      double precision xr,r,xm,L

      integer NCUT, i,N
      NCUT =27

      L=0.0
      i=0
      r=1.0

      xr=dreal(x)

      if (xr.ge.0.and.xr.le.0.5) then

       do while(i.le.NCUT)
            i=i+1
        r=r*xr

        L=L+r/i/i
       enddo
        Out=L
      elseif(xr.gt.0.5.and.xr.le.1.0) then
       xm = 1.0-xr
       L = pi*pi/6.0 - dlog(xr)*dlog(xm)

       do while(i.le.NCUT)
        i=i+1
        r=r*xm

        L=L-r/i/i
       enddo
      endif 
      
      Out=L  

      fLi2_2=Out      
      return
      end! fLi2_2 
cpn ------------------------------------------------------------------

      function rfLi2x(x)
cpn   Returns the real part of Li^2 divided by x
      implicit none
      double complex fLi2
      double precision x, rfLi2x
      rfLi2x = dreal(fLi2(dcmplx(x)))/x
      return
      end ! rfLi2x

cpn---------------------------------------------- Higher-order polylogs

      function fLi3(x)
cpn   Polylogarithm of the third order; for simplicity, implemented only
cpn   where it is purely real (i.e., for x<1)
cpn   Inputs: x  -- argument 
cpn           rlierr - relative error of integration (passed through the
cpn           common block lierr)
      double precision  rfLi2x, fLi3      
      double precision x, adzint, aerr, rerr, errest, tmpr
      integer iacta, iactb
      external rfLi2x

      double precision rlierr
      common/lierr/rlierr

      if (x.gt.1d0) then
        print *, 'STOP: Li(3,x) is implemented only for x<1'
        stop
      endif 
     
      aerr = 1e-16
      rerr = rlierr

cpn   Calculate the real part of Li[3,x]; since adzint cares about the
cpn   order of limits of integration, switch limits if x < 0
      tmpr = 0d0
      if (x.gt.0) then 
        iacta = 2
        iactb = 0
        tmpr = adzint(rfLi2x,0d0, x, aerr, rerr, errest, ier, iacta,
     >    iactb)
      else
        iacta = 0
        iactb = 2
        tmpr = -adzint(rfLi2x,x, 0d0, aerr, rerr, errest, ier, iacta,
     >    iactb)
      endif                     !x

c$$$      tmpi = 0d0
c$$$      if (x.gt.1d0) 
c$$$     >  tmpr = adzint(iLi2,0, x, aerr, rerr, errest, ier, iacta,
c$$$      iactb)
c$$$      
c$$$      fLi3 = tmpr + ie*tmpi

      fLi3 = tmpr

      return 
      end !fLi3

cpn ------------------------------------------------------------------
      function rfLi3x(x)
cpn   Returns the real part of Li^3 divided by x
      implicit none
      double precision x, fLi3, rfLi3x
!      rfLi3x = fLi3(dcmplx(x))/x   The argument of fLi3 should be a real number.  Feb 20225 Ibr
       rfLi3x = fLi3(x)/x
      return
      end ! rfLi3x

cpn---------------------------------------------- Higher-order polylogs

      function fLi4(x)
cpn   Polylogarithm of the fourth order; for simplicity, implemented
cpn   only where it is purely real (i.e., for x<1)
cpn   Inputs: x  -- argument 
cpn           rlierr - relative error of integration (passed through the
cpn           common block lierr)
      double precision  rfLi3x, fLi4      
      double precision x, adz2nt, aerr, rerr, errest, tmpr
      integer iacta, iactb
      external rfLi3x

      double precision rlierr
      common/lierr/rlierr

      if (x.gt.1d0) then
        print *, 'STOP: Li(3,x) is implemented only for x<1'
        stop
      endif 
     
      aerr = 1e-16
      rerr = rlierr

cpn   Calculate the real part of Li[4,x]; since adzint cares about the
cpn   order of limits of integration, switch limits if x<0
cpn   Use the second copy of adzint (adz2nt) to prevent interference
cpn   with integration in fLi3
      tmpr  =0d0
      if (x.gt.0) then
        iacta = 2
        iactb = 0
        tmpr = adz2nt(rfLi3x,0d0, x, aerr, rerr, errest, ier, iacta,
     >    iactb)
      else
        iacta = 0
        iactb = 2
        tmpr = -adz2nt(rfLi3x,x,0d0, aerr, rerr, errest, ier, iacta,
     >    iactb)
      endif                     !x

c$$$      tmpi = 0d0
c$$$      if (x.gt.1d0) 
c$$$     >  tmpr = adzint(iLi2,0, x, aerr, rerr, errest, ier, iacta,
c$$$      iactb)
c$$$      
c$$$      fLi4 = tmpr + ie*tmpi

      fLi4 = tmpr

      return 
      end !fLi4

cpn-------------------------------------------------------------------------
C                                                          =-=-= Adzint
      FUNCTION ADZINT (F, A, B, AERR, RERR, ERREST, IER, IACTA, IACTB)
C                                                   -=-=- adzint

C===========================================================================
C GroupName: Adzint
C Description: adaptive integration
C ListOfFiles: adzint adzcal adzspl intusz sglint totalz
C===========================================================================
C                                  Authors: Wu-Ki Tung and John C. Collins
C #Header: /Net/cteq06/users/wkt/1hep/1utl/RCS/Adzint.f,v 1.1 97/12/21 21:19:04 wkt Exp $
C #Log:     Adzint.f,v $
c Revision 1.1  97/12/21  21:19:04  wkt
c Initial revision
c 

C     FUNCTION   ADZINT (F, A, B, AERR, RERR, ERREST, IER, IACTA, IACTB)
C     SUBROUTINE ADZSPL (F, I, IER)
C     SUBROUTINE ADZCAL (F,I)
C     SUBROUTINE SGLINT (IACT, F1, F2, F3, DX, FINT, ESTER)
C     SUBROUTINE TOTALZ
C     FUNCTION   INTUSZ (X, FX)
C
C     COMMON / ADZWRK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES, 
C    > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), NUMINT,
C    > ICTA, ICTB, FA, FB, IB
C                        ------------------------

C     Adaptive integration routine which allows the integrand to be 
C     indeterminant at the lower and/or the upper ends of integration. 

C     Can self-adjust to any integrable singularity at the ends and compute 
C     the closest approximant, hence achieve the required accuracy efficiently
C     (provided the switch(s) IACTA (IACTB) are set to 2).
 
C     Input switches for end-treatment:
C        IACTA = 0 :   Use closed lower-end algorithm 
C                1 :   Open lower-end -- use open quadratic approximant
C                2 :   Open lower-end -- use adaptive singular approximant

C        IACTB = 0, 1, 2   (same as above, for the upper end)
 
C                Integral of F(X) from A to B, with error
C                less than ABS(AERR) + ABS(RERR*INTEGRAL)
C                Best estimate of error returned in ERREST.
CError code is IER: 0 :  o.k.
C                1 :  maximum calls to function reached before the 
C                     error criteria are met;
C                2 :  IACTA out of range, set to 1;
C                3 :  IACTB out of range, set to 1.
C                4 :  Error on Limits : B < A ; zero result returned.
C                5 :  Range of integration DX zero or close to roundoff
C                     returns DX * F(A+DX/2)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      EXTERNAL F
      PARAMETER (MAXINT = 1000)
C
C                   Work space:
      COMMON / ADZWRK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES, 
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
     > ICTA, ICTB, NUMINT, IB

      SAVE / ADZWRK /
      DATA SMLL, Sml / 1E-20, 1E-12 /
     
      IER = 0
      IF (AERR.LE.SMLL .AND. RERR.LE.SMLL)
     1 STOP 'Both Aerr and Rerr are zero in ADZINT!'
        
      IF (IACTA.LT.0 .OR. IACTA.GT.2) THEN
        PRINT '(A, I4/ A)', ' Illegal value of IACT in ADZINT call', 
     >  'IACTA =', IACTA, ' IACTA set for regular open-end option.'
        IACTA = 1
        IER = 2
      ENDIF 
      IF (IACTB.LT.0 .OR. IACTB.GT.2) THEN
        PRINT '(A, I4/ A)', ' Illegal value of IACT in ADZINT call', 
     >  'IACTB =', IACTB, ' IACTB set for regular open-end option.'
        IACTB = 1
        IER = 3
      ENDIF
      ICTA = IACTA
      ICTB = IACTB
 
      DDX = B - A
      If (DDX .Le. 0D0) Then
        AdzInt = 0D0
        Ier = 4
        If (DDX .Lt. 0D0) 
     >     Print '(/A/)', 'B < A in AdzInt; check limits!!'
        Return
      ElseIf (DDX .Le. Sml) Then
        AdzInt = F(A + DDX/2) * DDX
        Ier = 5
        Return
      EndIf

      NUMINT = 3
      DX = DDX/ NUMINT
      DO 10  I = 1, NUMINT
          IF (I .EQ. 1)  THEN
             U(1) = A 
             IF (IACTA .EQ. 0) THEN
               FU(1) = F(U(1))
             ELSE 
C                                   For the indeterminant end point, use the
C                                   midpoint as a substitue for the endpoint.
               FA = F(A+DX/2.)
             ENDIF
          ELSE
              U(I) = V(I-1)
              FU(I) = FV(I-1)
          ENDIF

          IF (I .EQ. NUMINT) THEN
             V(I) = B
             IF (IACTB .EQ. 0) THEN
               FV(I) = F(V(I))
             ELSE
               IB = I
               FB = F(B-DX/2.)
             ENDIF
          ELSE
              V(I) = A + DX * I
              FV(I) = F(V(I))
          ENDIF
          CALL ADZCAL(F,I)
   10     CONTINUE
       CALL TOTALZ
C                                                   Adaptive procedure:
   30     TARGET = ABS(AERR) + ABS(RERR * RES)
          IF (ERS .GT. TARGET)  THEN
              NUMOLD = NUMINT
              DO 40, I = 1, NUMINT
                  IF (ERR(I)*NUMOLD .GT. TARGET) CALL ADZSPL(F,I,IER)
   40         CONTINUE
              IF (IER.EQ.0 .AND. NUMINT.NE.NUMOLD)  GOTO 30
              ENDIF
      ADZINT = RES
      ERREST = ERS
      RETURN
C                        ****************************
      END

      SUBROUTINE ADZCAL (F,I)
C                                                   -=-=- adzcal
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D1 = 1.0, D2 = 2.0, HUGE = 1.E15)
C                        Fill in details of interval I given endpoints
      EXTERNAL F
      PARAMETER (MAXINT = 1000)
      COMMON / ADZWRK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES, 
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
     > ICTA, ICTB, NUMINT, IB
 
      SAVE / ADZWRK /

      DX =  V(I) - U(I)
      W  = (U(I) + V(I)) / 2.
     
      IF (I .EQ. 1 .AND. ICTA .GT. 0) THEN
C                                                                 Open LEFT end
        FW(I) = FA
        FA = F (U(I) + DX / 4.)

        CALL SGLINT (ICTA, FA, FW(I), FV(I), DX, TEM, ER)
      ELSEIF (I .EQ. IB .AND. ICTB .GT. 0) THEN
C                                                                open RIGHT end
        FW(I) = FB
        FB = F (V(I) - DX / 4.)
        CALL SGLINT (ICTB, FB, FW(I), FU(I), DX, TEM, ER)
      ELSE
C                                                                   Closed endS
        FW(I) = F(W)
        TEM = DX * (FU(I) + 4. * FW(I) + FV(I)) / 6.
C                                       Preliminary error Simpson - trapezoidal:
        ER  = DX * (FU(I) - 2. * FW(I) + FV(I)) / 12.
      ENDIF
 
      RESULT(I) = TEM         
      ERR   (I) = ABS (ER)
 
      RETURN
C                        ****************************
      END

      SUBROUTINE ADZSPL (F, I, IER)
C                                                   -=-=- adzspl
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C                                                      Split interval I
C                                                   And update RESULT & ERR
      EXTERNAL F
      PARAMETER (MAXINT = 1000)
      COMMON / ADZWRK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES, 
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
     > ICTA, ICTB, NUMINT, IB

      SAVE / ADZWRK /
      DATA TINY / 1.D-20 /
     
      IF (NUMINT .GE. MAXINT)  THEN
          IER = 1
          RETURN
          ENDIF
      NUMINT = NUMINT + 1
C                                                         New interval NUMINT
      IF (I .EQ. IB) IB = NUMINT
      U(NUMINT) = (U(I) + V(I)) / 2.
      V(NUMINT) = V(I)
 
      FU(NUMINT) = FW(I)
      FV(NUMINT) = FV(I)
C                                                             New interval I
       V(I) =  U(NUMINT)
      FV(I) = FU(NUMINT)
C                                                    Save old Result and Error
      OLDRES = RESULT(I)
      OLDERR = ERR(I)
     
      CALL ADZCAL (F, I)
      CALL ADZCAL (F, NUMINT)
C                                                               Update result
      DELRES = RESULT(I) + RESULT(NUMINT) - OLDRES
      RES = RES + DELRES
C                                  Good error estimate based on Simpson formula
      GODERR = ABS(DELRES) 
C                                                             Update new global 
      ERS = ERS + GODERR - OLDERR
C                                  Improve local error estimates proportionally
      SUMERR = ERR(I) + ERR(NUMINT)
      IF (SUMERR .GT. TINY) THEN
         FAC = GODERR / SUMERR 
      ELSE
         FAC = 1.
      ENDIF
      
      ERR(I)      = ERR(I) * FAC
      ERR(NUMINT) = ERR(NUMINT) * FAC
 
      RETURN
C                        ****************************
      END
 
      FUNCTION INTUSZ (X, FX)
C                                                   -=-=- intusz

C                    Return number of distinct points used in AdzInt
C                  Also returns the x- and F(x)-values (ie. the integrand)
C                       at these points (perhaps for plotting purposes).

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (MAXINT = 1000)
      COMMON / ADZWRK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES, 
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
     > ICTA, ICTB, NUMINT, IB

      Dimension X(MAXINT), FX(MAXINT)

      SAVE / ADZWRK /

       Call Sort2 (NumInt, U, Fu)
       Do 10 Ix = 1, NUMINT
         X(Ix)  =  U(Ix)
         FX(Ix) = Fu(Ix)
   10  Continue

       Call Sort2 (NumInt, V, Fv)
       X (NumInt + 1) = V (NumInt)
       FX(NumInt + 1) = Fv(NumInt)

      INTUSZ = NUMINT + 1

      RETURN
C                        ****************************
      END

      SUBROUTINE SGLINT (IACT, F1, F2, F3, DX, FINT, ESTER)
C                                                   -=-=- sglint

C     Calculate end-interval using open-end algorithm based on function values
C     at three points at (1/4, 1/2, 1)DX from the indeterminant endpoint (0).

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)

      DATA HUGE / 1.E20 /
C                                                         Use quadratic formula
      TEM = DX * (4.*F1 + 3.*F2 + 2.*F3) / 9.
C                 Error est based on Diff between quadratic and linear integrals
      ER  = DX * (4.*F1 - 6.*F2 + 2.*F3) / 9.

C                          Invoke adaptive singular parametrization if IACT = 2
C                      Algorithm is based on the formula F(x) = AA + BB * x **CC
C                 where AA, BB & CC are determined from F(Dx/4), F(Dx/2) & F(Dx)

      IF (IACT .EQ. 2) THEN
          T1 = F2 - F1
          T2 = F3 - F2
          IF (T1*T2 .LE. 0.) GOTO 7
          T3  = T2 - T1
          IF (ABS(T3)*HUGE .LT. T1**2) GOTO 7
          CC  = LOG (T2/T1) / LOG(D2)
          IF (CC .LE. -D1)  GOTO 7
          BB  = T1**2 / T3
          AA  = (F1*F3 - F2**2) / T3
C                                          Estimated integral based on A+Bx**C
          TMP = DX * (AA + BB* 4.**CC / (CC + 1.))
C                                       Error estimate based on the difference
          ER = TEM - TMP
C                                              Use the improved integral value
          TEM= TMP 
      ENDIF

    7 FINT = TEM
      ESTER= ER
      RETURN
C                        ****************************
      END
     
      SUBROUTINE TOTALZ
C                                                   -=-=- totalz
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (MAXINT = 1000)
      COMMON / ADZWRK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES, 
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
     > ICTA, ICTB, NUMINT, IB

      SAVE / ADZWRK /
      RES = 0.
      ERS = 0.
      DO 10  I = 1, NUMINT
          RES = RES + RESULT(I)
          ERS = ERS + ERR(I)
   10     CONTINUE
C                        ****************************
      END

C                                                          =-=-= Adz2nt
      FUNCTION ADZ2NT (F, A, B, AERR, RERR, ERREST, IER, IACTA, IACTB)
C                                                   -=-=- adz2nt
 
C===========================================================================
C GroupName: Adz2nt
C Description: second copy of adzint
C ListOfFiles: adz2nt adz2pl adz2al int2sz sgl2nt tot2lz
C=========================================================================== 
C #Header: /Net/cteq06/users/wkt/1hep/1utl/RCS/Adz2nt.f,v 1.1 97/12/21 21:19:00 wkt Exp $
C #Log:     Adz2nt.f,v $
c Revision 1.1  97/12/21  21:19:00  wkt
c Initial revision
c 

C List of GLOBAL Symbols

C     FUNCTION   ADZ2NT (F, A, B, AERR, RERR, ERREST, IER, IACTA, IACTB)
C     SUBROUTINE ADZ2PL (F, I, IER)
C     SUBROUTINE ADZ2AL (F,I)
C     SUBROUTINE SGL2NT (IACT, F1, F2, F3, DX, FINT, ESTER)
C     SUBROUTINE TOT2LZ
C     FUNCTION   INT2SZ (X, FX)
C
C     COMMON / ADZ2RK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES, 
C    > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
C    > ICTA, ICTB, NUMINT, IB
C                   ------------------------------------

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      EXTERNAL F
      PARAMETER (MAXINT = 1000)
C
C                   Work space:
      COMMON / ADZ2RK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES, 
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
     > ICTA, ICTB, NUMINT, IB

      SAVE / ADZ2RK /
      DATA SMLL / 1E-20 /
     
      IER = 0
      IF (AERR.LE.SMLL .AND. RERR.LE.SMLL)
     1 STOP 'Both Aerr and Rerr are zero in ADZ2NT!'
        
      IF (IACTA.LT.0 .OR. IACTA.GT.2) THEN
        PRINT '(A, I4/ A)', ' Illegal value of IACT in ADZ2NT call', 
     >  'IACTA =', IACTA, ' IACTA set for regular open-end option.'
        IACTA = 1
        IER = 2
      ENDIF 
      IF (IACTB.LT.0 .OR. IACTB.GT.2) THEN
        PRINT '(A, I4/ A)', ' Illegal value of IACT in ADZ2NT call', 
     >  'IACTB =', IACTB, ' IACTB set for regular open-end option.'
        IACTB = 1
        IER = 3
      ENDIF
      ICTA = IACTA
      ICTB = IACTB
 
      NUMINT = 3
      DX = (B-A)/ NUMINT
      DO 10  I = 1, NUMINT
          IF (I .EQ. 1)  THEN
             U(1) = A 
             IF (IACTA .EQ. 0) THEN
               FU(1) = F(U(1))
             ELSE 
C                                   For the indeterminant end point, use the
C                                   midpoint as a substitue for the endpoint.
               FA = F(A+DX/2.)
             ENDIF
          ELSE
              U(I) = V(I-1)
              FU(I) = FV(I-1)
          ENDIF

          IF (I .EQ. NUMINT) THEN
             V(I) = B
             IF (IACTB .EQ. 0) THEN
               FV(I) = F(V(I))
             ELSE
               IB = I
               FB = F(B-DX/2.)
             ENDIF
          ELSE
              V(I) = A + DX * I
              FV(I) = F(V(I))
          ENDIF
          CALL ADZ2AL(F,I)
   10     CONTINUE
       CALL TOT2LZ
C                                                   Adaptive procedure:
   30     TARGET = ABS(AERR) + ABS(RERR * RES)
          IF (ERS .GT. TARGET)  THEN
              NUMOLD = NUMINT
              DO 40, I = 1, NUMINT
                  IF (ERR(I)*NUMOLD .GT. TARGET) CALL ADZ2PL(F,I,IER)
   40             CONTINUE
              IF (IER.EQ.0 .AND. NUMINT.NE.NUMOLD)  GOTO 30
              ENDIF
      ADZ2NT = RES
      ERREST = ERS
      RETURN
C                        ****************************
      END

      SUBROUTINE ADZ2PL (F, I, IER)
C                                                   -=-=- adz2pl
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C                                                      Split interval I
C                                                   And update RESULT & ERR
      EXTERNAL F
      PARAMETER (MAXINT = 1000)
      COMMON / ADZ2RK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES, 
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
     > ICTA, ICTB, NUMINT, IB

      SAVE / ADZ2RK /
      DATA TINY / 1.D-20 /
     
      IF (NUMINT .GE. MAXINT)  THEN
          IER = 1
          RETURN
          ENDIF
      NUMINT = NUMINT + 1
C                                                         New interval NUMINT
      IF (I .EQ. IB) IB = NUMINT
      U(NUMINT) = (U(I) + V(I)) / 2.
      V(NUMINT) = V(I)
 
      FU(NUMINT) = FW(I)
      FV(NUMINT) = FV(I)
C                                                             New interval I
       V(I) =  U(NUMINT)
      FV(I) = FU(NUMINT)
C                                                    Save old Result and Error
      OLDRES = RESULT(I)
      OLDERR = ERR(I)
     
      CALL ADZ2AL (F, I)
      CALL ADZ2AL (F, NUMINT)
C                                                               Update result
      DELRES = RESULT(I) + RESULT(NUMINT) - OLDRES
      RES = RES + DELRES
C                                  Good error estimate based on Simpson formula
      GODERR = ABS(DELRES) 
C                                                             Update new global 
      ERS = ERS + GODERR - OLDERR
C                                  Improve local error estimates proportionally
      SUMERR = ERR(I) + ERR(NUMINT)
      IF (SUMERR .GT. TINY) THEN
         FAC = GODERR / SUMERR 
      ELSE
         FAC = 1.
      ENDIF
      
      ERR(I)      = ERR(I) * FAC
      ERR(NUMINT) = ERR(NUMINT) * FAC
 
      RETURN
C                        ****************************
      END
 
      SUBROUTINE ADZ2AL (F,I)
C                                                   -=-=- adz2al
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D1 = 1.0, D2 = 2.0, HUGE = 1.E15)
C                        Fill in details of interval I given endpoints
      EXTERNAL F
      PARAMETER (MAXINT = 1000)
      COMMON / ADZ2RK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES, 
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
     > ICTA, ICTB, NUMINT, IB
 
      SAVE / ADZ2RK /

      DX =  V(I) - U(I)
      W  = (U(I) + V(I)) / 2.
     
      IF (I .EQ. 1 .AND. ICTA .GT. 0) THEN
C                                                                 Open LEFT end
        FW(I) = FA
        FA = F (U(I) + DX / 4.)

        CALL SGL2NT (ICTA, FA, FW(I), FV(I), DX, TEM, ER)
      ELSEIF (I .EQ. IB .AND. ICTB .GT. 0) THEN
C                                                                open RIGHT end
        FW(I) = FB
        FB = F (V(I) - DX / 4.)
        CALL SGL2NT (ICTB, FB, FW(I), FU(I), DX, TEM, ER)
      ELSE
C                                                                   Closed endS
        FW(I) = F(W)
        TEM = DX * (FU(I) + 4. * FW(I) + FV(I)) / 6.
C                                       Preliminary error Simpson - trapezoidal:
        ER  = DX * (FU(I) - 2. * FW(I) + FV(I)) / 12.
      ENDIF
 
      RESULT(I) = TEM         
      ERR   (I) = ABS (ER)
 
      RETURN
C                        ****************************
      END

      FUNCTION INT2SZ (X, FX)
C                                                   -=-=- int2sz
C                    Return number of distinct points used in AdzInt
C                  Also returns the x- and F(x)-values (ie. the integrand)
C                       at these points (perhaps for plotting purposes).
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (MAXINT = 1000)

      COMMON / ADZ2RK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES, 
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
     > ICTA, ICTB, NUMINT, IB

      Dimension X(MAXINT), FX(MAXINT)

      SAVE / ADZ2RK /

       Call Sort2 (NumInt, U, Fu)
       Do 10 Ix = 1, NUMINT
         X(Ix)  =  U(Ix)
         FX(Ix) = Fu(Ix)
   10  Continue

       Call Sort2 (NumInt, V, Fv)
       X (NumInt + 1) = V (NumInt)
       FX(NumInt + 1) = Fv(NumInt)

      INT2SZ = NUMINT + 1
      RETURN
C                        ****************************
      END
C
      SUBROUTINE SGL2NT (IACT, F1, F2, F3, DX, FINT, ESTER)
C                                                   -=-=- sgl2nt

C     Calculate end-interval using open-end algorithm based on function values
C     at three points at (1/4, 1/2, 1)DX from the indeterminant endpoint (0).

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)

      DATA HUGE / 1.E20 /
C                                                         Use quadratic formula
      TEM = DX * (4.*F1 + 3.*F2 + 2.*F3) / 9.
C                 Error est based on Diff between quadratic and linear integrals
      ER  = DX * (4.*F1 - 6.*F2 + 2.*F3) / 9.

C                          Invoke adaptive singular parametrization if IACT = 2
C                      Algorithm is based on the formula F(x) = AA + BB * x **CC
C                 where AA, BB & CC are determined from F(Dx/4), F(Dx/2) & F(Dx)

      IF (IACT .EQ. 2) THEN
          T1 = F2 - F1
          T2 = F3 - F2
          IF (T1*T2 .LE. 0.) GOTO 7
          T3  = T2 - T1
          IF (ABS(T3)*HUGE .LT. T1**2) GOTO 7
          CC  = LOG (T2/T1) / LOG(D2)
          IF (CC .LE. -D1)  GOTO 7
          BB  = T1**2 / T3
          AA  = (F1*F3 - F2**2) / T3
C                                          Estimated integral based on A+Bx**C
          TMP = DX * (AA + BB* 4.**CC / (CC + 1.))
C                                       Error estimate based on the difference
          ER = TEM - TMP
C                                              Use the improved integral value
          TEM= TMP 
      ENDIF

    7 FINT = TEM
      ESTER= ER
      RETURN
C                        ****************************
      END
     
      SUBROUTINE TOT2LZ
C                                                   -=-=- tot2lz
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (MAXINT = 1000)
      COMMON / ADZ2RK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES, 
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
     > ICTA, ICTB, NUMINT, IB

      SAVE / ADZ2RK /
      RES = 0.
      ERS = 0.
      DO 10  I = 1, NUMINT
          RES = RES + RESULT(I)
          ERS = ERS + ERR(I)
   10     CONTINUE
C                        ****************************
      END

      SUBROUTINE SORT2(N,RA,RB)
C                                                   -=-=- sort2
C                                                  
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      DIMENSION RA(N),RB(N)
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA=RA(L)
          RRB=RB(L)
        ELSE
          RRA=RA(IR)
          RRB=RB(IR)
          RA(IR)=RA(1)
          RB(IR)=RB(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            RA(1)=RRA
            RB(1)=RRB
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(RA(J).LT.RA(J+1))J=J+1
          ENDIF
          IF(RRA.LT.RA(J))THEN
            RA(I)=RA(J)
            RB(I)=RB(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(I)=RRA
        RB(I)=RRB
      GO TO 10
C                    ************************
      END

CsB --------------------------------------------------------------------
C     Testing routines for gg -> AAg
CsB --------------------------------------------------------------------

      Subroutine TestggaagRoutines
      
      Implicit Double Precision (A-H,O-Z)

      Real*8 pLab,pCEM,wHat
      Common /Mom/ pLab(4,12),pCEM(4,12),wHat

CsB   A set of simple momenta to test the polarixation vectors
!     pLab(1,1) =  0.
!     pLab(2,1) =  0.
!     pLab(3,1) = -.5 ! corresponds to p(3,3) =  .5
!     pLab(4,1) = -.5 ! corresponds to p(4,3) =  .5
!
!     pLab(1,2) =  0.
!     pLab(2,2) =  .5 ! corresponds to p(3,4) = -.5
!     pLab(3,2) =  0.
!     pLab(4,2) = -.5 ! corresponds to p(4,4) =  .5
!
!     pLab(1,2) =  0.
!     pLab(2,2) =  0.
!     pLab(3,2) =  .5 ! corresponds to p(3,4) = -.5
!     pLab(4,2) = -.5 ! corresponds to p(4,4) =  .5
!
!     pLab(1,3) =  0.00005
!     pLab(2,3) =  .176776695
!     pLab(3,3) =  .176776695
!     pLab(4,3) = -.25
!
!     pLab(1,4) = -0.00005
!     pLab(2,4) =  .176776695
!     pLab(3,4) =  .176776695
!     pLab(4,4) = -.25
! 
!     pLab(1,5) =  0.
!     pLab(2,5) =  .35355339
!     pLab(3,5) =  .35355339
!     pLab(4,5) =  .5 ! corresponds to p(4,5) = -.5

      Call gggAAA2(amp2)
!     Print*, ' Amp2  =',amp2

CsB   Tests
!     Call TestEpsilon
!     Call TestPolVec
!     Call TestInvariants
!     Call TestPolarizationVectors
!     Call TestTraces
!     Call TestDilogs
!     Call TestI2n
!     Call TestF1
      Call TestAs

!     Stop
      Return
      End

CsB --------------------------------------------------------------------

      Subroutine TestEpsilon
      
      Implicit Double Precision (A-H,O-Z)

      External epsilon
      
      Do i1 = 1, 4
      Do i2 = 1, 4
      Do i3 = 1, 4
      Do i4 = 1, 4
        eps = epsilon(i1,i2,i3,i4)
        Print*, i1,i2,i3,i4, eps
      EndDo
      EndDo
      EndDo
      EndDo

!     Stop
      Return
      End

CsB --------------------------------------------------------------------

      Subroutine TestDilogs
      
      Implicit Double Precision (A-H,O-Z)
      Real*8 p,s
      Complex pc
      Common / sij / p(4,5),pc(4,5),s(5,5)
      Double Complex fLi2,arg,val

CsB   Checked fLi2 arguments and values against Mathematica.
C     In Mathematica the imaginary part of PolyLog[2,x] is always negative.
C     Here, the imaginary part of fLi2 is always positive.
      Do 100 i = 1, 5
      Do 100 j = 1, 5
      Do 100 k = 1, 5
      Do 100 l = 1, 5
        If (Abs(s(i,j)).GT.1e-9 .and. Abs(s(k,l)).GT.1e-9) then
          arg = 1.d0 - s(i,j)/s(k,l)
          val = fLi2(arg)
          Print*,'num,den=',s(i,j),s(k,l),'| Li2(',arg,')=', val
        End If
  100 Continue

!     Stop
      Return
      End

CsB --------------------------------------------------------------------

      Subroutine TestI2n
      
      Implicit Double Precision (A-H,O-Z)
      Real*8 p,s
      Complex pc
      Common / sij / p(4,5),pc(4,5),s(5,5)
      Complex I2n
      
      Do 100 i = 1, 5
        k = i + 1
        Do 100 j = k, 5
          Print*, 'I2n(',s(i,j),')=', I2n(i,j)
  100 Continue

      Stop
      Return
      End

CsB --------------------------------------------------------------------

      Subroutine TestF1
      
      Implicit Double Precision (A-H,O-Z)
      Real*8 p,s
      Complex pc
      Common / sij / p(4,5),pc(4,5),s(5,5)
      Complex F1,tmp
      
      Do 100 i = 1, 5
      Do 100 j = 1, 5
      Do 100 k = 1, 5
      Do 100 l = 1, 5
      Do 100 m = 1, 5
      Do 100 n = 1, 5
        If (i.Eq.j .or. k.Eq.l .or. m.Eq.n) Goto 100
        tmp = F1(i,j,k,l,m,n)
        Print*, 'F1(',s(i,j),',',s(k,l),',',s(m,n),')=', tmp
  100 Continue

      Stop
      Return
      End

CsB --------------------------------------------------------------------

      Subroutine TestInvariants

      Implicit Double Precision (A-H,O-Z)
      Dimension sij(10)
      Real*8 pLab,pCEM,wHat
      Common /Mom/ pLab(4,12),pCEM(4,12),wHat
      Real*8 p,s
      Complex pc
      Common / sij / p(4,5),pc(4,5),s(5,5)

CsB   A set of simple momenta for testing
!     Do i = 1,4
!       pLab(i,3) = -1./2.
!       pLab(i,4) = -2./2.
!       pLab(i,1) =  3./2.
!       pLab(i,2) =  4./2.
!       pLab(i,5) = -5./2.
!     EndDo

!     Call ConstructInvariants(iErr)

CsB   Look at the momenta
      Do i=1,5
        Write(*,'(A,i1,A,4f16.10)') 
     ,  'p(',i,')=',p(1,i),p(2,i),p(3,i),p(4,i)
      EndDo
!     Print*, 'q =', q(1),q(2),q(3),q(4)

CsB   Check momentum conservation
      Write(*,'(A,4f14.8)') 'Sum{i=1,5}p(i)=', 
     ,         p(1,1)+p(1,2)+p(1,3)+p(1,4)+p(1,5),
     ,         p(2,1)+p(2,2)+p(2,3)+p(2,4)+p(2,5),
     ,         p(3,1)+p(3,2)+p(3,3)+p(3,4)+p(3,5),
     ,         p(4,1)+p(4,2)+p(4,3)+p(4,4)+p(4,5)

CsB   Look at the invariants
      Do i = 1, 5
      Do j = i, 5
        Write(*,'(A,i1,A,i1,A,1f18.10)'),'2 p(',i,').p(',j,')=',s(i,j)        
      EndDo
      EndDo

CsB   Check s_{i,j} - s_{k,l}| { 8*1d-5 }
      k = 1
      Do i = 1, 5
        i1 = i + 1
        Do j = i1, 5
        sij(k) = s(i,j)
        k = k + 1
      EndDo
      EndDo
      Do k = 1, 10
        k1 = k + 1
        Do l = k1, 10
          Print*, '|sij(',k,')|-|sij(',l,')|/"+"=',
     ,         (Abs(sij(k))-Abs(sij(l)))/(Abs(sij(k))+Abs(sij(l))) 
      EndDo
      EndDo

!     Stop
      Return
      End

CsB --------------------------------------------------------------------

      Subroutine TestPolarizationVectors
      
      Implicit Double Precision (A-H,O-Z)
      Complex eps
      Common / PolarizationVectors / eps(4,5,2)
      Real*8 p,s
      Complex pc
      Common / sij / p(4,5),pc(4,5),s(5,5)
      Complex CDot
      Complex peps,epep
      Common / epsProducts / peps(5,5,2),epep(5,2,5,2)

!     Do i = 1, 5
!       Print*, 'p(',i,')=', p(1,i),p(2,i),p(3,i),p(4,i)
!     EndDo

CsB   Look at the polarization vectors
!     Do i = 1, 5
!       Print*, 'e+(',i,')=',eps(1,i,2),eps(2,i,2),eps(3,i,2),eps(4,i,2)
!       Print*, 'e-(',i,')=',eps(1,i,1),eps(2,i,1),eps(3,i,1),eps(4,i,1)
!     EndDo

CsB   Test e+ = (e-)*
      Do i = 1, 5
        Print*, 'e+(',i,')-(e-(',i,'))*=', 
     ,           eps(1,i,2)-Conjg(eps(1,i,1)),
     ,           eps(2,i,2)-Conjg(eps(2,i,1)),
     ,           eps(3,i,2)-Conjg(eps(3,i,1)),
     ,           eps(4,i,2)-Conjg(eps(4,i,1))
      EndDo

CsB   Test the Lorentz condition: e_i^+/-.p_i = 0
      Do i = 1, 5
        Print*, 'e+(',i,').p(',i,')=', peps(i,i,2)
        Print*, 'e-(',i,').p(',i,')=', peps(i,i,1)
      EndDo

CsB   Test e_i^+.e_i^+ = e_i^-.e_i^- = 0
      Do i = 1, 5
        Print*, 'e+(',i,').e+(',i,')=', epep(i,2,i,2)
        Print*, 'e-(',i,').e-(',i,')=', epep(i,1,i,1)
      EndDo

CsB   Test e_i^+.e_i^- = e_i^-.e_i^+ = -1
      Do i = 1, 5
        Print*, 'e+(',i,').e-(',i,')=', epep(i,2,i,1)
        Print*, 'e-(',i,').e+(',i,')=', epep(i,1,i,2)
      EndDo

!     Stop
      Return
      End

CsB --------------------------------------------------------------------

      Subroutine TestTraces
      
      Implicit Double Precision (A-H,O-Z)

      Complex TrFF,TrFFF,piFjlpk

      Do i = 1, 5
      Do j = i, 5
        Print*,'Tr(F',i,'+  F',j,'+)=',TrFF(i,j,2)
      End Do
      End Do

      Do i = 1, 5
      Do j = i, 5
        Print*,'Tr(F',i,'+  F',j,'+)-Tr(F',i,'-  F',j,'-)*=',
     ,          TrFF(i,j,2)-Conjg(TrFF(i,j,1))
      End Do
      End Do

      Do i = 1, 5
      Do j = i, 5
      Do k = j, 5
        Print*,'Tr(F',i,'+  F',j,'+  F',k,'+)=',TrFFF(i,j,k,2)
      End Do
      End Do
      End Do

      Do i = 1, 5
      Do j = i, 5
      Do k = j, 5
        Print*,'p(',i,').F',j,'^-.p(',k,')=',piFjlpk(i,j,1,k)
        Print*,'p(',i,').F',j,'^+.p(',k,')=',piFjlpk(i,j,2,k)
      End Do
      End Do
      End Do

!     Stop
      Return
      End

CsB --------------------------------------------------------------------

      Subroutine TestAs
      
      Implicit Double Precision (A-H,O-Z)

      Complex Appppp,Ampppp,Appppm,Apppmm,Ammppp,Ampppm,tmp
      Complex gggAAA
      Common / gggAAAmplitude / gggAAA(2,2,2,2,2,2)
      Data Pi / 3.1415926535897932 /

      Common / Test / iTest

      Call TestInvariants
      AmpSq=PentagonGGAAG(3,4,5,1,2)

      Sum1 = 0.d0
      Sum2 = 0.d0
      Do 100 i1 = 1, 2
      Do 100 i2 = 1, 2
      Do 100 i3 = 1, 2
      Do 100 i4 = 1, 2
      Do 100 i5 = 1, 2
        Sum1 = Sum1 +
     ,   Real(Conjg(gggAAA(i1,i2,i3,i4,i5,1))*gggAAA(i1,i2,i3,i4,i5,1))
        Sum2 = Sum2 +
     ,   Real(Conjg(gggAAA(i1,i2,i3,i4,i5,2))*gggAAA(i1,i2,i3,i4,i5,2))
        Write (*,'(x,5i1,2f14.6,1f12.6)') i1,i2,i3,i4,i5, (8.d0*Pi**4)*
     ,   Real(Conjg(gggAAA(i1,i2,i3,i4,i5,1))*gggAAA(i1,i2,i3,i4,i5,1)),
     ,   Real(Conjg(gggAAA(i1,i2,i3,i4,i5,2))*gggAAA(i1,i2,i3,i4,i5,2)),
     ,   (8.d0*Pi**4)*
     ,   Real(Conjg(gggAAA(i1,i2,i3,i4,i5,1))*gggAAA(i1,i2,i3,i4,i5,1))/
     ,   Real(Conjg(gggAAA(i1,i2,i3,i4,i5,2))*gggAAA(i1,i2,i3,i4,i5,2))
  100 Continue
      Print*,' --- Sum1,Sum2,Sum1/Sum2 =',(8.d0*Pi**4)*Sum1,Sum2,
     ,                                    (8.d0*Pi**4)*Sum1/Sum2, ' ---'

!     Print*,'|A--+++|^2=',
!    ,       Conjg(gggAAA(1,1,2,2,2,2))*gggAAA(1,1,2,2,2,2)
!     Print*,'|A-++++|^2=',
!    ,       Conjg(gggAAA(1,2,2,2,2,2))*gggAAA(1,2,2,2,2,2)
!     Print*,'|A++++-|^2=',
!    ,       Conjg(gggAAA(2,2,2,2,1,2))*gggAAA(2,2,2,2,1,2)
!     Print*,'|A+++++|^2=',
!    ,       Conjg(gggAAA(2,2,2,2,2,2))*gggAAA(2,2,2,2,2,2)

      Write(*,'(x,A,2f14.6)') '|A+++--|^2 ->', (8.d0*Pi**4)*Real(
     ,       Conjg(gggAAA(2,2,2,1,1,1))*gggAAA(2,2,2,1,1,1)/
     /       Conjg(gggAAA(2,2,2,1,1,2))/gggAAA(2,2,2,1,1,2)),
     ,                        (8.d0*Pi**4)*Real(
     ,       Conjg(gggAAA(1,1,1,2,2,1))*gggAAA(1,1,1,2,2,1)/
     /       Conjg(gggAAA(1,1,1,2,2,2))/gggAAA(1,1,1,2,2,2))
      Write(*,'(x,A,2f14.6)') '|A++++-|^2 ->', (8.d0*Pi**4)*Real(
     ,       Conjg(gggAAA(2,2,2,2,1,1))*gggAAA(2,2,2,2,1,1)/
     /       Conjg(gggAAA(2,2,2,2,1,2))/gggAAA(2,2,2,2,1,2)),
     ,                        (8.d0*Pi**4)*Real(
     ,       Conjg(gggAAA(1,1,1,1,2,1))*gggAAA(1,1,1,1,2,1)/
     /       Conjg(gggAAA(1,1,1,1,2,2))/gggAAA(1,1,1,1,2,2))
      Write(*,'(x,A,2f14.6)') '|A+++++|^2 ->', (8.d0*Pi**4)*Real(
     ,       Conjg(gggAAA(2,2,2,2,2,1))*gggAAA(2,2,2,2,2,1)/
     /       Conjg(gggAAA(2,2,2,2,2,2))/gggAAA(2,2,2,2,2,2)),
     ,                        (8.d0*Pi**4)*Real(
     ,       Conjg(gggAAA(1,1,1,1,1,1))*gggAAA(1,1,1,1,1,1)/
     /       Conjg(gggAAA(1,1,1,1,1,2))/gggAAA(1,1,1,1,1,2))
      Write(*,'(x,A,2f14.6)') '|A-++++|^2 ->', (8.d0*Pi**4)*Real(
     ,       Conjg(gggAAA(1,2,2,2,2,1))*gggAAA(1,2,2,2,2,1)/
     /       Conjg(gggAAA(1,2,2,2,2,2))/gggAAA(1,2,2,2,2,2)),
     ,                        (8.d0*Pi**4)*Real(
     ,       Conjg(gggAAA(2,1,1,1,1,1))*gggAAA(2,1,1,1,1,1)/
     /       Conjg(gggAAA(2,1,1,1,1,2))/gggAAA(2,1,1,1,1,2))
      Write(*,'(x,A,2f14.6)') '|A--+++|^2 ->', (8.d0*Pi**4)*Real(
     ,       Conjg(gggAAA(1,1,2,2,2,1))*gggAAA(1,1,2,2,2,1)/
     /       Conjg(gggAAA(1,1,2,2,2,2))/gggAAA(1,1,2,2,2,2)),
     ,                        (8.d0*Pi**4)*Real(
     ,       Conjg(gggAAA(2,2,1,1,1,1))*gggAAA(2,2,1,1,1,1)/
     /       Conjg(gggAAA(2,2,1,1,1,2))/gggAAA(2,2,1,1,1,2))
      Write(*,'(x,A,2f14.6)') '|A-+++-|^2 ->', (8.d0*Pi**4)*Real(
     ,       Conjg(gggAAA(1,2,2,2,1,1))*gggAAA(1,2,2,2,1,1)/
     /       Conjg(gggAAA(1,2,2,2,1,2))/gggAAA(1,2,2,2,1,2)),
     ,                        (8.d0*Pi**4)*Real(
     ,       Conjg(gggAAA(2,1,1,1,2,1))*gggAAA(2,1,1,1,2,1)/
     /       Conjg(gggAAA(2,1,1,1,2,2))/gggAAA(2,1,1,1,2,2))

!     iTest = 1
!     Print*, 'iTest =',iTest

!     Stop
      Return
      End

CsB --------------------------------------------------------------------
C     End routines for gg -> AAg
CsB --------------------------------------------------------------------


C --------------------------------------------------------------------------
      SUBROUTINE WGRIDPH(Y,WT)
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
      Real*8 XMZ_IN,XMZ_GRID,GAMZ_IN
      COMMON/ZMASS_IN/XMZ_IN,XMZ_GRID,GAMZ_IN

      REAL*8 YM12MIN,YM12MAX,QTMin,QTMax, TMmin, Tmmax, ETmin
      COMMON/MYYM12/YM12MIN,YM12MAX,QTMin,QTMax, TMmin, Tmmax, ETmin

      REAL*8 Y,WT
      DIMENSION Y(50)
      INTEGER I
      REAL*8 RN
      EXTERNAL RN

      INTEGER*4 IMONTE
      REAL*8 B_QT2,C_QT2,Q2_V,WT_Q2,QT2_V,WT_QT2,WT_YB,alpha1
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

      REAL*8 QT_V_MIN, QT_V_MAX,YB_V_MAX,Q_V_MIN,Q_V_MAX,QT_V_MXY
      COMMON/HARDCUT/ QT_V_MIN, QT_V_MAX,YB_V_MAX,
     >  Q_V_MIN,Q_V_MAX,QT_V_MXY

CCPY  USED FOR SEPERATING NLO AND NLO_SIG CONTRIBUTIONS
      REAL*8 QT_V_SEP
      COMMON/NLO/ QT_V_SEP

      REAL*8 T_NBAR,T_N
      DIMENSION T_NBAR(4),T_N(4)
      LOGICAL TESTING
      REAL*8 QT2_V_MIN,QT2_V_MAX,TEMP0,TEMP1,TEMP
      REAL*8 Y_EMAX
      Character*40 Boson
      Common / Boson / Boson
!      Integer N_SUD_AB,N_SUD_A,N_SUD_B,N_WIL_C,I_SUD
!      Common / Sudakov / N_SUD_AB,N_SUD_A,N_SUD_B,N_WIL_C,I_SUD
      REAL*8 PCM_L1,ECM_L1,ECM_L2,BETA
      Integer ICalls
      Common / Calls / ICalls
      Integer iYPiece, iReadIn
      Common / YPiece / iYPiece, iReadIn
      Data TESTING /.False./
      INTEGER*4 I_DECAY,I_HADRONIZE
      COMMON/DEC_AY/I_DECAY,I_HADRONIZE

!      Real*8 QT,MT
      
      Character*2 hDecayMode
      Common / hDM / hDecayMode

      Real*8 Reson,xNum
      External Reson
      
      Double Precision ScaPro
      External ScaPro

      TESTING=.false.
c      TESTING=.TRUE.
c      If (Testing) Print*, ' in WGRIDPH '

      PI=4.0*ATAN(1.0D0)
      EPS=1.D-2
      WT=-1.0

      Q_V_MIN=YM12MIN
      Q_V_MAX=YM12MAX

      IF(ABS(Q_V_MAX-Q_V_MIN).GT.EPS) THEN
        PRINT*,' THIS ONLY WORKS FOR ONE Q_V VALUE '
        CALL EXIT
      ENDIF
      Q_V=Q_V_MAX
      Q2_V=Q_V**2
      WT_Q2=1.D0

C FIND QT2_V
CsB The high qT end of the Y piece is read from the Y grid.
cC      If (QT_V_MXY.LE.1.d-6) QT_V_MXY = 100.d0
C ASSUME THE SHAPE AS d(sigma)/d(QT2_V) = c* EXP(-b*QT2_V)
C WHERE b IS FITTED WITH RESUMMED RESULT OF UNPOLARIZED W-BOSON CALCULATIONS
C QT2_V_MAX IS THE MAXIMUM QT2 OF V
CCPY INPUT FROM FIT TO W^+ AT 1.8 TEV, P-PBAR

        B_QT2=(5.9-1.2)/(QT_V_MAX**2-0.1)
        QT2_V_MAX=QT_V_MAX**2


c      print*,' qt_v_max, qtmax =',qt_v_max, qtmax 
      
cdump
c      print*,' y(1) =',y(1)
c      print*,' y(2) =',y(2)
c      print*,' Q_V_Min, Q_V_Max, q2_v, QT2_V_MAX'
c      print*,Q_V_Min, Q_V_Max, q2_v, QT2_V_MAX
      
        IMONTE=0
        If (iMonte.Eq.2) then ! assumes: d(sigma)/d(QT^2) ~ (1/QT2_V)
C_________HARD-WIRED CUT
          QT2_V=QT2_V_MIN*(QT2_V_MAX/QT2_V_MIN)**Y(1)
          QT_V=DSQRT(QT2_V)
          WT_QT2=QT2_V*DLOG(QT2_V_MAX/QT2_V_MIN)
        Else if (iMonte.Eq.3) then ! assumes: d(sigma)/d(QT^2) = 1/X * LOG(X)
C_________WITH X = QT^2/Q^2
C_________HARD-WIRED CUT
          TEMP0=LOG(QT2_V_MIN/Q2_V)**2
          TEMP=LOG(QT2_V_MAX/Q2_V)**2-TEMP0
          C_QT2=2.0/TEMP
          TEMP1=SQRT(TEMP0+TEMP*Y(1))
          QT2_V=Q2_V*EXP(-TEMP1)
          QT_V=DSQRT(QT2_V)
          WT_QT2=-1.0/C_QT2/TEMP1*QT2_V
        Else if (iMonte.Eq.0) then ! RESUMMED CASE
C HARD-WIRED CUT
          QT2_V_MIN=QT_V_MIN**2
          TEMP0=EXP(-B_QT2*QT2_V_MIN)
          TEMP=TEMP0-EXP(-B_QT2*QT2_V_MAX)
          C_QT2=B_QT2/TEMP
          QT2_V=-DLOG(TEMP0-TEMP*Y(1))/B_QT2
          QT_V=DSQRT(QT2_V)
          WT_QT2=1.0/C_QT2/(TEMP0-TEMP*Y(1))
        End if

c      print*,' Qt2_V_Min, Qt2_V_Max, qt2_v '
c      print*,Qt2_V_Min, Qt2_V_Max, qt2_v

C GET YBOOST
C      Print*, ' CALL YMAXIMUM '
      CALL YMAXIMUM(W,Q_V,QT_V,YBMAX)

cdump      print*,' W,Q_V,QT_V,YBMAX =',W,Q_V,QT_V,YBMAX
      
C HARD-WIRED CUT
c      IF(YBMAX.GT.YB_V_MAX) YBMAX=YB_V_MAX

      WHAT=Q_V
      TAU=(WHAT/W)**2
      YBOOST=YBMAX*(2.0D0*Y(2)-1.0D0)
      YB_V=YBOOST
      WT_YB=2.0D0*YBMAX
      X1=DSQRT(TAU)*DEXP(YBOOST)
      X2=TAU/X1

cdump      print*,' w,q_v,ybmax,yboost,YB_V_MAX '
c      print*,w,q_v,ybmax,yboost,YB_V_MAX
      
C ADD THE FACTOR FROM d^4(Q): PI/2
      WT=PI/2.0D0*WT_QT2*WT_Q2*WT_YB


C ADD THE OTHER FACTORS
      WT=WT*
     &   1.d0/(Pi*qT_V)                    ! <= d4Q = pi/2 dqT^2 dQ dy
 
      IF(TESTING) THEN
        PRINT*,'X1,X2,TAU,YBOOST'
        PRINT*,X1,X2,TAU,YBOOST
        PRINT*,'Q_V,QT_V,YB_V'
        PRINT*,Q_V,QT_V,YB_V
        PRINT*,'WT_QT2,WT_Q2,WT_YB,WT'
        PRINT*,WT_QT2,WT_Q2,WT_YB,WT
      ENDIF

C**************************************

C GIVE A RANDOM CHOICE OF AZIMUTHAL ANGLE FOR W-BOSON, IN LAB FRAME
      PHI_Q=2.0*PI*RN(2)
      TM2_Q=Q2_V+QT2_V
      TM_Q=SQRT(TM2_Q)

C (vector boson) 4 momentum in lab frame
      P_Q(1)=QT_V*COS(PHI_Q)
      P_Q(2)=QT_V*SIN(PHI_Q)
      P_Q(3)=TM_Q*(EXP(YB_V)-EXP(-YB_V))/2.0D0
      P_Q(4)=TM_Q*(EXP(YB_V)+EXP(-YB_V))/2.0D0

      P_Q(5)=SQRT(-P_Q(1)**2-P_Q(2)**2-P_Q(3)**2+
     >        P_Q(4)**2)

      IF(TESTING) THEN
        PRINT*,' P_Q, Q_V = ', P_Q,Q_V
      ENDIF

      IF(ABS(P_Q(5)-Q_V).GT.EPS) THEN
       PRINT*,' ERROR IN P_Q'
      ENDIF

CsB Fill the lab momenta PLAB for final particles 1,2 and 9=1&2
      DO 90 I=1,4
C NO INFORMATION ON TWO INCOMING PARTONS: SET THEM TO BE ZERO
        PLAB(I,1)=0.D0
        PLAB(I,2)=0.D0
C VECTOR BOSON
        PLAB(I,3)=P_Q(I)
90    CONTINUE

      IF(TESTING) THEN
        PRINT *,'P_Q  =  ',P_Q
      ENDIF

      RETURN
      END ! WGRIDPH

C ==========================================================================
C --------------------------------------------------------------------------
      SUBROUTINE WGRID(WT)
C --------------------------------------------------------------------------
      IMPLICIT NONE

      REAL*8 WT
      REAL*8 W,QCDL,X1X2MIN
      INTEGER*4 NJET,N_PRIMARY_PART,NF,I_STR,IBEAM,NPART,NEVENT,NORDER,
     ,          N_SECOND_PART
      INTEGER*4 IW
      REAL*8 PLAB,PCEM,WHAT,X1,X2,Q,Q1,Q2
      INTEGER*4 IERR_WT
      REAL*8 XMASS,X_DECAYMASS

      COMMON/AA/W,QCDL,X1X2MIN,NJET,N_PRIMARY_PART,NF,I_STR,IBEAM,NPART,
     &  NEVENT,NORDER,N_SECOND_PART
      COMMON/CHARGE/IW
      COMMON/MOM/PLAB(4,12),PCEM(4,12),WHAT
      COMMON/PMF/X1,X2
CsB 5/11/97 - Q is not used in this routine.
      COMMON/SCALE/Q
      integer nflmax
      parameter (nflmax=5)
      COMMON/STRUCTURE/q1(-nflmax:nflmax),q2(-nflmax:nflmax)
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
      Parameter (MD1 = 200, MD2 = 300, MD3 = 100)
      Dimension X1A(0:MD1), X2A(0:MD2), X3A(0:MD3)
      Dimension YAS(0:MD1,0:MD2,0:MD3), YAA(0:MD1,0:MD2,0:MD3)
      Dimension YDAS(0:MD1,0:MD2,0:MD3), YDAA(0:MD1,0:MD2,0:MD3)
      Dimension X1Y(0:MD1), X2Y(0:MD2), X3Y(0:MD3)
      Dimension YAY(0:MD1,0:MD2,0:MD3,0:7)
      Common / Arrays / X1A,X2A,X3A, YAS,YAA,  X1Y,X2Y,X3Y, YAY
      Common / Dims / iD1,iD2,iD3, iY1,iY2,iY3
      INTEGER N_QT_INTP,N_Y_INTP,N_Q_INTP

CJI Comment out, since they are not used anywhere
CCPY  The second set of grids
C      Real*8 X1A_C,X2A_C,X3A_C, YAS_C,YAA_C, 
C     > YDAS_C,YDAA_C, X1Y_C,X2Y_C,X3Y_C, YAY_C
C      Dimension X1A_C(0:MD1), X2A_C(0:MD2), X3A_C(0:MD3)
C      Dimension YAS_c(0:MD1,0:MD2,0:MD3), YAA_c(0:MD1,0:MD2,0:MD3)
C      Dimension X1Y_C(0:MD1), X2Y_C(0:MD2), X3Y_C(0:MD3)
C      Dimension YAY_C(0:MD1,0:MD2,0:MD3,0:6)
C      Dimension YDAS_C(0:MD1,0:MD2,0:MD3), YDAA_C(0:MD1,0:MD2,0:MD3)
C      Common / Arrays_C / X1A_C,X2A_C,X3A_C, YAS_c,YAA_c,  
C     >X1Y_C,X2Y_C,X3Y_C, YAY_C

      Character*40 Boson
      Common / Boson / Boson
      Integer N_SUD_AB,N_SUD_A,N_SUD_B,N_WIL_C,I_SUD
      Common / Sudakov / N_SUD_AB,N_SUD_A,N_SUD_B,N_WIL_C,I_SUD
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
      Dimension SigY (0:5)
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

      INTEGER iCSSPiece
      Common / CSSPiece / iCSSPiece

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
      
CCAO! Theorectical Soft-Cutoff
      REAL*8 DELTA_SOFT
      Common /TheoCutoff/ DELTA_SOFT
CCAO! Spercify the lepton flavor in W decay
      INTEGER*4 ID_LEPTON, ID_ELE,ID_MU,ID_TAU
      COMMON/LEPTONID/ID_LEPTON,ID_ELE,ID_MU,ID_TAU
CCAO! Input the lepton masses
      REAL*8 XM_ELE,XM_MUON,XM_TAU,XM_NU
      COMMON/LEPTONMASS/ XM_ELE,XM_MUON,XM_TAU,XM_NU
      
      Real*8 propagator
CCPY
      REAL*8 LE,MF,MFS,XNU_PH,XE_PH,FQED,S_QED,FQED1,FQED2,EE

      Data Pi / 3.141592653589793 /
      Data (SigY(ii), ii=0,4) / 0.d0,0.d0,0.d0,0.d0,0.d0 /
      Data FacN / 1.d0 /
      


C This routine is right only in the Collins-Soper frame.

CsB___Initialization
      SigS = 0.d0
      SigA = 0.d0
      SigY(0) = 0.d0
      SigY(1) = 0.d0
      SigY(2) = 0.d0
      SigY(3) = 0.d0
      SigY(4) = 0.d0
      SigY(5) = 0.d0

CCPY ORDER OF POLYNOMIAL USED IN INTERPOLATION
      N_QT_INTP=4
      N_Y_INTP=4
      N_Q_INTP=4

CsB_____Resummed case
          If (iCSSPiece.eq.1 .and. iYPiece.eq.1) then
CsB_________If iYPSw=1 then send the Pert piece back
CsB_________If iYPSw=0 then send the Y piece back
CsB_______For the matched resummed case the CSS piece is included in the Y grid
C         This way the CSS piece interpolation does not have to be called.
            SigA = 0.d0
            SigS = 0.d0
CsB_______Interpolate the Y piece
            iYPSw = 0
            Call ThreeYInterp ( N_QT_INTP,N_Y_INTP,N_Q_INTP,
     >        YDAS,YDAA, QT_V,yB_V,Q_V, SigY, ErrSgY )

          ElseIf (iCSSPiece.eq.1 .and. iYPiece.eq.0) then
CsB_______Interpolate the CSS piece
            Call ThreeDInterp ( N_QT_INTP,N_Y_INTP,N_Q_INTP,
     >      YDAS,YDAA, QT_V,yB_V,Q_V,SigS,SigA, ErrSgS,ErrSgA )

          ElseIf (iCSSPiece.eq.0 .and. iYPiece.eq.1) then
            SigA = 0.d0
            SigS = 0.d0
CsB_______Interpolate the Y piece
CsB_________If iYPSw=1 then send the Pert piece back
CsB_________If iYPSw=0 then send the Y piece back
            iYPSw = 0
            Call ThreeYInterp ( N_QT_INTP,N_Y_INTP,N_Q_INTP,
     >        YDAS,YDAA, QT_V,yB_V,Q_V, SigY, ErrSgY )

        End If
        
  
       WT = SigS+SigY(0)
       
cdump      print*,' QT_V,yB_V,Q_V =',QT_V,yB_V,Q_V
c      print*,' SigS, SigY(0) =',SigS, SigY(0)
      
 
      RETURN
      END ! WGRID

C=============================================================




c==================================================
