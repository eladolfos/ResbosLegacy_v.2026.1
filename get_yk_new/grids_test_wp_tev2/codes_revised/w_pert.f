cpn Sep 2007 Add a separate calculation for W- production, including the
c       option  to compute s+g -> W- + c for ido_cbar=1; updated
c      electroweak factors; introduced extrapolation of the K-factor 
c       outside of the kinematically allowed region 
CCPY Jan. 9, 2007: Add the option to return (c-bar + Z0 ) rate, which
C                   is coded through the flag "IDO_CBAR"  
C Dec. 20, 2006: Add the option to return (c-bar + W+ ) rate, which
C                   is coded through the flag "IDO_CBAR"  

      program w_pert
      IMPLICIT NONE        

CCPY===========================================

      REAL*8 ERREL
      INTEGER NINTVL
      COMMON /ERR/ ERREL,NINTVL
cpn2007 Keep errel < 0.01 to get good cancellation of PERT and ASY at
cpn2007 small QT
C      DATA NINTVL,ERREL /0,0.005/
      DATA NINTVL,ERREL /0,0.002/

      CHARACTER*40 PDF_FILE
      COMMON/FILINP/ PDF_FILE
      REAL*8 NC0,CF0,NU0,ND0,NUD0,CNP0,TR0
      REAL*8 CF,CNP,TR
      INTEGER NC,NU,ND,NUD
      COMMON /CLEBSH/ CF,CNP,TR,NC,NU,ND,NUD     
        PARAMETER (NC0=3, CF0=4./3., NU0=2, ND0=3, NUD0 = 2)
        PARAMETER (CNP0=CF0-NC0/2., TR0=(NU0+ND0)/2.)       
        DATA NC,CF,NU,ND,NUD /NC0,CF0,NU0,ND0,NUD0/         
        DATA CNP,TR /CNP0,TR0/

        double precision ss,tt,uu,qq,qm,qmm, qmu
      COMMON /STU/ SS,TT,UU,QQ,QM,QMM,qmu
        double precision xlim,qt,y
      COMMON /X1LIM/ XLIM(2),QT,Y    

      REAL*8 U1,UB1,D1,DB1,ST1,STB1,C1,CB1,B1,BB1,G1
      REAL*8 U2,UB2,D2,DB2,ST2,STB2,C2,CB2,B2,BB2,G2
      COMMON /STRUCT/         
     1         U1,UB1,D1,DB1,ST1,STB1,C1,CB1,B1,BB1,G1      
     2        ,U2,UB2,D2,DB2,ST2,STB2,C2,CB2,B2,BB2,G2      

      REAL*8 CC11,CC12,CC21,CC22
      COMMON /CABIBO/ CC11,CC12,CC21,CC22         
      DATA CC11,CC12 /.9506, .050/     
      DATA CC21,CC22 /.050, .9506/ !V_km^2 for km_l=0.224    

      REAL*8 GGU0,GGD0,GVU0,GVD0,GAU0,GAD0
      REAL*8 UL,UR,DL,DR,SIN2W,GZCOEF
      REAL*8 GF,ZMASS,WMASS,GGW,GGU,GGD,GVU,GVD,GAU,GAD
      COMMON /WPARM/ GF,ZMASS,WMASS,GGW,GGU,GGD,GVU,GVD,GAU,GAD       
CCPY
      REAL*8 GVUXGVD,GAUXGAD,GAUXGAU,GVUXGVU,GADXGAD,GVDXGVD
      COMMON/WPARM2/GVUXGVD,GAUXGAD,GAUXGAU,GVUXGVU,GADXGAD,GVDXGVD   


CsB        PARAMETER (SIN2W=0.227)
CZL        PARAMETER (SIN2W=0.223)         
        PARAMETER (SIN2W=0.23143)         
CZL        PARAMETER (GZCOEF=0.1848)  
        PARAMETER (GZCOEF=0.187523)  
C                                        ! SQRT(GF*ZMASS**2/2./SQRT(2))         
        PARAMETER (UR=-2.*(2./3.)*SIN2W, DR=-2.*(-1./3.)*SIN2W)       
        PARAMETER (UL=1.+UR, DL=-1.+DR) 
        PARAMETER (GVU0=GZCOEF*(UR+UL)  , GVD0=GZCOEF*(DR+DL))        
        PARAMETER (GAU0=GZCOEF*(UR-UL)  , GAD0=GZCOEF*(DR-DL))        
        PARAMETER (GGU0=GVU0**2+GAU0**2 , GGD0=GVD0**2+GAD0**2 )      
        DATA GF,ZMASS,WMASS /1.16637E-5, 91.1875, 80.385/        
CCPY        DATA GGW /0.1056/
C      GGW=SQRT(2.0)*GF*WMASS**2
        DATA GGW /0.1066/

        DATA GGU,GGD,GVU,GVD,GAU,GAD /GGU0,GGD0,GVU0,GVD0,GAU0,GAD0/  

      REAL*8 PI
      DATA PI /3.1415 92653 58979 32384 62643/  

          
C **  JWTYPE -- INDICATES THE TYPE OF PHOTON PRODUCED BY THE FOLLOWING
C **  CODES         
C **      
C **     JWP, JWM      W+, W- BOSON         
C **     JZ            Z BOSON
C **     JPH           PHOTON (NOT IMPLEMENTED)   
          
       integer jwtype,jwp,jwm,jz,jph,JHB
      COMMON /JCODES/ JWTYPE,jwp,jwm,jz,jph,JHB
          DATA jwm,jwp,jz,jph,JHB /-1,1,2,3,4/        

      Character*40 Boson
      integer lepasy
      Common / Boson / Boson,lepasy

      INTEGER MORDER
      COMMON /MCODES/ MORDER

      integer iset, norder
      COMMON/QCD1/ ISET,NORDER
  

C **  INSTAT -- INDICATES THE TYPE OF COLLIDING PARTONS BY THE FOLLOWING        
C **  CODES         
C **      
C ** 1    IQQB1         Q+QBAR: [-->GGQ] + [VIRTUAL] = (F1+F2)**2      
C ** 2    IQG           Q+G    
C ** 3    IGG           G+G    
C ** 4    IQQB          Q+QBAR: EVERYTHING ELSE    
C ** 5    IQQ           Q+Q    
C ** 6    IALL          SUM ALL OF THE ABOVE       
          
      integer INSTAT,  IQQB1,IQG,IGG,IQQB,IQQ,IALL 
      COMMON /ICODES/ INSTAT,  IQQB1,IQG,IGG,IQQB,IQQ,IALL  
        DATA IQQB1,IQG,IGG,IQQB,IQQ,IALL /1,2,3,4,5,6/      
          
C **    QQFLAG    -- IF TRUE, Q+QBAR CONTRIBUTIONS INCLUDED.
C **    QGFLAG    -- IF TRUE, Q+G AND QBAR+G CONTRIBUTIONS INCLUDED.
C **    JWTYPE    -- TYPE OF PHOTON PRODUCED (CODES = JWP,JWM,JZ,JPH)

      COMMON /PROCES/ QQFLAG,QGFLAG
        LOGICAL QQFLAG,QGFLAG
        DATA QQFLAG,QGFLAG /.TRUE.,.TRUE./

C use MS-bar scheme PDF
      COMMON /FLAGS/ DS2FLG,NRMFLG,S2AFLG,FACFLG  
        LOGICAL DS2FLG,NRMFLG,S2AFLG,FACFLG       
        DATA DS2FLG,NRMFLG,S2AFLG,FACFLG /.TRUE.,.TRUE.,.TRUE.,.FALSE./

C **  LTERM -- INDICATES WHICH TERM OF INTEGRAND IS CURRENTLY BEING EVALUATED IN
C **  PARTON MOMENTUM-FRACTION INTEGRALS FOR 2ND-ORDER CONTRIBUTIONS. 
C **      
C **     L2DIM         EVERYTHING BUT THE BELOW   
C **     LS20          COEFFICIENT OF [1/S2] AT S2=0
C **     LWS20         COEFFICIENT OF [LOG(S2)/S2] AT S2=0    
C **     LDS2          DELTA(S2) TERMS  
          
      integer LTERM, L2DIM,LS20,LWS20,LDS2
      COMMON /LCODES/ LTERM, L2DIM,LS20,LWS20,LDS2      
        DATA L2DIM,LS20,LWS20,LDS2 /1,2,3,4/    
C***
      Real*8 QQCsB
      Common /CsBQQ/ QQCsB

CCPY===========================================
      INTEGER IDO_CBAR
      COMMON /CBAR_ONLY/ IDO_CBAR      

      INTEGER NIN1, NOUT1

      CHARACTER*40 QGFN, QTGFN, YGFN
      REAL*8 QG(99),PT(200),RAPY(200)
      COMMON / GRIDFILE / QG,PT,RAPY,QGFN,QTGFN,YGFN
      INTEGER II, N_Q,N_QT,N_Y, LTO
      COMMON / NGRID / N_Q,N_QT,N_Y, LTO

      REAL*8 ECM
      INTEGER IBEAM,JWTYPE_IN
      common/collider/ecm,ibeam
      
      Real*8 pT_Stretch,YMAX
      REAL*8 TANB,VEV,COUP_MB,CR,COUPL

C=============================================

     
      REAL*8 ADYPT
      EXTERNAL ADYPT
      
      INTEGER I,J,K
      REAL*8 VPT,RAPIN,vmas
      INTEGER IQTMN, IQTMX, IQTST, IYMN, IYMX, IYST,
     &        IQMN, IQMX, IQST
      INTEGER IPTMIN, IPTMAX, IPTSTP, IYMIN, IYMAX, IYSTP,
     >        IQMIN, IQMAX, IQSTP

      REAL*8 pert_1,pert_2,factor_k,pert
cpn2007 For extrapolation of the K-factor outside of the kinematically
cpn2007 allowed region
      real*8 a_exterp, b_exterp, vpt_prev, factor_K_prev

      real*8 scale_fac_muR, scale_fac_muF
      
C=============================================

!ZL
	character*100 jobname
	integer length
       double precision ALPI
       CHARACTER*40 BOSON_IN
       integer JZ_TYPE

       REAL*8 PERTL0_1,PERTL0_2,PERTA3_1,PERTA3_2,TERML0,TERMA3,
     & FACTOR_KL0,FACTOR_KA3,DUMMY

CBY Oct26,18
       INTEGER iscale

	if(iargc().eq.0)then
		jobname='w_pert'
	else if(iargc().ge.1) then
		call getarg(1,jobname)
	endif
	call trmstr(jobname,length)
      
C I/O UNITS
      NIN1=23
      NOUT1=22
C
!ZL      OPEN(UNIT=NIN1,FILE='w_pert.in',STATUS='old')
!ZL      OPEN(UNIT=NOUT1,FILE='w_pert.out')
      OPEN(UNIT=NIN1,FILE=jobname(1:length)//'.in',STATUS='old')
      OPEN(UNIT=NOUT1,FILE=jobname(1:length)//'.out')


C=========================
CCPY Turn on all subprocesses
        INSTAT = IALL
C=========================

      READ(NIN1,*) ECM,IBEAM
      SS=ECM**2
      
C DATA jwm,jwp,jz,jph,JHB /-1,1,2,3,4/ 
CCPY August 2020
      IDO_CBAR=0
      JZ_TYPE=0
      lepasy=0
      BOSON='-'
      READ(NIN1,*)JWTYPE_IN,IDO_CBAR,JZ_TYPE
      JWTYPE=JWTYPE_IN
      IF (JWTYPE.eq.JZ) THEN 
        IF (JZ_TYPE.EQ.1) THEN 
          BOSON='ZU'
        ELSEIF (JZ_TYPE.EQ.-1) THEN 
          BOSON='ZD'
        ELSEIF (JZ_TYPE.EQ.0) THEN 
          BOSON='Z0'
        ELSE
          PRINT*,' THIS JZ_TYPE IS NOT YET IMPLEMENTED.'
          CALL EXIT
        ENDIF   
      ENDIF 
      
      IF(BOSON.EQ.'ZU') THEN
C replace GAU**2 by GAUXGAU,and GVU**2 by GVUXGVU, 
        GAUXGAU=GAU*GAU
        GVUXGVU=GVU*GVU
        GVUXGVD= GVU*GVD*GGU/(GGU+GGD)
        GAUXGAD= GAU*GAD*GGU/(GGU+GGD)
C replace GAD**2 by GADXGAD,and GVD**2 by GVDXGVD, 
        GADXGAD=0.D0
        GVDXGVD=0.D0
        GGD=0.D0
      ELSEIF(BOSON.EQ.'ZD') THEN
C replace GAD**2 by GADXGAD,and GVD**2 by GVDXGVD, 
        GADXGAD=GAD*GAD
        GVDXGVD=GVD*GVD
        GVUXGVD= GVU*GVD*GGD/(GGU+GGD)
        GAUXGAD= GAU*GAD*GGD/(GGU+GGD)
C replace GAU**2 by GAUXGAU,and GVU**2 by GVUXGVU, 
        GAUXGAU=0.D0
        GVUXGVU=0.D0
        GGU=0.D0
      ELSE
        GVUXGVD= GVU*GVD
        GAUXGAD= GAU*GAD
C replace GAU**2 by GAUXGAU,and GVU**2 by GVUXGVU, 
        GAUXGAU=GAU*GAU
        GVUXGVU=GVU*GVU
C replace GAD**2 by GADXGAD,and GVD**2 by GVDXGVD, 
        GADXGAD=GAD*GAD
        GVDXGVD=GVD*GVD
      ENDIF

      IF(JWTYPE.EQ.JHB) THEN
        TANB = 50.0
        VEV = 246.220569
CCPU HARD-WIRED RUNNING BOTTOM QUARK MASS at 100 GeV         
        COUP_MB=2.98D0
        CR = TANB*(COUP_MB/VEV)
        COUPL=0.5*CR**2
      ELSE
        COUPL=1.0D0
      ENDIF
      
CCPY
      IF(IDO_CBAR.EQ.1 .AND. (abs(JWTYPE).NE.JWP.AND.JWTYPE.NE.JZ)) THEN 
         WRITE(*,*) ' NOT IMPLEMENTED, EXCEPT FOR W/Z PRODUCTION '
         CALL EXIT
      ENDIF
      
      IF(IDO_CBAR.EQ.1 .AND. (abs(JWTYPE).EQ.JWP.OR.JWTYPE.EQ.JZ)) THEN 
        WRITE(*,*) ' ONLY CBAR FINAL STATE IS INCLUDED '
        WRITE(NOUT1,*) ' ONLY CBAR FINAL STATE IS INCLUDED '
C        PRINT*,' CC11,CC21,CC12,CC22 =', CC11,CC21,CC12,CC22 
C        PRINT*,' NUD =',NUD
        CC11=0.0
        CC12=0.0
        CC21=0.0
        CC22=1.0
        NUD=1
        ND=0
        NU=1
C        PRINT*,' CC11,CC21,CC12,CC22 =', CC11,CC21,CC12,CC22 
C        PRINT*,' NUD,ND,NU =',NUD,ND,NU

C=========================
CCPY Turn on only IQG
C        INSTAT = IQG
CCPY Turn on only IGG
C        INSTAT = IGG        
C=========================
        IF(INSTAT.EQ.IQG) THEN
          WRITE(NOUT1,*) ' FROM G+QBAR ONLY'
        ELSEIF(INSTAT.EQ.IGG) THEN
          WRITE(NOUT1,*) ' FROM G+G ONLY '
        ELSEIF(INSTAT.EQ.IALL) THEN
          WRITE(NOUT1,*) ' FROM G+Q AND G+G  '
        ELSE
          WRITE(NOUT1,*) 'NOT CONSISTENT '
          CALL EXIT
        ENDIF                
        
      ENDIF
 
cpn2006
C      READ(NIN1,*) iset, norder
ccpy July 2013
      READ(NIN1,*) iset, norder, scale_fac_muR, scale_fac_muF,iscale
      
      READ(NIN1,'(A40)') PDF_FILE

cpn2006
      if (iset.eq.902) then
        write (NOUT1,*) 'iset, pdf_file= ',iset,', ',pdf_file
      else
         write (NOUT1,*) 'iset = ', iset
       endif                    !iset
 

CCPY July 2013      WRITE(NOUT1,*) 'ECM, IBEAM, JWTYPE_IN, NORDER, IDO_CBAR = ',
C     >ECM, IBEAM, JWTYPE_IN, NORDER, IDO_CBAR 
    
      WRITE(NOUT1,*) 'ECM, IBEAM, JWTYPE_IN, NORDER, IDO_CBAR, JZ_TYPE,  
     >scale_fac_muR, scale_fac_muF, iscale= ',
     >ECM, IBEAM, JWTYPE_IN, NORDER, IDO_CBAR,JZ_TYPE,
     >scale_fac_muR, scale_fac_muF,iscale
    
      IF(JWTYPE.EQ.JHB) THEN   
        WRITE(NOUT1,25) 
     >    ' TANB,VEV,COUP_MB,CR,COUPL = ',
     >    TANB,VEV,COUP_MB,CR,COUPL
      else 
        write(NOUT1,*)
      ENDIF
 25   FORMAT (A30,5(1x,G10.3))

      READ(NIN1, '(A40)') QGFN
      READ(NIN1, '(A40)') QTGFN
      READ(NIN1, '(A40)') YGFN
      IF (QGFN.NE.'-') THEN
        PRINT*, ' Reading  Q values from ', QGFN
        OPEN(UNIT=31,FILE=QGFN,STATUS='OLD')
        DO II = 1,150
          READ(31,*,END=98) QG(II)
c          Print*, '>>>', ii,  QG(ii)
        END DO
   98   N_Q = II - 1
        CLOSE(31)
      END IF
      IF (QTGFN.NE.'-') THEN
        PRINT*, ' Reading qT values from ', QTGFN
        OPEN(UNIT=31,FILE=QTGFN,STATUS='OLD')
        DO II = 1,200
          READ(31,*,END=99) PT(II)
C          Print*, '>>>', ii, pt(ii)
        END DO
   99   N_QT = II - 1
        CLOSE(31)
      END IF
      IF (YGFN.NE.'-') THEN
        PRINT*, ' Reading  y values from ', YGFN
        OPEN(UNIT=31,FILE=YGFN,STATUS='OLD')
        DO II = 1,200
          READ(31,*,END=100) RAPY(II)
C          Print*, '>>>', ii, y(ii)
        END DO
  100   N_Y = II - 1
        CLOSE(31)
      END IF

      READ(NIN1, *) IQTMN, IQTMX, IQTST, IYMN, IYMX, IYST,
     &             IQMN, IQMX, IQST
      CLOSE(NIN1)

      IQMIN=IQMN
      IQMAX=IQMX
      IQSTP=IQST
      IPTMIN=IQTMN
      IPTMAX=IQTMX 
      IPTSTP=IQTST 
      IYMIN=IYMN
      IYMAX=IYMX
      IYSTP=IYST
 
      If (N_Q.LT.IQMAX) then
        Print*, ' Resetting IQMAX to N_Q = ', N_Q
        IQMAX = N_Q
      End If
      If (N_qT.LT.IpTMAX) then
        Print*, ' Resetting IpTMAX to N_qT = ', N_qT
        IpTMAX = N_qT
      End If
      If (N_y.LT.IyMAX) then
        Print*, ' Resetting IyMAX to N_y = ', N_y
        IyMAX = N_y
      End If

CCPY
       IF (JWTYPE.EQ.JZ) THEN
         WRITE(NOUT1,*) ' Q,qT,y, pertL0_1, pertA3_1,  pertL0_1+pertL0_2,
     > pertA3_1+pertA3_2 '
CCPY     > pertA3_1+pertA3_2, factor_KL0, FACTOR_KA3  '
       ELSE 
CCPY FEb 2021
         WRITE(NOUT1,*) ' Q,qT,y, pertL0_1, pertA3_1,  pertL0_1+pertL0_2,
     > pertA3_1+pertA3_2 '
C         WRITE(NOUT1,*) ' Q,qT,y, pert_1, pert_2, 
C     >  factor_K=(pert_2+pert_1)/pert_1  '
       ENDIF

        PRINT *, ' Q, qT, y, pert_1, pert_2, factor_K, LEPASY '

        pT_Stretch = 1.d0


      call Setthepdf
      
        pert=0.D0
        factor_K=0.D0
        
        DO K = IQMIN, IQMAX, IQSTP
          VMAS = QG(K)
c===================
CZL
!          QM=VMAS
!          QQ=QM**2
!          QQCSB=QQ
c===================   
CCPY July 2013; scale_factor
C          SCALE_FAC=1.0D0
C
c===================   

          DO J = IYMIN, IYMAX, IYSTP
            RAPIN = RAPY(J)
cpn2007     Initial values of parameters for extrapolation outside 
cpn2007     of the kinematically allowed region
            vpt_prev=0d0; factor_K_prev=0d0

            DO I = IPTMIN, IPTMAX, IPTSTP
CCPY
              LEPASY=0

111           CONTINUE
              VPT=PT(I)*pT_Stretch
CJI 2016:

CBY Oct26,18
              if(iscale.eq.0) then
                QM=scale_fac_muF*VMAS
                QMu=scale_fac_muR*VMAS
              else
                QM=scale_fac_muF*sqrt(VMAS**2+VPT**2)
                QMu=scale_fac_muR*sqrt(VMAS**2+VPT**2)
              endif
              QQ=VMAS**2
              QQCSB=QQ

C              PRINT*,k,j,i,vmas,VPT,RAPIN,norder

              CALL YMAXIMUM(ECM,VMAS,VPT,YMAX)
              
C              PRINT*,'VMAS,VPT,YMAX ,RAPIN ='
C              PRINT*,VMAS,VPT,YMAX ,RAPIN 

cpn2007  Compute the exact K-factor at kinematically allowed QT (0<QT<QT_kin);
c        otherwise (at QT > QT_kin) approximate the K-factor by
c        extrapolation from the region QT < QT_kin. This is needed for
c        reasonable interpolation of the K-factor inside the kinematically allowed region 
c        at QT close to QT_kin in ResBos, for which we have to have a
c        reasonable K-factor behavior at two QT values immediately above QT_kin. 

              IF(ABS(RAPIN).LT.YMAX) then
cpn2007       Compute the exact K-factor inside the kinematically allowed region
                IF(norder.EQ.1) THEN
                  call dypt1(vpt,rapin,pert_1)
                  pert_2=0.D0
                  factor_K=0.D0               
                ELSE
                  call dypt1(vpt,rapin,pert_1)
                  call dypt2(vpt,rapin,pert_2)
                  factor_K=(pert_2+pert_1)/pert_1
                ENDIF
CCPY              pert=2.0*VPT*(pert_1+pert_2)
C              pert=pert*COUPL
                
                pert_1=pert_1*2.0*VPT*COUPL
                pert_2=pert_2*2.0*VPT*COUPL
                pert=pert_1+pert_2
              
              else              !outside of the kinematical region
                pert_1=0d0; pert_2=0d0
                factor_K = a_exterp*vpt + b_exterp
              endif             !IF(ABS(RAPIN).LT.YMAX)

              if(pert_1 .lt. 1E-18 .and. pert_2 .lt. 1E-18) then
                  factor_K = 0
              endif

!ZL avoid too large k factor
		if(abs(factor_K).gt.5) then
			print *, 'warning: abs of k factor > 5 : ', 
     > factor_K, vmas,VPT,RAPIN
		endif
CZLdebug                        print *, "alpha_s(mz)=",ALPI(91.1876d0)*PI
C Output the results
CCPY123           PRINT 101,    vmas,VPT,RAPIN,pert_1,pert_1+pert_2,factor_K 
123           PRINT 101,    vmas,VPT,RAPIN,pert_1,pert_2,factor_K,LEPASY 

cpn2007 Compute coefficients of linear extrapolation for the K factor
              a_exterp=(factor_K-factor_K_prev)/(vpt - vpt_prev)
              b_exterp=(factor_K_prev*vpt-factor_K*vpt_prev)/
     >          (vpt - vpt_prev)
              vpt_prev = vpt
              factor_K_prev=factor_K

CCPY
            IF ((JWTYPE.eq.JZ).AND.(LEPASY.EQ.0)) THEN
              PERTL0_1=PERT_1
              PERTL0_2=PERT_2
CCPY Nov 30, 2020: This code cannot correctly predict the case for (LEPASY.EQ.1)
C              LEPASY=1
C              GOTO 111
C            ELSEIF ((JWTYPE.eq.JZ).AND.(LEPASY.EQ.1)) THEN
C              PERTA3_1=PERT_1
C              PERTA3_2=PERT_2
              TERML0=PERTL0_1+PERTL0_2
C              TERMA3=PERTA3_1+PERTA3_2
              FACTOR_KL0=TERML0/PERTL0_1
C              FACTOR_KA3=TERMA3/PERTA3_1
C Hence, we shall assign a strange number -999 here
              DUMMY=-999.0
              PERTA3_1=DUMMY
              TERMA3=DUMMY
              FACTOR_KA3=DUMMY
              WRITE(nout1,122) vmas,VPT,RAPIN,PERTL0_1,PERTA3_1,
     &                     TERML0,TERMA3
CCPY     &                     TERML0,TERMA3,factor_KL0,FACTOR_KA3
            ELSE
              DUMMY=-999.0
C              WRITE(nout1,102) vmas,VPT,RAPIN,
C     &                         pert_1,pert_1+pert_2,factor_K
              WRITE(nout1,122) vmas,VPT,RAPIN,
     &                         pert_1,DUMMY,pert_1+pert_2,DUMMY
            ENDIF

            END DO              ! I (QT)
          ENDDO       ! J (Y)
        ENDDO         ! K (Q)
        CLOSE(NOUT1)
        RETURN
cpn2007 Keep enough digits in pert_1, pert_2, factor_K in order not to
cpn2007 lose precision in the Y-piece K-factor despite large cancellations

  101 FORMAT(1X,3(G10.4,2x),3(G14.7,2x),I4)
  102 FORMAT(1X,3(G10.4,2x),3(G14.7,2x))
  122 FORMAT(1X,3(G10.4,2x),6(G14.7,2x))


      
      END                       !w_pert
      
C --------------------------------------------------------------------------
      SUBROUTINE YMAXIMUM(ECM,Q_V,QT_V,YMAX)
C --------------------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 ECM,Q_V,QT_V,YMAX
      REAL*8 X,ACOSH,T1,T2,T3

      ACOSH(X)=DLOG(X+DSQRT(X**2-1.0))

      T1=ECM**2+Q_V**2
      T2=DSQRT(QT_V**2+Q_V**2)
      T3=2.0*ECM*T2
      YMAX=ACOSH(T1/T3)

      if(YMAX.ne.YMAX) print*, "YMAX is NAN!"

      RETURN
      END


C --------------------------------------------------------------------------
      subroutine setthepdf
C --------------------------------------------------------------------------

      implicit none
      integer iset, norder
      COMMON/QCD1/ ISET,NORDER
      CHARACTER*40 PDF_FILE
      COMMON/FILINP/ PDF_FILE
      integer NIN, NOUT, NWRT
      COMMON / IOUNIT / NIN, NOUT, NWRT
      double precision  ECM
      INTEGER IBEAM

      common/collider/ecm,ibeam

cpn                   
      real*8 xminin, qmaxin
      integer jr, nxin, nqin,nfin,npdffile
      INTEGER II, N_Q,N_QT,N_Y, LTO
      COMMON / NGRID / N_Q,N_QT,N_Y, LTO

c transfer I/O units into QCD library common area from program common area
      NIN=23
      NOUT=71
      NWRT=13
C Set the values of xmin and qmax
      xminin = 1.e-5
      nxin=100

      qmaxin =max(ecm,1.d3)
      nqin=13
      if (lto.eq.0) then  
        qmaxin = 4*qmaxin
        nqin = 22
      endif

      nfin=5
      npdffile=67
      call parpdf(1,'XMIN',xminin,jr)
      call parpdf(1,'QMAX',qmaxin,jr)
      call parpdf(1,'NT',dble(nqin),jr)
      call parpdf(1,'NX',dble(nxin),jr)
      call parpdf(1,'NFL',dble(nfin),jr)
      call parpdf(1,'NFMX',dble(nfin),jr)
      call parpdf(1,'NuIni',dble(npdffile),jr)

      open(npdffile,file=pdf_file,status='unknown')

c Initialize PDFs
      Call SetQCD
      Call SetUtl

      close(npdffile)

      Call SetPDF(iset)

      return 
      end 
