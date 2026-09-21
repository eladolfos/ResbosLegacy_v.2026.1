cpn2007 Sep 2007 -- fixed a bug in the flavor dependence in the
cpn2007 subroutine charge; linked to the updated CTEQ package;
cpn    Add a separate calculation for W- production; updated
c      electroweak factors; introduced extrapolation of the K-factor
c       outside of the kinematically allowed region
CCPY Uisng the asym.f code to calculate the Y-piece up to alpha_s^2
C for  (b + bbar --> A) process


      program w_asym
      IMPLICIT NONE        

CCPY===========================================

      REAL*8 ERREL
      INTEGER NINTVL
      COMMON /ERR/ ERREL,NINTVL      
      DATA NINTVL,ERREL /1,0.005/
cpn2007 Keep errel < 0.01 to get good cancellation of PERT and ASY at
cpn2007 small QT

      CHARACTER*40 PDF_FILE
      integer LEN_PDF_FILE
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
        common /stu/ ss,tt,uu,qq,qm,qmm, qmu
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
      Character*40 Boson
      Common / Boson / Boson
   
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
CZL        PARAMETER (GGU0=2*(GZCOEF)**2*2/0.1473*0.5 , GGD0=0 )      
CZL        PARAMETER (GGU0=GVU0**2+GAU0**2 , GGD0=0 )      
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
C **     JWP, JWM            W+, W- BOSON         
C **     JZ            Z BOSON
C **     JPH           PHOTON (NOT IMPLEMENTED)   
          
       integer jwtype,jwp,jwm,jz,jph,JHB
      COMMON /JCODES/ JWTYPE,jwp,jwm,jz,jph,JHB
          DATA jwm,jwp,jz,jph,JHB /-1,1,2,3,4/        

      INTEGER MORDER
      COMMON /MCODES/ MORDER

      integer iset, norder
      COMMON/QCD1/ ISET,NORDER


CCPY===========================================
      INTEGER NIN1,NOUT1

      CHARACTER*40 QGFN, QTGFN, YGFN
      REAL*8 QG(100),PT(200),RAPY(200)
      COMMON / GRIDFILE / QG,PT,RAPY,QGFN,QTGFN,YGFN
      INTEGER II, N_Q,N_QT,N_Y, LTO
      COMMON / NGRID / N_Q,N_QT,N_Y, LTO

      REAL*8 ECM
      INTEGER IBEAM,JWTYPE_IN
      common/collider/ecm,ibeam
      
      Real*8 pT_Stretch,YMAX
      REAL*8 TANB,VEV,COUP_MB,CR,COUPL

      Integer KinCorr
      Common / CorrectKinematics / KinCorr
C=============================================

      REAL*8 ADYPT
      EXTERNAL ADYPT
      
      INTEGER I,J,K
      REAL*8 VPT,RAPIN,vmas
      INTEGER IQTMN, IQTMX, IQTST, IYMN, IYMX, IYST,
     &        IQMN, IQMX, IQST
      INTEGER IPTMIN, IPTMAX, IPTSTP, IYMIN, IYMAX, IYSTP,
     >        IQMIN, IQMAX, IQSTP

      REAL*8 asym_1,asym_2,factor_k,ASYM
CBY add asy part
      REAL*8 asymL0_1,asymA3_1,asymL0_2,asymA3_2
      real*8 scale_fac_muR,scale_fac_muF  
      integer JZ_TYPE

CJI Jan 2021: Add in LHAPDF
      integer initlhapdf, iset_save
      COMMON/LHAPDF/ initlhapdf, iset_save
      external LHAsub
   
C=============================================
cpn2007 For extrapolation of the K-factor outside of the kinematically
cpn2007 allowed region
      real*8 a_exterp, b_exterp, vpt_prev, factor_K_prev
!ZL
	character*100 jobname
	integer length
CBY Oct26,18
       INTEGER iscale

	if(iargc().eq.0)then
		jobname='w_asym'
	else if(iargc().ge.1) then
		call getarg(1,jobname)
	endif
	call trmstr(jobname,length)

      initlhapdf = 0
      
C I/O UNITS
      NIN1=23
      NOUT1=22
C
!ZL      OPEN(UNIT=NIN1,FILE='w_asym.in',STATUS='old')
!ZL      OPEN(UNIT=NOUT1,FILE='w_asym.out')
      OPEN(UNIT=NIN1,FILE=jobname(1:length)//'.in',STATUS='old')
      OPEN(UNIT=NOUT1,FILE=jobname(1:length)//'.out')

C            KinCorr -> 0: x_{1,2} =   Q/Sqrt[S] e^{+/-y} in Asy
C                       1: x_{1,2} = M_T/Sqrt[S] e^{+/-y} in Asy 

      READ(NIN1,*) ECM,IBEAM, KinCorr
      SS=ECM**2
      
C DATA jwm,jwp, jz,jph,JHB /-1,1,2,3,4/ 
      JZ_TYPE=0
      BOSON='-'
      READ(NIN1,*)JWTYPE_IN,JZ_TYPE
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
        GGD=0.D0
      ELSEIF(BOSON.EQ.'ZD') THEN
        GGU=0.D0
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
      
cpn2006
C      READ(NIN1,*) iset, norder
ccpy July 2013
      READ(NIN1,*) iset, norder, scale_fac_muR, scale_fac_muF, iscale
      
      READ(NIN1,'(A40)') PDF_FILE
      call TrmStr(PDF_FILE, LEN_PDF_FILE)
      PDF_FILE = PDF_FILE(1:LEN_PDF_FILE)

cpn2006
CJI Add in lhapdf
      if(PDF_FILE(1:3) .eq. "lha") then
          initlhapdf = 1
          call initpdfsetbyname(pdf_file(5:len_pdf_file)//".LHgrid")
          call initpdf(iset)
          iset_save = iset
      else if (iset.eq.902) then
        write (NOUT1,*) 'iset, pdf_file= ',iset,', ',pdf_file
      else
         write (NOUT1,*) 'iset = ', iset
       endif                    !iset
 
CCPY July 2013      WRITE(nout1,*) 'ECM, IBEAM, JWTYPE_IN, norder, KINCORR = ',
C     >ECM, IBEAM, JWTYPE_IN, norder,KINCORR 

      WRITE(NOUT1,*) 'ECM, IBEAM, JWTYPE_IN, NORDER, KINCORR, 
     >scale_fac_muR,scale_fac_muF= ',
     >ECM, IBEAM, JWTYPE_IN, NORDER, KINCORR, 
     >scale_fac_muR,scale_fac_muF
    
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
        DO II = 1,100
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

C       WRITE(nout1,*) 
C     >  '  Q, qT, y, asym_1, asym_2, factor_K=(asym_2+asym_1)/asym_1'
C        PRINT *, 
C     >   '  Q, qT, y, asym_1, asym_2, factor_K=(asym_2+asym_1)/asym_1'

CBY
       WRITE(nout1,*)
     >  '  Q, qT, y, asymL0_1, asymA3_1, asymL0_1+asymL0_2, 
     > asymA3_1+asymA3_2'
        PRINT *,
     >   '  Q, qT, y, asymL0_1, asymA3_1, asymL0_1+asymL0_2,
     > asymA3_1+ asymA3_2'

        pT_Stretch = 1.d0

        call setthepdf
        DO K = IQMIN, IQMAX, IQSTP
          VMAS = QG(K)
c===================
CZL
!          QM=VMAS
!          QQ=QM**2
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

              VPT=PT(I)*pT_Stretch
CBY Oct26,18
              if(iscale.eq.0) then
                QM=scale_fac_muF*VMAS
                QMu=scale_fac_muR*VMAS
              else
                QM=scale_fac_muF*sqrt(VMAS**2+VPT**2)
                QMu=scale_fac_muR*sqrt(VMAS**2+VPT**2)
              endif
                QQ=VMAS**2


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
C                  asym_1=ADYPT(vpt,rapin,1)
C                  asym_2=0.D0
                  factor_K=0.D0           
CBY,asyi=0 symmetric, asyi=1, asymmetric
                  asymL0_1=ADYPT(vpt,rapin,1,0)
                  asymA3_1=ADYPT(vpt,rapin,1,1)
                  asymL0_2=0.D0
                  asymA3_2=0.D0    
                ELSE
C                  asym_1=ADYPT(vpt,rapin,1)
C                  asym_2=ADYPT(vpt,rapin,2)
C                  factor_K=(asym_2+asym_1)/asym_1
CBY
                  asymL0_1=ADYPT(vpt,rapin,1,0)
                  asymA3_1=ADYPT(vpt,rapin,1,1)
                  asymL0_2=ADYPT(vpt,rapin,2,0)
                  asymA3_2=ADYPT(vpt,rapin,2,1)
                ENDIF

                asym_1=asym_1*2.0*VPT*COUPL
                asym_2=asym_2*2.0*VPT*COUPL

CBY
                asymL0_1=asymL0_1*2.0*VPT*COUPL 
                asymA3_1=asymA3_1*2.0*VPT*COUPL
                asymL0_2=asymL0_2*2.0*VPT*COUPL
                asymA3_2=asymA3_2*2.0*VPT*COUPL             
              else              !outside of the kinematical region
                asym_1=0d0; asym_2=0d0
                asymL0_1=0d0; asymA3_1=0d0
                asymL0_2=0d0; asymA3_2=0d0
CBY
                 factor_K=0d0
C                factor_K = a_exterp*vpt + b_exterp
              endif             !IF(ABS(RAPIN).LT.YMAX)
              
!ZL avoid too large k factor
		if(abs(factor_K).gt.5) then
			print *, 'warning: abs of k factor > 5 : ', 
     > factor_K, vmas,VPT,RAPIN
		endif
CZLDebug                        print *, "alpha_s(mz)=",ALPI(91.1876d0)*PI

C Output the results
C123           PRINT 102,    vmas,VPT,RAPIN,ASYM_1,ASYM_2+ASYM_1,factor_K 
C              WRITE(NOUT1,102) vmas,VPT,RAPIN,ASYM_1,ASYM_2,factor_K

CBY new ouput
123           PRINT 102,vmas,VPT,RAPIN,ASYML0_1,ASYMA3_1,
     > ASYML0_1+ASYML0_2,ASYMA3_1+ASYMA3_2
              WRITE(NOUT1,102) vmas,VPT,RAPIN,ASYML0_1,ASYMA3_1,
     > ASYML0_1+ASYML0_2,ASYMA3_1+ASYMA3_2
cpn2007 Compute coefficients of linear extrapolation for the K factor
              a_exterp=(factor_K-factor_K_prev)/(vpt - vpt_prev)
              b_exterp=(factor_K_prev*vpt-factor_K*vpt_prev)/
     >          (vpt - vpt_prev)
              vpt_prev = vpt
              factor_K_prev=factor_K

            END DO    ! I (QT)
          ENDDO       ! J (Y)
        ENDDO         ! K (Q)
        CLOSE(NOUT1)
        RETURN
cpn2007 Keep enough digits in asym_1, asym_2, factor_K in order not to
cpn2007 lose precision in the Y-piece K-factor despite large cancellations

  102 FORMAT(1X,3(G10.4,2x),4(G14.7,2x))


      
      END                       !w_asym
      
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
CJI 2014: To handle high pt region
      IF(ISNAN(YMAX)) THEN
          YMAX=-9999
      ENDIF

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

      if (iset.eq.902)
     >  open(npdffile,file=pdf_file,status='unknown')

c Initialize PDFs
      Call SetQCD
      Call SetUtl

      close(npdffile)

      Call SetPDF(iset)

      return 
      end 
