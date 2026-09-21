!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%.%%|%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%                                  %%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%              ResBos              %%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%                                  %%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%         Version 10.02.09cpy      %%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%         Version 06.09.08pn       %%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%                                  %%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
C----------------------------------------------------------------------
C- Common block to be used with the amplitudes calculated with the 
C- help of spinor-helicity methods. Can be used in the calculations
C- of the amplitudes with up to 10 legs.
C-
C-   Created   31-MAR-1998 by  Pavel Nadolsky, Michigan State University
C-
C----------------------------------------------------------------------
cpn 2004 Last modified: Jan. 16, 2004

      DOUBLE PRECISION HIGGS_MASS,
     &ZUP,ZDOWN,XMGW2,XMGZ2,XMGH2,XMH2,VSQ,
     &XMW,XMZ,XMW2,XMZ2,AEM,GAMW,GAMZ,GAMW2,GAMZ2,
     &SWS,CWS,XMBOT,XMTOP,XMB2,XMT2,GWEAK2,
     &XLU,XLD,XRU,XRD,QU,QD,GFERMI,XLE,XRE,
     &XLU2,XLD2,XRU2,XRD2,QU2,QD2,XLE2,XRE2,
     &XLU4,XLD4,XRU4,XRD4,QU4,QD4,PI,
     &XMCHARM,XMCHARM2,XMHIGGS,GAMH,GAMH2,
     &CURRENTW_MASS,CURRENTZ_MASS,
     &GAMT
      DOUBLE PRECISION VKM
CsB      DIMENSION VKM(6,6)

      COMMON/CKM/VKM(6,6)
      COMMON/HIGGS/HIGGS_MASS
      COMMON/HIGGS_STUFF/ZUP,ZDOWN,XMGW2,XMGZ2,XMGH2,XMH2,VSQ
      COMMON/STAND/XMW,XMZ,XMW2,XMZ2,AEM,GAMW,GAMZ,GAMW2,GAMZ2
      COMMON/STAND1/SWS,CWS,XMBOT,XMTOP,XMB2,XMT2,GWEAK2
      COMMON/STAND2/XLU,XLD,XRU,XRD,QU,QD,GFERMI,XLE,XRE
      COMMON/STAND3/XLU2,XLD2,XRU2,XRD2,QU2,QD2,XLE2,XRE2
      COMMON/STAND4/XLU4,XLD4,XRU4,XRD4,QU4,QD4,PI
      COMMON/STAND5/XMCHARM,XMCHARM2,XMHIGGS,GAMH,GAMH2
      COMMON/RUNNING_BOSON/CURRENTW_MASS,CURRENTZ_MASS
      COMMON/HTOP/GAMT

      COMMON/Tens/Sq(10,10),Tr(10,10),S(10,10),pp(10,3),h(10)
      double precision S        ! Relativistic invariants
      double precision pp       ! +,- components and polar angle of the particle momenta
      integer h        ! Helicities of the external particles
      double complex Sq,Tr ! Square & Dirac (triangular) spinor products


      COMMON /STD/const,cg,inout,Nleg,Nc,MomFlag
      COMMON /misc/ie,zero,epsil,pi2, zeta3, idebug,iwrite
      common /scale/Q
      Common / BuildInPDF / BIPDF, ISet
      common /pat_scale/itscale

      integer inout(10),Nleg,Nc,MomFlag,idebug,iwrite
      double complex ie
      double precision zero,const,cg,epsil, pi2, zeta3
      double precision Q
      Integer BIPDF, ISet
      integer itscale

cpn 2004 For h1 in g g -> AA
      integer Naux
      parameter (Naux = 12)
      double precision vvv1, vvv2
      common /vvvs/ vvv1(Naux), vvv2(Naux)

      double precision v, v2, dlogv, dlogv2, onepv, dlogv1mv, v1mv,
     >  fLi2v, fLi3v, fLi4v, fLi2v1mv, fLi3v1mv, fLi4v1mv, 
     >  onemv, onemv2, dlog1mv, dlog1mv2, twomv, dlog1mvv, onemvv, 
     >  fLi21mv, fLi31mv, fLi41mv,fLi21mvv, fLi31mvv, fLi41mvv 
      common/vs/ v, v2, dlogv, dlogv2, onepv, dlogv1mv, v1mv, 
     >  fLi2v, fLi3v, fLi4v, fLi2v1mv, fLi3v1mv, fLi4v1mv, 
     >  onemv, onemv2,dlog1mv, dlog1mv2, twomv, dlog1mvv, onemvv, 
     >  fLi21mv, fLi31mv, fLi41mv,fLi21mvv, fLi31mvv, fLi41mvv 

cpn Relative error of the integration in the polylogarithms
      double precision rlierr
      common/lierr/rlierr


cpn   for debugging
c      double precision tmpr, tmpi
c      double complex ctmp
c      integer id
