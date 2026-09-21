!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%.%%|%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%                                  %%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%              ResBos              %%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%                                  %%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%         Version 10.02.09cpy      %%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%         Version 06.09.08pn       %%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%                                  %%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
C----------------------------------------------------------------------------
C
C                                      *
C                                     *
C                                      *
C                                     *
C                                     PP
C                                  PP    PP
C                                PP        PP
C                              PP            PP
C                      *******PP   PENTAGON   PP*******
C                              PP            PP
C                               PP          PP
C                                PP        PP
C                                  PPPPPPPP
C                                 *        *
C                                *          *
C                               *            *
C                              *              *
C
C----------------------------------------------------------------------------
C-
C-   Purpose and Methods :Calculation of the squared matrix element for
C-   processes q(1)+qbar(2)+gamma(3)+gamma(4)+g(5) -> 0;
c              g(1)+g(2) -> gamma(3) + gamma(4) + g(5) at 1-loop level
c              q(1)+g(2) -> gamma(3)+gamma(4)+q(5) via the quark box  
C-
C-   Inputs  : PENTAGON.F  -- the main part of the code
C-             PNCOMMNON.F -- the common block
C-   Outputs :
C-   Controls:
C-
C-   Written by Pavel Nadolsky, Argonne National Laboratory
C-   Version November 28, 2006
C-

C==========================================================================
CPN Amplitude squared for  q Q or q G -> diphoton + parton at the tree level
C========================================================================
c
C -------------------------------------------------------------------------
      SUBROUTINE QQAAG(WT)
C -------------------------------------------------------------------------
      IMPLICIT NONE
      include 'pncommon.f'
      double precision pyalem,alphas, g_up4, g_dn4, Amp2, WT, WT1, WT2,
     >  WT3, Amp2Sub, WT1Sub, WT2Sub, WT3Sub, CF
      double precision qQBUL,qQBDL,QBqUL,QBqDL, qGUL,qGDL,GqUL,
     >  GqDL,GGL, Amp2qqVVV, qqVVVSubtractionIS, qqVVVSubtractionFS
      Integer iProc
      Common / PartonProcess / iProc
      double precision  Q_V,QT_V,YB_V,THE_STA,PHI_STA
      COMMON/CSFRAME/Q_V,QT_V,YB_V,THE_STA,PHI_STA
      double precision  QT_V_SEP
      COMMON/NLO/ QT_V_SEP

      integer isub
      double precision cutM45, SubtractionShape
      logical SubQGFlag
      Common / Subtraction / cutM45, SubtractionShape(3:4),iSub,
     >  SubQGFlag(3:4)

cpn2006 Set constants and overall normalization factors
      Pi= 3.141592653589793238d0
      nleg=5  !number of legs in the Feynman diagram
      Nc=3       !Number of colors
      CF = (Nc**2-1d0)/2/Nc

cpn2006 Overall normalization const in the qqbar channel already includes 
c      averaging factors for spin (1/4), color (1/Nc**2), and 
c      identical final-state particles (1/2)
c      const = 7937.6068301567539   !2*(4pi)**3
      const = 3968.803415 !2*(4pi)**3
      const= const*(CF/Nc)*alphas(Q)*(pyalem(Q**2))**2

cpn2006 Fermion-photon couplings 
      g_up4 = (2d0/3d0)**4
      g_dn4 = (-1d0/3d0)**4
cpn2006 Compute parton luminosities

      call PartonMomentumLumi(qQBUL,qQBDL,QBqUL,QBqDL,
     >  qGUL,qGDL,GqUL,GqDL,GGL)

      Amp2=0d0; Amp2Sub=0d0
      WT1=0d0; WT2=0d0; WT3=0d0; WT1Sub=0d0; WT2Sub=0d0; WT3Sub=0d0

Cpn The qqbar squared matrix element (following the
cpn2006 M2q-q-gamma-gamma-g notes)

cpn2006 Amp2qqVVV(i1,i2,i3,i4,i5) returns the squared matrix element for 
cpn     q(i1)+qB(i2)+gamma(i3)+gamma(i4)+g(i5) ->0
cpn     Note: In the diphoton case, Amp2qqVVV is symmetric with respect 
c       to the interchanges i1<->i2, i3 <-> i4
      Amp2 = Amp2qqVVV(1,2,3,4,5)
c      if (qt_v.lt.qt_v_sep)
c     >  Amp2Sub=qqVVVSubtractionIS()

      WT1 = ( (qQBUL+QBqUL)*g_up4 + (qQBDL+QBqDL)*g_dn4 )*(Amp2-Amp2Sub)
      WT1Sub=( (qQBUL+QBqUL)*g_up4 + (qQBDL+QBqDL)*g_dn4 )*Amp2Sub
Cpn The qG + Gq, qB G + G qB squared matrix elements

cpn2006 Obtain the overall normalization in the qg channel by adjusting
cpn2006 the color factor in the normalization for the qqbar channel
      const = const*Nc/(Nc**2-1)

Cpn q G + G q matrix amplitude squared and its final-state collinear subtraction
C   For the q G -> photon photon q process we define a subtraction term which
C   subtracts the collinear singular term when one of the final state photon is
C   emitted from the final state quark leg. 
C___iSub=0 no subtraction is done
C    iSub=1 a collinear piece is subtracted from the qG->AAq amplitude
C            unconditionally
      Amp2 = Amp2qqVVV(1,5,3,4,2)

      if (SubQGFlag(3).or.SubQGFlag(4))
     >  Amp2Sub = qqVVVSubtractionFS(1)

CPN   Note: qGUL and qGDL include both q G and G q contributions
      WT2 = ( qGUL*g_up4 + qGDL*g_dn4)*(Amp2-Amp2Sub) 
      WT2Sub = ( qGUL*g_up4 + qGDL*g_dn4)*Amp2Sub

Cpn qB G + G qB matrix amplitude squared and its collinear subtraction
C   For the q G -> photon photon q process we define a subtraction term which
C   subtracts the collinear singular term when one of the final state photon is
C   emitted from the final state quark leg. 
C___iSub=0 no subtraction is done
C    iSub=1 a collinear piece is subtracted from the qG->AAq amplitude
C            unconditionally
      Amp2 = Amp2qqVVV(5,2,3,4,1)

      if (SubQGFlag(3).or.SubQGFlag(4)) 
     >  Amp2Sub = qqVVVSubtractionFS(2)

CPN   Note: QBGUL and QBGDL include both qB G and G qB contributions
      WT3 = ( GqUL*g_up4 + GqDL*g_dn4)*(Amp2-Amp2Sub) 
      WT3Sub = ( GqUL*g_up4 + GqDL*g_dn4)*Amp2Sub

      If (iProc.Eq.0) then
        WT = WT1 + WT2 + WT3   ! Both Contributions
      ElseIf (iProc.Eq.1) then
        WT = WT1               ! q*qB + qB*q only
      ElseIf (iProc.Eq.2) then
        WT = WT2 + WT3         ! (q+qB)*G + G*(q+qB)  only
      ElseIf (iProc.Eq.12) then
        WT = WT2               ! (q+qB)*G only
      elseif (iproc.eq.13) then   
        WT = WT3               ! G*(q+qB) only  
      elseif (iproc.eq.-1) then 
        WT = WT1Sub            ! q*qB               (IS subtraction)
      elseif (iproc.eq.-2) then 
        WT = WT2Sub + WT3Sub   ! (q+qB)*G + G*(q+qB) (FS subtraction)
      elseif (iproc.eq.-12) then 
        WT = WT2Sub            ! (q+qB)*G (FS subtraction)
      elseif (iproc.eq.-13) then  
        WT = WT3Sub            ! G*(q+qB) (FS subtraction)
      Else
        Print*, ' iProc = ', iProc, ' value in not defined. '
        Stop
      EndIf

      RETURN
      END !QQAAG

cpn2006 ==============================================================

      double precision function Amp2qqVVV(i1, i2, i3, i4, i5)
cpn2006 Reduced q(1)q(2)V(3)V(4)V(5) amplitude squared, given for all
cpn2006 incoming massless particles 
cpn2006 Inputs: i1, i2, i3, i4, i5 -- indices of the external particle
cpn2006         momenta
      implicit none  
      include 'pncommon.f'
      double precision plab, pcem, what
      COMMON/MOM/PLAB(4,12),PCEM(4,12),WHAT
      integer i1, i2, i3, i4, i5, j, k
      
cpn2006 Scalar products 2*p_i*p_j (can differ by a sign from Mandelstam 
cpn2006 invariants (p_i+-p_j)^2 !)
      do j=1, nleg
        do k=1, nleg
          S(j,k) = 2d0*(plab(4,j)*plab(4,k) - 
     >    plab(3,j)*plab(3,k) -plab(2,j)*plab(2,k) -plab(1,j)*plab(1,k))  
        enddo 
      enddo 
      
      Amp2qqVVV =  const*(
     >    S(i1,i3)**3*S(i2,i3) + S(i1,i3)*S(i2,i3)**3 + 
     >    S(i1,i4)**3*S(i2,i4) + S(i1,i4)*S(i2,i4)**3 + 
     >    S(i1,i5)**3*S(i2,i5)+  S(i1,i5)*S(i2,i5)**3  ) *
     > S(i1,i2)/S(i1,i3)/S(i1,i4)/S(i1,i5)/S(i2,i3)/S(i2,i4)/S(i2,i5) 

      return 
      end !Amp2qqVV

c======================================================================

      double precision function qqVVVSubtractionIS()

cpn2006 Catani-Seymour initial-state (IS) collinear subtraction for 
c       the reduced q(1)q(2)V(3)V(4)V(5) amplitude squared. The 
c       phase-space integration of the 2->3 matrix element satisfies
c       the constraint (1/x1hat-1)*(1/x2hat-1)=qT^2/MT^2, corresponding
c       to a hyperbola in a plane of integration variables x1 and x2. 
c       The exact limits for collinear radiation off the initial-state
c       leg 1 and 2 correspond to x2hat=1 and x2hat=1, respectively.
c       We apply collinear subtraction for the leg 1 (leg 2) to the
c       lower (upper) branch of the hyperbola, defined by x2hat > x1hat
c       (x2hat < x1hat). The soft singularity, arising at x1hat=x2hat=1,
c       is subtracted implicitly as a part of the two collinear subtractions.
c
      implicit none  
      include 'pncommon.f'
      double precision x1hat, x2hat, Amp2QQBVV, Sub, xsplit
      COMMON/x1x2hat/X1HAT,x2hat
      double precision  Q_V,QT_V,YB_V,THE_STA,PHI_STA
      COMMON/CSFRAME/Q_V,QT_V,YB_V,THE_STA,PHI_STA

      xsplit=min(x1hat,x2hat)

cpn2006 Squared amplitude for Q+QBAR -> gamma + gamma
      Amp2QQBVV = const*(S(1,3)/S(1,4)+ S(1,4)/S(1,3))

      Sub=(xsplit**2+1d0)/xsplit*Amp2QQBVV/qt_v**2

      qqVVVSubtractionIS = Sub

      return 
      end ! qqVVVSubtractionIS

c======================================================================

      double precision function qqVVVSubtractionFS(isq)

cpn2006 Catani-Seymour final-state (FS) subtraction for 
c       the reduced q(1)q(2)V(3)V(4)V(5) amplitude squared
c Inputs: isq=1, 2 -- ID of the initial-state quark leg
      implicit none  
      include 'pncommon.f'
      integer ifragQ, isq, ifragA, ispecA, isub
      double precision cutM45, SubtractionShape, zCS, Sub, Sub45, 
     >  PgammaQ, Amp2QGVQ
      logical SubQGFlag
      Common / Subtraction / cutM45, SubtractionShape(3:4),iSub,
     >  SubQGFlag(3:4)

cpn2006 Catani-Seymour's momentum fraction z
      zCS()= 1d0 - S(5, ispecA)/
     >    (S(ifragA, ispecA)+S(5, ispecA) + S(5, ifragA))
c      zCS()= S(ifragA, ispecA)/                 !Alternative definition of z
c     >    (S(ifragA, ispecA)+S(5, ispecA))

cpn2006 Quark->photon splitting function
      PgammaQ()=(1d0 + (1d0-zCS())**2)/zCS()
cpn2006 Squared amplitude for Q+G -> gamma + Q
      Amp2QGVQ() = const*(S(1,2)/S(ispecA,isq)+ S(ispecA,isq)/S(1,2))

      if (isq.ne.1.and.isq.ne.2) then
        print *,
     >  'qqVVVSubtractionFS: wrong ID = ',
     >    isq, ' of the initial-state quark (must be 1 or 2)'
        call ExitGracefully
      endif                      

      Sub=0d0

      do ifragA=3,4 !conditionally perform subtraction in the vicinity 
                    !of the fragmentation photon ifragA (ifragA=3 or 4)

        if (SubQGFlag(ifragA)) then
          ispecA=7 - ifragA     !spectator photon ID (4 for ifragA=3, 
                                !                     3 for ifragA=4)
          Sub = Sub + SubtractionShape(ifragA)*PgammaQ()*
     >      Amp2QGVQ()/S(ifragA,5)
        endif                   !SubQGFlag
      enddo                     !ifragA

      qqVVVSubtractionFS = Sub

      return 
      end ! qqVVVSubtractionFS



C========================================================================
C Amplitude squared for G+G->diphoton+G via the box diagram
C========================================================================
C
C----------------------------------------------------------------------------
      double precision function PentagonGGAAG (i1,i2,i3,i4,i5)
C-   Returns the squared matrix element for
C-   the process g(i1) +g(i2) -> gamma(i3) + gamma(i4) + g(i5) at 1-loop level
C----------------------------------------------------------------------------
      IMPLICIT NONE
      include 'pncommon.f'
      integer i,j,iret,isum,sgn, i1, i2, i3, i4, i5
      double precision AmpSqdGG,xi(4)
      character*8 time1
      character*9 date1
      real*4 t1,t2                       ! Execution time controls
CsB___Info from ResBos to Pentagon
      DOUBLE PRECISION W,QCDL,X1X2MIN
      INTEGER*4 NJET,N_PRIMARY_PART,NF,I_STR,IBEAM,NPART,NEVENT,NORDER
      COMMON/AA/W,QCDL,X1X2MIN,NJET,N_PRIMARY_PART,NF,I_STR,IBEAM,NPART,
     &  NEVENT,NORDER
C      Integer BIPDF, ISet               ! already in PNCOMMON
C      Common / BuildInPDF / BIPDF, ISet ! already in PNCOMMON
      DOUBLE PRECISION PLAB,PCEM,WHAT
      COMMON/MOM/PLAB(4,12),PCEM(4,12),WHAT
CsB___Locally defined
      double precision tmp
CsB___Debug options
      Logical Debug
      Pi= 3.141592653589793238d0

CsB___Pentagon inputs
C      OPEN(UNIT=13,FORM='FORMATTED',file='pentagon.in',
C     >        STATUS='OLD',READONLY)
C      read(13,*) idebug,iwrite
C      READ(13,*) W,xi(1),xi(2),xi(3),xi(4),epsil,sgn
C      read(13,*) BiPDF,itscale
C      read(13,*) isum,h(1),h(2),h(3),h(4),h(5)      
      iDebug = 0      ! no print from Pentagon
      iWrite = 0      ! no write from Pentagon
      epsil  = 1.e-10 ! numeric zero
cpn2005 To produce a specific helicity amplitude
      iSum   = 1      ! sum over helicities

c$$$      do i=1,5
c$$$        h(i) = 1
c$$$      enddo

      h(1)=2; h(2)=2; h(3)=1; h(4)=1; h(5)=1

C_____The definitions of the four momenta which Pentagon uses in terms of the
C     four momenta which are used in ResBos (cf. Pavel's notes 3G2A)
C     Note:             in ResBos        in Pentagon
C           gluon 1 in     1                 4
C           gluon 2 in     2                 5
C           photon 1 out   3                 1
C           photon 2 out   4                 2
C           gluon out      5                 3
C      a = particle index
C      b = coordinate index
C        a b         b a
      pp(4,1) = pLab(4,i1) + pLab(3,i1)       ! p+  of gluon 4 in
      pp(4,2) = pLab(4,i1) - pLab(3,i1)       ! p-  of gluon 4 in
      pp(4,3) = 0.d0
      pp(5,1) = pLab(4,i2) + pLab(3,i2)       ! p+  of gluon 5 in
      pp(5,2) = pLab(4,i2) - pLab(3,i2)       ! p-  of gluon 5 in
      pp(5,3) = 0.d0
      pp(1,1) = pLab(4,i3) + pLab(3,i3)       ! p+  of photon 1
      pp(1,2) = pLab(4,i3) - pLab(3,i3)       ! p-  of photon 1
      pp(1,3) = ASin(pLab(2,i3)/Sqrt(pLab(1,i3)**2+pLab(2,i3)**2))
      If (pLab(1,i3).lt.0.0) pp(1,3) = Pi - pp(1,3)

cpn Alternative definition of the angle
c      pp(1,3) = ACos(pLab(1,i3)/Sqrt(pLab(1,i3)**2+pLab(2,i3)**2))
c      pp(1,3)=pp(1,3)*sign(1,pLab(2,i3))

      pp(2,1) = pLab(4,i4) + pLab(3,i4)       ! p+  of photon 2
      pp(2,2) = pLab(4,i4) - pLab(3,i4)       ! p-  of photon 2
      pp(2,3) = ASin(pLab(2,i4)/Sqrt(pLab(1,i4)**2+pLab(2,i4)**2))
      If (pLab(1,i4).lt.0.0) pp(2,3) = Pi - pp(2,3)
cpn 
c      pp(2,3) = ACos(pLab(1,i4)/Sqrt(pLab(1,i4)**2+pLab(2,i4)**2))
c      pp(2,3)=pp(2,3)*sign(1,pLab(2,i4))

      pp(3,1) = pLab(4,i5) + pLab(3,i5)       ! p+  of gluon out
      pp(3,2) = pLab(4,i5) - pLab(3,i5)       ! p-  of gluon out
      pp(3,3) = ASin(pLab(2,i5)/Sqrt(pLab(1,i5)**2+pLab(2,i5)**2))
      If (pLab(1,i5).lt.0.0) pp(3,3) = Pi - pp(3,3)
cpn
c      pp(3,3) = ACos(pLab(1,i5)/Sqrt(pLab(1,i5)**2+pLab(2,i5)**2))
c        pp(3,3)=pp(3,3)*sign(1,pLab(2,i5))


CsB___Output debug info headers
      if (iwrite.eq.1) then
        open(unit=14,file='ah.dat',status='new')
        open(unit=15,file='cop.dat',status='new')
        open(unit=17,file='ahcop.dat',status='new')
        call time(time1)
C        call date(date1)
C        write (14,*) 'AH  Created:',time1,'   ',date1
C        write (15,*) 'COP  Created:',time1,'   ',date1
C        write (17,*) 'AHCOP     Created:',time1,'   ',date1
      endif

CsB___Execute Pentagon subroutines
CPN xi2mom is only needed if pentagon is called standalone
C      Call xi2mom(W,xi,sgn)

      Call Standard          ! Get rid of this

      Call ConstInit1(1)     ! This part can go out of the loop 

      Call KinemInit      !Initialize kinematics 

      Call ConstInit2     !Set up the couplings

      Call TensorInit     !Fill in spinor products 

      PentagonGGAAG=AmpSqdGG(isum,iret)

C_____Output to the screen
      if (idebug.eq.1) then
        print *,'isum=',isum
        if (isum.eq.1) print *,'hel=',
     &    h(1)*10000+h(2)*1000+h(3)*100+h(4)*10+h(5)
        print *,'AmpSqdGG=',PentagonGGAAG
      endif
!     Print *,'AmpSqdGG=',PentagonGGAAG

CsB___Clean up
cpn      close(13)
      if (iwrite.eq.1) then
       close(14)
       close(15)
       close(17)
      endif

      RETURN
      END ! PentagonGGAAG

C----------------------------------------------------------------------------
      subroutine xi2mom(W,xi,sgn)
C----------------------------------------------------------------------------
C     Calculates all the particle momenta from the given values of
C     W=sqrt(s), p_{1+-},p_{2+-},sign(sin phi(p_1,p2))
CsB___NOT used in ResBos, only for testing purposes.

      implicit NONE
      include 'pncommon.f'
      integer i,j,sgn
      double precision W,xi(4),pi1
      double precision cosphi,sinphi

      pi1=3.14159265358979d0

      do i=1,5
       do j=1,3
        pp(i,j)=0d0
       enddo
      enddo

      pp(4,1)=W
      pp(5,2)=W

      pp(1,1)=W*xi(1)
      pp(1,2)=W*xi(2)

      pp(2,1)=W*xi(3)
      pp(2,2)=W*xi(4)

      cosphi=(1-xi(1)-xi(2)-xi(3)-xi(4)+xi(1)*xi(4)+xi(2)*xi(3))/
     &        2d0/dsqrt(xi(1)*xi(2)*xi(3)*xi(4))

      if (dabs(cosphi).gt.1) then
        print *,'STOP:cos(phi)=',cosphi
        STOP
      endif

      pp(2,3)=dacos(cosphi)+(1-sgn)/2d0*pi1   !+pi if sin(phi).lt.0
      sinphi=sgn*dsqrt(1d0-cosphi**2)

      pp(3,1)=W-pp(1,1)-pp(2,1)
      pp(3,2)=W-pp(1,2)-pp(2,2)

      pp(3,3)=dacos(-(dsqrt(pp(1,1)*pp(1,2))+dsqrt(pp(2,1)*pp(2,2))*
     &              cosphi)/dsqrt(pp(3,1)*pp(3,2)))
      if (sinphi.gt.0) pp(3,3)=2d0*pi1-pp(3,3)
      
      do 10 i=1,5
       do 10 j=1,2
       if (pp(i,j).lt.0) then
        print *,'STOP: pp(',i,',',j,')=',pp(i,j),' < 0'
        STOP
       endif
   10 continue

c$$$       if (idebug.eq.1) then
c$$$        print *,'pp(i,1:3):'
c$$$         do i=1,5
c$$$          print 20, i,': ',pp(i,1),pp(i,2),pp(i,3)
c$$$         enddo
c$$$       endif
 20    format (i2,':',3(g16.8,1x))

      End ! xi2mom

      
C---------------------------------------------------------------------
      subroutine ConstInit2
C---------------------------------------------------------------------
      implicit NONE
      include 'pncommon.f'
      double precision pyalem,alphas
cpn      external pyalem, alphas 
CsB >
!     Real*8 const1
CsB <
    
      const=const*(pyalem(Q**2)**2)*(alphas(Q)**3.0)
                              ! color averaging factor for 2 initial
                              ! gluons included
c$$$      if (idebug.eq.1) print *,'ConstInit2:Q,s(1,2),const='
c$$$     >       ,Q,s(1,2),const
      
C  Factor obtained by the summation over the quark flavors
CsB >
!     const1=const
CsB <

cpn2005 This is not correct
c$$$      if (Q.lt.xmcharm) then
c$$$       const=const*(6./9.)**2
c$$$      elseif(Q.lt.xmbot) then
c$$$       const=const*(10./9.)**2
c$$$      elseif (Q.lt.xmtop) then
c$$$        const=const*(11./9.)**2
c$$$      else 
c$$$        const=const*(15./9.)**2
c$$$      endif

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

      const=const*(11./9.)**2



CsB >
!     Print*, 'ConstInit: Q,const=',Q,const/const1
!     Print*, 'ConstInit: ',(pyalem(Q**2)**2)
!     Print*, 'ConstInit: ',(alphas(Q)**3)
CsB <

      End !ConstInit2

C----------------------------------------------------------------------------
      subroutine KinemInit
C----------------------------------------------------------------------------
C Calculate either the particle momenta from the relativistic invariants
C or vice versa, depending on the input mode
      implicit NONE
      include 'pncommon.f'
      integer i,j
      double precision sum1,sum2


      if (MomFlag.eq.1) then
C particle momenta ----> relativistic invariants
        do i=1,nleg
          S(i,i)=0
          do j=i+1,nleg
            S(i,j)=pp(i,1)*pp(j,2)+pp(i,2)*pp(j,1)-2*dsqrt(pp(i,1)*
     &        pp(i,2)*pp(j,1)*pp(j,2))*dcos(pp(i,3)-pp(j,3))

C Change the sign of S(i,j) in case one of the particles is ingoing and
C the other is outgoing

            S(i,j)=S(i,j)*inout(i)*inout(j)
          enddo
          do j=1,i-1
            S(i,j)=S(j,i)
          enddo
        enddo

      else
C invariants ---> particle momenta

        print *,'STOP: No relativistic invariants given'
        STOP
      endif

      if (idebug.eq.1) then
       do i=1,nleg
        sum1=0d0
        sum2=0d0
        do j=1,nleg
         sum1=sum1+s(i,j)
         sum2=sum2+s(j,i) 
        enddo

        if (dabs(sum1/s(4,5)).gt.1d-9) then
         print *,'Sum of kinematic invariants s(i,j) neq 0, i=',i
        endif

        if (dabs(sum2/s(4,5)).gt.1d-9) then
         print *,'Sum of kinematic invariants s(i,j) neq 0, j=',i
        endif
      enddo
     
      endif

      End ! KinemInit

C----------------------------------------------------------------------------
      subroutine TensorInit
C----------------------------------------------------------------------------
C     This subroutine calculates the square and Dirac bracket products
C     to be used in the primitive amplitude

      implicit NONE
      include 'pncommon.f'
      integer i,j,k
      double complex sum1


C  Calculation of the spinor products from the particle momenta
C  Conventions for the i-th partcile momenta are:
C    p_plus(i)=pp(i,1),p_minus(i)=pp(i,2), the polar angle
C    phi(i)=pp(i,3)

      do 20 i=1,nleg
       Sq(i,i)=0.0
       Tr(i,i)=0.0

       do 30 j=i+1,nleg
        Sq(i,j)=dsqrt(pp(i,1)*pp(j,2))*cdexp(ie*pp(j,3))-
     &          dsqrt(pp(i,2)*pp(j,1))*cdexp(ie*pp(i,3))
        Tr(i,j)=-dsign(1d0,S(i,j))*dconjg(Sq(i,j))
  30   continue

       do 40 j=1,i-1
        Sq(i,j)=-Sq(j,i)
        Tr(i,j)=-Tr(j,i)
  40   continue

  20  continue

      if (idebug.eq.1) then
       do i=1,nleg
        do j=1,nleg
         sum1=(0d0,0d0) 
         do k=1,nleg
           sum1=sum1+sq(i,k)*tr(k,j)
         enddo    

         if (cdabs(sum1).gt.1d-9) then
          print *,'Sum Sq(',i,',k)*(k,',j,')=',sum1,' neq 0'
         endif 
        enddo
       enddo 
    
      endif 

      End !TensorInit

C----------------------------------------------------------------------------
      function AmpSqdGG(isum,iret)
C----------------------------------------------------------------------------

      implicit NONE
      include 'pncommon.f'

      double precision AmpSqdGG
      integer k1,k2,k3,k4,k5
      integer iret,isum
      double complex Out,COP,AHCOP,AH,tmp1,tmp2,CAH
      external AHCOP,AH,CAH
CsB >
      Complex gggAAA
      Common / gggAAAmplitude / gggAAA(2,2,2,2,2,2)
      Data Pi / 3.1415926535897932 /
CsB <

      iret=0

      Out=(0.,0.)

      if (isum.eq.0) then
C no summation over helicities

cpn2005        print*,'Helicity summation disabled'

C        tmp1=COP(1,2,3,4,5,AHCOP)
C        tmp2=COP(1,2,3,5,4,AHCOP)
C
C        Out=tmp1+tmp2

      
        tmp1=COP(1,2,3,4,5,AH)
c        tmp2=COP(1,2,3,5,4,CAH) 
c        Out=-2d0*tmp1*tmp2
        Out=2d0*tmp1*dconjg(tmp1)

    
        if (idebug.eq.1) print 15,h(1),h(2),h(3),h(4),h(5),Out
   15   format (1x,'hel=',5i1,'  Out=',2g16.8)

  
      elseif (isum.eq.1) then
CsB >
        do  10 k1=1,2
!1      do  10 k1=1,1
         h(1)=k1
         do  10 k2=1,2
          h(2)=k2
          do  10 k3=1,2
           h(3)=k3
           do  10 k4=1,2
            h(4)=k4
            do  10 k5=1,2
             h(5)=k5

cpn             Out=Out-2d0*COP(1,2,3,4,5,AH)*COP(1,2,3,5,4,CAH)
             tmp1=COP(1,2,3,4,5,AH)                 

             gggAAA(k1,k2,k4,k5,k3,1) = tmp1

             Out=Out+2d0*tmp1*dconjg(tmp1) 
             
  10    continue
!1      Out = 2d0*Out
CsB <
      else
       print *,'isum=',isum,' not implemented yet!'
       STOP
      endif !isum

      if (dabs(DIMAG(Out)/DReal(Out)).ge.1d-6) then
        print *,
     &   'WARNING: The squared amplitude has a non-zero imaginary part'
        print *,'Im(AmpSqdGG)=',dimag(Out)*const
        iret=-1
      endif

      if (idebug.eq.1) then
       print *,'AmpSqdGG: const=',const
      endif
  
      AmpSqdGG=Dreal(Out)*const

      End !ampsqdGG

C----------------------------------------------------------------------------
      function AH(i1,i2,i3,i4,i5)
C----------------------------------------------------------------------------
C      Creates an amplitude for any set of helicities
C      using 8 primitive amplitudes AA1-AA8

      implicit NONE
      include 'pncommon.f'

      integer i1,i2,i3,i4,i5
      integer*4 hel
      double complex AH,AA1,AA2,AA3,AA4,AA5,AA6,AA7,AA8,temp

      hel=h(i5)+h(i4)*10+h(i3)*100+h(i2)*1000+h(i1)*10000

      if (hel.eq.11111) then
       TEMP=AA1(i1,i2,i3,i4,i5)       !1
      elseif (hel.eq.11112) then
       TEMP=AA2(i5,i1,i2,i3,i4)
      elseif (hel.eq.11121) then
       TEMP=AA2(i4,i5,i1,i2,i3)
      elseif (hel.eq.11122) then
       TEMP=AA3(i4,i5,i1,i2,i3)
      elseif (hel.eq.11211) then
       TEMP=AA2(i3,i4,i5,i1,i2)       !5
      elseif (hel.eq.11212) then
       TEMP=AA4(i3,i4,i5,i1,i2)
      elseif (hel.eq.11221) then
       TEMP=AA3(i3,i4,i5,i1,i2)
      elseif (hel.eq.11222) then
       TEMP=AA7(i1,i2,i3,i4,i5)
      elseif (hel.eq.12111) then
       TEMP=AA2(i2,i3,i4,i5,i1)
      elseif (hel.eq.12112) then
       TEMP=AA4(i5,i1,i2,i3,i4)       !10
      elseif (hel.eq.12121) then
       TEMP=AA4(i2,i3,i4,i5,i1)
      elseif (hel.eq.12122) then
       TEMP=AA8(i1,i2,i3,i4,i5)
      elseif (hel.eq.12211) then
       TEMP=AA3(i2,i3,i4,i5,i1)
      elseif (hel.eq.12212) then
       TEMP=AA8(i4,i5,i1,i2,i3)
      elseif (hel.eq.12221) then
       TEMP=AA7(i5,i1,i2,i3,i4)       !15
      elseif (hel.eq.12222) then
       TEMP=AA6(i1,i2,i3,i4,i5)
      elseif (hel.eq.21111) then
       TEMP=AA2(i1,i2,i3,i4,i5)
      elseif (hel.eq.21112) then
       TEMP=AA3(i5,i1,i2,i3,i4)
      elseif (hel.eq.21121) then
       TEMP=AA4(i4,i5,i1,i2,i3)
      elseif (hel.eq.21122) then
       TEMP=AA7(i2,i3,i4,i5,i1)       !20
      elseif (hel.eq.21211) then
       TEMP=AA4(i1,i2,i3,i4,i5)
      elseif (hel.eq.21212) then
       TEMP=AA8(i2,i3,i4,i5,i1)
      elseif (hel.eq.21221) then
       TEMP=AA8(i5,i1,i2,i3,i4)
      elseif (hel.eq.21222) then
       TEMP=AA6(i2,i3,i4,i5,i1)
      elseif (hel.eq.22111) then
       TEMP=AA3(i1,i2,i3,i4,i5)        !25
      elseif (hel.eq.22112) then
       TEMP=AA7(i3,i4,i5,i1,i2)
      elseif (hel.eq.22121) then
       TEMP=AA8(i3,i4,i5,i1,i2)
      elseif (hel.eq.22122) then
       TEMP=AA6(i3,i4,i5,i1,i2)
      elseif (hel.eq.22211) then
       TEMP=AA7(i4,i5,i1,i2,i3)
      elseif (hel.eq.22212) then
       TEMP=AA6(i4,i5,i1,i2,i3)        !30
      elseif (hel.eq.22221) then
       TEMP=AA6(i5,i1,i2,i3,i4)
      elseif (hel.eq.22222) then
       TEMP=AA5(i1,i2,i3,i4,i5)        !32
      endif

      if (iwrite.eq.1) write (14,*) 'AH(',hel,')=',temp

      AH=temp

      End !AH

C----------------------------------------------------------------------------
      function CAH(i1,i2,i3,i4,i5)
C----------------------------------------------------------------------------
      implicit NONE
      double complex CAH,AH
      integer i1,i2,i3,i4,i5

      CAH=dconjg(AH(i1,i2,i3,i4,i5))
      End! CAH
C--------------------------------------------------------------

      function COP(i1,i2,i3,i4,i5,AH)
C     Sum of the cyclically organized permutations of the amplitudes AH

      implicit NONE
      include 'pncommon.f'

      integer i1,i2,i3,i4,i5
      double complex COP,AH,Out
      external AH

      Out=AH(i1,i2,i3,i4,i5)+AH(i1,i3,i2,i4,i5)+AH(i1,i3,i4,i2,i5)+
     &    AH(i1,i3,i4,i5,i2)+AH(i1,i2,i5,i3,i4)+AH(i1,i5,i2,i3,i4)+
     &    AH(i1,i5,i3,i2,i4)+AH(i1,i5,i3,i4,i2)+AH(i1,i2,i4,i5,i3)+
     &    AH(i1,i4,i2,i5,i3)+AH(i1,i4,i5,i2,i3)+AH(i1,i4,i5,i3,i2)

      COP=Out

      End !COP
C--------------------------------------------------------------

      function AHCOP(i1,i2,i3,i4,i5)
C     product of the COP sum of the amplitudes and the conjugated amplitude

      implicit NONE
      include 'pncommon.f'

      integer i1,i2,i3,i4,i5
      double complex COP,AH,AHCOP,tmp1,tmp2,tmp3,tmp4
      external AH

      tmp1=dconjg(AH(i1,i2,i3,i4,i5))

      if (i2.eq.2) then
        tmp2=COP(i1,i2,i3,i4,i5,AH)
      elseif (i3.eq.2) then
        tmp2=COP(i1,i3,i2,i4,i5,AH)
      elseif (i4.eq.2) then
        tmp2=COP(i1,i4,i2,i3,i5,AH)
      elseif (i5.eq.2) then
        tmp2=COP(i1,i5,i2,i3,i4,AH)
      endif
     

      tmp3=tmp1*tmp2

      if (iwrite.eq.1) then
       write(15,15) i1,i2,i3,i4,i5,tmp2
       write(17,17) i1,i2,i3,i4,i5,tmp3
      endif
 
  15  format('COP(',5i1,')=',2g16.8)
  17  format('AHCOP(',5i1,')=',2g16.8)

      AHCOP=tmp3
      End !AHCOP

C----------------------------------------------------------------------------
      function AA1(i1,i2,i3,i4,i5)
C----------------------------------------------------------------------------
C     Primitive amplitude A(i1+,i2+,i3+,i4+,i5+)

      implicit NONE
      include 'pncommon.f'

      integer i1,i2,i3,i4,i5
      double complex AA1,feps,temp


      temp=-ie/6.*cg*(dcmplx(s(i1,i2)*s(i2,i3)+s(i2,i3)*s(i3,i4)
     &                       +s(i3,i4)*s(i4,i5)+s(i4,i5)*s(i5,i1)
     &                       +s(i5,i1)*s(i1,i2))+feps(i1,i2,i3,i4))/
     &           (tr(i1,i2)*tr(i2,i3)*tr(i3,i4)*tr(i4,i5)*tr(i5,i1))


      AA1=temp

      End !AA1

C----------------------------------------------------------------------------
      function AA2(i1,i2,i3,i4,i5)
C----------------------------------------------------------------------------
C     Primitive amplitude A(i1-,i2+,i3+,i4+,i5+)

      implicit NONE
      include 'pncommon.f'

      integer i1,i2,i3,i4,i5
      double complex AA2,num,denom,c23,c34,c45,temp

      c23=(tr(i1,i2)*tr(i1,i3))**2*sq(i2,i3)/tr(i2,i3)
      c34=(tr(i1,i3)*tr(i1,i4))**2*sq(i3,i4)/tr(i3,i4)
      c45=(tr(i1,i4)*tr(i1,i5))**2*sq(i4,i5)/tr(i4,i5)

      num=((s(i2,i3)+s(i3,i4)+s(i4,i5))*(sq(i2,i5)**2)-
     &     sq(i2,i4)*tr(i4,i3)*sq(i3,i5)*sq(i2,i5)-
     &     sq(i1,i2)*sq(i1,i5)/tr(i1,i2)/tr(i1,i5)*
     &     (c23+c34+c45))

      denom=sq(i1,i2)*tr(i2,i3)*tr(i3,i4)*tr(i4,i5)*sq(i5,i1)

      temp=-ie/3.*cg*num/denom

      AA2=temp

      End !AA2

C--------------------------------------------------------------


      function AA3(i1,i2,i3,i4,i5)
C     Primitive amplitude A(i1-,i2-,i3+,i4+,i5+)

      implicit NONE
      include 'pncommon.f'

      integer i1,i2,i3,i4,i5
      double complex AA3,Atree1,Ff,Fs
      double complex fL0,fL2,sum1,minuslog,cr,temp

C Corresponding tree amplitude

      Atree1=ie*tr(i1,i2)**4/
     &       (tr(i1,i2)*tr(i2,i3)*tr(i3,i4)*tr(i4,i5)*tr(i5,i1))


      sum1=tr(i2,i3)*sq(i3,i4)*tr(i4,i1)+tr(i2,i4)*sq(i4,i5)*tr(i5,i1)

      Ff=-0.5*fL0(cr(s(i2,i3),s(i5,i1)))/s(i5,i1)*(tr(i1,i2)**2)*
     &   sum1/(tr(i2,i3)*tr(i3,i4)*tr(i4,i5)*tr(i5,i1))

      Fs=(-Ff-fL2(cr(s(i2,i3),s(i5,i1)))/(s(i5,i1)**3)*sq(i3,i4)*
     &  tr(i4,i1)*tr(i2,i4)*sq(i4,i5)*sum1/tr(i3,i4)/tr(i4,i5)-
     3  tr(i3,i5)*(sq(i3,i5)**3)/sq(i1,i2)/sq(i2,i3)/tr(i3,i4)/
     &  tr(i4,i5)/sq(i5,i1)+
     4  tr(i1,i2)*(sq(i3,i5)**2)/sq(i2,i3)/tr(i3,i4)/tr(i4,i5)/sq(i5,i1)
     5  +0.5*tr(i1,i2)*sq(i3,i4)*tr(i4,i1)*tr(i2,i4)*sq(i4,i5)/s(i2,i3)
     &  /tr(i3,i4)/tr(i4,i5)/s(i5,i1))/3d0



      temp=-cg*(ie*(Ff+Fs)+Atree1/3d0*(minuslog(s(i2,i3))+
     &    minuslog(s(i5,i1))))


      AA3=temp
      End !AA3

C--------------------------------------------------------------

      function AA4(i1,i2,i3,i4,i5)
C     Primitive amplitude A(i1-,i2+,i3-,i4+,i5+)

      implicit NONE
      include 'pncommon.f'

      integer i1,i2,i3,i4,i5
      double complex AA4,Atree2,Ff,Fs,r23,r12,r00,d15,d12,d45,d54,temp
      double complex q23,q12,z23,z12,z34,fLs1,fL0,fL1,fL2,fLi2,minuslog
      double complex cr, s2351,s3451,s1234,s5134
      double precision s251,s234,s351,s334

C Various ratios and powers of scalar products
      s2351=cr(s(i2,i3),s(i5,i1))
      s3451=cr(s(i3,i4),s(i5,i1))
      s1234=cr(s(i1,i2),s(i3,i4))
      s5134=cr(s(i5,i1),s(i3,i4))
      s251=s(i5,i1)*s(i5,i1)
      s234=s(i3,i4)*s(i3,i4)
      s351=s251*s(i5,i1)
      s334=s234*s(i3,i4)

C Various spinor products
      d15=tr(i1,i5)*sq(i5,i2)*tr(i2,i3)-tr(i3,i4)*sq(i4,i2)*tr(i2,i1)
      d12=tr(i1,i2)*tr(i2,i3)*tr(i3,i4)*tr(i4,i5)*tr(i5,i1)
      d45=tr(i4,i5)*tr(i5,i1)*tr(i2,i4)
      d54=tr(i5,i4)*tr(i4,i3)*tr(i2,i5)

C Corresponding tree amplitude

      Atree2=ie*tr(i1,i3)**4/d12

C Logarithmic functions for Ff

      r23=fLs1(s2351,s3451)/s251
      r12=fLs1(s1234,s5134)/s234
      r00=0.5*fL0(s3451)/s(i5,i1)*(tr(i1,i3)**3)*d15/d12

      Ff=-r23*(tr(i1,i3)*sq(i2,i4))**2*tr(i4,i1)/tr(i4,i5)/tr(i5,i1)
     &   +r12*(tr(i1,i3)*sq(i2,i5))**2*tr(i5,i3)/tr(i3,i4)/tr(i4,i5)
     &   -r00

C Logarithmic functions for Fs
      q23=2*r23+(fL1(s2351)+fL1(s3451))/s251
      q12=2*r12+(fL1(s1234)+fL1(s5134))/s234
      z23=2./3.*fL2(s2351)/s351
      z12=2./3.*fL2(s1234)/s334
      z34=1./3.*fL2(s3451)/s351

      Fs=-q23*tr(i1,i2)*tr(i2,i3)*tr(i3,i4)*
     &  (tr(i4,i1)*sq(i2,i4))**2/d45/tr(i2,i4)
     2   +q12*tr(i3,i2)*tr(i2,i1)*tr(i1,i5)*
     &  (tr(i5,i3)*sq(i2,i5))**2/d54/tr(i2,i5)
     3   +z23*(tr(i2,i3)**2)*((tr(i4,i1)*sq(i2,i4))**3)/d45
     4   -z12*(tr(i2,i1)**2)*((tr(i5,i3)*sq(i2,i5))**3)/d54
     5   +z34*(tr(i1,i3)*sq(i2,i4)*sq(i2,i5)*d15/tr(i4,i5)
     6   +2.*(tr(i1,i2)*tr(i3,i4))**2*tr(i4,i1)*(sq(i2,i4)**3)/d45
     7   -2.*(tr(i3,i2)*tr(i1,i5))**2*tr(i5,i3)*(sq(i2,i5)**3)/d54)

      Fs=Fs+1./3.*(r00+(sq(i2,i4)*sq(i2,i5))**2/sq(i1,i2)/sq(i2,i3)
     &  /sq(i3,i4)/tr(i4,i5)/sq(i5,i1)
     3  -tr(i1,i2)*(tr(i4,i1)**2)*(sq(i2,i4)**3)/d45
     &   /sq(i2,i3)/sq(i3,i4)/s(i5,i1)
     4  +tr(i3,i2)*(tr(i5,i3)**2)*(sq(i2,i5)**3)/d54
     &   /sq(i2,i1)/sq(i1,i5)/s(i3,i4)
     5  +0.5*(tr(i1,i3)**2)*sq(i2,i4)*sq(i2,i5)/s(i3,i4)
     6   /tr(i4,i5)/s(i5,i1))


      temp=-cg*(ie*(Ff+Fs)+Atree2/3d0*(minuslog(s(i3,i4))+
     &    minuslog(s(i5,i1))))

 
      AA4=temp

      End !AA4


C-------------------------------------------------------------
C The following amplitudes are obtained from the amplitudes AA1-AA4
C by the replacement <> --> [] and the inversion of the overall sign
C
      function AA5(i1,i2,i3,i4,i5)
C     Primitive amplitude A(i1-,i2-,i3-,i4-,i5-)

      implicit NONE
      include 'pncommon.f'

      integer i1,i2,i3,i4,i5
      double complex AA5,feps,temp


      temp=ie/6.*cg*(dcmplx(s(i1,i2)*s(i2,i3)+s(i2,i3)*s(i3,i4)
     &                       +s(i3,i4)*s(i4,i5)+s(i4,i5)*s(i5,i1)
     &                       +s(i5,i1)*s(i1,i2))-feps(i1,i2,i3,i4))/
     &           (sq(i1,i2)*sq(i2,i3)*sq(i3,i4)*sq(i4,i5)*sq(i5,i1))


      AA5=temp

      End !AA5

C--------------------------------------------------------------

      function AA6(i1,i2,i3,i4,i5)
C     Primitive amplitude A(i1+,i2-,i3-,i4-,i5-)

      implicit NONE
      include 'pncommon.f'

      integer i1,i2,i3,i4,i5
      double complex AA6,num,denom,c23,c34,c45,temp

      c23=(sq(i1,i2)*sq(i1,i3))**2*tr(i2,i3)/sq(i2,i3)
      c34=(sq(i1,i3)*sq(i1,i4))**2*tr(i3,i4)/sq(i3,i4)
      c45=(sq(i1,i4)*sq(i1,i5))**2*tr(i4,i5)/sq(i4,i5)

      num=(dcmplx(s(i2,i3)+s(i3,i4)+s(i4,i5))*(tr(i2,i5)**2)-
     &     tr(i2,i4)*sq(i4,i3)*tr(i3,i5)*tr(i2,i5)-
     &     tr(i1,i2)*tr(i1,i5)/sq(i1,i2)/sq(i1,i5)*
     &     (c23+c34+c45))

      denom=tr(i1,i2)*sq(i2,i3)*sq(i3,i4)*sq(i4,i5)*tr(i5,i1)

      temp=ie/3.*cg*num/denom

      AA6=temp

      End !AA6

C--------------------------------------------------------------

      function AA7(i1,i2,i3,i4,i5)
C     Primitive amplitude A(i1+,i2+,i3-,i4-,i5-)

      implicit NONE
      include 'pncommon.f'

      integer i1,i2,i3,i4,i5
      double complex AA7,Atree1,Ff,Fs
      double complex fL0,fL2,sum1,minuslog,cr,temp

C Corresponding tree amplitude

      Atree1=ie*sq(i1,i2)**4/
     &       (sq(i1,i2)*sq(i2,i3)*sq(i3,i4)*sq(i4,i5)*sq(i5,i1))


      sum1=sq(i2,i3)*tr(i3,i4)*sq(i4,i1)+sq(i2,i4)*tr(i4,i5)*sq(i5,i1)

      Ff=-0.5*fL0(cr(s(i2,i3),s(i5,i1)))/s(i5,i1)*(sq(i1,i2)**2)*
     &   sum1/(sq(i2,i3)*sq(i3,i4)*sq(i4,i5)*sq(i5,i1))

      Fs=(-Ff-fL2(cr(s(i2,i3),s(i5,i1)))/(s(i5,i1)**3)*tr(i3,i4)*
     &  sq(i4,i1)*sq(i2,i4)*tr(i4,i5)*sum1/sq(i3,i4)/sq(i4,i5)-
     3  sq(i3,i5)*(tr(i3,i5)**3)/tr(i1,i2)/tr(i2,i3)/sq(i3,i4)/
     &  sq(i4,i5)/tr(i5,i1)+
     4  sq(i1,i2)*(tr(i3,i5)**2)/tr(i2,i3)/sq(i3,i4)/sq(i4,i5)/tr(i5,i1)
     5  +0.5*sq(i1,i2)*tr(i3,i4)*sq(i4,i1)*sq(i2,i4)*tr(i4,i5)/s(i2,i3)
     &  /sq(i3,i4)/sq(i4,i5)/s(i5,i1))/3.


      temp=cg*(ie*(Ff+Fs)+Atree1/3d0*(minuslog(s(i2,i3))+
     &    minuslog(s(i5,i1))))


      AA7=temp

      End !AA7

C--------------------------------------------------------------

      function AA8(i1,i2,i3,i4,i5)
C     Primitive amplitude A(i1+,i2-,i3+,i4-,i5-)

      implicit NONE
      include 'pncommon.f'

      integer i1,i2,i3,i4,i5
      double complex AA8,Atree2,Ff,Fs,r23,r12,r00,d15,d12,d45,d54,temp
      double complex q23,q12,z23,z12,z34,fLs1,fL0,fL1,fL2,fLi2,minuslog
      double complex cr, s2351,s3451,s1234,s5134
      double precision s251,s234,s351,s334

C Various ratios and powers of scalar products
      s2351=cr(s(i2,i3),s(i5,i1))
      s3451=cr(s(i3,i4),s(i5,i1))
      s1234=cr(s(i1,i2),s(i3,i4))
      s5134=cr(s(i5,i1),s(i3,i4))

      s251=s(i5,i1)*s(i5,i1)
      s234=s(i3,i4)*s(i3,i4)
      s351=s251*s(i5,i1)
      s334=s234*s(i3,i4)

C Various spinor products
      d15=sq(i1,i5)*tr(i5,i2)*sq(i2,i3)-sq(i3,i4)*tr(i4,i2)*sq(i2,i1)
      d12=sq(i1,i2)*sq(i2,i3)*sq(i3,i4)*sq(i4,i5)*sq(i5,i1)
      d45=sq(i4,i5)*sq(i5,i1)*sq(i2,i4)
      d54=sq(i5,i4)*sq(i4,i3)*sq(i2,i5)

C Corresponding tree amplitude

      Atree2=ie*sq(i1,i3)**4/d12

C Logarithmic functions for Ff

      r23=fLs1(s2351,s3451)/s251
      r12=fLs1(s1234,s5134)/s234
      r00=0.5*fL0(s3451)/s(i5,i1)*(sq(i1,i3)**3)*d15/d12

      Ff=-r23*(sq(i1,i3)*tr(i2,i4))**2*sq(i4,i1)/sq(i4,i5)/sq(i5,i1)
     &   +r12*(sq(i1,i3)*tr(i2,i5))**2*sq(i5,i3)/sq(i3,i4)/sq(i4,i5)
     &   -r00

C Logarithmic functions for Fs
      q23=2*r23+(fL1(s2351)+fL1(s3451))/s251
      q12=2*r12+(fL1(s1234)+fL1(s5134))/s234
      z23=2./3.*fL2(s2351)/s351
      z12=2./3.*fL2(s1234)/s334
      z34=1./3.*fL2(s3451)/s351

      Fs=-q23*sq(i1,i2)*sq(i2,i3)*sq(i3,i4)*
     &  (sq(i4,i1)*tr(i2,i4))**2/d45/sq(i2,i4)
     2   +q12*sq(i3,i2)*sq(i2,i1)*sq(i1,i5)*
     &  (sq(i5,i3)*tr(i2,i5))**2/d54/sq(i2,i5)
     3   +z23*(sq(i2,i3)**2)*((sq(i4,i1)*tr(i2,i4))**3)/d45
     4   -z12*(sq(i2,i1)**2)*((sq(i5,i3)*tr(i2,i5))**3)/d54
     5   +z34*(sq(i1,i3)*tr(i2,i4)*tr(i2,i5)*d15/sq(i4,i5)
     6   +2.*(sq(i1,i2)*sq(i3,i4))**2*sq(i4,i1)*(tr(i2,i4)**3)/d45
     7   -2.*(sq(i3,i2)*sq(i1,i5))**2*sq(i5,i3)*(tr(i2,i5)**3)/d54)

      Fs=Fs+1./3.*(r00+(tr(i2,i4)*tr(i2,i5))**2/tr(i1,i2)/tr(i2,i3)
     &  /tr(i3,i4)/sq(i4,i5)/tr(i5,i1)
     3  -sq(i1,i2)*(sq(i4,i1)**2)*(tr(i2,i4)**3)/d45
     &   /tr(i2,i3)/tr(i3,i4)/s(i5,i1)
     4  +sq(i3,i2)*(sq(i5,i3)**2)*(tr(i2,i5)**3)/d54
     &   /tr(i2,i1)/tr(i1,i5)/s(i3,i4)
     5  +0.5*(sq(i1,i3)**2)*tr(i2,i4)*tr(i2,i5)/s(i3,i4)
     6   /sq(i4,i5)/s(i5,i1))


      
      temp=cg*(ie*(Ff+Fs)+Atree2/3d0*(minuslog(s(i3,i4))+
     &    minuslog(s(i5,i1))))


      AA8=temp

      End !AA8


C-------------------------------------------------------------
      function feps(i1,i2,i3,i4)
C     Vector product of 4 4-momenta

      implicit NONE
      include 'pncommon.f'

      integer i1,i2,i3,i4,i5
      double complex feps

      feps=Sq(i1,i2)*Tr(i2,i3)*Sq(i3,i4)*Tr(i4,i1)-
     &    Tr(i1,i2)*Sq(i2,i3)*Tr(i3,i4)*Sq(i4,i1)

      End !feps

C#####################################################################
C---------------------------------------------------------------------

      function fL0(x)

      implicit NONE
      include 'pncommon.f'
      double precision zerom
      double complex fL0,x,temp


      zerom=1d-9      ! "Machine zero", is platform dependent

      if (cdabs(x).lt.zerom) then
       print *,'fL0: wrong argument, x=', x
       STOP
      endif

      if (cdabs(1d0-x).lt.zerom) then
       temp=-1d0
      else
       temp=cdlog(x)/dcmplx(1-x)
      endif   

      fL0=temp

      End !fL0

C---------------------------------------------------------------------
      function fL1(x)

      implicit NONE
      include 'pncommon.f'
      double precision zerom
      double complex fL1,x,temp

      zerom=1d-9     ! "Machine zero", is platform dependent

      if (cdabs(x).lt.zerom) then
       print *,'fL1: wrong argument, x=', x
       STOP
      endif

      if (cdabs(1d0-x).lt.zerom) then
        temp=-0.5
      else
        temp=cdlog(x)/(1-x)**2+1/(1-x)
      endif

      fL1=temp
 
      End !fL1
C---------------------------------------------------------------------
      function fL2(x)

      implicit NONE
      include 'pncommon.f'
      double precision zerom
      double complex fL2,x,temp,z

      zerom=1d-4    ! Cutoff,  platform dependent

      if (cdabs(x).lt.zero) then
       print *,'fL2: wrong argument, x=', x
       STOP
      endif

C This treats numeric complications around x=1

      if (cdabs(1d0-x).lt.zerom ) then
        z=dcmplx(x-1d0)
        temp=dcmplx(-1d0/3d0+0.25*z-0.2*z*z+0.5/x)
      else
        temp=(cdlog(x)-(x-1/x)/2)/(1-x)**3
      endif

      fL2=temp

      End !fL2


C----------------------------------------------------------------------
      function fLs1(x1,x2)

      implicit NONE
      include 'pncommon.f'
      double precision zero1,zero2 
      double complex fLs1,fLi2,x1,x2,dif,fL0,l1,l2,temp

      zero1=1d-5    
      zero2=1d-5
      Pi= 3.141592653589793238d0

      dif=1-x1-x2

      if (cdabs(x1).lt.zero.or.cdabs(x2).lt.zero) 
     &  then
       print *,'fLs1: wrong arguments, x1=', x1,' x2=', x2
       STOP
      endif

      if (cdabs(dif).lt.zero1) then
        if (cdabs(x1).lt.zero2) then
cpn2006          temp=-0.5*(cdlog(x1)+0.5+2d0*x1*cdlog(x1)) ! wrong 
          temp=-0.25d0-x1/3d0-cdlog(x1)/2d0-x1*cdlog(x1) 
        elseif (cdabs(x2).lt.zero2) then
cpn2006          temp=-0.5*(cdlog(x2)+0.5+2d0*x2*cdlog(x2)) ! wrong
          temp=-0.25d0-x2/3d0-cdlog(x2)/2d0-x2*cdlog(x2) 
        else !x1,x2>zero2
          temp=-0.5*(cdlog(x1)/x2**2+cdlog(x2)/x1**2+1d0/x1/x2)
        endif

      else ! dif>zero1
        temp=(fLi2(1-x1)+fLi2(1-x2)+cdlog(x1)*cdlog(x2)-pi*pi/6d0+
     &    dif*(fL0(x1)+fL0(x2)))/dif/dif
      endif
     
      fLs1=temp

      End !fLs1


C----------------------------------------------------------------------
      function fLs2(x1,x2)

      implicit NONE
      include 'pncommon.f'
      double precision zero1,zero2 
      double complex fLs2,fLs1,x1,x2,dif,fL1,temp

      zero1=1d-5    
      zero2=1d-5
      Pi= 3.141592653589793238d0

      dif=1-x1-x2

      if (cdabs(x1).lt.zero.or.cdabs(x2).lt.zero) 
     &  then
       print *,'fLs2: wrong arguments, x1=', x1,' x2=', x2
       STOP
      endif

      if (cdabs(dif).lt.zero1) then
        if (cdabs(x1).lt.zero2) then          
          temp=-7d0/18d0 - 3d0*x1/4. - cdlog(x1)/3.d0 - x1*cdlog(x1)
        elseif (cdabs(x2).lt.zero2) then
          temp=-7d0/18d0 - 3d0*x2/4. - cdlog(x2)/3.d0 - x2*cdlog(x2)
        else !x1,x2>zero2
          temp=-1/(4.*x1**2) - 1/(4.*x2**2) + 1/(24.*x1**2*x2**2)
     >      - 1/(8.*x1*x2**2) - 1/(8.*x1**2*x2)
     >      - cdlog(x1)/(3.*x2**3) - cdlog(x2)/(3.*x1**3)
        endif

      else ! dif>zero1
        temp = ( fLs1(x1,x2) + 0.5*(fL1(x1)+fL1(x2)) )/dif
      endif
     
      fLs2=temp

      End !fLs2

C----------------------------------------------------------------------
      function fLs3(x1,x2)

      implicit NONE
      include 'pncommon.f'
      double precision zero1,zero2 
      double complex fLs3,fLs2,x1,x2,dif,fL2,temp

      zero1=1d-5    
      zero2=1d-5
      Pi= 3.141592653589793238d0

      dif=1-x1-x2

      if (cdabs(x1).lt.zero.or.cdabs(x2).lt.zero) 
     &  then
       print *,'fLs3: wrong arguments, x1=', x1,' x2=', x2
       STOP
      endif

      if (cdabs(dif).lt.zero1) then
        if (cdabs(x1).lt.zero2) then          
          temp= -19d0/48d0 - 31d0*x1/30. - cdlog(x1)/4. - x1*cdlog(x1)

        elseif (cdabs(x2).lt.zero2) then
          temp= -19d0/48d0 - 31d0*x2/30. - cdlog(x2)/4. - x2*cdlog(x2)
        else !x1,x2>zero2
          temp= 1/(4.*x1**4) - 1/(24.*x1**3) - 1/(12.*x1**2) + 
     >      1/(384.*x1**3*x2**4) -  1/(192.*x1**2*x2**4) - 
     >      95/(384.*x1*x2**4) + x1/(4.*x2**4) + 5/(24.*x2**3) + 
     >      1/(192.*x1**3*x2**3) + 9/(128.*x1**2*x2**3) -
     >      13/(384.*x1*x2**3) - 1/(12.*x2**2) + 13/(192.*x1**3*x2**2) +
     >      1/(64.*x1**2*x2**2) - 1/(4.*x1**4*x2) - 13/(384.*x1**3*x2) - 
     >      Cdlog(x1)/(4.*x2**4) - Cdlog(x2)/(4.*x1**4)
        endif

      else ! dif>zero1
        temp = ( fLs2(x1,x2) + (fL2(x1)+fL2(x2))/3d0 
     >    - ( 1d0/x1 + 1d0/x2 )/6d0 )/dif
      endif
     
      fLs3=temp

      End !fLs3


C----------------------------------------------------------------------------
CsB   fLi2(x) and fLi2_2(x) were moved to resbosm.f 
!       function fLi2(x)
! C     Dilogarithm function for x+-i0. The precision is about 12 figures. 
! C     Translation of Carl Schmidt's C code. The imaginary part for real x>1
! C     is generated according to x -> x+i0
! 
!       implicit NONE
!       include 'pncommon.f'
! 
!       double complex fLi2,Out,x,fLi2_2
!       double precision czero,xr
!       external fLi2_2 
!       Data epsil / 1.e-10 / 
!   
!       czero=1d1*epsil
!       pi=3.14159265358979d0
! 
! 
! C     Check if the imaginary part of x is (almost) infinitesimally small
!       if (dabs(dimag(x)-epsil).gt.czero) then
!         print *,'fLi2: too large imaginary part of x=',dimag(x)
!         Print *,'epsil,czero:',epsil,czero
!         STOP
!       endif
! 
!       xr=dreal(x)
! 
! 
!       if (xr.ge.0.and.xr.le.1.0) then
!         Out=fLi2_2(x)
!       elseif (xr.lt.0) then
!         Out=fLi2_2(1/(1-x))+0.5*(dlog(1-xr))**2-
!      &      dlog(-xr)*dlog(1-xr)-pi*pi/6.
!       elseif (xr.gt.1.0) then
!         Out=-fLi2_2(1/x)+pi*pi/3.-0.5*(dlog(xr))**2
! 
!         if (dimag(x).ge.0) then
!           Out=Out+ie*pi*dlog(xr)
!         else
!           Out=Out-ie*pi*dlog(xr)
!         endif
! 
!       endif
! 
!       fLi2=Out
! 
!       End ! fLi2
! 
! C---------------------------------------------------------
!       double complex function fLi2_2(x)
!       implicit NONE
!       include 'pncommon.f'
!      
!       double complex x,Out
!       double precision xr,r,xm,L
! 
!       integer NCUT, i,N
!       NCUT =27
! 
!       L=0.0
!       i=0
!       r=1.0
! 
!       xr=dreal(x)
! 
!       if (xr.ge.0.and.xr.le.0.5) then
! 
!        do while(i.le.NCUT)
!        	i=i+1
!         r=r*xr
! 
!         L=L+r/i/i
!        enddo
!         Out=L
!       elseif(xr.gt.0.5.and.xr.le.1.0) then
!        xm = 1.0-xr
!        L = pi*pi/6.0 - dlog(xr)*dlog(xm)
! 
!        do while(i.le.NCUT)
!         i=i+1
!         r=r*xm
! 
!         L=L-r/i/i
!        enddo
!       endif 
!       
!       Out=L  
! 
! 
!       fLi2_2=Out      
!       return
!       end! fLi2_2 
!        
! 
C----------------------------------------------------------



      function minuslog(x)
C Calculate the complex log(-x), with the account for the imaginary part
      implicit NONE
      include 'pncommon.f'
      double precision x
      double precision x1
      double complex minuslog,temp
     

      if (x.lt.zero) then
       temp=dcmplx(dlog(-x))
      elseif (x.gt.zero) then
       temp=dcmplx(dlog(x))-ie*pi
      else 
       print *,'MINUSLOG: abs(x)=',dabs(x),' < ',zero
       STOP
      endif
     
      minuslog=temp

      End! minuslog
C----------------------------------------------------------


      function cr(x1,x2)
C A ratio of two numbers with infinitely small imaginary parts
      implicit NONE
      include 'pncommon.f'
      double precision x1,x2
      double complex cr,temp

      temp=dcmplx(x1/x2)
      
C Take care of the imaginary part if x1 and x2 have different signs

      if (dreal(temp).lt.0) then
 
        if (x1.gt.0.and.x2.lt.0) then
         temp=temp-ie*epsil
        elseif (x1.lt.0.and.x2.gt.0) then
         temp=temp+ie*epsil
        else  !this can be deleted in the future 
         print *,'cr:x1 and x2 have the same signs but cr<0. STOP'
         STOP
        endif

      endif 

      cr=temp 

      End !cr


C========================================================================
C Amplitude squared for Q+G->diphoton+Q via the box diagram
C========================================================================
C
C----------------------------------------------------------------------------
      double precision function PentagonQGAAQ (i1,i2,i3,i4,i5)
C-   Returns the squared matrix element for
C-   the process q(i4) +g(i5) -> gamma(i1) + gamma(i2) + q(i3) via the
cpn2006 quark box Feynman graph
C    The amplitude expression is taken from Y. Yasui, hep-ph/0203163 (2002)
C----------------------------------------------------------------------------
      IMPLICIT NONE
      include 'pncommon.f'
      integer i,j,iret,isum,sgn, i1, i2, i3, i4, i5
      double precision AmpSqdQG,xi(4)
CsB___Info from ResBos to Pentagon
      DOUBLE PRECISION W,QCDL,X1X2MIN
      INTEGER*4 NJET,N_PRIMARY_PART,NF,I_STR,IBEAM,NPART,NEVENT,NORDER
      COMMON/AA/W,QCDL,X1X2MIN,NJET,N_PRIMARY_PART,NF,I_STR,IBEAM,NPART,
     &  NEVENT,NORDER
      DOUBLE PRECISION PLAB,PCEM,WHAT
      COMMON/MOM/PLAB(4,12),PCEM(4,12),WHAT
CsB___Locally defined
      double precision tmp
CsB___Debug options
      Logical Debug
      Pi= 3.141592653589793238d0

      iDebug = 0      ! no print from Pentagon
      iWrite = 0      ! no write from Pentagon
      epsil  = 1.e-10 ! numeric zero
cpn2005 To produce a specific helicity amplitude
      iSum   = 1      ! sum over helicities

c$$$      do i=1,5
c$$$        h(i) = 1
c$$$      enddo


C_____The definitions of the four momenta which Pentagon uses in terms of the
C     four momenta which are used in ResBos 
C     Note:             in ResBos        in Pentagon/BDK paper
C           quark 1 in     1                 2
C           gluon 2 in     2                 5
C           photon 1 out   3                 3
C           photon 2 out   4                 4
C           quark out      5                 1
C      a = particle index
C      b = coordinate index
C        a b         b a
      pp(2,1) = pLab(4,i1) + pLab(3,i1)       ! p+  of quark 1 in
      pp(2,2) = pLab(4,i1) - pLab(3,i1)       ! p-  of quark 1 in
      pp(2,3) = 0.d0
      pp(5,1) = pLab(4,i2) + pLab(3,i2)       ! p+  of gluon 3 in
      pp(5,2) = pLab(4,i2) - pLab(3,i2)       ! p-  of gluon 3 in
      pp(5,3) = 0.d0
      pp(3,1) = pLab(4,i3) + pLab(3,i3)       ! p+  of photon 3
      pp(3,2) = pLab(4,i3) - pLab(3,i3)       ! p-  of photon 3
      pp(3,3) = ASin(pLab(2,i3)/Sqrt(pLab(1,i3)**2+pLab(2,i3)**2))
      If (pLab(1,i3).lt.0.0) pp(3,3) = Pi - pp(3,3)

cpn Alternative definition of the angle
c      pp(3,3) = ACos(pLab(4,i3)/Sqrt(pLab(1,i3)**2+pLab(2,i3)**2))
c      pp(3,3)=pp(3,3)*sign(1,pLab(2,i3))

      pp(4,1) = pLab(4,i4) + pLab(3,i4)       ! p+  of photon 4
      pp(4,2) = pLab(4,i4) - pLab(3,i4)       ! p-  of photon 4
      pp(4,3) = ASin(pLab(2,i4)/Sqrt(pLab(1,i4)**2+pLab(2,i4)**2))
      If (pLab(1,i4).lt.0.0) pp(4,3) = Pi - pp(4,3)
cpn 
c      pp(4,3) = ACos(pLab(1,i4)/Sqrt(pLab(1,i4)**2+pLab(2,i4)**2))
c      pp(4,3)=pp(4,3)*sign(1,pLab(2,i4))

      pp(1,1) = pLab(4,i5) + pLab(3,i5)       ! p+  of quark out
      pp(1,2) = pLab(4,i5) - pLab(3,i5)       ! p-  of quark out
      pp(1,3) = ASin(pLab(2,i5)/Sqrt(pLab(1,i5)**2+pLab(2,i5)**2))
      If (pLab(1,i5).lt.0.0) pp(1,3) = Pi - pp(1,3)
cpn
c      pp(1,3) = ACos(pLab(1,i5)/Sqrt(pLab(1,i5)**2+pLab(2,i5)**2))
c        pp(1,3)=pp(1,3)*sign(1,pLab(2,i5))



CsB___Output debug info headers
      if (iwrite.eq.1) then
        open(unit=14,file='ah.dat',status='new')
        open(unit=15,file='cop.dat',status='new')
        open(unit=17,file='ahcop.dat',status='new')
      endif

CsB___Execute Pentagon subroutines
CPN xi2mom is only needed if pentagon is called standalone
C      Call xi2mom(W,xi,sgn)

      Call Standard          ! Get rid of this

      Call ConstInit1(2)     ! This part can go out of the loop 

      Call KinemInit      !Initialize kinematics 

      Call ConstInit2     !Set up the couplings

      Call TensorInit     !Fill in spinor products 

      PentagonQGAAQ=AmpSqdQG(isum,iret)

C_____Output to the screen
      if (idebug.eq.1) then
        print *,'isum=',isum
        if (isum.eq.1) print *,'hel=',
     &    h(1)*10000+h(2)*1000+h(3)*100+h(4)*10+h(5)
        print *,'AmpSqdQG=',PentagonQGAAQ
      endif
!     Print *,'AmpSqdQG=',PentagonQGAAQ

CsB___Clean up
cpn      close(13)
      if (iwrite.eq.1) then
       close(14)
       close(15)
       close(17)
      endif

      RETURN
      END ! PentagonQGAAQ


cpn=====================================================================

C----------------------------------------------------------------------------
      function AmpSqdQG(isum,iret)
C----------------------------------------------------------------------------

      implicit NONE
      include 'pncommon.f'

      double precision AmpSqdQG
      integer k1,k2,k3,k4,k5
      integer iret,isum
      double complex Out,tmp1, AmpQGAAQ, S3_345
      external AmpQGAAQ

cpn2006 Sum over all permutations of 3, 4, 5
      S3_345()= AmpQGAAQ(1,2,3,4,5) + AmpQGAAQ(1,2,3,5,4)
     >        + AmpQGAAQ(1,2,4,3,5) + AmpQGAAQ(1,2,4,5,3)
     >        + AmpQGAAQ(1,2,5,3,4) + AmpQGAAQ(1,2,5,4,3)


      iret=0

      Out=(0.,0.)

      if (isum.eq.0) then
C no summation over helicities
   
        tmp1=S3_345()
        Out=tmp1*dconjg(tmp1)

    
        if (idebug.eq.1) print 15,h(1),h(2),h(3),h(4),h(5),Out
   15   format (1x,'hel=',5i1,'  Out=',2g16.8)
        
      elseif (isum.eq.1) then

        do  10 k1=1,2
          h(1)=k1
          h(2)=3-k1
          do  10 k3=1,2
            h(3)=k3
            do  10 k4=1,2
              h(4)=k4
              do  10 k5=1,2
                h(5)=k5

                tmp1=S3_345()
        
                Out=Out+tmp1*dconjg(tmp1) 
             
 10           continue
            else
              print *,'isum=',isum,' not implemented yet!'
              STOP
            endif               !isum

            AmpSqdQG=Dreal(Out)*const

      End !ampsqdGG


C----------------------------------------------------------------------------
      function AmpQGAAQ(i1,i2,i3,i4,i5)
C----------------------------------------------------------------------------
C      Creates an amplitude for any set of helicities

      implicit NONE
      include 'pncommon.f'

      integer i1,i2,i3,i4,i5
      integer*4 hel
      double complex AmpQGAAQ,AQ1, AQ2, AQ3, AQ4, AQ5, AQ6, AQ7, AQ8,
     >  temp

      hel=h(i5)+h(i4)*10+h(i3)*100+h(i2)*1000+h(i1)*10000

      if (hel.eq.12111) then
       TEMP =-AQ1(i2,i1,i5,i4,i3)       
      elseif (hel.eq.12112) then
       TEMP =-AQ2(i2,i1,i5,i4,i3)       
      elseif (hel.eq.12121) then
       TEMP =-AQ3(i2,i1,i5,i4,i3)
      elseif (hel.eq.12122) then
       TEMP =-AQ6(i1,i2,i3,i4,i5)

      elseif (hel.eq.12211) then
       TEMP =-AQ4(i2,i1,i5,i4,i3)
      elseif (hel.eq.12212) then
       TEMP =-AQ7(i1,i2,i3,i4,i5)
      elseif (hel.eq.12221) then
       TEMP =-AQ8(i1,i2,i3,i4,i5)       
      elseif (hel.eq.12222) then
       TEMP =-AQ5(i1,i2,i3,i4,i5)

      elseif (hel.eq.21111) then
       TEMP = AQ1(i1,i2,i3,i4,i5)
      elseif (hel.eq.21112) then
       TEMP = AQ4(i1,i2,i3,i4,i5)
      elseif (hel.eq.21121) then
       TEMP = AQ3(i1,i2,i3,i4,i5)
      elseif (hel.eq.21122) then
       TEMP = AQ8(i2,i1,i5,i4,i3)       

      elseif (hel.eq.21211) then
       TEMP = AQ2(i1,i2,i3,i4,i5)
      elseif (hel.eq.21212) then
       TEMP = AQ7(i2,i1,i5,i4,i3)
      elseif (hel.eq.21221) then
       TEMP = AQ6(i2,i1,i5,i4,i3)
      elseif (hel.eq.21222) then
       TEMP = AQ5(i2,i1,i5,i4,i3)
      endif

      if (iwrite.eq.1) write (14,*) 'AmpQGAAQ(',hel,') =  ',temp

      AmpQGAAQ =  temp

      End !AmpQGAAQ


C----------------------------------------------------------------------------
      function AQ1(i1,i2,i3,i4,i5)
C----------------------------------------------------------------------------
cpn  A primitive amplitude A_5(qbar(i1-),q(i2+),V(i3+_),V(i4+),V(i5+))
c    containing the quark box loop, where V=gluon or photon

      implicit NONE
      include 'pncommon.f'

      integer i1,i2,i3,i4,i5
      double complex AQ1, Fs

      
      Fs = - ie/3d0*(
     >          Tr(i1,i3)*Sq(i3,i4)*Tr(i4,i1)**2/
     >          Tr(i1,i2)/Tr(i3,i4)**2/Tr(i4,i5)/Tr(i5,i1)
     >  +       Tr(i1,i4)*Tr(i2,i4)*Sq(i4,i5)*Tr(i5,i1)/
     >          Tr(i1,i2)/Tr(i2,i3)/Tr(i3,i4)/Tr(i4,i5)**2
     >  +       Sq(i2,i3)*Sq(i2,i5)/Sq(i1,i2)/Tr(i3,i4)/Tr(i4,i5) )


      AQ1=-Fs/16d0/pi**2

      end !AQ1


C----------------------------------------------------------------------------
      function AQ2(i1,i2,i3,i4,i5)
C----------------------------------------------------------------------------
cpn A finite (free of infrared poles) part of the primitive amplitude 
c   A_5(qbar(i1-),q(i2+),V(i3-),V(i4+),V(i5+)) containing the quark 
c   box loop, where V=gluon or photon

      implicit NONE
      include 'pncommon.f'

      integer i1,i2,i3,i4,i5
      double complex AQ2, Fs, Ff, cr, fL0, fL1, fL2

      Fs = 1d0/3d0*(
     >     2d0*Tr(i1,i5)*Sq(i2,i5)*Tr(i3,i4)*Tr(i3,i5)*Sq(i4,i5)**2
     >          /Tr(i4,i5)*fL2(cr(S(i1,i2),S(i3,i4)))/ S(i3,i4)**3 
     >  -  Tr(i1,i3)*Tr(i1,i5)*Tr(i3,i4)*Sq(i4,i5)**2/Tr(i1,i2)
     >          /Tr(i4,i5)*fL1(cr(S(i1,i2),S(i3,i4)))/ S(i3,i4)**2 
     >  -  Tr(i1,i3)*Sq(i2,i4)*Sq(i4,i5)
     >          /Tr(i1,i2)/Sq(i1,i2)/Sq(i3,i4)/Tr(i4,i5) 
     >  +  Sq(i2,i4)**2*Sq(i2,i5)/Sq(i1,i2)/Sq(i2,i3)/Sq(i3,i4)
     >         /Tr(i4,i5) )

      Ff = - Tr(i1,i3)**2*Sq(i4,i5)/ Tr(i1,i2)/Tr(i4,i5)
     >  *fL0(cr(S(i1,i2),S(i3,i4))) / S(i3,i4) 

      AQ2=-(Fs+Ff)/(16d0*pi**2)
      
      end !AQ2



C----------------------------------------------------------------------------
      function AQ3(i1,i2,i3,i4,i5)
C----------------------------------------------------------------------------
cpn  A finite (free of infrared poles) part of the primitive amplitude
c    A_5(qbar(i1-),q(i2+),V(i3+),V(i4-),V(i5+)) containing the quark 
c   box loop, where V=gluon or photon

      implicit NONE
      include 'pncommon.f'

      integer i1,i2,i3,i4,i5
      double complex AQ3, Fs, Ff, fL2, fLs1, fLs3, cr

      Fs = 2d0*Tr(i1,i3)*Tr(i3,i4)*Tr(i4,i5)*Tr(i1,i5)*Sq(i3,i5)**4
     >         /Tr(i1,i2)
     >         *fLs3(cr(S(i3,i4),S(i1,i2)),cr(S(i4,i5),S(i1,i2)))
     >         / S(i1,i2)**4  
     >   + 2d0/3d0*Tr(i1,i4)*Tr(i1,i5)*Tr(i3,i4)*Sq(i3,i5)**3
     >         /Tr(i1,i2)*fL2(cr(S(i1,i2),S(i3,i4)))/S(i3,i4)**3
     >   + 2d0/3d0*Tr(i1,i3)*Tr(i1,i4)*Tr(i4,i5)*Sq(i3,i5)**3
     >         /Tr(i1,i2)*fL2(cr(S(i1,i2),S(i4,i5)))/S(i4,i5)**3
     >   - 1d0/3d0*Tr(i1,i3)*Tr(i1,i5)*Sq(i3,i5)**4
     >           /Tr(i1,i2)/Sq(i3,i4)/Sq(i4,i5)/S(i1,i2)**2

      Ff = - Tr(i1,i4)**2*Sq(i3,i5)**2/Tr(i1,i2)*
     >  fLs1(cr(S(i3,i4),S(i1,i2)),cr(S(i4,i5),S(i1,i2))) / S(i1,i2)**2

      AQ3=-(Fs+Ff)/(16d0*pi**2)

      end !AQ3

C----------------------------------------------------------------------------
      function AQ4(i1,i2,i3,i4,i5)
C----------------------------------------------------------------------------
cpn  A finite (free of infrared poles) part of the primitive amplitude
c    A_5(qbar(i1-),q(i2+),V(i3+),V(i4+),V(i5-)) containing the quark 
c    box loop, where V=gluon or photon

      implicit NONE
      include 'pncommon.f'

      integer i1,i2,i3,i4,i5
      double complex AQ4, Fs, Ff, fL0, fL1, fL2, cr

      Fs = 1d0/3d0*(
     >  -Tr(i1,i3)**2*Tr(i4,i5)**2*Sq(i3,i4)**3/Tr(i1,i2)/Tr(i3,i4)
     >       *2d0*fL2(cr(S(i4,i5),S(i1,i2)))/S(i1,i2)**3
     >  + Tr(i1,i3)*Tr(i4,i5)*Tr(i1,i5)*Sq(i3,i4)**2/Tr(i1,i2)
     >       /Tr(i3,i4)*3d0*fL1(cr(S(i4,i5),S(i1,i2)))/S(i1,i2)**2 
     >  + 2d0*Tr(i1,i5)**2*Sq(i3,i4)/Tr(i1,i2)**2
     >       /Sq(i1,i2)/Tr(i3,i4)
     >  + Tr(i1,i5)*Sq(i1,i3)*Sq(i2,i4)   
     >       /Tr(i1,i2)/Sq(i1,i2)/Tr(i3,i4)/Sq(i1,i5) 
     >  + Tr(i1,i4)*Sq(i1,i4)*Sq(i2,i4)*Sq(i3,i4)
     >      /Tr(i1,i2)/Sq(i1,i2)/Tr(i3,i4)/Sq(i4,i5)/Sq(i1,i5)  )

      Ff = -Tr(i1,i5)**2*Sq(i3,i4)/Tr(i1,i2)/Tr(i3,i4)
     >  *fL0(cr(S(i1,i2),S(i4,i5))) / S(i4,i5)    


      AQ4=-(Fs+Ff)/(16d0*pi**2)


      end !AQ4


C----------------------------------------------------------------------------
      function AQ5(i1,i2,i3,i4,i5)
C----------------------------------------------------------------------------
cpn  A primitive amplitude conjugate to AQ1(i1,i2,i3,i4,i5), obtained 
c by interchanging <ij> <-> [ji] = - [ij]

      implicit NONE
      include 'pncommon.f'

      integer i1,i2,i3,i4,i5
      double complex AQ5, Fs

      
      Fs = - ie/3d0*(
     >          Sq(i1,i3)*Tr(i3,i4)*Sq(i4,i1)**2/
     >          Sq(i1,i2)/Sq(i3,i4)**2/Sq(i4,i5)/Sq(i5,i1)
     >  +       Sq(i1,i4)*Sq(i2,i4)*Tr(i4,i5)*Sq(i5,i1)/
     >          Sq(i1,i2)/Sq(i2,i3)/Sq(i3,i4)/Sq(i4,i5)**2
     >  +       Tr(i2,i3)*Tr(i2,i5)/Tr(i1,i2)/Sq(i3,i4)/Sq(i4,i5) )


      AQ5 = + Fs/16d0/pi**2

      end !AQ5


C----------------------------------------------------------------------------
      function AQ6(i1,i2,i3,i4,i5)
C----------------------------------------------------------------------------
cpn  A primitive amplitude conjugate to AQ2(i1,i2,i3,i4,i5), obtained
c    by interchanging <ij> <-> [ji] = - [ij]

      implicit NONE
      include 'pncommon.f'

      integer i1,i2,i3,i4,i5
      double complex AQ6, Fs, Ff, cr, fL0, fL1, fL2

      Fs = 1d0/3d0*(
     >     2d0*Sq(i1,i5)*Tr(i2,i5)*Sq(i3,i4)*Sq(i3,i5)*Tr(i4,i5)**2
     >          /Sq(i4,i5)*fL2(cr(S(i1,i2),S(i3,i4)))/ S(i3,i4)**3 
     >  -  Sq(i1,i3)*Sq(i1,i5)*Sq(i3,i4)*Tr(i4,i5)**2/Sq(i1,i2)
     >          /Sq(i4,i5)*fL1(cr(S(i1,i2),S(i3,i4)))/ S(i3,i4)**2 
     >  -  Sq(i1,i3)*Tr(i2,i4)*Tr(i4,i5)
     >          /Sq(i1,i2)/Tr(i1,i2)/Tr(i3,i4)/Sq(i4,i5) 
     >  +  Tr(i2,i4)**2*Tr(i2,i5)/Tr(i1,i2)/Tr(i2,i3)/Tr(i3,i4)
     >         /Sq(i4,i5) )

      Ff = - Sq(i1,i3)**2*Tr(i4,i5)/ Sq(i1,i2)/Sq(i4,i5)
     >  *fL0(cr(S(i1,i2),S(i3,i4))) / S(i3,i4) 

      AQ6 = +(Fs+Ff)/(16d0*pi**2)
      
      end !AQ6



C----------------------------------------------------------------------------
      function AQ7(i1,i2,i3,i4,i5)
C----------------------------------------------------------------------------
cpn  A primitive amplitude conjugate to AQ3(i1,i2,i3,i4,i5), obtained
c by the replacement <ij> <-> [ji] = - [ij]


      implicit NONE
      include 'pncommon.f'

      integer i1,i2,i3,i4,i5
      double complex AQ7, Fs, Ff, fL2, fLs1, fLs3, cr

      Fs = 2d0*Sq(i1,i3)*Sq(i3,i4)*Sq(i4,i5)*Sq(i1,i5)*Tr(i3,i5)**4
     >         /Sq(i1,i2)
     >         *fLs3(cr(S(i3,i4),S(i1,i2)),cr(S(i4,i5),S(i1,i2)))
     >         / S(i1,i2)**4  
     >   + 2d0/3d0*Sq(i1,i4)*Sq(i1,i5)*Sq(i3,i4)*Tr(i3,i5)**3
     >         /Sq(i1,i2)*fL2(cr(S(i1,i2),S(i3,i4)))/S(i3,i4)**3
     >   + 2d0/3d0*Sq(i1,i3)*Sq(i1,i4)*Sq(i4,i5)*Tr(i3,i5)**3
     >         /Sq(i1,i2)*fL2(cr(S(i1,i2),S(i4,i5)))/S(i4,i5)**3
     >   - 1d0/3d0*Sq(i1,i3)*Sq(i1,i5)*Tr(i3,i5)**4
     >           /Sq(i1,i2)/Tr(i3,i4)/Tr(i4,i5)/S(i1,i2)**2

      Ff = - Sq(i1,i4)**2*Tr(i3,i5)**2/Sq(i1,i2)*
     >  fLs1(cr(S(i3,i4),S(i1,i2)),cr(S(i4,i5),S(i1,i2))) / S(i1,i2)**2


      AQ7 = +(Fs+Ff)/(16d0*pi**2)

      end !AQ7

C----------------------------------------------------------------------------
      function AQ8(i1,i2,i3,i4,i5)
C----------------------------------------------------------------------------
cpn  A primitive amplitude conjugate to AQ4(i1,i2,i3,i4,i5), obtained
c    by interchanging <ij> <-> [ji] = - [ij]

      implicit NONE
      include 'pncommon.f'

      integer i1,i2,i3,i4,i5
      double complex AQ8, Fs, Ff, fL0, fL1, fL2, cr

      Fs = 1d0/3d0*(
     >  -Sq(i1,i3)**2*Sq(i4,i5)**2*Tr(i3,i4)**3/Sq(i1,i2)/Sq(i3,i4)
     >       *2d0*fL2(cr(S(i4,i5),S(i1,i2)))/S(i1,i2)**3
     >  + Sq(i1,i3)*Sq(i4,i5)*Sq(i1,i5)*Tr(i3,i4)**2/Sq(i1,i2)
     >       /Sq(i3,i4)*3d0*fL1(cr(S(i4,i5),S(i1,i2)))/S(i1,i2)**2 
     >  + 2d0*Sq(i1,i5)**2*Tr(i3,i4)/Sq(i1,i2)**2
     >       /Tr(i1,i2)/Sq(i3,i4)
     >  + Sq(i1,i5)*Tr(i1,i3)*Tr(i2,i4)   
     >       /Sq(i1,i2)/Tr(i1,i2)/Sq(i3,i4)/Tr(i1,i5) 
     >  + Sq(i1,i4)*Tr(i1,i4)*Tr(i2,i4)*Tr(i3,i4)
     >      /Sq(i1,i2)/Tr(i1,i2)/Sq(i3,i4)/Tr(i4,i5)/Tr(i1,i5)  )

      Ff = -Sq(i1,i5)**2*Tr(i3,i4)/Sq(i1,i2)/Sq(i3,i4)
     >  *fL0(cr(S(i1,i2),S(i4,i5))) / S(i4,i5)    

      AQ8 = +(Fs+Ff)/(16d0*pi**2)


      end !AQ8




C The functions AAQG1 and AAQG2 are q+qbar+g+gamma+gamma amplitudes
c taken from 2002 Yasui's paper, and they may be wrong
c$$$C----------------------------------------------------------------------------
c$$$      function AAQG1(ipar,i1,i2,i3,i4,i5)
c$$$C----------------------------------------------------------------------------
c$$$C     Primitive amplitude for q(i1+)+qbar(i2-)+g(i3+)+gamma(i4+)
c$$$cpn2006 +gamma(i5+) -> 0
c$$$
c$$$      implicit NONE
c$$$      include 'pncommon.f'
c$$$
c$$$      integer i1,i2,i3,i4,i5, ipar, i, j
c$$$      double complex AAQG1,feps,temp, Tr1(10,10), Sq1(10,10)
c$$$
c$$$cpn2006 ipar=1 (2) -- +++++ or ----- (parity-conjugated amplitude)
c$$$      if (ipar.eq.1) then
c$$$        do i=1,nleg
c$$$          do j=1,nleg
c$$$             Tr1(i,j) = Tr(i,j)
c$$$             Sq1(i,j) = Sq(i,j)
c$$$          enddo 
c$$$        enddo 
c$$$      elseif (ipar.eq.2) then
c$$$        do i=1, nleg
c$$$          do j=1, nleg
c$$$             Tr1(i,j) = Sq(j,i)
c$$$             Sq1(i,j) = Tr(j,i)
c$$$           enddo 
c$$$         enddo                  
c$$$      endif                     !ipar
c$$$
c$$$c$$$      if (ipar.eq.2) then
c$$$c$$$        do i=1,nleg
c$$$c$$$          do j=i, nleg
c$$$c$$$            print *,i,j,Tr1(i,j)*Sq1(j,i),S(i,j)
c$$$c$$$          enddo 
c$$$c$$$        enddo 
c$$$c$$$      endif 
c$$$
c$$$      temp=-2d0*ie*( 
c$$$     >              Tr1(i4,i5)*Sq1(i5,i1)*Tr1(i1,i2)*Sq1(i2,i4)*S(i1,i4)
c$$$cpn2006 Original Yasui's expression
c$$$c     >             -Tr1(i4,i5)*Sq1(i5,i1)*Tr1(i1,i2)*Sq1(i2,i4)*S(i2,i4)
c$$$cpn2006 Expression with the fixed second term
c$$$     >             +Sq1(i4,i5)*Tr1(i5,i2)*Sq1(i2,i1)*Tr1(i1,i4)*S(i2,i4)
c$$$     >              -S(i2,i5)*S(i2,i4)*S(i1,i4) + S(i5,i1)*S(i2,i4)**2)/
c$$$     > (S(i1,i2)*Sq1(i2,i4)*Tr1(i4,i1)*Tr1(i3,i5)*Tr1(i5,i4)*Tr1(i4,i3))
c$$$
c$$$      AAQG1=temp
c$$$      
c$$$      return 
c$$$      End !AAQG1
c$$$
c$$$C---------------------------------------------------------------------------
c$$$      function AAQG2(ipar,i1,i2,i3,i4,i5)
c$$$C----------------------------------------------------------------------------
c$$$Cpn     Primitive amplitude for 
c$$$c       q(i1+)+qbar(i2-)+g(i3+)+gamma(i4-) +gamma(i5+) -> 0 or  
c$$$c       q(i1+)+qbar(i2-)+gamma(i3+)+g(i4-) +gamma(i5+) -> 0
c$$$
c$$$      implicit NONE
c$$$      include 'pncommon.f'
c$$$
c$$$      integer i1,i2,i3,i4,i5, ipar, i,j
c$$$      double complex AAQG2,feps,temp, Rln23412,regL, cr, fLi2,s1,s2, s3, 
c$$$     >  temp01, temp02,temp03,temp04, temp05, ru, Tr1(10,10), Sq1(10,10)
c$$$      data ru /(1d0,0d0)/
c$$$      save
c$$$
c$$$
c$$$cpn2006 ipar=1 (2) -- +++++ or ----- (parity-conjugated amplitude)
c$$$      if (ipar.eq.1) then
c$$$        do i=1,nleg
c$$$          do j=1,nleg
c$$$             Tr1(i,j) = Tr(i,j)
c$$$             Sq1(i,j) = Sq(i,j)
c$$$          enddo 
c$$$        enddo 
c$$$      elseif (ipar.eq.2) then
c$$$        do i=1, nleg
c$$$          do j=1, nleg
c$$$             Tr1(i,j) = Sq(j,i)
c$$$             Sq1(i,j) = Tr(j,i)
c$$$           enddo 
c$$$         enddo                  
c$$$      endif                     !ipar
c$$$
c$$$      Rln23412=cdlog(cr(S(i3,i4),S(i1,i2)))/(S(i3,i4)-S(i1,i2))**2
c$$$      regL=fLi2(ru-cr(S(i3,i4),S(i1,i2)))+
c$$$     >  fLi2(ru-cr(S(i4,i5),S(i1,i2)))+ 
c$$$     >  cdlog(cr(S(i3,i4),S(i1,i2)))*cdlog(cr(S(i4,i5),S(i1,i2)))
c$$$     >  -pi**2/6d0
c$$$
c$$$      temp01=-2d0*S(i2,i4)**2 * (-S(i1,i2)-S(i3,i5))*S(i3,i5)
c$$$     >  *Rln23412
c$$$
c$$$      temp02= 2d0*( -Tr1(i5,i1)*Sq1(i1,i2)*Tr1(i2,i3)*Sq1(i3,i5)
c$$$     >              -Tr1(i4,i5)*Sq1(i1,i4)*Tr1(i1,i3)*Sq1(i3,i5)
c$$$     >              +S(i2,i3)*S(i1,i2) )*S(i1,i4)/S(i3,i5)             +
c$$$     >        2d0*( Tr1(i1,i2)*Sq1(i2,i3)*Tr1(i3,i4)*Sq1(i1,i4)
c$$$     >              -Tr1(i5,i1)*Sq1(i1,i2)*Tr1(i2,i3)*Sq1(i3,i5) )
c$$$     >             *S(i1,i2)*S(i2,i3)/(-S(i1,i2)+S(i4,i5))/S(i3,i5)    +
c$$$     >       2d0*(-S(i2,i3)**2+S(i1,i2)*S(i5,i1)+S(i2,i3)*S(i5,i1))
c$$$     >           *S(i3,i4)**2/(-S(i1,i2)+S(i4,i5))**2                  -
c$$$     >       ( Tr1(i1,i2)*Sq1(i1,i3)*Tr1(i3,i4)*Sq1(i2,i4)
c$$$     >        -Tr1(i5,i1)*Sq1(i1,i2)*Tr1(i2,i3)*Sq1(i3,i5)
c$$$     >        -S(i3,i4)*S(i2,i3) + 2d0*S(i3,i4)*S(i5,i1) )
c$$$     >       *(S(i3,i4)*S(i1,i2)+S(i3,i4)*S(i2,i3)+S(i2,i3)*S(i1,i2) )
c$$$     >       /(-S(i1,i2)+S(i4,i5) )**2                                 -
c$$$     >      ( Tr1(i2,i3)*Sq1(i3,i4)*Tr1(i4,i5)*Sq1(i2,i5) + 
c$$$     >        Tr1(i1,i2)*Sq1(i2,i3)*Tr1(i3,i4)*Sq1(i1,i4) ) 
c$$$     >      *(-S(i3,i4)+S(i5,i1)+2d0*S(i2,i3))/(-S(i1,i2)+S(i4,i5))    -
c$$$     >      2d0*Sq1(i1,i3)*S(i1,i2)*(Tr1(i1,i2)*Tr1(i3,i4)*Sq1(i2,i4)
c$$$     >           -Tr1(i2,i3)*Sq1(i2,i5)*Tr1(i5,i1) )
c$$$     >      /(-S(i1,i2)+S(i4,i5)) -
c$$$     >      (-S(i3,i4)**2+S(i3,i4)*S(i5,i1)+2d0*S(i1,i2)**2+
c$$$     >        4d0*S(i1,i2)*S(i5,i1) )*S(i2,i3)/(-S(i1,i2)+S(i4,i5))    +
c$$$     >     2d0*(Sq1(i3,i4)*Tr1(i3,i5)*Sq1(i5,i1)*Tr1(i1,i4)+
c$$$     >          Tr1(i4,i5)*Sq1(i1,i4)*Tr1(i1,i3)*Sq1(i3,i5)+
c$$$     >          Sq1(i4,i5)*Tr1(i1,i4)*Sq1(i1,i2)*Tr1(i2,i5)-
c$$$     >          Tr1(i2,i3)*Sq1(i2,i4)*Tr1(i4,i5)*Sq1(i3,i5)  )         -
c$$$     >    4d0*(S(i2,i3)*S(i5,i1)+S(i3,i4)*S(i2,i3)-S(i2,i3)*S(i1,i2))  +
c$$$     >    2d0*S(i5,i1)*S(i4,i5) +S(i5,i1)**2+S(i3,i4)**2
c$$$     >   -2d0*S(i3,i4)*S(i1,i2)
c$$$
c$$$      temp03=cdlog(cr(S(i4,i5),S(i1,i2)))-
c$$$     >    ((S(i4,i5)-S(i1,i2))/(S(i3,i4)-S(i1,i2)))**2
c$$$     >    *cdlog(cr(S(i3,i4),S(i1,i2)))
c$$$
c$$$
c$$$      temp04=-2d0*Sq1(i3,i5)**2*Sq1(i2,i4)**2*
c$$$     >  (Tr1(i3,i4)**2*Tr1(i2,i5)**2+Tr1(i2,i3)**2*Tr1(i4,i5)**2 ) *
c$$$     >  ( regL/S(i3,i5)**2 
c$$$     >   -(2d0*S(i4,i5)-3d0*S(i1,i2)+2d0*S(i3,i4))*Rln23412/S(i3,i5))
c$$$     >    *( Tr1(i2,i3)*Sq1(i3,i4)*Tr1(i4,i5)*Sq1(i2,i5)
c$$$     >      +Tr1(i1,i2)*Sq1(i2,i3)*Tr1(i3,i4)*Sq1(i1,i4) )
c$$$     >   *S(i2,i3)*S(i3,i4)/S(i1,i2)/(-S(i1,i2)+S(i4,i5))
c$$$
c$$$      temp05=-( Tr1(i3,i4)*Sq1(i4,i5)*Tr1(i5,i1)*Sq1(i1,i3)
c$$$     >         +Tr1(i4,i5)*Sq1(i5,i1)*Tr1(i1,i2)*Sq1(i2,i4) )
c$$$     >        *S(i4,i5)*S(i5,i1)/S(i1,i2)/(S(i3,i4)-S(i1,i2))
c$$$     >       *( Tr1(i1,i2)*Sq1(i1,i3)*Tr1(i3,i4)*Sq1(i2,i4)
c$$$     >         -Tr1(i5,i1)*Sq1(i1,i2)*Tr1(i2,i3)*Sq1(i3,i5) )
c$$$     >        *(S(i3,i4)+S(i2,i3))/(-S(i1,i2)+S(i4,i5))
c$$$     >       +2d0*S(i2,i3)*S(i5,i1)*S(i3,i4)/(-S(i1,i2)+S(i4,i5))
c$$$     >       +(-Tr1(i4,i5)*Sq1(i5,i1)*Tr1(i1,i2)*Sq1(i2,i4)
c$$$     >         -Tr1(i5,i1)*Sq1(i1,i2)*Tr1(i2,i3)*Sq1(i3,i5)
c$$$     >         -S(i1,i2)*S(i4,i5) ) 
c$$$     >        *(-S(i4,i5)+S(i5,i1))/(S(i3,i4)-S(i1,i2))
c$$$     >       +2d0*S(i5,i1)*S(i4,i5)*(-S(i4,i5)+S(i2,i3)+S(i1,i2))
c$$$     >          /(S(i3,i4)-S(i1,i2))
c$$$     >       +(S(i4,i5)+S(i3,i4))*(-S(i4,i5)+S(i2,i3))
c$$$     >          *(-S(i3,i4)+S(i5,i1))/S(i1,i2)
c$$$     >       +( Tr1(i1,i2)*Sq1(i2,i3)*Tr1(i3,i4)*Sq1(i1,i4)
c$$$     >         -Tr1(i4,i5)*Sq1(i5,i1)*Tr1(i1,i2)*Sq1(i2,i4)
c$$$     >         +S(i4,i5)*S(i3,i4) )
c$$$     >        *(-S(i4,i5)-S(i3,i4)+S(i2,i3)+S(i5,i1))/S(i1,i2)
c$$$     >       +S(i3,i4)**2-2d0*S(i3,i4)*S(i5,i1)+S(i4,i5)**2
c$$$     >       +4d0*S(i3,i4)*S(i2,i3)+S(i5,i1)**2+S(i2,i3)**2
c$$$
c$$$c      temp=temp01 + temp02*temp03 + temp04 + temp05
c$$$      temp=temp05
c$$$
c$$$      temp=-ie*temp/Sq1(i2,i4)**2/Tr1(i3,i5)**2/Tr1(i1,i2)
c$$$
c$$$
c$$$      AAQG2=temp
c$$$      
c$$$      return 
c$$$      End !AAQG2


C========================================================================
C Miscalleneous 
C
C======================================================================

      FUNCTION GausInt(F,XL,XR,AERR,RERR,ERR)
C                                                   -=-=- gausint

C                                           Adptive Gaussian integration
C     
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
        DIMENSION XLIMS(100), R(93), W(93)
        INTEGER PTR(4),NORD(4)

        external f

        DATA PTR,NORD/4,10,22,46,  6,12,24,48/
        DATA R/.2386191860,.6612093865,.9324695142,
     1 .1252334085,.3678314990,.5873179543,.7699026742,.9041172563,
     1 .9815606342,.0640568929,.1911188675,.3150426797,.4337935076,
     1 .5454214714,.6480936519,.7401241916,.8200019860,.8864155270,
     1 .9382745520,.9747285560,.9951872200,.0323801710,.0970046992,
     1 .1612223561,.2247637903,.2873624873,.3487558863,.4086864820,
     1 .4669029048,.5231609747,.5772247261,.6288673968,.6778723796,
     1 .7240341309,.7671590325,.8070662040,.8435882616,.8765720203,
     1 .9058791367,.9313866907,.9529877032,.9705915925,.9841245837,
     1 .9935301723,.9987710073,.0162767488,.0488129851,.0812974955,
     1 .1136958501,.1459737146,.1780968824,.2100313105,.2417431561,
     1 .2731988126,.3043649444,.3352085229,.3656968614,.3957976498,
     1 .4254789884,.4547094222,.4834579739,.5116941772,.5393881083,
     1 .5665104186,.5930323648,.6189258401,.6441634037,.6687183100,
     1 .6925645366,.7156768123,.7380306437,.7596023411,.7803690438,
     1 .8003087441,.8194003107,.8376235112,.8549590334,.8713885059,
     1 .8868945174,.9014606353,.9150714231,.9277124567,.9393703398,
     1 .9500327178,.9596882914,.9683268285,.9759391746,.9825172636,
     1 .9880541263,.9925439003,.9959818430,.9983643759,.9996895039/
        DATA W/.4679139346,.3607615730,.1713244924,
     1 .2491470458,.2334925365,.2031674267,.1600783285,.1069393260,
     1 .0471753364,.1279381953,.1258374563,.1216704729,.1155056681,
     1 .1074442701,.0976186521,.0861901615,.0733464814,.0592985849,
     1 .0442774388,.0285313886,.0123412298,.0647376968,.0644661644,
     1 .0639242386,.0631141923,.0620394232,.0607044392,.0591148397,
     1 .0572772921,.0551995037,.0528901894,.0503590356,.0476166585,
     1 .0446745609,.0415450829,.0382413511,.0347772226,.0311672278,
     1 .0274265097,.0235707608,.0196161605,.0155793157,.0114772346,
     1 .0073275539,.0031533461,.0325506145,.0325161187,.0324471637,
     1 .0323438226,.0322062048,.0320344562,.0318287589,.0315893308,
     1 .0313164256,.0310103326,.0306713761,.0302999154,.0298963441,
     1 .0294610900,.0289946142,.0284974111,.0279700076,.0274129627,
     1 .0268268667,.0262123407,.0255700360,.0249006332,.0242048418,
     1 .0234833991,.0227370697,.0219666444,.0211729399,.0203567972,
     1 .0195190811,.0186606796,.0177825023,.0168854799,.0159705629,
     1 .0150387210,.0140909418,.0131282296,.0121516047,.0111621020,
     1 .0101607705,.0091486712,.0081268769,.0070964708,.0060585455,
     1 .0050142027,.0039645543,.0029107318,.0018539608,.0007967921/
        DATA TOLABS,TOLREL,NMAX/1.E-35,5.E-4,100/
C
C
        TOLABS=AERR
        TOLREL=RERR
     
        GausInt=0.
        NLIMS=2
        XLIMS(1)=XL
        XLIMS(2)=XR
C
10      AA=(XLIMS(NLIMS)-XLIMS(NLIMS-1))/2D0
        BB=(XLIMS(NLIMS)+XLIMS(NLIMS-1))/2D0
        TVAL=0.
        DO 15 I=1,3
15      TVAL=TVAL+W(I)*(F(BB+AA*R(I))+F(BB-AA*R(I)))
        TVAL=TVAL*AA
        DO 25 J=1,4
        VAL=0.
        DO 20 I=PTR(J),PTR(J)-1+NORD(J)
20      VAL=VAL+W(I)*(F(BB+AA*R(I))+F(BB-AA*R(I)))
        VAL=VAL*AA
        TOL=MAX(TOLABS,TOLREL*ABS(VAL))
        IF (ABS(TVAL-VAL).LT.TOL) THEN
                GausInt=GausInt+VAL
                NLIMS=NLIMS-2
                IF (NLIMS.NE.0) GO TO 10
                RETURN
                END IF
25      TVAL=VAL
        IF (NMAX.EQ.2) THEN
                GausInt=VAL
                RETURN
                END IF
        IF (NLIMS.GT.(NMAX-2)) THEN
                WRITE(*,50) GausInt,NMAX,BB-AA,BB+AA
                RETURN
                END IF
        XLIMS(NLIMS+1)=BB
        XLIMS(NLIMS+2)=BB+AA
        XLIMS(NLIMS)=BB
        NLIMS=NLIMS+2
        GO TO 10
C
50      FORMAT (' GausInt FAILS, GausInt,NMAX,XL,XR=',G15.7,I5,2G15.7)
C                        ****************************
        END


CPN  -----------------------------------------------------------------------
