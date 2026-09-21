CCPY March 2, 2005: Add the process (b+bbar --> HB)
C keyword: JHB
C Both C1 and B2 for (b+bbar --> HB) are implemented.

C **************************************************************************    
C         
C  MODULE: RASYM                                      Version 3.0
C  -------------       
C                   by Peter Arnold
C                   Ref:  P. Arnold and R. Kauffman, ANL-HEP-PR-90-70,
C                         to appear in Nucl. Phys. B
C
C      This Module computes the expansion of the small-pT resummation to
C  second order in alpha(QCD).  The result is returned by the main
C  routine, ADYPT.  The result is computed for W or Z production.
C  The expansion is done in the MS-bar factorization scheme for the
C  structure functions.  The MS-bar convention is also used for
C  the running of alpha(QCD).
C
C  References
C  ----------
C     C. Davies and W. Stirling, Nucl. Phys. B244 (1984) 337.
C     C. Davies, B. Webber and W. Stirling, Nucl. Phys. B256 (1985) 413.
C     C. Davies, Ph.D. Dissertation, Churchill College (1984).
C     J. Collins, D. Soper and G. Sternman, Nucl. Phys. B250 (1985) 199.
C     G. Curci, W. Furmanski and R. Petronzio, Nucl. Phys. B175 (1980) 27.
C     W. Furmanski and R. Petronzio, Phys. Lett. 97B (1980) 437.     
C
C **************************************************************************    


C      subroutine error
C      pause 'R_ASYM: programming error.'
C      end


C **
C **  Common: Ncodes, Nconv
C **  ---------------------
C **
C **    In order that all of the routines will work just as well with
C **  any flavor of quark, they pass along a code indicating the
C **  flavor.  Then, at the lowest level, when structure functions are
C **  evaluated, that code is used to decide which structure function
C **  to use.  The codes used are nu1, nub1, etc.  nu1 is the up
C **  content of the particle 1, nub1 the down content of particle 1,
C **  nu2 the up content of particle 2, etc.
C **
C **         u  : up       d  : down      st : strange
C **         c  : charm    b : bottom     g  : glue
C **         all: sum of all quarks and anti-qaurks
C **
C **  There are some arrays provided for conversions:
C **
C **       ng(n)   : returns ng1 if n belongs to particle 1 and ng2
C **                 if n belongs to particle 2.
C **
C **       nall(n) : returns nall1 if n belongs to particle 1 and nall2
C **                 if n belongs to particle 2.
C **
C **       nbar(n) : returns the anti-particle corresponding to n;
C **                 e.g. n(u1) = ub1.
C **

      block data
      
      IMPLICIT NONE
        integer nu1,nub1,nd1,ndb1,nst1,nstb1,nc1,ncb1,nb1,nbb1,ng1,nall1
     1        , nu2,nub2,nd2,ndb2,nst2,nstb2,nc2,ncb2,nb2,nbb2,ng2,nall2
        common /ncodes/
     1          nu1,nub1,nd1,ndb1,nst1,nstb1,nc1,ncb1,nb1,nbb1,ng1,nall1
     2        , nu2,nub2,nd2,ndb2,nst2,nstb2,nc2,ncb2,nb2,nbb2,ng2,nall2
        integer ng(24),nall(24),nbar(24)
        common /nconv/ nall,ng,nbar
        data nu1,nub1,nd1,ndb1,nst1,nstb1,nc1,ncb1,nb1,nbb1,ng1,nall1
     1            / 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12/
        data nu2,nub2,nd2,ndb2,nst2,nstb2,nc2,ncb2,nb2,nbb2,ng2,nall2
     1            /13,14,15,16,17,18,19,20,21,22,23,24/
        data ng   /11,11,11,11,11,11,11,11,11,11,11,11 
     1           , 23,23,23,23,23,23,23,23,23,23,23,23/
        data nall /12,12,12,12,12,12,12,12,12,12,12,12
     1           , 24,24,24,24,24,24,24,24,24,24,24,24/
        data nbar / 2, 1, 4, 3, 6, 5, 8, 7,10, 9,11,12
     1           , 14,13,16,15,18,17,20,19,22,21,23,24/
      end


C **
C **  FQ -- Returns the structure function f(x) corresponding to the
C **  flavor code n.  This is the only routine which assumes that
C **  particle 1 is a proton and particle 2 is either a proton or
C **  anti-proton (determined by Kcolld).  Structure functions are
C **  evaluated at scale passed by argument Scale.
C **    The m-codes used below indicate what type of *proton* structure
C **  function is used.
C **
C **  F -- As FQ but the scale is fixed at QM (from common).
C **

      function fq(n,x,qm)
      IMPLICIT NONE
      
      double precision fq,x,qm
      double precision xfu,xfd,xfsea,xfs,xfc,xfb,xfg
      integer n,m
      integer nu1,nub1,nd1,ndb1,nst1,nstb1,nc1,ncb1,nb1,nbb1,ng1,nall1
     >  , nu2,nub2,nd2,ndb2,nst2,nstb2,nc2,ncb2,nb2,nbb2,ng2,nall2
      common /ncodes/
     >  nu1,nub1,nd1,ndb1,nst1,nstb1,nc1,ncb1,nb1,nbb1,ng1,nall1
     >  , nu2,nub2,nd2,ndb2,nst2,nstb2,nc2,ncb2,nb2,nbb2,ng2,nall2
      integer ng(24),nall(24),nbar(24)
      common /nconv/ nall,ng,nbar
      integer mall,mtype(12), ihadron, ipdf
      integer jwtype,jwp,jwm,jz,jph,JHB
      COMMON /JCODES/ JWTYPE,jwp,jwm,jz,jph,JHB
      real*8 ecm, pdfh
      integer ibeam
      common/collider/ecm,ibeam


      external pdfh
      parameter(mall=12)
      

cpn2007 A hash between the local and CTEQ notations for PDF indices
c                 u  ub   d  db   s    sb  c  cb  d  db  g  all
      data mtype /1, -1,  2, -2,  3,   -3, 4, -4, 5, -5, 0, mall/

cpn2007 A hotfix to include W- production; the cross section for pp->W-
cpn2007 is the same as pbar+pbar->W+; hence we incorporate W- by
cpn2007 exchanging protons and antiprotons

      if (n.lt.nu2) then        !hadron 1
        if (jwtype.eq.jwm) then
          ihadron=-1
        else
          ihadron=1
        endif                   !jwtype
        m=mtype(n)   
      else                      !hadron 2
        if (jwtype.eq.jwm) then
          ihadron=-ibeam
        else
          ihadron=ibeam
        endif                   !jwtype
        m=mtype(n-nu2+1)  
      endif                     !n.lt.nu2
      
      if (m.ne.mall) then !PDF for the flavor m
        fq=pdfh(ihadron, m,x,qm)
      else                      !singlet quark PDF
        fq=0d0
        do ipdf=1,5
          fq=fq+pdfh(ihadron, ipdf,x,qm)+pdfh(ihadron, -ipdf,x,qm)
        enddo                   !ipdf
      endif

      end


      function f(n,x)
      IMPLICIT NONE

        double precision f,x,fq
        integer n
        double precision ss,tt,uu,qq,qm,qmm,qmu
      COMMON /STU/ SS,TT,UU,QQ,QM,QMM,qmu
      f = fq(n,x,qm)
      end



C **
C **  INTEG, INTEGC, INTGW1, INTGW0 -- These routines integrate Func
C **  between a (or x0) and b.  Integ and Integc use a straightforward
C **  Gaussian integration.  Intgw1 provides faster convergence for
C **  integrands that diverge like a power of log(1-x) near x=1 by
C **  changing variables to y = -log(1-x).  Intgw0 provides faster
C **  convergence for integrands that diverge like a power of log(x)
C **  near x=0 by changing to y = -log(x).  Integc and Integw1 both
C **  store the lower integration limit as x in common and take as
C **  an extra parameter a parton-type n which is also stored in common
C **  for use by the integrand.
C **    For Intgw1 and Intgw0, the upper limit of integration is
C **  infinity.  The routines cut it off based on the desired relative
C **  error Errel.  Note that if one is using non-adaptive gaussian
C **  integration (Nintvl not zero), one still needs to decrease Errel, 
C **  as well as increasing Nintvl, to get better accuracy.
C **

C !!!!!!!!!!!  Change above description. !!!!!!!!!!!

      function integ(func,a,b)
       IMPLICIT NONE

       double precision integ,func,a,b,agauss
        external func
        integer nintvl
        double precision errel
        common /err/ errel,nintvl       
      integ = agauss(func,a,b,nintvl,errel)
      end
      
      function integc(func,x0,b,n0)
      IMPLICIT NONE

        double precision integc,func,x0,b,agauss
        external func
        integer n0
        integer nintvl
        double precision errel
        common /err/ errel,nintvl       
        double precision x
        integer n
        common /convol/ x,n
      x = x0
      n = n0
      integc = agauss(func,x,b,nintvl,errel)
      end
     

      function intgw1(func,x0,b,n0)
      IMPLICIT NONE

       double precision intgw1,func,x0,b,ya,yb,zgauss,intw1b
        external func,intw1b
        integer n0
        integer nintvl
        double precision errel
        common /err/ errel,nintvl       
        double precision x
        integer n
        common /convol/ x,n
      x = x0
      n = n0
      ya = sqrt(-log(1.-x))/(1.+sqrt(-log(1.-x)))
      if (b.lt.1.d0) then
        yb = sqrt(-log(1.-b))/(1.+sqrt(-log(1.-b)))
      else
        yb = 1.d0
      endif
      intgw1 = zgauss(intw1b,ya,yb,nintvl,errel,func)
      end

      function intw1b(y,func)
      IMPLICIT NONE

        double precision intw1b,y,func,z
        external func
      z = 1 - exp(-(y/(1.-y))**2)
      if (z.lt.1.) then
        intw1b = 2.*y*(1.-z)/(1.-y)**3 * func(z)
      else
        intw1b = 0
      endif
      end

     
      function intgw0(func,a,b)
      IMPLICIT NONE

        double precision intgw0,func,a,b,ya,yb,zgauss,intw0b
        external func,intw0b
        integer nintvl
        double precision errel
        common /err/ errel,nintvl
      ya = sqrt(-log(b))/(1.+sqrt(-log(b)))
      if (a.gt.0) then
        yb = sqrt(-log(a))/(1.+sqrt(-log(a)))
      else
        yb = 1.d0
      endif       
      intgw0 = zgauss(intw0b,ya,yb,nintvl,errel,func)
      end

      function intw0b(y,func)
      IMPLICIT NONE

        double precision intw0b,y,func,z,eps
        external func
        parameter (eps=1.e-10)
      z = exp(-(y/(1.-y))**2)
      if (z.gt.eps) then
        intw0b = 2.*y*z/(1.-y)**3 * func(z)
      else
        intw0b = 0
      endif
      end

     

C **************************************************************************    
C         
C  Convolutions
C  ------------
C    The following routines convolute structure functions with
C  Altarelli-Parisi kernels and with the function C(x).
C
C  Naming conventions:
C
C         g1  :  leading order A-P function
C         g2  :  next-to-leading order A-P function; MS-bar
C         c1  :  order alpha(s)/2/pi piece of C(x); MS-bar
C
C  The structure functions to be convoluted is determined by the argument
C  n, which is the parton-code of the particle that produces the photon.
C  For example, if G1 were convoluted with structure fns. and n was
C  NU1, then the result would be the sum of PQQ convoluted with the
C  up structure fn. and PQG convoluted with the gluon structure fn.
C
C  References for c1: Davies and Stirling
C                     Collins, Soper and Sterman
C               Note: Davies' thesis has mistake on c1(qg).
C
C  References for g2: Curci, Furmanski and Petronzio
C                     Furmanski and Petronzio
C
C **************************************************************************    


C **
C **  C1F -- Convolution of C1 with structure fns.
C **

      function c1f(n,x)
      IMPLICIT NONE

        double precision c1f,c1fx,x,f,pi,qqdel,integc
        external c1fx
        integer n
        parameter (pi=3.141592654)
        integer nc,nu,nd,nud
        double precision cf,cnp,tr
        common /clebsh/ cf,cnp,tr,nc,nu,nd,nud      
        integer jwtype,jwp,jwm,jz,jph,JHB	
        common /jcodes/ jwtype,jwp,jwm,jz,jph,JHB
CCPY
      IF(JWTYPE.EQ.JHB) THEN        
        qqdel = cf*(pi**2/2.-1.)
      ELSE
        qqdel = cf*(pi**2/2.-4.)
      ENDIF
      
      c1f = integc(c1fx,x,1.d0,n) + qqdel*f(n,x)
      end
     
      function c1fx(z)
      IMPLICIT NONE

        double precision c1fx,z,xqq,xqg,f
        double precision x
        integer n
        integer ng(24),nall(24),nbar(24)
        common /nconv/ nall,ng,nbar
        common /convol/ x,n
        integer nc,nu,nd,nud
        double precision cf,cnp,tr
        common /clebsh/ cf,cnp,tr,nc,nu,nd,nud      
      xqq = cf*(1.-z) * f(n,x/z)/z
      xqg = z*(1.-z) * f(ng(n),x/z)/z
      c1fx = xqq + xqg
      end


C **
C **  G1F -- Convolution of G1 with structure fns.
C **

      function g1f(n,x)
      IMPLICIT NONE

        double precision x,g1f,g1fx,g1fy,f,qqdel,integ,integc
        external g1fx,g1fy
        integer n
      qqdel = - integ(g1fy,0.d0,x)
      g1f = integc(g1fx,x,1.d0,n) + qqdel*f(n,x)
      end

      function g1fx(z)
      IMPLICIT NONE

        double precision g1fx,g1fy,z,xqq,xqg,f
        integer ng(24),nall(24),nbar(24)
        common /nconv/ nall,ng,nbar
        double precision x
        integer n
        common /convol/ x,n
      xqq = g1fy(z) * (f(n,x/z)/z-f(n,x))
      xqg = (z**2+(1.-z)**2)/2. * f(ng(n),x/z)/z
      g1fx = xqq + xqg
      end

      function g1fy(z)
      IMPLICIT NONE
        double precision g1fy,z
        integer nc,nu,nd,nud
        double precision cf,cnp,tr
        common /clebsh/ cf,cnp,tr,nc,nu,nd,nud      
      g1fy = cf*(1.+z**2)/(1.-z)
      end


C **
C **  G1G1F -- Convolution of G1 with G1 with sturcture fns.
C **

      function g1g1f(n,x)
      IMPLICIT NONE
      double precision x,g1g1f,g1g1fx,g1g1fy,f,qqqdel,integ,intgw1,pi
      external g1g1fx,g1g1fy
      integer n
      parameter (pi=3.141592654)
      integer nc,nu,nd,nud
      double precision cf,cnp,tr
      common /clebsh/ cf,cnp,tr,nc,nu,nd,nud      
      qqqdel = cf**2*(9./4.-2.*pi**2/3.) - integ(g1g1fy,0.d0,x)
      g1g1f = intgw1(g1g1fx,x,1.d0,n) + qqqdel*f(n,x)
      end

      function g1g1fx(z)
      IMPLICIT NONE
        double precision g1g1fx,g1g1fy,z,beta0,xqqq,xqqg,xqgq,xqgg,f
        integer ng(24),nall(24),nbar(24)
        common /nconv/ nall,ng,nbar
        double precision x
        integer n
        common /convol/ x,n
        integer nc,nu,nd,nud
        double precision cf,cnp,tr
        common /clebsh/ cf,cnp,tr,nc,nu,nd,nud      
      beta0 = (11.*nc-4.*tr)/6.
      xqqq = g1g1fy(z) * (f(n,x/z)/z-f(n,x))
     1     + cf**2 * ( -4.*log(z)/(1.-z) - 2.*(1.-z)
     2         + (1.+z)*(3.*log(z)-4.*log(1.-z)-3.) ) * f(n,x/z)/z
      xqqg = cf * ( (z**2+(1.-z)**2)*log((1.-z)/z) - (z-1./2.)*log(z)
     1         + z - 1./4. ) * f(ng(n),x/z)/z
      xqgq = cf * ( 2./3./z + (1.+z)*log(z) - 2.*z**2/3. - z/2.
     1         + 1./2. ) * f(nall(n),x/z)/z
      xqgg = ( 2.*nc * ( 1./3./z + (z**2-z+1./2.)*log(1.-z)
     1         + (2.*z+1./2.)*log(z) + 1./4. + 2.*z - 31.*z**2/12. )
     2       + beta0 * (z**2+(1.-z)**2)/2. )  *  f(ng(n),x/z)/z
      g1g1fx = xqqq + xqqg + xqgq + xqgg
      end

      function g1g1fy(z)
      IMPLICIT NONE
        double precision g1g1fy,z
        integer nc,nu,nd,nud
        double precision cf,cnp,tr
        common /clebsh/ cf,cnp,tr,nc,nu,nd,nud      
      g1g1fy = cf**2 * (8.*log(1.-z)/(1.-z) + 6./(1-z))
      end


C **
C **  C1G1F -- Convolution of C1 with G1 with structure fns.
C **

      function c1g1f(n,x)
      IMPLICIT NONE
        double precision x,c1g1f,c1g1fx,c1g1fy,f,qqqdel,integ,intgw1
        external c1g1fx,c1g1fy
        integer n
      qqqdel = - integ(c1g1fy,0.d0,x)
      c1g1f = intgw1(c1g1fx,x,1.d0,n) + qqqdel*f(n,x)
      end

      function c1g1fx(z)
      IMPLICIT NONE
        double precision c1g1fx,c1g1fy,z,beta0,xqqq,xqqg,xqgq,xqgg,f,pi
        parameter (pi=3.141592654)
        integer ng(24),nall(24),nbar(24)
        common /nconv/ nall,ng,nbar
        double precision x
        integer n
        common /convol/ x,n
        integer nc,nu,nd,nud
        double precision cf,cnp,tr
        common /clebsh/ cf,cnp,tr,nc,nu,nd,nud      
      beta0 = (11.*nc-4.*tr)/6.
      xqqq = c1g1fy(z) * (f(n,x/z)/z-f(n,x))
     1     + cf**2 * (1.-z)*(2.*log(1.-z)-log(z)-1./2.) * f(n,x/z)/z
      xqqg = cf * ( (pi**2/2.-4.)*(z**2+(1.-z)**2)/2.
     1         - 1. + z/2. + z**2/2. - (z+1./2.)*log(z) )
     2       * f(ng(n),x/z)/z
      xqgq = cf * (1./3./z - 1. + 2.*z**2/3. - z*log(z))
     1       * f(nall(n),x/z)/z
      xqgg = ( nc * (2.*z*(1.-z)*log(1.-z) - 4.*z*log(z) + 1./3./z
     1         - 1. - 5.*z + 17.*z**2/3.) + beta0*z*(1.-z) )
     2         * f(ng(n),x/z)/z
      c1g1fx = xqqq + xqqg + xqgq + xqgg
      end

      function c1g1fy(z)
      IMPLICIT NONE
        double precision c1g1fy,z,pi
        parameter (pi=3.141592654)
        integer nc,nu,nd,nud
        double precision cf,cnp,tr
        common /clebsh/ cf,cnp,tr,nc,nu,nd,nud      
        integer jwtype,jwp,jwm,jz,jph,JHB	
        common /jcodes/ jwtype,jwp,jwm,jz,jph,JHB
CCPY
      IF(JWTYPE.EQ.JHB) THEN        
        c1g1fy = cf*(pi**2/2.-1.) * cf*(1.+z**2)/(1.-z)
      else
        c1g1fy = cf*(pi**2/2.-4.) * cf*(1.+z**2)/(1.-z)
      endif
      
      end


C **
C **  G2F -- Convolution of G2 with structure fns.
C **

      function g2f(n,x)
      IMPLICIT NONE
        double precision g2f,x,g2pqq,g2pqqb,g2pff,g2pfg,f
        integer n
        integer ng(24),nall(24),nbar(24)
        common /nconv/ nall,ng,nbar
        integer nc,nu,nd,nud
        double precision cf,cnp,tr
        common /clebsh/ cf,cnp,tr,nc,nu,nd,nud      
      g2f = g2pqq(n,x) + g2pqqb(nbar(n),x)
     1    + (g2pff(nall(n),x)-g2pqq(nall(n),x)-g2pqqb(nall(n),x))
     2       / (4.*tr)
     2    + g2pfg(ng(n),x) / (4.*tr)
      end

C * Convolution of next-order, non-singlet PQQ with structure fns.

      function g2pqq(n,x)
      IMPLICIT NONE
        double precision g2pqq,x,qqdel,intgw0,intgw1,p2qq,p2qqb,p2qqx,f
        external p2qq,p2qqb,p2qqx
        integer n
      qqdel = intgw0(p2qq,0.d0,x) - intgw0(p2qqb,0.d0,1.d0)
      g2pqq = intgw1(p2qqx,x,1.d0,n) - qqdel*f(n,x)
      end

      function p2qqx(z)
        double precision p2qqx,p2qq,z,f
        double precision x
        integer n
        common /convol/ x,n
      p2qqx = p2qq(z) * (f(n,x/z)/z-f(n,x))
      end      

C * Convolution of next-order, non-singlet PQQB with structure fns.

      function g2pqqb(n,x)
      IMPLICIT NONE
        double precision g2pqqb,x,qqbdel,integ,intgw1,p2qqb,p2qqbx,f
        external p2qqb,p2qqbx
        integer n
      qqbdel = - integ(p2qqb,x,1.d0)
      g2pqqb = intgw1(p2qqbx,x,1.d0,n) - qqbdel*f(n,x)
      end

      function p2qqbx(z)
      IMPLICIT NONE
        double precision p2qqbx,p2qqb,z,f
        double precision x
        integer n
        common /convol/ x,n
      p2qqbx = p2qqb(z) * (f(n,x/z)/z-f(n,x))
      end      

C * Convolution of next-order, singlet PFF with structure fns.

      function g2pff(nnall,x)
      IMPLICIT NONE
        double precision g2pff,x,ffdel,intgw0,intgw1,p2zgf,p2zffx,f
        double precision integ,p2zff
        external p2zgf,p2zff,p2zffx
        integer nnall
      ffdel = intgw0(p2zgf,0.d0,0.5d0) + intgw1(p2zgf,0.5d0,1.0d0,-1)
     1      + integ(p2zff,0.d0,x)
      g2pff = intgw1(p2zffx,x,1.d0,nnall) - ffdel*f(nnall,x)
      end

      function p2zffx(z)
      IMPLICIT NONE
        double precision p2zffx,p2ff,z,f
        double precision x
        integer n
        common /convol/ x,n
      p2zffx = z*p2ff(z) * (f(n,x/z)/z**2-f(n,x))
      end

      function p2zff(z)
      IMPLICIT NONE
        double precision p2zff,p2ff,z
      p2zff = z*p2ff(z)
      end

C *  z times next-order, singlet PGF.

      function p2zgf(z)
      IMPLICIT NONE
        double precision p2zgf,z,p2gf
      p2zgf = z*p2gf(z)
      end

C * Convolution of next-order, singlet PFG with structure fns.

      function g2pfg(nng,x)
      IMPLICIT NONE
        double precision g2pfg,x,intgw1,p2fgx
        external p2fgx
        integer nng
      g2pfg = intgw1(p2fgx,x,1.d0,nng)
      end

      function p2fgx(z)
      IMPLICIT NONE
        double precision p2fgx,p2fg,z,f
        double precision x
        integer n
        common /convol/ x,n
      p2fgx = p2fg(z) * f(n,x/z)/z
      end      


C **
C **  P2QQ, P2QQB, P2FF, P2FG, P2GF -- Second order Altarelli-Parisi
C **  kernels for non-singlet (P2QQ, P2QQB) and singlet (P2FF, P2FG,
C **  P2GF) cases.  The factorization convention is MS-bar.
C **

      function p2qq(z)
      IMPLICIT NONE
        double precision p2qq,z,pi,wz,w1mz,pf,pg,pnf
        parameter (pi=3.141592654)
        integer nc,nu,nd,nud
        double precision cf,cnp,tr
        common /clebsh/ cf,cnp,tr,nc,nu,nd,nud      
      wz = log(z)
      w1mz = log(1.-z)
      pf = - 2.*(1.+z**2)/(1.-z)*wz*w1mz - (3./(1.-z)+2.*z)*wz
     1     - (1.+z)*wz**2/2. - 5.*(1.-z)
      pg = (1.+z**2)/(1.-z) * (wz**2+11.*wz/3.+67./9.-pi**2/3.)
     1     + 2.*(1.+z)*wz + 40.*(1.-z)/3.
      pnf = 2./3. * ( (1.+z**2)/(1.-z)*(-wz-5./3.) - 2.*(1.-z) )
      p2qq = cf**2*pf + cf*nc*pg/2. + cf*tr*pnf
      end

      function p2qqb(z)
      IMPLICIT NONE
        double precision p2qqb,z,pi,wz,pa,ddilog
        parameter (pi=3.141592654)
        integer nc,nu,nd,nud
        double precision cf,cnp,tr
        common /clebsh/ cf,cnp,tr,nc,nu,nd,nud      
      wz = log(z)
      pa = 2.*(1.+z)*wz + 4.*(1.-z)  +  (1.+z**2)/(1.+z)
     1     * (wz**2-4.*wz*log(1.+z)-4.*ddilog(-z)-pi**2/3.)
      p2qqb = (cf**2-cf*nc/2.) * pa
      end

      function p2ff(z)
      IMPLICIT NONE
        double precision p2ff,z,pi,wz,w1mz,pf,pg2,pnf,pff,pffmx
        double precision s2,ddilog
        parameter (pi=3.141592654)
        integer nc,nu,nd,nud
        double precision cf,cnp,tr
        common /clebsh/ cf,cnp,tr,nc,nu,nd,nud      
      wz = log(z)
      w1mz = log(1.-z)
      pff = (1.+z**2)/(1.-z)
      pffmx = (1.+z**2)/(1.+z)
      s2 = log(z)**2/2. - 2.*log(z)*log(1.+z) - 2*ddilog(-z) - pi**2/6.
      pf = - 1. + z + (1./2.-3.*z/2.)*wz - (1+z)*wz**2/2.
     1     - (3.*wz/2.+2.*wz*w1mz)*pff + 2.*pffmx*s2
      pg2 = 14.*(1.-z)/3. + (11.*wz/6.+wz**2/2.+67./18.-pi**2/6.)*pff
     1     - pffmx*s2
      pnf = - 16./3. + 40.*z/3. + (10.*z+16.*z**2/3.+2.)*wz
     1      - 112.*z**2/9. + 40./9./z - 2.*(1.+z)*wz**2
     2      - (10./9. + 2.*wz/3.)*pff
      p2ff = cf**2*pf + cf*nc*pg2 + cf*tr*pnf
      end

      function p2gf(z)
      IMPLICIT NONE
        double precision p2gf,z,pi,wz,w1mz,pf,pg2,pnf,pgf,pgfmx
        double precision s2,ddilog
        parameter (pi=3.141592654)
        integer nc,nu,nd,nud
        double precision cf,cnp,tr
        common /clebsh/ cf,cnp,tr,nc,nu,nd,nud      
      wz = log(z)
      w1mz = log(1.-z)
      pgf = (1.+(1.-z)**2)/z
      pgfmx = (1.+(1.+z)**2)/(-z)
      s2 = log(z)**2/2. - 2.*log(z)*log(1.+z) - 2*ddilog(-z) - pi**2/6.
      pf = - 5./2. - 7.*z/2. + (2.+7.*z/2.)*wz + (-1.+z/2.)*wz**2
     1     - 2.*z*w1mz + (-3.*w1mz-w1mz**2)*pgf
      pg2 = 28./9. + 65.*z/18. + 44.*z**2/9. + (-12.-5.*z-8.*z**2/3.)*wz
     1     + (4.+z)*wz**2 + 2.*z*w1mz
     2     + (-2.*wz*w1mz+wz**2/2.+11.*w1mz/3.+w1mz**2-pi**2/6.
     3        +1./2.)*pgf
     4     + pgfmx*s2
      pnf = - 4.*z/3. - (20./9.+4.*w1mz/3.)*pgf
      p2gf = cf**2*pf + cf*nc*pg2 + cf*tr*pnf
      end
      
      function p2fg(z)
      IMPLICIT NONE
        double precision p2fg,z,pi,wz,w1mz,s2,pfg,qg1,qg2,ddilog
        parameter (pi=3.141592654)
        integer nc,nu,nd,nud
        double precision cf,cnp,tr
        common /clebsh/ cf,cnp,tr,nc,nu,nd,nud      
      s2 = log(z)**2/2. - 2.*log(z)*log(1.+z) - 2*ddilog(-z) - pi**2/6.
      wz = log(z)
      w1mz = log(1.-z)
      pfg = z**2+(1.-z)**2
      qg1 = 4. - 9.*z + (-1.+4.*z)*wz + (-1.+2.*z)*wz**2 + 4.*w1mz
     1      + ( - 4.*wz*w1mz + 4.*wz + 2.*wz**2 - 4.*w1mz + 2.*w1mz**2
     2        - 2.*pi**2/3. + 10.) * pfg
      qg2 = 182./9. + 14.*z/9. + 40./9./z + (136.*z/3.-38./3.)*wz
     1      - 4.*w1mz - (2.+8.*z)*wz**2
     2      + ( - wz**2 + 44.*wz/3. - 2.*w1mz**2 + 4.*w1mz + pi**2/3.
     3        - 218./9. ) * pfg
     4      + 2.*(z**2+(1.+z)**2)*s2
      p2fg = cf*tr*qg1 + nc*tr*qg2
      return
      end



C **************************************************************************    
C         
C  High-Level Routines
C  -------------------
C
C **************************************************************************    


C **
C **  ASIG -- Returns the sum of
C **
C **        C(n,m) * (alphas(mu)/2/pi)**N * log(Q**2/qt**2)**M
C **
C **  convolved with structure fns.  The C(n,m) here are Davies'
C **  perturbative coefficients, not to be confused with the C(x)
C **  structure fn. corrections implemented earlier.
C **    The result is returned for a photon vertex involving
C **  parton types n1 and n2.
C **

      function asig(n1,n2)
      IMPLICIT NONE
        double precision asig
        double precision f,f1,f2,g1f,g1f1,g1f2,g2f,g2f1,g2f2
        double precision g1g1f,g1g1f1,g1g1f2,c1g1f,c1g1f1,c1g1f2
        double precision c1f,c1f1,c1f2,asig1,asig2
        double precision c11,c10,c23,c22,c21,c20
        double precision d12,d11,d23,d22,d21,a1,a2,b1,b2
        double precision wq2qt2,wm2q2,wmu2q2,wqt2q2,beta0,zeta3,pi
       
        double precision alpi,al2pi
        external alpi
       
        integer n1,n2
        parameter (pi=3.141592654)
        parameter (zeta3=1.20205 69031 59594 28540)
        double precision ss,tt,uu,qq,qm,qmm,qmu
      COMMON /STU/ SS,TT,UU,QQ,QM,QMM,qmu
        double precision xlim,qt,y
        common /x1lim/ xlim(2),qt,y
        integer morder
        common /mcodes/ morder
        integer nc,nu,nd,nud
        double precision cf,cnp,tr
        common /clebsh/ cf,cnp,tr,nc,nu,nd,nud      
        double precision qcdlam
        common /qcdlam/ qcdlam
        integer mord
        double precision xa,xb
        common /aparam/ xa,xb,mord
        integer jwtype,jwp,jwm,jz,jph,JHB	
        common /jcodes/ jwtype,jwp,jwm,jz,jph,JHB

!      qmu = qm
C !!! Note me !!!!

      beta0 = (11.*nc-4.*tr)/6.

      wq2qt2 = log(qq/qt**2)
      wqt2q2 = - wq2qt2
      wm2q2  = log(qm**2/qq)
      wmu2q2 = log(qmu**2/qq)

      f1 = f(n1,xa)
      f2 = f(n2,xb)
      c1f1 = c1f(n1,xa)
      c1f2 = c1f(n2,xb)
      g1f1 = g1f(n1,xa)
      g1f2 = g1f(n2,xb)

      a1 = 2.*cf
      b1 = -3.*cf
      a2 = 2.*cf * ( (67./18.-pi**2/6.)*nc - 10.*tr/9. )
      b2 = cf**2 * ( pi**2 - 3./4. - 12.*zeta3 )
     1   + cf*nc * ( 11.*pi**2/9. - 193./12. + 6.*zeta3 )
     2   + cf*tr * ( -4.*pi**2/9. + 17./3. )

CCPY
      IF(JWTYPE.EQ.JHB) B2 = B2 + BETA0*CF*6.0        

      d12 = - a1/2.
      d11 = - b1
      d23 = - beta0*a1/3.
      d22 = - a2/2. - beta0/2.*(a1*wmu2q2+b1)
      d21 = - b2 - beta0*b1*wmu2q2

      asig1 = 0
      asig2 = 0

      if (mord.eq.1) then
        c11 = - 2.*d12*f1*f2
        c10 = - d11*f1*f2 + g1f1*f2 + g1f2*f1
        asig1 = c11*wq2qt2 + c10
      endif

      if (mord.eq.2) then
        g2f1 = g2f(n1,xa)
        g2f2 = g2f(n2,xb)
        g1g1f1 = g1g1f(n1,xa)
        g1g1f2 = g1g1f(n2,xb)
        c1g1f1 = c1g1f(n1,xa)
        c1g1f2 = c1g1f(n2,xb)
        c23 = - 2.*(d12**2)*f1*f2
        c22 = 3.*d12*(g1f1*f2 + g1f2*f1) - 3.*(d12*d11+d23)*f1*f2
        c21 = - 2.*d12 * ( (c1f1-g1f1*wm2q2)*f2 + (c1f2-g1f2*wm2q2)*f1 )
     1        + 2.*d11 * (g1f1*f2 + g1f2*f1)
     2        - 2.*g1f1*g1f2 - (d11**2+2.*d22)*f1*f2
     3        - (g1g1f1-beta0*g1f1)*f2 - (g1g1f2-beta0*g1f2)*f1
        c20 = c1f1*g1f2 + c1f2*g1f1
     1        - d11 * ( (c1f1-g1f1*wm2q2)*f2 + (c1f2-g1f2*wm2q2)*f1 )
     2        - 2.*g1f1*g1f2*wm2q2
     3        + (8.*zeta3*(d12**2)-d21)*f1*f2
     4        - ( beta0*c1f1 - c1g1f1 + g1g1f1*wm2q2
     5          - beta0*g1f1*wmu2q2 - g2f1 ) * f2
     6        - ( beta0*c1f2 - c1g1f2 + g1g1f2*wm2q2
     7          - beta0*g1f2*wmu2q2 - g2f2 ) * f1
        asig2 = c23*wq2qt2**3 + c22*wq2qt2**2 + c21*wq2qt2 + c20
      endif

      al2pi = alpi(qmu)/2.0
CZL
      asig = al2pi*asig1 + al2pi**2*asig2 
      
      end


C **  CHARGE -- Takes a function SIG(n1,n2), which should return
C **  a differential cross-section for a process where the
C **  photon/W/Z couples to quark species n1/n2.  This routine
C **  sums SIG over all the appropriate species, putting in
C **  coupling constants and mixing angles.

      function charge(sig)
      IMPLICIT NONE
        double precision charge,sig,udbk,uub1,ddb1
        external sig
        integer nu1,nub1,nd1,ndb1,nst1,nstb1,nc1,ncb1,nb1,nbb1,ng1,nall1
     1        , nu2,nub2,nd2,ndb2,nst2,nstb2,nc2,ncb2,nb2,nbb2,ng2,nall2
        common /ncodes/
     1          nu1,nub1,nd1,ndb1,nst1,nstb1,nc1,ncb1,nb1,nbb1,ng1,nall1
     2        , nu2,nub2,nd2,ndb2,nst2,nstb2,nc2,ncb2,nb2,nbb2,ng2,nall2
        double precision gf,zmass,wmass,ggw,ggu,ggd,gvu,gvd,gau,gad
        common /wparm/ gf,zmass,wmass,ggw,ggu,ggd,gvu,gvd,gau,gad
        double precision cc11,cc12,cc21,cc22
        common /cabibo/ cc11,cc12,cc21,cc22
        integer jwtype,jwp,jwm,jz,jph,JHB	
        common /jcodes/ jwtype,jwp,jwm,jz,jph,JHB
	
	real*8 p1,p2
cpn2007 The flavor calls in the commented block are wrong.
cpn2007 Below is the corrected block
c$$$      if (jwtype.eq.jwp.or.jwtype.eq.jwm) then
c$$$        udbk = cc11*sig(nu1,ndb2) + cc12*sig(nu1,nstb2)
c$$$     1       + cc21*sig(nc1,ndb2) + cc22*sig(nc1,nstb2)
c$$$     2       + cc11*sig(nu2,ndb1) + cc12*sig(nu2,nstb1)
c$$$     3       + cc21*sig(nc2,ndb1) + cc22*sig(nc2,nstb1)
c$$$        charge = ggw*udbk
c$$$      else if (jwtype.eq.jz) then
c$$$        uub1 = sig(nu1,nub2) + sig(nc1,ncb2)
c$$$     1       + sig(nu2,nub1) + sig(nc2,ncb1)
c$$$        ddb1 = sig(nd1,ndb2) + sig(nst1,nstb2) + sig(nb1,nbb2)
c$$$     1       + sig(nd2,ndb1) + sig(nst2,nstb1) + sig(nb2,nbb1)
c$$$        charge = ggu*uub1 + ggd*ddb1
c$$$      else if (jwtype.eq.JHB) then
c$$$        ddb1 = sig(nb1,nbb2) + sig(nb2,nbb1)
c$$$        charge = ddb1
c$$$      else
c$$$        call error
c$$$      endif
c$$$
      if (jwtype.eq.jwp.or.jwtype.eq.jwm) then
        udbk = cc11*sig(nu1,ndb2) + cc12*sig(nu1,nstb2)
     1       + cc21*sig(nc1,ndb2) + cc22*sig(nc1,nstb2)
     2       + cc11*sig(ndb1,nu2) + cc12*sig(nstb1,nu2)
     3       + cc21*sig(ndb1,nc2) + cc22*sig(nstb1,nc2)
        charge = ggw*udbk
      else if (jwtype.eq.jz) then
        uub1 = sig(nu1,nub2) + sig(nc1,ncb2)
     1       + sig(nub1,nu2) + sig(ncb1,nc2)
        ddb1 = sig(nd1,ndb2) + sig(nst1,nstb2) + sig(nb1,nbb2)
     1       + sig(ndb1,nd2) + sig(nstb1,nst2) + sig(nbb1,nb2)
CZL
        charge = ggu*uub1 !+ ggd*ddb1
      else if (jwtype.eq.JHB) then
        ddb1 = sig(nb1,nbb2) + sig(nbb1,nb2)
        charge = ddb1
      else
        call error
      endif
  
      end


C **  ADYPT -- Returns the contribution to d(sigma)/d(y)/d(qt**2)
C **  due to the 2nd-order expansion of the small-pT resummation.
C **  If Mord is 1, the 1st-order piece is returned; if Mord is 2,
C **  the 2nd-order piece is returned.  d(sigma) is computed for
C **  the specified values of qt and y.  The result is returned
C **  in units of pb/(GeV**2).
C **    Parameters that need to be passed in common: SS, QM,
C **  QQ, JWTYPE.

      function adypt(qt0,y0,mord0)
      IMPLICIT NONE
      
        double precision adypt,qt0,y0,pi,gevxpb,charge,asig,taup
        external asig
        parameter (pi=3.141592654)
        parameter (gevxpb=3.8938568e8)
        double precision ss,tt,uu,qq,qm,qmm,qmu,ch
      COMMON /STU/ SS,TT,UU,QQ,QM,QMM,qmu
        double precision xlim,qt,y
        common /x1lim/ xlim(2),qt,y
        integer mord0,mord
        double precision xa,xb
        common /aparam/ xa,xb,mord
CCPY
      Integer KinCorr
      Common / CorrectKinematics / KinCorr

      qt = qt0
      y = y0
      mord = mord0

c      taup = sqrt((qq+qt**2)/ss) + sqrt(qt**2/ss)
C      xa = exp( y) * taup
C      xb = exp(-y) * taup
ccpy
      IF(KINCORR.EQ.1) THEN
        xa = dexp( y)*dsqrt((qq+qt**2)/ss)
        xb = dexp(-y)*dsqrt((qq+qt**2)/ss)
      ELSE
        xa = dexp( y)*dsqrt(qq/ss)
        xb = dexp(-y)*dsqrt(qq/ss)
      ENDIF
      
      if(abs(xa).ge.1.0 .or.abs(xb).ge.1.0) then    
        print*,' abs(xa).ge.1.0 .or.abs(xa).ge.1.0 !!!! '
        print*,' qt,y,mord =',qt,y,mord
        print*,' y,xa,xb,qq,ss =',y,xa,xb,qq ,ss
        adypt=0.d0
      else
        adypt = gevxpb * pi/3./ss/qt**2 * charge(asig)
      endif
          
      END
