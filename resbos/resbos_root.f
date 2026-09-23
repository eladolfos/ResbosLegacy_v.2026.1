!----------------------------------------------------------------------!
!                                                                      !
!                                             /////                    !
!           _/_/_/   _/_/_/   _/_/_/   _/_/_/    _/_/     _/_/_/       !
!          _/   _/  _/      _/        _/   _/  _/   _/  _/             !
!         _/_/_/   _/_/     _/_/_/   _/_/_/   _/    _/  _/_/_/         !
!        _/   _/  _/             _/ _/    _/ _/     _/       _/        !
!       _/    _/ _/_/_/_/ _/_/_/_/ _/_/_/_/  _/_/_/_/ _/_/_/_/         !
!                                     ///                              !
!                                   //                                !
!                                   /                                  !
!                  ------->--------*--------<-------                   !
!                                 /                                    !
!                                .|                                    !
!                               . |                                    !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%.%%|%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%                                  %%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%             ResBos_root          %%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%                                  %%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%         Version 10.05.20cpy      %%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%         Version 06.19.20cpy      %%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%         Version 11.22.16cpy      %%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%         Version 12.19.14cpy      %%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%         Version 01.16.11cpy      %%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%         Version 06.16.10cpy      %%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%         Version 06.09.08pn       %%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%                                  %%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
CCPY Dec 19, 2014: This version does not require CERNLIB, nor LHAPDF
C
C Authors                   
C -------             C. Balazs, and C.-P. Yuan
C                  balazs@hep.anl.gov & yuan@pa.msu.edu
C             with substantial modifications by Pavel Nadolsky
C                      (nadolsky@pa.msu.edu)
C
C  The input grid files are generated using a code called Legacy, which
C  was written by C. Balazs, G. Ladinsky, P. Nadolsky and C.-P. Yuan. 
C  Legacy also uses various libraries of the CTEQ Collaboration, and it is 
C  developed and maintained by P. Nadolsky (nadolsky@pa.msu.edu). 
C  For more information check http://hep.pa.msu.edu/wwwlegacy/.
C
C  The event generator was created by modifying the popular Monte Carlo 
C  program Papageno (version 3.07) written by Ian Hinchliffe.
C  With questions or comments contact: balazs@hep.anl.gov
C  For the latest version check http://www.pa.msu.edu/~balazs/ResBos
C
C Description
C -----------
C     This program computes the fully differential cross section
C     for processes 
C                          h1 h2 -> B(-> l1 l2) X
C     Here
C          h1 and h2 are hadrons (pion, proton, pBar or a nucleon),
C          B is a single vector boson (a virtual photon, W+-, Z0 or Higgs),
C            or a pair of vector bosons (presently photon- or Z-pair).
C          l1 and l2 are decay products (when applicable).
C
C     The cross section is calculated in either LO (mostly O(alpha_s^0)) or with
C     NLO QCD corrections or with resummed QCD corrections included,
C     depending on the type of the input grid.
C     The code implements the CSS resummation formalism extended to the
C     decay lepton level (PLB355 (1995) 548 and hep-ph/9704258).
C     It contains the exact matrix elements for the B production 
C     (including initial state soft gluon resummation effects) 
C     and for the decay. It contains the B width (when applicable) and full 
C     spin-correlation effects.
C
C     The structure of the Monte Carlo event
C     generator is the same as popular MC program PAPAGENO, written by
C     Ian Hinchliffe. The program uses VEGAS as integration routine
C     to evaluate the phase space integrals. It produces the relevant
C     distributions in forms of histograms.
C
C     The code is written such that the file resbosm.f contain the modules for
C     the MC calculation and the file resbos.f is a front end to resbosm.f.
C     The code is ready for use as it is, in principle only the input file
C     should be changed.
C
C     In case of an effort to customize ResBos it is desirable to change only
C     the routines in resbos.f.
C
C
C Inputs
C ------
C     The input file is assumed to be called resbos.in.
C ----------------------------------------------------------------------------
C  Format of the input file        Description of input
C ----------------------------------------------------------------------------
C  nCall,ITMX1,NCALL1,ITMX2,NCALL2, ISEED   Vegas parameters
C  GRID_FILE1                             Main data grid file name
C  GRID_FILE2                             High QT data grid file name
C  IWGT,iKFacP,iKFacY,iYG, iMG, iMM      Weighting, K-factor and matching info
C  PTl1CUT,Yl1CUT,DELR_CUT,PTl2CUT,Yl2CUT Cuts on final state particles
C  YM12MIN,YM12MAX,QTMin,QTMax, yWmin, yWmax            Cuts on V
C  MTmin, MTmax, ETmin      Minimal and maximal transverse mass, missing ET
C  YLUMI, mt, mH                          Luminosity, top and Higgs masses
C  Output                                 Type of output desired
C  qTSep,PDF,Proc,Sc,ACut,iSub,HQA        NLO PS separator & others for diphotons
C  MV, GamV                               W/Z boson mass and width
C ----------------------------------------------------------------------------
C 1,10,20000,10,100000,182019                          > # runs; test iterations (50*20000); MC events (50*100000); seed for random #
C ./grids/legacy_w321.out                              > Main data grid
C ./grids/legacy_yk.out                                > Y piece grid (same PDF or -)
C 0, 0, 0, 2, 1, 0                                     > Un/reweight,kKFacP&Y,iYG,iMG,iMM
C 0.0, -4.5, 0.0, 10000.0, 4.5                         > Cuts(1): pTlmin,ymin,DelR,pTlmax,ymax
C 66.0, 116.0, 0.0, 10000., -99.0, 99.0                > Cuts(2): QMn,QMx,QTMn,QTMx,yWmn,yWmx
C 0.0, 10000.,0.0                                      > Cuts(3):min,max transvers mass of boson, min missing ET
C 1.0, 173.5,125.0                                     > Luminosity (pb^-1), mt,mh
C GBook                                                > Output fromat
C 1.5, 903, 0, 5, 0, 0, 0, 0.0                         > qT_Sep,PDF,Proc,Sc,ACut,ACut,iSub,iHQApp,cutM45
C WW,0                                                 > hDecay, i_decay
C 80.385, 2.108                                        > W mass and width
C 91.1876, 2.498353                                    > Z mass and width
C 1, 1, 0, 0						      > EW flag, qed flag, CKM_flag, KFAC_ANGFUNC
C ----------------------------------------------------------------------------
C
C    Description and typical values for input parameters
C
C    iGen
C    Number of loops over the whole code. Use just 1 since as of 7/1/97 the
C    total cross section is not corrected by this value.
C
C    ITMX1
C    Number of preliminary MC iterations (E.g.: 5-10)
C
C    NCALL1
C    Number of MC events in preliminary iterations (E.g.: 5000-10000)
C    In the firts round (ITMX1,NCALL1) Vegas just samples the cross section
C    but does not contribute to the final MC events.
C    NCALL1/ITMX1 should be kept around 1000.
C
C    ITMX2
C    Number of main MC iterations (E.g.: 10-100)
C
C    NCALL2
C    Number of MC events in main iterations (E.g.: 10000-50000)
C    The trick is to keep the NCALL/ITMX ratio around 1000 and preferably
C    run with low NCALL2 (~10-20k) and many iterations iTMX2 (20-100), 
C    depending on how much MC events (how much precision) is desired.
C
C    ISEED  an integer seed for the random number generator
 
C    GRID_FILEx
C    To run this program, the appropriate input data grid files are needed,
C    which is the output of the vector boson production code Legacy. The
C    grids contain the resummed or the Y-piece portions of
C       d(sigma)/d(Q_T)/dy/dQ^2
C    for _polarized_ V boson.
C    The Y-piece grid is not neccesary, the code can be run with only the
C    resummed grid. In that case in the place of the Y-piece grid name the
C    user has to write "-" (a negative sign without quotes), indicating that
C    only the resummed result is desired.
C    The presently available grids are listed below.
C
C    IWGT
C    0 (default): return weighted (conventional) cross sections 
C    1: returns unweighted events
C   -1: return weighted events and save weights for the matrix element
c       and phase space in an external file weights.dat
C   -2: read the weight for the matrix element from weights.dat and
c       substitute for the current ME weight in the Monte-Carlo
c       integration. IWGT=-1 and -2 are useful for generation of events 
c       N error PDF sets using kinematically the same events for all 
c       PDF sets. See REWEIGHTING below
C
C
C    iKFacP, iKFacY
cpn2008 NOTE: this feature is not needed with some latest grids which already
C     include O(alpha_s^2) contributions from the beginning; 
c     just keep iKFacP=iKFacY=0 with such grids
C    Original description (circa 1997): The Y piece grids contain 
C    the O(alpha_s) Y regular and Perturbative pieces.
C    If the user wants to use O(alpha_s^2) Y and Pert pieces then iKFacP and
C    iKFacY should be set to 1. This feature is implemented for W and Z
C    production as of 7/1/97. (Typical: 1, 1)
C
C    iYGrid
CCPY iYGrid=1 for reading in Y-grid files with old format
CCPY iYGrid=2 for reading in Y-grid files with new format
C iYGrid=1
C    ' Q,qT,y, Singular (K*L0,K*A3), Pert. (KP*L0,KP*A3,A1,A2,A4)'
C iYGrid=2
C    ' Q,qT,y, Singular (L0,A3), Pert. (L0,A3), R_A1, R_A2, R_A4, R_A0'  
C (Typical: 2)
C
CCPY
c    iMG (iMatchOut)
c    deterimines if the matched resummed+Y grid is saved in a file (1)
c    or not (0, default).
c
C    iMM (iMatchMethod)
c    determines the procedure for the matching of the resummed 
c    and Y grids (default: 0)
c
C    PTl1CUT
C    Cut on the pT of the first final state particle
C    ( anti-lepton for W/Z/photon production ).
C    E.g.: for W+ the positron pT > PTl1CUT.
C
C    Yl1CUT
C    Cut on the rapidity of the first final state particle.
C    E.g.: for W+ the positron |y| < Yl1CUT.
C
C    PTNCUT
C    Cut on P_T of the second final state particle
C    ( lepton for W/Z/photon production ).
C    E.g.: for W+ the missing E_T (P_T of neutrino) < PTl2CUT.
C
C    DELR_CUT
C    Delta R between the two final state particles.
C    Used for double boson production only (typically: 0.4).
C    For sinle V production should be set it to zero.
C
C    YM12MIN, YM12MAX
C    Minimum and maximum invariant mass of V (limited by the input grid).
C    (For W production typically: 50, 110 [GeV]).

C    YM12MIN, YM12MAX
C    Minimum and maximum invariant mass of V (limited by the input grid).
C    (For W production typically: 40, 500 [GeV]).
C
C    QT_Min, QT_Max
C    Minimun and maximum transverse momentum of V (limited by the input grid).
C    (For W production typically: 0, 300 [GeV]).
C
C    yWmin, yWMax
C    Minimun and maximum rapidity of V in the laboratory frame
C    
C    MTmin, MTMax, ETmin: minimal and maximal transverse mass of the
c    lepton pair, minimum missing transverse energy associated with the neutrino
c    in the case of W boson production
C
C    YLUMI
C    Luminosity of the collider. E.g.: 1000 [pb^(-1)].
C    If set to 1 the code returns cross sections in pb.
C    
C    mt, mH: top quark and Higgs boson masses
C
C    Output
C    Possible choices:
C    GBook  -> Produces built in coarse histograms
C              Results in resbos.out (ASCII) can be viewed
C              without any further software.
C    HBook  -> Produces built in histograms in RESBOS.HB (binary).
C              User needs PAW to view these histograms.
C    NT -> Produces NTuple with built in variables.
C              User needs PAW to view RESBOS.NT (binary).
C    Custom -> Produces NTuple with custom variables and cuts.
C              For more on this option see the
C              usr_resbos_... routines in this file.
C              ( The user is recommended to change only these
C                routines to preserve the consistency of the
C                physics in this ResBos. )
C    qTSep
C    Used only by the NLO option. Gives the scale where the
C    LO and NLO phase spaces are separated in QT.
C
C    PDF,Proc,Sc,ACut
C    These parameters are active only for calculating perturbative contributions
C    for diphoton production.  Except ACut which is also active for the W piece.
C
C    BIPDF 
C    Selects the PDF for the O(alpha_s) h1 h2 -> photon photon process.
C    ResBos has the built in amplitude squared for the Perturbative
C    calculation, in contrast to the CSS and Asymptotic cases where
C    it uses the grids. These grids were generated for a particular
C    parton distribution, and the PDF input (which is only active
C    for the Perturbative calculation) has to match the PDF for which
C    the grid was generated. (For the specific PDF choices look into
C    FUNCTION ALPHAS).
C
C    Proc
C    Selects the partonic initial state for the above process
C    = 1 generate q qBar -> A A X only,
C    = 2 generate q g -> A A X only,   
C    = 0 generate all the above.
C    Only active for the Perturbative process.
C 
C    Sc
C    ITScale is the scale at which the above matrix-elements and PDF is
C    evaluated.
C    = 1 -> mu = wHat
C    Specific renormalization and factorization scale choices
C    for diphoton production:
C    = 5 -> mu = invariant mass of the photon pair,
C    = 6 -> mu = half invariant mass of the photon pair,
C    = 7 -> mu = (pT(photon_1)+pT(photon_2))/2.
C
C    Acut   
C    A theory cut that prevents the above matrix element from blowing up.
C    Look into SUBROUTINE PRIMARY_CUTS if you want to see what it
C    does. Just use it as I did in my samples above. It should not
C    effect the physical results.
C
C    iSub   
C    Controls the calculation of the fragmentation piece.
C    = 0 no subtraction is done
C    = 1 a collinear piece is subtracted from the qG->AAq amplitude
C        unconditionally
C    = 2 a collinear piece is subtracted from the qG->AAq amplitude
C        if the invariant mass of the collinear pair < 5 GeV
C    = 3 a collinear piece is subtracted from the qG->AAq amplitude
C        and a fragmentation piece is added 
C    = 4 a collinear piece is subtracted from the qG->AAq amplitude
C            if the invariant mass of the collinear pair < cutM45 GeV
C            and a fragmenation piece is added
C    =-1 experimental isolation at Q_T > E_T^{iso}
C    =-2 quasi-experimental isolation at Q_T > E_T^{iso} combined
C             with smooth-cone isolation at Q_T < E_T^{iso}
C     isub=12 quasi-experimental isolation at Q_T > E_T^{iso} combined
C             with smoothly turned-on subtraction at Q_T < E_T^{iso}
C     isub=13 quasi-experimental isolation for Q_T > E_T_gamma*E_T^{iso}
c             (where E_T^{iso} denotes a dimensionless ratio of hadronic and
c              photonic energy in the isolation cone),
c             combined with smoothly turned on subtracton at
c             Q_T < E_T_gamma*E_T^{iso}      
cpn2006 isub=3,4 temporarily disabled as of August 2006
C isub=-1 or -2: defined in function ConeETmax
C     iSub note: iSub=2, 12, 13 achieve Y = P - A ~ 0 at QT ~ 0 both at the 
C            Tevatron and the LHC.
C
CsB___iHQApprox = 1 and 2 for finite m_t and infinite m_t respectively
C
C    hDecayMode
C    = AA higgs is decayed into a photon pair
C    = ZZ higgs is decayed into a Z-pair
C
C    I_DECAY
C    = 0 If hDecayMode=AA then the higgs is not decayed into AA
C        If hDecayMode=ZZ then the Z's are not decayed to leptons
C    = 1 If hDecayMode=AA then the higgs is decayed into AA
C        If hDecayMode=ZZ then the Z's are decayed to leptons
C
C
C    MV, GamV: mass and width of V (GamV=-1 for the built-in width calculation)
C
C Example of an input file for the resummed option
C
c 1,5,10000,5,10000,52321                   > # runs; iterations; MC events
c res.out                                   > Main data grid
c y.out                                     > Y piece grid (same PDF or -)
c -2, 0,0,0,1,0                             >Un/reweight,kKFacP&Y,iYG,iMG,iMM
c 0, -6, 0.0,1000,6                   > Cuts(1): pTlmin,ymin,DelR,pTlmax,ymax
c 40,500., 0.1, 300.0, -5, 5          > Cuts(2):QMn,QMx,QTMn,QTMx,yWmin,yWmax
c 0, 500,0                      > Cuts(3):min,max trans. mass, min missing ET
c 1.0, 175.0, 125.0                      > Luminosity (pb^-1), mt,mH
c GBook                                 > Output fromat
c 1.4, 31, 0, 5, 1, 12, 0, 5              > qT_Sep, PDF,Proc,Sc,ACut,Sub,iHQA,cutM45
c AA,0,0                                 >hDecay, i_decay
c 91.1875,-1                             >W/Z mass and width
c
C
C Important Remarks about the Y piece and the O(alpha_s) options
C
C     The resummed grids contain only the resummed part of the resummation
C     formula (see PLB355 (1995) 549). When comparing with experiment
C     one has to include the Y piece. The O(alpha_s) calculation is impelented
C     such that the "NLO" grids contain only the finite parts of the Q_T
C     singular pieces. Again, when comparing with experiment one has to
C     include the Y piece.
C
C     The O(alpha_s) calculation is somewhat tricky. One has to run ResBos
C     twice and add the resulted histograms by hand. This is because the
C     phase space is split into 2 parts
C       2 -> 2 with V having Q_T<Q_T^Separation and
C       2 -> 3 with V having Q_T>Q_T^Separation.
C     It is done to regulate the singularity at Q_T^W = 0.
C
C     Operationally, the user runs the asymptotic piece first with the Y-piece
C     (the sample below is for p pBar -> Z(-> e+ e-) X @ Tevatron, no cuts):
C
C 1, 5,1000, 5, 20000                               > Vegas parameters
C Z0_ASY_CTEQ4M.DAT                                 > Low QT Grid
C Z0_Y_CTEQ4M.DAT                                   > Y piece
C 0                                                 > 0=Weighted 1=Unveighted
C  0.00,99.0,0.0, 0.00,99.0                         > Cuts: pT1,y1,R,pT2,y2
C 60., 110.0, 30.0                                  > Cuts: QMn, QMx, QTMax
C 1.0d3                                             > Luminosity (pb^-1)
C HBook                                             > Output
C 1.2, 50.0                                         > qT_Sep for NLO, QT match
C
C     Then the user runs the "DSI" piece without the Y-piece:
C
C 1, 5,1000, 5, 20000                               > Vegas parameters
C Z0_DSI_CTEQ4M.DAT                                 > Low QT Grid
C -                                                 > Y piece
C 0                                                 > 0=Weighted 1=Unveighted
C  0.00,99.0,0.0, 0.00,99.0                         > Cuts: pT1,y1,R,pT2,y2
C 60., 120.0, 30.0                                  > Cuts: QMn, QMx, QTMax
C 1.0d3                                             > Luminosity (pb^-1)
C HBook                                             > Output
C 1.2, 50.0                                         > qT_Sep for NLO, QT match
C
C     Then the HBooks has to be added together histogram by histogram.
C
C
C Outputs
C
C    If the output switches HBook or NTuple are used
C    a file called resbos.hb or resbos.nt is created which are binary files.
C    resbos.hb contains numerous HBooks. The format of HBooks is compatible
C    with PAW. Another file, resbos.out, summarizes the inputs and the total
C    rate for the given run.
C
C
C Files
C -----
C     Program files
C     -------------
C     resbos.f      :  this file, interface for RESBOSM.FOR.
C     resbosm.f     :  contains all the modules needed for the
C                      calculations.
C
C     Data files
C     ----------
C     resbos.in       : input file
C     Input data grids: Most grids are for current Tevatron parameters.
C     See Note in Purpose and the header of resbos.f for more details.
C
C     K-factor grids:
C     vv_Kx_yyyy.DAT  : O(alpha_s^2) K-factor grid for vv production
C    (vv = WP,WM,Z0) for the x piece (x=Y or P [for Pert]) and for PDF YYYY.
C
C     Input data grids:
C
C     see http://www.pa.msu.edu/~balazs/ResBos/ for the full list of grid files.
C
C
C Compilation
C -----------
C    This program was developed and tested on VAX, Alpha, Sun, NeXT and PC
C    machines using various FORTRAN compilers.
C    The source code is self-contained, for the purpose of the
C    calculation it is not required to link with any libraries.
C    However, the user has to link this code with the CERN libraries
C       PACKLIB and KERNLIB
C    to generate the HBook or NTuple Output.
C    An example of a linux makefile is included below:
C    --------------------------------------------------------------------------
C    # ResBos makefile for CygWin with g77 - C. Balazs, Apr. 21, 2001
C    
C    FC = g77
C    
C    FCFlags = -c -O2 -fno-automatic -fvxt 
C    LCFlags = '-oresbos.x'
C    
C    Libs = -L/cygdrive/d/Csabi/Codes/CERNLib/2003/lib -lpacklib
C    
C    resbos.x: resbos.o resbosm.o pentagon.o
C    	$(FC) $(LCFlags) resbos.o resbosm.o pentagon.o $(Libs)
C    
C    %.o: %.for
C    	$(FC) $(FCFlags) $<
C    --------------------------------------------------------------------------
C
C Reweighting
C -----------
C It is often necessary to compute resummed cross sections for
c the same kinematics of particles in each Monte-Carlo event and a 
c large number N of special PDF sets in order to estimate the
c PDF uncertainty in the cross section.  Since the kinematical 
c variables are kept the same in each event for all the PDF sets, 
c the random Monte-Carlo fluctuations cancel in the computation of the PDF
c uncertainty, and accurate results for the PDF errors can be
c obtained with a smaller number of events. The same procedure can
c be applied to reduce the effect of Monte-Carlo fluctuations on
c various cross section ratios, provided the kinematical shape of all
c cross sections in such ratios is roughly the same.
c
c In this calculation, one needs to obtain a set of N cross section
c weights for each Monte-Carlo event. The reweighting is carried out 
c in two steps. In the first step, run ResBos to compute cross sections 
c for a reference (e.g., best-fit) PDF set using iwgt=-1 option (save event
c information in an auxiliary file weights.dat). Rename the
c produced HBOOK or NTUPLE (e.g., move resbos.hbook to resbos_1.hbook). 
c
c In the second step, set iwgt=-2 (read event information 
c from weights.dat) and  run ResBos N-1 times for each of 
c the remaining PDF sets, renaming resulting resbos.hbook
c files each time (e.g., resbos.2.hbook, resbos.3.hbook, ..., 
c resbos.N.hbook). Upon the completion of the calculation, events 
c in all .hbook files will have exactly the same kinematics, 
c but the matrix element weights will be different due to the 
c variations in the PDF sets.
c
c If the output format is an ntuple, all ntuples in the second
c step contain only the cross section weights, while the
c kinematical variables are the same as in the first step. This is
c done to save disk space while generating ntuples for many PDF sets.

C
C Program History
C ---------------
C      4/1994 Created by C. Balazs and C.-P. Yuan.
C     12/1994 The code was cleaned up from unnecersary routines,
C             consistency of variables, common blocks and modules
C             has been checked.
C      1/1995 Interpolation routine has been replaced by 3DInterp.
C             Now the code can take any (for upper bound see
C             declarations of arrays) size input data grid.
C             A primitive documentation on the structure of the code
C             is now available.
C      2/1995 Photon-lepton and Z-lepton coupling has been set.
C             The output of the histograms has been 'modulized',
C             so that custom output device(s) can be implemented easily
C             and without disturbing the physics content of the code.
C             Now the code runs on both NeXT and VAX platforms.
C      3/1995 Both HBook and NTuple type of output works.
C      4/1995 GBook was installed.
C      8/1995 Custom output routine is attached.
C             Photon-lepton and Z-lepton coupling has been rewised.
C      9/1995 Couplings and constants set to be compatible with Legacy.
C     10/1995 It's checked to correctly generate the kinematics of the
C             decay products of W+, W-, Z0 and A0.
C             Legacy and ResBos out/input files are compatible.
C     12/1995 NLO process works.
C      3/1996 Double photon production implemented.
C      4/1996 Renormalization scale dependence for the gluon-gluon
C             initiated processes installed.
C      6/1996 pp initial state option is added.
C      7/1996 Interpolation at the boundary of y phase space is improved.
C      9/1996 Grids with just one Q value are accepted.
C             In its effect it is equivalent with the narrow
C             width approximation for W and Z.
C             Analytic approximation for the low qT behavior of the
C             Asymptotic piece is implemented.
C             Coupled with this the interpolation of the CSS, ASY and Y
C             pieces was changed. For these the order of interpolation is:
C             y first, then Q, finally qT. A small bug in the Y piece
C             interpolation was fixed.
C     10/1996 ReadIn routine fixed: negative values for the grid are allowed.
C             Negative weights are allowed. Negative histograms can be
C             generated for the delta-sigma piece.
C             Breit-Wigner shape (W/Z line-shapes) was changed according
C             the LEP convention.
C     11/1996 The setup of interpolation of the Delta Sigma piece (nPart = 13)
C             was changed from N_QT_INTP=1 to N_QT_INTP=4. In other words
C             now we do interpolate in the pT direction for Delta Sigma.
C      1/1997 Matching is done using the grid. In particular: for W+/- and Z0
C             Q = M_V and y = 0 are fixed. Then the approximate crossing point
C             of the CSS + Y and the K-factor * Perturbative qT distributions
C             is determined. This point is used as matching point. With a very
C             good approximation this is the same as the crossing points of
C             the Q and y integrated cross sections.
C             For photon, higgs, gluon and double gauge bosons the matching
C             point is still has to be given in the input file.
C      2/1997 Introduced iGen. The generator routine is in a loop now, so
C             that several run can be done and written into the same HBook.
C             This way the distributions are smoother since the random number
C             generation starts from a different value every time.
C             Cut on minimum of QT (of produced [lepton/photon] pair)
C             from input file.
C             The W,Z and top masses are read from the grid.
C      3/1997 Crossing points of CSS and Asymptotic pieces are determined.
C             For W/Z this is done for Q = MV and y = 0. For A0, AA and AG is
C             done for the highest Q in the grid and y = 0.
C      4/1997 D0 cuts (for W mass measurement) are included as an option.
C      5/1997 Renormalization and factorization scale dependences for the
C             built in NLO AA process (for both q qB and q G initial states)
C             are corrected. Scale choices implemented:
C             - invariant mass of the photon pair (ITSCALE = 5),
C             - half invariant mass of the photon pair (ITSCALE = 6),
C             - (pT(photon_1)+pT(photon_2))/2 (ITSCALE = 7).
C      6/1997 Q, QT and Y dependent O(alpha_s^2) K-factors for the
C             Perturbative and Y L0 (alias 'symmetric') pieces are implemented.
C             To include these K-factors two K-factor grid files are required
C             which are the output of Arnold's O(alpha_s^2) perturbative and
C             resummed code.
C             Matching is done event by event. First MatPt is determined.
C             MatPt is the lowest second crossing point of the Total and Pert
C             cross sections in the whole Q and y range. When the weight
C             is calculate for QT > MatPt the Total QT cross section is
C             tested vs. the Pert and the higher one is used in the weight.
C      7/1997 Matching is rewised. It is done at the grid level. After the CSS
C             and Y grids are read in a unified (matched grid) is created
C             in the following fashion. The Q and y coordinates
C             of the Y grid are taken as the coordinates of the matched grid and
C             the CSS piece is interpolated at those Q and y points. Once the
C             CSS + Y and Pert pieces are known for these Q and y values the
C             crossing points of the triple differential cross sections (CSS+Y
C             and Pert) are determined for each Q and y. In the knowledge of the
C             second crossing pionts between the CSS + Y and the Pert curves
C             for qT values lower than the crossing point the CSS + Y grid is
C             loaded into the L0 and A3 values of the Y grid. For qT's larger
C             than the second crossing point the Pert piece is copied there,
C             for each Q and y.
C             When determining on the weight of an event, only this matched
C             gris is interpolated.
C      9/1997 A collinear subtraction term for the q G -> A A q process is
C             included which subtracts the collinear singularity when one of
C             the photons and the final state fermion are collinear. The
C             calculation of the subtraction term is done along the line of
C             Catani-Seymour hep-ph/9607318.
C             HBooks are reorganized.
C     10/1997 Fragmentation is included for the q G -> A A q type process.
C     11/1997 Running of the EW coupling is implemented.
C     12/1997 Fixed target q Q -> A A, q G -> A A implemented.
C      1/1998 q Q, q G, Q G -> Z Z resummed, asymptotic and fixed order
C             perturbative (real emissions) process is installed and tested
C             against the NLO results.
C      2/1998 G G -> H0 G -> Z Z G process is installed.
C             The LO process was tested
C             against Papageno and results in the literature (Dicus &
C             Glover-Bi). The CSS and Perturbative
C             processes are also working. The G G -> H0 G -> Z Z G process
C             checks out against Baur-Glover.
C      4/1998 W+/- h production is installed.
C      5/1998 Kappa coefficient for G G -> H0 is added.
C      5/1999 Formal modifications are done to make ResBos compile and run
C             on Linux platforms.
C      2/2000 New type matching implemented based on the 2001 les Houches 
C             proceedings (cf. Balazs, Collins, Soper).
C      3/2000 Z-decays included for h1 h2 -> H X -> Z Z X -> l L l' L' X
C      1/2004 New Legacy grid is recognized and can be read. This allows 
C             the use of kinematics independent B^(2) coefficients for 
C             for both of the qQ and gg initiated resummed diphoton processes. 
C             The kinematics independent C^(1) coefficient for the gg->AA process
C             is also implemented. The CFG resummation scheme is installed for
C             the diphoton processes. 
C             The 3D interpolation was also somewhat improved.
C      4/2004 C(1) and B(2) coefficients for g g -> A A X are implemented.
C             B(2) coefficient for q Q -> A A X is implemented.
C      5/2004 Most of ResBosA code and modifications are fused into ResBos.
C             Especially, subroutine Standard was updated.
C      9/2004 Z mass and width are made optional inputs.
C             The photon-Z interference terms was improved.
cpn   10/2005 Modifications for AA production: an alternative phase
c             space integration AARESPH for PERT piece;
c             improved sampling of qT space; revised isolation cuts  
cpn   06/2006 various improvements related to diphoton production;
c             modified treatment of the hard-vertex correction in the
c             Catani-deFlorian-Grazzini resummation scheme    
cpn   07/2006 small modifications to read in the grids
c             with the value of the A3 coefficient recorded in the header
cpn  08-10/06 A bug in the parton luminosity fo qg contribution 
cpn           to the diphoton production corrected; 
c             Subtraction of the diphoton fragmentation corrected, 
c             options for quasi-experimental isolation and subtraction 
c             (isub=12 and 13) implemented; qg->gamma gamma q
c              contribution (via the quark box diagram) implemented;
c             an optional auxiliary Catani-Seymour subtraction for
c             initial-state singularities in qqbar->gamma gamma implemented
c   11/28/06  minor adjustments
c   01/03/07  Pavel's and C.-P.'s changes: adjustments in the GBOOK
c              /NTUPLE output, updated NNLO K-factor
cpn 09/26/07  Fixed a problem with the premature matching of the
c             resummed and PERT cross sections due to the negative
c             small-Q_T asymptotic behavior of PERT at NNLO
cpn 10/18/07  Fixed a typo in the coupling of the ZA interference term 
c             (per Csaba's email)
cpn 05/19/08  Rudimentary event reweighting is introduced, with the purpose
c             to compute cross sections for N PDF's using the same 
c             Monte-Carlo events as for the central PDF set. The output into
c             ROOT ntuples, StdHEP, and plain format is implemented
c 

CCPY MARCH 2008: ADD H-->ZZ->4 FERMION CHANNEL
CCPY June 2009:  Improve gg -> H0 calculation with HWIDTH(); adjust QTV_LOW_CUT.
C To calcaulte NLO rate of gg -> H0, we could use the 'H0' process
C with (hDecayMode.Eq.'AA') and (I_DECAY.EQ.1).
C We need to make two separate runs. 
C One run for 
C ./grids/h0_delsig111.out           > Main data grid
C -                                  > Y piece grid (same PDF or -)
C another run for 
C./grids/h0_asy111.out                 > Main data grid
C./grids/h0_y.out                      > Y piece grid (same PDF or -)
C The NLO rate is the sum of the above two runs.
C
CCPY September 2009: Clean up the codes about the flag ACUT 
C                    for diphoton processes; Correct the codes for 
C                    the case IBEAM.EQ.1 with one-side (in y) grid files.
C
CCPY June 2010: refine SW2_EFF in SUBROUTINE DY_COUPL
C
CCPY Jan 16, 2011: Restore 'GL' and 'GG' processes, which are NOT 
C exact calculations of  top quark pair productions.  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      Program ResBos

      IMPLICIT NONE

      INTEGER*4 NCALL1,ITMX1,NCALL2,ITMX2
      COMMON/VEGASI/NCALL1,ITMX1,NCALL2,ITMX2

      INTEGER I,IRUN
      Character*78 dummy

      Integer iGen
      Common / iGen / iGen

      INTEGER*4 I_SEED
      COMMON/SEED/I_SEED

      Integer lPr
      Common / PrintLevel / lPr
c htj
      character*100 jobname
      integer jlen

      if(iargc().eq.0)then
        jobname='resbos'
      else if(iargc().ge.1) then
        call getarg(1,jobname)
      endif      
      jlen=len_trim(jobname)

c      print*,jobname(1:jlen)
c      stop
CsB * I/O
CsB On Alpha: if resbos.out exists and STATUS='Unknown' then the output
C   will be appended to the old file.
C   On NeXT : if resbos.out exists and STATUS='New' then the code crashes
C   with 'existing file specified new' message.
C Open output file
      OPEN(UNIT=16,FILE=jobname(1:jlen)//'.out',FORM='FORMATTED'
     >)
!    >,    Action='Write')
!    >,    STATUS='New' )  

C Dump input file into output
      Print*, ' Read from input file --------------------------------- '
      Write(16,*) ' Input file --------------------------------------- '
      OPEN( UNIT=5,FILE=jobname(1:jlen)//'.in',FORM='FORMATTED',
c      OPEN( UNIT=5,FILE='resbos.in',FORM='FORMATTED',
     >      STATUS='OLD') !,READONLY )
      Do i=1,20
        Read(5,'(A)',Err=199,End=199) dummy
        Print*, dummy
        Write(16,*) dummy
      End Do
 199  Close(5)
      Print*, ' ------------------------------------------------------ '
      Write(16,*) ' -------------------------------------------------- '

      OPEN( UNIT=5,FILE=jobname(1:jlen)//'.in',FORM='FORMATTED',
c      OPEN( UNIT=5,FILE='resbos.in',FORM='FORMATTED',
     >      STATUS='OLD') !,READONLY )

CsB___Set Print level to highest (2 = every warning is printed)
      lPr = 2

CsB___iRun = 1:   hardwired
      IRUN = 1
      READ(5,*) iGen, ITMX1,NCALL1,ITMX2,NCALL2,i_Seed
      DO 10 I=1,IRUN
CsB *   Main process:
        CALL YUAN_MAIN
10    CONTINUE
      CLOSE(5)
      CLOSE(16)
      STOP

      End

C -------------------------------------------------------------------------
      SUBROUTINE YUAN_MAIN
C -------------------------------------------------------------------------
      IMPLICIT NONE
C     I_GTMARK: strips the trailing "> comment" off a filename read with
C     '(A)' (whole-line, unlike the rest of this codebase's list-directed
C     reads, which stop at the first blank on their own)
      INTEGER*4 I_GTMARK

      REAL*8 W,QCDL,X1X2MIN
      INTEGER*4 NJET,N_PRIMARY_PART,NF,I_STR,NPART,NEVENT,NORDER,iBeam,
     &IW,IAPPROX,I_DECAY,N_SECOND_PART
      LOGICAL I_HADRONIZE
      REAL*8 ETOTMN,PTMIN,ETAMAX,ETASEP
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
      LOGICAL SMEARING
      REAL*8 JET_SM_RANGE
      REAL*8 MAX_WGT,NEW_MAX_WGT,WGT_FACTOR
      INTEGER*4 NBINIT, iwgt
      LOGICAL CHARM_ONLY,BOTTOM_ONLY
      REAL*8 XLU,XLD,XRU,XRD,QU,QD,GFERMI,XLE,XRE,
     &XLU2,XLD2,XRU2,XRD2,QU2,QD2,XLE2,XRE2,
     &XLU4,XLD4,XRU4,XRD4,QU4,QD4,PI,
     &CURRENTW_MASS,CURRENTZ_MASS,
     &ZUP,ZDOWN,XMGW2,XMGZ2,XMGH2,XMH2,VSQ,
     &GAMT
      REAL*8 XL,XU,ACC
      INTEGER*4 NDIM,NCALL,ITMX,NPRN
      REAL*8 XI,SI,SI2,SWGT,SCHI,ALPH
      INTEGER*4 NDO,IT,
     &I_RESONANCE
      REAL*8 Q
      INTEGER*4 NCALL1,ITMX1,NCALL2,ITMX2
      REAL*8 A_SUSY,B_SUSY,XMHM,GAMHM,XMH0,GAMH0,
     &FINAL_P
      INTEGER*4 NP_TYPE,N_PARTS
      LOGICAL IBINNER
      INTEGER*4 NIN,NOUT

      COMMON/AA/W,QCDL,X1X2MIN,NJET,N_PRIMARY_PART,NF,I_STR,IBEAM,NPART,
     &  NEVENT,NORDER,N_SECOND_PART
      COMMON/CHARGE/IW
      COMMON/CUTVAL/ETOTMN,PTMIN,ETAMAX,ETASEP,IAPPROX
      COMMON/DEC_AY/I_DECAY,I_HADRONIZE
      COMMON/DETAILS/IVEGOUT,I_READIN,I_WRITEOUT,SAVE_FILE,RESTORE_FILE
      COMMON/HIGGS/HIGGS_MASS
      COMMON/PAT_SCALE/ITSCALE
      COMMON/STAND0/XMW,XMZ,XMW2,XMZ2,AEM,GAMW,GAMZ,GAMW2,GAMZ2
      COMMON/STAND1/SWS,CWS,XMBOT,XMTOP,XMB2,XMT2,GWEAK2
      COMMON/STAND5/XMCHARM,XMCHARM2,XMHIGGS,GAMH,GAMH2
      COMMON/XPART_MASS/XMASS(100),X_DECAYMASS(100)
      COMMON/MYX12/Y_M1,Y_M2,Y_EMIN
      COMMON/JET_SMEAR/JET_SM_RANGE,SMEARING
      COMMON/UNWEIGHT/MAX_WGT,NEW_MAX_WGT,WGT_FACTOR,iwgt,NBINIT
      COMMON/HEAVY_FLAGS/CHARM_ONLY,BOTTOM_ONLY
      COMMON/STAND2/XLU,XLD,XRU,XRD,QU,QD,GFERMI,XLE,XRE
      COMMON/STAND3/XLU2,XLD2,XRU2,XRD2,QU2,QD2,XLE2,XRE2
      COMMON/STAND4/XLU4,XLD4,XRU4,XRD4,QU4,QD4,PI
      COMMON/RUNNING_BOSON/CURRENTW_MASS,CURRENTZ_MASS
      COMMON/HIGGS_STUFF/ZUP,ZDOWN,XMGW2,XMGZ2,XMGH2,XMH2,VSQ
      COMMON/HTOP/GAMT
      COMMON/BVEG1/XL(50),XU(50),ACC,NDIM,NCALL,ITMX,NPRN
      COMMON/BVEG2/XI(50,50),SI,SI2,SWGT,SCHI,NDO,IT,ALPH
      COMMON/RESONANCE/I_RESONANCE
      COMMON/SCALE/Q
      COMMON/VEGASI/NCALL1,ITMX1,NCALL2,ITMX2
      COMMON/MY_SUSY/A_SUSY,B_SUSY,XMHM,GAMHM,XMH0,GAMH0
      COMMON/OUT_MOMENTA/FINAL_P(4,50),NP_TYPE(50),N_PARTS
      COMMON/NININF/IBINNER
      COMMON/EFFIC/NIN,NOUT

      INTEGER*4 I_SEED
      COMMON/SEED/I_SEED
      REAL*8 PTECUT,YECUT,DELR_CUT,PTNCUT,yNCut,YWLCUT,YWUCUT
      COMMON/MYCUT/PTECUT,YECUT,DELR_CUT,PTNCUT,yNCut,YWLCUT,YWUCUT

      REAL*8 YM12MIN,YM12MAX,QTMin,QTMax, TMmin, Tmmax, ETmin
      COMMON/MYYM12/YM12MIN,YM12MAX,QTMin,QTMax, TMmin, Tmmax, ETmin
      REAL*8 YLUMI
      COMMON/CHIST/ YLUMI

      INTEGER*4 MXHIST,NHSHFT
      LOGICAL HEXIST
      COMMON/MYHIST/MXHIST,NHSHFT,HEXIST(200)

      REAL*8 YMW_AMP,YGW_AMP
      COMMON/VMASS/YMW_AMP,YGW_AMP

CCPY USED FOR SEPERATING NLO AND NLO_SIG CONTRIBUTIONS
      REAL*8 QT_V_SEP
      COMMON/NLO/ QT_V_SEP

      CHARACTER*145 Filename, YFilename, KPFileName, KYFileName
      Common/GridName/ Filename, YFilename
      Real*8 ECMC
      Common /ECMC/ ECMC
      Character*40 Boson
      Common / Boson / Boson
      Integer N_SUD_AB,N_SUD_A,N_SUD_B,N_WIL_C,I_SUD
      Common / Sudakov / N_SUD_AB,N_SUD_A,N_SUD_B,N_WIL_C,I_SUD
      Integer nyEdge
      Common / nyEdge / nyEdge

CsB___Interpolation
C     3DInterp:
      Integer MD1,MD2,MD3, ID1,ID2,ID3, IY1,IY2,IY3
      Integer iYPiece, iReadIn
CsB___Cooridinates for the Y piece stored in XiY.
C     Values for the Y piece stored in YAY (See ReadIn).
C     The perturbative L0&A3 pieces are stored in YAY(,,,5&6)
      Real*8 X1A,X2A,X3A,YAS,YAA, X1Y,X2Y,X3Y,YAY
      Parameter (MD1 = 200, MD2 = 300, MD3 = 100)
      Dimension X1A(0:MD1), X2A(0:MD2), X3A(0:MD3)
      Dimension YAS(0:MD1,0:MD2,0:MD3), YAA(0:MD1,0:MD2,0:MD3)
      Dimension X1Y(0:MD1), X2Y(0:MD2), X3Y(0:MD3)
      Dimension YAY(0:MD1,0:MD2,0:MD3,0:7)
      Common / Arrays / X1A,X2A,X3A, YAS,YAA,  X1Y,X2Y,X3Y, YAY
      Common / Dims / iD1,iD2,iD3, iY1,iY2,iY3

      Character*8 Output
      Common / Output / Output
      Common / YPiece / iYPiece, iReadIn
      Integer iReadInK
      Common / cReadInK / iReadInK
      Integer LTOpt
      Common / Process / LTOpt
      REAL*8 QT_V_MIN, QT_V_MAX,YB_V_MAX,Q_V_MIN,Q_V_MAX,QT_V_MXY
      COMMON/HARDCUT/ QT_V_MIN, QT_V_MAX,YB_V_MAX,Q_V_MIN,
     >  Q_V_MAX,QT_V_MXY
      Real*8 QT_V_Mat
      Common / Matching / QT_V_Mat
      Integer iOneQ
      Real*8 eps
      Common / OneQ / iOneQ
      Real*8 MatPt
      Integer iGen,ii
      Common / iGen / iGen
CsB___BIPDF: Selects the built in PDF
C     ISet : Within the CTEQ PDF's selects M, L, D, etc.
      Integer BIPDF, ISet, iProc, iDummy
      Common / BuildInPDF / BIPDF, ISet
      Common / PartonProcess / iProc
      Integer ACut
      Common / AAGCut / ACut
      Integer iKFacP,iKFacY,iYGrid,iMatchOut,iMatchMethod
      Common / KFactor / iKFacP,iKFacY,iYGrid
      Common / MatchGridCom / iMatchOut
      Integer iSub
      Real*8 cutM45, SubtractionShape, dummy
      logical SubQGFlag
      Common / Subtraction / cutM45,SubtractionShape(3:4),iSub,
     >  SubQGFlag(3:4)
      Integer iHQApprox
      Common / HQApprox / iHQApprox

      Character*2 hDecayMode
      Common / hDM / hDecayMode

!     Real*8 Wwidth
!     Common / InputWwidth / Wwidth
      Real*8 XMW_IN,XMW_GRID,GAMW_IN
      COMMON/WMASS_IN/XMW_IN,XMW_GRID,GAMW_IN
      Real*8 XMZ_IN,XMZ_GRID,GAMZ_IN
      COMMON/ZMASS_IN/XMZ_IN,XMZ_GRID,GAMZ_IN
CCPY
      INTEGER IFLAV
      REAL*8 CAA_E,CAA_O,CZA_E,CZA_O,CZZ_E,CZZ_O
      COMMON/DY_IFLAV/ CAA_E,CAA_O,CZA_E,CZA_O,CZZ_E,CZZ_O,IFLAV

cpn2005 Isolation cone parameters
      integer iqqt
      real*8 DelRIso, ETmaxIso
      common/isolation/DelRIso, ETMaxIso, iqqt
CCPY
      LOGICAL Y_ONE_SIDE
      COMMON/GRID_FILE_FORM/ Y_ONE_SIDE

CCPY jAN 16, 2011
      REAL*8 ASY_TC
      COMMON/GL_ASY/ ASY_TC
CBIN Sept, 2016
      InTEGER*4 EW_FLAG,qed_FLAG,CKM_FLAG 
      COMMON/FLAG/EW_FLAG,qed_FLAG,CKM_FLAG    
CJI NOV 2020
      INTEGER*4 KFAC_ANGFUNC, KFAC_ANGFUNC_dummy
      COMMON/ANGFUNC/ KFAC_ANGFUNC

CsB___INITILIZE ALL THE VARIABLES
      CALL STARTHERE

CsB___Read the name of the Input data grid(s)
      Read(5,'(A)') Filename
      I_GTMARK = INDEX(Filename, '>')
      IF (I_GTMARK.GT.0) Filename = Filename(1:I_GTMARK-1)
      If (Filename.Eq.'AA' .Or. Filename.Eq.'AG' .Or.
     .    Filename.Eq.'ZZ' .Or. Filename.Eq.'H0') then
CsB_____This is the case when we use the built in matrix elements.
        iYPiece = 1
        Print*, ' !! Warning: Built in amplitude and PDF is used. !! '
        LTOpt = 9
CsB_____ECM is read in
C        Read(5,*) W, iBeam, xMtop, xMHiggs
        Read(5,*) W, iBeam
C        Print*, ' Ecm, iBeam =', W, iBeam
        YFilename = '-'
        Goto 300
      End If

      Read(5,'(A)') YFilename
      I_GTMARK = INDEX(YFilename, '>')
      IF (I_GTMARK.GT.0) YFilename = YFilename(1:I_GTMARK-1)
      If (YFilename.Eq.'-') then
        iYPiece = 0
        Print*, ' +---------------------------------------+ '
        Print*, ' ! Y piece contribution is NOT included. ! '
        Print*, ' +---------------------------------------+ '
      Else
        iYPiece = 1
      End If

  300 Continue
CsB___For processes with grids nPart is set acording to LTOpt which is read from
C     the grid.
C     NPART = 11   FOR RESUMMED
C     NPart = 12   for "NLO"     (Asymptotic piece)
C     NPart = 13   for "NLO_Sig" (Delta sigma piece)
C     NPART = 14   FOR LEADING ORDER
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
CsB___K factors:
C     iKFacP = 0: no K-factor applied: Pert(1,L0) is used
C     iKFacP = 1: K-factor applied: Pert(2,L0)=K(Pert)*Pert(1,L0) is used
C     iKFacY = 0: no K-factor applied: Y(1,L0) is used
C     iKFacY = 1: K-factor applied: Y(2,L0)=K(Y)*Y(1,L0) is used
C     iYGrid = 1 or 2:  read in Y-grid files with old (1) or new (2) format
C iYGrid=1
C    ' Q,qT,y, Singular (K*L0,K*A3), Pert. (KP*L0,KP*A3,A1,A2,A4)'
C iYGrid=2
C    ' Q,qT,y, Singular (L0,A3), Pert. (L0,A3), R_A1, R_A2, R_A4, R_A0'  
C
C================================
      Read(5,*) 
     > IWGT, iKFacP,iKFacY,iYGrid,iMatchOut,iMatchMethod
CCPY
      IF(iMatchMethod.NE.0) THEN
        WRITE(*,*) ' USE ONLY iMatchMethod=0 '
        iMatchMethod=0
      ENDIF

      NJET=0
      I_DECAY=0
      IAPPROX=0
      I_HADRONIZE=.FALSE.
      SMEARING=.FALSE.
      WGT_FACTOR=1.

CsB   I_SEED=733372319 ! 5/1/2004 read from input
!     Print*, ' Seed: ', i_Seed
      IVEGOUT=1

CsB___Renormalization and factorization scale choice
      ITSCALE=4

      If (LTOpt.EQ.9) then ! Processes with built-in amplitudes =============
C================================
C     NPART = 20    G  G -> photon photon X NLO, WPRESPH-like phase space
C     NPART = 21    q qB -> photon photon X NEXT-TO-LEADING ORDER
C     NPART = 22    G  G -> photon photon X NLO, standard phase space routine
C     NPART = 23    q qB -> Z Z NLO (O(alpha_s), 1 loop)
C     NPART = 24    G  G -> H0 -> Z Z NLO (O(alpha_s^3), 2 loops)

C     NPART = 50    MANIPULATING LEGACY GRIDS
C================================
        If (Filename.Eq.'AA') then
          Boson = 'AA'
          NPart = 21            !AARESPH phase space
c          NPart = 22     !AAJPH phase space
        Else If (Filename.Eq.'AG') then
          Boson = 'AG'
          NPart = 21        !AARESPH phase space
c          NPart = 22 !AAJPH phase space
        Else If (Filename.Eq.'ZZ') then
          Boson = 'ZZ'
          NPart = 23
        Else If (Filename.Eq.'H0') then
          Boson = 'H0'
          NPart = 24
        End If
      Else ! Processes requiring grids --------------------------------------
C_______Read the CSS grid.
        IReadIn = 0
        Print*, ' Opening main data grid '
        Call CheckSum(Filename)
        Call PreReadIn(Filename, 1)
        Call ReadIn(Filename)
C        Print*, ' nPart, iProc =', npart, iProc
C        PRINT*, ' Low QT grid read in succesfully.'
        QT_V_MXY = 500.0
        If (iYPiece.Eq.1) then
C_______If Y piece is used
          If (iKFacP.Eq.1) then
C___________If O(a_s^2) K factors wanted for the L0 Pert(1) piece:
            If (Boson.Eq.'W+' .or. Boson.Eq.'W-' .or.
     .          Boson.Eq.'HP') then
C              KPFileName = 'wp_prk_ct4m.dat'
               READ(5,'(A)') KPFileName
               I_GTMARK = INDEX(KPFileName, '>')
               IF (I_GTMARK.GT.0) KPFileName = KPFileName(1:I_GTMARK-1)
            Else If (Boson.Eq.'Z0' .or. Boson.Eq.'HZ' ) then
ccpy              KPFileName = 'z0_pk_ct61m.dat'
C              KPFileName = 'pert_kf_z0_tev2_ct61m_kin.dat'
               READ(5,'(A)') KPFileName
               I_GTMARK = INDEX(KPFileName, '>')
               IF (I_GTMARK.GT.0) KPFileName = KPFileName(1:I_GTMARK-1)
            Else If (Boson.Eq.'H0') then
C              KPFileName = './grids/h0_prk_ct66m.dat'
C              KPFileName = 'h0_prk_ct66m.dat'
               READ(5,'(A)') KPFileName
               I_GTMARK = INDEX(KPFileName, '>')
               IF (I_GTMARK.GT.0) KPFileName = KPFileName(1:I_GTMARK-1)
            Else
              Print*, ' K-factor grid is not available for Pert piece. '
              iKFacP = 0
              Goto 210
C              Print*, ' K-factor grids for W-bosons are used. '
C              KPFileName = 'wp_prk_ct4m.dat'
            End If

            Print*, ' K-factor grid for Pert piece is used: ',
     >      KPFileName

CCPY!CCPY  No KPFile for 'H0' yet
CCPY!            If (Boson.NE.'H0') then

              Print*, ' Opening K-factor grid for Pert(1) '
              iReadInK = 0
              Call CheckSum(KPFileName)
              Call PreReadIn(KPFileName, 3)
              Call ReadInK(KPFileName)

CCPY!            End If  ! Boson.NE.'H0'

          End If ! iKFacP.Eq.1
  210     Continue
          If (iKFacY.Eq.1) then
C___________If O(a_s^2) K factor wanted for Y(1) piece:
            If (Boson.Eq.'W+' .or. Boson.Eq.'W-' .or.
     .          Boson.Eq.'HP') then
C              KYFileName = 'wp_yk_ct4m.dat'
               READ(5,'(A)') KYFileName
               I_GTMARK = INDEX(KYFileName, '>')
               IF (I_GTMARK.GT.0) KYFileName = KYFileName(1:I_GTMARK-1)
            Else If (Boson.Eq.'Z0' .or. Boson.Eq.'HZ' ) then
C              KYFileName = 'z0_yk_ct61m.dat'
               READ(5,'(A)') KYFileName
               I_GTMARK = INDEX(KYFileName, '>')
               IF (I_GTMARK.GT.0) KYFileName = KYFileName(1:I_GTMARK-1)
            Else If (Boson.Eq.'H0') then
C              KYFileName = './grids/h0_yk_ct66m.dat'
C              KYFileName = 'h0_yk_ct66m.dat'
               READ(5,'(A)') KYFileName
               I_GTMARK = INDEX(KYFileName, '>')
               IF (I_GTMARK.GT.0) KYFileName = KYFileName(1:I_GTMARK-1)
            Else
              Print*, ' K-factor grid is not available for Y piece. '
              iKFacY = 0
              Goto 220
C              Print*, ' K-factor grids for W-bosons are used. '
C              KYFileName = 'wp_Yk_ct4m.dat'
            End If

            Print*, ' K-factor grid for Y piece is used: ',
     >      KYFileName

CCPY!CCPY  No KYFile for 'H0' yet
CCPY!            If (Boson.NE.'H0') then

              Print*, ' Opening K-factor grid for Y(1) '
              iReadInK = 1
              Call CheckSum(KYFilename)
              Call PreReadIn(KYFileName, 4)
              Call ReadInK(KYFileName)

CCPY!            End If  ! Boson.NE.'H0'

          End If ! iKFacY.Eq.1
          
  220     Continue
C_________Read in Y piece grid
          Print*, ' Opening Y grid '
          IReadIn = 1
          Call CheckSum(YFileName)
          Call PreReadIn(YFilename, 2)
          Call ReadIn(YFilename)
!          PRINT*, ' Y grid', YFilename, ' read in succesfully.'
          If (nPart.Eq.11) then
C___________For resummed case: create a unified, CSS+Y matched grid with Pert
            Print*, ' Matching CSS+Y and Pert grids '
            If (iMatchMethod.Eq.0) Call MatchGrids1
            If (iMatchMethod.Eq.1) Call MatchGrids2
          End If
        End If ! iYPiece.Eq.1
        W = ECMC
      End If ! LTOpt.EQ.9 ===================================================

CsB___Some hard-wired stuff
      PTMIN=0.0
      ETAMAX=10.0 ! ETAMAX is not in use as of 2/9/98

CsB___Read cuts on primary final-state particles
      READ(5,*)PTECUT,YECUT,DELR_CUT,PTNCUT,YnCut

CsB___Read the invariant mass (Q) and transverse momentum ranges
      READ(5,*)YM12MIN,YM12MAX,QTMin,QTMax,YWLCUT,YWUCUT

      iqqt=0; DelRiso=0; ETmaxIso=0
      if (Boson.eq.'AA'.or.Boson.eq.'AG') then
cpn2005 Read the angular size (Delta R) of the isolation cone around the
cpn2005 photon and maximal hadronic energy    
        read(5,*) iqqt, DelRiso, ETmaxIso
      else
cpn   Read the "leptonic" transverse mass and minimal missing 
cpn   transverse energy
        read(5,*) TMmin, TMmax, ETmin
      endif

CsB___Read the luminosity of the collider, top mass, higgs mass
      READ(5,*) YLUMI, xMtop, xMHiggs

CsB___Read the type of the output
      Read(5,'(A)') Output

CCPY QT_V_SEP USED FOR SEPERATING NLO AND NLO_SIG CONTRIBUTIONS
CsB  QT_V_Mat was the QT at which the switch is done from the CSS piece
C    to the perturbative piece (a simple matching). Inactive as of 1/15/97.
CsB BIPDF: Selects the built in PDF
CsB iProc=0/1/2<->all/qqB/qG contribution
CsB ITScale: type of renormalization scale = factorization scale for processes
C   with built in matrix elements and PDF's.
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
C     iSub=-1 experimental isolation at Q_T > E_T^{iso}
C     iSub=-2 quasi-experimental isolation at Q_T > E_T^{iso} combined
C             with smooth-cone isolation at Q_T < E_T^{iso}
C     isub=12 quasi-experimental isolation at Q_T > E_T^{iso} combined
C             with smoothly turned-on subtraction at Q_T < E_T^{iso}
C     isub=13 quasi-experimental isolation for Q_T > E_T_gamma*E_T^{iso}
c             (where E_T^{iso} denotes a dimensionless ratio of hadronic and
c              photonic energy in the isolation cone),
c             combined with smoothly turned on subtracton at
c             Q_T < E_T_gamma*E_T^{iso}      
cpn2006 isub=3,4 temporarily disabled as of August 2006
C isub=-1 or -2: defined in function ConeETmax
C     iSub note: iSub=2, 12, 13 achieve Y = P - A ~ 0 at QT ~ 0 both at the 
C            Tevatron and the LHC.
C
CsB___iHQApprox = 1 and 2 for finite m_t and infinite m_t respectively
      If (nPart.LE.19) then ! Processes with grids
        READ(5,*) QT_V_SEP,BIPDF,iDummy,ITScale,ACut,iSub,iHQApprox,
     ,            cutM45
      Else
        READ(5,*) QT_V_SEP,BIPDF,iProc,ITScale,ACut,iSub,iHQApprox,
     ,            cutM45
      End If

CCPY Sept 2009
      If (BOSON.EQ.'ZZ') ISUB=0
      
      if (Boson.eq.'AA') then !set up isolation and subtraction options
        if (isub.eq.3.or.isub.eq.4) then
cpn2006 computation of the fragmentation contribution temporary disabled
          print *, 'isub = 3, 4 temporarily disabled;',
     >      ' choose another value'
          call ExitGracefully
        else if (iSub.eq.12.or.iSub.eq.13) then
cpn2006  for smoothly turned-on subtraction, cutM45 is a part of
cpn      the minimal squared invariant mass of the photon and 
cpn      final-state quark excluded by isolation
          cutM45 = ETMaxIso*DelRIso**2
        endif                   !isub
      endif                     !Boson.eq.'AA'

      If (Boson.Eq.'H0' .and. iHQApprox.NE.1) then
        Print*, ' The m_top = infinity option has to be debugged! '
        Print*, ' m_top =',xMtop,' is used for gg (top loop) -> HX. '
        iHQApprox = 1
      EndIf

CCPY Jan 16, 2011
C In leading order QCD, the asymmetry of top: ASY_TC = 0.d0.
C This initialization is needed for 'GL' calculation.       
      ASY_TC=0.D0
      If (Boson.Eq.'GL') Read(5,*) ASY_TC      
      
      Read(5,'(A2,XI1)') hDecayMode,I_DECAY
      If (Boson.Eq.'H0') then
        Print*, ' hDecayMode, I_DECAY = ', hDecayMode,I_DECAY
      EndIf

CCPY SET UP COUPLINGS FOR 'DY' 
      IF(BOSON.EQ.'ZU') THEN
C THIS IS FOR UP-TYPE QUARKS IN THE INITIAL STATE
        IFLAV=1
        BOSON='DY'
      ELSEIF(BOSON.EQ.'ZD') THEN
C THIS IS FOR DOWN-TYPE QUARKS IN THE NITIAL STATE
        IFLAV=2
        BOSON='DY'
      ENDIF

CCPY!*************************************************************************
CCPY!
CCPY!CsB   For W or Z production read optional W or Z mass & width
CCPY!      IF(BOSON.EQ.'Z0' .or. BOSON.EQ.'DY')THEN
CCPY!        Read(5,*) XMZ_IN,GAMZ_IN
CCPY!        WRITE(*,*) ' XMZ_IN,GAMZ_IN =',XMZ_IN,GAMZ_IN
CCPY!        IF(XMZ_IN.NE.XMZ_GRID .and. XMZ_IN.GT.0.) THEN
CCPY!           XMZ=XMZ_IN
CCPY!           PRINT*,
CCPY!     >     'Using input Z mass, instead of the value from the grid.'     
CCPY!        ELSE
CCPY!           XMZ=XMZ_GRID
CCPY!        ENDIF        
CCPY!        XMW=XMW_GRID
CCPY!      ELSEIF(BOSON.EQ.'W+' .OR. BOSON.EQ.'W-')THEN
CCPY!        Read(5,*) XMW_IN,GAMW_IN
CCPY!        WRITE(*,*) ' XMW_IN,GAMW_IN =',XMW_IN,GAMW_IN
CCPY!CCAO!   Check the w-boson mass and use the input w-boson mass to do the
CCPY!C       following calculation
CCPY!        IF(XMW_IN.NE.XMW_GRID .and. XMW_IN.GT.0.) THEN
CCPY!           XMW=XMW_IN
CCPY!           PRINT*,
CCPY!     >     'Using input W mass, instead of the value from the grid.'     
CCPY!        ELSE
CCPY!           XMW=XMW_GRID
CCPY!        ENDIF        
CCPY!        XMZ=XMZ_GRID
CCPY!      ELSE
CCPY!        XMW=XMW_GRID
CCPY!        XMZ=XMZ_GRID
CCPY!!       Print*, ' XMW,XMZ: ', XMW,XMZ
CCPY!      ENDIF
CCPY!
CCPY!
CCPY!
C*************************************************************************
CCPY Sept , 2016  Do not use massses from Grid files
        Read(5,*) XMW_IN,GAMW_IN
        Read(5,*) XMZ_IN,GAMZ_IN
        IF (XMW_IN*GAMW_IN .LT. 0.0) THEN
           PRINT*,' XMW_IN AND GAMW_IN NEED TO BE SET IN resbos,in file'
C           CALL EXIT
        ENDIF 
        IF (XMZ_IN*GAMZ_IN .LT. 0.0) THEN
           PRINT*,' XMZ_IN AND GAMZ_IN NEED TO BE SET IN resbos,in file'
           CALL EXIT
        ENDIF 

        XMW=XMW_IN
        XMZ=XMZ_IN

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

        Read(5,*) EW_FLAG, qed_FLAG, CKM_FLAG, KFAC_ANGFUNC_dummy
        if(KFAC_ANGFUNC.ne.2) KFAC_ANGFUNC = KFAC_ANGFUNC_dummy

C EW_FLAG=3, all the parameters are read from FEWZ3.1 input file
        IF(EW_FLAG.eq.3) THEN
          PRINT*,' XMW, XMZ, GAMW_IN, GAMZ_IN ARE ALL FIXED 
     >TO BE THE SAME AS FEWZ3.1 INPUT'
           XMW=80.403d0
           XMZ=91.1876d0
           GAMW_IN=2.141d0
      	    GAMZ_IN=2.4952d0
         ENDIF

CJI Nov 2020
        if(KFAC_ANGFUNC.gt.2 .or. KFAC_ANGFUNC.lt.0) then
            print*, ' Invalid angular k-factor.'
            print*, ' Allowed values are 0, 1, and 2'
            call exit
        endif

CCPY!*********************************************************************************
CCPY!
CCPY!      IF(LTOpt.EQ.9) THEN
CCPY!CCPY JAN 2016; INPUT VALUES OF XMW AND XMZ FOR PROCESSES WITH BUILD-IN AMPLITUDES
CCPY!CsB_____This is the case when we use the built in matrix elements.
CCPY!CsB_____Since the EW parameters are not read from the grid:
CCPY!CZL        XMW = 80.36
CCPY!CZL        XMZ = 91.187
CCPY!        XMW = 80.419
CCPY!        XMZ = 91.188
CCPY!C        Print*, ' mt, mH = ', XMTop, XMHiggs
CCPY!        Print*, ' Hard wired MW, MZ are used. '
CCPY!cpn2005 Npart for built-in matrix elements is set in resbos.for
CCPY!      ENDIF
CCPY!
CCPY!*********************************************************************************

      Close(5)
CsB___End of reading input file

CCPY SEPT 2009: I TURN OFF THIS OPTION.
      If (Boson.Eq.'ZG') then
        PRINT*,' gg -> ZZ has not been implemented in this version '
        Call ExitGracefully
      ENDIF

CsB___Mass of primary final state particles
      If (Boson.Eq.'GL') then
        XMASS(1)=XMTOP
        XMASS(2)=XMTOP
      Else If (Boson.Eq.'ZZ' .or. Boson.Eq.'ZG') then
        XMASS(1) = xMZ
        XMASS(2) = xMZ
      Else If (Boson.Eq.'H0') then
CCPY MARCH 2008: TO INCLDUE THE CASE FOR OFF-SHELL W AND Z, WE SET 
C XMASS(1)=XMASS(2)=0 FOR 'ZZ' AND 'WW' DECAY MODES OF 'H0'.
        If (hDecayMode.Eq.'ZZ') then
          XMASS(1) = 0.D0
          XMASS(2) = 0.D0
        ElseIf (hdecayMode.Eq.'WW')then
          XMASS(1) = 0.D0
          XMASS(2) = 0.D0
        ElseIf (hDecayMode.Eq.'AA') then
          XMASS(1) = 0.d0
          XMASS(2) = 0.d0
        End If
      Else If (Boson.Eq.'HP' .or. Boson.Eq.'HM') then
        XMASS(1) = xMW
        XMASS(2) = xMHiggs
      Else If (Boson.Eq.'HZ') then
        XMASS(1) = xMZ
        XMASS(2) = xMHiggs
      ELSE
        XMASS(1)=0.0D0
        XMASS(2)=0.0D0
      ENDIF

CsB___Determine invariant mass (Q) range
      IF(YM12MIN.LT.X3A(0)) THEN
        Print*, ' !! Warning: Adjusting minimal Q to', X3A(0), ' !!'
        YM12MIN=X3A(0)
      ENDIF
      If (Boson.Eq.'ZZ' .or. Boson.Eq.'ZG') ! No Z Z production if Q < 2 Mz
     >    YM12MIN = Max(2.d0*xMz+.1, YM12MIN)
      If (Boson.Eq.'HP' .or. Boson.Eq.'HM') ! No W h production if Q < MW+Mh
     >    YM12MIN = Max(xMW+xMHiggs+.1, YM12MIN)
      If (Boson.Eq.'HZ') ! No Z h production if Q < MZ+Mh
     >    YM12MIN = Max(xMZ+xMHiggs+.1, YM12MIN)
      WRITE(*,*) ' Minimal invariant mass (Q):       ', YM12MIN
      IF(YM12MAX.GT.X3A(iD3) .and. LTOpt.NE.9 ) THEN
        YM12MAX = X3A(iD3)
        Print*, ' !! Warning: Adjusting maximal Q to', X3A(iD3), ' !!'
      ENDIF
      WRITE(*,*) ' Maximal invariant mass (Q):       ', YM12MAX
      IF(YM12MIN.GT.YM12MAX) THEN
        WRITE(*,*) ' !! Warning: Fatal input YM12MIN > YM12MAX !!'
        CALL ExitGracefully
      ENDIF

CsB___Determine the transverse momentum (QT) range
      If (QTMax.GT.X1A(iD1) .and. LTOpt.NE.9) then
        QTMax    = X1A(iD1)
C        QT_V_Max = X1A(iD1)
        Print*, ' !! Warning: Adjusting maximal QT to', X1A(iD1), ' !!'
      End If
      If (nPart.Eq.11)
     .  WRITE(*,*) ' Minimum transverse momentum (QT): ', QTMin
      WRITE(*,*) ' Maximum transverse momentum (QT): ', QTMax
      If (iYPiece.Eq.1) then
        IF(YM12MAX.GT.X3Y(iY3) .and. LTOpt.NE.9 ) THEN
          WRITE(*,*) ' !! Warning: YM12MAX IS RESET TO BE: ', X3Y(iY3)
          YM12MAX = X3Y(iY3)
        ENDIF
      End If
      ETOTMN=YM12MIN

CsB___For processes with built-in amplitudes:
      If (LTOpt.EQ.9) then
        Q_V_Max  = YM12MAX
        QT_V_MIN = QT_V_SEP
        QT_V_MAX = QTMax
        QT_V_MXY = QTMax
        YB_V_MAX =   5.0
      End If

CsB___If grid is used then scale = invariant mass (Q)
      If (NPart.LE.19) ITSCALE = 4

CsB___Set PDF for processes with built-in amplitudes:
      If (LTOpt.Eq.9 .Or.
     .   (nPart.Eq.11.And.Boson.Eq.'AA') .Or.
     .   (nPart.Eq.13.And.Boson.Eq.'AA') .Or.
     .   (nPart.Eq.13.And.Boson.Eq.'AG')) then
        If (BIPDF.GT.40 .and. BIPDF.LT.50) then
          ISet = BIPDF - 40   ! CTEQ4M
          Print*, ' Initializing CTEQ4 '
          Call SetCtq4(ISet)
        Else If (BIPDF.GT.50 .and. BIPDF.LT.60) then
          ISet = BIPDF - 50   ! CTEQ5M
          Print*, ' Initializing CTEQ5 '
          Call SetCtq5(Iset)
CCPY        Else If (BIPDF.GT.60 .and. BIPDF.LT.70) then
CCPY          ISet = BIPDF - 60   ! CTEQ6M
        Else If (BIPDF.GE.200 .AND. BIPDF.LT.500) then ! CTEQ6
          ISet = BIPDF   ! CTEQ6
          Print*, ' Initializing CTEQ6 '
          Call SetCtq6(Iset)
CCPY
        Else If (BIPDF.EQ.903) then ! USING pdf.pds file
          ISet = BIPDF   ! CTEQ6
          Print*, ' The file pdf.pds should be in the same directory '
          Print*, ' Initializing CTEQ6 '
          Call SetCtq6(Iset)
        End If
      End If

CsB   nYEdge = 0 do nothing about the edge of the y phase space
C              when interpolating.
      nYEdge = 0

CsB___This call was moved into SUBROUTINE GENERATOR
C      CALL STANDARD

CsB___As of 7/1/97 matching is done at the gird level.

CsB___Initialize histograms
      CALL HistInit

CsB___Main calculation part
      Do ii = 1,iGen
        CALL GENERATOR
      End Do

CsB___Write header of output file
      WRITE(16,*) ' W = ',W, '   TYPE_V = ', Boson, '   iProc =', iProc
      write(16,*) ' XMASS(1), XMASS(2) =',XMASS(1),XMASS(2)
      WRITE(16,*) ' PTMIN, ETAMAX =',PTMIN,ETAMAX
      WRITE(16,*) ' PTECUT,YECUT,PTNCUT,yNCut =',
     >PTECUT,YECUT,PTNCUT,yNCut
      WRITE(16,*) ' YLUMI =',YLUMI, '    Output = ',Output

CsB___Output Histograms
      CALL HistOutp

      RETURN
      END ! YUAN_MAIN

C -------------------------------------------------------------------------
      Subroutine HistOutp
C -------------------------------------------------------------------------
	Implicit None
	Character*8 Output
	Common / Output / Output

	If (Output.Eq.'ROOTNT1'.or.Output.Eq.'ROOTNT2') then
          Print *,'Closing the ROOT ntuple'
          call ROOTNTOutp
	ElseIf (Output.Eq.'StdHEP') then
          Print*, ' Closing StdHEP...'
          call StdHEPOutp
	ElseIf (Output.Eq.'Text1'.or.Output.Eq.'Text2') then
          Print*, ' Closing ASCII files with events'
          close(18)
	ElseIf (Output.Eq.'GBook ') then
          Print*, ' Closing GBook...'
	  CALL GCLEAR
	Else
	  Print*, ' Problem in HistOutP. Output: ', Output
	  Stop
	EndIf

	Return
	End ! HistOutP

CsB * *********************************************************************

      subroutine StdHEPOutp
cpn2008 Close StdHEP output file
      implicit none
      integer photoscor, event_num, decaymode, wpm
      real*4 mymaxwgt
      common/stdhep/mymaxwgt, photoscor, event_num, decaymode, wpm

      write (18,*) 0,0
      close (18) 

      PRINT*, 'Maximum event weight = ', mymaxwGT

      return 
      end                       !StdHEPOutp


CsB * *********************************************************************

      SUBROUTINE HISTINIT
	Implicit Real*4 (A-H,O-Z)
	Character*8 Output
	Common / Output / Output

        INTEGER*4 NBINIT, iwgt
        REAL*8 MAX_WGT,NEW_MAX_WGT,WGT_FACTOR
        COMMON/UNWEIGHT/MAX_WGT,NEW_MAX_WGT,WGT_FACTOR,iwgt,NBINIT

c htj 
      character*100 jobname
      integer jlen

      if(iargc().eq.0)then
        jobname='resbos'
      else if(iargc().ge.1) then
        call getarg(1,jobname)
      endif
      jlen=len_trim(jobname)

	If (Output.Eq.'ROOTNT1'.or.Output.Eq.'ROOTNT2') then
cpn2008   Initialize a ROOT ntuple
          if (iwgt.eq.-2) then
            Call InitROOTNT(jobname(1:jlen)//".root",'UPDATE')
c            Call InitROOTNT('resbos.root','UPDATE')
          else 
            Call InitROOTNT(jobname(1:jlen)//".root",'RECREATE')
c            Call InitROOTNT('resbos.root','RECREATE')
          endif                 !iwgt

          call set_element_tags !must be called after the NT initialization 
          call SetRootNTBranches
	ElseIf (Output.Eq.'NTuple') then
cpn20008
          print *,
     >      'STOP: the NTuple output type in resbos.in is phased out'
          print *,'Use ROOTNT1 or PAWNT1 to record kinematics in terms',
     >      'of Q, y, Q_T, cos(theta), phi, etc.; or ROOTNT2 or PAWNT2',
     >      'to record kinematics in terms of 4-momenta of secondary',
     >      'particles'
          call Quit
	ElseIf (Output.Eq.'StdHEP') then
cpn2008   Initialize the event output in StdHEP format
          call set_element_tags
          call InitStdHEP
	ElseIf (Output.Eq.'Text1'.or.Output.Eq.'Text2') then
cpn2008   Initialize the event output in ASCII format
          call set_element_tags
          call InitText(Output)
	ElseIf (Output.Eq.'GBook ') then
CsB   Initialize GBook
C         Print*, ' Initializing GBook...'
	  Call InitGB
	Else
	  Print*, 'HistInit: wrong output type ', Output
	  Stop
	EndIf

C       Print*, ' Histograms are initialized.'

      RETURN
      END ! HistInit


cpn2008 ********************************************************************
      subroutine SetRootNTBranches
      implicit none
      integer max_siz,i
      parameter(max_siz=50)
      Character*80 chtitl
      Character*8  chtags(max_siz)
      common /NT_TAGS/chtitl, chtags

      real*4 event_vector(max_siz)
      common/evec/event_vector

      Integer lunit, nt_id, icycle, nt_elements
      Character*80 nt_file_name
      Character*8  nt_subdir
      Common / NTPL / nt_id, lunit, nt_subdir, nt_file_name, nt_elements


      do i=1, nt_elements
        print *,'chtags(',i,') =',chtags(i),'=========='
        call AddNTBranch( event_vector(i), chtags(i)) 
      enddo                     !i

      return 
      end !SetRootNTBranches

CsB * *********************************************************************

      subroutine InitText(Output)
cpn2008 Initialize the output of the events in the text format
      implicit none
      character*5 Output

      INTEGER*4 NBINIT, iwgt
      real*8 MAX_WGT,NEW_MAX_WGT,WGT_FACTOR
      COMMON/UNWEIGHT/MAX_WGT,NEW_MAX_WGT,WGT_FACTOR,iwgt,NBINIT

      open(unit=18,file='resbos_events.dat',form='formatted')

      if (iwgt.eq.-2) then
        write (18,*) '#Weight'
      elseif (Output.Eq.'Text1') then
       write (18,*)
     >    '# px_ele,  py_ele, pz_ele, E_ele, px_nu, py_nu, pz_nu, E_nu,
     >    px_V, py_V, pz_V, E_V, Weight'
      else if (Output.Eq.'Text2') then 
        write (18,*) 
CCPY     >    '#pt_ele, pt_nu, eta_ele, eta_nu, pt_V, y_V, M_V, cos(theta_*)
C     >    phi_*, Del_phi34, DelR34, Weight, 0' 
     >    '#pt_ele, pt_nu, eta_ele, eta_nu, pt_V, y_V, M_V, phi_ele, 
     >    cos(theta_*), phi_nu, DelR34, Weight, 0' 
      endif                     !Output
      
      return 
      end

cpn2008 =================================================================

      Subroutine InitGB
      IMPLICIT REAL*4 (A-H,O-Z)
      REAL*4 A
      COMMON /GBOOK/ A(10000)
      Character*40 Boson
      Common / Boson / Boson
      Real*8 W,QCDL,X1X2MIN
      COMMON/AA/W,QCDL,X1X2MIN,NJET,N_PRIMARY_PART,NF,I_STR,IBEAM,NPART,
     &  NEVENT,NORDER,N_SECOND_PART

C-----------------------------------------------------------------------
CsB Book Gbook
C-
      If (Boson.Eq.'Z0' .Or. Boson.Eq.'W+' .Or.
     .    Boson.Eq.'W-'.Or.Boson.Eq.'A0'.Or.Boson.Eq.'DY') then
        CALL GBOOK1(1,'p_T^l+',50,0.0,60.0)
        CALL GBOOK1(2,'p_T^n',50,0.0,60.0)
C        CALL GBOOK1(3,'y^l+',40,-3.0,3.0)
C        CALL GBOOK1(4,'y^n',50,-5.0,5.0)
        CALL GBOOK1(3,'y^l+',50,0.0,5.0)
        CALL GBOOK1(4,'y^n',50,0.0,5.0)
C        CALL GBOOK1(5,'M^V',50,50.0,120.0)
        CALL GBOOK1(5,'M^V',50,100.0,1100.0)
C        CALL GBOOK1(6,'p_T^V',50,0.0,120.0)
        CALL GBOOK1(6,'p_T^V',50,0.0,25.0)
        CALL GBOOK1(7,'y^V',50,-3.0,3.0)
        CALL GBOOK1(8,'M_T^V',30,2800.0,3100.0)
        CALL GBOOK1(9,'CTS_exp',50,-1.0,1.0)
        CALL GBOOK1(10,'Del phi',50,0.0,3.1416)
        CALL GBOOK1(11,'E^l+' ,50,0.0,180.0)
        CALL GBOOK1(12,'E^n',50,0.0,180.0)
        CALL GBOOK1(13,'p_Z^V',50,-700.0,700.0)
        CALL GBOOK1(14,'CTS_PHASE',50,-1.0,1.0)
        CALL GBOOK1(15,'phi_sta',50,0.0,3.1416)
CCPY        CALL GBOOK1(16,'Z_Corr',50,-1.0,1.0)
        CALL GBOOK1(17,'Del_R',50,0.0,6.2)
        CALL GBOOK1(18,'p_T^V',50,0.0,30.0)
        CALL GBOOK1(19,'M_T^V',50,70.0,85.0)
        CALL GBOOK1(20,'p_T^l+',50,35.0,45.0)
        CALL GBOOK1(21,'p_T^n',50,35.0,45.0)
CCPY
        CALL GBOOK1(22,'y^l+',50,-5.0,5.0)
        CALL GBOOK1(23,'y^n',50,-5.0,5.0)
        CALL GBOOK1(24,'y^V',50,-5.0,5.0)
        
CCPY      ElseIf (Boson.Eq.'H0') then
      ElseIf (Boson.Eq.'H0'
     >.OR.Boson.Eq.'AA'.OR.Boson.Eq.'AG') then
        CALL GBOOK1(1,'p_T^1',50, 0.0,150.)
        CALL GBOOK1(2,'p_T^2',50, 0.0,150.)
        CALL GBOOK1(3,  'y^1',50,-5.0, 5.0)
        CALL GBOOK1(4,  'y^2',50,-5.0, 5.0)
        CALL GBOOK1(5,  'M^V',50,225.,275.)
        CALL GBOOK1(6,'p_T^V',50, 0.0,50.0)
        CALL GBOOK1(7,  'y^V',50,-5.0, 5.0)
        CALL GBOOK1(8,'M_T^V',50,225.,275.)
        CALL GBOOK1(9,'M_T^V',50,225.,275.)
        CALL GBOOK1(10,'Del phi',50,0.0,3.1416)
        CALL GBOOK1(11,'E^1' ,50,150.0,950.0)
        CALL GBOOK1(12,'E^2',50,150.0,950.0)
        CALL GBOOK1(13,'p_Z^V',50,-700.0,700.0)
        CALL GBOOK1(14,'Cos(the_sta)',50,-1.0,1.0)
        CALL GBOOK1(15,'phi_sta',50,0.0,3.1416)
c        CALL GBOOK1(16,'Z_Corr',50,-1.0,1.0)
        CALL GBOOK1(17,'Del_R',50,0.0,6.2)
        CALL GBOOK1(18,'p_T^V',50,0.0,30.0)
        CALL GBOOK1(19,'M_T^V',50,200.0,300.0)
        CALL GBOOK1(20,'p_T^1',50,0.0,100.0)
        CALL GBOOK1(21,'p_T^2',50,0.0,100.0)
CCPY
        CALL GBOOK1(22,'y^1',50,-5.0,5.0)
        CALL GBOOK1(23,'y^2',50,-5.0,5.0)
        CALL GBOOK1(24,'y^V',50,-5.0,5.0)
 
      Else
        If (iBeam.NE.-2) then
          CALL GBOOK1(1,   'p_T^1', 50, 0.0, 12.0)
          CALL GBOOK1(2,   'p_T^2', 50, 0.0, 12.0)
          CALL GBOOK1(3,     'y^1', 50,-1.3,  1.3)
          CALL GBOOK1(4,     'y^2', 50,-1.3,  1.3)
          CALL GBOOK1(5,    'M^12', 50, 0.0, 75.0)
          CALL GBOOK1(6,  'p_T^12', 50, 0.0, 40.0)
          CALL GBOOK1(7,    'y^12', 50,-1.0,  1.0)
          CALL GBOOK1(8,   'M?T!^[gg]$', 50,79.0, 81.5)
          CALL GBOOK1(9,   'M?T!^[gg]$', 50, 0.0, 90.0)
          CALL GBOOK1(10,  '[Df]^[gg]$', 50, 0.0, 3.15)
          CALL GBOOK1(11,    'E?^[g]1$', 50, 0.0, 60.0)
          CALL GBOOK1(12,    'E?^[g]2$', 50, 0.0, 60.0)
          CALL GBOOK1(13,  'p?Z!^[gg]$', 50,-700.0,700.0)
          CALL GBOOK1(14,   'CosThe^*$', 50,-1.0,  1.0)
          CALL GBOOK1(15,      'phi^*$', 50, 0.0, 6.28)
          CALL GBOOK1(16,     'Z?Corr$', 50, -1.0, 1.0)
          CALL GBOOK1(17,     'Del_R34',50,0.0,6.2)
          CALL GBOOK1(18,     'Del_R35',50,0.0,6.2)
          CALL GBOOK1(19,     'Del_R45',50,0.0,6.2)
          CALL GBOOK1(20,         'M35',50,0.0,80.)
          CALL GBOOK1(21,         'M45',50,0.0,80.)
        Else
          CALL GBOOK1(1,  'p_T^l+',50, 0.0, 10.0)
          CALL GBOOK1(2,   'p_T^l',50, 0.0, 10.0)
          CALL GBOOK1(3,    'y^l+',50,-1.25,1.25)
          CALL GBOOK1(4,     'y^l',50,-1.25,1.25)
          CALL GBOOK1(5,     'M^V',50, 0.0,20.0)
          CALL GBOOK1(6,   'p_T^V',50, 0.0,40.0)
          CALL GBOOK1(7,     'y^V',50,-1.0, 1.0)
          CALL GBOOK1(8,   'M_T^V',50,79.0,81.5)
          CALL GBOOK1(9,   'M_T^V',50, 0.0,120.0)
          CALL GBOOK1(10,'Del phi',50, 0.0,3.1416)
          CALL GBOOK1(11,   'E^l+',50, 0.0,50.0)
          CALL GBOOK1(12,    'E^l',50, 0.0,50.0)
          CALL GBOOK1(13,  'p_Z^V',50,-700.0,700.0)
          CALL GBOOK1(14,'Cos(the_sta)',50,-1.0,1.0)
          CALL GBOOK1(15,'phi_sta',50,0.0,3.1416)
c          CALL GBOOK1(16, 'Z_Corr',50,-1.0,1.0)
          CALL GBOOK1(17,  'Del_R34',50,0.0,6.2)
          CALL GBOOK1(18,  'Del_R35',50,0.0,6.2)
          CALL GBOOK1(19,  'Del_R45',50,0.0,6.2)
        End If !ibeam
      End If

      If (NPart.Eq.21.or.npart.eq.22) then
          CALL GBOOK1(88,  'Del_R35',50,0.0,6.2)
          CALL GBOOK1(89,  'Del_R45',50,0.0,6.2)
          CALL GBOOK1(90,  'M_35',50,0.0,200.)
          CALL GBOOK1(91,  'M_45',50,0.0,200.)
        End If                  !Npart

C-
C-----------------------------------------------------------------------
	Return
	End ! InitGB


      subroutine  InitStdHEP
cpn2008 Initialize StdHEP (text) output
      implicit none
      integer photoscor, event_num, decaymode, wpm
      real*4 mymaxwgt
      common/stdhep/mymaxwgt, photoscor, event_num, decaymode, wpm
      
      open (unit=18,file='resbos.hep',FORM='FORMATTED')
      mymaxwgt=0d0; event_num=0; decaymode=1; wpm=0

      photoscor=0
CZHU  IF WE WANT TO USE PHOTOS FOR FSR, NEED TO INITIALIZE PHOTOS FIRST
c      IF(PHOTOSCOR.EQ.1) THEN
c        CALL PHOINI
c      ENDIF
      
      return 
      end !InitStdHEP

C -------------------------------------------------------------------------
      SUBROUTINE PRIMARY_CUTS(PASS_CUTS)
C -------------------------------------------------------------------------
      IMPLICIT NONE
      LOGICAL PASS_CUTS
      REAL*8 PBOOST(4)
      REAL*8 PLAB3(4),PLAB4(4),PCMS3(4),PCMS4(4)

      REAL*8 W,QCDL,X1X2MIN
      INTEGER*4 NJET,N_PRIMARY_PART,NF,I_STR,IBEAM,NPART,NEVENT,NORDER,
     ,          N_SECOND_PART
      COMMON/AA/W,QCDL,X1X2MIN,NJET,N_PRIMARY_PART,NF,I_STR,IBEAM,NPART,
     &  NEVENT,NORDER,N_SECOND_PART
      REAL*8 SEC_PLAB,PLAB,PCEM,WHAT
      COMMON/DECAY_MOMENTA/SEC_PLAB(4,50)
      COMMON/MOM/PLAB(4,12),PCEM(4,12),WHAT
      REAL*8 ETOTMN,PTMIN,ETAMAX,ETASEP
      INTEGER*4 IAPPROX
      COMMON/CUTVAL/ETOTMN,PTMIN,ETAMAX,ETASEP,IAPPROX

      REAL*8 PTECUT,YECUT,DELR_CUT,PTNCUT,yNCut,YWLCUT,YWUCUT
      COMMON/MYCUT/PTECUT,YECUT,DELR_CUT,PTNCUT,yNCut,YWLCUT,YWUCUT
      REAL*8 YM12,YPT12,YP12,ymij, ymij2, tranm
      COMMON/MY1/YM12,YPT12,YP12(4)
      REAL*8 YM12MIN,YM12MAX,QTMin,QTMax, TMmin, TMMax, Etmin
      COMMON/MYYM12/YM12MIN,YM12MAX,QTMin,QTMax, TMMin, TMMax, Etmin
      REAL*8 YXP,YXETA,YXPT,YXDELR,YXPP,yxy(30)
      COMMON/SECCONV/ YXP(4,30),YXETA(30),YXPT(30),YXDELR(30,30),
     >YXPP(30)
      Character*40 Boson
      Common / Boson / Boson
      Integer iYPiece, iReadIn
      Common / YPiece / iYPiece, iReadIn

      Character*2 hDecayMode
      Common / hDM / hDecayMode

      Real*8 YXPT3,YXPT4,PTXCUT

      REAL*8 QT_V_MIN, QT_V_MAX,YB_V_MAX,Q_V_MIN,Q_V_MAX,QT_V_MXY
      COMMON/HARDCUT/ QT_V_MIN, QT_V_MAX,YB_V_MAX,
     >  Q_V_MIN,Q_V_MAX,QT_V_MXY

      REAL*8 Q_V,QT_V,YB_V,THE_STA,PHI_STA
      COMMON/CSFRAME/Q_V,QT_V,YB_V,THE_STA,PHI_STA

      Integer LTOpt
      Common / Process / LTOpt

      INTEGER I,J,NFINAL,ID_ELE,ID_NU
      REAL*8 YXDELETA,YXCOS1,YXCOS2,YXSIN1,YXSIN2,YXCOSDPH,YXDELPHI
CsB   LOGICAL NOJET - Removed on 05/10 by Pavel
      REAL*8 ZCORR,ZCORR_CUT

      Integer ACut
      Common / AAGCut / ACut
      Real*8 CosTheta
      REAL*8 XMW,XMZ,XMW2,XMZ2,AEM,GAMW,GAMZ,GAMW2,GAMZ2
      COMMON/STAND0/XMW,XMZ,XMW2,XMZ2,AEM,GAMW,GAMZ,GAMW2,GAMZ2
 
cpn2005 Isolation cone parameters
      integer iqqt
      real*8 DelRIso, ETmaxIso
      common/isolation/DelRIso, ETMaxIso, iqqt
      data iqqt /0/

      integer isub
      double precision cutM45, SubtractionShape, cutSIJ, rmij2, Aexp
      logical SubQGFlag
      Common / Subtraction / cutM45, SubtractionShape(3:4),iSub,
     >  SubQGFlag(3:4)

cpn2005 For smooth cone isolation
      real*8 ConeEtmax
      external ConeEtmax

      Logical Debug
      Data Debug /.False./

       real*8 ycut_1,ycut_2,ycut_3

cdump
C      Debug = .True.

      If (Debug) Print*, ' Primary cuts -------------------------- '

      PASS_CUTS=.false.

CsB___Rename lab-momenta
      DO 517 J=1,12
       DO 517 I=1,4
        YXP(I,J)=PLAB(I,J)
517   CONTINUE

CsB___Construct pT's, rapidities, and pseudorapidities 
C     for the produced particles
      DO 519 J=3,12
        YXPT(J)=DSQRT(YXP(1,J)**2+YXP(2,J)**2)
        IF((YXP(4,J)-YXP(3,J)).GT.0.0 .AND.
     >  (YXP(4,J)+YXP(3,J)).GT.0.0)THEN
          YXY(J)=0.5*DLOG((YXP(4,J)+YXP(3,J))
     >    /(YXP(4,J)-YXP(3,J)))
          costheta=yxp(3,j)/sqrt(yxp(1,j)**2+yxp(2,j)**2+yxp(3,j)**2)
          yxeta(j)=0.5*dlog((1d0+costheta)/(1d0-costheta))
        ELSEIF((YXP(4,J)-YXP(3,J)).LE.0.)THEN
C          Print*, ' (YXP(4,J)-YXP(3,J)).LE.0. '
          YXY(J)=20.
          YXETA(J)=20.
        ELSEIF((YXP(4,J)+YXP(3,J)).LE.0.)THEN
C          Print*, ' (YXP(4,J)+YXP(3,J)).LE.0. '
          YXY(J)=-20.
          YXETA(J)=-20.
        END IF
519   CONTINUE

CCPY This is for LHCb Z production only. 

CCPY ID_NU refers to fermion and ID_ELE to anti-fermion.
C For example, in Z -> e- e+, ID_NU refers to e-.
C In W^+ -> \nu e^+ , ID_ELE refers to e^+.
C In W^- -> e^- {\bar\nu}  , due to the variable "sign=-1" included in 
C resbosm.f , ID_ELE refers to e^-.
C
      ID_NU=3
      ID_ELE=4

CsB___Construct pT3&4
      DO 513 I=1,4
        YP12(I)=0.D0
        DO 515 J=3,4
          YP12(I)=YP12(I)+PLAB(I,J)
515     CONTINUE
513   CONTINUE

      If (Debug) Print*,' PLAB(I,3) = ',(PLAB(I,3),I=1,4) 
      If (Debug) Print*,' PLAB(I,4) = ',(PLAB(I,4),I=1,4) 

CsB Construct pT3&4 (magnitude), M3&4 and MT3&4
      YPT12 = DSQRT(YP12(1)**2+YP12(2)**2)
      YM12  = DSQRT(YP12(4)**2-YP12(3)**2-YPT12**2)
C      YMT12 = DSQRT(YPT12 + YM12)
C      Print*, ' YPT12,YM12 = ', YPT12,YM12

CsB___VB level cuts read from the input file
C   --- INVARIANT-MASS3&4 range ---
      If (Debug) Print*, YM12MIN, '< M12=', YM12, '<', YM12MAX
      IF(YM12.GT.YM12MAX .OR. YM12.LT.YM12MIN) RETURN
C   --- QT 3&4 Maximum ---
      If (Debug) Print*, ' QT34 = ', YpT12, ' < ', QTMax
      IF(YpT12.GT.QTMax) RETURN

CsB___Among others: to study the high pT perturbative contribution
C   to the low pT individual photon distribution
C   in AA & AG case.
C   Cut on
C   --- QT 3&4 Minimunm ---
      If (Debug) Print*, ' QT34 = ', YpT12, ' > ', QTMin
      If (YPT12.LT.QTMin) Return

C___HARD-WIRED CUTS ON VECTOR-BOSON:
CsB QT_V_MXY is read from the Y grid.
C   The last point of the Y piece grid in Q_T is read from the grid.
C     --- PT12 maximum ---
      If (Debug) Print*, ' YXPT(9) =', YXPT(9),'<', QT_V_MXY
      IF(YXPT(9).GT.QT_V_MXY)  RETURN
C     --- Y12 maximum ---
      If (Debug) Print*, ' ABS(YXY(9)) < ', YB_V_MAX
      IF(ABS(YXY(9)).GT.YB_V_MAX)  RETURN
C     --- INVARIANT-MASS12 range ---
      If (Debug) Print*, Q_V_MIN, ' < YM12 < ', Q_V_MAX
      IF(YM12.LT.Q_V_MIN .OR. YM12.GT.Q_V_MAX)  RETURN

CCPY May 2015
CPN   Cuts on the rapidity of the vector boson 
      if (yxy(9).lt.YWLCUT.or.yxy(9).gt.YWUCUT) return

        YXPT4 = YXPT(ID_ELE)
        YXPT3 = YXPT(ID_NU)
 
C_____CUTS:


CCPY ID_NU refers to fermion and ID_ELE to anti-fermion.
C For example, in Z -> e- e+, ID_NU refers to e-.

CCPY Feb 2013
        If (Boson. Eq. 'W+' .Or. Boson. Eq. 'W-') then

CCPY ID_NU refers to fermion and ID_ELE to anti-fermion.
C For example, in Z -> e- e+, ID_NU refers to e-.
C In W^+ -> \nu e^+ , ID_ELE refers to e^+.
C In W^- -> e^- {\bar\nu}  , due to the variable "sign=-1" included in 
C resbosm.f , ID_ELE refers to e^-.
C
CPN 31.01.2001 PTNCUT and YNCUT are rewired to serve as upper cuts 
CPN on the electron's momentum and rapidity

C For Cuts(1): PTECUT,YECUT,DELR_CUT,PTNCUT,YnCut
C For Cuts(1): pTlmin,ymin,DelR,pTlmax,ymax

          IF(YXPT(ID_ELE) .LT. PTECUT.OR.YXPT(ID_ELE) .GT. PTNCUT)
     >      RETURN

CCPY May 2015
C         IF(YXETA(ID_ELE) .LT. YECUT.OR.YXETA(ID_ELE) .GT. YNCUT)
C     >      RETURN

          IF(ABS(YXETA(ID_ELE)) .LT. YECUT .OR.
     > ABS(YXETA(ID_ELE)).GT. YNCUT) RETURN

          IF(YXPT(ID_NU) .LT. ETMin) return



CBIN inclusive cross section
c         YCUT_1=1.4442
c         YCUT_2=1.566
c         YCUT_3=2.5

c         IF(ABS(YXETA(ID_ELE)).GT.YCUT_1.AND.
c     >        ABS(YXETA(ID_ELE)).LT.YCUT_2 )   RETURN

c         IF(ABS(YXETA(ID_ELE)).GT.YCUT_3)  RETURN





        Else If (Boson. Eq. 'Z0' .Or. Boson. Eq. 'A0' .Or. 
     > Boson. Eq. 'DY') then

CPN 31.01.2001 PTNCUT and YNCUT are rewired to serve as upper cuts 
CPN on the electron's momentum and rapidity

          IF(YXPT(ID_ELE) .LT. PTECUT.OR.YXPT(ID_ELE) .GT. PTNCUT)
     >      RETURN
          IF(YXPT(ID_NU) .LT. PTECUT.OR.YXPT(ID_NU) .GT. PTNCUT)
     >      RETURN

CCPY THIS IS FOR GENERAL PURPOSE CUTS:
          IF(YXETA(ID_ELE) .LT. YECUT.OR.YXETA(ID_ELE) .GT. YNCUT)
     >      RETURN
          IF(YXETA(ID_NU) .LT. YECUT.OR.YXETA(ID_NU) .GT. YNCUT)
     >      RETURN

CCPY IMPOSE HARD-WIRED CUTS HERE FOR COMPARING TO CMS DATA:


CBIN inclusive cross section
c         YCUT_1=1.4442
c         YCUT_2=1.566
c         YCUT_3=2.5

c         IF(ABS(YXETA(ID_ELE)).LE.YCUT_1) THEN 
c            IF(ABS(YXETA(ID_NU)).GT.YCUT_3) RETURN
c            IF(ABS(YXETA(ID_NU)).LT.YCUT_2.AND.
c     >         ABS(YXETA(ID_NU)).GT.YCUT_1) RETURN

c         ELSEIF(ABS(YXETA(ID_NU)).LE.YCUT_1) THEN 
c            IF(ABS(YXETA(ID_ELE)).GT.YCUT_3) RETURN
c            IF(ABS(YXETA(ID_ELE)).LT.YCUT_2.AND.
c     >         ABS(YXETA(ID_ELE)).GT.YCUT_1) RETURN
c         ELSE
c            RETURN
c         ENDIF


        else 
          PRINT*, ' THIS BOSON IS NOT YET IMPLEMENTED '
          call ExitGracefully

        end if  !Boson
      

CsB___Construct DeltaPhi and DeltaR
        YXDELETA=YXETA(3)-YXETA(4)
        YXCOS1=YXP(1,3)/YXPT(3)
        YXSIN1=YXP(2,3)/YXPT(3)
        YXCOS2=YXP(1,4)/YXPT(4)
        YXSIN2=YXP(2,4)/YXPT(4)
        YXCOSDPH=YXCOS1*YXCOS2+YXSIN1*YXSIN2
        If(ABS(YXCOSDPH).LE.1.0D0) THEN
          YXDELPHI=ACOS(YXCOSDPH)
        ENDIF

        YXDELR(3,4)=DSQRT(YXDELETA**2+YXDELPHI**2)
        YXDELR(4,3)=YXDELR(3,4)

C    --- DELTA R minimum ---
          If (Debug) Print*, ' YXDELR(3,4) > ', DELR_CUT
          IF(YXDELR(3,4) .LT. DELR_CUT) RETURN


CPN Cut on the "leptonic" transverse mass, which is defined in J. Smith,
C   W. L. van Neerven, and J.A.M. Vermaseren, PRL 50, 1738 (1983) and 
C   is not the same as sqrt(Q^2 + Q_T^2)

          TRANM=2.0*YXPT(ID_ELE)*YXPT(ID_NU)*(1.0D0-COS(YXDELPHI))
          TRANM=SQRT(TRANM)

CPN   Cuts on the transverse mass of the vector boson
          if (TRANM.lt.TMMin .or. TRANM.gt.TMMax) return

C===============

      PASS_CUTS=.TRUE.
      If (Debug) Print*, ' Passed primary-cuts ... '
      
      RETURN
      END                       ! PRIMARY_CUTS 



C -------------------------------------------------------------------------
      SUBROUTINE PRIMARY_CUTS_ALL(PASS_CUTS)
C -------------------------------------------------------------------------
      IMPLICIT NONE
      LOGICAL PASS_CUTS

      REAL*8 W,QCDL,X1X2MIN
      INTEGER*4 NJET,N_PRIMARY_PART,NF,I_STR,IBEAM,NPART,NEVENT,NORDER,
     ,          N_SECOND_PART
      COMMON/AA/W,QCDL,X1X2MIN,NJET,N_PRIMARY_PART,NF,I_STR,IBEAM,NPART,
     &  NEVENT,NORDER,N_SECOND_PART
      REAL*8 SEC_PLAB,PLAB,PCEM,WHAT
      COMMON/DECAY_MOMENTA/SEC_PLAB(4,50)
      COMMON/MOM/PLAB(4,12),PCEM(4,12),WHAT
      REAL*8 ETOTMN,PTMIN,ETAMAX,ETASEP
      INTEGER*4 IAPPROX
      COMMON/CUTVAL/ETOTMN,PTMIN,ETAMAX,ETASEP,IAPPROX

      REAL*8 PTECUT,YECUT,DELR_CUT,PTNCUT,yNCut,YWLCUT,YWUCUT
      COMMON/MYCUT/PTECUT,YECUT,DELR_CUT,PTNCUT,yNCut,YWLCUT,YWUCUT
      REAL*8 YM12,YPT12,YP12,ymij, ymij2, tranm
      COMMON/MY1/YM12,YPT12,YP12(4)
      REAL*8 YM12MIN,YM12MAX,QTMin,QTMax, TMmin, TMMax, Etmin
      COMMON/MYYM12/YM12MIN,YM12MAX,QTMin,QTMax, TMMin, TMMax, Etmin
      REAL*8 YXP,YXETA,YXPT,YXDELR,YXPP,yxy(30)
      COMMON/SECCONV/ YXP(4,30),YXETA(30),YXPT(30),YXDELR(30,30),
     >YXPP(30)
      Character*40 Boson
      Common / Boson / Boson
      Integer iYPiece, iReadIn
      Common / YPiece / iYPiece, iReadIn

      Character*2 hDecayMode
      Common / hDM / hDecayMode

      Real*8 YXPT3,YXPT4,PTXCUT

      REAL*8 QT_V_MIN, QT_V_MAX,YB_V_MAX,Q_V_MIN,Q_V_MAX,QT_V_MXY
      COMMON/HARDCUT/ QT_V_MIN, QT_V_MAX,YB_V_MAX,
     >  Q_V_MIN,Q_V_MAX,QT_V_MXY

      REAL*8 Q_V,QT_V,YB_V,THE_STA,PHI_STA
      COMMON/CSFRAME/Q_V,QT_V,YB_V,THE_STA,PHI_STA

      Integer LTOpt
      Common / Process / LTOpt

      INTEGER I,J,NFINAL,ID_ELE,ID_NU
      REAL*8 YXDELETA,YXCOS1,YXCOS2,YXSIN1,YXSIN2,YXCOSDPH,YXDELPHI
CsB   LOGICAL NOJET - Removed on 05/10 by Pavel
      REAL*8 ZCORR,ZCORR_CUT

      Integer ACut
      Common / AAGCut / ACut
      Real*8 CosTheta
      REAL*8 XMW,XMZ,XMW2,XMZ2,AEM,GAMW,GAMZ,GAMW2,GAMZ2
      COMMON/STAND0/XMW,XMZ,XMW2,XMZ2,AEM,GAMW,GAMZ,GAMW2,GAMZ2
 
cpn2005 Isolation cone parameters
      integer iqqt
      real*8 DelRIso, ETmaxIso
      common/isolation/DelRIso, ETMaxIso, iqqt
      data iqqt /0/

      integer isub
      double precision cutM45, SubtractionShape, cutSIJ, rmij2, Aexp
      logical SubQGFlag
      Common / Subtraction / cutM45, SubtractionShape(3:4),iSub,
     >  SubQGFlag(3:4)

cpn2005 For smooth cone isolation
      real*8 ConeEtmax
      external ConeEtmax

      Logical Debug
      Data Debug /.False./

cdump
C      Debug = .True.

      If (Debug) Print*, ' Primary cuts -------------------------- '

      PASS_CUTS=.false.

CsB___Rename lab-momenta
      DO 517 J=1,12
       DO 517 I=1,4
        YXP(I,J)=PLAB(I,J)
517   CONTINUE

CsB___Construct pT's, rapidities, and pseudorapidities 
C     for the produced particles
      DO 519 J=3,12
        YXPT(J)=DSQRT(YXP(1,J)**2+YXP(2,J)**2)
        IF((YXP(4,J)-YXP(3,J)).GT.0.0 .AND.
     >  (YXP(4,J)+YXP(3,J)).GT.0.0)THEN
          YXY(J)=0.5*DLOG((YXP(4,J)+YXP(3,J))
     >    /(YXP(4,J)-YXP(3,J)))
          costheta=yxp(3,j)/sqrt(yxp(1,j)**2+yxp(2,j)**2+yxp(3,j)**2)
          yxeta(j)=0.5*dlog((1d0+costheta)/(1d0-costheta))
        ELSEIF((YXP(4,J)-YXP(3,J)).LE.0.)THEN
C          Print*, ' (YXP(4,J)-YXP(3,J)).LE.0. '
          YXY(J)=20.
          YXETA(J)=20.
        ELSEIF((YXP(4,J)+YXP(3,J)).LE.0.)THEN
C          Print*, ' (YXP(4,J)+YXP(3,J)).LE.0. '
          YXY(J)=-20.
          YXETA(J)=-20.
        END IF
519   CONTINUE

CCPY
      IF(BOSON.EQ.'HP') GOTO 333
      IF(BOSON.EQ.'HM') GOTO 333
      IF(BOSON.EQ.'HZ') GOTO 333

      IF(BOSON.EQ.'H0') THEN
        IF(hDecayMode.EQ.'ZZ' .OR. hDecayMode.EQ.'WW') THEN
          CALL PRIMARY_H0(PASS_CUTS)
          RETURN
        ELSEIF(hDecayMode.EQ.'AA') THEN
          GOTO 10
        ENDIF
      ENDIF

10    CONTINUE

C_____INCLUSIVE KINEMATICS - Removed on 05/10 by Pavel
C     IF( nPart.Eq.21 .Or. nPart.Eq.22 .Or. nPart.Eq.23 .Or.
C    .    nPart.Eq.24) THEN
CsB     When comparing the real perturbative piece to the asymptotic piece
C       we do not want to make different cuts on them.
C       NOJET=.False.
C       NFINAL=5
C     ELSE
C       NOJET=.TRUE.
C       NFINAL=4
C     ENDIF

CCPY ID_NU refers to fermion and ID_ELE to anti-fermion.
C For example, in Z -> e- e+, ID_NU refers to e-.
C In W^+ -> \nu e^+ , ID_ELE refers to e^+.
C In W^- -> e^- {\bar\nu}  , due to the variable "sign=-1" included in 
C resbosm.f , ID_ELE refers to e^-.
C
      ID_NU=3
      ID_ELE=4

CsB___Construct pT3&4
      DO 513 I=1,4
        YP12(I)=0.D0
        DO 515 J=3,4
          YP12(I)=YP12(I)+PLAB(I,J)
515     CONTINUE
513   CONTINUE

      If (Debug) Print*,' PLAB(I,3) = ',(PLAB(I,3),I=1,4) 
      If (Debug) Print*,' PLAB(I,4) = ',(PLAB(I,4),I=1,4) 

CsB Construct pT3&4 (magnitude), M3&4 and MT3&4
      YPT12 = DSQRT(YP12(1)**2+YP12(2)**2)
      YM12  = DSQRT(YP12(4)**2-YP12(3)**2-YPT12**2)
C      YMT12 = DSQRT(YPT12 + YM12)
C      Print*, ' YPT12,YM12 = ', YPT12,YM12

CsB___VB level cuts read from the input file
C   --- INVARIANT-MASS3&4 range ---
      If (Debug) Print*, YM12MIN, '< M12=', YM12, '<', YM12MAX
      IF(YM12.GT.YM12MAX .OR. YM12.LT.YM12MIN) RETURN
C   --- QT 3&4 Maximum ---
      If (Debug) Print*, ' QT34 = ', YpT12, ' < ', QTMax
      IF(YpT12.GT.QTMax) RETURN

CsB___Among others: to study the high pT perturbative contribution
C   to the low pT individual photon distribution
C   in AA & AG case.
C   Cut on
C   --- QT 3&4 Minimunm ---
      If (Debug) Print*, ' QT34 = ', YpT12, ' > ', QTMin
      If (YPT12.LT.QTMin) Return

      if (iqqt.eq.1.and.ym12.lt.ypt12) return
      if (iqqt.eq.2.and.ym12.gt.ypt12) return

C___HARD-WIRED CUTS ON VECTOR-BOSON:
CsB QT_V_MXY is read from the grid.
C   The last point of the Y piece grid in Q_T is read from the grid.
C     --- PT12 maximum ---
      If (Debug) Print*, ' YXPT(9) =', YXPT(9),'<', QT_V_MXY
      IF(YXPT(9).GT.QT_V_MXY)  RETURN
C     --- Y12 maximum ---
      If (Debug) Print*, ' ABS(YXY(9)) < ', YB_V_MAX
      IF(ABS(YXY(9)).GT.YB_V_MAX)  RETURN
C     --- INVARIANT-MASS12 range ---
      If (Debug) Print*, Q_V_MIN, ' < YM12 < ', Q_V_MAX
      IF(YM12.LT.Q_V_MIN .OR. YM12.GT.Q_V_MAX)  RETURN

CCPY May 2015
CPN   Cuts on the rapidity of the vector boson 
      if (yxy(9).lt.YWLCUT.or.yxy(9).gt.YWUCUT) return

CsB   pT's of 3&4
      If (BOSON.EQ.'AA' .OR. BOSON.EQ.'AG') THEN
      	If (PTECUT.LT.PTNCUT) THEN ! make sure that PTECUT > or = PTNCUT
          PTXCUT = PTECUT
          PTECUT = PTNCUT
          PTNCUT = PTXCUT
        End If
        YXPT3 = Max(YXPT(ID_ELE),YXPT(ID_NU))
        YXPT4 = Min(YXPT(ID_ELE),YXPT(ID_NU))
      ELSE
CCPY     YXPT3 = YXPT(ID_ELE)
C        YXPT4 = YXPT(ID_NU)
        YXPT4 = YXPT(ID_ELE)
        YXPT3 = YXPT(ID_NU)
      END IF

C_____CUTS:

      IF( iBeam.Eq.-2 ) THEN

CsB Cuts on individial photons...
CsB ... pT's:
C FOR WA70 data:
C        IF(YXPT(3).LT.MIN(PTECUT,PTNCUT)) RETURN
C        IF(YXPT(4).LT.MIN(PTECUT,PTNCUT)) RETURN
C        If (Debug)
C     >     Print*, ' MAX(YXPT(3),YXPT(4) > ', MIN(PTECUT,PTNCUT)
C        IF(MAX(YXPT(3),YXPT(4)).LT.MAX(PTECUT,PTNCUT)) RETURN
C For E706 Data
        If (Debug) Print*, ' YXPT(3) > ', PTECUT
        If (Debug) Print*, ' YXPT(4) > ', PTNCUT
        IF(YXPT(3).LT.PTECUT) RETURN
        IF(YXPT(4).LT.PTNCUT) RETURN
C FOR WA70 DATA
CsB ... rapidities for asymmetric y cuts (same on both photons):
C For WA70 and E706 data
C Reminder: set YECUT in data file to lowest rapidity value (e.g. -1.0)
C and YNCUT to highest rapidity value (e.g. +1.25)
        IF((YXETA(3).LT.YECUT).OR.(YXETA(3).GT.YNCUT))  RETURN
        IF((YXETA(4).LT.YECUT).OR.(YXETA(4).GT.YNCUT))  RETURN

      ELSE ! iBeam.Eq.-1,0,1,2

        if (Boson.eq.'AA'.or.Boson.eq.'AG'.or.Boson.Eq.'ZZ'
     >    .or.Boson.eq.'ZG'.or.Boson.eq.'H0') then 
cpn2006   Asymmetric cuts on pT and rapidities of the photons
          IF(YXPT3 .LT. PTECUT)  RETURN
          IF(YXPT4 .LT. PTNCUT)  RETURN
          IF(ABS(YXETA(ID_ELE)) .GT. YECUT)  RETURN
          IF(ABS(YXETA(ID_NU)) .GT. ynCut)  RETURN

CCPY ID_NU refers to fermion and ID_ELE to anti-fermion.
C For example, in Z -> e- e+, ID_NU refers to e-.

CCPY Feb 2013
        else if (Boson. Eq. 'W+' .Or. Boson. Eq. 'W-') then
CCPY ID_NU refers to fermion and ID_ELE to anti-fermion.
C For example, in Z -> e- e+, ID_NU refers to e-.
C In W^+ -> \nu e^+ , ID_ELE refers to e^+.
C In W^- -> e^- {\bar\nu}  , due to the variable "sign=-1" included in 
C resbosm.f , ID_ELE refers to e^-.
C
CPN 31.01.2001 PTNCUT and YNCUT are rewired to serve as upper cuts 
CPN on the electron's momentum and rapidity

C For Cuts(1): PTECUT,YECUT,DELR_CUT,PTNCUT,YnCut
C For Cuts(1): pTlmin,ymin,DelR,pTlmax,ymax

          IF(YXPT(ID_ELE) .LT. PTECUT.OR.YXPT(ID_ELE) .GT. PTNCUT)
     >      RETURN

CCPY May 2015
C          IF(YXETA(ID_ELE) .LT. YECUT.OR.YXETA(ID_ELE) .GT. YNCUT)
C     >      RETURN
          IF(ABS(YXETA(ID_ELE)) .LT. YECUT .OR.
     > ABS(YXETA(ID_ELE)).GT. YNCUT) RETURN

          IF(YXPT(ID_NU) .LT. ETMin) return

        else                    !for Z0, A0 (and DY) Bosons    

CPN 31.01.2001 PTNCUT and YNCUT are rewired to serve as upper cuts 
CPN on the electron's momentum and rapidity
          IF(YXPT(ID_ELE) .LT. PTECUT.OR.YXPT(ID_ELE) .GT. PTNCUT)
     >      RETURN
          IF(YXETA(ID_ELE) .LT. YECUT.OR.YXETA(ID_ELE) .GT. YNCUT)
     >      RETURN
          IF(YXPT(ID_NU) .LT. PTECUT.OR.YXPT(ID_NU) .GT. PTNCUT)
     >      RETURN
          IF(YXETA(ID_NU) .LT. YECUT.OR.YXETA(ID_NU) .GT. YNCUT)
     >      RETURN

        endif                   !Boson
      ENDIF ! iBeam.Eq.-2

CCPY Feb 2010: LTOPT = -1 (for Leading order calcaultion) 
      If (LTOpt.Eq.-1) GOTO 333


CsB___Construct DeltaPhi and DeltaR
        YXDELETA=YXETA(3)-YXETA(4)
        YXCOS1=YXP(1,3)/YXPT(3)
        YXSIN1=YXP(2,3)/YXPT(3)
        YXCOS2=YXP(1,4)/YXPT(4)
        YXSIN2=YXP(2,4)/YXPT(4)
        YXCOSDPH=YXCOS1*YXCOS2+YXSIN1*YXSIN2
        If(ABS(YXCOSDPH).LE.1.0D0) THEN
          YXDELPHI=ACOS(YXCOSDPH)
        ELSE
          If ((nPart.Eq.13 .or. nPart.Eq.14)
     .        .and. ABS(YXCOSDPH).LE.1.000001D0) then
CsB For the delta sigma piece ABS(YXCOSDPH) frequently (in ~20% of the events)
C   larger than 1.d0, by a small numerical offset. This value can not be
C   fed into ACOS, because it will crash. To avoid loosing these
C   events, we relax the numerical accuracy.
C            Print*, ' ABS(YXCOSDPH) = ', ABS(YXCOSDPH)
            YXDELPHI=ACOS(YXCOSDPH/1.000001d0)
            Goto 200
          End If
          RETURN
        ENDIF
 200    Continue
        YXDELR(3,4)=DSQRT(YXDELETA**2+YXDELPHI**2)
        YXDELR(4,3)=YXDELR(3,4)

CPN Cut on the "leptonic" transverse mass, which is defined in J. Smith,
C   W. L. van Neerven, and J.A.M. Vermaseren, PRL 50, 1738 (1983) and 
C   is not the same as sqrt(Q^2 + Q_T^2)
        if (Boson.Eq.'W+'.or.Boson.Eq.'W-'.or.Boson.Eq.'Z0'
     >    .or.Boson.eq.'A0') then
          TRANM=2.0*YXPT(ID_ELE)*YXPT(ID_NU)*(1.0D0-COS(YXDELPHI))
          TRANM=SQRT(TRANM)

CPN   Cuts on the transverse mass of the vector boson
          if (TRANM.lt.TMMin .or. TRANM.gt.TMMax) return
        endif ! Boson

CsB___For AA & AG from data file
CsB______Checking small the_sta region for GG -> gamma gamma
c        if(the_sta.lt.0.05) then
c          print*,'  Q_V,QT_V,YB_V,THE_STA,PHI_STA '
c          print*,  Q_V,QT_V,YB_V,THE_STA,PHI_STA
c          Print*, 'PLAB3,4,9'
c          Print*, (PLAB(I,3), i=1,4)
c          Print*, (PLAB(I,4), i=1,4)
c          Print*, (PLAB(I,9), i=1,4)
c          Print*, 'YXDelR, YXDELETA, YXDELPHI'
c          Print*,  YXDelR(3,4), YXDELETA, YXDELPHI
c        endif

      IF (BOSON.EQ.'AA' .OR. BOSON.EQ.'AG') Then ! ==========================
        IF(IBEAM.EQ.-2) THEN ! Fixed target, like WA70 ----------------------
C         Calculate the z-distribution
          ZCORR_CUT=2.75/MAX(YXPT(3),YXPT(4))
          ZCORR=YXP(1,3)*YXP(1,4)+YXP(2,3)*YXP(2,4)
          ZCORR=-ZCORR/MAX(YXPT(3),YXPT(4))**2
C Hard-wired z-cut
C    --- Z_CORRELATION minimum ---
          If (Debug) Print*, ' ZCORR > ', ZCORR_CUT
          IF(ZCORR.LT.ZCORR_CUT) RETURN
        ELSE ! Colliders ----------------------------------------------------
C    --- DELTA R minimum ---
          If (Debug) Print*, ' YXDELR(3,4) > ', DELR_CUT
          IF(YXDELR(3,4) .LT. DELR_CUT) RETURN
CsB_______Delta phi cut
C         CDF can measure events with any Delta phi
C         separation, provided they are separated in rapidity.
C         This is designed to the CDF detector having an EM calorimeter
C         with the Delta phi = 15 degrees = 0.2617 radians.
C    --- DELTA PHI and Delta Eta minimum ---
C         If(YXDELPHI.LT.0.2617 .and. YXDelEta.LT.0.2) Return
C         If(YXDELPHI.LT.0.06) Return
C
CCPY: SEPT 2009: HARD-WIRED THE FOLLOWING TWO CUTS 
C TO AVOID SINGULARITIES IN THE QQ->AA AND GG->AA CALCAULTION
C
            If (ACut.GT.0) then
C           If (Debug) Print*, ' YXPT(9) < ', YM12
CsB         M > QT cut for resummed diphoton production: Q(3&4) > QT(3&4)
C    --- M(photon pair) < QT(photon pair) ---
CsB         !!! This hard wired cut is needed for the resummed AA and AG !!!
!           If (YM12.LE.YXPT(9) .And. nPart.Eq.11 ) Return
C    --- |Cos(theta*)| < .99 ---
!           Print*, ' THE_STA = ', THE_STA
CsB         This cut is neccessary at the LHC (on W and A), because the 
C           u/t+t/u = (1 + Cos(th*)^2)/(1 - Cos(th*)^2)
C           term is occasionally large. Vegas takes this as an indication of a
C           large integrand and groups events into the |Cos(th*)|~1 region.
C           It is verified that changing this cut to 0.98 effects the
C           cross section only in the fifth digit with experimental pT cuts 
C           in place.
cpn2005     This cut is only needed for the resummed piece
cpn            If (Abs(Cos(THE_STA)).GT.0.99) Return
            If ( Abs(Cos(THE_STA)).GT.0.99 .and.
     >        (nPart.Eq.11 .or. nPart.Eq.12) ) Return 
cpn2006     This cut is needed to avoid large cancellations in the gg box
cpn2006 diagram happening when the photons are exactly back-to-back. It
cpn2006 should be applied both to the PERT and ASY piece to ensure
cpn2006 smooth matching at Q_T -> 0
            if (yxdelphi.gt.3.14) return

CCPY: BELOW ARE THE ADDITIONAL CUTS IMPLEMENTED FOR EACH NON-ZERO ACUT
          If (ACut.Eq.1) then
CsB M > QT cut for resummed diphoton production
C    --- M(2 photon) > QT(2 photon) ---
C           Q(3&4) > QT(3&4)
C            If (Debug) Print*, ' YXPT(9) < ', YM12
C            If (YM12.LE.YXPT(9) .And. nPart.Eq.11 ) Return
C            If (YXPT(9).GT.YM12 .and. nPart.Eq.12 ) Return
             If (YXPT(9).GT.YM12) Return
          Else If (ACut.Eq.2) then
C           D(phi) > Pi/2
C           If (YXDELPHI.LT.1.57 .And. nPart.Eq.11 ) Return
CsB         If this cut is lowered (without any others active in this paragraph then
C           the integration becomes unstable.  This is because at the Tevatron
C           the resummed delta_phi dist'n becomes negative below about 1.5.
C    --- DelPhi > .9 for resummed case
            If (YXDELPHI.LT.0.9 .And. nPart.Eq.11 ) Return
          Else If (ACut.Eq.3) then
C    --- QT(3&4) < 30 GeV for resummed case
            If (YXPT(9).GT.30. .And. nPart.Eq.11 ) Return
          Else If (ACut.Eq.4) then
CCPY: THIS IS FOR THOSE CUTS ALREADY IMPLEMENTED UNDER "ACUT.GT.0)
         Else If (ACut.Eq.5) then
C    --- |Cos(phi*)| < .99 ---
CsB         This cut is here to suppress the final state collinear singularities
C           for q g -> A A q.  
            If (Abs(Cos(Phi_Sta)).GT.0.99) Return
          else
            print *,'Primary_cuts: no such acut = ',acut
          End If

            ENDIF   ! ACut.GT.0

        ENDIF ! IBEAM -------------------------------------------------------
      ENDIF ! BOSON.EQ.'AA' .OR. BOSON.EQ.'AG' ==============================

      If (BOSON.EQ.'ZZ' .OR. BOSON.EQ.'ZG') then
CCPY Sept 2009        If (ACut.Eq.1) then
        If (ACut.GT.0) then
cpn2005     This cut is only needed for the resummed piece
cpn            If (Abs(Cos(THE_STA)).GT.0.99) Return
            If ( Abs(Cos(THE_STA)).GT.0.99 .and.
     >        (nPart.Eq.11 .or. nPart.Eq.12) ) Return 
cpn2006     This cut is needed to avoid large cancellations in the gg box
cpn2006 diagram happening when the photons are exactly back-to-back. It
cpn2006 should be applied both to the PERT and ASY piece to ensure
cpn2006 smooth matching at Q_T -> 0
            if (yxdelphi.gt.3.14) return
        Else
CCPY
        End If

        Goto 333

      End If

CsB   IF(NOJET) GOTO 333

CPN2005 Introduce a cut on the final-state jet when
C       applicable (e.g. in diphoton processes) 

      If (Debug) Print*, ' npart = ',npart
      if (npart.lt.21.or.npart.gt.24) goto 333

CPN2005 Check if the jet is resolved, i.e. its invariant masses with the
C       other final-state particles are larger than the selection cut
C       cutM45. If yes, demand the jet to be angular separated from the
C       other final-state particles [Delta R(i,5)> delr_cut]
                                !perturbative piece
      nfinal = 5

CsB___Cuts on jets
      If (Debug) Print*, ' Cuts on jets '
      If (Debug) Print*, ' (YXP(J,5),J=1,4)= ',(YXP(J,5),J=1,4)
C      DO 10 I=3,NFINAL
C       IF(I.EQ.ID_ELE .OR. I.EQ.ID_NU) GOTO 10
C       IF(YXPT(I) .LT. PTMIN)  RETURN
C       IF(ABS(YXETA(I)) .GT. ETAMAX)  RETURN
C 10    CONTINUE

CsB___These two cuts were used to reproduce results of NPB339('90)38.
CCsB___Rapidity of the jet
C      If (nPart.Eq.24) then ! for parton parton -> H0 X -> Z Z X
C        If (Abs(yxEta(5)).GT.yeCut) Return
C      End If
CCsB___Z-jet angle
C      If (nPart.Eq.24) then ! for parton parton -> H0 X -> Z Z X
C        Do i = 3,4
CC         Agle between Z(pi) and jet(p5)
C          CosTheta =
C     =    (pLab(1,i)*pLab(1,5)+pLab(2,i)*pLab(2,5)+pLab(3,i)*pLab(3,5))/
C     /     Sqrt(pLab(4,i)**2-xMz**2)/pLab(4,5)
C          If (CosTheta.GT.0.7) Return
C        End Do
C      Goto 333
C      End If

CsB___Delta phi
      Do I = 3, 4
        Do J = 5, NFINAL
CPN2005   Reject events with the invariant mass M(I,J)=YMIJ outside 
C         of the region of numerical stability of the 
C         perturbative matrix element 
C         if (Boson.eq.'AG') then
C           YMIJ = DSQRT( (YXP(4,I)+YXP(4,J))**2 - (YXP(3,I)+YXP(3,J))
C    >        **2-(YXP(2,I)+YXP(2,J))**2 - (YXP(1,I)+YXP(1,J))**2 )
C           if (YMIJ.lt.0.1d0) return
C         endif                 ! Boson.eq.'AG'
CPN2005 Reject events if the jet falls in the angular cone Delta R around
C       one of the photons and E_T of the jet is larger than ETmin.
C       ETmin can depend on Delta_R and the photon transverse energy; it is
C       returned by the function ConeETmin          
          YXDELETA=YXETA(I)-YXETA(J)
          YXCOS1=YXP(1,I)/YXPT(I)
          YXSIN1=YXP(2,I)/YXPT(I)
          YXCOS2=YXP(1,J)/YXPT(J)
          YXSIN2=YXP(2,J)/YXPT(J)
          YXCOSDPH=YXCOS1*YXCOS2+YXSIN1*YXSIN2
          If (Abs(YXCOSDPH).GT.1.d0) then
            Print*, ' Unphysical event: YXCOSDPH =', YXCOSDPH
            Return
          End If
          YXDELPHI=ACOS(YXCOSDPH)
          YXDELR(I,J)=DSQRT(YXDELETA**2+YXDELPHI**2)
          YXDELR(J,I)=YXDELR(I,J)
cpn2005   If the jet is within the delta R cone for the photon, decide 
cpn2005   if it should be accepted or rejected depending on the
cpn2005   isolation condition
cpn2005   Note: the present setup only works for one final parton (nfinal=5)
          if (YXDELR(I,J) .LT. DelRIso.and.
     >      Yxpt(J).gt.ConeETmax(yxpt(i),yxdelr(i,j))) return 
          
cpn2006 For events surviving isolation, conditionally apply fragmentation
cpn2006 subtraction. For isub=12 or 13, the subtraction is smoothly turned on
cpn2006 for sufficiently light photon-quark pairs. In this case, the
cpn2006 subtraction is multiplied by a "subtraction shape function",
cpn2006 which smoothly increases from 0 to 1 with decreasing invariant
cpn2006 mass of the pair.
Cpn2006 The current implementation is for nfinal=5.

          if (Boson.eq.'AA') then
cpn2006   Defaults: no subtraction performed; the subtraction shape
cpn2006 function is a unity
            SubQGFlag(i)=.false.
            SubtractionShape(i)=1d0

            if (isub.eq.1.or.isub.eq.3) 
     >        SubQGFlag(i)=.true. !subtract uncondtionally
   
            if (isub.eq.2.or.isub.eq.4) then
                               ! subtract only if the invariant mass 
                               !of the collinear pair < constant cutM45 GeV
                               !supplied from resbos.in
              YMIJ = DSQRT( (YXP(4,I)+YXP(4,J))**2 - (YXP(3,I)+YXP(3,J))
     >          **2-(YXP(2,I)+YXP(2,J))**2 - (YXP(1,I)+YXP(1,J))**2 )
              SubQGFlag(i)=YMIJ.lt.cutM45
            endif               !isub=2 or 4

            if (isub.eq.12.or.isub.eq.13) then
                               ! subtract if the squared invariant mass 
                               !of the collinear pair < cutSIJ, with
                               ! cutSIJ determined by the isolation condition
              YMIJ2 = (YXP(4,I)+YXP(4,J))**2 - (YXP(3,I)+YXP(3,J))
     >          **2-(YXP(2,I)+YXP(2,J))**2 - (YXP(1,I)+YXP(1,J))**2

              cutSIJ=Yxpt(I)*cutM45 !minimal squared inv. mass of the collinear
                                !pair allowed by the experimental isolation
c isub=13: ETMaxIso in cutM45 denotes the ratio of max(ET_hadron)/ET_gamma
c          in the cone, not max(ET_hadron) itself
              if (isub.eq.13) cutSIJ=cutSIJ*Yxpt(I)
c              SubQGFlag(i)=yxpt(j).lt.ConeETmax(yxpt(i),yxdelr(i,j))
c     >                     .and.YMIJ2.lt.cutSIJ
                  SubQGFlag(i)=YMIJ2.lt.cutSIJ
c              SubtractionShape(i)=(abs(1d0 - YMIJ2/cutSIJ))
cpn2006 Subtraction form inspired by the DDT approximation for the
cpn2006 resummed cross section; Aexp is a free parameter determining the 
cpn2006 magnitude of the pseudo-resummed fragmentation contribution
c              rmij2 = ymij2/cutSIJ
c              Aexp  = 1d0
c              SubtractionShape(i)=1d0 + 
c     >          Aexp*log(rmij2)*exp(-0.5*Aexp*log(rmij2)**2)
            endif               !isub=12 or 13
          endif                 !Boson.eq.'AA'

        enddo                   !j=5,nfinal
      enddo                     !i=3,4
        
 333  CONTINUE
      
      PASS_CUTS=.TRUE.
      If (Debug) Print*, ' Passed primary-cuts ... '
      
      RETURN
      END                       ! PRIMARY_CUTS_ALL 


C -------------------------------------------------------------------------
      SUBROUTINE PRIMARY_H0(PASS_CUTS)
C -------------------------------------------------------------------------
      IMPLICIT NONE
      LOGICAL PASS_CUTS

      REAL*8 W,QCDL,X1X2MIN
      INTEGER*4 NJET,N_PRIMARY_PART,NF,I_STR,IBEAM,NPART,NEVENT,NORDER,
     ,          N_SECOND_PART
      COMMON/AA/W,QCDL,X1X2MIN,NJET,N_PRIMARY_PART,NF,I_STR,IBEAM,NPART,
     &  NEVENT,NORDER,N_SECOND_PART
      REAL*8 SEC_PLAB,PLAB,PCEM,WHAT
      COMMON/DECAY_MOMENTA/SEC_PLAB(4,50)
      COMMON/MOM/PLAB(4,12),PCEM(4,12),WHAT
      REAL*8 ETOTMN,PTMIN,ETAMAX,ETASEP
      INTEGER*4 IAPPROX
      COMMON/CUTVAL/ETOTMN,PTMIN,ETAMAX,ETASEP,IAPPROX

      REAL*8 PTECUT,YECUT,DELR_CUT,PTNCUT,yNCut,YWLCUT,YWUCUT
      COMMON/MYCUT/PTECUT,YECUT,DELR_CUT,PTNCUT,yNCut,YWLCUT,YWUCUT
      REAL*8 YM12,YPT12,YP12,ymij, ymij2, tranm
      COMMON/MY1/YM12,YPT12,YP12(4)
      REAL*8 YM12MIN,YM12MAX,QTMin,QTMax, TMmin, TMMax, Etmin
      COMMON/MYYM12/YM12MIN,YM12MAX,QTMin,QTMax, TMMin, TMMax, Etmin
      REAL*8 YXP,YXETA,YXPT,YXDELR,YXPP,yxy(30)
      COMMON/SECCONV/ YXP(4,30),YXETA(30),YXPT(30),YXDELR(30,30),
     >YXPP(30)
      Character*40 Boson
      Common / Boson / Boson
      Integer iYPiece, iReadIn
      Common / YPiece / iYPiece, iReadIn

      Character*2 hDecayMode
      Common / hDM / hDecayMode

      Real*8 YXPT3,YXPT4,PTXCUT

      REAL*8 QT_V_MIN, QT_V_MAX,YB_V_MAX,Q_V_MIN,Q_V_MAX,QT_V_MXY
      COMMON/HARDCUT/ QT_V_MIN, QT_V_MAX,YB_V_MAX,
     >  Q_V_MIN,Q_V_MAX,QT_V_MXY

      REAL*8 Q_V,QT_V,YB_V,THE_STA,PHI_STA
      COMMON/CSFRAME/Q_V,QT_V,YB_V,THE_STA,PHI_STA

      INTEGER I,J,NFINAL,ID_ELE,ID_NU
      REAL*8 YXDELETA,YXCOS1,YXCOS2,YXSIN1,YXSIN2,YXCOSDPH,YXDELPHI
CsB   LOGICAL NOJET - Removed on 05/10 by Pavel
      REAL*8 ZCORR,ZCORR_CUT

      Integer ACut
      Common / AAGCut / ACut
      Real*8 CosTheta
      REAL*8 XMW,XMZ,XMW2,XMZ2,AEM,GAMW,GAMZ,GAMW2,GAMZ2
      COMMON/STAND0/XMW,XMZ,XMW2,XMZ2,AEM,GAMW,GAMZ,GAMW2,GAMZ2
 
      Logical Debug
      Data Debug /.False./
C      Debug = .True.

C================================
      REAL*8 YP_H
      DIMENSION YP_H(4)
      REAL*8 WPT2,W_RAP,W_PT,W_MASS,YP_W
      DIMENSION YP_W(4)
      COMMON/W_CUT/ WPT2,W_RAP,W_PT,W_MASS,YP_W
CCHEN! LOCAL VARIABLES FOR ACCEPTANCE STUDY IN HIGGS PROJECT MAY 2008
      REAL*8 PTL1,PTL2,MLL,MISET,COSLL,PHI_L1,PHI_L2,PHI_MIS,MT1,MT2,
     >       MT_MIN,YL1,YL2,PTSUM,COSPHILL,DPHILL,PLLTOT,PHI_LLT,
     >       PTLL,MT_H,COS1,COS2,SIN1,SIN2,CMIS,SMIS,CLL,SLL,THE_LL,
     >       DELTA_R,DY,PNNTOT
      DIMENSION PLLTOT(4),PNNTOT(4)
      COMMON/H0_CUT2/ PTL1,PTL2,MLL,MISET,COSLL,PHI_L1,PHI_L2,PHI_MIS,
     >       MT1,MT2,PNNTOT,
     >       MT_MIN,YL1,YL2,PTSUM,COSPHILL,DPHILL,PLLTOT,PHI_LLT,
     >       PTLL,MT_H,COS1,COS2,SIN1,SIN2,CMIS,SMIS,CLL,SLL,THE_LL,
     >       DELTA_R,DY
      INTEGER*4 MACHINE
      REAL*4 CHEN1,CHEN2
      LOGICAL USECUT
C================================

c THIS IS FOR (hDecayMode.EQ.'ZZ' .OR. hDecayMode.EQ.'WW')

      If (Debug) Print*, ' Primary cuts for H0 ----------'

      PASS_CUTS=.false.

CCPY FOR H -> VV -> 4*F      
CCPY____Mass of 3,4,5,6
        DO I=1,4
          YP_H(I)=YXP(I,3)+YXP(I,4)+YXP(I,5)+YXP(I,6)
        ENDDO
        WPT2=YP_H(1)**2+YP_H(2)**2
        W_PT=SQRT(WPT2)
        W_MASS=SQRT(YP_H(4)**2-YP_H(3)**2-WPT2)
C_______RAPIDITY OF H
        IF((YP_H(4)-YP_H(3)).GT.0.0 .AND.
     >  (YP_H(4)+YP_H(3)).GT.0.0)THEN
          W_RAP=0.5*DLOG((YP_H(4)+YP_H(3))
     >    /(YP_H(4)-YP_H(3)))
        ELSEIF((YP_H(4)-YP_H(3)).LE.0.)THEN
          W_RAP=20.
        ELSEIF((YP_H(4)+YP_H(3)).LE.0.)THEN
          W_RAP=-20.
        ENDIF
CCHEN! INPUT CUTS FOR ACCEPTANCE STUDY
CCHEN! LET'S FIX ID FOR CHARGED LEPTONS TO BE 4,5
C_______RAPIDITY OF LEPTONS
        IF((YXP(4,4)-YXP(3,4)).GT.0.0 .AND.
     >  (YXP(4,4)+YXP(3,4)).GT.0.0)THEN
          YL1=0.5*DLOG((YXP(4,4)+YXP(3,4))
     >    /(YXP(4,4)-YXP(3,4)))
        ELSEIF((YXP(4,4)-YXP(3,4)).LE.0.)THEN
          YL1=20.
        ELSEIF((YXP(4,4)+YXP(3,4)).LE.0.)THEN
          YL1=-20.
        ENDIF
         IF((YXP(4,5)-YXP(3,5)).GT.0.0 .AND.
     >  (YXP(4,5)+YXP(3,5)).GT.0.0)THEN
          YL2=0.5*DLOG((YXP(4,5)+YXP(3,5))
     >    /(YXP(4,5)-YXP(3,5)))
        ELSEIF((YXP(4,5)-YXP(3,5)).LE.0.)THEN
          YL2=20.
        ELSEIF((YXP(4,5)+YXP(3,5)).LE.0.)THEN
          YL2=-20.
        ENDIF
       PTL1=SQRT(YXP(1,4)**2 + YXP(2,4)**2)
       PTL2=SQRT(YXP(1,5)**2 + YXP(2,5)**2)
       MLL=(YXP(4,5)+YXP(4,4))**2-(YXP(3,5)+YXP(3,4))**2
     >     -(YXP(2,5)+YXP(2,4))**2-(YXP(1,5)+YXP(1,4))**2
       MLL = SQRT(ABS(MLL))
       MISET=SQRT(((YXP(2,3)+YXP(2,6))**2+(YXP(1,3)+YXP(1,6))**2))
       COSLL=(YXP(1,4)*YXP(1,5)+YXP(2,4)*YXP(2,5)+YXP(3,4)*YXP(3,5))/
     >       (SQRT(YXP(1,4)**2 + YXP(2,4)**2+YXP(3,4)**2)*
     >        SQRT(YXP(1,5)**2 + YXP(2,5)**2+YXP(3,5)**2))
       THE_LL = ACOS(COSLL)
       COS1=YXP(1,4)/PTL1
       SIN1=YXP(2,4)/PTL1
       COS2=YXP(1,5)/PTL2
       SIN2=YXP(2,5)/PTL2
       CMIS=(YXP(1,3)+YXP(1,6))/MISET
       SMIS=(YXP(2,3)+YXP(2,6))/MISET
       MT1=2.0*PTL1*MISET*(1-(COS1*CMIS+SIN1*SMIS))
       MT1=SQRT(MT1)
       
       MT2=2.0*PTL2*MISET*(1-(COS2*CMIS+SIN2*SMIS))
       MT2=SQRT(MT2)
       MT_MIN=MIN(MT1,MT2)
       DO I = 1,4
       PLLTOT(I)=YXP(I,4)+YXP(I,5)
       ENDDO
       DO I = 1,4
       PNNTOT(I)=YXP(I,3)+YXP(I,6)
       ENDDO
       PTLL= SQRT(PLLTOT(1)**2+PLLTOT(2)**2)     
       CLL=PLLTOT(1)/PTLL 
       SLL=PLLTOT(2)/PTLL 
       MT_H=2.0*PTLL*MISET*(1-(CLL*CMIS+SLL*SMIS))
       MT_H=SQRT(MT_H)
       PTSUM=PTL1+PTL2+MISET
       COSPHILL=COS1*COS2+SIN1*SIN2
       DPHILL=ABS(ACOS(COSPHILL))
       DY=ABS(YL1-YL2)
       DELTA_R=SQRT(DPHILL**2+DY**2)
       
CCHEN!********* SWITCH ON AND OFF THE CUTS,USING DZERO'S PAPER FOR TEVATRON*****
       USECUT = .FALSE.
CCPY       USECUT = .TRUE.       
       MACHINE= 1 ! Tevatron Run-2
CCPY       MACHINE=2  !LHC
       IF(USECUT .AND. MACHINE .EQ. 1) THEN
         IF(PTL1.LE.10.0.OR. PTL2.LE. 10.0) RETURN     
         IF(ABS(YL1).GE. 3.0 .OR. ABS(YL2).GE. 3.0) RETURN
         IF(MISET.LE. 20.0) RETURN
         IF(MLL .LE.15.0) RETURN
         IF(MT_MIN .LE. 20.0) RETURN
         IF(DPHILL .GE. 2.5) RETURN
         IF(DELTA_R .LE.0.5) RETURN
       ELSEIF(USECUT .AND. MACHINE .EQ. 2) THEN
         IF(PTL1.LE.20.0.OR. PTL2.LE. 20.0) RETURN     
         IF(ABS(YL1).GE. 2.5 .OR. ABS(YL2).GE. 2.5) RETURN
         IF(MISET.LE. 40.0) RETURN
         IF(MLL .GE. 80.0) RETURN
         IF(MT_H .LE. 130.0 .OR. MT_H .GE. 160.0) RETURN
         IF(THE_LL .GE. 0.9) RETURN
         IF(ABS(YL1-YL2).GE.1.5) RETURN
         IF(DPHILL .GE. 1.0) RETURN
         IF(DELTA_R .LE.0.5) RETURN
       ENDIF

      PASS_CUTS=.TRUE.
      If (Debug) Print*, ' Passed primary-cuts for H0 ...'
      
      RETURN
      END                       ! PRIMARY_H0



C -------------------------------------------------------------------------
      SUBROUTINE SECONDARY_CUTS(PASS_CUTS)
C -------------------------------------------------------------------------
      IMPLICIT NONE
      LOGICAL PASS_CUTS
      PASS_CUTS=.TRUE.
      RETURN
      END

C -----------------------------------------------------------------------------
      SUBROUTINE PASS_Output(WT,WT1)
C -----------------------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 WT,WT1

      REAL*8 PTECUT,YECUT,DELR_CUT,PTNCUT,yNCut,YWLCUT,YWUCUT
      COMMON/MYCUT/PTECUT,YECUT,DELR_CUT,PTNCUT,yNCut,YWLCUT,YWUCUT
      REAL*8 YM12,YPT12,YP12
      COMMON/MY1/YM12,YPT12,YP12(4)
      REAL*8 YM12MIN,YM12MAX,QTMin,QTMax, TMmin, Tmmax, ETmin
      COMMON/MYYM12/YM12MIN,YM12MAX,QTMin,QTMax, TMmin, Tmmax, ETmin
      REAL*8 YXP,YXETA,YXPT,YXDELR,YXPP
      COMMON/SECCONV/ YXP(4,30),YXETA(30),YXPT(30),YXDELR(30,30),
     >YXPP(30)
      REAL*8 YLUMI
      COMMON/CHIST/ YLUMI

      REAL*4 A
      COMMON /GBOOK/ A(10000)
      REAL*8 ZERO,WEIGHT,Weight1
      REAL*4 ZERO4,WEIGHT4,WT4,Dummy

      Real*8 v4e, v4n
      Dimension v4e(4), v4n(4)

      INTEGER*4 ID_NU,ID_ELE,I
      REAL*8 PI,PHI_E,PHI_N,DELPHI,TRANM
      REAL*8 WPT2,W_RAP,W_PT,W_MASS,YP_W
      DIMENSION YP_W(4)
CsB      NT_ELEMENTS  number of elements to fill
C-       EVENT_VECTOR returned values to enter into the ntuple
C-                 the order should match the order given for the tags
      Integer nt_id, lunit, nt_elements
      Character*8 nt_subdir
      Character*80 nt_file_name
      Common / NTPL / nt_id, lunit, nt_subdir, nt_file_name, nt_elements
 

      integer max_siz
      parameter(max_siz=50)
      real*4 event_vector(max_siz)
      common/evec/event_vector
      Character*80 chtitl
      Character*8  chtags(max_siz)
      common /NT_TAGS/chtitl, chtags


      Character*8 Output
      Common / Output / Output

      Real*8 CTS
      REAL*8 Q_V,QT_V,YB_V,THE_STA,PHI_STA
      COMMON/CSFRAME/Q_V,QT_V,YB_V,THE_STA,PHI_STA

      Real*8 ZCorr, Del_y

      Character*40 Boson
      Common / Boson / Boson
      Character*2 hDecayMode
      Common / hDM / hDecayMode

      Real*8 Mass35, Mass45

      REAL*8 W,QCDL,X1X2MIN
      INTEGER*4 NJET,N_PRIMARY_PART,NF,I_STR,IBEAM,NPART,NEVENT,NORDER,
     ,          N_SECOND_PART
      COMMON/AA/W,QCDL,X1X2MIN,NJET,N_PRIMARY_PART,NF,I_STR,IBEAM,NPART,
     &  NEVENT,NORDER,N_SECOND_PART

      Real*8 x_1,x_2, pT1,pT2

      Real*8 SEC_PLAB
      COMMON/DECAY_MOMENTA/SEC_PLAB(4,50)
      
      Real*8 SEC_PT,SEC_ETA
      COMMON/SEC_VARS/SEC_PT(10),SEC_ETA(10)

      Data Pi / 3.1415926535897932 /, Zero/ 0.d0 /, Zero4/ 0.0 /

      REAL*8 MAX_WGT,NEW_MAX_WGT,WGT_FACTOR
      INTEGER*4 NBINIT, iwgt
      COMMON/UNWEIGHT/MAX_WGT,NEW_MAX_WGT,WGT_FACTOR,iwgt,NBINIT

      real*8 weight_new, weight_old, xweight_old
      integer icut_mismatch
C yfu add a tag for qqbar or qbarq event
      integer myEventType
      CHARACTER*40 envEventType

      common/reweight/weight_new, weight_old, icut_mismatch

      REAL*8 cts_PHASE,E1E2_MASS,E1E2_TMASS,SIGN_E1E2_Z,
     > E1_PLUS,E1_MINUS,E2_PLUS,E2_MINUS,COS_THETA_STAR,
     > EPS
      REAL*8 SIN_THETA_STAR,COS_PHI_STAR,CPS_PHASE,CPS,SPS_PHASE,SPS

CCPY
      IF(BOSON.EQ.'H0') THEN
        IF(hDecayMode.EQ.'ZZ' .OR. hDecayMode.EQ.'WW') THEN
         CALL PASS_H0(WT,WT1)
         RETURN
        ENDIF
      ENDIF

CCPY
      WEIGHT=WT*YLUMI
      Weight1 = WT1*YLumi

cpn2008 for reweighting
      if (iwgt.eq.-2) then
        xweight_old =  sign(max(abs(weight_old),1e-14),weight_old)
        weight=weight*weight_new/xweight_old
        weight1=weight1*weight_new/xweight_old
      endif                     !iwgt
      

CsB___Convert to *4
      WT4 = WT
      Weight4 = Weight

CsB___For any process: fermion=3, anti-fermion=4
C     E.g. for W+ production: nu=3, e+=4
      ID_NU=3
      ID_ELE=4

      If (nPart.Eq.21.or.npart.eq.22) then
CsB_____Mass of 3&5
        DO I=1,4
          YP_W(I)=YXP(I,3)+YXP(I,5)
        ENDDO
        WPT2=YP_W(1)**2+YP_W(2)**2
        W_PT=SQRT(WPT2)
        MASS35=SQRT(YP_W(4)**2-YP_W(3)**2-WPT2)

CsB_____Mass of 4&5
        DO I=1,4
          YP_W(I)=YXP(I,4)+YXP(I,5)
        ENDDO
        WPT2=YP_W(1)**2+YP_W(2)**2
        W_PT=SQRT(WPT2)
        MASS45=SQRT(YP_W(4)**2-YP_W(3)**2-WPT2)
      End if

CsB___Mass of 3&4
      DO I=1,4
        YP_W(I)=YXP(I,ID_ELE)+YXP(I,ID_NU)
      ENDDO
      WPT2=YP_W(1)**2+YP_W(2)**2
      W_PT=SQRT(WPT2)
      W_MASS=SQRT(YP_W(4)**2-YP_W(3)**2-WPT2)

C_______RAPIDITY OF 3&4
        IF((YP_W(4)-YP_W(3)).GT.0.0 .AND.
     >  (YP_W(4)+YP_W(3)).GT.0.0)THEN
          W_RAP=0.5*DLOG((YP_W(4)+YP_W(3))
     >    /(YP_W(4)-YP_W(3)))
        ELSEIF((YP_W(4)-YP_W(3)).LE.0.)THEN
          W_RAP=20.
        ELSEIF((YP_W(4)+YP_W(3)).LE.0.)THEN
          W_RAP=-20.
        ENDIF

CsB_____Rapidity difference: y3-y4
C        Del_y = YXETA(ID_ELE)-YXETA(ID_NU)
        If (YXPT(ID_ELE).GT.YXPT(ID_NU)) then
          Del_y =  YXETA(ID_ELE)-YXETA(ID_NU)
        Else
          Del_y = -YXETA(ID_ELE)+YXETA(ID_NU)
        End If

CsB___Angles of decay products
      PHI_E=DATAN2(YXP(2,ID_ELE),YXP(1,ID_ELE))
      IF(PHI_E.LT.0.D0) PHI_E=PHI_E+2.0D0*PI
      PHI_N=DATAN2(YXP(2,ID_NU),YXP(1,ID_NU))
      IF(PHI_N.LT.0.D0) PHI_N=PHI_N+2.0D0*PI

CsB___Delta phi: angle between 3 and 4 (in the Lab frame)
C      DELPHI=PHI_E-PHI_N
      DELPHI=ACos(Cos(PHI_E-PHI_N))

CsB___Folding the Delta Phi distribution
      If (Boson.EQ.'AA' .OR. Boson.EQ.'AG' .OR.
     .    Boson.EQ.'ZZ' .OR. Boson.EQ.'ZG' .Or.
     .    Boson.Eq.'H0') then ! Identical particles in the final state -----
         If(DELPHI.GT.PI) DELPHI=DELPHI-PI
      End If

C_____TRANSVERSE MASS OF 3&4 
Cpn Defined in J. Smith, W. L. van Neerven, and J.A.M. Vermaseren, and
Cpn not the same as sqrt(Q^2 + Q_T^2)
      TRANM=2.0*YXPT(ID_ELE)*YXPT(ID_NU)*(1.0D0-COS(DELPHI))
      TRANM=SQRT(TRANM)

C_____Z-distribution for 3 and 4
      ZCORR=YXP(1,ID_ELE)*YXP(1,ID_NU)+YXP(2,ID_ELE)*YXP(2,ID_NU)
      ZCORR=-ZCORR/MAX(YXPT(ID_NU),YXPT(ID_ELE))**2

C_yfu____Cos(the_sta) in CS frame (14)
      cts_PHASE = Cos(the_sta)

CCPY June 19, 2020
CCPY calculate cos(theta_sta) and cso(phi_sta) in the Collin-Soper frame

C MASS
       E1E2_MASS=W_MASS
C TRANSVERSE MASS
       E1E2_TMASS= SQRT(W_MASS**2+W_PT**2)
C Z-COMPONET
       SIGN_E1E2_Z= YP_W(3)/ABS(YP_W(3))
       
       E1_PLUS= (YXP(4,ID_NU)+YXP(3,ID_NU))/SQRT(2.0)
       E1_MINUS= (YXP(4,ID_NU)-YXP(3,ID_NU))/SQRT(2.0)

       E2_PLUS= (YXP(4,ID_ELE)+YXP(3,ID_ELE))/SQRT(2.0)
       E2_MINUS= (YXP(4,ID_ELE)-YXP(3,ID_ELE))/SQRT(2.0)

       COS_THETA_STAR= 2.0*(E1_PLUS*E2_MINUS-E1_MINUS*E2_PLUS)/
     > E1E2_MASS/E1E2_TMASS

C yfu note that CTS is the cos theta defined by Z boost (13)
       CTS = SIGN_E1E2_Z*COS_THETA_STAR
C SIN_THETA_STAR > 0 for (the_Sta) between 0 adn pi.
       EPS=1.0D-8
C       SIN_THETA_STAR=DSQRT(1.0D0-COS_THETA_STAR**2)
C       IF(ABS(SIN_THETA_STAR).GT.EPS) THEN         
C         COS_PHI_STAR= 2.0*(E1_PLUS*E1_MINUS-E2_MINUS*E2_PLUS)/
C     > W_PT/E1E2_TMASS/SIN_THETA_STAR
C       ELSE
C         PRINT*,' ABS(SIN_THETA_STAR).LT.EPS '
C         CALL ExitGracefully
C       ENDIF

C_____Cos(the_sta) in CS frame, from subroutine WRESPH
       CTS_PHASE = Cos(the_sta)
       CPS_PHASE = Cos(phi_Sta)
       SPS_PHASE = Sin(phi_Sta)      
       
C=========================
CCPY June 19, 2020
       IF (IBEAM.EQ.1) THEN
          CTS = SIGN_E1E2_Z*COS_THETA_STAR
          CPS = COS_PHI_STAR
       ELSE
          CTS = COS_THETA_STAR
          CPS = COS_PHI_STAR

CCPY This SPS may be wrong. It is only correct up to an overall sign. 
C We need to use Eq.(A1) of the ResBos paper, with (p^\mu \dot Y_\mu), 
C to solve for SPS=SIN_PHI_STAR  in order to get the correct sign.
          SPS= DSQRT(1.0D0-CPS**2)           

       ENDIF 

C yfu add a tag for qqbar or qbarq event

         call getenv('myEventType',envEventType)
         IF ( envEventType(1:1).eq.' ') then
           myEventType = 0
         ELSE IF ( envEventType(1:5).eq.'uubar') then
           myEventType = 1
         ELSE IF ( envEventType(1:5).eq.'ddbar') then
           myEventType = 2
         ELSE IF ( envEventType(1:5).eq.'ssbar') then
           myEventType = 3
         ELSE IF ( envEventType(1:5).eq.'ccbar') then
           myEventType = 4
         ELSE IF ( envEventType(1:5).eq.'bbbar') then
           myEventType = 5
         ELSE IF ( envEventType(1:5).eq.'ubaru') then
           myEventType = -1
         ELSE IF ( envEventType(1:5).eq.'dbard') then
           myEventType = -2
         ELSE IF ( envEventType(1:5).eq.'sbars') then
           myEventType = -3
         ELSE IF ( envEventType(1:5).eq.'cbarc') then
           myEventType = -4
         ELSE IF ( envEventType(1:5).eq.'bbarb') then
           myEventType = -5
         ELSE IF ( envEventType(1:5).eq.'udbar') then
           myEventType = 112
         ELSE IF ( envEventType(1:5).eq.'usbar') then
           myEventType = 113
         ELSE IF ( envEventType(1:5).eq.'ubbar') then
           myEventType = 115
         ELSE IF ( envEventType(1:5).eq.'cdbar') then
           myEventType = 142
         ELSE IF ( envEventType(1:5).eq.'csbar') then
           myEventType = 143
         ELSE IF ( envEventType(1:5).eq.'cbbar') then
           myEventType = 145
         ELSE IF ( envEventType(1:5).eq.'dbaru') then
           myEventType = -112
         ELSE IF ( envEventType(1:5).eq.'sbaru') then
           myEventType = -113
         ELSE IF ( envEventType(1:5).eq.'bbaru') then
           myEventType = -115
         ELSE IF ( envEventType(1:5).eq.'dbarc') then
           myEventType = -142
         ELSE IF ( envEventType(1:5).eq.'sbarc') then
           myEventType = -143
         ELSE IF ( envEventType(1:5).eq.'bbarc') then
           myEventType = -145
         ELSE IF ( envEventType(1:5).eq.'dubar') then
           myEventType = 212
         ELSE IF ( envEventType(1:5).eq.'subar') then
           myEventType = 213
         ELSE IF ( envEventType(1:5).eq.'bubar') then
           myEventType = 215
         ELSE IF ( envEventType(1:5).eq.'dcbar') then
           myEventType = 242
         ELSE IF ( envEventType(1:5).eq.'scbar') then
           myEventType = 243
         ELSE IF ( envEventType(1:5).eq.'bcbar') then
           myEventType = 245
         ELSE IF ( envEventType(1:5).eq.'ubard') then
           myEventType = -212
         ELSE IF ( envEventType(1:5).eq.'ubars') then
           myEventType = -213
         ELSE IF ( envEventType(1:5).eq.'ubarb') then
           myEventType = -215
         ELSE IF ( envEventType(1:5).eq.'cbard') then
           myEventType = -242
         ELSE IF ( envEventType(1:5).eq.'cbars') then
           myEventType = -243
         ELSE IF ( envEventType(1:5).eq.'cbarb') then
           myEventType = -245
         END IF

C CTS_PHASE and COS_THETA_STAR should be equal (checked).
C the above CTS maybe different with CTS_PHASE and COS_THETA_STAR due to
C the boost effects, SIGN_E1E2_Z

C FOR TESTING:
C       EPS=1.0D-4
C      IF( ABS(CTS - CTS_PHASE).GT.EPS) THEN 
C        PRINT*,' IBEAM, SIGN_E1E2_Z =',IBEAM,SIGN_E1E2_Z
C        PRINT*,' CTS, CTS_PHASE =',CTS, CTS_PHASE
C        PRINT*,' ABS(CTS - CTS_PHASE).GT.EPS '
C        CALL ExitGracefully
C      ENDIF
C       IF( ABS(CPS - CPS_PHASE).GT.EPS) THEN 
C         PRINT*,' '
C         PRINT*,' IBEAM, SIGN_E1E2_Z =',IBEAM,SIGN_E1E2_Z
C         PRINT*,' CTS, CTS_PHASE =',CTS, CTS_PHASE
C         PRINT*,' CPS, CPS_PHASE =',CPS, CPS_PHASE
C         PRINT*,' ABS(CTP - CPS_PHASE).GT.EPS '
C         CALL ExitGracefully
C       ENDIF
C       IF( ABS(SPS - SPS_PHASE).GT.EPS) THEN 
C         PRINT*,' '
C         PRINT*,' IBEAM, SIGN_E1E2_Z =',IBEAM,SIGN_E1E2_Z
C         PRINT*,' CTS, CTS_PHASE =',CTS, CTS_PHASE
C         PRINT*,' CPS, CPS_PHASE =',CPS, CPS_PHASE
C         PRINT*,' SPS, SPS_PHASE =',SPS, SPS_PHASE
C         PRINT*,' ABS(STP - SPS_PHASE).GT.EPS '
C         CALL ExitGracefully
C       ENDIF


c      If (W_Mass.GT.36.8 .AND.  W_Mass.LT.38.4) then
c        Print*, 'Q,qT,y, y_l1,y_l2, Cos(the_sta), pT_1,pT_2, ZCorr'
c        Print*, W_Mass,W_PT,W_RAP,YXETA(ID_ELE),YXETA(ID_NU), CTS
c        Print*, YXPT(ID_ELE),YXPT(ID_NU), ZCorr
c      End if

C_____Momentum fractions x_i
C     x_1 = exp[ y] Q/Sqrt[S]
C     x_2 = exp[-y] Q/Sqrt[S]
      x_1 = DExp( W_Rap)*W_Mass/W
      x_2 = DExp(-W_Rap)*W_Mass/W

      If (iwgt.eq.-2) then
        event_vector(1) = Weight
      Else If (Output.Eq.'ROOTNT1'.or.Output.Eq.'PAWNT1'
     >    .or.Output.Eq.'Text1') then
cpn2008 Prepare the event variables to fill
        event_vector(1)  =  YXP(1,ID_ELE)
        event_vector(2)  =  YXP(2,ID_ELE)
        event_vector(3)  =  YXP(3,ID_ELE)
        event_vector(4)  =  YXP(4,ID_ELE)
        event_vector(5)  =  YXP(1,ID_NU)
        event_vector(6)  =  YXP(2,ID_NU)
        event_vector(7)  =  YXP(3,ID_NU)
        event_vector(8)  =  YXP(4,ID_NU)
        event_vector(9)  =  YP_W(1)
        event_vector(10) =  YP_W(2)
        event_vector(11) =  YP_W(3)
        event_vector(12) =  YP_W(4)
C yfu add cos theta of hadron level
        event_vector(13) = CTS
C yfu add cos theta of quark level
        event_vector(14) = cts_PHASE
C yfu add a tag for qqbar or qbarq event
C qqbar=1 qbarq=-1 all=0
        event_vector(15) = myEventType
        event_vector(16) = Weight
      ElseIf (Output.Eq.'ROOTNT2'.or.Output.Eq.'PAWNT2') then
        event_vector(1) = YXPT(ID_ELE)
        event_vector(2) = YXPT(ID_NU)
        event_vector(3) = YXETA(ID_ELE)
        event_vector(4) = YXETA(ID_NU)
        event_vector(5) = W_PT
        event_vector(6) = W_RAP
        event_vector(7) = W_MASS
        event_vector(8) = DELPHI
        event_vector(9) = CTS
        event_vector(10) = phi_sta
        event_vector(11) = yxdelr(3,4)
        event_vector(12) = Weight
        event_vector(13) = 0
cpn Alternative choices for the output
cpn       Invariant masses of 3 and 5, 4 and 5
c         event_vector(11) =  DSQRT( (YXP(4,3)+YXP(4,5))**2 
c     >      - (YXP(3,3)+YXP(3,5))**2-(YXP(2,3)+YXP(2,5))**2 
c     >      - (YXP(1,3)+YXP(1,5))**2 )
c          event_vector(13) =  DSQRT( (YXP(4,4)+YXP(4,5))**2 
c     >      - (YXP(3,4)+YXP(3,5))**2-(YXP(2,4)+YXP(2,5))**2 
c     >      - (YXP(1,4)+YXP(1,5))**2 )
cpn       Cosines of the azimuthal angles between 3 and 5, 4 and 5
c          event_vector(11) =  (YXP(1,3)*YXP(1,5)+YXP(2,3)*YXP(2,5))/
c     >      yxpt(3)/w_pt
c          event_vector(13) =  (YXP(1,4)*YXP(1,5)+YXP(2,4)*YXP(2,5))/
c     >      yxpt(4)/yxpt(5)

cpn Catani-Seymour's z for 3||5
c          event_vector(13)=1d0-mass45/(W_mass+mass35+mass45)
CCPY Oct 2008
      ElseIf (Output.Eq.'Text2') then
        event_vector(1) = YXPT(ID_ELE)
        event_vector(2) = YXPT(ID_NU)
        event_vector(3) = YXETA(ID_ELE)
        event_vector(4) = YXETA(ID_NU)
        event_vector(5) = W_PT
        event_vector(6) = W_RAP
        event_vector(7) = W_MASS
        event_vector(8) = DELPHI + phi_sta
        event_vector(9) = CTS
        event_vector(10) = phi_sta
        event_vector(11) = yxdelr(3,4)
        event_vector(12) = Weight
        event_vector(13) = 0
      Endif                     !iwgt.ne.-2


cpn2008 Fill the output file
      if (Output.Eq.'GBook ') then
	  Dummy = YXPT(ID_ELE)
          CALL GFill1(1,Dummy,WEIGHT4)
	  Dummy = YXPT(ID_NU)
          CALL GFill1(2,Dummy,WEIGHT4)
	  Dummy = YXETA(ID_ELE)
          CALL GFill1(3,Dummy,WEIGHT4)
	  Dummy = YXETA(ID_NU)
          CALL GFill1(4,Dummy,WEIGHT4)
	  Dummy = W_MASS
          CALL GFill1(5,Dummy,WEIGHT4)
	  Dummy = W_PT
          CALL GFill1(6,Dummy,WEIGHT4)
	  Dummy = W_RAP
          CALL GFill1(7,Dummy,WEIGHT4)
	  Dummy = TRANM
          CALL GFill1(8,Dummy,WEIGHT4)
	  Dummy = CTS
          CALL GFill1(9,Dummy,WEIGHT4)
	  Dummy = DELPHI
          CALL GFill1(10,Dummy,WEIGHT4)
	  Dummy = YXP(4,ID_ELE)
          CALL GFill1(11,Dummy,WEIGHT4)
 	  Dummy = YXP(4,ID_NU)
          CALL GFill1(12,Dummy,WEIGHT4)
	  Dummy = YP_W(3)
          CALL GFill1(13,Dummy,WEIGHT4)
	  Dummy = CTS_PHASE
          CALL GFill1(14,Dummy,WEIGHT4)
	  Dummy = phi_sta
          CALL GFill1(15,Dummy,WEIGHT4)
c          Dummy = ZCorr
c        CALL GFill1(16,Dummy,WEIGHT4)
cc        If (Abs(YXETA(ID_ELE)).GT.1.d0) then
cc 	    Dummy = YXPT(ID_ELE)
cc          CALL GFill1(17,Dummy,WEIGHT4)
cc        EndIf
          Dummy=YXDELR(3,4)
          CALL GFill1(17,Dummy,WEIGHT4)
ccpy
	  Dummy = W_PT
          CALL GFill1(18,Dummy,WEIGHT4)
          Dummy = TRANM
          CALL GFill1(19,Dummy,WEIGHT4)
	  Dummy = YXPT(ID_ELE)
          CALL GFill1(20,Dummy,WEIGHT4)
	  Dummy = YXPT(ID_NU)
          CALL GFill1(21,Dummy,WEIGHT4)
CCPY
	  Dummy = YXETA(ID_ELE)
          CALL GFill1(22,Dummy,WEIGHT4)          
	  Dummy = YXETA(ID_NU)
          CALL GFill1(23,Dummy,WEIGHT4)
	  Dummy = W_RAP
          CALL GFill1(24,Dummy,WEIGHT4)

          If (NPart.Eq.21.or.npart.eq.22) then
            Dummy=YXDELR(3,5)
            CALL GFill1(88,Dummy,WEIGHT4)
            Dummy=YXDELR(4,5)
            CALL GFill1(89,Dummy,WEIGHT4)
            Dummy=Mass35
            CALL GFill1(90,Dummy,WEIGHT4)
            Dummy=Mass45
            CALL GFill1(91,Dummy,WEIGHT4)
          End If                !Npart
        
        ElseIf (Output.Eq.'ROOTNT1'.or.Output.Eq.'ROOTNT2') then

cpn2008   Fill ROOT ntuple
          if (iwgt.eq.-2) then
            call FillNTBranch(chtags(1))
          else 
            call FillNTBranch('all')
          endif                 !iwgt

      elseif (Output.eq.'Text1'.or.Output.Eq.'Text2') then
        if(iwgt.eq.-2) then
          write (18, '(g12.5)') Weight
        else
          write (18,20) (event_vector(i),i=1,nt_elements)
 20       format (11(2x,f10.4),2(2x,g12.5))
        endif                   !iwgt
      ElseIf (Output.Eq.'StdHEP') then
        call StdHEPFill(Weight4)
      Else
        Print*, ' Problem in PASS_Output, Output: ', Output
        Stop
      EndIf

      RETURN
      END ! PASS_Output

C -----------------------------------------------------------------------------
      SUBROUTINE PASS_H0(WT,WT1)
C -----------------------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 WT,WT1

      REAL*8 PTECUT,YECUT,DELR_CUT,PTNCUT,yNCut,YWLCUT,YWUCUT
      COMMON/MYCUT/PTECUT,YECUT,DELR_CUT,PTNCUT,yNCut,YWLCUT,YWUCUT
      REAL*8 YM12,YPT12,YP12
      COMMON/MY1/YM12,YPT12,YP12(4)
      REAL*8 YM12MIN,YM12MAX,QTMin,QTMax, TMmin, Tmmax, ETmin
      COMMON/MYYM12/YM12MIN,YM12MAX,QTMin,QTMax, TMmin, Tmmax, ETmin
      REAL*8 YXP,YXETA,YXPT,YXDELR,YXPP
      COMMON/SECCONV/ YXP(4,30),YXETA(30),YXPT(30),YXDELR(30,30),
     >YXPP(30)
      REAL*8 YLUMI
      COMMON/CHIST/ YLUMI

      REAL*4 A
      COMMON /GBOOK/ A(10000)
      REAL*8 ZERO,WEIGHT,Weight1
      REAL*4 ZERO4,WEIGHT4,WT4,Dummy

      Real*8 v4e, v4n
      Dimension v4e(4), v4n(4)

      INTEGER*4 ID_NU,ID_ELE,I
      REAL*8 PI,PHI_E,PHI_N,DELPHI,TRANM

CsB      NT_ELEMENTS  number of elements to fill
C-       EVENT_VECTOR returned values to enter into the ntuple
C-                 the order should match the order given for the tags
      Integer nt_id, lunit, nt_elements
      Character*8 nt_subdir
      Character*80 nt_file_name
      Common / NTPL / nt_id, lunit, nt_subdir, nt_file_name, nt_elements
 

      integer max_siz
      parameter(max_siz=50)
      real*4 event_vector(max_siz)
      common/evec/event_vector
      Character*80 chtitl
      Character*8  chtags(max_siz)
      common /NT_TAGS/chtitl, chtags


      Character*8 Output
      Common / Output / Output

      Real*8 CTS
      REAL*8 Q_V,QT_V,YB_V,THE_STA,PHI_STA
      COMMON/CSFRAME/Q_V,QT_V,YB_V,THE_STA,PHI_STA

      Real*8 ZCorr, Del_y

      Character*40 Boson
      Common / Boson / Boson
      Character*2 hDecayMode
      Common / hDM / hDecayMode

      Real*8 Mass35, Mass45

      REAL*8 W,QCDL,X1X2MIN
      INTEGER*4 NJET,N_PRIMARY_PART,NF,I_STR,IBEAM,NPART,NEVENT,NORDER,
     ,          N_SECOND_PART
      COMMON/AA/W,QCDL,X1X2MIN,NJET,N_PRIMARY_PART,NF,I_STR,IBEAM,NPART,
     &  NEVENT,NORDER,N_SECOND_PART

      Real*8 x_1,x_2, pT1,pT2

      Real*8 SEC_PLAB
      COMMON/DECAY_MOMENTA/SEC_PLAB(4,50)
      
      Real*8 SEC_PT,SEC_ETA
      COMMON/SEC_VARS/SEC_PT(10),SEC_ETA(10)

      Data Pi / 3.1415926535897932 /, Zero/ 0.d0 /, Zero4/ 0.0 /

      REAL*8 MAX_WGT,NEW_MAX_WGT,WGT_FACTOR
      INTEGER*4 NBINIT, iwgt
      COMMON/UNWEIGHT/MAX_WGT,NEW_MAX_WGT,WGT_FACTOR,iwgt,NBINIT

      real*8 weight_new, weight_old, xweight_old
      integer icut_mismatch
      common/reweight/weight_new, weight_old, icut_mismatch

C================================
      REAL*8 YP_H
      DIMENSION YP_H(4)
      REAL*8 WPT2,W_RAP,W_PT,W_MASS,YP_W
      DIMENSION YP_W(4)
      COMMON/W_CUT/ WPT2,W_RAP,W_PT,W_MASS,YP_W
CCHEN! LOCAL VARIABLES FOR ACCEPTANCE STUDY IN HIGGS PROJECT MAY 2008
      REAL*8 PTL1,PTL2,MLL,MISET,COSLL,PHI_L1,PHI_L2,PHI_MIS,MT1,MT2,
     >       MT_MIN,YL1,YL2,PTSUM,COSPHILL,DPHILL,PLLTOT,PHI_LLT,
     >       PTLL,MT_H,COS1,COS2,SIN1,SIN2,CMIS,SMIS,CLL,SLL,THE_LL,
     >       DELTA_R,DY,PNNTOT
      DIMENSION PLLTOT(4),PNNTOT(4)
      COMMON/H0_CUT2/ PTL1,PTL2,MLL,MISET,COSLL,PHI_L1,PHI_L2,PHI_MIS,
     >       MT1,MT2,PNNTOT,
     >       MT_MIN,YL1,YL2,PTSUM,COSPHILL,DPHILL,PLLTOT,PHI_LLT,
     >       PTLL,MT_H,COS1,COS2,SIN1,SIN2,CMIS,SMIS,CLL,SLL,THE_LL,
     >       DELTA_R,DY
      INTEGER*4 MACHINE
      REAL*4 CHEN1,CHEN2
      LOGICAL USECUT
C================================


CCPY
      WEIGHT=WT*YLUMI
      Weight1 = WT1*YLumi

cpn2008 for reweighting
      if (iwgt.eq.-2) then
        xweight_old =  sign(max(abs(weight_old),1e-14),weight_old)
        weight=weight*weight_new/xweight_old
        weight1=weight1*weight_new/xweight_old
      endif                     !iwgt
      

      If (iwgt.eq.-2) then
        event_vector(1) = Weight
      Else If (Output.Eq.'ROOTNT1'.or.Output.Eq.'PAWNT1'
     >    .or.Output.Eq.'Text1') then
        event_vector(1)  =  YXP(1,4)
        event_vector(2)  =  YXP(2,4)
        event_vector(3)  =  YXP(3,4)
        event_vector(4)  =  YXP(4,4)
        event_vector(5)  =  YXP(1,5)
        event_vector(6)  =  YXP(2,5)
        event_vector(7)  =  YXP(3,5)
        event_vector(8)  =  YXP(4,5)
        event_vector(9)  =  PNNTOT(1)
        event_vector(10) =  PNNTOT(2)
        event_vector(11) =  PNNTOT(3)
        event_vector(12) =  PNNTOT(4)
        event_vector(13) =  YP_H(1)
        event_vector(14) =  YP_H(2)
        event_vector(15) =  YP_H(3)
        event_vector(16) =  YP_H(4)
        event_vector(17) =  Weight
      ElseIf (Output.Eq.'ROOTNT2'.or.Output.Eq.'PAWNT2'
     >    .or.Output.Eq.'Text2') then
        event_vector(1) = W_MASS
        event_vector(2) = W_PT
        event_vector(3) = W_RAP
        event_vector(4) = PTL1
        event_vector(5) = YL1
        event_vector(6) = PTL2
        event_vector(7) = YL2
        event_vector(8) = COSLL
        event_vector(9) = PTLL
        event_vector(10) = DPHILL
        event_vector(11) = DY
        event_vector(12) = DELTA_R
        event_vector(13) =  MISET
        event_vector(14) =  MT1
        event_vector(15) =  MT2
        event_vector(16) =  MT_H
        event_vector(17) =  Weight
      Endif

      WT4 = WT
      Weight4 = Weight

       If (Output.Eq.'ROOTNT1'.or.Output.Eq.'ROOTNT2') then
cpn2008   Fill ROOT ntuple
          if (iwgt.eq.-2) then
            call FillNTBranch(chtags(1))
          else 
            call FillNTBranch('all')
          endif                 !iwgt

      elseif (Output.eq.'Text1'.or.Output.Eq.'Text2') then
        if(iwgt.eq.-2) then
          write (18, '(g12.5)') Weight
        else
          write (18,20) (event_vector(i),i=1,nt_elements)
 20       format (11(2x,f10.4),2(2x,g12.5))
        endif                   !iwgt
      Else
        Print*, ' Problem in PASS_Output, Output: ', Output
        Stop
      Endif                    


      RETURN
      END ! PASS_H0


CsB * *******************************************************************

      Subroutine ReadIn(Filename)
C
C     The following notation apply:
C     For the CSS (in case of NLO run: Asymptotic) piece
C     X1A (X2A, X3A) contains the qT (y, Q) vector
C     iD1+1 (iD2+1, iD3+1) is the number of qT (y, Q) values in the grid
C     The 1st index of YAS and YAA belong to qT, the 2nd to y, the 3rd to Q.
C     For the Y piece:
C     X1Y (X2Y, X3Y) contains the qT (y, Q) vector
C     iY1+1 (iY2+1, iY3+1) is the number of qT (y, Q) values in the grid
C     The 1st index of YAY belong to qT, the 2nd to y, the 3rd to Q.
C
C     iDD1 (iDD2, iDD3) is the number of qT (y, Q)
C          values used at the interpolation
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      Real*8 ECMC
      Character*78 dummy
      Character*40 Boson
      Character*145 Filename
CsB   The dimensions below give the maximal possible gridsize.
      Parameter (MD1 = 200, MD2 = 300, MD3 = 100)
      Dimension X1A(0:MD1), X2A(0:MD2), X3A(0:MD3)
      Dimension YAS(0:MD1,0:MD2,0:MD3), YAA(0:MD1,0:MD2,0:MD3)
      Dimension X1Y(0:MD1), X2Y(0:MD2), X3Y(0:MD3)
      Dimension YAY(0:MD1,0:MD2,0:MD3,0:7)
      Dimension YDAS(0:MD1,0:MD2,0:MD3), YDAA(0:MD1,0:MD2,0:MD3)
      Common / Arrays / X1A,X2A,X3A, YAS,YAA,  X1Y,X2Y,X3Y, YAY
      Common / Dims / iD1,iD2,iD3, iY1,iY2,iY3
CsB___K-Factor stuff
      Dimension X1P(0:MD1), X2P(0:MD2), X3P(0:MD3)
      Dimension YKP(0:MD1,0:MD2,0:MD3)
      Dimension X1K(0:MD1), X2K(0:MD2), X3K(0:MD3)
      Dimension YKY(0:MD1,0:MD2,0:MD3)
      Integer iKFacP,iKFacY
      Common / KFactor / iKFacP,iKFacY,iYGrid
      Common / ArraysK / X1P,X2P,X3P,YKP, X1K,X2K,X3K,YKY
      Common / DimsK / iP1,iP2,iP3, iK1,iK2,iK3
      Real*8 aKP,aKY
CsB___
      Common /ECMC/ ECMC
      Common / Boson / Boson
      Integer N_SUD_AB,N_SUD_A,N_SUD_B,N_WIL_C,I_SUD
      Common / Sudakov / N_SUD_AB,N_SUD_A,N_SUD_B,N_WIL_C,I_SUD
!      CHARACTER*10 PDF_EVL_NAME
!      COMMON/PDFFILE/ PDF_EVL_NAME
      Common / Process / LTOpt
CsB NCMx = # of columns in the grid
      Parameter (NCMx = 20)
      Dimension A(NCMx)
      Common / NumCol / NC
      Common / YPiece / iYPiece, iReadIn
      REAL*8 QT_V_MIN, QT_V_MAX,YB_V_MAX,Q_V_MIN,Q_V_MAX,QT_V_MXY
      COMMON/HARDCUT/ QT_V_MIN, QT_V_MAX,YB_V_MAX,
     >  Q_V_MIN,Q_V_MAX,QT_V_MXY
      COMMON/AA/W,QCDL,X1X2MIN,NJET,N_PRIMARY_PART,NF,I_STR,IBEAM,NPART,
     &  NEVENT,NORDER,N_SECOND_PART
      Common / OneQ / iOneQ

      REAL*8 XMW,XMZ,XMW2,XMZ2,AEM,GAMW,GAMZ,GAMW2,GAMZ2
     &       SWS,CWS,XMBOT,XMTOP,XMB2,XMT2,GWEAK2
     &       XMCHARM,XMCHARM2,XMHIGGS,GAMH,GAMH2
      COMMON/STAND0/XMW,XMZ,XMW2,XMZ2,AEM,GAMW,GAMZ,GAMW2,GAMZ2
      COMMON/STAND1/SWS,CWS,XMBOT,XMTOP,XMB2,XMT2,GWEAK2
      COMMON/STAND5/XMCHARM,XMCHARM2,XMHIGGS,GAMH,GAMH2
      Integer iProc
      Common / PartonProcess / iProc
      Integer iResScheme, iLV
      Real*8 h1v,h1L
      Common/ ReSumScheme / iResScheme,iLV,h1v,h1L
      Integer lHead
      Integer BIPDF,iSet, iNonpert,iFlag_c3,KinCorr
      Common / BuildInPDF / BIPDF,iSet
      Real*8 C1,C2,C3,C4
      Common / LegacyCi / C1,C2,C3,C4
CCAO! The read W-boson mass and width
      Real*8 XMW_IN,XMW_GRID,GAMW_IN
      COMMON/WMASS_IN/XMW_IN,XMW_GRID,GAMW_IN
      Real*8 XMZ_IN,XMZ_GRID,GAMZ_IN
      COMMON/ZMASS_IN/XMZ_IN,XMZ_GRID,GAMZ_IN

CCPY
      LOGICAL Y_ONE_SIDE
      COMMON/GRID_FILE_FORM/ Y_ONE_SIDE

      Logical Testing
      Data Testing / .False. /

CJI NOV 2020
      INTEGER*4 KFAC_ANGFUNC
      COMMON/ANGFUNC/ KFAC_ANGFUNC

C      Testing = .True.

      If (Testing) Print*, ' iD1, iD2, iD3 = ', iD1, iD2, iD3
      If (Testing) Print*, ' iY1, iY2, iY3 = ', iY1, iY2, iY3

ccc      Print*, ' iD1, iD2, iD3 = ', iD1, iD2, iD3
ccc      Print*, ' iY1, iY2, iY3 = ', iY1, iY2, iY3

CsB___Find out the length of the header
      lHead = iHeadLen(2,Filename)
      Print*, ' The length of the header is ',lhead
      
CsB___Dump the header
      If (iReadIn.Eq.0) then ! W,Asy,Prt... pieces
CsB     Length of old header is 16 lines
        OPEN(UNIT=2,FILE=Filename,STATUS='OLD') !,ReadOnly)
C       Read in additional two or three lines for identifying the grid file.
        Do I = 1, 18
          Read(2,'(A)',Err=99,End=99) dummy
          Write(16,'(A)') dummy
        End Do
        Write(16,*)
     >'#############################################################'
        Close (2)

CsB_____Read the header. Length of new (old=FORTRAN) header is 15 (16) lines.
        OPEN(UNIT=2,FILE=Filename,STATUS='OLD') !,ReadOnly)
        Read(2,'(A)',Err=99,End=99) dummy
        If (lHead.Eq.16) then ! Grid from FORTRAN version of Legacy 
          Read(2,*) ECMC, iBeam, nDummy, LTOpt, iProc
          iResScheme = 1 
          iLV = 1  ! version <= 2.2
        Else                  ! Grid from C++ version of Legacy 
          Read(2,*,Err=87,End=87)             ! version > 2.2
     ,              ECMC, iBeam, LTOpt, iProc, iResScheme, 
     ,              nDummy, nDummy, nDummy, nDummy
!       Print*,     ECMC, iBeam, LTOpt, iProc, iResScheme, 
!    ,              nDummy, nDummy, nDummy, nDummy
   87     Print*, ' Legacy v2.2 or earlier header '
CCPY   88     Print*, ' Legacy v2.3 or later header '
          iResScheme = 1 
          iLV = 2  
        End If
	  Read (2,'(2X,A)',Err=99,End=99) dummy
          Read (2,*) iSet,iNonpert,iFlag_c3,C1,C2,C3,C4,KinCorr
!       Print*, iset,inonpert,iflag_c3,C1,C2,C3,C4,KinCorr
C_______Skip 1 line on old grids
        Read (2,'(2X,A)',Err=99,End=99) dummy
C_______Read  MT,MW,MZ,MH,MA
          Read (2,*) xMTop, xMW_Grid, xMZ_Grid 
          Read (2,*) xMHiggs,Dummy         
!       Print*, ' mt,mW,mZ,mH,mA =', xMTop,xMW_Grid,xMZ_Grid,xMHiggs,nDummy
C_______Read Boson type
        Read (2,'(12X,A)',Err=99,End=99) Boson
!       Print*, Boson
C_______Skip 2 lines
        Do I = 1, 2
          Read (2,'(2X,A)',Err=99,End=99) dummy
        EndDo
        If (iLv.eq.1) then
          h1L=0
          Read (2,'(2X,A)',Err=99,End=99) dummy
          Read (2,*) nDummy,nDummy,N_SUD_AB,N_WIL_C,I_FSR
        Else
          Read (2,*) rDummy,   h1L,rDummy
          Read (2,'(2X,A)',Err=99,End=99) dummy
          Read (2,*) nDummy,nDummy,N_SUD_A,N_SUD_B,N_WIL_C,I_SUD
        EndIf
!       Print*, N_SUD_A,N_SUD_B,N_WIL_C,I_SUD
        If (iLv.eq.1) then
          nnn = 3
        Else
          nnn = 1
        End If
        Do I = 1, nnn
          Read (2,'(2X,A)',Err=99,End=99) dummy
        End Do

        Print*, ' Input parameters from the grid header:'
        Print*, '  ECMC = ', ECMC
        Print*, '  iBeam,nDummy,LTOpt,iProc = ',iBeam,nDummy,LTOpt,iProc
        Print*, '  iResScheme,h1L = ', iResScheme,h1L
        Print*, '  iSet,iNonpert,iFlag_c3,KinCorr = ', 
     ,             iSet,iNonpert,iFlag_c3,KinCorr
        Print*, '  C1,C2,C3,C4 = ', C1,C2,C3,C4
        Print*, '  Boson type = ', Boson
        If (lHead.Eq.15) Print*, '  N_SUD_A,N_SUD_B,N_WIL_C,I_SUD = ', 
     ,                              N_SUD_A,N_SUD_B,N_WIL_C,I_SUD
        If (lHead.Eq.16) Print*, '  N_SUD_AB,N_WIL_C,I_SUD = ', 
     ,                              N_SUD_AB,N_WIL_C,I_SUD
      Else If (iReadIn.Eq.1) then ! Y piece
CCPY        If (iYGrid.Eq.0) then ! no Y grid (with O(alpha_s^2) K-factors) output
CsB_______Skip the header (16 lines)
C          OPEN(UNIT=2,FILE=Filename,STATUS='OLD') !,ReadOnly)
C          Do ii = 1,lHead
C            Read(2,'(A)',Err=99,End=99) dummy
C          End Do
CCPY        Else If (iYGrid.Eq.1) then ! Y grid (with O(alpha_s^2) K-factors output)
        If (iYGrid.Eq.1) then ! Y grid file with old format
            IF(KFAC_ANGFUNC.eq.2) THEN
              PRINT*, 'KFAC_ANGFUNC=2 requires Angular k-factor grid'
              Call Exit
            END IF 
C iYGrid=1
C    ' Q,qT,y, Singular (K*L0,K*A3), Pert. (KP*L0,KP*A3,A1,A2,A4)' 
C
CsB_______Dump the header (15 lines)
          OPEN( UNIT=2,FILE=Filename,STATUS='OLD') !,ReadOnly)
          OPEN(UNIT=15,FILE='Y_Grid.out') !, STATUS='NEW')
          IF(lHead.ne.15) THEN
            PRINT*, 'This is for iYGrid=1'
            PRINT*, 'Must have 15 comment lines in the Y-Grid file'
            Call Exit
          END IF 
          Do ii = 1,lHead
            Read(2,'(A)',Err=99,End=99) dummy
          End Do
          Write(15,'(A)') 'This is for iYGrid=1'
          Write(15,'(A)')
     >    '  Q,qT,y, Singular (K*L0,K*A3), Pert. (KP*L0,KP*A3,A1,A2,A4)'
        Else If (iYGrid.Eq.2) then ! Y grid file with new format
CJI Nov 2020
          PRINT*, 'Forcing KFAC_ANGFUNC=2 for this grid format'
          KFAC_ANGFUNC=2
C iYGrid=2
C    ' Q,qT,y, Singular (L0,A3), Pert. (L0,A3), R_A1, R_A2, R_A4, R_A0'
C
CsB_______Dump the header (17 lines)
          OPEN( UNIT=2,FILE=Filename,STATUS='OLD') !,ReadOnly)
          OPEN(UNIT=15,FILE='Y_Grid.out') !, STATUS='NEW')
          IF(lHead.ne.17) THEN
            PRINT*, 'This is for iYGrid=2'
            PRINT*, 'Must have 17 comment lines in the Y-Grid file'
            Call Exit
          END IF 
          Do ii = 1,lHead
            Read(2,'(A)',Err=99,End=99) dummy
          End Do
          Write(15,'(A)') 'This is for iYGrid=2'
          Write(15,'(A)')
     >    ' Q,qT,y, Singular (L0,A3), Pert. (L0,A3), R_A1, R_A2, R_A4, 
     > R_A0'
        End If ! iYgrid
      End If ! iReadIn

CsB___Set nPart by LTOpt (which is read from grid)
C     NPart -> Switches between LO, NLO and resummed processes.
C              11 does Resummed, 12 Asymptotic, 13 Delta-sigma, 14 LO.
      If (LTOpt.Eq.-1) then
        NPart = 14
      Else If (LTOpt.Eq.0 .or. LTOpt.Eq.4) then
        NPart = 11
      Else If (LTOpt.Eq.2) then
        NPart = 12
      Else If (LTOpt.Eq.1) then
        NPart = 13
      Else
        Print*, ' No such LTOpt: ', LTOpt
      End If

CsB___Read the data
      If (iReadIn.Eq.0) then !
        iDD1 = iD1
        iDD2 = iD2
        iDD3 = iD3
      Else
        iDD1 = iY1
        iDD2 = iY2
        iDD3 = iY3
      End if
      Print*, ' Reading content of ', Filename
      Q_V = -10.d0
      IQ  = -1
!      Print*, ' iDD3+1, iDD2-1, iDD1 ', iDD3+1, iDD2-1, iDD1
      Do 160 K = 1, iDD3+1
        Do 140 J = 1, iDD2-1
        Do 140 I = 1, iDD1
          Read(2,*,ERR=99,End=99) (A(ii), ii=1,NC)
          If (Q_V.NE.A(1)) then
            Q_V = A(1)
            IQ = IQ + 1
          EndIf
          If (iReadIn.Eq.0) then ! for other then the Y piece ===============
            X3A(IQ) = Q_V        ! Print*, Q_V, IQ
            X1A(I)  = A(2)       ! qT             
            X2A(J)  = A(3)       ! y              
            YAS(I,J,IQ) = A(4)   ! CSS L0         
            YAA(I,J,IQ) = A(5)   ! CSS A3         
            If (Testing) then
            If (A(4).lt.0.d0 .or. A(5).lt.0.d0) then
              Print*, X3A(IQ),X1A(I),X2A(J), YAS(I,J,IQ),YAA(I,J,IQ)
!              Stop
            End If
            End If
          Else                 ! for Y piece --------------------------------
            If (iYGrid.Eq.1) then ! Y grid file with old format
C    ' Q,qT,y, Singular (K*L0,K*A3), Pert. (KP*L0,KP*A3,A1,A2,A4)' 
C
CsB YAY has the Y piece and the perturbative pieces in it.
C   The last index of YAY is I_Pert.
C   I_Pert Assignments. I_Pert:XSect -> 0:L0, 1:A1, 2:(A2+A0), 3:2*A3, 4:A4.
C   YAY(I,J,IQ,5&6) has the perturbative L0&A3 pieces
            X3Y(IQ) = Q_V               ! Print*, Q_V, IQ
            X1Y(I)  = A(2)              ! qT                        
            X2Y(J)  = A(3)              ! y                         
C Y pieces
            YAY(I,J,IQ,0) = A(6) - A(4) ! regular L0 piece
            YAY(I,J,IQ,1) = A(8)        !         A1 piece
            YAY(I,J,IQ,2) = A(9)        !         A2 = A0 piece
            YAY(I,J,IQ,3) = A(7) - A(5) ! regular A3 piece
            YAY(I,J,IQ,4) = A(10)       !         A4 piece
C Perturbative pieces
            YAY(I,J,IQ,5) = A(6)        ! perturbative L0 piece
            YAY(I,J,IQ,6) = A(7)        !              A3 piece
            YAY(I,J,IQ,7) = 0.D0        ! dummy 
!            Print*, X3Y(IQ), X1Y(I), X2Y(J)
!            Print*, YAY(I,J,IQ,5), YAY(I,J,IQ,6)
!            Print*, YAY(I,J,IQ,1), YAY(I,J,IQ,2), YAY(I,J,IQ,4)
!            Stop
            Else If (iYGrid.Eq.2) then ! Y grid fiel with new format
C iYGrid=2
C    ' Q,qT,y, Singular (L0,A3), Pert. (L0,A3), R_A1, R_A2, R_A4, R_A0'
C
C   I_Pert Assignments. I_Pert:XSect -> 0:L0, 1:A1, 2:(A2+A0), 3:2*A3, 4:A4.
C   YAY(I,J,IQ,5&6) has the perturbative L0&A3 pieces
            X3Y(IQ) = Q_V               ! Print*, Q_V, IQ
            X1Y(I)  = A(2)              ! qT                        
            X2Y(J)  = A(3)              ! y                         
C Y pieces
            YAY(I,J,IQ,0) = A(6) - A(4) ! regular L0 piece
            YAY(I,J,IQ,1) = A(8)        !         R_A1 piece
            YAY(I,J,IQ,2) = A(9)        !         R_A2 piece
            YAY(I,J,IQ,3) = A(7) - A(5) ! regular A3 piece
            YAY(I,J,IQ,4) = A(10)       !         R_A4 piece
            YAY(I,J,IQ,7) = A(11)       !         R_A0 piece
C Perturbative pieces
            YAY(I,J,IQ,5) = A(6)        ! perturbative L0 piece
            YAY(I,J,IQ,6) = A(7)        !              A3 piece
C            Print*, X3Y(IQ), X1Y(I), X2Y(J)
C            Print*, YAY(I,J,IQ,5), YAY(I,J,IQ,6)
C            Print*, YAY(I,J,IQ,1), YAY(I,J,IQ,4)
C            Print*, YAY(I,J,IQ,2), YAY(I,J,IQ,7)
C            Stop

            End If ! iYGrid

CsB_________K-Factor Stuff
            If (iKfacP.Eq.1 .or. iKfacY.Eq.1 ) then ! = = = = = = = = = = = =

CCPY!              If (Boson.Eq.'H0') then ! :::::::::::::::::::::::::::::::::::::
CCPY!CsB_____________Pert. K-factor from Florian,Grazzini,Kunszt ph/9902483
CCPY!                aKP = 1.5 ! for CTEQ4M
CCPY!                YAY(I,J,IQ,5) = aKP*YAY(I,J,IQ,5)
CCPY!              Else ! ........................................................

cCsB_________This is to avoid to hit the edge of the K-factor grids:
c              If (X3Y(IQ).Eq.X3P(0)) then
c                QInt = 1.0001d0*X3Y(IQ)
cC                Print*, ' X3Y(IQ), X3P(0): ', X3Y(IQ), X3P(0)
c              Else If (X3Y(IQ).Eq.X3P(iP3)) then
c                QInt = 0.9999d0*X3Y(IQ)
cC                Print*, ' X3Y(IQ), X3P(iP3): ', X3Y(IQ), X3P(iP3)
c              Else
                QInt = X3Y(IQ)
c              End If
cC              If (X1Y(I).Eq.X1P(0)) then
cC                Print*, ' X1Y(I), X1P(0): ', X1Y(I), X1P(0)
C                QTInt = 1.0001d0*X1Y(I)
CsB___________Note: As of 6/9/97
C             The Y piece K-Factor input makes not much sense for 
C             QT < 10-15 GeV
              If (X1Y(I).LE.10) then
                QTInt = 10.d0
c              Else If (X1Y(I).Eq.X1P(iP1)) then
c                QTInt = 0.9999d0*X1Y(I)
              Else
                QTInt = X1Y(I)
              End If
c              If (X2Y(J).Eq.X2P(0)) then
c                yInt = 1.0001d0*X2Y(J)
c              Else If (X2Y(J).Eq.X2P(iP2)) then
c                yInt = 0.9999d0*X2Y(J)
c              Else
                yInt = X2Y(J)
c              End If
              Call ThreeKInterp ( 4,4,4, YDAS,YDAA
     >        ,QTInt,yInt,QInt, aKP,aKY,aErrP,aErrY )
C     >        ,YKP,YKY, QTInt,yInt,QInt, aKP,aKY,aErrP,aErrY )
C     >        ,YKP,YKY, X1Y(I),X2Y(J),X3Y(IQ), aKP,aKY,aErrP,aErrY )
              If (Testing) Print*, QInt,yInt,QTInt, aKP,aKY
CsB___________Apply K-factors
              ! for L0 pieces
              YAY(I,J,IQ,0) = aKY*YAY(I,J,IQ,0)
              YAY(I,J,IQ,5) = aKP*YAY(I,J,IQ,5)
              ! for A3 piece - added on 10/22/98
              YAY(I,J,IQ,3) = aKY*YAY(I,J,IQ,3)
              YAY(I,J,IQ,6) = aKP*YAY(I,J,IQ,6)
              YAY(I,J,IQ,7) = 0.D0   ! dummy

CCPY Oct 2008
CCPY JAN 8, 2007 APPLY aKP TO ALL OTHER ANGULAR FUNCTIONS
              YAY(I,J,IQ,1) = aKP*YAY(I,J,IQ,1)
              YAY(I,J,IQ,2) = aKP*YAY(I,J,IQ,2)
              YAY(I,J,IQ,4) = aKP*YAY(I,J,IQ,4)
C
CsB___________Write out new Y-grid
C               -> Q, Q_T, y:
                WRITE(15,103) X3Y(IQ), X1Y(I), X2Y(J)
C               -> Asymptotic L0, A3 and perturbative L0, A3:
CsB             Note: Asy2 = Prt2 - Y2 = KP*Prt1 - KY*(Prt1 - Asy1)
C                PRINT 107,   aKP*A(6)-aKY*(A(6)-A(4)),
C     ,                       aKP*A(7)-aKY*(A(7)-A(5)),aKP*A(6),aKP*A(7)
                WRITE(15,107) aKP*A(6)-aKY*(A(6)-A(4)),
     ,                        aKP*A(7)-aKY*(A(7)-A(5)),aKP*A(6),aKP*A(7)
C               -> Perturbative A1, A2, A4:
C                PRINT 105,   A(8),A(9),A(10)
                WRITE(15,105) A(8),A(9),A(10)
  103           Format(1X,3(G10.4,2X))
  105           Format(1X,3(G18.12,2X))
  107           Format(1X,4(G18.12,2X))

            End If ! iKfacP.Eq.1 .or. iKfacY.Eq.1 = = = = = = = = = = = = = =
CsB_________
          End If ! ==========================================================


        Y_ONE_SIDE=.FALSE.
CCPY Sept 2009: In case that one-side y-grid file is used, we set the flag 
C to keep track on it. Here we assume the one-side y-grid starts with 
C y=0.0, and positive y-grids are given.
C 
      If (iReadIn.Eq.0) then ! W,Asy,Prt... pieces
C The first y value in the y-grid file
        IF(DABS(X2A(1)).LT.1.0D-5) Y_ONE_SIDE=.TRUE.
      Else If (iReadIn.Eq.1) then ! Y piece
        IF(DABS(X2Y(1)).LT.1.0D-5) Y_ONE_SIDE=.TRUE.
      End If

      IF(Y_ONE_SIDE) THEN
        If (Boson.NE.'A0' 
     ..AND.Boson.NE.'AA' .AND. Boson.NE.'AG'
     ..AND.Boson.NE.'ZZ' .AND. Boson.NE.'ZG'
     ..AND.Boson.NE.'H0'
     ..AND.iBeam.NE.1) then
          PRINT*,'Can y_one_side grid file be used for this process?'
          call ExitGracefully
        ENDIF
      EndIf

      If (Testing) then
C          Print*, IQ, I, J, X3A(IQ), X1A(I), X2A(J), YAS(I,J,IQ), YAA(I,J,IQ)
        If (IQ.Eq.3) then
          If (iReadIn.Eq.1) then
C            Print*, X3A(IQ), X1A(I), X2A(J), YAS(I,J,IQ), YAA(I,J,IQ)
            Print*, X3A(IQ),X1A(I),X2A(J), YAY(I,J,IQ,1),YAY(I,J,IQ,2)
          End If
        End If
      End If
140   Continue

CCPY March 2010
      IF(BOSON.EQ.'AA'.OR.BOSON.EQ.'AG'.OR.
     >  BOSON.EQ.'ZZ'.OR.
     >  BOSON.EQ.'WW_UUB'.OR.BOSON.EQ.'WW_DDB') THEN
C DO NOTHING
      ELSE
CsB * Xsect(qT=0) = 0 for both resummed and regular (Y) pieces.
C     For perturbative pieces it does not really matter here.
        X1A(0) = 0.d0
        X1Y(0) = 0.d0
        Do jj = 0, iDD2
          If (iReadIn.Eq.0) then
            YAS(0,jj,iQ) = 0.0d0
            YAA(0,jj,iQ) = 0.0d0
          Else
CCPY Sept 2009
            Do ii = 0,7
              YAY(0,jj,iQ,ii) = 0.0d0
CCPY              YAY(0,jj,iQ,ii) = YAY(1,jj,iQ,ii)
            End Do
            Do ii = 5,6
              YAY(0,jj,iQ,ii) = YAY(1,jj,iQ,ii)
            End Do
          End If
        End Do
      ENDIF

CCPY!CsB * Xsect(y=ymin,ymax) = 0 for both resummed and Y pieces
CCPY!        Do 160 jj = 0, iDD1
CCPY!          If (iReadIn.Eq.0) then
CCPY!            YAS(jj,0,IQ) = 0.0d0
CCPY!            YAA(jj,0,IQ) = 0.0d0
CCPY!            YAS(jj,iD2,iQ) = 0.0d0
CCPY!            YAA(jj,iD2,iQ) = 0.0d0
CCPY!          Else
CCPY!            Do ii = 0,7
CCPY!              YAY(jj,  0,iQ,ii) = 0.d0
CCPY!              YAY(jj,iY2,iQ,ii) = 0.d0
CCPY!            End Do
CCPY!            Do ii = 5,6
CCPY!              YAY(jj,  0,iQ,ii) = 0.d0
CCPY!              YAY(jj,iY2,iQ,ii) = 0.d0
CCPY!            End Do
CCPY!          End If
CCPY!

CCPY Sept 2009 Deal with Xsect(y=ymin,ymax) for both resummed and Y pieces
CCPY: Assuming a large enough y range is covered in the y-grid file.
        Do 160 jj = 0, iDD1
          If (iReadIn.Eq.0) then
            YAS(jj,0,IQ) = YAS(jj,1,IQ)
            YAA(jj,0,IQ) = YAA(jj,1,IQ)
            YAS(jj,iD2,iQ) = YAS(jj,iD2-1,iQ)
            YAA(jj,iD2,iQ) = YAA(jj,iD2-1,iQ)
          Else
            Do ii = 0,7
              YAY(jj,  0,iQ,ii) = YAY(jj,  1,iQ,ii)
              YAY(jj,iY2,iQ,ii) = YAY(jj,iY2-1,iQ,ii)
            End Do
            Do ii = 5,6
              YAY(jj,  0,iQ,ii) = YAY(jj,  1,iQ,ii)
              YAY(jj,iY2,iQ,ii) = YAY(jj,iY2-1,iQ,ii)
            End Do
          End If


CsB   Except for photon, or p p collider where:
C     (p N collision is an exception)
          If (Boson.Eq.'AA' .and. iBeam.Eq.0) Goto 160

          IF(Y_ONE_SIDE) THEN

            If (iReadIn.Eq.0) then
C             Assuming that grid is generated only for positive rapidites.
              YAS(jj,0,IQ) = YAS(jj,2,IQ)
              YAA(jj,0,IQ) = YAA(jj,2,IQ)
C              Assuming that grid is generated only for negative rapidites.
c              YAS(jj,iD2,iQ) = YAS(jj,iD2-2,iQ)
c              YAA(jj,iD2,iQ) = YAA(jj,iD2-2,iQ)
            Else
              Do ii = 0,7
                YAY(jj,0,iQ,ii) = YAY(jj,2,iQ,ii)
              End Do
              Do ii = 5,6
                YAY(jj,0,iQ,ii) = YAY(jj,2,iQ,ii)
              End Do
            End If
          End If
160   Continue

      Close(2)
      If ((iYGrid.Eq.1).OR.(iYGrid.Eq.2)) Close(15)

CCPY!CsB * The values of yMax, yMin are dummies at this moment.
CCPY!CsB * The actual values will be determined before the interpolation.
CCPY!      If (iReadIn.Eq.0) then
CCPY!        yStep = ( X2A(iDD2-1)-X2A(1) )/(iDD2 - 2)
CCPY!        X2A(  0) = X2A(1)     - yStep
CCPY!        X2A(iD2) = X2A(iD2-1) + yStep
CCPY!      Else
CCPY!        yStep = ( X2Y(iDD2-1)-X2Y(1) )/(iDD2 - 2)
CCPY!        X2Y(  0) = X2Y(1)     - yStep
CCPY!        X2Y(iY2) = X2Y(iY2-1) + yStep
CCPY!      End If
CCPY!

CCPY Sept 2009: Assuming the y-grid is arranged so that the value 
C of y increases. This setup is also suitable for the case with 
C only one-side (in y) grid file.
      If (iReadIn.Eq.0) then
        X2A(  0) = X2A(1)- DABS( X2A(2)-X2A(1) )
        X2A(iD2) = X2A(iD2-1) + DABS( X2A(iD2-1)-X2A(iD2-2) )
      Else
        X2Y(  0) = X2Y(1)- DABS( X2Y(2)-X2Y(1) )
        X2Y(iY2) = X2Y(iY2-1) + DABS( X2Y(iY2-1)-X2Y(iY2-2) )
      End If

cdump
c      print*,' idd2,id2,iy2 =',idd2,id2,iy2
c      Print*, ' yStep = ', yStep
c      Print*, X2A
c      Print*, '...'
c      Print*, X2Y

CsB For the special case when there's only one Q in the grid.
      If (iDD3.Eq.0) then
        iOneQ = 1
C   Pretend that there are 4 Q's, with small separation
        If (iReadIn.Eq.0) then
          iD3 = iD3 + 3
        Else If (iReadIn.Eq.1) then
          iY3 = iY3 + 3
        End If
        iDD3 = iDD3+3
C        Print*, ' iD3 is increased by ', 3
        Print*, ' +---------------------------------+ '
        Print*, ' ! The V width is taken to be ZERO ! '
        Print*, ' +---------------------------------+ '
        ep  = .001d0
        If (iReadIn.Eq.0) then
          Do ii = 1,3
            X3A(ii) = X3A(0) + ii*ep
          End Do
        Else If (iReadIn.Eq.1) then
          Do ii = 1,3
            X3Y(ii) = X3Y(0) + ii*ep
          End Do
        End If

CsB Fill up the ghost arrays, so that they contain the same values as the
C   single Q array.
        If (iReadIn.Eq.0) then
          Do I = 0,iDD1+1
          Do J = 0,iDD2+1
          Do K = 1,iDD3+1
            YAS(I,J,K) = YAS(I,J,0)
            YAA(I,J,K) = YAA(I,J,0)
          End Do
          End Do
          End Do
        Else If (iReadIn.Eq.1) then
          Do I = 0,iDD1+1
          Do J = 0,iDD2+1
          Do K = 1,iDD3+1
            Do ii = 0,7
              YAY(I,J,K,ii) = YAY(I,J,0,ii)
            End Do
          End Do
          End Do
          End Do
        End If
      End If

CsB Determine QT_V_MIN, QT_V_MAX,YB_V_MAX,Q_V_MIN,Q_V_MAX
C   from the resummed grid.
C     X1A (X2A, X3A) contains the qT (y, Q) vector
C     iD1+1 (iD2+1, iD3+1) is the number of qT (y, Q) values in the grid
      QT_V_MIN = X1A(0)
      QT_V_MAX = X1A(iD1)
      YB_V_MAX = X2A(iD2) - yStep
      Q_V_MIN  = X3A(0)
      Q_V_MAX  = X3A(iD3)
CsB Determine QT_V_MXY
C   from the Y grid.
      QT_V_MXY = X1Y(iDD1)
C      Print*, ' Q_V_MIN,Q_V_MAX, QT_V_MIN,QT_V_MIN, QT_V_MAX, YB_V_MAX, QT_V_MXY:'
C      Print*,   Q_V_MIN,Q_V_MAX, QT_V_MIN,QT_V_MIN, QT_V_MAX, YB_V_MAX, QT_V_MXY

      Return

 99   Print*, ' Error at reading data in ReadIn. '
      Print*, ' i = ', i, ', j = ', j, ', k = ', k
      Print*, ' X1A(i) = ',X1A(i),', X2A(j) =',X2A(j),', Q_V = ',Q_V
      Stop
      End ! ReadIn

CsB * *******************************************************************
C
      Subroutine PreReadIn(Filename, Mode)
C
C_____Mode: 1 = main grid
C           2 = Y piece
C           3 = K-factor for Pert piece
C           4 = K-factor for Y piece
C
      Implicit Double Precision (A-H,O-Z)
      Character*145 Filename
      Common / Dims  / iD1,iD2,iD3, iY1,iY2,iY3
      Common / DimsK / iP1,iP2,iP3, iK1,iK2,iK3
      Logical Testing
      Common / NumCol / NC
c      Common / YPiece / iYPiece, iReadIn

C      Testing = .True.
      Testing = .False.

CsB___Find out the length of the header
      lHead = iHeadLen(2,Filename)

CsB   Sample the file
CsB   Determine the number of data columns in the file
      Call Column (2, lHead, Filename, NC)
C      If (Testing) Print*, 'NC = ', NC
CsB   Determine the length of the file,
C     and the number of Q points in the grid.
C     The maximum length of file is hard-wired here to be 500000 lines.
      Length = 10000000
      Call Sample (2, lHead, Filename, 1, Length, NC, NQ)
C      If (Testing) Print*, ' Length, NQ = ', Length, NQ
CsB   Determine the number of y and qT points in the grid
      Length = Length/NQ
      Call Sample (2, lHead, Filename, 3, Length, NC, Ny)
C      If (Testing) Print*, ' Length, Ny = ', Length, Ny
      If (Mode.Eq.1) then
        iD1 = Length/Ny
        iD2 = Ny + 1
        iD3 = NQ - 1
      Else If (Mode.Eq.2) then
        iY1 = Length/Ny
        iY2 = Ny + 1
        iY3 = NQ - 1
      Else If (Mode.Eq.3) then
        iP1 = Length/Ny
        iP2 = Ny + 1
        iP3 = NQ - 1
      Else If (Mode.Eq.4) then
        iK1 = Length/Ny
        iK2 = Ny + 1
        iK3 = NQ - 1
      End If
      If (Testing) Print*, ' iD1, iD2, iD3 = ', iD1,iD2,iD3
      If (Testing) Print*, ' iY1, iY2, iY3 = ', iY1,iY2,iY3
      If (Testing) Print*, ' iP1, iP2, iP3 = ', iP1,iP2,iP3
      If (Testing) Print*, ' iK1, iK2, iK3 = ', iK1,iK2,iK3

      Return
      End

CPN * *******************************************************************

      function iHeadLen(iUnit,FileName)
CPN Returns the length of the comment header in the file by looking for the
C   line marking the end of the block. 
CsB > This is not right anymore: The first 5 symbols in this line must
C   be #####.

      implicit NONE

      integer iHead,iUnit,iHeadLen
      character*80 dummy1
      character*8  dummy2
      equivalence (dummy1,dummy2)

      character*145 FileName

      iHead=0
      dummy1=' '

      OPEN(UNIT=iUnit,FILE=Filename,STATUS='OLD')

      do while (.not.(dummy2.eq.'  Q,qT,y'))
       read (iUnit,'(A)') dummy1
!       Print*, dummy2
       iHead=iHead+1
      enddo

      close(iUnit)

      iHeadLen=iHead
      End! iHeadLen


CsB * *******************************************************************

      Subroutine Sample (iUnit,iHead,Filename,iComp,Length,NC, Number)

      Implicit Double Precision (A-H,O-Z)
      Character*145 Filename
      Logical Testing
      Parameter (NCMx = 20)
      Dimension A(NCMx)

      Testing = .True.
      Testing = .False.

      If (Testing) Print*, 'in Sample'

      OPEN(UNIT=iUnit,FILE=Filename,STATUS='OLD') !,ReadOnly)
CsB   Skip the header
      Do i = 1, iHead
        Read(iUnit,'(A)',Err=99,End=99) dummy
      EndDo
CsB   Count different A(iComp)'s until Length/EOF
      D = -100.0
      Number = 0
      Do i = 1,Length
        Read(iUnit,*,ERR=98,End=98) (A(j), j=1,NC)
        Q = A(1)
        If (A(iComp).NE.D) then
           D = A(iComp)
           Number = Number + 1
        EndIf
      EndDo
CsB   When counting y values (double) check the line after the "last" line
C     to avoid coincidental agreements in checksum
      If (iComp.EQ.3) then
        Read(iUnit,*,ERR=98,End=98) (A(j), j=1,NC)
        If (A(1).Eq.Q) then
          Print*, ' File checksum BAD. The grid file is corrupt ! '
          Stop
        End If
      End If
      Close(iUnit)
      Return

  98  Continue
      Length = i - 1
      Close(iUnit)
      If (Testing) Print*, Length, Number
      Return

  99  Print*, ' Error at reading the header of the data. '
      Stop
      End

CsB * *******************************************************************

      Subroutine Column (iUnit, iHead, Filename, NC)
      Implicit Double Precision (A-H,O-Z)
      Character*145 Filename
      Logical Testing
      Parameter (NCMx = 20)
      Dimension A(NCMx)
      Common / YPiece / iYPiece, iReadIn

      Testing = .True.
      Testing = .False.

      If (Testing) Print*, ' in Column '
C      If (iReadIn.Eq.0) Print*, ' Opening low  qT grid... '
C      If (iReadIn.Eq.1) Print*, ' Opening Y piece grid... '
      OPEN(UNIT=iUnit,FILE=Filename,STATUS='OLD') !,ReadOnly)
CsB   Skip the header
      Do i = 1, iHead
c        Read(2,'(A)',Err=99,End=99) dummy
        Read(iUnit,'(A)',Err=99,End=99) dummy
      EndDo
CsB   Count different entries
c      Read(2,*,ERR=98,End=98) (A(i), i = 1,NCMx)
      Read(iUnit,*,ERR=98,End=98) (A(i), i = 1,NCMx)
      Do i = 2,NCMx
        If (A(i) .eq. A(1)) then
          NC = i - 1
          Close(iUnit)
          Return
        End if
      End Do
CsB   If same entry was not found in NCMx steps.
      Print*, ' Error finding same Q. '
      Stop

  99  Print*, ' Error at reading the header of datafile. '
  98  Print*, ' Error at reading datafile. '
      Stop
      End

CsB * *******************************************************************

      Subroutine ReadInK(FileName)
C
C     This routine reads K factors for Pert(1) and Y(1) pieces.
C
C     The following notation apply for the perturbative K factor:
C     X1P (X2P, X3P) contains the qT (y, Q) vector
C     iP1+1 (iP2+1, iP3+1) is the number of qT (y, Q) values in the grid
C     iPP1 (iPP2, iPP3) is the number of qT (y, Q)
C          values used at the interpolation
C     The following notation apply for the Y-piece factor:
C     X1K (X2K, X3K) contains the qT (y, Q) vector
C     iK1+1 (iK2+1, iK3+1) is the number of qT (y, Q) values in the grid
C     iKK1 (iKK2, iKK3) is the number of qT (y, Q)
C          values used at the interpolation
C     The 1st index of YKP and YKY belong to qT, the 2nd to y, the 3rd to Q.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      Real*8 ECMC
      Character*78 Dummy
      Character*40 Boson
      Character*145 FileName
CsB   The dimensions below give the maximal possible gridsize.
      Parameter (MD1 = 200, MD2 = 300, MD3 = 100)
      Dimension X1P(0:MD1), X2P(0:MD2), X3P(0:MD3)
      Dimension YKP(0:MD1,0:MD2,0:MD3)
      Dimension X1K(0:MD1), X2K(0:MD2), X3K(0:MD3)
      Dimension YKY(0:MD1,0:MD2,0:MD3)
      Common / ArraysK / X1P,X2P,X3P,YKP, X1K,X2K,X3K,YKY
      Common / DimsK / iP1,iP2,iP3, iK1,iK2,iK3
C      Common / Dims / iD1,iD2,iD3, iY1,iY2,iY3
C      Common /ECMC/ ECMC
C      Common / Boson / Boson
      Integer N_SUD_AB,N_SUD_A,N_SUD_B,N_WIL_C,I_SUD
C      Common / Sudakov / N_SUD_AB,N_SUD_A,N_SUD_B,N_WIL_C,I_SUD
      CHARACTER*20 PDF_EVL_NAME
C      COMMON/PDFFILE/ PDF_EVL_NAME
CsB NCMx = Max # of columns in the grid
      Parameter (NCMx = 20)
      Dimension A(NCMx)
      Common / NumCol / NC
c      Common / YPiece / iYPiece, iReadIn
      Common / cReadInK / iReadInK
C      REAL*8 QT_V_MIN, QT_V_MAX,YB_V_MAX,Q_V_MIN,Q_V_MAX,QT_V_MXY
C      COMMON/HARDCUT/ QT_V_MIN, QT_V_MAX,YB_V_MAX,Q_V_MIN,Q_V_MAX,QT_V_MXY
C      COMMON/AA/W,QCDL,X1X2MIN,NJET,N_PRIMARY_PART,NF,I_STR,IBEAM,NPART,
C     &  NEVENT,NORDER,N_SECOND_PART
      Common / OneQ / iOneQ

CCPY
      LOGICAL Y_ONE_SIDE
      COMMON/GRID_FILE_FORM/ Y_ONE_SIDE

C      REAL*8 XMW,XMZ,XMW2,XMZ2,AEM,GAMW,GAMZ,GAMW2,GAMZ2
C      COMMON/STAND0/XMW,XMZ,XMW2,XMZ2,AEM,GAMW,GAMZ,GAMW2,GAMZ2

      Logical Testing
      Data Testing / .False. /

C      Testing = .True.

      If (iReadInK.Eq.0) then
        iD1 = iP1
        iD2 = iP2
        iD3 = iP3
        iDD1 = iP1
        iDD2 = iP2
        iDD3 = iP3
      Else
        iY1 = iK1
        iY2 = iK2
        iY3 = iK3
        iDD1 = iK1
        iDD2 = iK2
        iDD3 = iK3
      End If
      If (Testing) Print*, ' iP1, iP2, iP3 = ', iP1,iP2,iP3
      If (Testing) Print*, ' iK1, iK2, iK3 = ', iK1,iK2,iK3

CsB___Read the header
CsB   Length of header is 16 lines
      OPEN(UNIT=2,FILE=FileName,STATUS='OLD') !,ReadOnly)
      Read(2,'(A)',Err=99,End=99) Dummy
      Read(2,*) ECMC, iBeam, nDummy, LTOpt
CCPY
      Read (2,'(2X,A)',Err=99,End=99) Dummy
CCPY      READ(2,'(5X,A)') PDF_EVL_NAME
      READ(2,'(A)') PDF_EVL_NAME
      Print*, ' K-factor calculated for PDF: ', PDF_EVL_NAME
C_____Skip 2 lines
      Do I = 1, 2
        Read (2,'(2X,A)',Err=99,End=99) Dummy
      EndDo
C_____Read  MT,MW,MZ,MH,MA
      Read (2,*) xMTop, XMW, XMZ
      Read (2,'(2X,A)',Err=99,End=99) dummy
C_____Read Boson type
        Read (2,'(12X,A)',Err=99,End=99) Boson
C_____Skip 4 lines
      Do I = 1, 4
        Read (2,'(2X,A)',Err=99,End=99) dummy
      EndDo
      Read (2,*) NDUMMY,NDUMMY,N_SUD_AB,N_WIL_C,I_SUD
      Do I = 1, 2
        Read (2,'(2X,A)',Err=99,End=99) dummy
      EndDo
      If (Testing) Print*, ' ECMC = ', ECMC, '  Boson : ', Boson

cCsB___Skip the header
cCsB   Length of header is 16 lines
c      OPEN(UNIT=2,FILE=Filename,STATUS='OLD') !,ReadOnly)
c      Do ii = 1,16
c        Read(2,'(A)',Err=99,End=99) dummy
c      End Do

CsB___Read the data
      Print*, ' Reading content of ', Filename
      Q_V = -10.d0
      IQ  = -1
c      Init = 0
      Do 160 K = 1, iDD3+1
        Do 140 J = 1, iDD2-1
        Do 140 I = 1, iDD1
          Read(2,*,ERR=99,End=99) (A(ii), ii=1,NC)
          If (Q_V.NE.A(1)) then
            Q_V = A(1)
            IQ = IQ + 1
          EndIf
          If (iReadInK.Eq.0) then ! for other then the Y piece
            X3P(IQ) = Q_V
C            Print*, Q_V, IQ
            X1P(I)  = A(2)
            X2P(J)  = A(3)
            YKP(I,J,IQ) = A(6)
            If (Testing) then
              Print*, X3P(IQ),X1P(I),X2P(J), YKP(I,J,IQ)
            End If
          Else                   ! for Y piece
            X3K(IQ) = Q_V
            X1K(I)  = A(2)
            X2K(J)  = A(3)
            YKY(I,J,IQ) = A(6)
            If (Testing) then
C              Print*, X3K(IQ),X1K(I),X2K(J), YKY(I,J,IQ)
            End If
          End If
      If (Testing) then
C          Print*, IQ,I,J, X3K(IQ),X1K(I), X2A(J), YAS(I,J,IQ), YAA(I,J,IQ)
        If (IQ.Eq.3 .and. J.Eq.4) then
          If (iReadInK.Eq.0) then
            Print*, X3P(IQ),X1P(I),X2P(J), YKP(I,J,IQ)
          Else If (iReadInK.Eq.1) then
            Print*, X3K(IQ),X1K(I),X2K(J), YKY(I,J,IQ)
          End If
        End If
      End If
140   Continue

CsB * K-Factor(qT=0) ~ K-Factor(qT=qT[first grid point])
C     for both Pert and Y pieces.
        X1P(0) = 0.d0
        X1K(0) = 0.d0
        Do jj = 0, iDD2
          If (iReadInK.Eq.0) then
            YKP(0,jj,iQ) = YKP(1,jj,iQ)
          Else
            YKY(0,jj,iQ) = YKY(1,jj,iQ)
          End If
        End Do

ccpy Oct 2008
CCPY ERROR: OCT 2007, THE FOLLOWING IS NOT CORRECT. 
CsB * K-Factor(y=ymin,ymax) = 1.0 for both Pert and Y pieces
C     (Both O(alpha_s) and O(alpha_s^2) cross sections are zero there,
C      so the value of K should not matter.)
C        Do 160 jj = 0, iDD1
C          If (iReadInK.Eq.0) then
C            YKP(jj,0,IQ)   = 1.0d0
C            YKP(jj,iP2,iQ) = 1.0d0
C          Else
C            YKY(jj,  0,iQ) = 1.d0
C            YKY(jj,iK2,iQ) = 1.d0
C          End If

CsB   Except for photon, or p p collider where:
C     (p N collision is an exception)
          If (Boson.Eq.'AA' .and. iBeam.Eq.0) Goto 160

          IF(Y_ONE_SIDE) THEN

            If (iReadInK.Eq.0) then
C             Assuming that grid is generated only for positive rapidites.
              YKP(jj,0,IQ) = YKP(jj,2,IQ)
C              Assuming that grid is generated only for negative rapidites.
c              YAS(jj,iD2,iQ) = YAS(jj,iD2-2,iQ)
c              YAA(jj,iD2,iQ) = YAA(jj,iD2-2,iQ)
            Else
              YKY(jj,0,iQ) = YKY(jj,2,iQ)
            End If
          End If
160   Continue

      Close(2)

CsB___The values of yMax, yMin are dummies at this moment.
CsB___The actual values will be determined before the interpolation.
      If (iReadInK.Eq.0) then
        yStep = ( X2P(iDD2-1)-X2P(1) )/(iDD2 - 2)
        X2P(  0) = X2P(1)     - yStep
        X2P(iD2) = X2P(iD2-1) + yStep
C        Print*, ' yStepP = ', yStep
      Else
        yStep = ( X2K(iDD2-1)-X2K(1) )/(iDD2 - 2)
        X2K(  0) = X2K(1)     - yStep
        X2K(iY2) = X2K(iY2-1) + yStep
C        Print*, ' yStepK = ', yStep
      End If

CsB For the special case when there's only one Q in the grid.
      If (iDD3.Eq.0) then
        iOneQ = 1
C   Pretend that there are 4 Q's, with small separation
        If (iReadInK.Eq.0) then
          iP3 = iP3 + 3
        Else If (iReadInK.Eq.1) then
          iK3 = iK3 + 3
        End If
        iDD3 = iDD3+3
C        Print*, ' iD3 is increased by ', 3
        Print*, ' +----------------------------------------+ '
        Print*, ' ! No V width when calculating K-Factor   ! '
        Print*, ' +----------------------------------------+ '
        ep  = .001d0
        If (iReadInK.Eq.0) then
          Do ii = 1,3
            X3P(ii) = X3P(0) + ii*ep
          End Do
        Else If (iReadInK.Eq.1) then
          Do ii = 1,3
            X3K(ii) = X3K(0) + ii*ep
          End Do
        End If

CsB Fill up the ghost arrays, so that they contain the same values as the
C   single Q array.
        If (iReadInK.Eq.0) then
          Do I = 0,iDD1+1
          Do J = 0,iDD2+1
          Do K = 1,iDD3+1
            YKP(I,J,K) = YKP(I,J,0)
          End Do
          End Do
          End Do
        Else If (iReadInK.Eq.1) then
          Do I = 0,iDD1+1
          Do J = 0,iDD2+1
          Do K = 1,iDD3+1
            YKY(I,J,K) = YKY(I,J,0)
          End Do
          End Do
          End Do
        End If
      End If

      Return

 99   Print*, ' Error at reading data in ReadInK. '
      Print*, ' i = ', i, ', j = ', j, ', k = ', k
      Stop
      End ! ReadInK

CsB -------------------------------------------------------------------------
      Subroutine CheckSum (Filename)
CsB -------------------------------------------------------------------------

      Implicit Double Precision (A-H,O-Z)
      Character*145 Filename

      Print*, ' Checking consistency of ', Filename

CsB___Find out the length of the header
!      Print*, ' Finding length of header of ', Filename
      lHead = iHeadLen(2,Filename)
!      Print*, ' Length of header  = ', lHead

CsB___Determine the number of data columns in the file
      Call Column (2, lHead, Filename, NC)
!      Print*, ' # of columns = ', NC

CsB___Determine the length of the file, and the number of Q points in the grid.
C     The maximum length of file is hard-wired here to be 500000 lines.
      Length = 10000000
      Call Sample (2, lHead, Filename, 1, Length, NC, NQ)
      LengthFile = Length
!      Print*, ' # of  Q points = ', NQ

CsB___Determine the number of y points in the grid
      Length = Length/NQ
      Call Sample (2, lHead, Filename, 3, Length, NC, Ny)
!      Print*, ' # of  y points = ', Ny

CsB___Determine the number of qT points in the grid
      Length = Length/Ny
      Call Sample (2, lHead, Filename, 2, Length, NC, NqT)
!      Print*, ' # of qT points = ', NqT
!      Print*, ' Length of file, NQ*Ny*NqT =', LengthFile, NQ*Ny*NqT
      If (LengthFile.Eq.NQ*Ny*NqT) then
        Print*, ' File checksum OK '
      Else
        Print*, ' File checksum BAD. The grid file is corrupt! '
        Stop
      EndIf
CCPY      If (NQ.GT.50 .or. Ny.GT.75 .or. NqT.GT.150) then
CZL      If (NQ.GT.50 .or. Ny.GT.85 .or. NqT.GT.150) then
c htj    If (NQ.GT.50 .or. Ny.GT.100 .or. NqT.GT.200) then
      If (NQ.GT.100 .or. Ny.GT.500 .or. NqT.GT.200) then
c        print*,"NQ=",NQ,"Ny=",NY,"Nqt=",Nqt
        Print*, ' Grid file too large.  See comments in the code.'
CsB     These dimensions are limited to 200, 100 and 50.  
C       To change search for (in both resbos.f and resbosm.f):
C "C      The dimensions below give the maximal possible gridsize"
C "        Parameter (MD1 = 200, MD2 = 100, MD3 = 50)"

        Stop
      EndIf

      Return
      End ! CheckSum

CsB -------------------------------------------------------------------------
      Subroutine MatchGrids1
CsB -------------------------------------------------------------------------
C
C     The following notation apply:
C     For the CSS (in case of NLO run: Asymptotic) piece
C     X1A (X2A, X3A) contains the qT (y, Q) vector
C     iD1+1 (iD2+1, iD3+1) is the number of qT (y, Q) values in the grid
C     The 1st index of YAS and YAA belong to qT, the 2nd to y, the 3rd to Q.
C     For the Y piece:
C     X1Y (X2Y, X3Y) contains the qT (y, Q) vector
C     iY1+1 (iY2+1, iY3+1) is the number of qT (y, Q) values in the grid
C     The 1st index of YAY belong to qT, the 2nd to y, the 3rd to Q.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

CsB   The dimensions below give the maximal possible gridsize.
      Parameter (MD1 = 200, MD2 = 300, MD3 = 100)
      Dimension X1A(0:MD1), X2A(0:MD2), X3A(0:MD3)
      Dimension YAS(0:MD1,0:MD2,0:MD3), YAA(0:MD1,0:MD2,0:MD3)
      Dimension X1Y(0:MD1), X2Y(0:MD2), X3Y(0:MD3)
      Dimension YAY(0:MD1,0:MD2,0:MD3,0:7)
      Dimension YDAS(0:MD1,0:MD2,0:MD3), YDAA(0:MD1,0:MD2,0:MD3)
      Common / Arrays / X1A,X2A,X3A, YAS,YAA,  X1Y,X2Y,X3Y, YAY
      Common / Dims / iD1,iD2,iD3, iY1,iY2,iY3

      Integer NPart
      COMMON/AA/W,QCDL,X1X2MIN,NJET,N_PRIMARY_PART,NF,I_STR,IBEAM,NPART,
     &  NEVENT,NORDER,N_SECOND_PART
      Character*40 Boson
      Common / Boson / Boson

CCPY
      LOGICAL Y_ONE_SIDE
      COMMON/GRID_FILE_FORM/ Y_ONE_SIDE
      Common / MatchGridCom / iMatchOut

      CHARACTER*145 Filename,YFilename, KPFileName, KYFileName
      Common/GridName/ Filename, YFilename
      Character*78 dummy

      Logical Testing
      Data Testing /.False./

CCPY  DO NOT PERFORM MATCHING FOR QTV SMALLER THAN QTV_LOW_CUT
      REAL*8 QTV_LOW_CUT

cdump
C      Testing=.true.
      
CsB___Use 4 points in each direction in the interpolation.
      N_QT_INTP=4
      N_Y_INTP=4
      N_Q_INTP=4

      If (iMatchOut.Eq.1) then
        Print*, ' Writing matched grid '
        Open(Unit=9,File='Matched.out') !,Status='New')
      End If

CsB___Dump the header of the CSS grid
CsB   Find out the length of its header
      lHead = iHeadLen(2,Filename)
CsB   Dump the header
        OPEN(UNIT=2,FILE=Filename,STATUS='OLD') !,ReadOnly)
        Do i = 1, lHead-1
          Read(2,'(A)',Err=99,End=99) dummy
          If (iMatchOut.Eq.1) Write(9,'(A)') dummy
        End Do
        Close (2)
CsB   Modify the last line
      Write(9,'(A)') '  Q,qT,y, CSS,CSS+Y,Pert,CSS+Y matched to Pert'

CsB___Loop over Q and y
      Do 230 iQ = 0,iY3

      Do 231 iy = 0,iY2 ! write full y range
CCPY Nov 2009:      Do 231 iy = 1,iY2 ! don't write first y value
        Qv = X3Y(iQ)
        yv = X2Y(iy)
        If (Qv.GT.X3A(iD3)) Goto 230
        If (yv.GT.X2A(iD2)) Goto 231

CCPY Sept 2009: Assuming X2A(1)=0.0d0, i.e. y=0.0 is the first y-bin in 
C one-side y-grid file.

      IF(Y_ONE_SIDE) THEN
CCPY Use 2nd polynomial to do the interpolation.
        IF(DABS(YV).LE.X2A(2)) N_Y_INTP=4
C        print*,' YV,X2A(1),X2A(2) =', YV,X2A(1),X2A(2)
      EndIf

CsB_____Refill the L0 and A3 parts of the Y array:
        YAY(0,iy,iQ,0) = 0.d0
        YAY(0,iy,iQ,3) = 0.d0

        iCrossS = 0
        iCrossA = 0

cpn2005 This is to avoid the fake crossing count at the lowest q_T on
cpn2005 the grid
        Tot0Prev = 0d0
        Tot3Prev = 0d0

CsB_____Loop over qT
        Do 232 iqT = 1,iY1
          qTv = X1Y(iqT)
          If (qTv.GT.X1A(iD1)) Goto 232
CsB_______Interpolate the resummed piece
          Call ThreeDInterp ( N_QT_INTP,N_Y_INTP,N_Q_INTP,
     >      YDAS,YDAA, qTv,yv,Qv, SigS,SigA, ErrSgS,ErrSgA )

CsB_______Define: Tot=CSS+Y and Pert
          Tot0 = SigS + YAY(iqT,iy,iQ,0)
          Tot3 = SigA + YAY(iqT,iy,iQ,3)
          Per0 = YAY(iqT,iy,iQ,5)
          Per3 = YAY(iqT,iy,iQ,6)
          Per0Prev = YAY(iqT-1,iy,iQ,5)
          Per3Prev = YAY(iqT-1,iy,iQ,6)

          IF(Boson.EQ.'GL' .OR. Boson.EQ.'GG') THEN
CCPY Jan 16, 2011
C Do not include Y contribution for qTv < QTV_LOW_CUT          
            QTV_LOW_CUT=10.0
            IF(qTv.LT.QTV_LOW_CUT) THEN
              Tot0 = SigS
              Tot3 = SigA
            ENDIF
          ENDIF              

CCPY  DO NOT PERFORM MATCHING FOR QTV SMALLER THAN QTV_LOW_CUT
          IF(Boson.EQ.'H0') THEN
            QTV_LOW_CUT=2.5
          ELSE
CCPY Sept 2009            QTV_LOW_CUT=1.0
            QTV_LOW_CUT=1.5

          ENDIF
CCPY April 2018: At the NNLO, when the magnitude of Y-piece is large, we could try this:
          IF(qTv.LT.QTV_LOW_CUT) THEN
            Tot0 = SigS
            Tot3 = SigA
          ENDIF

CsB_______Check whether CSS+Y and Pert cross
C         L0 piece
          If(Per0Prev-Tot0Prev.Eq.0.) Goto 220
cpn2007   Condition for the crossing in the "matching region", defined
cpn2007 that Q_T should not too small
          If(qTv.gt.QTV_LOW_CUT.and.
     >       (Per0-Tot0)*(Per0Prev-Tot0Prev).le.0)
     >      iCrossS = iCrossS + 1
C         A3 piece
          If (Per3Prev-Tot3Prev.Eq.0.) Goto 220
          If(qTv.gt.QTV_LOW_CUT.and.
     >       (Per3-Tot3)*(Per3Prev-Tot3Prev).le.0)
     >      iCrossA = iCrossA + 1
  220     Continue

CsB_______This temporary fix is for W production at the 1.8 TeV Tevatron
C         to force matching for Q ~ 80 GeV
C         It is not needed when the correct kinematics 
C         (x_{1,2}=M_T^2/S e^(+-y)) is used.
c          If (Boson.Eq.'W+' .and. W.LT.1.9e3) then
c            If (qTv.GT.80.) then 
c              iCrossS = 2
c              iCrossA = 2
c            EndIf
c          EndIf

CsB_______Reload the L0 (symmetric) piece
          If (iCrossS.LT.2) then
            YAY(iqT,iy,iQ,0) = Tot0
          Else If (iCrossS.GE.2) then
            YAY(iqT,iy,iQ,0) = Per0
          End If

CsB_______Reload the A3 (anti-symmetric) piece
          If (iCrossA.LT.2) then
            YAY(iqT,iy,iQ,3) = Tot3
          Else If (iCrossA.GE.2) then
            YAY(iqT,iy,iQ,3) = Per3
          End If

CsB_______Inspect the matched grid
          If (iMatchOut.Eq.1) then
            Write (9,101)
C            Print 101,
     ,      X3Y(iQ),X1Y(iqT),X2Y(iy),
C     ,      YAY(iqT,iy,iQ,0),YAY(iqT,iy,iQ,5),YAY(iqT,iy,iQ,3),
C     ,      YAY(iqT,iy,iQ,0),YAY(iqT,iy,iQ,3),YAY(iqT,iy,iQ,6)
C           CSS,CSS+Y,Pert,Matched
     ,      SigS,Tot0,Per0,YAY(iqT,iy,iQ,0),YAY(iqT,iy,iQ,7)
     ,      ,iCrossS
C     ,      ,iCrossA
  101     Format(1X,3(G9.3,1X),5(G11.5,1X),1(I6))
C  101       Format(1X,3(G9.3,1X),3(G9.3,1X),2(I6))
          End If

          Tot0Prev = Tot0
          Tot3Prev = Tot3

  232   Continue ! End loop on qT
  231 Continue ! End loop on y
  230 Continue ! End loop on Q

      If (iMatchOut.Eq.1) Close (9)
C        Print*, ' Stopping in MatchGrids1 '
C        Stop

      Return

 99   Print*, ' Error at reading data in MatchGrids1. '
      Stop

      End ! MatchGrids1

CsB * ***********************************************************************
CsB -------------------------------------------------------------------------
      Subroutine MatchGrids2
CsB -------------------------------------------------------------------------
C
C     The following notation apply:
C     For the CSS (in case of NLO run: Asymptotic) piece
C     X1A (X2A, X3A) contains the qT (y, Q) vector
C     iD1+1 (iD2+1, iD3+1) is the number of qT (y, Q) values in the grid
C     The 1st index of YAS and YAA belong to qT, the 2nd to y, the 3rd to Q.
C     For the Y piece:
C     X1Y (X2Y, X3Y) contains the qT (y, Q) vector
C     iY1+1 (iY2+1, iY3+1) is the number of qT (y, Q) values in the grid
C     The 1st index of YAY belong to qT, the 2nd to y, the 3rd to Q.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

CsB   The dimensions below give the maximal possible gridsize.
      Parameter (MD1 = 200, MD2 = 300, MD3 = 100)
      Dimension X1A(0:MD1), X2A(0:MD2), X3A(0:MD3)
      Dimension YAS(0:MD1,0:MD2,0:MD3), YAA(0:MD1,0:MD2,0:MD3)
      Dimension X1Y(0:MD1), X2Y(0:MD2), X3Y(0:MD3)
      Dimension YAY(0:MD1,0:MD2,0:MD3,0:7)
      Dimension YDAS(0:MD1,0:MD2,0:MD3), YDAA(0:MD1,0:MD2,0:MD3)
      Common / Arrays / X1A,X2A,X3A, YAS,YAA,  X1Y,X2Y,X3Y, YAY
      Common / Dims / iD1,iD2,iD3, iY1,iY2,iY3

      Integer NPart
      COMMON/AA/W,QCDL,X1X2MIN,NJET,N_PRIMARY_PART,NF,I_STR,IBEAM,NPART,
     &  NEVENT,NORDER,N_SECOND_PART
      Character*40 Boson
      Common / Boson / Boson

CCPY
      LOGICAL Y_ONE_SIDE
      COMMON/GRID_FILE_FORM/ Y_ONE_SIDE

      Common / MatchGridCom / iMatchOut

      CHARACTER*145 Filename,YFilename, KPFileName, KYFileName
      Common/GridName/ Filename, YFilename
      Character*78 dummy

      Logical Testing
      Data Testing /.False./

CsB___Use 4 points in each direction in the interpolation.
      N_QT_INTP=4
      N_Y_INTP=4
      N_Q_INTP=4

      If (iMatchOut.Eq.1) then
        Print*, ' Writing matched grid '
        Open(Unit=9,File='Matched.out') !,Status='New')
      End If

CsB___Dump the header of the CSS grid
CsB   Find out the length of its header
      lHead = iHeadLen(2,Filename)
CsB   Dump the header
        OPEN(UNIT=2,FILE=Filename,STATUS='OLD') !,ReadOnly)
        Do i = 1, lHead-1
          Read(2,'(A)',Err=99,End=99) dummy
          Write(9,'(A)') dummy
        End Do
        Close (2)
CsB   Modify the last line
      Write(9,'(A)') '  Q,qT,y, CSS,CSS+Y,Pert,CSS+Y matched to Pert'

      noCrossFlag = 0

CsB___Loop over Q and y
      Do 230 iQ = 0,iY3
!      Do 231 iy = 0,iY2 ! write full y range
      Do 231 iy = 1,iY2 ! don't write first y value
        Qv = X3Y(iQ)
        yv = X2Y(iy)
        If (Qv.GT.X3A(iD3)) Goto 230
        If (yv.GT.X2A(iD2)) Goto 231

CCPY Sept 2009: Assuming X2A(1)=0.0d0, i.e. y=0.0 is the first y-bin in 
C one-side y-grid file.

      IF(Y_ONE_SIDE) THEN
CCPY Use 2nd polynomial to do the interpolation.
        IF(DABS(YV).LE.X2A(2)) N_Y_INTP=4
c        print*,' YB_V,X2A(1),X2A(2) =', YB_V,X2A(1),X2A(2)
      EndIf

CsB_____Refill the L0 and A3 parts of the Y array:
        YAY(0,iy,iQ,0) = 0.d0
        YAY(0,iy,iQ,3) = 0.d0

        iCrossS = 0
        iCrossA = 0

cpn2005 This is to avoid the fake crossing count at the lowest q_T on
cpn2005 the grid
        Tot0Prev = 0d0
        Tot3Prev = 0d0


CsB_____Loop over qT to find the second crossing point
        Do 233 iqT = 1,iY1
          qTv = X1Y(iqT)
          If (qTv.GT.X1A(iD1)) Goto 233
CsB_______Interpolate the resummed piece
          Call ThreeDInterp ( N_QT_INTP,N_Y_INTP,N_Q_INTP,
     >      YDAS,YDAA, qTv,yv,Qv, SigS,SigA, ErrSgS,ErrSgA )

CsB_______Define: Tot=CSS+Y and Pert
          Tot0 = SigS + YAY(iqT,iy,iQ,0)
          Tot3 = SigA + YAY(iqT,iy,iQ,3)
          Per0 = YAY(iqT,iy,iQ,5)
          Per3 = YAY(iqT,iy,iQ,6)
          Per0Prev = YAY(iqT-1,iy,iQ,5)
          Per3Prev = YAY(iqT-1,iy,iQ,6)

CsB_______Check whether CSS+Y and Pert cross
C         L0 piece
          If(Per0Prev-Tot0Prev.Eq.0.) Goto 221
          Cross = (Per0-Tot0)/(Per0Prev-Tot0Prev)
          If (Cross.LE.0.) iCrossS = iCrossS + 1
C         A3 piece
          If (Per3Prev-Tot3Prev.Eq.0.) Goto 221
          Cross = (Per3-Tot3)/(Per3Prev-Tot3Prev)
          If (Cross.LE.0.) iCrossA = iCrossA + 1
  221     Continue

CsB_______Store the crossing point
          If (iCrossS.Eq.1) iqTCPtS = iqT
          If (iCrossA.Eq.1) iqTCPtA = iqT
C          If (iCrossS.Eq.2) Print*, Qv, yv, qTCPtS
C          If (iCrossA.Eq.2) Print*, Qv, yv, qTCPtA
C          If (iCrossS*iCrossA.Eq.4) Stop

          Tot0Prev = Tot0
          Tot3Prev = Tot3

  233   Continue ! End 1st loop on qT

!	If (iCrossS.Eq.0) then
!	  Print*, ' Q,y, qTCPtS, qTCPtA = ', 
!     ,            Qv, yv, X1Y(iqTCPtS+1), X1Y(iqTCPtA+1)
!	  Stop
!	End If
  
CsB_____If crossing doesn't occur, find the point where CSS+Y and Pert 
C       are the closest.
        If (iCrossS.Eq.0 .or. iCrossA.Eq.0) then
        noCrossFlag = noCrossFlag + 1
        Do 244 iqT = 1,iY1
          qTv = X1Y(iqT)
          If (qTv.GT.X1A(iD1)) Goto 244
CsB_______Interpolate the resummed piece
          Call ThreeDInterp ( N_QT_INTP,N_Y_INTP,N_Q_INTP,
     >      YDAS,YDAA, qTv,yv,Qv, SigS,SigA, ErrSgS,ErrSgA )

CsB_______Define: Tot=CSS+Y and Pert
          Tot0 = SigS + YAY(iqT,iy,iQ,0)
          Tot3 = SigA + YAY(iqT,iy,iQ,3)
          Per0 = YAY(iqT,iy,iQ,5)
          Per3 = YAY(iqT,iy,iQ,6)
          Per0Prev = YAY(iqT-1,iy,iQ,5)
          Per3Prev = YAY(iqT-1,iy,iQ,6)

CsB_______Check whether CSS+Y and Pert approach each other
C         L0 piece
          If (Per0Prev-Tot0Prev.LE.Per0-Tot0 .and. iCrossS.Eq.0) 
     .       iCrossS = 1
C         A3 piece
          If (Per3Prev-Tot3Prev.LE.Per3-Tot3 .and. iCrossA.Eq.0) 
     .       iCrossA = 1
  222     Continue

CsB_______Store the crossing point
          If (iCrossS.Eq.1) iqTCPtS = iqT - 1
          If (iCrossA.Eq.1) iqTCPtA = iqT - 1
C          If (iCrossS.Eq.2) Print*, Qv, yv, qTCPtS
C          If (iCrossA.Eq.2) Print*, Qv, yv, qTCPtA
C          If (iCrossS*iCrossA.Eq.4) Stop

          Tot0Prev = Tot0
          Tot3Prev = Tot3

  244   Continue ! End loop on qT
        End If ! iCrossS.Eq.0 .or. iCrossA.Eq.0

CsB_____Loop over qT to modify values
        Do 232 iqT = 1,iY1
          qTv = X1Y(iqT)
          If (qTv.GT.X1A(iD1)) Goto 232
CsB_______Interpolate the resummed piece
          Call ThreeDInterp ( N_QT_INTP,N_Y_INTP,N_Q_INTP,
     >      YDAS,YDAA, qTv,yv,Qv, SigS,SigA, ErrSgS,ErrSgA )

CsB_______Define the matching functions F for Symm and Asym pieces
          If (iqT.LT.iqTCPtS+1 .And. X1Y(iqTCPtS+1).NE.0.) then
            FS = (1.d0 - qTv**2/X1Y(iqTCPtS+1)**2)**2
          Else
            FS = 0.d0
          End If
          If (iqT.LT.iqTCPtA+1 .And. X1Y(iqTCPtA+1).NE.0.) then
            FA = (1.d0 - qTv**2/X1Y(iqTCPtA+1)**2)**2
          Else
            FA = 0.d0
          End If

CsB_______Define: Tot = F*CSS + Pert - F*Asy
          ASY0 = YAY(iqT,iy,iQ,5)-YAY(iqT,iy,iQ,0)
          Tot0 = FS*SigS + YAY(iqT,iy,iQ,5) - FS*ASY0
          ASY3 = YAY(iqT,iy,iQ,6)-YAY(iqT,iy,iQ,3)
          Tot3 = FA*SigA + YAY(iqT,iy,iQ,6) - FA*ASY3
c          Per0 = YAY(iqT,iy,iQ,5)
c          Per3 = YAY(iqT,iy,iQ,6)

CsB_______Reload the L0 (symmetric) piece
            YAY(iqT,iy,iQ,0) = Tot0

CsB_______Reload the A3 (anti-symmetric) piece
            YAY(iqT,iy,iQ,3) = Tot3

CsB_______Inspect the matched grid
          If (iMatchOut.Eq.1) then
            Write (9,101)
C            Print 101,
     ,      X3Y(iQ),X1Y(iqT),X2Y(iy),
C     ,      YAY(iqT,iy,iQ,0),YAY(iqT,iy,iQ,5),YAY(iqT,iy,iQ,3),
C     ,      YAY(iqT,iy,iQ,0),YAY(iqT,iy,iQ,3),YAY(iqT,iy,iQ,6)
C           CSS,CSS+Y,Pert,Matched
     ,      SigS,SigS+YAY(iqT,iy,iQ,5)-ASY0,YAY(iqT,iy,iQ,5),Tot0,FS
C     ,      ,iCrossS
C     ,      ,iCrossA
  101       Format(1X,3(G9.3,1X),6(G11.5,1X))
C  101       Format(1X,3(G9.3,1X),4(G9.3,2X),1(I6))
C  101       Format(1X,3(G9.3,1X),3(G9.3,1X),2(I6))
          End If
          
  232   Continue ! End 2nd loop on qT

  231 Continue ! End loop on y
  230 Continue ! End loop on Q

      If (noCrossFlag.GE.1) 
     .   Print*, ' Xing doesn`t happen for some Q and y !'

      If (iMatchOut.Eq.1) Close (9)
C      Print*, ' Stopping in MatchGrids2 '
C      Stop

      Return

 99   Print*, ' Error at reading data in MatchGrids2. '
      Stop

      End ! MatchGrids2

C----------------------------------------------------------------------
      SUBROUTINE SET_ELEMENT_TAGS
C----------------------------------------------------------------------
C-
C    Define the number of elements and their names (tags) for events
C    written into ntuples
C----------------------------------------------------------------------
      IMPLICIT NONE
      Integer lunit, istatus, icycle, nt_id, nt_elements,
     >  nt_iwgt0_elements
      Character*80 nt_file_name
      Character*8  nt_subdir
      Common / NTPL / nt_id, lunit, nt_subdir, nt_file_name, nt_elements

      integer max_siz
      parameter(max_siz=50)
      Character*80 chtitl
      Character*8  chtags(max_siz)
      common /NT_TAGS/chtitl, chtags

      Character*40 Boson
      Common / Boson / Boson
      Character*2 hDecayMode
      Common / hDM / hDecayMode

      INTEGER*4 NBINIT, iwgt, GetNumBranches
      external GetNumBranches
      REAL*8 MAX_WGT,NEW_MAX_WGT,WGT_FACTOR
      COMMON/UNWEIGHT/MAX_WGT,NEW_MAX_WGT,WGT_FACTOR,iwgt,NBINIT

      Character*8 Output
      Common / Output / Output

c   CHTITL      is the character title, a human readable label associated
c        with it (but is NOT looked at by any of the HBOOK routines)
c   CHTAGS      is an array of character "tags" associated with each element
c        (they must be defined for each column in order to access it via PAW)
c
c       additional "tag" rules:
c         a. up to 8 characters long (but no longer)
c         b. for compatibility reasons it is useful for them to
c             be all uppercase characters, so FORTRAN routines
c             distinquished between them
c         c. start with alphabetic character
c         d. include only alphanumerics [A-Z][0-9] and "_"

      IF(BOSON.NE.'H0') GOTO 10


      IF (hDecayMode.EQ.'ZZ' .OR. hDecayMode.EQ.'WW') THEN
         nt_iwgt0_elements=17 !number of elements if iwgt=0 or -1 is chosen
      ELSE
         GO TO 10
      ENDIF

      If (HDecayMode.Eq.'ZZ') then
        chtitl        = 'h1 h2 -> Z Z X -> l L l` L` X'
      ELSEIf (HDecayMode.Eq.'WW') then
        chtitl        = 'h1 h2 -> W W X -> l N L n X'
      endif                     

      if (iwgt.eq.-2) then  !Save only the weight
        nt_elements=1
        if (Output.Eq.'ROOTNT1'.or.Output.Eq.'ROOTNT2') then
          write(chtags(1),'(''WT'',i2.2)') 
     >      GetNumBranches()-nt_iwgt0_elements+1
          print *,GetNumBranches(),'----------',chtags(1),'==========='
       else
          chtags(1)='WTXX'
        endif                   !Output
      else              !iwgt.ne.2; write both kinematical variables and weights
        nt_elements=nt_iwgt0_elements
        If (Output.Eq.'ROOTNT1'.or.Output.Eq.'PAWNT1'
     >    .or.Output.Eq.'Text1') then !4-momentum output
          chtags(1)  = 'Px_d1'
          chtags(2)  = 'Py_d1'
          chtags(3)  = 'Pz_d1'
          chtags(4)  = 'E_d1'
          chtags(5)  = 'Px_d2'
          chtags(6)  = 'Py_d2'
          chtags(7)  = 'Pz_d2'
          chtags(8)  = 'E_d2'
          chtags(9)  = 'Px_nn'
          chtags(10)  = 'Py_nn'
          chtags(11)  = 'Pz_nn'
          chtags(12)  = 'E_nn'
          chtags(13)  = 'Px_H'
          chtags(14) = 'Py_H'
          chtags(15) = 'Pz_H'
          chtags(16) = 'E_H'
          chtags(17) = 'WT00'
        ElseIf (Output.Eq.'ROOTNT2'.or.Output.Eq.'PAWNT2'
     >      .or.Output.Eq.'Text2') then !"physical" momentum variables 
          chtags(1)  = 'M_H'
          chtags(2)  = 'pT_H'
          chtags(3)  = 'y_H'
          chtags(4)  = 'pT_d1'
          chtags(5)  = 'y_d1'
          chtags(6)  = 'pT_d2'
          chtags(7)  = 'y_d2'
          chtags(8)  = 'cos_the_dd'
          chtags(9)  = 'pT_dd'
          chtags(10) = 'Del_phi_dd'
          chtags(11) = 'Dy_dd'
          chtags(12) = 'DelR_dd'
          chtags(13) = 'mis_ET'
          chtags(14) = 'MT_d1'
          chtags(15) = 'MT_d2'
          chtags(16) = 'MT_H'
          chtags(17) = 'WT00'

        Endif                   !Output
      endif                     ! iwgt

      GOTO 999

10	CONTINUE
C================================
C yfu: you need to change the maximum number of this parameter
      nt_iwgt0_elements=16 !number of elements if iwgt=0 or -1 is chosen

      chtitl        = 'h1 h2 -> B(-> d1[e] d2[n]) X'

      if (iwgt.eq.-2) then  !Save only the weight
        nt_elements=1
        if (Output.Eq.'ROOTNT1'.or.Output.Eq.'ROOTNT2') then
          write(chtags(1),'(''WT'',i2.2)') 
     >      GetNumBranches()-nt_iwgt0_elements+1
          print *,GetNumBranches(),'----------',chtags(1),'==========='
       else
          chtags(1)='WTXX'
        endif                   !Output
      else              !iwgt.ne.2; write both kinematical variables and weights
        nt_elements=nt_iwgt0_elements
        If (Output.Eq.'ROOTNT1'.or.Output.Eq.'PAWNT1'
     >    .or.Output.Eq.'Text1') then !4-momentum output
          chtags(1)  = 'Px_d1'
          chtags(2)  = 'Py_d1'
          chtags(3)  = 'Pz_d1'
          chtags(4)  = 'E_d1'
          chtags(5)  = 'Px_d2'
          chtags(6)  = 'Py_d2'
          chtags(7)  = 'Pz_d2'
          chtags(8)  = 'E_d2'
          chtags(9)  = 'Px_V'
          chtags(10) = 'Py_V'
          chtags(11) = 'Pz_V'
          chtags(12) = 'E_V'
C yfu add cos theta of hadron level
          chtags(13) = 'costhe_h'
C yfu add cos theta of quark level
          chtags(14) = 'costhe_q'
C yfu add a tag for qqbar or qbarq event
C qqbar=1 qbarq=-1 all=0
          chtags(15) = 'Type'
          chtags(16) = 'WT00'
        ElseIf (Output.Eq.'ROOTNT2'.or.Output.Eq.'PAWNT2'
     >      .or.Output.Eq.'Text2') then !"physical" momentum variables 
          chtags(1)  = 'pT_d1'
          chtags(2)  = 'pT_d2'
          chtags(3)  = 'y_d1'
          chtags(4)  = 'y_d2'
          chtags(5)  = 'pT_B'
          chtags(6)  = 'y_B'
          chtags(7)  = 'M_B'
          chtags(8)  = 'D_phi'
          chtags(9)  = 'cos_the_sta'
          chtags(10) = 'phi_sta'
          chtags(11) = 'DelR34'
          chtags(12) = 'WT00'
          chtags(13) = 'dummy'
cpn2006 Alternative choices for the output
c        chtags(11) = 'M35'
c        chtags(13) = 'M45'
c        chtags(13) = 'zCS' !Catani-Seymour's fragmentation fraction z
        Endif                   !Output
      endif                     ! iwgt

 999  RETURN
      END ! set_element_tags


C***************************************************************************

C----------------------------------------

cpn2008 Fill the output in the StdHEP format
      SUBROUTINE StdHEPFill(WEIGHT)
      IMPLICIT NONE

      REAL*4 WEIGHT
      REAL*8 W,QCDL,X1X2MIN
      INTEGER*4 NJET,N_PRIMARY_PART,NF,I_STR,IBEAM,NPART,NEVENT,NORDER,
     ,          N_SECOND_PART
      COMMON/AA/W,QCDL,X1X2MIN,NJET,N_PRIMARY_PART,NF,I_STR,IBEAM,NPART,
     &  NEVENT,NORDER,N_SECOND_PART
      REAL*8 SEC_PLAB,PLAB,PCEM,WHAT
      COMMON/DECAY_MOMENTA/SEC_PLAB(4,50)
      COMMON/MOM/PLAB(4,12),PCEM(4,12),WHAT
      REAL*8 ETOTMN,PTMIN,ETAMAX,ETASEP
      INTEGER*4 IAPPROX
      COMMON/CUTVAL/ETOTMN,PTMIN,ETAMAX,ETASEP,IAPPROX

      REAL*8 PTECUT,YECUT,DELR_CUT,PTNCUT,yNCut,YWLCUT,YWUCUT
      COMMON/MYCUT/PTECUT,YECUT,DELR_CUT,PTNCUT,yNCut,YWLCUT,YWUCUT
      REAL*8 YM12,YPT12,YP12
      REAL*8 YM12MIN,YM12MAX,QTMin,QTMax, TMmin, TMMax, Etmin
      COMMON/MYYM12/YM12MIN,YM12MAX,QTMin,QTMax, TMMin, TMMax, Etmin
 
      COMMON/MY1/YM12,YPT12,YP12(4)
       REAL*8 YXP,YXETA,YXPT,YXDELR,YXPP
      COMMON/SECCONV/ YXP(4,30),YXETA(30),YXPT(30),YXDELR(30,30),
     >YXPP(30)
      Character*40 Boson
      Common / Boson / Boson
      Integer iYPiece, iReadIn
      Common / YPiece / iYPiece, iReadIn

      REAL*8 QT_V_MIN, QT_V_MAX,YB_V_MAX,Q_V_MIN,Q_V_MAX,QT_V_MXY
      COMMON/HARDCUT/ QT_V_MIN, QT_V_MAX,YB_V_MAX,Q_V_MIN,
     >  Q_V_MAX,QT_V_MXY

      REAL*8 Q_V,QT_V,YB_V,THE_STA,PHI_STA
      COMMON/CSFRAME/Q_V,QT_V,YB_V,THE_STA,PHI_STA

      INTEGER*4 I,J,NFINAL
      REAL*8 YXDELETA,YXCOS1,YXCOS2,YXSIN1,YXSIN2,YXCOSDPH,YXDELPHI
      LOGICAL NOJET
      REAL*8 ZCORR,ZCORR_CUT

      Integer ACut
      Common / AAGCut / ACut
      Real*8 CosTheta
      REAL*8 XMW,XMZ,XMW2,XMZ2,AEM,GAMW,GAMZ,GAMW2,GAMZ2
      COMMON/STAND0/XMW,XMZ,XMW2,XMZ2,AEM,GAMW,GAMZ,GAMW2,GAMZ2

      REAL*8 DELTA_SOFT
      Common /TheoCutoff/ DELTA_SOFT

CZHU  ID AND MASS OF THE PARENT AND DAUGHTER PARTICLES
      REAL*8 BOSON_VECT(4), BOSON_MASS
      INTEGER ID_PARENT, ID_PART1, ID_PART2, ID_PART3
      REAL*8 MASS_PART1, MASS_PART2, MASS_PART3

      INTEGER*4 NCALL1,ITMX1,NCALL2,ITMX2
      COMMON/VEGASI/NCALL1,ITMX1,NCALL2,ITMX2

      INTEGER EVENT_THROW, EVENT_ALL

      INTEGER*4 ID_LEPTON,ID_ELE,ID_MU,ID_TAU
      COMMON/LEPTONID/ID_LEPTON,ID_ELE,ID_MU,ID_TAU

      REAL*8 RANDOM
      

C     4 VECTORS BEFORE AND AFTER PHOTOS
C     AFTER PHOTOS, AT MOST 5 PARTICLES
      INTEGER INDEX, NPHOTOS, ISSTABLE
      REAL*8 VECT1(4,3),VECT2(4,5)
      INTEGER IDPARTS(5)

cpn2008
      integer photoscor, event_num, decaymode, wpm
      real*4 mymaxwgt
      common/stdhep/mymaxwgt, photoscor, event_num, decaymode, wpm
 
cpn2008 Temporary plug
      print *,'Photoscor, Decaymode and WPM temporarily not defined'
      call Quit

      
CZHU --------------------------------------------------------
CZHU     START TO WRITE OUT INFORMATION
CZHU --------------------------------------------------------
      
CZHU   ONLY WRITE OUT THE LAST NCALL2 EVENTS
      EVENT_NUM = EVENT_NUM + 1
      
      IF(MYMAXWGT.LE.WEIGHT) THEN
        MYMAXWGT = WEIGHT
      ENDIF
      
cpn2008
c$$$      IF(MOD(EVENT_NUM,10000).EQ.0) THEN
c$$$        PRINT*, '# OF EVENTS GENERATED: ', EVENT_NUM
c$$$      ENDIF
      
CZHU  WRITE OUT EVENT NUMBER AND EVENT WEIGHT
      WRITE(18, *) EVENT_NUM, WEIGHT
      
      NFINAL = 4
      IF(NPART.EQ.15.OR.NPART.EQ.16)THEN
        NFINAL=5
      ENDIF

CZHU  Construct pT for the W/Z boson (Generator-level)
      DO 602 I=1, 4
         YP12(I)=0.D0
         
CZHU  add the px, py components of leptons and photons to get the pT of W/Z boson
CZHU  note: PLAB(, 1) and PLAB(, 2) are not filled
         DO 600 J=3, NFINAL
           YP12(I)=YP12(I)+PLAB(I,J)
 600     CONTINUE
 602   CONTINUE
       
CZHU  Construct pT (magnitude) of the W/Z boson and the mass
       YPT12 = DSQRT(YP12(1)**2+YP12(2)**2)
       YM12  = DSQRT(YP12(4)**2-YP12(3)**2-YPT12**2)
       
CZHU   add the four-momenta of the W/Z boson
       DO 515 I=1, 4
         BOSON_VECT(I)=YP12(I)
 515   CONTINUE
       
C     add the mass of the W/Z boson
       BOSON_MASS=YM12

CZHU  *********************************
CZHU  SETUP PARTICLE ID AND MASS
CZHU  *********************************

CZHU  Use ISAJET particle ID convertion
CZHU  http://www-d0.fnal.gov/~serban/isajet/ptable.txt

CZHU  FOR W-
       IF ((BOSON(1:2).EQ.'W-')) THEN
         ID_PARENT = -80
         
         IF(DECAYMODE.EQ.1) THEN
           ID_PART1 = 12
           ID_PART2 = -11
           MASS_PART1 = 0.000511
         ELSE IF(DECAYMODE.EQ.2) THEN
           ID_PART1 = 14
           ID_PART2 = -13
           MASS_PART1 = 0.105658
         ELSE IF(DECAYMODE.EQ.3) THEN
           ID_PART1 = 16
           ID_PART2 = -15
           MASS_PART1 = 1.77699
         ENDIF
         MASS_PART2 = 0.        ! MASS OF ANTI-ELECTRON, MUON AND TAU NEUTRINOS
         
CZHU  FOR W+
       ELSE IF (BOSON(1:2).EQ.'W+') THEN
         ID_PARENT = 80
         
         IF(DECAYMODE.EQ.1) THEN
           ID_PART1 = 11
           ID_PART2 = -12
           MASS_PART2 = 0.000511
         ELSE IF (DECAYMODE.EQ.2) THEN
           ID_PART1 = 13
           ID_PART2 = -14
           MASS_PART2 = 0.105658
         ELSE IF (DECAYMODE.EQ.2) THEN
           ID_PART1 = 15
           ID_PART2 = -16
           MASS_PART2 = 1.77699
         ENDIF
         MASS_PART1 = 0.        ! MASS OF ELECTRON, MUON AND TAU NEUTRINOS
         
CZHU  FOR Z0
         
       ELSE IF (BOSON(1:2).EQ.'Z0' .OR.  BOSON(1:2).EQ.'A0' .OR.
     >     BOSON(1:2).EQ.'DY') THEN
         ID_PARENT = 90
         IF (DECAYMODE.EQ.1) THEN
           ID_PART1 = 12
           ID_PART2 = -12
           MASS_PART1 = 0.000511
           MASS_PART2 = 0.000511
         ELSE IF (DECAYMODE.EQ.2) THEN
           ID_PART1 = 14
           ID_PART2 = -14
           MASS_PART1 = 0.105658
           MASS_PART2 = 0.105658
         ELSE IF (DECAYMODE.EQ.3) THEN
           ID_PART1 = 16
           ID_PART2 = -16
           MASS_PART1 = 1.77699
           MASS_PART2 = 1.77699
         ENDIF                  ! ID_LEPTON
         
       ENDIF                    !W+, W-, Z0

CZHU  *************************************************
CZHU  IF WE WANT TO GENERATE W+, W- or W+ and W- EVENTS
CZHU  **************************************************
       IF ((BOSON(1:2).EQ.'W+').OR.(BOSON(1:2).EQ.'W-')) THEN

CZHU  PROVIDE EITHER W+ OR W- GRID FILE AND WANT TO GENERATE W+ AND W- EVENTS TOGETHER
         IF(WPM.EQ.0) THEN
CZHU  USE A RANDOM NUMBER TO DETERMINE WHETHER FLIP SIGN
           RANDOM=RAND()
           IF(RANDOM.LE.0.5D0) THEN
             ID_PARENT=-ID_PARENT
             ID_PART1=-ID_PART1
             ID_PART2=-ID_PART2
             
             BOSON_VECT(1)=-BOSON_VECT(1)
             BOSON_VECT(2)=-BOSON_VECT(2)
             BOSON_VECT(3)=-BOSON_VECT(3)
             
             PLAB(1,3)=-PLAB(1,3)
             PLAB(2,3)=-PLAB(2,3)
             PLAB(3,3)=-PLAB(3,3)
             
             PLAB(1,4)=-PLAB(1,4)
             PLAB(2,4)=-PLAB(2,4)
             PLAB(3,4)=-PLAB(3,4)
           ENDIF
         ENDIF

CZHU  PROVIDE W+ GRID FILE AND WANT TO GENERATE W- EVENTS OR
CZHU  PROVIDE W- GRID FILE AND WANT TO GENERATE W+ EVENTS
         IF( ((BOSON(1:2).EQ.'W+').AND.(WPM.EQ.2)).OR.
     *     ((BOSON(1:2).EQ.'W-').AND.(WPM.EQ.1)) ) THEN
           ID_PARENT=-ID_PARENT
           ID_PART1=-ID_PART1
           ID_PART2=-ID_PART2
           
           BOSON_VECT(1)=-BOSON_VECT(1)
           BOSON_VECT(2)=-BOSON_VECT(2)
           BOSON_VECT(3)=-BOSON_VECT(3)
           
           PLAB(1,3)=-PLAB(1,3)
           PLAB(2,3)=-PLAB(2,3)
           PLAB(3,3)=-PLAB(3,3)
           
           PLAB(1,4)=-PLAB(1,4)
           PLAB(2,4)=-PLAB(2,4)
           PLAB(3,4)=-PLAB(3,4)
         ENDIF
       ENDIF

CZHU  *************************************************
CZHU  IF WE WANT TO DO FSR USING PHOTOS
CZHU  NEED TO MODIFY 4-VECTORS
CZHU  *************************************************
       IF(PHOTOSCOR.EQ.1) THEN
CZHU  SAVE WZ AND LEPTONS 4 VECTORS
         DO 18  INDEX=1,4
           VECT1(INDEX,1)=BOSON_VECT(INDEX)
           VECT1(INDEX,2)=PLAB(INDEX,3)
           VECT1(INDEX,3)=PLAB(INDEX,4)
 18      CONTINUE

cpn2008         CALL PHOTOS_RAD(EVENT_NUM, ID_PARENT, ID_PART1,
c     &     ID_PART2, VECT1, VECT2, IDPARTS, NPHOTOS) 
       ENDIF
       
CZHU  *************************************************
CZHU  WRITE OUT EVENT INFORMATION INTO A TEXT FILE
CZHU  **************************************************
CZHU  WRITE OUT THE PRIMARY VERTEX
       WRITE(18, *) 0., 0., 0.

CZHU  *************************************************
CZHU  FOR NORMAL RESBOS+PHOTOS OUTPUT
CZHU  *************************************************
       IF(PHOTOSCOR.EQ.1) THEN
         DO 19 I=1,NPHOTOS

CZHU  STATUS OF THE PARTICLE
CZHU  W/Z UNSTABLE, OTHER PARTICLES STABLE
           IF((ABS(IDPARTS(I)).EQ.80).OR.
     &       (ABS(IDPARTS(I)).EQ.90)) THEN
             ISSTABLE=0
           ELSE
             ISSTABLE=1
        ENDIF
        
        WRITE(18, *) IDPARTS(I),VECT2(1,I),VECT2(2,I),
     &    VECT2(3,I), VECT2(4,I), ISSTABLE, 0
 19   CONTINUE 
      
CZHU  *************************************************
CZHU  FOR NORMAL RESBOS OUTPUT
CZHU  *************************************************
      ELSE 
        WRITE(18, *) ID_PARENT,BOSON_VECT(1),BOSON_VECT(2),
     &    BOSON_VECT(3), BOSON_VECT(4), 0, 0
        WRITE(18, *) ID_PART1,PLAB(1, 3),PLAB(2, 3),
     &    PLAB(3, 3),PLAB(4, 3), 1, 0
        WRITE(18, *) ID_PART2,PLAB(1, 4),PLAB(2, 4),
     &    PLAB(3, 4),PLAB(4, 4), 1, 0
      ENDIF
      
CZHU  WRITE THE PHOTON INFORMATION
CZHU  MOST TIME THIS SHOULD NOT HAPPEN
       IF(NFINAL.EQ.5) THEN
        WRITE(18, *) 10,PLAB(1, 5),PLAB(2, 5),PLAB(3, 5),
     &       PLAB(4, 5), 1, 0
       ENDIF

CZHU  MARK THE END OF AN EVENT
       WRITE(18, *), 0

      RETURN
      END                       ! StdHEPFill


C######################################################################

