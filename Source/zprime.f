! ======================================================================
      program zprime
! ----------------------------------------------------------------------
! Header
  ! Authors: Declan Millar, Stefano Moretti.

  ! Calculates the cross section and generates distributions for
  !   pp -> tt,
  !   pp -> tt -> bW^+bbarW^- -> bbbare^+nue^-nubarc
  ! (Future:) pp -> bW^+bbarW^- -> b bbar e^+ nu qqbar'
  ! Uses adapted Madgraph functions.
  ! Uses cteq6 and mrs99 PDF subroutines.
! ----------------------------------------------------------------------
! Declarations
  ! implicit
      implicit real*8 (a-h,p-z)
      implicit integer (l-o)

  ! Global variables     
  !   vegas
      common/bveg1/ncall,itmx,nprn,ndev,xl(100),xu(100),acc
      common/bveg2/it,ndo,si,swgt,schi,xi(50,100)
      common/rndm/iseed
      common/reslocal/resl(20),standdevl(20)
  !   Kinematics
      common/par/rm3,rm4,rm5,rm6,rm7,rm8,s
      common/limfac/fac
      common/EW/a_em,s2w
      common/final/o_final,ipmax   
      common/stat/npoints
      common/symint/ixmax,jxmax
      common/coll/ecm_coll
      common/cuts/ytmax,yttmin
  !   Debugging
      common/debug/o_M_eq_1
  !   Polarised/Spatial cross sections
      integer nasym
      parameter (nasym=9)
      integer nspat
      parameter (nspat=6) !nasym-3
      common/polarised/polcross(20,-1:1,-1:1),polerror(20,-1:1,-1:1)
      common/spatial/spatcross(nspat,20,-1:1),spaterror(nspat,20,-1:1)
  !   Permitted gauge sectors
      common/igauge/o_QCD,o_EW,o_BSM
  !   Interference       
      common/interference/o_int      
  !   Z' masses and VA/LR couplings
      common/fermions/ fmass,     fwidth
      dimension fmass(12), fwidth(12)
      common/vmass1/rm_W,Gamma_W,rm_Z,Gamma_Z
      common/vmass2/rm_A,Gamma_a,rm_h,Gamma_h
      common/Zp/rmZp(5),gamZp(5)
      common/Zpparam/paramZp(5)
      common/coupZpVA/gp(5),gV_d(5),gA_d(5),gV_u(5),gA_u(5)
      common/coupZp/gZpd(2,5),gZpu(2,5)
  !   Narrow width approximation (NWA)
      common/NWA/o_NWA
  !   Structure functions
      common/partdist/o_structure
      common/QCD/rlambdaQCD4,nloops
      common/collider/o_coll
  !   Switch for all distributions
      common/distros/o_distros
  !   Distributions in pTs of external particles
      common/ext_pT/pTmax(8),pTmin(8),pTw(8)
      common/dist_pT/xpT(8,500),fxpT(8,500,20),fxpTtot(8,500)
      common/inp_pT/o_pT(8)
      common/div_pT/ndiv_pT(8)
  !   Distributions in etas of external particles
      common/ext_eta/etamax(8),etamin(8),etaw(8)
      common/dist_eta/xeta(8,500),fxeta(8,500,20),fxetatot(8,500)
      common/inp_eta/o_eta(8)
      common/div_eta/ndiv_eta(8)
  !   Distributions in phis of external particles
      common/ext_phi/phimax(8),phimin(8),phiw(8)
      common/dist_phi/xphi(8,500),fxphi(8,500,20),fxphitot(8,500)
      common/inp_phi/o_phi(8)
      common/div_phi/ndiv_phi(8)
  !   Distribution in ETmiss
      common/ext_ETmiss/ETmissmax,ETmissmin,ETmissw
      common/dist_ETmiss/xETmiss(500),fxETmiss(500,20),fxETmisstot(500)
      common/inp_ETmiss/o_ETmiss
      common/div_ETmiss/ndiv_ETmiss
  !   Distribution in pT of the top
      common/ext_pT356/pT356max,pT356min,pT356w
      common/dist_pT356/xpT356(500),fxpT356(500,20),fxpT356tot(500)
      common/inp_pT356/o_pT356
      common/div_pT356/ndiv_pT356
  !   Distribution in eta of the top
      common/ext_eta356/eta356max,eta356min,eta356w
      common/dist_eta356/xeta356(500),fxeta356(500,20),fxeta356tot(500)
      common/inp_eta356/o_eta356
      common/div_eta356/ndiv_eta356
  !   Distribution in phi of the top
      common/ext_phi356/phi356max,phi356min,phi356w
      common/dist_phi356/xphi356(500),fxphi356(500,20),fxphi356tot(500)
      common/inp_phi356/o_phi356
      common/div_phi356/ndiv_phi356
  !   Distribution in pT of the anti-top
      common/ext_pT478/pT478max,pT478min,pT478w
      common/dist_pT478/xpT478(500),fxpT478(500,20),fxpT478tot(500)
      common/inp_pT478/o_pT478
      common/div_pT478/ndiv_pT478
  !   Distribution in eta of the anti-top
      common/ext_eta478/eta478max,eta478min,eta478w
      common/dist_eta478/xeta478(500),fxeta478(500,20),fxeta478tot(500)
      common/inp_eta478/o_eta478
      common/div_eta478/ndiv_eta478
  !   Distribution in phi of the anti-top
      common/ext_phi478/phi478max,phi478min,phi478w
      common/dist_phi478/xphi478(500),fxphi478(500,20),fxphi478tot(500)
      common/inp_phi478/o_phi478
      common/div_phi478/ndiv_phi478     
  !   Distribution in invarient mass of the top pair
      common/ext_rMtt/rMttmax,rMttmin,rMttw
      common/dist_rMtt/xrMtt(500),fxrMtt(500,20),fxrMtttot(500)
      common/inp_rMtt/o_rMtt
      common/div_rMtt/ndiv_rMtt
  !   Distribution in boost of top pair centre of mass frame
      common/ext_beta/betamax,betamin,betaw
      common/dist_beta/xbeta(500),fxbeta(500,20),fxbetatot(500)
      common/inp_beta/o_beta
      common/div_beta/ndiv_beta
  !   Distribution in cos(theta_t)
      common/ext_cost/costmax,costmin,costw
      common/dist_cost/xcost(500),fxcost(500,20),fxcosttot(500)
      common/inp_cost/o_cost
      common/div_cost/ndiv_cost
  !   Distribution in top energy
      common/ext_Et/Etmax,Etmin,Etw
      common/dist_Et/xEt(500),fxEt(500,20),fxEttot(500)
      common/inp_Et/o_Et
      common/div_Et/ndiv_Et
  !   Distributions in transverse variables
      integer ntrans
      parameter (ntrans=10)
      common/ext_trans/transmax(ntrans),transmin(ntrans),transw(ntrans)
      common/dist_trans/xtrans(ntrans,500),fxtrans(ntrans,500,20)
     &                                           ,fxtranstot(ntrans,500)
      common/inp_trans/o_tran(ntrans)
      common/div_trans/ndiv_trans(ntrans)
      dimension sfxtranstot(ntrans)
  !   Distribution in phi_l (lepton azimuthal angle)
      common/ext_fl/flmax,flmin,flw
      common/dist_fl/xfl(500),fxfl(500,20),fxfltot(500)
      common/inp_fl/o_fl
      common/div_fl/ndiv_fl
  !   Distribution in cos_phi_l
      common/ext_cosfl/cosflmax,cosflmin,cosflw
      common/dist_cosfl/xcosfl(500),fxcosfl(500,20),fxcosfltot(500)
      common/inp_cosfl/o_cosfl
      common/div_cosfl/ndiv_cosfl
  !   Distribution in delta phi
      common/ext_dphi/dphimax,dphimin,dphiw
      common/dist_dphi/xdphi(500),fxdphi(500,20),fxdphitot(500)
      common/inp_dphi/o_dphi
      common/div_dphi/ndiv_dphi
  !   Distribution in cost5
      common/ext_cost5/cost5max,cost5min,cost5w
      common/dist_cost5/xcost5(500),fxcost5(500,20),fxcost5tot(500)
      common/inp_cost5/o_cost5
      common/div_cost5/ndiv_cost5
  !   Distribution in cost7
      common/ext_cost7/cost7max,cost7min,cost7w
      common/dist_cost7/xcost7(500),fxcost7(500,20),fxcost7tot(500)
      common/inp_cost7/o_cost7
      common/div_cost7/ndiv_cost7
  !   Distribution in cost7
      common/ext_ct7ct5/ct7ct5max,ct7ct5min,ct7ct5w
      common/dist_ct7ct5/xct7ct5(500),fxct7ct5(500,20),fxct7ct5tot(500)
      common/inp_ct7ct5/o_ct7ct5
      common/div_ct7ct5/ndiv_ct7ct5  
  !   Distributions in asymmetries
      common/inp_asym/o_asym(nasym)    
  !   Distribution in sigp
      common/ext_sigp/sigpmax,sigpmin,sigpw
      common/dist_sigp/xsigp(1000),fxsigp(nasym,1000,20)
     &,fxsigptot(nasym,1000)
      common/inp_sigp/o_sigp
      common/div_sigp/ndiv_sigp
  !   Distribution in sigm
      common/ext_sigm/sigmmax,sigmmin,sigmw
      common/dist_sigm/xsigm(1000),fxsigm(nasym,1000,20)
     &,fxsigmtot(nasym,1000)
      common/inp_sigm/o_sigm
      common/div_sigm/ndiv_sigm

  ! Local variables
  !   Flag for Zp width specification
      dimension o_width(5)
  !   Model name
      character*50 model
  !   Polarised/hemispherised cross sections
      dimension cnorm(20)
      !dimension snorm(6)  !,ave(4) 
      dimension poltot(-1:1,-1:1),polchi(-1:1,-1:1)
      dimension spattot(nspat,-1:1),spatchi(nspat,-1:1)
      dimension sfxpTtot(8),sfxetatot(8),sfxphitot(8)
      dimension sfxsigptot(nasym),sfxsigmtot(nasym)
      dimension asym_int(nasym)
      dimension Atot(nasym),Atoterr(nasym)
  ! test particular matrix element
      dimension testp1(0:3),testp2(0:3),
     &          testp3(0:3),testp4(0:3)

  ! Local constants
  !   pi
      parameter (pi=3.14159265358979323846d0)
  !   Unit conversion GeV -> nb
      parameter (conv=0.38937966d9)
  !   Date and time
      integer today(3), now(3)
  !   Branching ratio for t->bev=bmv=btv (with QCD corrections?)     
  !    real*8 BRtbln/0.10779733d0/
  !   Branching ratio for t->bev=bmv=btv=1/9 (tree level)      
      real*8 BRtbln/0.11111111d0/ 
  !   Branching ratio for t->beq=bmq=btq=6/9 (tree level)      
      real*8 BRtbeq/0.66666666d0/

  ! External procedures
      external dxsec
! ----------------------------------------------------------------------
! Read input files
  ! Read config file
  !   Collider flag (o_coll=0: pp; o_coll=1: ppbar)
      read(5,*) o_coll
  !   Collider energy     
      read(5,*) ecm_coll
  !   PDFs
      read(5,*) o_structure  
  !   Name of model file
      read(5,*) model      
  !   Permitted gauge sector options
      read(5,*) o_QCD
      read(5,*) o_EW
      read(5,*) o_BSM
  !   Interference options
      read(5,*) o_int
  !   Final state option (o_final=ttbar: no decay; o_final=1: bbllnn)      
      read(5,*) o_final
  !   NWA flag (o_NWA = 0: Actual top widths; o_NWA = 1: tops in NWA)
      read(5,*) o_NWA
  !   Branching ratio flag
      read(5,*) o_BR
  !   Transverse mass variables flag
      read(5,*) o_trans
  !   Asymmetry observable flag
      read(5,*) o_asyms
  !   Cut on top rapidity 
      read(5,*) ytmax
  !   Cut on top pair boost
      read(5,*) yttmin
  !   Random number seed
      read(5,*) iseed
  !   Maximum number of Vegas iterations
      read(5,*) itmx      
  !   Number of Vegas calls per iteration
      read(5,*) ncall
  !   Desired Accuracy (If negative, run maximum iterations.)
      read(5,*) acc  
  !   Symmatrise of x1 and x2
      read(5,*) o_symx1x2
  !   Symmatrise of x1 and x2
      read(5,*) o_symcost      
  !   Standard distributions flag
      read(5,*) o_distros   
  !   set |M|^2=1
      read(5,*) o_M_eq_1
        
  ! Interpret config
    ! Number of external lines
      if(o_final.eq.0)then 
        ipmax=4
      else if(o_final.eq.1)then 
        ipmax=8
      else 
        write(*,*)'Invalid final state identifier!'
        stop
      end if
    ! NWA only for six-body final state
      if(o_final.eq.0) o_NWA=0
    ! itmx no more than 20.      
      if(itmx.gt.20)then
        write(*,*)'itmx does not have to exceed 20!'
        stop
      end if
    ! For every point in phase space with x1 and x2, include the point
    ! in phase space with x1<->x2
      if(o_symx1x2.eq.1)then
        ixmax=2
      else
        ixmax=1
      end if
    ! in phase space with cost->-cost
      if(o_symcost.eq.1)then
        jxmax=2
      else
        jxmax=1
      end if
    ! in phase space with cost->-cost
      if(o_M_eq_1.eq.1)then
        o_QCD=0
        o_EW=0
        o_BSM=0
      else
        jxmax=1
      end if
    ! Extract model filename (Remove white space.)
      imodel = len(model)
      do while(model(imodel:imodel).eq.'') 
        imodel = imodel-1
      end do

  ! Read model file
      open(unit=42,file='Models/'//model(1:imodel)//'.mdl',status='old')
      read(42,*) rmZp
      read(42,*) gamZp      
      read(42,*) gp 
      read(42,*) paramZp
      read(42,*) gV_u
      read(42,*) gA_u
      read(42,*) gV_d
      read(42,*) gA_d      
  !   Check whether width has been specified
  !   (If gamZp is zero, the function widthZp is used instead.)
      do i=1,5
        if ((gamZp(i).eq.0d0).and.(rmZp(i).gt.0d0)) then
          o_width(i) = 0
        else
          o_width(i) = 1
        end if
      enddo
! ----------------------------------------------------------------------
! Distributions Setup
  ! (Set flags, binning range and divisions.)
  ! pT distributions
      do ip=1,ipmax
        o_pT(ip)=o_distros
        pTmax(ip)=7000.d0/(1+o_coll*6)
        pTmin(ip)=0.d0
        ndiv_pT(ip)=70
  ! eta distributions
        o_eta(ip)=o_distros
        etamax(ip)=+10
        etamin(ip)=-10
        ndiv_eta(ip)=50
  ! phi distributions
        o_phi(ip)=o_distros
        phimax(ip)=+pi
        phimin(ip)=-pi
        ndiv_phi(ip)=50
      end do
  !   missing transverse momentum
      o_ETmiss=o_distros
      ETmissmax=7000.d0/(1+o_coll*6)
      ETmissmin=0.d0
      ndiv_ETmiss=70      
  !   top transverse momentum
      o_pT356=o_distros
      pT356max=7000.d0/(1+o_coll*6)
      pT356min=0.d0
      ndiv_pT356=70
  !   2to6 top pseudorapidity
      o_eta356=o_distros
      eta356max=+10
      eta356min=-10
      ndiv_eta356=50 
  !   2to6 top pseudorapidity
      o_phi356=o_distros
      phi356max=+pi
      phi356min=-pi
      ndiv_phi356=50        
  !   anti-top transverse momentum
      o_pT478=o_distros
      pT478max=7000.d0/(1+o_coll*6)
      pT478min=0.d0
      ndiv_pT478=70
  !   2to6 anti-top pseudorapidity
      o_eta478=o_distros
      eta478max=+10
      eta478min=-10
      ndiv_eta478=50
  !   2to6 top pseudorapidity
      o_phi478=o_distros
      phi478max=+pi
      phi478min=-pi
      ndiv_phi478=50
  !   invarient mass of tt pair (always on)
      o_rMtt=1
      rMttmax=14000.d0/(1+o_coll*6)
      rMttmin=0.d0
      ndiv_rMtt=140
  !   boost of parton CoM
      o_beta=o_distros
      betamax=1000.d0
      betamin=0.d0
      ndiv_beta=100
  !   costheta
      o_cost=o_distros
      costmax=+1.d0
      costmin=-1.d0
      ndiv_cost=50
  !   top energy
      o_Et=o_distros
      Etmax=7000.d0/(1+o_coll*6)
      Etmin=0.d0
      ndiv_Et=70

  !   transverse variables
      do itrans=1,ntrans
        if(o_final.eq.0)then
          o_tran(itrans)=0
        else 
          o_tran(itrans)=o_trans
        end if
      end do
  !   invarient mass of the visible decay products of the tt pair
      transmax(1)=4000
      transmin(1)=0.d0
      ndiv_trans(1)=40
  !   sum of tranvserse energy  
      transmax(2)=4000
      transmin(2)=0.d0
      ndiv_trans(2)=40   
  !   transverse mass 1
      transmax(3)=4000
      transmin(3)=0.d0
      ndiv_trans(3)=40
  !   transverse mass 2
      transmax(4)=4000
      transmin(4)=0.d0
      ndiv_trans(4)=40
  !   transverse mass 3
      transmax(5)=4000
      transmin(5)=0.d0
      ndiv_trans(5)=40
   !  lepton transverse mass
      transmax(6)=500
      transmin(6)=0.d0
      ndiv_trans(6)=40     
  !   contransverse mass 1
      transmax(7)=4000
      transmin(7)=0.d0
      ndiv_trans(7)=40    
  !   contransverse mass 2
      transmax(8)=4000
      transmin(8)=0.d0
      ndiv_trans(8)=40            
  !   contransverse mass 3
      transmax(9)=4000
      transmin(9)=0.d0
      ndiv_trans(9)=40
  !   lepton contransverse mass 
      transmax(10)=500
      transmin(10)=0.d0
      ndiv_trans(10)=50     

  !   phi_l
      o_fl=o_asyms
      flmax=+2*pi
      flmin=0
      ndiv_fl=100
  !   cosphi_l
      o_cosfl=o_asyms
      cosflmax=+1.d0
      cosflmin=-1.d0
      ndiv_cosfl=100
  !   delta phi
      o_dphi=o_asyms
      dphimax=+pi
      dphimin=0
      ndiv_dphi=10
  !   cost5
      o_cost5=o_asyms
      cost5max=+1
      cost5min=-1
      ndiv_cost5=10
  !   cost7
      o_cost7=o_asyms
      cost7max=+1
      cost7min=-1
      ndiv_cost7=10
   !  ct7ct5
      o_ct7ct5=o_asyms
      ct7ct5max=+1
      ct7ct5min=-1
      ndiv_ct7ct5=10   
  !   sigp
      o_sigp=o_asyms
      sigpmax=rMttmax
      sigpmin=rMttmin
      ndiv_sigp=ndiv_rMtt/5
  !   sigm
      o_sigm=o_asyms
      sigmmax=rMttmax
      sigmmin=rMttmin
      ndiv_sigm=ndiv_rMtt/5      

  !   asymmetries
      do iasy=1,nasym
        o_asym(iasy)=o_asyms
      end do

  !   Turn off 2->6 only distributions
      if (o_final.eq.0)then
        do ip=5,8
          o_pT(i)   = 0
          o_eta(i)  = 0
          o_phi(i)  = 0          
        end do
        do itrans=1,ntrans
          o_tran(itrans)=0
        end do
        o_tran  = 0
        o_pT356  = 0
        o_eta356 = 0
        o_phi356 = 0
        o_pT478  = 0
        o_eta478 = 0
        o_phi478 = 0
        o_ETmiss = 0
        o_HT     = 0
        o_fl     = 0
        o_dphi   = 0
        o_cosfl  = 0
        o_cost7  = 0
        o_cost5  = 0
        o_ct7ct5 = 0
        o_asym(9) = 0  ! turn off A_l
      end if
  !   Turn off 2->2 only distributions
      if (o_final.eq.1)then
        o_asym(1) = 0 ! turn off A_LL
        o_asym(2) = 0 ! turn off A_L
        o_asym(3) = 0 ! turn off A_PV
      end if
! ----------------------------------------------------------------------
! Set-up physics

  ! Collider CM energy squared.      
      s=ecm_coll*ecm_coll

  ! Factor outside integration
  !   Conversion GeV^-2 -> pb
      fac=conv
  !   Azimuthal angle integrated out (No initial transverse polarisation.)
      fac=fac*2.d0*pi      

  ! QCDL4 is QCD LAMBDA4 (to match PDF fits).
  ! (PDFs are intrinsically linked to the value of lamda_QCD; alpha_QCD) 
      if(o_structure.eq.1)qcdl4=0.326d0
      if(o_structure.eq.2)qcdl4=0.326d0
      if(o_structure.eq.3)qcdl4=0.326d0
      if(o_structure.eq.4)qcdl4=0.215d0
      if(o_structure.eq.5)qcdl4=0.300d0
      if(o_structure.eq.6)qcdl4=0.300d0
      if(o_structure.eq.7)qcdl4=0.300d0
      if(o_structure.eq.8)qcdl4=0.229d0
      if(o_structure.eq.9)qcdl4=0.383d0
      rlambdaQCD4=QCDL4
  ! Initialise CTEQ grids.
      if(o_structure.le.4)then
        icteq=o_structure
        call SetCtq6(ICTEQ)
      end if

  ! Use appropriately evolved alphas.
      if(o_structure.eq.1)nloops=2
      if(o_structure.eq.2)nloops=2
      if(o_structure.eq.3)nloops=1
      if(o_structure.eq.4)nloops=1
      if(o_structure.eq.5)nloops=1
      if(o_structure.eq.6)nloops=1
      if(o_structure.eq.7)nloops=1
      if(o_structure.eq.8)nloops=1
      if(o_structure.eq.9)nloops=1

  ! initialise MadGraph - masses and coupling constants of particles
      call initialise_madGraph(o_NWA,model)

  ! Calculate Zp couplings
      call coupZpx

  ! Calculate sequential Zp widths
      do i=1,5
        if (o_width(i).eq.0) gamZp(i)=
     &             widthZp(rm_W,rm_Z,rmZp(i),a_em,s2w,rlambdaQCD4,nloop)
      end do
! ---------------------------------------------------------------------- 
! VEGAS parameters
  ! Dimensions of integration
      if(o_final.eq.0)then
        ndim=3
      else if(o_final.eq.1)then
        ndim=15
      end if
  !   (If nprn<0 no print-out.)
      nprn=0
      if(o_final.eq.0)then
  !   Final state masses     
        rm3=fmass(11)
        rm4=rm3
        rm5=0.d0
        rm6=0.d0
        rm7=0.d0
        rm8=0.d0
  !   Integrates on:
  !   x(3)=(x1-tau)/(1-tau),
  !   x(2)=(ecm-rm3-rm4)/(ecm_max-rm3-rm4),
  !   x(1)=cos(theta3_cm)
  !   Limits:
        do i=3,2,-1
          xl(i)=0.d0
          xu(i)=1.d0
        end do
  !         if(isycost.eq.1)then  ! might not work
  !           nctpoints = 200
  !           do i=1,1
  !             xl(i)=0.d0
  !             xu(i)=0.d0
  !           end do
  !         else
  !           nctpoints = 0
        do i=1,1
          xl(i)=-1.d0
          xu(i)=1.d0
        end do
  !         end if

      else if(o_final.eq.1)then
  !   Final state masses
        rm3=fmass(12)
        rm4=rm3
        rm5=0.d0
        rm6=0.d0
        rm7=0.d0
        rm8=0.d0
  !   Integrates on:
 
  !   x(15)=(x1-tau)/(1-tau),
  !   x(14)=(ecm-rm3-rm4-rm5-rm6-rm7-rm8)
  !        /(ecm_max-rm3-rm4-rm5-rm6-rm7-rm8),
  !   x(13)=(XX356-XX356min)/(XX356max-XX356min),
  !   where XX356=arctg((rm356**2-rm3**2)/rm3/gamt),
  !   x(12)=(XX478-XX478min)/(XX478max-XX478min),
  !   where XX478=arctg((rm478**2-rm3**2)/rm3/gamt),
  !   x(11)=(XX56-XX56min)/(XX56max-XX56min),
  !   where XX56=arctg((rm56**2-rm_W**2)/rm_W/gamW),
  !   x(10)=(XX78-XX78min)/(XX78max-XX78min),
  !   where XX78=arctg((rm78**2-rm_W**2)/rm_W/gamW),
  !   x(9)=cos(theta_cm_356)=-cos(theta_cm_478)  
  !   x(8)=cos(theta56_cm_356),
  !   x(7)=cos(theta78_cm_478),
  !   x(6)=cos(theta5_cm_56),
  !   x(5)=cos(theta7_cm_78),
  !   x(4)=fi56_cm_356,
  !   x(3)=fi78_cm_478,
  !   x(2)=fi5_cm_56,
  !   x(1)=fi8_cm_78;
  !   Limits:
        do i=15,14,-1
          xl(i)=0.d0
          xu(i)=1.d0
        end do
        do i=13,10,-1
          xl(i)=0.d0
          xu(i)=1.d0
        end do
        do i=9,5,-1
          xl(i)=-1.d0
          xu(i)=1.d0
        end do
        do i=4,1,-1
          xl(i)=0.d0
          xu(i)=2.d0*pi
        end do
      end if
! ----------------------------------------------------------------------
! Generate bins
  ! (Finds bin width, finds midpoints.)

      do ip=3,ipmax
        if(o_pT(ip).eq.1)then
          pTw(ip)=(pTmax(ip)-pTmin(ip))/ndiv_pT(ip)
          do j=1,ndiv_pT(ip)
            xpT(ip,j)=pTmin(ip)+pTw(ip)*(j-1)+pTw(ip)/2.d0
          end do
        end if
        if(o_eta(ip).eq.1)then
          etaw(ip)=(etamax(ip)-etamin(ip))/ndiv_eta(ip)
          do j=1,ndiv_eta(ip)
            xeta(ip,j)=etamin(ip)+etaw(ip)*(j-1)+etaw(ip)/2.d0
          end do
        end if
        if(o_phi(ip).eq.1)then
          phiw(ip)=(phimax(ip)-phimin(ip))/ndiv_phi(ip)
          do j=1,ndiv_phi(ip)
            xphi(ip,j)=phimin(ip)+phiw(ip)*(j-1)+phiw(ip)/2.d0
          end do
        end if
      end do

      if(o_ETmiss.eq.1)then
        ETmissw=(ETmissmax-ETmissmin)/ndiv_ETmiss
        do i=1,ndiv_ETmiss
          xETmiss(i)=ETmissmin+ETmissw*(i-1)+ETmissw/2.d0
        end do
      end if

      if(o_pT356.eq.1)then
        pT356w=(pT356max-pT356min)/ndiv_pT356
        do i=1,ndiv_pT356
          xpT356(i)=pT356min+pT356w*(i-1)+pT356w/2.d0
        end do
      end if

      if(o_eta356.eq.1)then
        eta356w=(eta356max-eta356min)/ndiv_eta356
        do i=1,ndiv_eta356
          xeta356(i)=eta356min+eta356w*(i-1)+eta356w/2.d0
        end do
      end if

      if(o_phi356.eq.1)then
        phi356w=(phi356max-phi356min)/ndiv_phi356
        do i=1,ndiv_phi356
          xphi356(i)=phi356min+phi356w*(i-1)+phi356w/2.d0
        end do
      end if

      if(o_pT478.eq.1)then
        pT478w=(pT478max-pT478min)/ndiv_pT478
        do i=1,ndiv_pT478
          xpT478(i)=pT478min+pT478w*(i-1)+pT478w/2.d0
        end do
      end if      

      if(o_eta478.eq.1)then
        eta478w=(eta478max-eta478min)/ndiv_eta478
        do i=1,ndiv_eta478
          xeta478(i)=eta478min+eta478w*(i-1)+eta478w/2.d0
        end do
      end if

      if(o_phi478.eq.1)then
        phi478w=(phi478max-phi478min)/ndiv_phi478
        do i=1,ndiv_phi478
          xphi478(i)=phi478min+phi478w*(i-1)+phi478w/2.d0
        end do
      end if

      if(o_rMtt.eq.1)then
        rMttw=(rMttmax-rMttmin)/ndiv_rMtt
        do i=1,ndiv_rMtt
          xrMtt(i)=rMttmin+rMttw*(i-1)+rMttw/2.d0
        end do
      end if

      if(o_beta.eq.1)then
        betaw=(betamax-betamin)/ndiv_beta
        do i=1,ndiv_beta
          xbeta(i)=betamin+betaw*(i-1)+betaw/2.d0
        end do
      end if
      if(o_cost.eq.1)then
        costw=(costmax-costmin)/ndiv_cost
        do i=1,ndiv_cost
          xcost(i)=costmin+costw*(i-1)+costw/2.d0
        end do
      end if

      if(o_Et.eq.1)then
        Etw=(Etmax-Etmin)/ndiv_Et
        do i=1,ndiv_Et
          xEt(i)=Etmin+Etw*(i-1)+Etw/2.d0
        end do
      end if

      do itrans=1,ntrans
        if(o_tran(itrans).eq.1)then
          transw(itrans)=(transmax(itrans)-transmin(itrans))
     &                                               /ndiv_trans(itrans)
          do i=1,ndiv_trans(itrans)
            xtrans(itrans,i)=transmin(itrans)+transw(itrans)*(i-1)
     &                                              +transw(itrans)/2.d0
          end do
        end if
      end do

      if(o_fl.eq.1)then
        flw=(flmax-flmin)/ndiv_fl
        do i=1,ndiv_fl
          xfl(i)=flmin+flw*(i-1)+flw/2.d0
        end do
      end if

      if(o_cosfl.eq.1)then
        cosflw=(cosflmax-cosflmin)/ndiv_cosfl
        do i=1,ndiv_cosfl
          xcosfl(i)=cosflmin+cosflw*(i-1)+cosflw/2.d0
        end do
      end if

      if(o_dphi.eq.1)then
        dphiw=(dphimax-dphimin)/ndiv_dphi
        do i=1,ndiv_dphi
          xdphi(i)=dphimin+dphiw*(i-1)+dphiw/2.d0
        end do
      end if

      if(o_cost5.eq.1)then
        cost5w=(cost5max-cost5min)/ndiv_cost5
        do i=1,ndiv_cost5
          xcost5(i)=cost5min+cost5w*(i-1)+cost5w/2.d0
        end do
      end if

      if(o_cost7.eq.1)then
        cost7w=(cost7max-cost7min)/ndiv_cost7
        do i=1,ndiv_cost7
          xcost7(i)=cost7min+cost7w*(i-1)+cost7w/2.d0
        end do
      end if

      if(o_ct7ct5.eq.1)then
        ct7ct5w=(ct7ct5max-ct7ct5min)/ndiv_ct7ct5
        do i=1,ndiv_ct7ct5
          xct7ct5(i)=ct7ct5min+ct7ct5w*(i-1)+ct7ct5w/2.d0
        end do
      end if

      if(o_sigp.eq.1)then
        sigpw=(sigpmax-sigpmin)/ndiv_sigp
        do i=1,ndiv_sigp
          xsigp(i)=sigpmin+sigpw*(i-1)+sigpw/2.d0
        end do
      end if

      if(o_sigm.eq.1)then
        sigmw=(sigmmax-sigmmin)/ndiv_sigm
        do i=1,ndiv_sigm
          xsigm(i)=sigmmin+sigmw*(i-1)+sigmw/2.d0
        end do
      end if
! ----------------------------------------------------------------------
! Output information before integration
      write(*,*)'====================================================='
      call idate(today)   ! today(1)=day, (2)=month, (3)=year
      call itime(now)     ! now(1)=hour, (2)=minute, (3)=second
      write(*,*)'DATE ',today(3),today(2),today(1)
      write(*,*)'TIME ',now(1),now(2),now(3)
      write(*,*)'-----------------------------------------------------'
      write(*,*)'PROCESS'
      if(o_coll.eq.0)then
        if(o_final.eq.0)
     &    write(*,*)'pp #rightarrow t#bar{t}',
     &               ' #times BR(t#rightarrow bl#nu)^{2}'
        if(o_final.eq.1)
     &    write(*,*)'pp #rightarrow t#bar{t}',
     &               '#rightarrow b#bar{b} W^{+}W^{-}',
     &               '#rightarrow b#bar{b} l^{+}l^{-} #nu#bar{#nu}'
      else if(o_coll.eq.1)then
        if(o_final.eq.0)
     &    write(*,*)'p#bar{p} #rightarrow t#bar{t}',
     &               ' #times BR(t#rightarrow bl#nu)^{2}'
        if(o_final.eq.1) 
     &    write(*,*)'p#bar{p} #rightarrow t#bar{t}',
     &               '#rightarrow b#bar{b} W^{+}W^{-}',
     &               '#rightarrow b#bar{b} l^{+}l^{-} #nu#bar{#nu}'
      end if
      write(*,*)'-----------------------------------------------------'
      write(*,*)'NOTES'            
      write(*,*)'Units: GeV'
      write(*,*)'Quarks: all massless except t, b.'           
      if(o_structure.eq.1)write(*,*)'PDFs: cteq6m.'
      if(o_structure.eq.2)write(*,*)'PDFs: cteq6d.'
      if(o_structure.eq.3)write(*,*)'PDFs: cteq6l.'
      if(o_structure.eq.4)write(*,*)'PDFs: cteq6l1.'
      if(o_structure.eq.5)write(*,*)'PDFs: mrs99 (cor01).'
      if(o_structure.eq.6)write(*,*)'PDFs: mrs99 (cor02).'
      if(o_structure.eq.7)write(*,*)'PDFs: mrs99 (cor03).'
      if(o_structure.eq.8)write(*,*)'PDFs: mrs99 (cor04).'
      if(o_structure.eq.9)write(*,*)'PDFs: mrs99 (cor05).'
      if((o_final.eq.1).and.(o_NWA.eq.0))write(*,*)'Tops: off-shell.'
      if((o_final.eq.1).and.(o_NWA.eq.1))write(*,*)'Tops: NWA.'
      write(*,*)'BSM model: ',model
      if(o_QCD.eq.1)write(*,*)'QCD: On '
      if(o_QCD.eq.0)write(*,*)'QCD: Off'
      if(o_EW.eq.1) write(*,*)'EW:  On '
      if(o_EW.eq.0) write(*,*)'EW:  Off'
      if(o_BSM.eq.1)write(*,*)'BSM: On '
      if(o_BSM.eq.0)write(*,*)'BSM: Off'   
      if(o_int.eq.0)write(*,*)'Interference: None'
      if(o_int.eq.1)write(*,*)'Interference: SM'
      if(o_int.eq.2)write(*,*)'Interference: Full'
      if(o_int.eq.3)write(*,*)'Interference: No square terms.'
      if(o_M_eq_1.eq.1)write(*,*)'Phase space only'
      if(o_symx1x2.eq.1)write(*,*)'Symmetrical integration over x1<->x2'
      if(o_symcost.eq.1)write(*,*)'Symmetrical integration over cost'
      write(*,*)'iseed: ',iseed 
      write(*,*)'-----------------------------------------------------'      
      write(*,*)'PARAMETERS'
      write(*,*)'#sqrt{s}              ',ecm_coll
      write(*,*)'at |y| <              ',abs(ytmax)
      write(*,*)'Loops a_s evaluated at',nloops
      write(*,*)'a_{s}(M_{Z})          ',alfas(rm_Z,rLambdaQCD4,nloops)
      write(*,*)'#Lambda_{QCD}(4)      ',QCDL4
      write(*,*)'m_{b}                 ',fmass(12)    
      write(*,*)'#Gamma_{b}            ',fwidth(12) 
      write(*,*)'m_{t}                 ',fmass(11)  
      write(*,*)'#Gamma_{t}            ',fwidth(11)   
      write(*,*)'m_{Z}                 ',rm_Z    
      write(*,*)'#Gamma_{Z}            ',gamma_Z   
      write(*,*)'m_{W}                 ',rm_W    
      write(*,*)'#Gamma_{W}            ',gamma_W   
      write(*,*)'m_{H}                 ',rm_h    
      write(*,*)'#Gamma_{H}            ',gamma_h
      write(*,*)'-----------------------------------------------------' 
      write(*,*)'ZPRIME PARAMETERS'
      do i=1,5
        if(rmZp(i).gt.0)then
          write(*,*)'Z#prime               ',i 
          write(*,*)'m_{Z#prime}           ',rmZp(i) 
          write(*,*)'o_width:              ',o_width(i)   
          write(*,*)'#Gamma_{Z#prime}      ',gamZp(i)
          write(*,*)'g_{p}                 ',gp(i)
          write(*,*)'g_{V}^{u}             ',gV_u(i)
          write(*,*)'g_{A}^{u}             ',gA_u(i)
          write(*,*)'g_{V}^{d}             ',gV_d(i)
          write(*,*)'g_{A}^{d}             ',gA_d(i)
          write(*,*)         
        end if
      end do
      write(*,*)'-----------------------------------------------------' 
      write(*,*)'CUTS'      
      write(*,*)'-----------------------------------------------------'
! ----------------------------------------------------------------------
! Integration
  ! Section header
      write(*,*)'INTEGRATION'
  !   Reset counter
      npoints=0  
  !   Reset various iterative quantities
      if(o_final.eq.0)then
        do i=1,20
          resl(i)=0.d0
          standdevl(i)=0.d0
          cnorm(i)=0.d0
          do iphel=-1,+1,2
            do jphel=-1,+1,2
              polcross(i,iphel,jphel)=0.d0
              polerror(i,iphel,jphel)=0.d0
            end do
          end do
          do ispat=1,nspat
            do iasy=-1,+1,2
              spatcross(ispat,i,iasy)=0.d0
              spaterror(ispat,i,iasy)=0.d0
            end do
          end do 
        end do
      end if
  !   Integrate
      it=0
      call vegas(ndim,dxsec,avgi,sd,chi2a)
      if(o_final.eq.0)then
  !   Multiply by branching ratios (if o_final = 0)
        if(o_BR.eq.1)then     
          avgi=avgi*(BRtbln)**2
          sd=sd*(BRtbln)**2
        else 
          continue
        end if
      end if
  !  Collect total cross-section
      cross=avgi
      error=sd

    ! Print integrated cross section
      write(*,*)'-----------------------------------------------------'
      write(*,*)'INTEGRATED CROSS-SECTION'
      if(cross.eq.0d0)then
        write(*,*)'sigma = 0! Check permitted gauge sectors.'
        stop
      else      
        write(*,*)'sigma (pb)','error (same units)'
        write(*,*)cross,error
        write(*,*)'(using ',npoints,' points)'
      end if
  ! Re-weight distributions for different iterations     
      stantot=0.d0
        do i=1,it
          stantot=stantot+1.d0/standdevl(i)/standdevl(i)
        end do
        do i=1,it
          standdevl(i)=standdevl(i)*standdevl(i)*stantot
        end do
        do i=1,it
          cnorm(i)=resl(i)*standdevl(i)
        end do
! ----------------------------------------------------------------------
! Total asymmetries
  ! Collect polarised cross sections.
      if(o_asyms.eq.1)then
        if(o_final.eq.0)then  
          do iphel=-1,+1,2
            do jphel=-1,+1,2
              do i=1,it
                polcross(i,iphel,jphel)=polcross(i,iphel,jphel)
     &                               *avgi/cnorm(i)
                polerror(i,iphel,jphel)=polcross(i,iphel,jphel)
     &                               *sd/cnorm(i)
              end do
              poltot(iphel,jphel)=0.d0
              polchi(iphel,jphel)=0.d0
              do i=1,it
                poltot(iphel,jphel)=poltot(iphel,jphel)
     &                       +polcross(i,iphel,jphel)
                polchi(iphel,jphel)=polchi(iphel,jphel)
     &                       +polerror(i,iphel,jphel)
              end do
              polchi(iphel,jphel)=polchi(iphel,jphel)
     &                         /poltot(iphel,jphel)
    !          polchi(iphel,jphel)=
    !   & sqrt(abs(polchi(iphel,jphel)
    !   &         -poltot(iphel,jphel)**2*dfloat(ncall)))
    !   & /dfloat(ncall)
            end do
          end do
        end if

    ! Collect unpolarised spatial asymmetry
        do ispat=1,nspat
          if(o_asym(ispat+3).eq.0)then
            continue
          else
            do iAB=-1,+1,2
              do i=1,it
                spatcross(ispat,i,iAB)=spatcross(ispat,i,iAB)
     &                            *avgi/cnorm(i)
                spaterror(ispat,i,iAB)=spatcross(ispat,i,iAB)
     &                            *sd/cnorm(i)
              end do        
              spattot(ispat,iAB)=0.d0
              spatchi(ispat,iAB)=0.d0
              do i=1,it   ! add up each iteration
                spattot(ispat,iAB)=spattot(ispat,iAB)
     &                      +spatcross(ispat,i,iAB)
                spatchi(ispat,iAB)=spatchi(ispat,iAB)
     &                      +spaterror(ispat,i,iAB)
              end do
              spatchi(ispat,iAB)=spatchi(ispat,iAB)
     &                        /spattot(ispat,iAB)
    !           spatchi(iasy)=
    !      & sqrt(abs(spatchi(iasy)
    !      &         -spattot(iasy)**2*dfloat(ncall)))
    !      & /dfloat(ncall)
            end do
          end if
        end do

    ! Define asymmetries
        if(o_final.eq.0)then
    ! ALL
          Atot(1)=
     &          +(poltot(+1,+1)-poltot(+1,-1)
     &           -poltot(-1,+1)+poltot(-1,-1))
     &          /cross
          Atoterr(1)=
     &             +(polchi(+1,+1)+polchi(+1,-1)
     &              +polchi(-1,+1)+polchi(-1,-1))
     &             /4.d0*Atot(1)
    ! AL
          Atot(2)=
     &          +(poltot(-1,-1)-poltot(+1,-1) 
     &           +poltot(-1,+1)-poltot(+1,+1))
     &          /cross
          Atoterr(2)=
     &            +(polchi(-1,-1)+polchi(+1,-1) 
     &             +polchi(-1,+1)+polchi(+1,+1))
     &            /4.d0*Atot(2)
    ! APV
          Atot(3)=
     &          +(poltot(-1,-1)-poltot(+1,+1))
     &          /cross/2.d0
          Atoterr(3)=
     &             +(polchi(-1,-1)+polchi(+1,+1))
     &             /2.d0*Atot(3)
        end if

        do iasy=4,nasym
          ispat=iasy-3
          if(o_asym(iasy).gt.0)then
            Atot(iasy)=
     &             +(spattot(ispat,+1)-spattot(ispat,-1))
     &             /cross
            Atoterr(iasy)=
     &               +sd/avgi*Atot(iasy)
           end if
        end do

  ! Print Asymmetries     
        write(*,*)'TOTAL ASYMMETRIES'
        if(o_final.eq.0)then
          write(*,*)'ALL:                  uncertainty (same units):'
          write(*,*)Atot(1),Atoterr(1) 
          write(*,*)'AL:                   uncertainty (same units):'
          write(*,*)Atot(2),Atoterr(2) 
          write(*,*)'APV:                  uncertainty (same units):'
          write(*,*)Atot(3),Atoterr(3) 
          write(*,*)'AFB:                 uncertainty (same units):'
          write(*,*)Atot(4),Atoterr(4)
          write(*,*)'AFB*:                   uncertainty (same units):'
          write(*,*)Atot(5),Atoterr(5)
          write(*,*)'AtRFB:                  uncertainty (same units):'
          write(*,*)Atot(6),Atoterr(6)
          write(*,*)"AttbRFB/A:              uncertainty (same units):"
          write(*,*)Atot(7),Atoterr(7)
          write(*,*)"ARFB/A':              uncertainty (same units):"
          write(*,*)Atot(8),Atoterr(8)
        else if(o_final.gt.0)then
          write(*,*)'A_l:                  uncertainty (same units):'
          write(*,*)Atot(9),Atoterr(9) 
        end if
      end if 
! ----------------------------------------------------------------------
! Plot Distributions
  ! Section header
      write(*,*)'-----------------------------------------------------'
      write(*,*)'HISTOGRAMS'
      do ip=3,8
  ! Plot distributions in pT
        if(o_pT(ip).eq.1)then        
          sfxpTtot(ip)=0d0
          do j=1,ndiv_pT(ip)
            fxpTtot(ip,j)=0.d0
            do i=1,it
              fxpT(ip,j,i)=fxpT(ip,j,i)*avgi/cnorm(i)/pTw(ip)
              fxpTtot(ip,j)=fxpTtot(ip,j)+fxpT(ip,j,i)
            end do
            sfxpTtot(ip)=sfxpTtot(ip)+fxpTtot(ip,j)*pTw(ip)
          end do
          write(*,*)'DISTRIBUTION'
          write(*,'(A,I1)')'pT',ip
          write(*,'(A,I1,A)')'d#sigma-/dp_{T}(',ip,')--[pb/GeV]'
          write(*,'(A,I1,A)')'p_{T}(',ip,')--[GeV]'
          do i=1,ndiv_pT(ip)
            write(*,*)xpT(ip,i),fxpTtot(ip,i)
          end do
          write(*,*)'END'
        end if
  ! Plot distributions in eta
        if(o_eta(ip).eq.1)then
          sfxetatot(ip)=0d0
          do j=1,ndiv_eta(ip)
            fxetatot(ip,j)=0.d0
            do i=1,it
              fxeta(ip,j,i)=fxeta(ip,j,i)*avgi/cnorm(i)/etaw(ip)
              fxetatot(ip,j)=fxetatot(ip,j)+fxeta(ip,j,i)            
            end do
            sfxetatot(ip)=sfxetatot(ip)+fxetatot(ip,j)*etaw(ip)
          end do
          write(*,*)'DISTRIBUTION'
          write(*,'(A,I1)')'eta',ip
          write(*,'(A,I1,A)')'d#sigma-/d#eta(',ip,')--[pb]'
          write(*,'(A,I1,A)')'#eta(',ip,')'          
          do i=1,ndiv_eta(ip)
            write(*,*)xeta(ip,i),fxetatot(ip,i)
          end do
          write(*,*)'END'
        end if
  ! Plot distributions in phi        
        if(o_phi(ip).eq.1)then
          sfxphitot(ip)=0d0
          do j=1,ndiv_phi(ip)
            fxphitot(ip,j)=0.d0
            do i=1,it
              fxphi(ip,j,i)=fxphi(ip,j,i)*avgi/cnorm(i)/phiw(ip)
              fxphitot(ip,j)=fxphitot(ip,j)+fxphi(ip,j,i)            
            end do
            sfxphitot(ip)=sfxphitot(ip)+fxphitot(ip,j)*phiw(ip)
          end do
          write(*,*)'DISTRIBUTION'
          write(*,'(A,I1)')'phi',ip
          write(*,'(A,I1,A)')'d#sigma-/d#phi(',ip,')--[pb/rad]'
          write(*,'(A,I1,A)')'#phi(',ip,')--[rad]'          
          do i=1,ndiv_phi(ip)
            write(*,*)xphi(ip,i),fxphitot(ip,i)
          end do
          write(*,*)'END'
        end if
      end do
  
  ! Plot distribution in ETmiss
      if(o_ETmiss.eq.1)then
        sfxETmisstot=0d0
        do j=1,ndiv_ETmiss
          fxETmisstot(j)=0.d0
          do i=1,it
            fxETmiss(j,i)=fxETmiss(j,i)*avgi/cnorm(i)/ETmissw
            fxETmisstot(j)=fxETmisstot(j)+fxETmiss(j,i)            
          end do
          sfxETmisstot=sfxETmisstot+fxETmisstot(j)*ETmissw
        end do
        write(*,*)'DISTRIBUTION'
        write(*,*)'ETmiss'
        write(*,*)'d#sigma-/dp_{Tmiss}--[pb/GeV]'
        write(*,*)'p_{T}(miss)--[GeV]'
        do i=1,ndiv_ETmiss
          write(*,*)xETmiss(i),fxETmisstot(i)
        end do
        write(*,*)'END'
      end if

  ! Plot distribution in pT356
      if(o_pT356.eq.1)then
        sfxpT356tot=0d0
        do j=1,ndiv_pT356
          fxpT356tot(j)=0.d0
          do i=1,it
            fxpT356(j,i)=fxpT356(j,i)*avgi/cnorm(i)/pT356w
            fxpT356tot(j)=fxpT356tot(j)+fxpT356(j,i)            
          end do
          sfxpT356tot=sfxpT356tot+fxpT356tot(j)*pT356w
        end do
        write(*,*)'DISTRIBUTION'
        write(*,*)'pT356'
        write(*,*)'d#sigma-/dp_{T}--[pb/GeV]'
        write(*,*)'p_{T}(t)--[GeV]'
        do i=1,ndiv_pT356
          write(*,*)xpT356(i),fxpT356tot(i)
        end do
        write(*,*)'END'
      end if
  ! Plot distribution in eta356
      if(o_eta356.eq.1)then
        sfxeta356tot=0d0
        do j=1,ndiv_eta356
          fxeta356tot(j)=0.d0
          do i=1,it
            fxeta356(j,i)=fxeta356(j,i)*avgi/cnorm(i)/eta356w
            fxeta356tot(j)=fxeta356tot(j)+fxeta356(j,i)            
          end do
          sfxeta356tot=sfxeta356tot+fxeta356tot(j)*eta356w
        end do
        write(*,*)'DISTRIBUTION'
        write(*,*)'eta356'
        write(*,*)'d#sigma-/d#eta--[pb]'
        write(*,*)'#eta(#bar{t})'
        do i=1,ndiv_eta356
          write(*,*)xeta356(i),fxeta356tot(i)
        end do
        write(*,*)'END'
      end if  
  ! Plot distribution in phi356
      if(o_phi356.eq.1)then
        sfxphi356tot=0d0
        do j=1,ndiv_phi356
          fxphi356tot(j)=0.d0
          do i=1,it
            fxphi356(j,i)=fxphi356(j,i)*avgi/cnorm(i)/phi356w
            fxphi356tot(j)=fxphi356tot(j)+fxphi356(j,i)            
          end do
          sfxphi356tot=sfxphi356tot+fxphi356tot(j)*phi356w
        end do
        write(*,*)'DISTRIBUTION'
        write(*,*)'phi356'
        write(*,*)'d#sigma-/d#phi--[pb]'
        write(*,*)'#phi(#bar{t})'
        do i=1,ndiv_phi356
          write(*,*)xphi356(i),fxphi356tot(i)
        end do
        write(*,*)'END'
      end if       
  ! Plot distribution in pT478
      if(o_pT478.eq.1)then  
        sfxpT478tot=0d0
        do j=1,ndiv_pT478
          fxpT478tot(j)=0.d0
          do i=1,it
            fxpT478(j,i)=fxpT478(j,i)*avgi/cnorm(i)/pT478w
            fxpT478tot(j)=fxpT478tot(j)+fxpT478(j,i)            
          end do
          sfxpT478tot=sfxpT478tot+fxpT478tot(j)*pT478w
        end do
        write(*,*)'DISTRIBUTION'
        write(*,*)'pT478'
        write(*,*)'d#sigma-/dp_{T}--[pb/GeV]'
        write(*,*)'p_{T}(#bar{t})--[GeV]'
        do i=1,ndiv_pT478
          write(*,*)xpT478(i),fxpT478tot(i)
        end do
        write(*,*)'END'
      end if
  ! Plot distribution in eta478
      if(o_eta478.eq.1)then
        sfxeta478tot=0d0
        do j=1,ndiv_eta478
          fxeta478tot(j)=0.d0
          do i=1,it
            fxeta478(j,i)=fxeta478(j,i)*avgi/cnorm(i)/eta478w
            fxeta478tot(j)=fxeta478tot(j)+fxeta478(j,i)            
          end do
          sfxeta478tot=sfxeta478tot+fxeta478tot(j)*eta478w
        end do
        write(*,*)'DISTRIBUTION'
        write(*,*)'eta478'
        write(*,*)'d#sigma-/d#eta--[pb]'
        write(*,*)'#eta(#bar{t})'
        do i=1,ndiv_eta478
          write(*,*)xeta478(i),fxeta478tot(i)
        end do
        write(*,*)'END'
      end if      
  ! Plot distribution in phi478
      if(o_phi478.eq.1)then
        sfxphi478tot=0d0
        do j=1,ndiv_phi478
          fxphi478tot(j)=0.d0
          do i=1,it
            fxphi478(j,i)=fxphi478(j,i)*avgi/cnorm(i)/phi478w
            fxphi478tot(j)=fxphi478tot(j)+fxphi478(j,i)            
          end do
          sfxphi478tot=sfxphi478tot+fxphi478tot(j)*phi478w
        end do
        write(*,*)'DISTRIBUTION'
        write(*,*)'phi478'
        write(*,*)'d#sigma-/d#phi--[pb]'
        write(*,*)'#phi(#bar{t})'
        do i=1,ndiv_phi478
          write(*,*)xphi478(i),fxphi478tot(i)
        end do
        write(*,*)'END'
      end if    
  ! Plot distribution in Mtt
      if(o_rMtt.eq.1)then
        sfxrMtttot=0d0
        do j=1,ndiv_rMtt
          fxrMtttot(j)=0.d0
          do i=1,it
            fxrMtt(j,i)=fxrMtt(j,i)*avgi/cnorm(i)/rMttw
            fxrMtttot(j)=fxrMtttot(j)+fxrMtt(j,i)            
          end do
          sfxrMtttot=sfxrMtttot+fxrMtttot(j)*rMttw
        end do
        write(*,*)'DISTRIBUTION'
        write(*,*)'Mtt'
        write(*,*)'d#sigma-/dM_{tt}--[pb/GeV]'
        write(*,*)'M_{tt}--[GeV]'
        do i=1,ndiv_rMtt
          write(*,*)xrMtt(i),fxrMtttot(i)
        end do
        write(*,*)'END'
      end if
  ! Plot distribution in beta.
      if(o_beta.eq.1)then
        sfxbetatot=0d0
        do j=1,ndiv_beta
          fxbetatot(j)=0.d0
          do i=1,it
            fxbeta(j,i)=fxbeta(j,i)*avgi/cnorm(i)/betaw
            fxbetatot(j)=fxbetatot(j)+fxbeta(j,i)            
          end do
          sfxbetatot=sfxbetatot+fxbetatot(j)*betaw
        end do
        write(*,*)'DISTRIBUTION'
        write(*,*)'Beta'
        write(*,*)'d#sigma-/d#Beta_{t}--[pb]'
        write(*,*)'#Beta_{t}'
        do i=1,ndiv_beta
          write(*,*)xbeta(i),fxbetatot(i)
        end do
        write(*,*)'END'
      end if  
  ! Plot distribution in cost
      if(o_cost.eq.1)then
        sfxcosttot=0d0
        do j=1,ndiv_cost
          fxcosttot(j)=0.d0
          do i=1,it
            fxcost(j,i)=fxcost(j,i)*avgi/cnorm(i)/costw
            fxcosttot(j)=fxcosttot(j)+fxcost(j,i)            
          end do
          sfxcosttot=sfxcosttot+fxcosttot(j)*costw
        end do
        write(*,*)'DISTRIBUTION'
        write(*,*)'cost'
        write(*,*)'d#sigma-/dcos#theta--[pb]'
        write(*,*)'cos#theta'
        do i=1,ndiv_cost
          write(*,*)xcost(i),fxcosttot(i)
        end do
        write(*,*)'END'
      end if 
  ! Plot distribution in Et
      if(o_Et.eq.1)then

        sfxEttot=0d0
        do j=1,ndiv_Et
          fxEttot(j)=0.d0
          do i=1,it
            fxEt(j,i)=fxEt(j,i)*avgi/cnorm(i)/Etw
            fxEttot(j)=fxEttot(j)+fxEt(j,i)            
          end do
          sfxEttot=sfxEttot+fxEttot(j)*Etw
        end do
        write(*,*)'DISTRIBUTION'
        write(*,*)'Et'
        write(*,*)'d#sigma-/dE_{t}--[pb/GeV]'
        write(*,*)'E_{t}--[GeV]'
        do i=1,ndiv_Et
          write(*,*)xEt(i),fxEttot(i)
        end do
        write(*,*)'END'
      end if
  ! Plot distributions in all transverse variables
      do itrans=1,ntrans
        if(o_tran(itrans).eq.1)then
          sfxtranstot(itrans)=0d0
          do j=1,ndiv_trans(itrans)
            fxtranstot(itrans,j)=0.d0
            do i=1,it
              fxtrans(itrans,j,i)=fxtrans(itrans,j,i)
     &                           *avgi/cnorm(i)/transw(itrans)
              fxtranstot(itrans,j)=fxtranstot(itrans,j)
     &                            +fxtrans(itrans,j,i)            
            end do
            sfxtranstot(itrans)=sfxtranstot(itrans)+
     &                               fxtranstot(itrans,j)*transw(itrans)
          end do
          write(*,*)'DISTRIBUTION'
          if (itrans.eq.1)then
            write(*,*)'Mvis'
            write(*,*)'d#sigma-/dM_{vis}--[pb/GeV]'
            write(*,*)'M_{vis}--[GeV]'
          else if (itrans.eq.2)then
            write(*,*)'HT'
            write(*,*)'d#sigma-/dH_{T}--[pb/GeV]'
            write(*,*)'H_{T}--[GeV]'          
          else if (itrans.eq.3)then
            write(*,*)'M_T1'
            write(*,*)'d#sigma-/dM_{T1}--[pb/GeV]'
            write(*,*)'M_{T}--[GeV]'
          else if (itrans.eq.4)then
            write(*,*)'M_T2'
            write(*,*)'d#sigma-/dM_{T2}--[pb/GeV]'
            write(*,*)'M_{T}--[GeV]'
          else if (itrans.eq.5)then
            write(*,*)'M_T3'
            write(*,*)'d#sigma-/dM_{T3}--[pb/GeV]'
            write(*,*)'M_{T}--[GeV]'
          else if (itrans.eq.6)then
            write(*,*)'MlT'
            write(*,*)'d#sigma-/dM^{l}_{T}--[pb/GeV]'
            write(*,*)'M_{T}--[GeV]'
          else if (itrans.eq.7)then
            write(*,*)'M_CT1'
            write(*,*)'d#sigma-/dM_{T1}--[pb/GeV]'
            write(*,*)'M_{T}--[GeV]'
          else if (itrans.eq.8)then
            write(*,*)'M_CT2'
            write(*,*)'d#sigma-/dM_{T2}--[pb/GeV]'
            write(*,*)'M_{T}--[GeV]'
          else if (itrans.eq.9)then
            write(*,*)'M_CT3'
            write(*,*)'d#sigma-/dM_{T3}--[pb/GeV]'
            write(*,*)'M_{T}--[GeV]'
          else if (itrans.eq.10)then
            write(*,*)'MlCT'
            write(*,*)'d#sigma-/dM^{l}_{CT}--[pb/GeV]'
            write(*,*)'M^{l}_{CT}--[GeV]'
          else
            continue
          end if
          do i=1,ndiv_trans(itrans)
            write(*,*)xtrans(itrans,i),fxtranstot(itrans,i)
          end do
          write(*,*)'END'
        end if
      end do
  ! Plot distribution in fl
      if(o_fl.eq.1)then
        sfxfltot=0d0
        do j=1,ndiv_fl
          fxfltot(j)=0.d0
          do i=1,it
            fxfl(j,i)=fxfl(j,i)*avgi/cnorm(i)/flw
            fxfltot(j)=fxfltot(j)+fxfl(j,i)            
          end do
          sfxfltot=sfxfltot+fxfltot(j)*flw
        end do
        write(*,*)'DISTRIBUTION'
        write(*,*)'fl'
        write(*,*)'d#sigma-/d#phi_{l}--[pb]'
        write(*,*)'#phi_{l}--[-rad-]'
        do i=1,ndiv_fl
          write(*,*)xfl(i),fxfltot(i)
        end do
        write(*,*)'END'
      end if 
  ! Plot distribution in cosfl
      if(o_cosfl.eq.1)then
        sfxcosfltot=0d0
        do j=1,ndiv_cosfl
          fxcosfltot(j)=0.d0
          do i=1,it
            fxcosfl(j,i)=fxcosfl(j,i)*avgi/cnorm(i)/cosflw
            fxcosfltot(j)=fxcosfltot(j)+fxcosfl(j,i)            
          end do
          sfxcosfltot=sfxcosfltot+fxcosfltot(j)*cosflw
        end do
        write(*,*)'DISTRIBUTION'
        write(*,*)'cosfl'
        write(*,*)'d#sigma-/dcos#phi_{l}--[pb]'
        write(*,*)'cos#phi_{l}'
        do i=1,ndiv_cosfl
          write(*,*)xcosfl(i),fxcosfltot(i)
        end do
        write(*,*)'END'
      end if
  ! Plot distribution in delta phi
      if(o_dphi.eq.1)then
        sfxdphitot=0d0
        do j=1,ndiv_dphi
          fxdphitot(j)=0.d0
          do i=1,it
            fxdphi(j,i)=fxdphi(j,i)*avgi/cnorm(i)/dphiw
            fxdphitot(j)=fxdphitot(j)+fxdphi(j,i)            
          end do
          sfxdphitot=sfxdphitot+fxdphitot(j)*dphiw
        end do
        write(*,*)'DISTRIBUTION'
        write(*,*)'dphi'
        write(*,*)'d#sigma-/d#Delta#phi--[pb]'
        write(*,*)'#Delta#phi'
        do i=1,ndiv_dphi
          write(*,*)xdphi(i),fxdphitot(i)
        end do
        write(*,*)'END'
      end if
  ! Plot distribution in cost5
      if(o_cost5.eq.1)then
        sfxcost5tot=0d0
        do j=1,ndiv_cost5
          fxcost5tot(j)=0.d0
          do i=1,it
            fxcost5(j,i)=fxcost5(j,i)*avgi/cnorm(i)/cost5w
            fxcost5tot(j)=fxcost5tot(j)+fxcost5(j,i)            
          end do
          sfxcost5tot=sfxcost5tot+fxcost5tot(j)*cost5w
        end do
        write(*,*)'DISTRIBUTION'
        write(*,*)'cost5'
        write(*,*)'d#sigma-/dcos#theta_{+}--[pb]'
        write(*,*)'cos#theta_{+}'
        do i=1,ndiv_cost5
          write(*,*)xcost5(i),fxcost5tot(i)
        end do
        write(*,*)'END'
      end if
  ! Plot distribution in cost7
      if(o_cost7.eq.1)then
        sfxcost7tot=0d0
        do j=1,ndiv_cost7
          fxcost7tot(j)=0.d0
          do i=1,it
            fxcost7(j,i)=fxcost7(j,i)*avgi/cnorm(i)/cost7w
            fxcost7tot(j)=fxcost7tot(j)+fxcost7(j,i)            
          end do
          sfxcost7tot=sfxcost7tot+fxcost7tot(j)*cost7w
        end do
        write(*,*)'DISTRIBUTION'
        write(*,*)'cost7'
        write(*,*)'d#sigma-/dcos#theta_{-}--[pb]'
        write(*,*)'cos#theta_{-}'
        do i=1,ndiv_cost7
          write(*,*)xcost7(i),fxcost7tot(i)
        end do
        write(*,*)'END'
      end if
  ! Plot distribution in ct7ct5
      if(o_ct7ct5.eq.1)then
        sfxct7ct5tot=0d0
        do j=1,ndiv_ct7ct5
          fxct7ct5tot(j)=0.d0
          do i=1,it
            fxct7ct5(j,i)=fxct7ct5(j,i)*avgi/cnorm(i)/ct7ct5w
            fxct7ct5tot(j)=fxct7ct5tot(j)+fxct7ct5(j,i)            
          end do
          sfxct7ct5tot=sfxct7ct5tot+fxct7ct5tot(j)*ct7ct5w
        end do
        write(*,*)'DISTRIBUTION'
        write(*,*)'ct7ct5'
        write(*,*)
     &      'd^{2}#sigma-/d(cos#theta^{*}_{+}cos#theta^{*}_{-})--[pb]'
        write(*,*)'cos#theta_{+}cos#theta_{-}'
        do i=1,ndiv_ct7ct5
          write(*,*)xct7ct5(i),fxct7ct5tot(i)
        end do
        write(*,*)'END'
      end if
  ! Plot distributions in all asymmetries
      if((o_sigp.eq.1).and.(o_sigm.eq.1))then
        do jasy=1,nasym
          if(o_asym(jasy).eq.0)then
            continue
          else
            ! snorm(jasy)=0.d0
            sfxsigptot(jasy)=0d0
            do j=1,ndiv_sigp
              fxsigptot(jasy,j)=0.d0
              do i=1,it
                fxsigp(jasy,j,i)=fxsigp(jasy,j,i)*avgi/cnorm(i)/sigpw
                fxsigptot(jasy,j)=fxsigptot(jasy,j)+fxsigp(jasy,j,i)
              end do
              sfxsigptot(jasy)=sfxsigptot(jasy)+fxsigptot(jasy,j)*sigpw
            end do
            sfxsigmtot(jasy)=0d0           
            do j=1,ndiv_sigm
              do i=1,it
                fxsigm(jasy,j,i)=fxsigm(jasy,j,i)*avgi/cnorm(i)/sigmw
                fxsigmtot(jasy,j)=fxsigmtot(jasy,j)+fxsigm(jasy,j,i)
              end do
              sfxsigmtot(jasy)=sfxsigmtot(jasy)+fxsigmtot(jasy,j)*sigmw
            end do
            write(*,*)'ASYMMETRY'
            if(jasy.eq.1)then
                write(*,*)'ALL'
                write(*,*)'A_{LL}'
            else if(jasy.eq.2)then
                write(*,*)'AL'
                write(*,*)'A_{L}'
            else if(jasy.eq.3)then
                write(*,*)'APV'
                write(*,*)'A_{PV}'
            else if(jasy.eq.4)then
                write(*,*)'AFB'
                write(*,*)'A_{FB}'
            else if(jasy.eq.5)then
                write(*,*)'AFBst'
                write(*,*)'A_{FB^{*}}'
            else if(jasy.eq.6)then
                write(*,*)'AtRFB'
                write(*,*)'A^{t}_{RFB}'
            else if(jasy.eq.7)then
                write(*,*)'AttbRFB'
                write(*,*)"A^{b\bar{b}}_{RFB}"
            else if(jasy.eq.8)then
                write(*,*)'ARFB'
                write(*,*)"A_{RFB}"  
            else if(jasy.eq.9)then
                write(*,*)'A_l'
                write(*,*)'A_{l}'
            end if            
            write(*,*)'M_{tt}'
            ndiv_sig=(ndiv_sigm+ndiv_sigp)/2
            do i=1,ndiv_sig
              if(fxsigptot(jasy,i)+fxsigmtot(jasy,i).eq.0.d0)then
                write(*,*)(xsigm(i)+xsigp(i))/2.d0,0.d0
      !             snorm(jasy)=snorm(jasy)+0.d0
              else  
                write(*,*)(xsigm(i)+xsigp(i))/2.d0,
     &                   (fxsigptot(jasy,i)-fxsigmtot(jasy,i))/
     &                   (fxsigptot(jasy,i)+fxsigmtot(jasy,i)),
     &                    fxsigptot(jasy,i),fxsigmtot(jasy,i)
      !               snorm(jasy)=snorm(jasy)+
      !      &               (fxsigptot(jasy,i)-fxsigmtot(jasy,i))/
      !      &               (fxsigptot(jasy,i)+fxsigmtot(jasy,i))
      !      &               *fxrMtttot(i)*rMttw/avgi      
              end if
            end do
            asym_int(jasy)=(sfxsigptot(jasy)-sfxsigmtot(jasy))/
     &                     (sfxsigptot(jasy)+sfxsigmtot(jasy))          
            write(*,*)'END'
      !           write(*,*)'(Total Asymmetry:',asym_int(jasy),')'
      !           write(*,*)'(Integrated Asymmetry:',snorm(jasy),' )'
          end if
        end do
      end if
! ----------------------------------------------------------------------
! Check distributions
      diff_max=1E-12
      n_error=0
      do ip=3,ipmax
        if(o_pT(ip).eq.1)then
          if(abs(cross-sfxpTtot(ip))>diff_max)then
            write(*,*)'pT',ip,' Error:',sfxpTtot(ip)
            n_error=n_error+1
          end if
        end if
        if(o_eta(ip).eq.1)then
          if(abs(cross-sfxetatot(ip))>diff_max)then
            write(*,*)'eta',ip,' Error:',sfxetatot(ip)
            n_error=n_error+1
          end if
        end if
        if(o_phi(ip).eq.1)then
          if(abs(cross-sfxphitot(ip))>diff_max)then
            write(*,*)'phi',ip,' Error:',sfxphitot(ip)
            n_error=n_error+1
          end if
        end if
      end do      
      if(o_ETmiss.eq.1)then
        if(abs(cross-sfxETmisstot)>diff_max)then
          write(*,*)'ETmiss Error:',sfxETmisstot
          n_error=n_error+1
        end if
      end if
      if(o_pT356.eq.1)then
        if(abs(cross-sfxpT356tot)>diff_max)then
          write(*,*)'pT356 Error:',sfxpT356tot
          n_error=n_error+1
        end if
      end if
      if(o_eta356.eq.1)then
        if(abs(cross-sfxeta356tot)>diff_max)then
          write(*,*)'eta356 Error:',sfxeta356tot
          n_error=n_error+1
        end if
      end if
      if(o_phi356.eq.1)then
        if(abs(cross-sfxphi356tot)>diff_max)then
          write(*,*)'phi356 Error:',sfxphi356tot
          n_error=n_error+1
        end if
      end if
      if(o_pT478.eq.1)then
        if(abs(cross-sfxpT478tot)>diff_max)then
          write(*,*)'pT478 Error:',sfxpT478tot
          n_error=n_error+1
        end if
      end if
      if(o_eta478.eq.1)then
        if(abs(cross-sfxeta478tot)>diff_max)then
          write(*,*)'eta478 Error:',sfxeta478tot
          n_error=n_error+1
        end if
      end if
      if(o_phi478.eq.1)then
        if(abs(cross-sfxphi478tot)>diff_max)then
          write(*,*)'phi478 Error:',sfxphi478tot
          n_error=n_error+1
        end if
      end if
      if(o_rMtt.eq.1)then
        if(abs(cross-sfxrMtttot)>diff_max)then
          write(*,*)'rMtt Error:',sfxrMtttot
          n_error=n_error+1
        end if
      end if
      if(o_beta.eq.1)then
        if(abs(cross-sfxbetatot)>diff_max*10)then
          write(*,*)'beta Error:',sfxbetatot
          n_error=n_error+1
        end if
      end if
      if(o_cost.eq.1)then
        if(abs(cross-sfxcosttot)>diff_max)then
          write(*,*)'cost Error:',sfxcosttot
          n_error=n_error+1
        end if
      end if
      if(o_Et.eq.1)then
        if(abs(cross-sfxEttot)>diff_max)then
          write(*,*)'Et Error:',sfxEttot
          n_error=n_error+1
        end if
      end if
      do iasym=1,nasym
        if(o_asym(iasy).eq.0)then
          continue
        else
          if(abs(Atot(iasy)-asym_int(iasy))>diff_max)then
            write(*,*)'A Error:',iasy,asym_int(iasy)
            n_error=n_error+1
          end if
        end if
      end do
      if(o_fl.eq.1)then
        if(abs(cross-sfxfltot)>diff_max)then
          write(*,*)'fl Error:',sfxfltot
          n_error=n_error+1
        end if
      end if
      if(o_cosfl.eq.1)then
        if(abs(cross-sfxcosfltot)>diff_max)then
          write(*,*)'cosfl Error:',sfxcosfltot
          n_error=n_error+1
        end if
      end if
      if(o_dphi.eq.1)then
        if(abs(cross-sfxdphitot)>diff_max)then
          write(*,*)'dphi Error:',sfxdphitot
          n_error=n_error+1
        end if
      end if
      do iasy=1,nasym
        if(o_asym(iasy).eq.0)then
          continue
        else
          if(abs(Atot(iasy)-asym_int(iasy))>diff_max)then
            write(*,*)'A Error:',iasy,asym_int(iasy)
            n_error=n_error+1
          end if
        end if
      end do
      write(*,*)'INTEGRATION ERRORS:',n_error
! ----------------------------------------------------------------------
! End program
      write(*,*)'CLOSE'
      write(*,*)'======================================================'
      stop
      end
! ======================================================================