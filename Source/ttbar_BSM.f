! ======================================================================
      program ppbbllnn_BSM
! Authors: Declan Millar, Stefano Moretti.

! Calculates the cross section and generates distributions for
!   pp -> tt,
!   pp -> tt -> bW^+bbarW^- -> bbbare^+nue^-nubarc
! (Future:) pp -> bW^+bbarW^- -> b bbar e^+ nu qqbar'
! Uses adapted Madgraph functions.
! Uses cteq6 and mrs99 PDF subroutines.
! ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)

!  LHE format
!   Initialization information
      integer maxpup
      parameter (maxpup=100)
      integer idbmup,pdfgup,pdfsup,idwtup,nprup,lprup 
      double precision ebmup,xsecup,xerrup,xmaxup 
      common/heprup/idbmup(2),ebmup(2),pdfgup(2),pdfsup(2),
     & idwtup,nprup,xsecup(maxpup),xerrup(maxpup),
     & xmaxup(maxpup),lprup(maxpup)
!   Information on each separate event
      integer maxnup
      parameter (maxnup=500)
      integer nup,idprup,idup,istup,mothup,icolup
      double precision xwgtup,scalup,aqedup,aqcdup,pup,vtimup,spinup
      common/hepeup/nup,idprup,xwgtup,scalup,aqedup,aqcdup, 
     & idup(maxnup),istup(maxnup),mothup(2,maxnup),icolup(2,maxnup),
     & pup(5,maxnup),vtimup(maxnup),spinup(maxnup)

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
      common/final/ifinal
      common/top/rmt,gamt
      common/W/rmW,gamW
      common/Z/rmZ,gamZ
      common/H/rmH,gamH      
      common/stat/npoints
      common/coll/ecm_coll
      common/cuts/ytcut,yttcut
!   Polarised/Spatial cross sections
      common/polarised/polcross(20,-1:1,-1:1),polerror(20,-1:1,-1:1)
      common/spatial/asycross(5,20,-1:1),asyerror(5,20,-1:1) !nasy -3
!   Permitted gauge sectors
      common/igauge/iQCD,iEW,iBSM
!   Interference       
      common/interference/iint      
!   Z' masses and VA/LR couplings
      common/Zp/rmZp(5),gamZp(5)
      common/Zpparam/paramZp(5)
      common/ZpAVcoup/gp(5),gV_d(5),gA_d(5),gV_u(5),gA_u(5)
      common/ZpLRcoup/gZpd(2,5),gZpu(2,5)
!   Narrow width approximation (NWA)
      common/NWA/gNWA
!   Structure functions
      common/partdist/istructure
      common/ALFASTRONG/rlambdaQCD4,nloops
      common/collider/icoll
!   Switch for all distributions
      common/distros/idist
!   Distributions in pTs of external particles
      common/ext_pT/pTmax(8),pTmin(8),pTw(8)
      common/dist_pT/xpT(8,500),fxpT(8,500,20),fxpTtot(8,500)
      common/inp_pT/m_pT(8)
      common/div_pT/ndiv_pT(8)
!   Distributions in phis of external particles
      common/ext_phi/phimax(8),phimin(8),phiw(8)
      common/dist_phi/xphi(8,500),fxphi(8,500,20),fxphitot(8,500)
      common/inp_phi/m_phi(8)
      common/div_phi/ndiv_phi(8)
!   Distribution in pT of the top
      common/ext_pT356/pT356max,pT356min,pT356w
      common/dist_pT356/xpT356(500),fxpT356(500,20),fxpT356tot(500)
      common/inp_pT356/m_pT356
      common/div_pT356/ndiv_pT356
!   Distribution in pT of the anti-top
      common/ext_pT478/pT478max,pT478min,pT478w
      common/dist_pT478/xpT478(500),fxpT478(500,20),fxpT478tot(500)
      common/inp_pT478/m_pT478
      common/div_pT478/ndiv_pT478
!   Distribution in ETmiss
      common/ext_ETmiss/ETmissmax,ETmissmin,ETmissw
      common/dist_ETmiss/xETmiss(500),fxETmiss(500,20),fxETmisstot(500)
      common/inp_ETmiss/m_ETmiss
      common/div_ETmiss/ndiv_ETmiss
!   Distribution in eta of the top / bottom
      common/ext_eta3/eta3max,eta3min,eta3w
      common/dist_eta3/xeta3(500),fxeta3(500,20),fxeta3tot(500)
      common/inp_eta3/m_eta3
      common/div_eta3/ndiv_eta3
!   Distribution in eta of the anti-top / bottom
      common/ext_eta4/eta4max,eta4min,eta4w
      common/dist_eta4/xeta4(500),fxeta4(500,20),fxeta4tot(500)
      common/inp_eta4/m_eta4
      common/div_eta4/ndiv_eta4
!   Distribution in eta of the anti-leetaon
      common/ext_eta5/eta5max,eta5min,eta5w
      common/dist_eta5/xeta5(500),fxeta5(500,20),fxeta5tot(500)
      common/inp_eta5/m_eta5
      common/div_eta5/ndiv_eta5
!   Distribution in eta of the neutrino
      common/ext_eta6/eta6max,eta6min,eta6w
      common/dist_eta6/xeta6(500),fxeta6(500,20),fxeta6tot(500)
      common/inp_eta6/m_eta6
      common/div_eta6/ndiv_eta6
!   Distribution in eta of the leetaon
      common/ext_eta7/eta7max,eta7min,eta7w
      common/dist_eta7/xeta7(500),fxeta7(500,20),fxeta7tot(500)
      common/inp_eta7/m_eta7
      common/div_eta7/ndiv_eta7
!   Distribution in eta of the anti-neutrino
      common/ext_eta8/eta8max,eta8min,eta8w
      common/dist_eta8/xeta8(500),fxeta8(500,20),fxeta8tot(500)
      common/inp_eta8/m_eta8
      common/div_eta8/ndiv_eta8
!   Distribution in eta of the top
      common/ext_eta356/eta356max,eta356min,eta356w
      common/dist_eta356/xeta356(500),fxeta356(500,20),fxeta356tot(500)
      common/inp_eta356/m_eta356
      common/div_eta356/ndiv_eta356
!   Distribution in eta of the anti-top
      common/ext_eta478/eta478max,eta478min,eta478w
      common/dist_eta478/xeta478(500),fxeta478(500,20),fxeta478tot(500)
      common/inp_eta478/m_eta478
      common/div_eta478/ndiv_eta478
!   Distribution in invarient mass of the top pair
      common/ext_rmass/rmassmax,rmassmin,rmassw
      common/dist_rmass/xrmass(500),fxrmass(500,20),fxrmasstot(500)
      common/inp_rmass/m_rmass
      common/div_rmass/ndiv_rmass
!   Distribution in invarient mass of all visible decay products
      common/ext_rmvis/rmvismax,rmvismin,rmvisw
      common/dist_rmvis/xrmvis(500),fxrmvis(500,20),fxrmvistot(500)
      common/inp_rmvis/m_rmvis
      common/div_rmvis/ndiv_rmvis  
!   Distribution in  boost of top pair centre of mass frame
      common/ext_beta/betamax,betamin,betaw
      common/dist_beta/xbeta(500),fxbeta(500,20),fxbetatot(500)
      common/inp_beta/m_beta
      common/div_beta/ndiv_beta
!   Distribution in cos(theta_t)
      common/ext_cost/costmax,costmin,costw
      common/dist_cost/xcost(500),fxcost(500,20),fxcosttot(500)
      common/inp_cost/m_cost
      common/div_cost/ndiv_cost
!   Distribution in top energy
      common/ext_Et/Etmax,Etmin,Etw
      common/dist_Et/xEt(500),fxEt(500,20),fxEttot(500)
      common/inp_Et/m_Et
      common/div_Et/ndiv_Et
!   Distribution in sum of ET
      common/ext_HT/HTmax,HTmin,HTw
      common/dist_HT/xHT(500),fxHT(500,20),fxHTtot(500)
      common/inp_HT/m_HT
      common/div_HT/ndiv_HT
!   Distribution in transverse mass
      common/ext_rM_T/rM_Tmax,rM_Tmin,rM_Tw
      common/dist_rM_T/xrM_T(500),fxrM_T(500,20),fxrM_Ttot(500)
      common/inp_rM_T/m_rM_T
      common/div_rM_T/ndiv_rM_T
!   Distribution in *total* contransverse mass
      common/ext_rM_CT/rM_CTmax,rM_CTmin,rM_CTw
      common/dist_rM_CT/xrM_CT(500),fxrM_CT(500,20),fxrM_CTtot(500)
      common/inp_rM_CT/m_rM_CT
      common/div_rM_CT/ndiv_rM_CT
!   Distribution in lepton contransverse mass
      common/ext_rMlCT/rMlCTmax,rMlCTmin,rMlCTw
      common/dist_rMlCT/xrMlCT(500),fxrMlCT(500,20),fxrMlCTtot(500)
      common/inp_rMlCT/m_rMlCT
      common/div_rMlCT/ndiv_rMlCT
!   Distribution in phi_l (lepton azimuthal angle)
      common/ext_fl/flmax,flmin,flw
      common/dist_fl/xfl(500),fxfl(500,20),fxfltot(500)
      common/inp_fl/m_fl
      common/div_fl/ndiv_fl
!   Distribution in cos_phi_l
      common/ext_cosfl/cosflmax,cosflmin,cosflw
      common/dist_cosfl/xcosfl(500),fxcosfl(500,20),fxcosfltot(500)
      common/inp_cosfl/m_cosfl
      common/div_cosfl/ndiv_cosfl
!   Distribution in sigp
      common/ext_sigp/sigpmax,sigpmin,sigpw
      common/dist_sigp/xsigp(1000),fxsigp(8,1000,20),fxsigptot(8,1000)
      common/inp_sigp/m_sigp
      common/div_sigp/ndiv_sigp
!   Distribution in sigm
      common/ext_sigm/sigmmax,sigmmin,sigmw
      common/dist_sigm/xsigm(1000),fxsigm(8,1000,20),fxsigmtot(8,1000)
      common/inp_sigm/m_sigm
      common/div_sigm/ndiv_sigm
!   Distributions in transverse variables
      common/inp_tran/m_tran(8) ! nasy
!   Distributions in asymmetries
      common/inp_asym/m_asy(8) ! nasy

! Local variables
!   Flag for Zp width specification
      dimension iwidth(5)
!   Model name
      character*50 model
!   Polarised/hemispherised cross sections
      dimension cnorm(20)
      dimension snorm(6) !,ave(4) 
      dimension poltot(-1:1,-1:1),polchi(-1:1,-1:1)
      dimension asytot(5,-1:1),asychi(5,-1:1)
      dimension sfxpTtot(8)
      dimension sfxsigptot(8),sfxsigmtot(8)
      dimension asym_int(8)
      dimension Atot(8),Atoterr(8)

! Local constants
!   pi
      parameter (pi=3.14159265358979323846d0)
!   Unit conversion GeV -> nb
      parameter (conv=0.38937966d9)
!   Date and time
      integer today(3), now(3)      
!   Quark masses
      data rmu/0.00d0/,
     &     rmd/0.00d0/,
     &     rms/0.00d0/,
     &     rmc/0.00d0/,
     &     rmb/4.25d0/
!   Higgs/gauge boson masses and widths
      data              gamW/2.08d0/
      data rmZ/91.19d0/,gamZ/2.50d0/
      data rmH/125.0d0/,gamH/0.31278d-2/
!   Branching ratio for t->bln      
      real*8 BRtbln/0.10779733d0/
!       real*8 BRtbln/1d0/

! External procedures
      external fxn
! ----------------------------------------------------------------------
! Read config file
!   Permitted gauge sector flags
      read(5,*) iQCD
      read(5,*) iEW
      read(5,*) iBSM
!   Interference flag
      read(5,*) iint
!   Final state flag (ifinal=0: no top decay; ifinal=1: dileptonic top decay)      
      read(5,*) ifinal
!   NWA flag (iNWA = 0: Actual top widths; iNWA = 1: tops in NWA)
      read(5,*) iNWA
!   Name of model file
      read(5,*) model
!   Collider energy     
      read(5,*) ecm_coll
!   Collider flag (icoll=0: pp; icoll=1: ppbar)
      read(5,*) icoll
!   PDFs
      read(5,*) istructure
!   Cut on top rapidity 
      read(5,*) ytcut
!   Cut on top pair boost
      read(5,*) yttcut
!   Number of Vegas calls per iteration
      read(5,*) ncall
!   Maximum number of Vegas iterations
      read(5,*) itmx
!   Desired Accuracy (If negative, run maximum iterations.)
      read(5,*) acc
!   Random number seed
      read(5,*) iseed
!   Distributions flag
      read(5,*) idist
!   Outout in lhe format
      read(5,*) ilhe
! !   Manually sum over costheta
!       read(5,*) isycost

! Modify config
!   NWA only for six-body final state.
      if(ifinal.eq.0) iNWA=0
!   itmx no more than 20.      
      if(itmx.gt.20)then
        write(*,*)'itmx does not have to exceed 20!'
        stop
      end if
!   No event weighting for *true* event generation
      if(ilhe.eq.1)then
        itmx=1
      end if
!   Extract model filename (Remove white space.)
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
          iwidth(i) = 0
        else
          iwidth(i) = 1
        end if
      enddo

! Distributions Setup
! (Set flags, binning range and divisions.)
! pT distributions
      do i=1,8
        m_pT(i)=idist
        pTmax(i)=7000.d0/(1+icoll*6)
        pTmin(i)=0.d0
        ndiv_pT(i)=175
! phi distributions
        m_phi(i)=idist
        phimax(i)=2*pi
        phimin(i)=0.d0
        ndiv_phi(i)=100
      end do  
!   top transverse momentum
      m_pT356=idist
      pT356max=7000.d0/(1+icoll*6)
      pT356min=0.d0
      ndiv_pT356=175
!   anti-top transverse momentum
      m_pT478=idist
      pT478max=7000.d0/(1+icoll*6)
      pT478min=0.d0
      ndiv_pT478=175
!   missing transverse momentum
      m_ETmiss=idist
      ETmissmax=7000.d0/(1+icoll*6)
      ETmissmin=0.d0
      ndiv_ETmiss=175
!   top/bottom pseudorapidity
      m_eta3=idist
      eta3max=+10
      eta3min=-10
      ndiv_eta3=50
!   anti-top/bottom pseudorapidity
      m_eta4=idist
      eta4max=+10
      eta4min=-10
      ndiv_eta4=50
!   anti-leetaon pseudorapidity
      m_eta5=idist
      eta5max=+10
      eta5min=-10
      ndiv_eta5=50
!   neutrino pseudorapidity
      m_eta6=idist
      eta6max=+10
      eta6min=-10
      ndiv_eta6=50
!   lepton pseudorapidity
      m_eta7=idist
      eta7max=+10
      eta7min=-10
      ndiv_eta7=50
!   neutrino pseudorapidity
      m_eta8=idist
      eta8max=+10
      eta8min=-10
      ndiv_eta8=50
!   top pseudorapidity
      m_eta356=idist
      eta356max=+10
      eta356min=-10
      ndiv_eta356=50
!   anti-top pseudorapidity
      m_eta478=idist
      eta478max=+10
      eta478min=-10
      ndiv_eta478=50
!   invarient mass of tt pair
      m_rmass=idist
      rmassmax=14000.d0/(1+icoll*6)
      rmassmin=0.d0
      ndiv_rmass=500
!   invarient mass of the visible decay products of the tt pair
      m_rmvis=idist
      rmvismax=14000.d0/(1+icoll*6)
      rmvismin=0.d0
      ndiv_rmvis=500
!   pseudorapidity
      m_eta=idist
      etamax=1000.d0
      etamin=0.d0
      ndiv_eta=100
!   boost of parton CoM
      m_beta=idist
      betamax=1000.d0
      betamin=0.d0
      ndiv_beta=100
!   costheta
      m_cost=idist
      costmax=+1.d0
      costmin=-1.d0
      ndiv_cost=50
!   top energy
      m_Et=idist
      Etmax=7000.d0/(1+icoll*6)
      Etmin=0.d0
      ndiv_Et=175
!   sum of tranvserse energy
      m_HT=idist
      HTmax=7000.d0/(1+icoll*6)
      HTmin=0.d0
      ndiv_HT=175   
!   transverse mass
      m_rM_T=idist
      rM_Tmax=14000.d0/(1+icoll*6)
      rM_Tmin=0.d0
      ndiv_rM_T=175
!   contransverse mass 1
      m_rM_CT=idist
      rM_CTmax=14000.d0/(1+icoll*6)
      rM_CTmin=0.d0
      ndiv_rM_CT=175    
!   contransverse mass 2
      m_rMlCT=idist
      rMlCTmax=500.d0/(1+icoll*6)
      rMlCTmin=0.d0
      ndiv_rMlCT=175      

!   phi_l
      m_fl=idist
      flmax=+2*pi
      flmin=0
      ndiv_fl=100
!   cosphi_l
      m_cosfl=idist
      cosflmax=+1.d0
      cosflmin=-1.d0
      ndiv_cosfl=100
!   sigp
      m_sigp=m_rmass
      sigpmax=rmassmax
      sigpmin=rmassmin
      ndiv_sigp=ndiv_rmass/10
!   sigm
      m_sigm=m_rmass
      sigmmax=rmassmax
      sigmmin=rmassmin
      ndiv_sigm=ndiv_rmass/10      

!   asymmetries
      do i_asym=1,8 ! N_asym
        m_asy(i_asym)=idist
      end do

!   Turn off 2->6 only distributions
      if (ifinal.eq.0)then
        do i=5,8
          m_pT(i)   = 0
          m_phi(i)  = 0
        end do
        m_eta5   = 0
        m_eta6   = 0
        m_eta7   = 0
        m_eta8   = 0
        m_pT356  = 0
        m_eta356 = 0
        m_pT478  = 0
        m_eta478 = 0
        m_ETmiss = 0
        m_HT     = 0
        m_rM_T   = 0
        m_rM_CT  = 0
        m_rMlCT  = 0
        m_rMvis  = 0
        m_fl     = 0
        m_cosfl  = 0
        m_asy(8) = 0
      end if
!   Turn off 2->2 only distributions
      if (ifinal.eq.1)then
        m_asy(1) = 0
        m_asy(2) = 0
        m_asy(3) = 0
      end if    

! Initialize MadGraph for MEs
      rmt=175.d0
      gamt=1.55d0
      gNWA=gamt
      if(iNWA.eq.1)gamt=1.d-5
      call initialize(rmt,gamt)

! Collider CM energy squared.      
      s=ecm_coll*ecm_coll

! QCDL4 is QCD LAMBDA4 (to match PDF fits).
! (PDFs are intrinsically linked to the value of lamda_QCD; alpha_QCD) 
      if(istructure.eq.1)qcdl4=0.326d0
      if(istructure.eq.2)qcdl4=0.326d0
      if(istructure.eq.3)qcdl4=0.326d0
      if(istructure.eq.4)qcdl4=0.215d0
      if(istructure.eq.5)qcdl4=0.300d0
      if(istructure.eq.6)qcdl4=0.300d0
      if(istructure.eq.7)qcdl4=0.300d0
      if(istructure.eq.8)qcdl4=0.229d0
      if(istructure.eq.9)qcdl4=0.383d0
      rlambdaQCD4=QCDL4

! Initialise CTEQ grids.
      if(istructure.le.4)then
        icteq=istructure
        call SetCtq6(ICTEQ)
      end if

! Use appropriately evolved alphas.
      if(istructure.eq.1)nloops=2
      if(istructure.eq.2)nloops=2
      if(istructure.eq.3)nloops=1
      if(istructure.eq.4)nloops=1
      if(istructure.eq.5)nloops=1
      if(istructure.eq.6)nloops=1
      if(istructure.eq.7)nloops=1
      if(istructure.eq.8)nloops=1
      if(istructure.eq.9)nloops=1

! Factor outside integration
!   Conversion GeV^-2 -> pb
      fac=conv
!   Azimuthal angle integrated out (No initial transverse polarisation.)
      fac=fac*2.d0*pi

! some EW parameters.
      a_em=1.d0/128.d0
      s2w=.2320d0
      rmW=rmZ*sqrt(1.d0-s2w)
      sw=sqrt(s2w)
      c2w=1.d0-s2w
      cw=sqrt(c2w)

! Calculate Zp couplings
      call coupZp

! Calculate sequential Zp widths
      do i=1,5
        if (iwidth(i).eq.0) gamZp(i)=
     &               widthZp(rmW,rmZ,rmZp(i),a_em,s2w,rlambdaQCD4,nloop)
      end do

! VEGAS parameters
      if(ifinal.eq.0)then
        ndim=3
      else if(ifinal.eq.1)then
        ndim=15
      end if
!   (If nprn<0 no print-out.)
      nprn=0
      if(ifinal.eq.0)then
!   Final state masses     
        rm3=rmt
        rm4=rmt
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
!         if(isycost.eq.1)then ! might not work
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

      else if(ifinal.eq.1)then
!   Final state masses
        rm3=rmb
        rm4=rmb
        rm5=0.d0
        rm6=0.d0
        rm7=0.d0
        rm8=0.d0
!   Integrates on:
!   x(9)=cos(theta_cm_356)=-cos(theta_cm_478),
!   x(15)=(x1-tau)/(1-tau),
!   x(14)=(ecm-rm3-rm4-rm5-rm6-rm7-rm8)
!        /(ecm_max-rm3-rm4-rm5-rm6-rm7-rm8),
!   x(13)=(XX356-XX356min)/(XX356max-XX356min),
!   where XX356=arctg((rm356**2-rmt**2)/rmt/gamt),
!   x(12)=(XX478-XX478min)/(XX478max-XX478min),
!   where XX478=arctg((rm478**2-rmt**2)/rmt/gamt),
!   x(11)=(XX56-XX56min)/(XX56max-XX56min),
!   where XX56=arctg((rm56**2-rmW**2)/rmW/gamW),
!   x(10)=(XX78-XX78min)/(XX78max-XX78min),
!   where XX78=arctg((rm78**2-rmW**2)/rmW/gamW),
!   
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
! Generate bins
! (Finds bin width, finds midpoints.)

      do i=1,8
        if(m_pT(i).eq.1)then
          pTw(i)=(pTmax(i)-pTmin(i))/ndiv_pT(i)
          do j=1,ndiv_pT(i)
            xpT(i,j)=pTmin(i)+pTw(i)*(j-1)+pTw(i)/2.d0
          end do
        end if
      end do

      if(m_pT356.eq.1)then
        pT356w=(pT356max-pT356min)/ndiv_pT356
        do i=1,ndiv_pT356
          xpT356(i)=pT356min+pT356w*(i-1)+pT356w/2.d0
        end do
      end if

      if(m_pT478.eq.1)then
        pT478w=(pT478max-pT478min)/ndiv_pT478
        do i=1,ndiv_pT478
          xpT478(i)=pT478min+pT478w*(i-1)+pT478w/2.d0
        end do
      end if

      if(m_ETmiss.eq.1)then
        ETmissw=(ETmissmax-ETmissmin)/ndiv_ETmiss
        do i=1,ndiv_ETmiss
          xETmiss(i)=ETmissmin+ETmissw*(i-1)+ETmissw/2.d0
        end do
      end if

      if(m_eta3.eq.1)then
        eta3w=(eta3max-eta3min)/ndiv_eta3
        do i=1,ndiv_eta3
          xeta3(i)=eta3min+eta3w*(i-1)+eta3w/2.d0
        end do
      end if

      if(m_eta4.eq.1)then
        eta4w=(eta4max-eta4min)/ndiv_eta4
        do i=1,ndiv_eta4
          xeta4(i)=eta4min+eta4w*(i-1)+eta4w/2.d0
        end do
      end if

      if(m_eta5.eq.1)then
        eta5w=(eta5max-eta5min)/ndiv_eta5
        do i=1,ndiv_eta5
          xeta5(i)=eta5min+eta5w*(i-1)+eta5w/2.d0
        end do
      end if

      if(m_eta6.eq.1)then
        eta6w=(eta6max-eta6min)/ndiv_eta6
        do i=1,ndiv_eta6
          xeta6(i)=eta6min+eta6w*(i-1)+eta6w/2.d0
        end do
      end if

      if(m_eta7.eq.1)then
        eta7w=(eta7max-eta7min)/ndiv_eta7
        do i=1,ndiv_eta7
          xeta7(i)=eta7min+eta7w*(i-1)+eta7w/2.d0
        end do
      end if

      if(m_eta8.eq.1)then
        eta8w=(eta8max-eta8min)/ndiv_eta8
        do i=1,ndiv_eta8
          xeta8(i)=eta8min+eta8w*(i-1)+eta7w/2.d0
        end do
      end if

      if(m_eta356.eq.1)then
        eta356w=(eta356max-eta356min)/ndiv_eta356
        do i=1,ndiv_eta356
          xeta356(i)=eta356min+eta356w*(i-1)+eta356w/2.d0
        end do
      end if

      if(m_eta478.eq.1)then
        eta478w=(eta478max-eta478min)/ndiv_eta478
        do i=1,ndiv_eta478
          xeta478(i)=eta478min+eta478w*(i-1)+eta478w/2.d0
        end do
      end if

      if(m_rmass.eq.1)then
        rmassw=(rmassmax-rmassmin)/ndiv_rmass
        do i=1,ndiv_rmass
          xrmass(i)=rmassmin+rmassw*(i-1)+rmassw/2.d0
        end do
      end if

      if(m_rMvis.eq.1)then
        rMvisw=(rMvismax-rMvismin)/ndiv_rMvis
        do i=1,ndiv_rMvis
          xrMvis(i)=rMvismin+rMvisw*(i-1)+rMvisw/2.d0
        end do
      end if

      if(m_beta.eq.1)then
        betaw=(betamax-betamin)/ndiv_beta
        do i=1,ndiv_beta
          xbeta(i)=betamin+betaw*(i-1)+betaw/2.d0
        end do
      end if
      if(m_cost.eq.1)then
        costw=(costmax-costmin)/ndiv_cost
        do i=1,ndiv_cost
          xcost(i)=costmin+costw*(i-1)+costw/2.d0
        end do
      end if

      if(m_Et.eq.1)then
        Etw=(Etmax-Etmin)/ndiv_Et
        do i=1,ndiv_Et
          xEt(i)=Etmin+Etw*(i-1)+Etw/2.d0
        end do
      end if

      if(m_HT.eq.1)then
        HTw=(HTmax-HTmin)/ndiv_HT
        do i=1,ndiv_HT
          xHT(i)=HTmin+HTw*(i-1)+HTw/2.d0
        end do
      end if

      if(m_rM_T.eq.1)then
        rM_Tw=(rM_Tmax-rM_Tmin)/ndiv_rM_T
        do i=1,ndiv_rM_T
          xrM_T(i)=rM_Tmin+rM_Tw*(i-1)+rM_Tw/2.d0
        end do
      end if

      if(m_rM_CT.eq.1)then
        rM_CTw=(rM_CTmax-rM_CTmin)/ndiv_rM_CT
        do i=1,ndiv_rM_CT
          xrM_CT(i)=rM_CTmin+rM_CTw*(i-1)+rM_CTw/2.d0
        end do
      end if

      if(m_rMlCT.eq.1)then
        rMlCTw=(rMlCTmax-rMlCTmin)/ndiv_rMlCT
        do i=1,ndiv_rMlCT
          xrMlCT(i)=rMlCTmin+rMlCTw*(i-1)+rMlCTw/2.d0
        end do
      end if

      if(m_fl.eq.1)then
        flw=(flmax-flmin)/ndiv_fl
        do i=1,ndiv_fl
          xfl(i)=flmin+flw*(i-1)+flw/2.d0
        end do
      end if

      if(m_cosfl.eq.1)then
        cosflw=(cosflmax-cosflmin)/ndiv_cosfl
        do i=1,ndiv_cosfl
          xcosfl(i)=cosflmin+cosflw*(i-1)+cosflw/2.d0
        end do
      end if

      if(m_sigp.eq.1)then
        sigpw=(sigpmax-sigpmin)/ndiv_sigp
        do i=1,ndiv_sigp
          xsigp(i)=sigpmin+sigpw*(i-1)+sigpw/2.d0
        end do
      end if

      if(m_sigm.eq.1)then
        sigmw=(sigmmax-sigmmin)/ndiv_sigm
        do i=1,ndiv_sigm
          xsigm(i)=sigmmin+sigmw*(i-1)+sigmw/2.d0
        end do
      end if
! Output information before integration
      write(*,*)'====================================================='
      call idate(today)   ! today(1)=day, (2)=month, (3)=year
      call itime(now)     ! now(1)=hour, (2)=minute, (3)=second
      write(*,*)'DATE ',today(3),today(2),today(1)
      write(*,*)'TIME ',now(1),now(2),now(3)
      write(*,*)'-----------------------------------------------------'
      write(*,*)'PROCESS'
      if(icoll.eq.0)then
        if(ifinal.eq.0)
     &    write(*,*)'pp #rightarrow t#bar{t}',
     &               ' #times BR(t#rightarrow bl#nu)^{2}'
        if(ifinal.eq.1)
     &    write(*,*)'pp #rightarrow t#bar{t}',
     &               '#rightarrow b#bar{b} W^{+}W^{-}',
     &               '#rightarrow b#bar{b} l^{+}l^{-} #nu#bar{#nu}'
      else if(icoll.eq.1)then
        if(ifinal.eq.0)
     &    write(*,*)'p#bar{p} #rightarrow t#bar{t}',
     &               ' #times BR(t#rightarrow bl#nu)^{2}'
        if(ifinal.eq.1) 
     &    write(*,*)'p#bar{p} #rightarrow t#bar{t}',
     &               '#rightarrow b#bar{b} W^{+}W^{-}',
     &               '#rightarrow b#bar{b} l^{+}l^{-} #nu#bar{#nu}'
      end if
      write(*,*)'-----------------------------------------------------'
      write(*,*)'NOTES'            
      write(*,*)'Units: GeV'
      write(*,*)'Quarks: all massless except t, b.'           
      if(istructure.eq.1)write(*,*)'PDFs: cteq6m.'
      if(istructure.eq.2)write(*,*)'PDFs: cteq6d.'
      if(istructure.eq.3)write(*,*)'PDFs: cteq6l.'
      if(istructure.eq.4)write(*,*)'PDFs: cteq6l1.'
      if(istructure.eq.5)write(*,*)'PDFs: mrs99 (cor01).'
      if(istructure.eq.6)write(*,*)'PDFs: mrs99 (cor02).'
      if(istructure.eq.7)write(*,*)'PDFs: mrs99 (cor03).'
      if(istructure.eq.8)write(*,*)'PDFs: mrs99 (cor04).'
      if(istructure.eq.9)write(*,*)'PDFs: mrs99 (cor05).'
      if((ifinal.eq.1).and.(iNWA.eq.0))write(*,*)'Tops: off-shell.'
      if((ifinal.eq.1).and.(iNWA.eq.1))write(*,*)'Tops: NWA.'
      write(*,*)'BSM model: ',model
      if(iQCD.eq.1)write(*,*)'QCD: On '
      if(iQCD.eq.0)write(*,*)'QCD: Off'
      if(iEW.eq.1) write(*,*)'EW:  On '
      if(iEW.eq.0) write(*,*)'EW:  Off'
      if(iBSM.eq.1)write(*,*)'BSM: On '
      if(iBSM.eq.0)write(*,*)'BSM: Off'   
      if(iint.eq.0)write(*,*)'Interference: None'
      if(iint.eq.1)write(*,*)'Interference: SM'
      if(iint.eq.2)write(*,*)'Interference: Full'
      if(iint.eq.3)write(*,*)'Interference: No square terms.'      
      write(*,*)'-----------------------------------------------------'      
      write(*,*)'PARAMETERS'
      write(*,*)'#sqrt{s}              ',ecm_coll
      write(*,*)'at |y| <              ',abs(ytcut)
      write(*,*)'Loops a_s evaluated at',nloops
      write(*,*)'a_{s}(M_{Z})          ',alfas(rmZ,rLambdaQCD4,nloops)
      write(*,*)'#Lambda_{QCD}(4)      ',QCDL4
      write(*,*)'m_{b}                 ',rmb    
      write(*,*)'#Gamma_{b}            ',0.d0   
      write(*,*)'m_{t}                 ',rmt    
      write(*,*)'#Gamma_{t}            ',gamt   
      write(*,*)'m_{Z}                 ',rmZ    
      write(*,*)'#Gamma_{Z}            ',gamZ   
      write(*,*)'m_{W}                 ',rmW    
      write(*,*)'#Gamma_{W}            ',gamW   
      write(*,*)'m_{H}                 ',rmH    
      write(*,*)'#Gamma_{H}            ',gamH
      write(*,*)'-----------------------------------------------------' 
      write(*,*)'ZPRIME PARAMETERS'
      do i=1,5
        if(rmZp(i).gt.0)then
          write(*,*)'Z#prime',i 
          write(*,*)'m_{Z#prime}           ',rmZp(i) 
          write(*,*)'iwidth:               ',iwidth(i)   
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
      write(*,*)'INTEGRATION'

! Integration     
!   Reset counter
      npoints=0  
!   Reset various iterative quantities
      if(ifinal.eq.0)then
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
          do jasy=1,5 ! nasy-3
            do iasy=-1,+1,2
              asycross(jasy,i,iasy)=0.d0
              asyerror(jasy,i,iasy)=0.d0
            end do
          end do 
        end do
      end if
!   Integrate
      it=0
      call vegas(ndim,fxn,avgi,sd,chi2a)
      if(ifinal.eq.0)then
!   Multiply by branching ratios (if ifinal = 0)      
        avgi=avgi*(BRtbln)**2
        sd=sd*(BRtbln)**2
      end if
!  Collect total cross-section
      cross=avgi
      error=sd

! Integration output
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

      if(ifinal.eq.0)then  
! Collect polarised cross sections.
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
!     & sqrt(abs(polchi(iphel,jphel)
!     &         -poltot(iphel,jphel)**2*dfloat(ncall)))
!     & /dfloat(ncall)
          end do
        end do
      end if

! Collect unpolarised spatial asymmetry
      do iasy=1,5 ! nasy-3
!         write(*,*)iasy,m_asy(iasy+3)
        if(m_asy(iasy+3).eq.0)then
          continue
        else
          do iAB=-1,+1,2
            do i=1,it
              asycross(iasy,i,iAB)=asycross(iasy,i,iAB)
     &                            *avgi/cnorm(i)
              asyerror(iasy,i,iAB)=asycross(iasy,i,iAB)
     &                            *sd/cnorm(i)
              write(*,*)asycross(iasy,i,iAB)
            end do        
            asytot(iasy,iAB)=0.d0
            asychi(iasy,iAB)=0.d0
            do i=1,it ! add up each iteration
              asytot(iasy,iAB)=asytot(iasy,iAB)
     &                      +asycross(iasy,i,iAB)
              asychi(iasy,iAB)=asychi(iasy,iAB)
     &                      +asyerror(iasy,i,iAB)
            end do
!             write(*,*)'bork1011'
            asychi(iasy,iAB)=asychi(iasy,iAB)
     &                        /asytot(iasy,iAB)
!             write(*,*)'bork1014'
!           asychi(iasy)=
!        & sqrt(abs(asychi(iasy)
!        &         -asytot(iasy)**2*dfloat(ncall)))
!        & /dfloat(ncall)
          end do
        end if
      end do

! Define asymmetries
      if(ifinal.eq.0)then
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
! AFBcm
      Atot(4)=
     &          +(asytot(1,+1)-asytot(1,-1))
     &          /cross
        Atoterr(4)=
     &            +sd/avgi*Atot(4)
! AtFB
      Atot(5)=
     &          +(asytot(2,+1)-asytot(2,-1))
     &          /cross
      Atoterr(5)=
     &            +sd/avgi*Atot(5)

      if(m_asy(6).gt.0)then
! A
        Atot(6)=
     &          +(asytot(3,+1)-asytot(3,-1))
     &          /cross
        Atoterr(6)=
     &            +sd/avgi*Atot(6)
      end if

      if(m_asy(7).gt.0)then
! A'
        Atot(7)=
     &          +(asytot(4,+1)-asytot(4,-1))
     &          /cross
        Atoterr(7)=
     &            +sd/avgi*Atot(7)
      end if

      
      if(m_asy(8).gt.0)then
! A_l
        Atot(8)=
     &          +(asytot(5,+1)-asytot(5,-1))
     &          /cross
        Atoterr(8)=
     &            +sd/avgi*Atot(8)
      end if

! Print Asymmetries     
      write(*,*)'TOTAL ASYMMETRIES'
      if(ifinal.eq.0)then
        write(*,*)'ALL:                  uncertainty (same units):'
        write(*,*)Atot(1),Atoterr(1) 
        write(*,*)'AL:                   uncertainty (same units):'
        write(*,*)Atot(2),Atoterr(2) 
        write(*,*)'APV:                  uncertainty (same units):'
        write(*,*)Atot(3),Atoterr(3) 
        write(*,*)'AFB*:                 uncertainty (same units):'
        write(*,*)Atot(4),Atoterr(4)
        write(*,*)'AtFB:                    uncertainty (same units):'
        write(*,*)Atot(5),Atoterr(5)
        write(*,*)'A:                    uncertainty (same units):'
        write(*,*)Atot(6),Atoterr(6)
        write(*,*)"A':                    uncertainty (same units):"
        write(*,*)Atot(7),Atoterr(7)
      else if(ifinal.gt.0)then
        write(*,*)'A_l:                  uncertainty (same units):'
        write(*,*)Atot(6),Atoterr(6) 
      end if

! Plot Distributions
      write(*,*)'-----------------------------------------------------'
      write(*,*)'HISTOGRAMS'

!  plot distributions in transverse momentum
      do ip=1,8
        if(m_pT(ip).eq.0)then
          continue
        else
          sfxpTtot(ip)=0d0
          do j=1,ndiv_sigp
            fxpTtot(ip,j)=0.d0
            do i=1,it
              fxpT(ip,j,i)=fxpT(ip,j,i)*avgi/cnorm(i)/pTw(ip)
              fxpTtot(ip,j)=fxpTtot(ip,j)+fxpT(ip,j,i)
            end do
            sfxpTtot(ip)=sfxpTtot(ip)+fxpTtot(ip,j)*pTw(ip)
          end do          
          write(*,*)'HISTOGRAM'
          write(*,'(A,I1)')'pT',ip
          write(*,'(A,I1,A)')'d#sigma-/dp_{T}(',ip,')--[pb/GeV]'
          write(*,'(A,I1,A)')'p_{T}(',ip,')--[GeV]'
          do i=1,ndiv_pT(ip)
            write(*,*)xpT(ip,i),fxpTtot(ip,i)
          end do
          write(*,*)'END'
          write(*,*)'(Integrated cross-section:',sfxpTtot(ip),')'
        end if
      end do

      if(m_pT356.eq.1)then
!   plot distribution in pT356
        sfxpT356tot=0d0
        do j=1,ndiv_pT356
          fxpT356tot(j)=0.d0
          do i=1,it
            fxpT356(j,i)=fxpT356(j,i)*avgi/cnorm(i)/pT356w
            fxpT356tot(j)=fxpT356tot(j)+fxpT356(j,i)            
          end do
          sfxpT356tot=sfxpT356tot+fxpT356tot(j)*pT356w
        end do
        write(*,*)'HISTOGRAM'
        write(*,*)'pT356'
        write(*,*)'d#sigma-/dp_{T}--[pb/GeV]'
        write(*,*)'p_{T}(t)--[GeV]'
        do i=1,ndiv_pT356
          write(*,*)xpT356(i),fxpT356tot(i)
        end do
        write(*,*)'END'
      end if

      if(m_pT478.eq.1)then
!   plot distribution in pT478
        sfxpT478tot=0d0
        do j=1,ndiv_pT478
          fxpT478tot(j)=0.d0
          do i=1,it
            fxpT478(j,i)=fxpT478(j,i)*avgi/cnorm(i)/pT478w
            fxpT478tot(j)=fxpT478tot(j)+fxpT478(j,i)            
          end do
          sfxpT478tot=sfxpT478tot+fxpT478tot(j)*pT478w
        end do
        write(*,*)'HISTOGRAM'
        write(*,*)'pT478'
        write(*,*)'d#sigma-/dp_{T}--[pb/GeV]'
        write(*,*)'p_{T}(#bar{t})--[GeV]'
        do i=1,ndiv_pT478
          write(*,*)xpT478(i),fxpT478tot(i)
        end do
        write(*,*)'END'
      end if

      if(m_ETmiss.eq.1)then
!   plot distribution in ETmiss
        sfxETmisstot=0d0
        do j=1,ndiv_ETmiss
          fxETmisstot(j)=0.d0
          do i=1,it
            fxETmiss(j,i)=fxETmiss(j,i)*avgi/cnorm(i)/ETmissw
            fxETmisstot(j)=fxETmisstot(j)+fxETmiss(j,i)            
          end do
          sfxETmisstot=sfxETmisstot+fxETmisstot(j)*ETmissw
        end do
        write(*,*)'HISTOGRAM'
        write(*,*)'ETmiss'
        write(*,*)'d#sigma-/dp_{Tmiss}--[pb/GeV]'
        write(*,*)'p_{T}(miss)--[GeV]'
        do i=1,ndiv_ETmiss
          write(*,*)xETmiss(i),fxETmisstot(i)
        end do
        write(*,*)'END'
      end if

      if(m_eta3.eq.1)then
!   plot distdistribution in eta3
        sfxeta3tot=0d0
        do j=1,ndiv_eta3
          fxeta3tot(j)=0.d0
          do i=1,it
            fxeta3(j,i)=fxeta3(j,i)*avgi/cnorm(i)/eta3w
            fxeta3tot(j)=fxeta3tot(j)+fxeta3(j,i)            
          end do
          sfxeta3tot=sfxeta3tot+fxeta3tot(j)*eta3w
        end do
        write(*,*)'HISTOGRAM'
        write(*,*)'eta3'
        write(*,*)'d#sigma-/d#eta--[pb]'
        if(ifinal.eq.0) write(*,*)'#eta(t)'
        if(ifinal.eq.1) write(*,*)'#eta(b)'
        do i=1,ndiv_eta3
          write(*,*)xeta3(i),fxeta3tot(i)
        end do
        write(*,*)'END'
      end if

      if(m_eta4.eq.1)then
!   plot distribution in eta4
        sfxeta4tot=0d0
        do j=1,ndiv_eta4
          fxeta4tot(j)=0.d0
          do i=1,it
            fxeta4(j,i)=fxeta4(j,i)*avgi/cnorm(i)/eta4w
            fxeta4tot(j)=fxeta4tot(j)+fxeta4(j,i)            
          end do
          sfxeta4tot=sfxeta4tot+fxeta4tot(j)*eta4w
        end do
        write(*,*)'HISTOGRAM'
        write(*,*)'eta4'
        write(*,*)'d#sigma-/d#eta--[pb]'
        if(ifinal.eq.0)write(*,*)'#eta(#bar{t})'
        if(ifinal.eq.1)write(*,*)'#eta(#bar{b})'
        do i=1,ndiv_eta4
          write(*,*)xeta4(i),fxeta4tot(i)
        end do
        write(*,*)'END'
      end if

      if(m_eta5.eq.1)then
!   plot distribution in eta5
        sfxeta5tot=0d0
        do j=1,ndiv_eta5
          fxeta5tot(j)=0.d0
          do i=1,it
            fxeta5(j,i)=fxeta5(j,i)*avgi/cnorm(i)/eta5w
            fxeta5tot(j)=fxeta5tot(j)+fxeta5(j,i)            
          end do
          sfxeta5tot=sfxeta5tot+fxeta5tot(j)*eta5w
        end do
        write(*,*)'HISTOGRAM'
        write(*,*)'eta5'
        write(*,*)'d#sigma-/d#eta--[pb]'
        write(*,*)'#eta(l^{+})'
        do i=1,ndiv_eta5
          write(*,*)xeta5(i),fxeta5tot(i)
        end do
        write(*,*)'END'
      end if

      if(m_eta6.eq.1)then
!   plot distribution in eta6
        sfxeta6tot=0d0
        do j=1,ndiv_eta6
          fxeta6tot(j)=0.d0
          do i=1,it
            fxeta6(j,i)=fxeta6(j,i)*avgi/cnorm(i)/eta6w
            fxeta6tot(j)=fxeta6tot(j)+fxeta6(j,i)            
          end do
          sfxeta6tot=sfxeta6tot+fxeta6tot(j)*eta6w
        end do
        write(*,*)'HISTOGRAM'
        write(*,*)'eta6'
        write(*,*)'d#sigma-/d#eta--[pb]'
        write(*,*)'#eta(#bar{#nu})'
        do i=1,ndiv_eta6
          write(*,*)xeta6(i),fxeta6tot(i)
        end do
        write(*,*)'END'
      end if

      if(m_eta7.eq.1)then
!   plot distribution in eta7
        sfxeta7tot=0d0
        do j=1,ndiv_eta7
          fxeta7tot(j)=0.d0
          do i=1,it
            fxeta7(j,i)=fxeta7(j,i)*avgi/cnorm(i)/eta7w
            fxeta7tot(j)=fxeta7tot(j)+fxeta7(j,i)            
          end do
          sfxeta7tot=sfxeta7tot+fxeta7tot(j)*eta7w
        end do
        write(*,*)'HISTOGRAM'
        write(*,*)'eta7'
        write(*,*)'d#sigma-/d#eta--[pb]'
        write(*,*)'#eta(l^{-})'
        do i=1,ndiv_eta7
          write(*,*)xeta7(i),fxeta7tot(i)
        end do
        write(*,*)'END'
      end if

      if(m_eta8.eq.1)then
!   plot distribution in eta8
        sfxeta8tot=0d0
        do j=1,ndiv_eta8
          fxeta8tot(j)=0.d0
          do i=1,it
            fxeta8(j,i)=fxeta8(j,i)*avgi/cnorm(i)/eta8w
            fxeta8tot(j)=fxeta8tot(j)+fxeta8(j,i)            
          end do
          sfxeta8tot=sfxeta8tot+fxeta8tot(j)*eta8w
        end do
        write(*,*)'HISTOGRAM'
        write(*,*)'eta8'
        write(*,*)'d#sigma-/d#eta--[pb]'
        write(*,*)'#eta(#bar{#nu})'
        do i=1,ndiv_eta8
          write(*,*)xeta8(i),fxeta8tot(i)
        end do
        write(*,*)'END'
      end if

      if(m_eta356.eq.1)then
!   plot distribution in eta356
        do j=1,ndiv_eta356
          do i=1,it
            fxeta356(j,i)=fxeta356(j,i)*avgi/cnorm(i)/eta356w
          end do
        end do
        do j=1,ndiv_eta356
          fxeta356tot(j)=0.d0 
          do i=1,it
            fxeta356tot(j)=fxeta356tot(j)+fxeta356(j,i)
          end do
        end do
        write(*,*)'HISTOGRAM'
        write(*,*)'eta356'
        write(*,*)'d#sigma-/d#eta--[pb]'
        write(*,*)'#eta(t)'
        do i=1,ndiv_eta356
          write(*,*)xeta356(i),fxeta356tot(i)
        end do
        write(*,*)'END'
      end if

      if(m_eta478.eq.1)then
!   plot distribution in eta478
        sfxeta478tot=0d0
        do j=1,ndiv_eta478
          fxeta478tot(j)=0.d0
          do i=1,it
            fxeta478(j,i)=fxeta478(j,i)*avgi/cnorm(i)/eta478w
            fxeta478tot(j)=fxeta478tot(j)+fxeta478(j,i)            
          end do
          sfxeta478tot=sfxeta478tot+fxeta478tot(j)*eta478w
        end do
        write(*,*)'HISTOGRAM'
        write(*,*)'eta478'
        write(*,*)'d#sigma-/d#eta--[pb]'
        write(*,*)'#eta(#bar{t})'
        do i=1,ndiv_eta478
          write(*,*)xeta478(i),fxeta478tot(i)
        end do
        write(*,*)'END'
      end if      

      if(m_rmass.eq.1)then
!   plot distribution in rmass.
        sfxrmasstot=0d0
        do j=1,ndiv_rmass
          fxrmasstot(j)=0.d0
          do i=1,it
            fxrmass(j,i)=fxrmass(j,i)*avgi/cnorm(i)/rmassw
            fxrmasstot(j)=fxrmasstot(j)+fxrmass(j,i)            
          end do
          sfxrmasstot=sfxrmasstot+fxrmasstot(j)*rmassw
        end do
        write(*,*)'HISTOGRAM'
        write(*,*)'Mtt'
        write(*,*)'d#sigma-/dM_{tt}--[pb/GeV]'
        write(*,*)'M_{tt}--[GeV]'
        do i=1,ndiv_rmass
          write(*,*)xrmass(i),fxrmasstot(i)
        end do
        write(*,*)'END'
      end if

      if(m_rMvis.eq.1)then
!   plot distribution in rMvis.
        sfxrMvistot=0d0
        do j=1,ndiv_rMvis
          fxrMvistot(j)=0.d0
          do i=1,it
            fxrMvis(j,i)=fxrMvis(j,i)*avgi/cnorm(i)/rMvisw
            fxrMvistot(j)=fxrMvistot(j)+fxrMvis(j,i)            
          end do
          sfxrMvistot=sfxrMvistot+fxrMvistot(j)*rMvisw
        end do
        write(*,*)'HISTOGRAM'
        write(*,*)'Mvis'
        write(*,*)'d#sigma-/dM_{vis}--[pb/GeV]'
        write(*,*)'M_{vis}--[GeV]'
        do i=1,ndiv_rMvis
          write(*,*)xrMvis(i),fxrMvistot(i)
        end do
        write(*,*)'END'
      end if

      if(m_beta.eq.1)then
!   plot distribution in beta.
        sfxbetatot=0d0
        do j=1,ndiv_beta
          fxbetatot(j)=0.d0
          do i=1,it
            fxbeta(j,i)=fxbeta(j,i)*avgi/cnorm(i)/betaw
            fxbetatot(j)=fxbetatot(j)+fxbeta(j,i)            
          end do
          sfxbetatot=sfxbetatot+fxbetatot(j)*betaw
        end do
        write(*,*)'HISTOGRAM'
        write(*,*)'Beta'
        write(*,*)'d#sigma-/d#Beta_{t}--[pb]'
        write(*,*)'#Beta_{t}'
        do i=1,ndiv_beta
          write(*,*)xbeta(i),fxbetatot(i)
        end do
        write(*,*)'END'
      end if  

      if(m_cost.eq.1)then
!   plot distribution in cost.
        sfxcosttot=0d0
        do j=1,ndiv_cost
          fxcosttot(j)=0.d0
          do i=1,it
            fxcost(j,i)=fxcost(j,i)*avgi/cnorm(i)/costw
            fxcosttot(j)=fxcosttot(j)+fxcost(j,i)            
          end do
          sfxcosttot=sfxcosttot+fxcosttot(j)*costw
        end do
        write(*,*)'HISTOGRAM'
        write(*,*)'cost'
        write(*,*)'d#sigma-/dcos#theta--[pb]'
        write(*,*)'cos#theta'
        do i=1,ndiv_cost
          write(*,*)xcost(i),fxcosttot(i)
        end do
        write(*,*)'END'
      end if

      if(m_fl.eq.1)then
!   plot distribution in fl.
        sfxfltot=0d0
        do j=1,ndiv_fl
          fxfltot(j)=0.d0
          do i=1,it
            fxfl(j,i)=fxfl(j,i)*avgi/cnorm(i)/flw
            fxfltot(j)=fxfltot(j)+fxfl(j,i)            
          end do
          sfxfltot=sfxfltot+fxfltot(j)*flw
        end do
        write(*,*)'HISTOGRAM'
        write(*,*)'fl'
        write(*,*)'d#sigma-/d#phi_{l}--[pb]'
        write(*,*)'#phi_{l}--[-rad-]'
        do i=1,ndiv_fl
          write(*,*)xfl(i),fxfltot(i)
        end do
        write(*,*)'END'
      end if 

      if(m_cosfl.eq.1)then
!   plot distribution in cosfl.
        sfxcosfltot=0d0
        do j=1,ndiv_cosfl
          fxcosfltot(j)=0.d0
          do i=1,it
            fxcosfl(j,i)=fxcosfl(j,i)*avgi/cnorm(i)/cosflw
            fxcosfltot(j)=fxcosfltot(j)+fxcosfl(j,i)            
          end do
          sfxcosfltot=sfxcosfltot+fxcosfltot(j)*cosflw
        end do
        write(*,*)'HISTOGRAM'
        write(*,*)'cosfl'
        write(*,*)'d#sigma-/dcos#phi_{l}--[pb]'
        write(*,*)'cos#phi_{l}'
        do i=1,ndiv_cosfl
          write(*,*)xcosfl(i),fxcosfltot(i)
        end do
        write(*,*)'END'
      end if

      if(m_Et.eq.1)then
!   plot distribution in Et.
        sfxEttot=0d0
        do j=1,ndiv_Et
          fxEttot(j)=0.d0
          do i=1,it
            fxEt(j,i)=fxEt(j,i)*avgi/cnorm(i)/Etw
            fxEttot(j)=fxEttot(j)+fxEt(j,i)            
          end do
          sfxEttot=sfxEttot+fxEttot(j)*Etw
        end do
        write(*,*)'HISTOGRAM'
        write(*,*)'Et'
        write(*,*)'d#sigma-/dE_{t}--[pb/GeV]'
        write(*,*)'E_{t}--[GeV]'
        do i=1,ndiv_Et
          write(*,*)xEt(i),fxEttot(i)
        end do
        write(*,*)'END'
      end if

      if(m_HT.eq.1)then
!   plot distribution in scalar sum of ET: HT.
        sfxHTtot=0d0
        do j=1,ndiv_HT
          fxHTtot(j)=0.d0
          do i=1,it
            fxHT(j,i)=fxHT(j,i)*avgi/cnorm(i)/HTw
            fxHTtot(j)=fxHTtot(j)+fxHT(j,i)            
          end do
          sfxHTtot=sfxHTtot+fxHTtot(j)*HTw
        end do
        write(*,*)'HISTOGRAM'
        write(*,*)'HT'
        write(*,*)'d#sigma-/dH_{T}--[pb/GeV]'
        write(*,*)'H_{T}--[GeV]'
        do i=1,ndiv_HT
          write(*,*)xHT(i),fxHTtot(i)
        end do
        write(*,*)'END'
      end if

      if(m_rM_T.eq.1)then
!   plot distribution in M_T.
        sfxrM_Ttot=0d0
        do j=1,ndiv_rM_T
          fxrM_Ttot(j)=0.d0
          do i=1,it
            fxrM_T(j,i)=fxrM_T(j,i)*avgi/cnorm(i)/rM_Tw
            fxrM_Ttot(j)=fxrM_Ttot(j)+fxrM_T(j,i)            
          end do
          sfxrM_Ttot=sfxrM_Ttot+fxrM_Ttot(j)*rM_Tw
        end do
        write(*,*)'HISTOGRAM'
        write(*,*)'M_T'
        write(*,*)'d#sigma-/dM_{T}--[pb/GeV]'
        write(*,*)'M_{T}--[GeV]'
        do i=1,ndiv_rM_T
          write(*,*)xrM_T(i),fxrM_Ttot(i)
        end do
        write(*,*)'END'
      end if

      if(m_rM_CT.eq.1)then
!   plot distribution in Mcon.
        sfxrM_CTtot=0d0
        do j=1,ndiv_rM_CT
          fxrM_CTtot(j)=0.d0
          do i=1,it
            fxrM_CT(j,i)=fxrM_CT(j,i)*avgi/cnorm(i)/rM_CTw
            fxrM_CTtot(j)=fxrM_CTtot(j)+fxrM_CT(j,i)            
          end do
          sfxrM_CTtot=sfxrM_CTtot+fxrM_CTtot(j)*rM_CTw
        end do
        write(*,*)'HISTOGRAM'
        write(*,*)'M_CT'
        write(*,*)'d#sigma-/dM_{CT}--[pb/GeV]'
        write(*,*)'M_{CT}--[GeV]'
        do i=1,ndiv_rM_CT
          write(*,*)xrM_CT(i),fxrM_CTtot(i)
        end do
        write(*,*)'END'
      end if

      if(m_rMlCT.eq.1)then
!   plot distribution in Mcon.
        sfxrMlCTtot=0d0
        do j=1,ndiv_rMlCT
          fxrMlCTtot(j)=0.d0
          do i=1,it
            fxrMlCT(j,i)=fxrMlCT(j,i)*avgi/cnorm(i)/rMlCTw
            fxrMlCTtot(j)=fxrMlCTtot(j)+fxrMlCT(j,i)            
          end do
          sfxrMlCTtot=sfxrMlCTtot+fxrMlCTtot(j)*rMlCTw
        end do
        write(*,*)'HISTOGRAM'
        write(*,*)'MlCT'
        write(*,*)'d#sigma-/dM^{l}_{CT}--[pb/GeV]'
        write(*,*)'M^{l}_{CT}--[GeV]'
        do i=1,ndiv_rMlCT
          write(*,*)xrMlCT(i),fxrMlCTtot(i)
        end do
        write(*,*)'END'
      end if


!   plot distributions in all asymmetries.
      if((m_sigp.eq.1).and.(m_sigm.eq.1))then
        do jasy=1,8
          if(m_asy(jasy).eq.0)then
            continue
          else
  !           snorm(jasy)=0.d0
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
            write(*,*)'HISTOGRAM'
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
                write(*,*)'AFBcm'
                write(*,*)'A_{FB^{*}}'
            else if(jasy.eq.5)then
                write(*,*)'AtFB'
                write(*,*)'A^{t}_{FB}'
            else if(jasy.eq.6)then
                write(*,*)'A'
                write(*,*)'A'
            else if(jasy.eq.7)then
                write(*,*)'Ap'
                write(*,*)"A'"
            else if(jasy.eq.8)then
                write(*,*)'A_l'
                write(*,*)'A_{l}'
            end if            
            write(*,*)'M_{tt}'
            ndiv_sig=(ndiv_sigm+ndiv_sigp)/2
            do i=1,ndiv_sig
              if(fxsigptot(jasy,i)+fxsigmtot(jasy,i).eq.0.d0)then
                write(*,*)(xsigm(i)+xsigp(i))/2.d0,0.d0
  !               snorm(jasy)=snorm(jasy)+0.d0
              else  
                write(*,*)(xsigm(i)+xsigp(i))/2.d0,
     &                   (fxsigptot(jasy,i)-fxsigmtot(jasy,i))/
     &                   (fxsigptot(jasy,i)+fxsigmtot(jasy,i))
  !               snorm(jasy)=snorm(jasy)+
  !      &               (fxsigptot(jasy,i)-fxsigmtot(jasy,i))/
  !      &               (fxsigptot(jasy,i)+fxsigmtot(jasy,i))
  !      &               *fxrmasstot(i)*rmassw/avgi      
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

! Check distributions
      diff_max=1E-12
      n_error=0
!       if(m_pT3.eq.1)then
!         if(abs(cross-sfxpT3tot)>diff_max)then
!           write(*,*)'pT3 Integration Error:',sfxpT3tot
!           n_error=n_error+1
!         end if
!       end if
!       if(m_pT4.eq.1)then
!         if(abs(cross-sfxpT4tot)>diff_max)then
!           write(*,*)'pT4 Integration Error:',sfxpT4tot
!           n_error=n_error+1
!         end if
!       end if
!       if(m_pT5.eq.1)then
!         if(abs(cross-sfxpT5tot)>diff_max)then
!           write(*,*)'pT5 Integration Error:',sfxpT5tot
!           n_error=n_error+1
!         end if
!       end if
!       if(m_pT6.eq.1)then
!         if(abs(cross-sfxpT6tot)>diff_max)then
!           write(*,*)'pT6 Integration Error:',sfxpT6tot
!           n_error=n_error+1
!         end if
!       end if
!       if(m_pT7.eq.1)then
!         if(abs(cross-sfxpT7tot)>diff_max)then
!           write(*,*)'pT7 Integration Error:',sfxpT7tot
!           n_error=n_error+1
!         end if
!       end if
!       if(m_pT8.eq.1)then
!         if(abs(cross-sfxpT8tot)>diff_max)then
!           write(*,*)'pT8 Integration Error:',sfxpT8tot
!           n_error=n_error+1
!         end if
!       end if
      if(m_pT356.eq.1)then
        if(abs(cross-sfxpT356tot)>diff_max)then
          write(*,*)'pT356 Integration Error:',sfxpT356tot
          n_error=n_error+1
        end if
      end if
      if(m_pT478.eq.1)then
        if(abs(cross-sfxpT478tot)>diff_max)then
          write(*,*)'pT478 Integration Error:',sfxpT478tot
          n_error=n_error+1
        end if
      end if
      if(m_ETmiss.eq.1)then
        if(abs(cross-sfxETmisstot)>diff_max)then
          write(*,*)'ETmiss Integration Error:',sfxETmisstot
          n_error=n_error+1
        end if
      end if
      if(m_eta3.eq.1)then
        if(abs(cross-sfxeta3tot)>diff_max)then
          write(*,*)'eta3 Integration Error:',sfxeta3tot
          n_error=n_error+1
        end if
      end if
      if(m_eta4.eq.1)then
        if(abs(cross-sfxeta4tot)>diff_max)then
          write(*,*)'eta4 Integration Error:',sfxeta4tot
          n_error=n_error+1
        end if
      end if
      if(m_eta5.eq.1)then
        if(abs(cross-sfxeta5tot)>diff_max)then
          write(*,*)'eta5 Integration Error:',sfxeta5tot
          n_error=n_error+1
        end if
      end if
      if(m_eta6.eq.1)then
        if(abs(cross-sfxeta6tot)>diff_max)then
          write(*,*)'eta6 Integration Error:',sfxeta6tot
          n_error=n_error+1
        end if
      end if
      if(m_eta7.eq.1)then
        if(abs(cross-sfxeta7tot)>diff_max)then
          write(*,*)'eta7 Integration Error:',sfxeta7tot
          n_error=n_error+1
        end if
      end if
      if(m_eta8.eq.1)then
        if(abs(cross-sfxeta8tot)>diff_max)then
          write(*,*)'eta8 Integration Error:',sfxeta8tot
          n_error=n_error+1
        end if
      end if
      if(m_rmass.eq.1)then
        if(abs(cross-sfxrmasstot)>diff_max)then
          write(*,*)'rmass Integration Error:',sfxrmasstot
          n_error=n_error+1
        end if
      end if
      if(m_rMvis.eq.1)then
        if(abs(cross-sfxrMvistot)>diff_max)then
          write(*,*)'rMvis Integration Error:',sfxrMvistot
          n_error=n_error+1
        end if
      end if
      if(m_beta.eq.1)then
        if(abs(cross-sfxbetatot)>diff_max*10)then
          write(*,*)'beta Integration Error:',sfxbetatot
          n_error=n_error+1
        end if
      end if
      if(m_cost.eq.1)then
        if(abs(cross-sfxcosttot)>diff_max)then
          write(*,*)'cost Integration Error:',sfxcosttot
          n_error=n_error+1
        end if
      end if
      if(m_Et.eq.1)then
        if(abs(cross-sfxEttot)>diff_max)then
          write(*,*)'Et Integration Error:',sfxEttot
          n_error=n_error+1
        end if
      end if
      if(m_rM_T.eq.1)then
        if(abs(cross-sfxrM_Ttot)>diff_max)then
          write(*,*)'M_T Integration Error:',sfxrM_Ttot
          n_error=n_error+1
        end if
      end if
      if(m_rM_CT.eq.1)then
        if(abs(cross-sfxrM_CTtot)>diff_max)then
          write(*,*)'M_CT1 Integration Error:',sfxrM_CTtot
          n_error=n_error+1
        end if
      end if
      if(m_rMlCT.eq.1)then
        if(abs(cross-sfxrMlCTtot)>diff_max)then
          write(*,*)'M_CT2 Integration Error:',sfxrMlCTtot
          n_error=n_error+1
        end if
      end if
      if(m_fl.eq.1)then
        if(abs(cross-sfxfltot)>diff_max)then
          write(*,*)'fl Integration Error:',sfxfltot
          n_error=n_error+1
        end if
      end if
      if(m_cosfl.eq.1)then
        if(abs(cross-sfxcosfltot)>diff_max)then
          write(*,*)'cosfl Integration Error:',sfxcosfltot
          n_error=n_error+1
        end if
      end if
      do iasy=1,8 !nasy
        if(m_asy(iasy).eq.0)then
          continue
        else
          if(abs(Atot(iasy)-asym_int(iasy))>diff_max)then
            write(*,*)'A Integration Error:',iasy,asym_int(iasy)
            n_error=n_error+1
          end if
        end if
      end do
      write(*,*)'INTEGRATION ERRORS:',n_error
      write(*,*)'CLOSE'
      write(*,*)'======================================================'
      stop
      end
! ======================================================================