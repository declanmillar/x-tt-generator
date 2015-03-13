module Distributions

  ! use Configuration, only: o_distros, o_dist2d, o_trans, o_asyms, ifinal_state

  ! switch for all distributions
  ! integer :: o_distros

  ! distributions in pTs of asymptotic particles
  real :: pTmax(8),pTmin(8),pTw(8)
  real :: xpT(8,500),fxpT(8,500,20),fxpTtot(8,500)
  integer :: o_pT(8)
  integer :: ndiv_pT(8)

  ! distributions in etas of asymptotic particles
  real :: etamax(8),etamin(8),etaw(8)
  real :: xeta(8,500),fxeta(8,500,20),fxetatot(8,500)
  integer :: o_eta(8)
  integer :: ndiv_eta(8)

  ! distributions in phis of asymptotic particles
  real :: phimax(8),phimin(8),phiw(8)
  real :: xphi(8,500),fxphi(8,500,20),fxphitot(8,500)
  integer :: o_phi(8)
  integer :: ndiv_phi(8)

  ! distributions in ycol of asymptotic particles
  real :: ycolmax(8),ycolmin(8),ycolw(8)
  real :: xycol(8,500),fxycol(8,500,20),fxycoltot(8,500)
  integer :: o_ycol(8)
  integer :: ndiv_ycol(8)

  ! distribution in ETmiss
  real :: ETmissmax,ETmissmin,ETmissw
  real :: xETmiss(500),fxETmiss(500,20),fxETmisstot(500)
  integer :: o_ETmiss
  integer :: ndiv_ETmiss

  ! distribution in pT of the top
  real :: pT356max,pT356min,pT356w
  real :: xpT356(500),fxpT356(500,20),fxpT356tot(500)
  integer :: o_pT356
  integer :: ndiv_pT356

  ! distribution in eta of the top
  real :: eta356max,eta356min,eta356w
  real :: xeta356(500),fxeta356(500,20),fxeta356tot(500)
  integer :: o_eta356
  integer :: ndiv_eta356

  ! distribution in phi of the top
  real :: phi356max,phi356min,phi356w
  real :: xphi356(500),fxphi356(500,20),fxphi356tot(500)
  integer :: o_phi356
  integer :: ndiv_phi356

  ! distribution in pT of the anti-top
  real :: pT478max,pT478min,pT478w
  real :: xpT478(500),fxpT478(500,20),fxpT478tot(500)
  integer :: o_pT478
  integer :: ndiv_pT478

  ! distribution in eta of the anti-top
  real :: eta478max,eta478min,eta478w
  real :: xeta478(500),fxeta478(500,20),fxeta478tot(500)
  integer :: o_eta478
  integer :: ndiv_eta478

  ! distribution in phi of the anti-top
  real :: phi478max,phi478min,phi478w
  real :: xphi478(500),fxphi478(500,20),fxphi478tot(500)
  integer :: o_phi478
  integer :: ndiv_phi478

  ! distribution in invarient mass of the top pair
  real :: rMttmax,rMttmin,rMttw
  real :: xrMtt(500),fxrMtt(500,20),fxrMtttot(500)
  integer :: o_rMtt
  integer :: ndiv_rMtt

  ! distribution in boost of top pair centre of mass frame
  real :: betamax,betamin,betaw
  real :: xbeta(500),fxbeta(500,20),fxbetatot(500)
  integer :: o_beta
  integer :: ndiv_beta

  ! distribution in cos(theta_t)
  real :: costmax,costmin,costw
  real :: xcost(500),fxcost(500,20),fxcosttot(500)
  integer :: o_cost
  integer :: ndiv_cost

  ! distribution in top energy
  real :: Etmax,Etmin,Etw
  real :: xEt(500),fxEt(500,20),fxEttot(500)
  integer :: o_Et
  integer :: ndiv_Et

  ! distribution in Delta_y
  real :: Delta_ymax,Delta_ymin,Delta_yw
  real :: xDelta_y(500),fxDelta_y(500,20), fxDelta_ytot(500)
  integer :: o_Delta_y
  integer :: ndiv_Delta_y

  ! distributions in transverse variables
  integer :: ntrans
  parameter (ntrans=10)
  real :: transmax(ntrans),transmin(ntrans),transw(ntrans)
  real :: xtrans(ntrans,500),fxtrans(ntrans,500,20) ,fxtranstot(ntrans,500)
  integer :: o_tran(ntrans)
  integer :: ndiv_trans(ntrans)
  real :: sfxtranstot(ntrans),sfxtransdptot(ntrans)

  ! distribution in phi_l (lepton azimuthal angle)
  real :: flmax,flmin,flw
  real :: xfl(500),fxfl(500,20),fxfltot(500)
  integer :: o_fl
  integer :: ndiv_fl

  ! distribution in cos_phi_l
  real :: cosflmax,cosflmin,cosflw
  real :: xcosfl(500),fxcosfl(500,20),fxcosfltot(500)
  integer :: o_cosfl
  integer :: ndiv_cosfl

  ! distribution in delta phi
  real :: dphimax,dphimin,dphiw
  real :: xdphi(500),fxdphi(500,20),fxdphitot(500)
  integer :: o_dphi
  integer :: ndiv_dphi

  ! distribution in cost5
  real :: cost5max,cost5min,cost5w
  real :: xcost5(500),fxcost5(500,20),fxcost5tot(500)
  integer :: o_cost5
  integer :: ndiv_cost5

  ! distribution in cost7
  real :: cost7max,cost7min,cost7w
  real :: xcost7(500),fxcost7(500,20),fxcost7tot(500)
  integer :: o_cost7
  integer :: ndiv_cost7

  ! distribution in cost7
  real :: ct7ct5max,ct7ct5min,ct7ct5w
  real :: xct7ct5(500),fxct7ct5(500,20),fxct7ct5tot(500)
  integer :: o_ct7ct5
  integer :: ndiv_ct7ct5

  ! distributions in asymmetries
  integer :: nasym
  parameter (nasym=9)
  integer :: nspat
  parameter (nspat=6) ! nasym-3
  integer :: o_asym(nasym)

  ! distribution in sigp
  real :: sigpmax,sigpmin,sigpw
  real :: xsigp(1000),fxsigp(nasym,1000,20) ,fxsigptot(nasym,1000)
  integer :: o_sigp
  integer :: ndiv_sigp

  ! distribution in sigm
  real :: sigmmax,sigmmin,sigmw
  real :: xsigm(1000),fxsigm(nasym,1000,20) ,fxsigmtot(nasym,1000)
  integer :: o_sigm
  integer :: ndiv_sigm

  ! 2d dist
  real :: xdphi2d(500,500),fxdphi2d(500,500,20), fxdphi2dtot(500,500)
  integer :: o_dphi2d

  ! trans 2d dist
  real :: xtransdp(ntrans,500,500), fxtransdp(ntrans,500,500,20),fxtransdptot(ntrans,500,500)
  integer :: o_transdp(ntrans)

  ! public :: initialise_distributions
  ! public :: 

  ! contains

  !   subroutine initialise_distributions

  !     ! Set flags, binning range and divisions

  !     ! pT distributions
  !     do ip=1,nfinal
  !       o_pT(ip)=o_distros
  !       pTmax(ip)=7000.d0/(1+initial_state*6)
  !       pTmin(ip)=0.d0
  !       ndiv_pT(ip)=70
  !       ! eta distributions
  !       o_eta(ip)=o_distros
  !       etamax(ip)=+10
  !       etamin(ip)=-10
  !       ndiv_eta(ip)=50
  !       ! phi distributions
  !       o_phi(ip)=o_distros
  !       phimax(ip)=+pi
  !       phimin(ip)=-pi
  !       ndiv_phi(ip)=50
  !       ! ycol distributions
  !       o_ycol(ip)=o_distros
  !       ycolmax(ip)=+4.d0
  !       ycolmin(ip)=-4.d0
  !       ndiv_ycol(ip)=100
  !     end do
  !     !   missing transverse momentum
  !     o_ETmiss=o_distros
  !     ETmissmax=7000.d0/(1+initial_state*6)
  !     ETmissmin=0.d0
  !     ndiv_ETmiss=70
  !     !   top transverse momentum
  !     o_pT356=o_distros
  !     pT356max=7000.d0/(1+initial_state*6)
  !     pT356min=0.d0
  !     ndiv_pT356=70
  !     !   2to6 top pseudorapidity
  !     o_eta356=o_distros
  !     eta356max=+10
  !     eta356min=-10
  !     ndiv_eta356=50
  !     !   2to6 top pseudorapidity
  !     o_phi356=o_distros
  !     phi356max=+pi
  !     phi356min=-pi
  !     ndiv_phi356=50
  !     !   anti-top transverse momentum
  !     o_pT478=o_distros
  !     pT478max=7000.d0/(1+initial_state*6)
  !     pT478min=0.d0
  !     ndiv_pT478=70
  !     !   2to6 anti-top pseudorapidity
  !     o_eta478=o_distros
  !     eta478max=+10
  !     eta478min=-10
  !     ndiv_eta478=50
  !     !   2to6 top pseudorapidity
  !     o_phi478=o_distros
  !     phi478max=+pi
  !     phi478min=-pi
  !     ndiv_phi478=50
  !     !   invarient mass of tt pair (always on)
  !     o_rMtt=1
  !     rMttmax=14000.d0/(1+initial_state*6)
  !     rMttmin=0.d0
  !     ndiv_rMtt=140
  !     !   boost of parton CoM
  !     o_beta=o_distros
  !     betamax=1000.d0
  !     betamin=0.d0
  !     ndiv_beta=100
  !     !   costheta
  !     o_cost=o_distros
  !     costmax=+1.d0
  !     costmin=-1.d0
  !     ndiv_cost=50
  !     !   top energy
  !     o_Et=o_distros
  !     Etmax=7000.d0/(1+initial_state*6)
  !     Etmin=0.d0
  !     ndiv_Et=70
  !     !   delta_y
  !     o_Delta_y=o_distros
  !     Delta_ymax=4.d0
  !     Delta_ymin=-4.d0
  !     ndiv_Delta_y=100
  !     !   transverse variables
  !     do itrans=1,ntrans
  !       if(ifinal_state == 0)then
  !         o_tran(itrans)=0
  !       else
  !         o_tran(itrans)=o_trans
  !       end if
  !     end do
  !     !   invarient mass of the visible decay products of the tt pair
  !     transmax(1)=4000
  !     transmin(1)=0.d0
  !     ndiv_trans(1)=40
  !     !   sum of tranvserse energy
  !     transmax(2)=4000
  !     transmin(2)=0.d0
  !     ndiv_trans(2)=40
  !     !   transverse mass 1
  !     transmax(3)=4000
  !     transmin(3)=0.d0
  !     ndiv_trans(3)=40
  !     !   transverse mass 2
  !     transmax(4)=4000
  !     transmin(4)=0.d0
  !     ndiv_trans(4)=40
  !     !   transverse mass 3
  !     transmax(5)=4000
  !     transmin(5)=0.d0
  !     ndiv_trans(5)=40
  !     !  lepton transverse mass
  !     transmax(6)=500
  !     transmin(6)=0.d0
  !     ndiv_trans(6)=40
  !     !   contransverse mass 1
  !     transmax(7)=4000
  !     transmin(7)=0.d0
  !     ndiv_trans(7)=40
  !     !   contransverse mass 2
  !     transmax(8)=4000
  !     transmin(8)=0.d0
  !     ndiv_trans(8)=40
  !     !   contransverse mass 3
  !     transmax(9)=4000
  !     transmin(9)=0.d0
  !     ndiv_trans(9)=40
  !     !   lepton contransverse mass
  !     transmax(10)=500
  !     transmin(10)=0.d0
  !     ndiv_trans(10)=50

  !     !   phi_l
  !     o_fl=o_asyms
  !     flmax=+2*pi
  !     flmin=0
  !     ndiv_fl=100
  !     !   cosphi_l
  !     o_cosfl=o_asyms
  !     cosflmax=+1.d0
  !     cosflmin=-1.d0
  !     ndiv_cosfl=100
  !     !   delta phi
  !     o_dphi=o_asyms
  !     dphimax=+pi
  !     dphimin=0
  !     ndiv_dphi=10
  !     !   cost5
  !     o_cost5=o_asyms
  !     cost5max=+1
  !     cost5min=-1
  !     ndiv_cost5=10
  !     !   cost7
  !     o_cost7=o_asyms
  !     cost7max=+1
  !     cost7min=-1
  !     ndiv_cost7=10
  !     !  ct7ct5
  !     o_ct7ct5=o_asyms
  !     ct7ct5max=+1
  !     ct7ct5min=-1
  !     ndiv_ct7ct5=10
  !     !   sigp
  !     o_sigp=o_asyms
  !     sigpmax=rMttmax
  !     sigpmin=rMttmin
  !     ndiv_sigp=ndiv_rMtt/5
  !     !   sigm
  !     o_sigm=o_asyms
  !     sigmmax=rMttmax
  !     sigmmin=rMttmin
  !     ndiv_sigm=ndiv_rMtt/5
  !     !   dphi2d
  !     if((o_dphi == 1) .AND. (o_rMtt == 1))then
  !       o_dphi2d=O_dist2d
  !     else
  !       o_dphi2d=0
  !     end if
  !     !   dtransph
  !     do itrans=1, ntrans
  !       if((o_dphi == 1) .AND. (o_tran(itrans) == 1))then
  !         o_transdp(itrans)=o_dist2d
  !       else
  !         o_transdp(itrans)=0
  !       end if
  !     end do
  !     !   asymmetries
  !     do iasy=1,nasym
  !       o_asym(iasy)=o_asyms
  !     end do

  !     !   Turn off 2->6 only distributions
  !     if (ifinal_state == 0)then
  !       do ip=5,8
  !         o_pT(i)   = 0
  !         o_eta(i)  = 0
  !         o_phi(i)  = 0
  !       end do
  !       do itrans=1,ntrans
  !         o_tran(itrans)=0
  !       end do
  !       o_tran  = 0
  !       o_pT356  = 0
  !       o_eta356 = 0
  !       o_phi356 = 0
  !       o_pT478  = 0
  !       o_eta478 = 0
  !       o_phi478 = 0
  !       o_ETmiss = 0
  !       o_HT     = 0
  !       o_fl     = 0
  !       o_dphi   = 0
  !       o_cosfl  = 0
  !       o_cost7  = 0
  !       o_cost5  = 0
  !       o_ct7ct5 = 0
  !       o_dphi2d = 0
  !       o_asym(9) = 0    ! turn off A_l
  !     end if
  !     !   Turn off 2->2 only distributions
  !     if (ifinal_state >= 1)then
  !       o_asym(1) = 0   ! turn off A_LL
  !       o_asym(2) = 0   ! turn off A_L
  !       o_asym(3) = 0   ! turn off A_PV
  !     end if
  !   end subroutine initialise_distributions
end module Distributions
