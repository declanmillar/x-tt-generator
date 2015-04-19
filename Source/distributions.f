module distributions

  use configuration, only: print_all_distributions, print_2d_distributions, include_transverse, include_asymmetries, o_asym, &
                           n_fb_asymmetries, n_asymmetries, final_state, initial_state, n_final, tops_decay, include_errors
  use mathematics, only: pi
  use kinematics, only: sigma
  use integration, only: it
  use class_histogram

  implicit none

  public

  integer :: o_ptb = 1
  integer :: o_ptbb = 1
  integer :: o_ptlp = 1
  integer :: o_ptlm = 1
  integer :: o_ptnu = 1
  integer :: o_ptnub = 1
  integer :: o_ptt = 1
  integer :: o_pttb = 1
  integer :: o_etab = 1
  integer :: o_etabb = 1
  integer :: o_etalp = 1
  integer :: o_etalm = 1
  integer :: o_etanu = 1
  integer :: o_etanub = 1
  integer :: o_etat = 1
  integer :: o_etatb = 1
  integer :: o_phib = 1
  integer :: o_phibb = 1
  integer :: o_philp = 1
  integer :: o_philm = 1
  integer :: o_phinu = 1
  integer :: o_phinub = 1
  integer :: o_phit   = 1
  integer :: o_phitb = 1
  integer :: o_ycolb = 1
  integer :: o_ycolbb = 1
  integer :: o_ycollp = 1
  integer :: o_ycollm = 1
  integer :: o_ycolnu = 1
  integer :: o_ycolnub = 1
  integer :: o_ycolt   = 1
  integer :: o_ycoltb = 1
  integer :: o_mtt = 1
  integer :: o_mtt_reco = 1
  integer :: o_mtb = 1
  integer :: o_mt_reco = 1
  integer :: o_etmiss = 1
  integer :: o_beta = 1
  integer :: o_cost = 1
  integer :: o_et = 1
  integer :: o_delta_y = 1
  integer :: o_fl = 1
  integer :: o_cosfl = 1
  integer :: o_dphi = 1
  integer :: o_cost5 = 1
  integer :: o_cost7 = 1
  integer :: o_ct7ct5 = 1
  integer :: o_ht = 1
  integer :: o_mttvis = 1
  integer :: o_mt1 = 1
  integer :: o_mt2 = 1
  integer :: o_mt3 = 1
  integer :: o_mct1 = 1
  integer :: o_mct2 = 1
  integer :: o_mct3 = 1
  integer :: o_mlt = 1
  integer :: o_mlct = 1

  type(histogram) :: h_ptb
  type(histogram) :: h_ptbb
  type(histogram) :: h_ptlp
  type(histogram) :: h_ptlm
  type(histogram) :: h_ptnu
  type(histogram) :: h_ptnub
  type(histogram) :: h_ptt
  type(histogram) :: h_pttb
  type(histogram) :: h_etab
  type(histogram) :: h_etabb
  type(histogram) :: h_etalp
  type(histogram) :: h_etalm
  type(histogram) :: h_etanu
  type(histogram) :: h_etanub
  type(histogram) :: h_etat
  type(histogram) :: h_etatb
  type(histogram) :: h_phib
  type(histogram) :: h_phibb
  type(histogram) :: h_philp
  type(histogram) :: h_philm
  type(histogram) :: h_phinu
  type(histogram) :: h_phinub
  type(histogram) :: h_phit  
  type(histogram) :: h_phitb
  type(histogram) :: h_ycolb
  type(histogram) :: h_ycolbb
  type(histogram) :: h_ycollp
  type(histogram) :: h_ycollm
  type(histogram) :: h_ycolnu
  type(histogram) :: h_ycolnub
  type(histogram) :: h_ycolt  
  type(histogram) :: h_ycoltb
  type(histogram) :: h_mtt
  type(histogram) :: h_mtt_reco
  type(histogram) :: h_mtb
  type(histogram) :: h_mt_reco
  type(histogram) :: h_etmiss
  type(histogram) :: h_beta
  type(histogram) :: h_cost
  type(histogram) :: h_et
  type(histogram) :: h_delta_y
  type(histogram) :: h_fl
  type(histogram) :: h_cosfl
  type(histogram) :: h_dphi
  type(histogram) :: h_cost5
  type(histogram) :: h_cost7
  type(histogram) :: h_ct7ct5
  type(histogram) :: h_ht
  type(histogram) :: h_mttvis
  type(histogram) :: h_mt1
  type(histogram) :: h_mt2
  type(histogram) :: h_mt3
  type(histogram) :: h_mct1
  type(histogram) :: h_mct2
  type(histogram) :: h_mct3
  type(histogram) :: h_mlt
  type(histogram) :: h_mlct

  ! distribution in sigp
  real :: sigpmax,sigpmin,sigpw
  real :: xsigp(1000),fxsigp(n_asymmetries,1000,20) ,fxsigptot(n_asymmetries,1000)
  integer :: o_sigp
  integer :: ndiv_sigp

  ! distribution in sigm
  real :: sigmmax,sigmmin,sigmw
  real :: xsigm(1000),fxsigm(n_asymmetries,1000,20) ,fxsigmtot(n_asymmetries,1000)
  integer :: o_sigm
  integer :: ndiv_sigm

  ! 2d dist
  real :: xdphi2d(500,500),fxdphi2d(500,500,20), fxdphi2dtot(500,500)
  integer :: o_dphi2d

  ! 2d dist
  real :: xct7ct52d(500,500),fxct7ct52d(500,500,20), fxct7ct52dtot(500,500)
  integer :: o_ct7ct52d

  ! 2d dist
  real :: xcost72d(500,500),fxcost72d(500,500,20), fxcost72dtot(500,500)
  integer :: o_cost72d

  ! 2d dist
  real :: xcost52d(500,500),fxcost52d(500,500,20), fxcost52dtot(500,500)
  integer :: o_cost52d

!   integer :: include_transversedp(ntrans)

  public :: create_distributions
  public :: initialise_distributions
  public :: generate_bins
  public :: finalise_distributions

  
!   real :: cnorm(20)
  real :: atot(n_asymmetries),atoterr(n_asymmetries)
  integer, private :: i, j, k, ip, iasy, jasy
  real, private :: sfxpttot(8), sfxetatot(8), sfxphitot(8), sfxycoltot(8)
  real, private :: sfxsigptot(n_asymmetries), sfxsigmtot(n_asymmetries), asym_int(n_asymmetries)

contains

subroutine create_distributions

  implicit none 

  print*, "Initialising distributions"

  h_ptb = histogram("pTb", "d#sigma-/dp_{T}^{b}", "p_{T}^{b}", 0.d0, 1000.d0, 100)
  h_ptbb = histogram("pTbb", "d#sigma-/dp_{T}^{#bar{b}}", "p_{T}^{#bar{b}}", 0.d0, 1000.d0, 100)
  h_ptlp = histogram("pTlp", "d#sigma-/dp_{T}^{l^{+}}", "p_{T}^{l^{+}}", 0.d0, 1000.d0, 100)
  h_ptlm = histogram("pTlm", "d#sigma-/dp_{T}^{l^{-}}", "p_{T}^{l^{-}}", 0.d0, 1000.d0, 100)
  h_ptnu = histogram("pTnu", "d#sigma-/dp_{T}^{#nu}", "p_{T}^{#nu}", 0.d0, 1000.d0, 100)
  h_ptnub = histogram("pTnub", "d#sigma-/dp_{T}^{#bar{#nu}}", "p_{T}^{#bar{#nu}}", 0.d0, 1000.d0, 100)
  h_ptt = histogram("pTt", "d#sigma-/dp^{t}_{T}--[pb]", "p^{t}_{T}", 0.d0, 1000.d0, 100)
  h_pttb = histogram("pTtb", "d#sigma-/dp^{#bar{t}}_{T}--[pb]", "p^{#bar{t}}_{T}", 0.d0, 1000.d0, 100)

  h_etab = histogram("etab", "d#sigma-/d#eta_{b}", "#eta_{b}", -10.d0, 10.d0, 100)
  h_etabb = histogram("etabb", "d#sigma-/d#eta_{#bar{b}}", "#eta_{#bar{b}}", -10.d0, 10.d0, 100)
  h_etalp = histogram("etalp", "d#sigma-/d#eta_{l^{+}}", "#eta_{l^{+}}", -10.d0, 10.d0, 100)
  h_etalm = histogram("etalm", "d#sigma-/d#eta_{l^{-}}", "#eta_{l^{-}}", -10.d0, 10.d0, 100)
  h_etanu = histogram("etanu", "d#sigma-/d#eta_{#nu}", "#eta_{#nu}", -10.d0, 10.d0, 100)
  h_etanub = histogram("etanub", "d#sigma-/d#eta_{#bar{#nu}}", "#eta_{#bar{#nu}}", -10.d0, 10.d0, 100)
  h_etat = histogram("etat", "d#sigma-/dd#eta_{t}--[pb]", "#eta_{t}", -10.d0, 10.d0, 100)
  h_etatb = histogram("etatb", "d#sigma-/d#eta_{#bar{t}}--[pb]", "#eta_{#bar{t}}", -10.d0, 10.d0, 100)

  h_phib = histogram("phib", "d#sigma-/d#phi_{b}", "#phi_{b}", -pi, pi, 100)
  h_phibb = histogram("phibb", "d#sigma-/d#phi_{#bar{b}}", "#phi_{#bar{b}}", -pi, pi, 100)
  h_philp = histogram("philp", "d#sigma-/d#phi_{l^{+}}", "#phi_{l^{+}}", -pi, pi, 100)
  h_philm = histogram("philm", "d#sigma-/d#phi_{l^{-}}", "#phi_{l^{-}}", -pi, pi, 100)
  h_phinu = histogram("phinu", "d#sigma-/d#phi_{#nu}", "#phi_{#nu}", -pi, pi, 100)
  h_phinub = histogram("phinub", "d#sigma-/d#phi_{#bar{#nu}}", "#phi_{#bar{#nu}}", -pi, pi, 100)
  h_phit = histogram("phit", "d#sigma-/dd#phi_{t}--[pb]", "#phi_{t}", -pi, pi, 100)
  h_phitb = histogram("phitb", "d#sigma-/d#phi_{#bar{t}}--[pb]", "#phi_{#bar{t}}", -pi, pi, 100)

  h_ycolb = histogram("ycolb", "d#sigma-/dy^b}", "#phi_{b}", -10.d0, 10.d0, 100)
  h_ycolbb = histogram("ycolbb", "d#sigma-/dy^{#bar{b}}_{col}--[pb]", "#phi_{#bar{b}}", -10.d0, 10.d0, 100)
  h_ycollp = histogram("ycollp", "d#sigma-/dy^{l^{+}}_{col}--[pb]", "#phi_{l^{+}}", -10.d0, 10.d0, 100)
  h_ycollm = histogram("ycollm", "d#sigma-/dy^{l^{-}}_{col}--[pb]", "#phi_{l^{-}}", -10.d0, 10.d0, 100)
  h_ycolnu = histogram("ycolnu", "d#sigma-/dy^{#nu}_{col}--[pb]", "#phi_{#nu}", -10.d0, 10.d0, 100)
  h_ycolnub = histogram("ycolnub", "d#sigma-/y^{#bar{#nu}}_{col}--[pb]", "#phi_{#bar{#nu}}", -10.d0, 10.d0, 100)
  h_ycolt = histogram("ycolt", "d#sigma-/dy^{t}_{col}--[pb]", "#phi_{t}", -10.d0, 10.d0, 100)
  h_ycoltb = histogram("ycoltb", "d#sigma-/dy^{#bar{t}}_{col}--[pb]", "#phi_{#bar{t}}", -10.d0, 10.d0, 100)

  h_etmiss = histogram("ETmiss", "d#sigma-/dE_{tt}--[pb/GeV]", 'M_{tt}--[GeV]', 0.d0, 1000.d0, 140)
  h_mtt = histogram("Mtt", "d#sigma-/dM_{tt}--[pb/GeV]", 'M_{tt}--[GeV]', 0.d0, 14000.d0, 140)
  h_mtt_reco = histogram("Mtt_reco", "d#sigma-/dM_{tt}--[pb/GeV]", 'M_{tt}--[GeV]', 0.d0, 14000.d0, 140)
  h_beta = histogram("beta", "d#sigma-/d#beta_{tt}--[pb/GeV]", 'beta_{t}', 0.d0, 1000.d0, 100)
  h_mtb = histogram("mtb", "d#sigma-/d#m_{#bar{t}}--[pb]", "#m_{#bar{t}}", 0.d0, 1000.d0, 100)
  h_mt_reco = histogram("mt_reco", "d#sigma-/dm^{reco}_{t}--[pb]", "m^{reco}_{t}", 0.d0, 1000.d0, 100)
  h_beta = histogram("beta", "d#sigma-/d#beta--[pb]", "#beta", 0.d0, 1000.d0, 100)
  h_cost = histogram("cost", "d#sigma-/dcos#theta_{t}--[pb]", "cos#theta_{t}", -1.d0, 1000.d0, 100)
  h_et = histogram("Et","d#sigma-/dE_{t}--[pb]", "E_{t}", 0.d0, 1000.d0, 100)
  h_delta_y = histogram("delta_y", "d#sigma-/d#Delta-y--[pb]", "#Delta-y", -4.d0, 4.d0, 100)
  h_fl = histogram("phil", "d#sigma-/d#phi_l--[pb]", "#phi_l", 0.d0, 1000.d0, 100)
  h_cosfl = histogram("cosphil", "d#sigma-/dcos#phi_l--[pb]", "cos#phi_l", -1.d0, 1.d0, 100)
  h_dphi = histogram("dphi", "d#sigma-/d#delta #phi_l--[pb]", "#delta #phi_l", 0.d0, 2*pi, 100)
  h_cost5 = histogram("cost5", "d#sigma-/dcos#theta_{l^+}--[pb]", "cos#theta_{l^+}", -1.d0, 1.d0, 100)
  h_cost7 = histogram("cost7", "d#sigma-/dcos#theta_{l^-}--[pb]", "cos#theta_{l^-}", -1.d0, 1.d0, 100)
  h_ct7ct5 = histogram("ct7ct5","d#sigma-/dcos#theta_{l^+}cos#theta_{l^-}--[pb]","cos#theta_{l^+}cos#theta_{l^-}",-1.d0,1.d0,100)
  h_ht = histogram("HT", "d#sigma-/dH_{T}--[pb/GeV]", "H_{T}", 0.d0, 4000.d0, 40)
  h_mttvis = histogram("mttvis", "d#sigma-/dM_{tt}^{vis}--[pb/GeV]", "M_{tt}", 0.d0, 4000.d0, 40)
  h_mt1 = histogram("MT1", "d#sigma-/dM_{T1}--[pb/GeV]", "M_{T1}", 0.d0, 4000.d0, 40)
  h_mt2 = histogram("MT2", "d#sigma-/dM_{T2}--[pb/GeV]", "M_{T2}", 0.d0, 4000.d0, 40)
  h_mt3 = histogram("MT3", "d#sigma-/dM_{T3}--[pb/GeV]", "M_{T3}", 0.d0, 4000.d0, 40)
  h_mct1 = histogram("MCT1", "d#sigma-/dM_{CT1}--[pb/GeV]", "M_{CT}", 0.d0, 4000.d0, 40)
  h_mct2 = histogram("MCT2", "d#sigma-/dM_{CT2}--[pb/GeV]", "M_{CT}", 0.d0, 4000.d0, 40)
  h_mct3 = histogram("MCT3", "d#sigma-/dM_{CT3}--[pb/GeV]", "M_{CT}", 0.d0, 4000.d0, 40)
  h_mlt = histogram("MlT", "d#sigma-/dM^{l}_{T}--[pb/GeV]", "M^{l}_{T}}", 0.d0, 500.d0, 50)
  h_mlct = histogram("MlCT", "d#sigma-/dM^{l}_{CT}--[pb/GeV]", "M^{l}_{CT}}", 0.d0, 500.d0, 50)

end subroutine create_distributions

subroutine initialise_distributions

  implicit none

  if (o_ptb == 1) call h_ptb%initialise()
  if (o_ptbb == 1) call h_ptbb%initialise()
  if (o_ptlp == 1) call h_ptlp%initialise()
  if (o_ptlm == 1) call h_ptlm%initialise()
  if (o_ptnu == 1) call h_ptnu%initialise()
  if (o_ptnub == 1) call h_ptnub%initialise()
  if (o_ptt == 1) call h_ptt%initialise()
  if (o_pttb == 1) call h_pttb%initialise()
  if (o_etab == 1) call h_etab%initialise()
  if (o_etabb == 1) call h_etabb%initialise()
  if (o_etalp == 1) call h_etalp%initialise()
  if (o_etalm == 1) call h_etalm%initialise()
  if (o_etanu == 1) call h_etanu%initialise()
  if (o_etanub == 1) call h_etanub%initialise()
  if (o_etat == 1) call h_etat%initialise()
  if (o_etatb == 1) call h_etatb%initialise()
  if (o_phib == 1) call h_phib%initialise()
  if (o_phibb == 1) call h_phibb%initialise()
  if (o_philp == 1) call h_philp%initialise()
  if (o_philm == 1) call h_philm%initialise()
  if (o_phinu == 1) call h_phinu%initialise()
  if (o_phinub == 1) call h_phinub%initialise()
  if (o_phit == 1) call h_phit%initialise()
  if (o_phitb == 1) call h_phitb%initialise()
  if (o_ycolb == 1) call h_ycolb%initialise()
  if (o_ycolbb == 1) call h_ycolbb%initialise()
  if (o_ycollp == 1) call h_ycollp%initialise()
  if (o_ycollm == 1) call h_ycollm%initialise()
  if (o_ycolnu == 1) call h_ycolnu%initialise()
  if (o_ycolnub == 1) call h_ycolnub%initialise()
  if (o_ycolt == 1) call h_ycolt%initialise()
  if (o_ycoltb == 1) call h_ycoltb%initialise()
  if (o_etmiss == 1) call h_etmiss%initialise()
  if (o_mtb == 1) call h_mtb%initialise()
  if (o_mt_reco == 1) call h_mt_reco%initialise()
  if (o_beta == 1) call h_beta%initialise()
  if (o_mtt == 1) call h_mtt%initialise
  if (o_mtt_reco == 1) call h_mtt_reco%initialise()
  if (o_beta == 1) call h_beta%initialise()
  if (o_cost == 1) call h_cost%initialise()
  if (o_et == 1) call h_et%initialise()
  if (o_delta_y == 1) call h_delta_y%initialise()
  if (o_fl == 1) call h_fl%initialise()
  if (o_cosfl == 1) call h_cosfl%initialise()
  if (o_dphi == 1) call h_dphi%initialise()
  if (o_cost5 == 1) call h_cost5%initialise()
  if (o_cost7 == 1) call h_cost7%initialise()
  if (o_ct7ct5 == 1) call h_ct7ct5%initialise()
  if (o_ht == 1) call h_ht%initialise()
  if (o_mttvis == 1) call h_mttvis%initialise()
  if (o_mt1 == 1) call h_mt1%initialise()
  if (o_mt2 == 1) call h_mt2%initialise()
  if (o_mt3 == 1) call h_mt3%initialise()
  if (o_mct1 == 1) call h_mct1%initialise()
  if (o_mct2 == 1) call h_mct2%initialise()
  if (o_mct3 == 1) call h_mct3%initialise()
  if (o_mlt == 1) call h_mlt%initialise()
  if (o_mlct == 1) call h_mlct%initialise()

  ! sigp
  o_sigp=include_asymmetries
  sigpmax=14000
  sigpmin=0
  ndiv_sigp=100/5
  ! sigm
  o_sigm=include_asymmetries
  sigmmax=14000
  sigmmin=0
  ndiv_sigm=100/5
!     ! dphi2d
!     if((o_dphi == 1) .and. (o_mtt == 1))then
!       o_dphi2d=print_2d_distributions
!     else
!       o_dphi2d=0
!     end if
!     ! dtransph
!     do itrans=1, ntrans
!       if((o_dphi == 1) .and. (o_tran(itrans) == 1))then
!         include_transversedp(itrans)=print_2d_distributions
!       else
!         include_transversedp(itrans)=0
!       end if
!     end do
!     ! asymmetries
  do iasy=1,n_asymmetries
    o_asym(iasy)=include_asymmetries
  end do

end subroutine initialise_distributions

subroutine disable_distributions

  if (initial_state == 0) then
    ! disable non-useful variables for pp
    o_asym(4) = 0
    o_asym(7) = 0
    o_asym(8) = 0 
  end if

  if (initial_state == 1) then
    ! disable non-useful variables for ppbar
    o_asym(5) = 0
    o_asym(6) = 0
  end if

  if (final_state == 0) then
    ! disable 2to6 variables 
    o_ptb = 0
    o_ptbb = 0
    o_ptlp = 0
    o_ptlm = 0
    o_ptnu = 0
    o_ptnub = 0
    o_etab = 0
    o_etabb = 0
    o_etalp = 0
    o_etalm = 0
    o_etanu = 0
    o_etanub = 0
    o_phib = 0
    o_phibb = 0
    o_philp = 0
    o_philm = 0
    o_phinu = 0
    o_phinub = 0
    o_ycolb = 0
    o_ycolbb = 0
    o_ycollp = 0
    o_ycollm = 0
    o_ycolnu = 0
    o_ycolnub = 0
    o_etmiss = 0
    o_fl = 0
    o_dphi = 0
    o_cosfl = 0
    o_cost7 = 0
    o_cost5 = 0
    o_ct7ct5 = 0
    o_dphi2d = 0
    o_ct7ct52d = 0
    o_cost72d = 0
    o_cost52d = 0
    o_asym(6) = 0
    o_asym(10) = 0
    o_asym(11) = 0
    o_asym(12) = 0
    o_mtt_reco = 0
    o_mt_reco = 0
    o_mtb = 0
    o_ht = 0
    o_mttvis = 0
    o_mt1 = 0
    o_mt2 = 0
    o_mt3 = 0
    o_mct1 = 0
    o_mct2 = 0
    o_mct3 = 0
    o_mlt = 0
    o_mlct = 0
  end if

  if (final_state > 0) then
    ! disable non 2to6 variables
    o_asym(1) = 0
    o_asym(2) = 0
    o_asym(3) = 0
  end if

  if (final_state == 1) then 
    ! disable non-useful variables in dileptonic
    o_mtt_reco = 0 
    o_mt_reco = 0
    o_mtb = 0
    do i = 4, 10
      o_asym(i) = 0
    end do
  end if

  if (final_state == 2) then
    ! disable non-useful variables in semi-leptonic
    o_cost7 = 0
    o_ct7ct5 = 0
    o_dphi = 0
    o_dphi2d = 0
    o_ct7ct52d = 0
    o_cost72d = 0
    o_cost52d = 0
    o_etmiss = 0
    o_ht = 0
    o_mttvis = 0
    o_mt1 = 0
    o_mt2 = 0
    o_mt3 = 0
    o_mct1 = 0
    o_mct2 = 0
    o_mct3 = 0
    o_mlt = 0
    o_mlct = 0
    o_asym(11) = 0
  end if

  if (final_state == 3) then
    ! disable non-useful variables in fully hadronic
    o_mtt_reco = 0
    o_cost5 = 0
    o_cost7 = 0
    o_ct7ct5 = 0
    o_dphi = 0
    o_dphi2d = 0
    o_ct7ct52d = 0
    o_cost72d = 0
    o_cost52d = 0
    o_etmiss = 0
    o_mt_reco = 0
    o_mtb = 0
    o_ht = 0
    o_mttvis = 0
    o_mt1 = 0
    o_mt2 = 0
    o_mt3 = 0
    o_mct1 = 0
    o_mct2 = 0
    o_mct3 = 0
    o_mlt = 0
    o_mlct = 0

    o_asym(6) = 0
    o_asym(10) = 0
    o_asym(11) = 0
    o_asym(12) = 0
  end if

  ! disable A_PV
  o_asym(3) = 0

end subroutine disable_distributions

subroutine generate_bins
  ! (finds bin width, finds midpoints.)
  implicit none

  print*, "Generating bins"

  if(o_sigp == 1)then
    sigpw=(sigpmax-sigpmin)/ndiv_sigp
    do i=1,ndiv_sigp
      xsigp(i)=sigpmin+sigpw*(i-1)+sigpw/2.d0
    end do
  end if

  if(o_sigm == 1)then
    sigmw=(sigmmax-sigmmin)/ndiv_sigm
    do i=1,ndiv_sigm
      xsigm(i)=sigmmin+sigmw*(i-1)+sigmw/2.d0
    end do
  end if

end subroutine generate_bins  

subroutine finalise_distributions

  implicit none

  integer :: ndiv_sig

  print*, "Printing histograms"

  write(10,*) '-----------------------------------------------------'
  write(10,*)'HISTOGRAMS'
    
  if (o_ptb == 1) call h_ptb%finalise()
  if (o_ptbb == 1) call h_ptbb%finalise()
  if (o_ptlp == 1) call h_ptlp%finalise()
  if (o_ptlm == 1) call h_ptlm%finalise()
  if (o_ptnu == 1) call h_ptnu%finalise()
  if (o_ptnub == 1) call h_ptnub%finalise()
  if (o_ptt == 1) call h_ptt%finalise()
  if (o_pttb == 1) call h_pttb%finalise()
  if (o_etab == 1) call h_etab%finalise()
  if (o_etabb == 1) call h_etabb%finalise()
  if (o_etalp == 1) call h_etalp%finalise()
  if (o_etalm == 1) call h_etalm%finalise()
  if (o_etanu == 1) call h_etanu%finalise()
  if (o_etanub == 1) call h_etanub%finalise()
  if (o_etat == 1) call h_etat%finalise()
  if (o_etatb == 1) call h_etatb%finalise()
  if (o_phib == 1) call h_phib%finalise()
  if (o_phibb == 1) call h_phibb%finalise()
  if (o_philp == 1) call h_philp%finalise()
  if (o_philm == 1) call h_philm%finalise()
  if (o_phinu == 1) call h_phinu%finalise()
  if (o_phinub == 1) call h_phinub%finalise()
  if (o_phit == 1) call h_phit%finalise()
  if (o_phitb == 1) call h_phitb%finalise()
  if (o_ycolb == 1) call h_ycolb%finalise()
  if (o_ycolbb == 1) call h_ycolbb%finalise()
  if (o_ycollp == 1) call h_ycollp%finalise()
  if (o_ycollm == 1) call h_ycollm%finalise()
  if (o_ycolnu == 1) call h_ycolnu%finalise()
  if (o_ycolnub == 1) call h_ycolnub%finalise()
  if (o_ycolt == 1) call h_ycolt%finalise()
  if (o_ycoltb == 1) call h_ycoltb%finalise()
  if (o_etmiss == 1) call h_etmiss % finalise()
  if (o_mtt == 1) call h_mtt % finalise()
  if (o_mtt_reco == 1) call h_mtt_reco % finalise()
  if (o_mtb == 1) call h_mtb % finalise()
  if (o_mt_reco == 1) call h_mt_reco % finalise()
  if (o_beta == 1) call h_beta % finalise()
  if (o_cost == 1) call h_cost % finalise()
  if (o_et == 1) call h_et % finalise()
  if (o_delta_y == 1) call h_delta_y % finalise()
  if (o_fl == 1) call h_fl%finalise()
  if (o_cosfl == 1) call h_cosfl%finalise()
  if (o_dphi == 1) call h_dphi%finalise()
  if (o_cost5 == 1) call h_cost5%finalise()
  if (o_cost7 == 1) call h_cost7%finalise()
  if (o_ct7ct5 == 1) call h_ct7ct5%finalise()
  if (o_ht == 1) call h_ht%finalise()
  if (o_mttvis == 1) call h_mttvis%finalise()
  if (o_mt1 == 1) call h_mt1%finalise()
  if (o_mt2 == 1) call h_mt2%finalise()
  if (o_mt3 == 1) call h_mt3%finalise()
  if (o_mct1 == 1) call h_mct1%finalise()
  if (o_mct2 == 1) call h_mct2%finalise()
  if (o_mct3 == 1) call h_mct3%finalise()
  if (o_mlt == 1) call h_mlt%finalise()
  if (o_mlct == 1) call h_mlct%finalise()

!     ! plot distributions in all asymmetries
!     if((o_sigp == 1) .and. (o_sigm == 1))then
!       do jasy=1,n_asymmetries
!         if(o_asym(jasy) == 0)then
!           continue
!         else
!           ! snorm(jasy)=0.d0
!           sfxsigptot(jasy)=0d0
!           do j=1,ndiv_sigp
!             fxsigptot(jasy,j)=0.d0
!             do i=1,it
!               fxsigp(jasy,j,i)=fxsigp(jasy,j,i)*sigma/cnorm(i)/sigpw
!               fxsigptot(jasy,j)=fxsigptot(jasy,j)+fxsigp(jasy,j,i)
!             end do
!             sfxsigptot(jasy)=sfxsigptot(jasy)+fxsigptot(jasy,j)*sigpw
!           end do
!           sfxsigmtot(jasy)=0d0
!           do j=1,ndiv_sigm
!             do i=1,it
!               fxsigm(jasy,j,i)=fxsigm(jasy,j,i)*sigma/cnorm(i)/sigmw
!               fxsigmtot(jasy,j)=fxsigmtot(jasy,j)+fxsigm(jasy,j,i)
!             end do
!             sfxsigmtot(jasy)=sfxsigmtot(jasy)+fxsigmtot(jasy,j)*sigmw
!           end do
!           write(10,*)'ASYMMETRY'
!           if(jasy == 1)then
!             write(10,*)'ALL'
!             write(10,*)'A_{ll}'
!           else if(jasy == 2)then
!             write(10,*)'AL'
!             write(10,*)'A_{l}'
!           else if(jasy == 3)then
!             write(10,*)'APV'
!             write(10,*)'A_{pv}'
!           else if(jasy == 4)then
!             write(10,*)'AFB'
!             write(10,*)'A_{fb}'
!           else if(jasy == 5)then
!             write(10,*)'AFBSTAR'
!             write(10,*)'A_{fb^{*}}'
!           else if(jasy == 6)then
!             write(10,*)'AFBSTAR_reco'
!             write(10,*)'A_{fb^{*}}(reco)'
!           else if(jasy == 7)then
!             write(10,*)'AtRFB'
!             write(10,*)'a^{t}_{rfb}'
!           else if(jasy == 8)then
!             write(10,*)'AttbRFB'
!             write(10,*)"a^{b\bar{b}}_{rfb}"
!           else if(jasy == 9)then
!             write(10,*)'ARFB'
!             write(10,*)"A_{rfb}"
!            else if(jasy == 10)then
!             write(10,*)'ARFB_reco'
!             write(10,*)"A_{rfb}(reco)"
!           else if(jasy == 11)then
!             write(10,*)'A_l'
!             write(10,*)'A_{l^+}'
!           else if(jasy == 12)then
!             write(10,*)'AlFB'
!             write(10,*)'A^{l^{+}}_{FB}'
!           end if
!           write(10,*)'M_{tt}'
!           ndiv_sig=(ndiv_sigm+ndiv_sigp)/2
!           do i=1,ndiv_sig
!             if(fxsigptot(jasy,i)+fxsigmtot(jasy,i) == 0.d0)then
!               write(10,*)(xsigm(i)+xsigp(i))/2.d0,0.d0,0.d0,0.d0
!               !           snorm(jasy)=snorm(jasy)+0.d0
!             else
!               write(10,*)(xsigm(i)+xsigp(i))/2.d0, &
!               (fxsigptot(jasy,i)-fxsigmtot(jasy,i))/ &
!               (fxsigptot(jasy,i)+fxsigmtot(jasy,i)), &
!               fxsigptot(jasy,i),fxsigmtot(jasy,i)
!               !             snorm(jasy)=snorm(jasy)+
!               !    &               (fxsigptot(jasy,i)-fxsigmtot(jasy,i))/
!               !    &               (fxsigptot(jasy,i)+fxsigmtot(jasy,i))
!               !    &               *fxmtttot(i)*mttw/sigma
!             end if
!           end do
!           asym_int(jasy)=(sfxsigptot(jasy)-sfxsigmtot(jasy))/(sfxsigptot(jasy)+sfxsigmtot(jasy))
!           write(10,*)'END'
!           !         write(10,*)'(total asymmetry:',asym_int(jasy),')'
!           !         write(10,*)'(integrated asymmetry:',snorm(jasy),' )'
!         end if
!       end do
!     end if

!     ! plot 2d-distribution in delta phi
!     if(o_dphi2d == 1)then
! !       sfxdphi2dtot=0d0
!       do i=1,ndiv_dphi
!         do j=1,ndiv_mtt
!           fxdphi2dtot(i,j)=0.d0
!           do k=1,it
!             fxdphi2d(i,j,k)=fxdphi2d(i,j,k)*sigma/cnorm(k)/dphiw/mttw
!             fxdphi2dtot(i,j)=fxdphi2dtot(i,j)+fxdphi2d(i,j,k)
!           end do
! !           sfxdphi2dtot=sfxdphitot+fxdphitot(j)*dphiw
!         end do
!       end do
!       write(10,*)'2D-DISTRIBUTION'
!       write(10,*)'dphi2d'
!       write(10,*)'d^2#sigma-/d#delta#phi-dM_{tt}--[pb/GeV]'
!       write(10,*)'#delta#phi--[rad]'
!       write(10,*) ndiv_dphi
!       write(10,*) dphimin
!       write(10,*) dphimax
!       write(10,*)'M_{tt}--[GeV]'
!       write(10,*) ndiv_mtt
!       write(10,*) mttmin
!       write(10,*) mttmax
!       do i=1,ndiv_dphi
!         do j=1,ndiv_mtt
!           write(10,*)xdphi(i),xmtt(j),fxdphi2dtot(i,j)
!         end do
!       end do
!       write(10,*)'END'
!     end if

!     ! plot 2d-distribution in ct7ct5
!     if(o_ct7ct52d == 1)then
! !       sfxct7ct52dtot=0d0
!       do i=1,ndiv_ct7ct5
!         do j=1,ndiv_mtt
!           fxct7ct52dtot(i,j)=0.d0
!           do k=1,it
!             fxct7ct52d(i,j,k)=fxct7ct52d(i,j,k)*sigma/cnorm(k)/ct7ct5w/mttw
!             fxct7ct52dtot(i,j)=fxct7ct52dtot(i,j)+fxct7ct52d(i,j,k)
!           end do
! !           sfxdphitot=sfxdphitot+fxdphitot(j)*dphiw
!         end do
!       end do
!       write(10,*)'2D-DISTRIBUTION'
!       write(10,*)'ct7ct52d'
!       write(10,*)'d^{2}#sigma-/d(cos#theta^{*}_{+}cos#theta^{*}_{-})--[pb]'
!       write(10,*)'cos#theta_{+}cos#theta_{-}'
!       write(10,*) ndiv_ct7ct5
!       write(10,*) ct7ct5min
!       write(10,*) ct7ct5max
!       write(10,*)'M_{tt}--[GeV]'
!       write(10,*) ndiv_mtt
!       write(10,*) mttmin
!       write(10,*) mttmax
!       do i=1,ndiv_ct7ct5
!         do j=1,ndiv_mtt
!           write(10,*)xct7ct5(i),xmtt(j),fxct7ct52dtot(i,j)
!         end do
!       end do
!       write(10,*)'END'
!     end if

!     ! plot 2d-distribution in cost7
!     if(o_cost72d == 1)then
! !       sfxcost72dtot=0d0
!       do i=1,ndiv_cost7
!         do j=1,ndiv_mtt
!           fxcost72dtot(i,j)=0.d0
!           do k=1,it
!             fxcost72d(i,j,k)=fxcost72d(i,j,k)*sigma/cnorm(k)/cost7w/mttw
!             fxcost72dtot(i,j)=fxcost72dtot(i,j)+fxcost72d(i,j,k)
!           end do
! !           sfxdphitot=sfxdphitot+fxdphitot(j)*dphiw
!         end do
!       end do
!       write(10,*)'2D-DISTRIBUTION'
!       write(10,*)'cost72d'
!       write(10,*)'d^{2}#sigma-/d(cos#theta^{*}_{-})--[pb]'
!       write(10,*)'cos#theta_{-}'
!       write(10,*) ndiv_cost7
!       write(10,*) cost7min
!       write(10,*) cost7max
!       write(10,*)'M_{tt}--[GeV]'
!       write(10,*) ndiv_mtt
!       write(10,*) mttmin
!       write(10,*) mttmax
!       do i=1,ndiv_cost7
!         do j=1,ndiv_mtt
!           write(10,*)xcost7(i),xmtt(j),fxcost72dtot(i,j)
!         end do
!       end do
!       write(10,*)'END'
!     end if

!     ! plot 2d-distribution in cost5
!     if(o_cost52d == 1)then
! !       sfxcost52dtot=0d0
!       do i=1,ndiv_cost5
!         do j=1,ndiv_mtt
!           fxcost52dtot(i,j)=0.d0
!           do k=1,it
!             fxcost52d(i,j,k)=fxcost52d(i,j,k)*sigma/cnorm(k)/cost5w/mttw
!             fxcost52dtot(i,j)=fxcost52dtot(i,j)+fxcost52d(i,j,k)
!           end do
! !           sfxdphitot=sfxdphitot+fxdphitot(j)*dphiw
!         end do
!       end do
!       write(10,*)'2D-DISTRIBUTION'
!       write(10,*)'cost52d'
!       write(10,*)'d^{2}#sigma-/d(cos#theta^{*}_{+}cos#theta^{*}_{-})--[pb]'
!       write(10,*)'cos#theta_{+}cos#theta_{-}'
!       write(10,*) ndiv_cost5
!       write(10,*) cost5min
!       write(10,*) cost5max
!       write(10,*)'M_{tt}--[GeV]'
!       write(10,*) ndiv_mtt
!       write(10,*) mttmin
!       write(10,*) mttmax
!       do i=1,ndiv_cost5
!         do j=1,ndiv_mtt
!           write(10,*)xcost5(i),xmtt(j),fxcost52dtot(i,j)
!         end do
!       end do
!       write(10,*)'END'
!     end if

!     ! plot 2d-distributions in delta_phi and all transverse variables
!     do itrans=1,ntrans
!       if(include_transversedp(itrans) == 1)then
!         sfxtransdptot(itrans)=0d0
!         do i=1,ndiv_ct7ct5
!           do j=1,ndiv_trans(itrans)
!             fxtransdptot(itrans,i,j)=0.d0
!             do k=1,it
!               fxtransdp(itrans,i,j,k)=fxtransdp(itrans,i,j,k) &
!               *sigma/cnorm(k)/transw(itrans)/dphiw
!               fxtransdptot(itrans,i,j)=fxtransdptot(itrans,i,j) &
!               +fxtransdp(itrans,i,j,k)
!             end do
!             sfxtransdptot(itrans)=sfxtransdptot(itrans)+ &
!             fxtransdptot(itrans,i,j)*transw(itrans)
!           end do
!         end do
!         write(10,*)'2D-DISTRIBUTION'
!         if (itrans == 1)then
!           write(10,*)'dphiMvis'
!           write(10,*)'d#sigma-/dM_{vis}--[pb/GeV]'
!         else if (itrans == 2)then
!           write(10,*)'dphiHT'
!           write(10,*)'d#sigma-/dH_{T}--[pb/GeV]'
!         else if (itrans == 3)then
!           write(10,*)'dphiM_T1'
!           write(10,*)'d#sigma-/dM_{CT1}--[pb/GeV]'
!         else if (itrans == 4)then
!           write(10,*)'dphiM_T2'
!           write(10,*)'d#sigma-/dM_{CT2}--[pb/GeV]'
!         else if (itrans == 5)then
!           write(10,*)'dphiM_T3'
!           write(10,*)'d#sigma-/dM_{CT3}--[pb/GeV]'
!         else if (itrans == 6)then
!           write(10,*)'dphiMlT'
!           write(10,*)'d#sigma-/dM^{l}_{T}--[pb/GeV]'
!         else if (itrans == 7)then
!           write(10,*)'dphiM_CT1'
!           write(10,*)'d#sigma-/dM_{CT1}--[pb/GeV]'
!         else if (itrans == 8)then
!           write(10,*)'dphiM_CT2'
!           write(10,*)'d#sigma-/dM_{CT2}--[pb/GeV]'
!         else if (itrans == 9)then
!           write(10,*)'dphiM_CT3'
!           write(10,*)'d#sigma-/dM_{CT3}--[pb/GeV]'
!         else if (itrans == 10)then
!           write(10,*)'dphiMlCT'
!           write(10,*)'d#sigma-/dM^{l}_{CT}--[pb/GeV]'
!         else
!           continue
!         end if
!         write(10,*)'#delta#phi'
!         write(10,*) ndiv_dphi
!         write(10,*) dphimin
!         write(10,*) dphimax
!         if (itrans == 1)then
!           write(10,*)'M_{vis}--[GeV]'
!         else if (itrans == 2)then
!           write(10,*)'H_{T}--[GeV]'
!         else if (itrans == 3)then
!           write(10,*)'M_{T1}--[GeV]'
!         else if (itrans == 4)then
!           write(10,*)'M_{T2}--[GeV]'
!         else if (itrans == 5)then
!           write(10,*)'M_{T3}--[GeV]'
!         else if (itrans == 6)then
!           write(10,*)'M^{l}_{T}--[GeV]'
!         else if (itrans == 7)then
!           write(10,*)'M_{T1}--[GeV]'
!         else if (itrans == 8)then
!           write(10,*)'M_{T2}--[GeV]'
!         else if (itrans == 9)then
!           write(10,*)'M_{T3}--[GeV]'
!         else if (itrans == 10)then
!           write(10,*)'M^{l}_{CT}--[GeV]'
!         else
!           continue
!         end if
!         write(10,*)ndiv_trans(itrans)
!         write(10,*)transmin(itrans)
!         write(10,*)transmax(itrans)
!         do i=1,ndiv_dphi
!           do j=1,ndiv_trans(itrans)
!             write(10,*)xdphi(i),xtrans(itrans,j) &
!             ,fxtransdptot(itrans,i,j)
!           end do
!         end do
!         write(10,*)'END'
!       end if
!     end do

end subroutine finalise_distributions

end module distributions
