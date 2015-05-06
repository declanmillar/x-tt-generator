module distributions

  use configuration
  use mathematics, only: pi
  use scattering, only: sigma
  use integration, only: it, cnorm
  use class_histogram
  use class_histogram2d

  implicit none

  public

! 1d distributions
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

  public :: disable_distributions
  public :: create_distributions
  public :: initialise_distributions
  public :: finalise_distributions

  integer, private :: i, j, k, iasy, jasy
  real, private :: sfxsigptot(n_asymmetries), sfxsigmtot(n_asymmetries), asym_int(n_asymmetries)

contains

subroutine create_distributions

  implicit none 

  print *, "Creating distributions..."

  if (o_ptb == 1) h_ptb = histogram("pTb", "d#sigma-/dp_{T}^{b}", "p_{T}^{b}", 100, 0.d0, 1000.d0)
  if (o_ptbb == 1) h_ptbb = histogram("pTbb", "d#sigma-/dp_{T}^{#bar{b}}", "p_{T}^{#bar{b}}", 100, 0.d0, 1000.d0)
  if (o_ptlp == 1) h_ptlp = histogram("pTlp", "d#sigma-/dp_{T}^{l^{+}}", "p_{T}^{l^{+}}", 100, 0.d0, 1000.d0)
  if (o_ptlm == 1) h_ptlm = histogram("pTlm", "d#sigma-/dp_{T}^{l^{-}}", "p_{T}^{l^{-}}", 100, 0.d0, 1000.d0)
  if (o_ptnu == 1) h_ptnu = histogram("pTnu", "d#sigma-/dp_{T}^{#nu}", "p_{T}^{#nu}", 100, 0.d0, 1000.d0)
  if (o_ptnub == 1) h_ptnub = histogram("pTnub", "d#sigma-/dp_{T}^{#bar{#nu}}", "p_{T}^{#bar{#nu}}", 100, 0.d0, 1000.d0)
  if (o_ptt == 1) h_ptt = histogram("pTt", "d#sigma-/dp^{t}_{T}--[pb]", "p^{t}_{T}", 100, 0.d0, 1000.d0)
  if (o_pttb == 1) h_pttb = histogram("pTtb", "d#sigma-/dp^{#bar{t}}_{T}--[pb]", "p^{#bar{t}}_{T}", 100, 0.d0, 1000.d0)
  if (o_etab == 1) h_etab = histogram("etab", "d#sigma-/d#eta_{b}", "#eta_{b}", 100, -10.d0, 10.d0)
  if (o_etabb == 1) h_etabb = histogram("etabb", "d#sigma-/d#eta_{#bar{b}}", "#eta_{#bar{b}}", 100, -10.d0, 10.d0)
  if (o_etalp == 1) h_etalp = histogram("etalp", "d#sigma-/d#eta_{l^{+}}", "#eta_{l^{+}}", 100, -10.d0, 10.d0)
  if (o_etalm == 1) h_etalm = histogram("etalm", "d#sigma-/d#eta_{l^{-}}", "#eta_{l^{-}}", 100, -10.d0, 10.d0)
  if (o_etanu == 1) h_etanu = histogram("etanu", "d#sigma-/d#eta_{#nu}", "#eta_{#nu}", 100, -10.d0, 10.d0)
  if (o_etanub == 1) h_etanub = histogram("etanub", "d#sigma-/d#eta_{#bar{#nu}}", "#eta_{#bar{#nu}}", 100, -10.d0, 10.d0)
  if (o_etat == 1) h_etat = histogram("etat", "d#sigma-/dd#eta_{t}--[pb]", "#eta_{t}", 100, -10.d0, 10.d0)
  if (o_etatb == 1) h_etatb = histogram("etatb", "d#sigma-/d#eta_{#bar{t}}--[pb]", "#eta_{#bar{t}}", 100, -10.d0, 10.d0)
  if (o_phib == 1) h_phib = histogram("phib", "d#sigma-/d#phi_{b}", "#phi_{b}", 100, -pi, pi)
  if (o_phibb == 1) h_phibb = histogram("phibb", "d#sigma-/d#phi_{#bar{b}}", "#phi_{#bar{b}}", 100,  -pi, pi)
  if (o_philp == 1) h_philp = histogram("philp", "d#sigma-/d#phi_{l^{+}}", "#phi_{l^{+}}", 100, -pi, pi)
  if (o_philm == 1) h_philm = histogram("philm", "d#sigma-/d#phi_{l^{-}}", "#phi_{l^{-}}", 100, -pi, pi)
  if (o_phinu == 1) h_phinu = histogram("phinu", "d#sigma-/d#phi_{#nu}", "#phi_{#nu}", 100, -pi, pi)
  if (o_phinub == 1) h_phinub = histogram("phinub", "d#sigma-/d#phi_{#bar{#nu}}", "#phi_{#bar{#nu}}", 100, -pi, pi)
  if (o_phit == 1) h_phit = histogram("phit", "d#sigma-/dd#phi_{t}--[pb]", "#phi_{t}", 100, -pi, pi)
  if (o_phitb == 1) h_phitb = histogram("phitb", "d#sigma-/d#phi_{#bar{t}}--[pb]", "#phi_{#bar{t}}", 100, -pi, pi)
  if (o_ycolb == 1) h_ycolb = histogram("ycolb", "d#sigma-/dy^b}", "#phi_{b}", 100, -10.d0, 10.d0)
  if (o_ycolbb == 1) h_ycolbb = histogram("ycolbb", "d#sigma-/dy^{#bar{b}}_{col}--[pb]", "#phi_{#bar{b}}", 100, -10.d0, 10.d0)
  if (o_ycollp == 1) h_ycollp = histogram("ycollp", "d#sigma-/dy^{l^{+}}_{col}--[pb]", "#phi_{l^{+}}", 100, -10.d0, 10.d0)
  if (o_ycollm == 1) h_ycollm = histogram("ycollm", "d#sigma-/dy^{l^{-}}_{col}--[pb]", "#phi_{l^{-}}", 100, -10.d0, 10.d0)
  if (o_ycolnu == 1) h_ycolnu = histogram("ycolnu", "d#sigma-/dy^{#nu}_{col}--[pb]", "#phi_{#nu}", 100, -10.d0, 10.d0)
  if (o_ycolnub == 1) h_ycolnub = histogram("ycolnub", "d#sigma-/y^{#bar{#nu}}_{col}--[pb]", "#phi_{#bar{#nu}}",100, -10.d0, 10.d0)
  if (o_ycolt == 1) h_ycolt = histogram("ycolt", "d#sigma-/dy^{t}_{col}--[pb]", "#phi_{t}", 100, -10.d0, 10.d0)
  if (o_ycoltb == 1) h_ycoltb = histogram("ycoltb", "d#sigma-/dy^{#bar{t}}_{col}--[pb]", "#phi_{#bar{t}}", 100, -10.d0, 10.d0)
  if (o_etmiss == 1) h_etmiss = histogram("ETmiss", "d#sigma-/dE_{tt}--[pb/GeV]", 'M_{tt}--[GeV]', 140, 0.d0, 1000.d0)
  if (o_mtb == 1) h_mtt = histogram("Mtt", "d#sigma-/dM_{tt}--[pb/GeV]", 'M_{tt}--[GeV]', 140, 0.d0, 14000.d0)
  if (o_mt_reco == 1) h_mtt_reco = histogram("Mtt_reco", "d#sigma-/dM_{tt}--[pb/GeV]", 'M_{tt}--[GeV]', 140, 0.d0, 14000.d0)
  if (o_beta == 1) h_beta = histogram("beta", "d#sigma-/d#beta_{tt}--[pb/GeV]", 'beta_{t}', 100, 0.d0, 1000.d0)
  if (o_mtb == 1) h_mtb = histogram("mtb", "d#sigma-/d#m_{#bar{t}}--[pb]", "#m_{#bar{t}}", 100, 0.d0, 1000.d0)
  if (o_mtt == 1) h_mtt = histogram("mtt", "d#sigma-/dM_{tt}--[pb]", "M_{tt}", 100, 0.d0, 14000.d0)
  if (o_mtt_reco == 1) h_mt_reco = histogram("mt_reco", "d#sigma-/dm^{reco}_{t}--[pb]", "m^{reco}_{t}", 100, 0.d0, 1000.d0)
  if (o_beta == 1) h_beta = histogram("beta", "d#sigma-/d#beta--[pb]", "#beta", 100, 0.d0, 50.d0)
  if (o_cost == 1) h_cost = histogram("cost", "d#sigma-/dcos#theta_{t}--[pb]", "cos#theta_{t}", 100, -1.d0, 1.d0)
  if (o_et == 1) h_et = histogram("Et","d#sigma-/dE_{t}--[pb]", "E_{t}", 100, 0.d0, 1000.d0)
  if (o_delta_y == 1) h_delta_y = histogram("delta_y", "d#sigma-/d#Delta-y--[pb]", "#Delta-y", 100, -4.d0, 4.d0)
  if (o_fl == 1) h_fl = histogram("phil", "d#sigma-/d#phi_l--[pb]", "#phi_{l}", 100, 0.d0, 1000.d0)
  if (o_cosfl == 1) h_cosfl = histogram("cosphil", "d#sigma-/dcos#phi_{l}--[pb]", "cos#phi_{l}", 100, -1.d0, 1.d0)
  if (o_dphi == 1) h_dphi = histogram("dphi", "d#sigma-/d#delta-#phi_{l}--[pb]", "#delta-#phi_{l}", 100, 0.d0, 2*pi)
  if (o_cost5 == 1) h_cost5 = histogram("cost5", "d#sigma-/dcos#theta_{l^{+}}--[pb]", "cos#theta_{l^{+}}", 100, -1.d0, 1.d0)
  if (o_cost7 == 1) h_cost7 = histogram("cost7", "d#sigma-/dcos#theta_{l^{-}}--[pb]", "cos#theta_{l^{-}}", 100, -1.d0, 1.d0)
  if (o_ct7ct5 == 1) h_ct7ct5 = histogram("ct7ct5","d#sigma-/dcos#theta_{l^{+}}cos#theta_{l^{-}}--[pb]", &
     "cos#theta_{l^+}cos#theta_{l^-}", 100,-1.d0,1.d0)

  ! sigp
  o_sigp = 1
  sigpmax = 14000
  sigpmin = 0
  ndiv_sigp = 100/5
  ! sigm
  o_sigm = 1
  sigmmax = 14000
  sigmmin = 0
  ndiv_sigm = 100/5

  print *, "done."
end subroutine create_distributions

subroutine initialise_distributions

  implicit none

  print*, "Initialising distributions..."

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
  if (o_mtt == 1) call h_mtt%initialise()
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

  print *, "done."

end subroutine initialise_distributions

subroutine disable_distributions

  print *, "Disabling irrelevent distributions..."

  do iasy = 1, n_asymmetries
    o_asym(iasy) = 1
  end do

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
    o_ptt = 0
    o_pttb = 0
    o_ptlp = 0
    o_ptlm = 0
    o_ptnu = 0
    o_ptnub = 0
    o_etat = 0
    o_etatb = 0
    o_etalp = 0
    o_etalm = 0
    o_etanu = 0
    o_etanub = 0
    o_phit = 0
    o_phitb = 0
    o_philp = 0
    o_philm = 0
    o_phinu = 0
    o_phinub = 0
    o_ycolt = 0
    o_ycoltb = 0
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
    o_asym(6) = 0
    o_asym(10) = 0
    o_asym(11) = 0
    o_asym(12) = 0
    o_mtt_reco = 0
    o_mt_reco = 0
    o_mtb = 0
  end if

  if (final_state > 0) then
    ! disable non 2to6 variables
    o_asym(1) = 0
    o_asym(2) = 0
    o_asym(3) = 0
  end if

  if (final_state == 1) then 
    ! disable non-useful variables in dileptonic
    o_mll = 1
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
    o_etmiss = 0
    o_asym(5) = 0
    o_asym(9) = 0
  end if

  if (final_state == 3) then
    ! disable non-useful variables in fully hadronic
    o_mtt_reco = 0
    o_cost5 = 0
    o_cost7 = 0
    o_ct7ct5 = 0
    o_dphi = 0
    o_etmiss = 0
    o_mt_reco = 0
    o_mtb = 0

    o_asym(6) = 0
    o_asym(10) = 0
    o_asym(11) = 0
    o_asym(12) = 0
  end if

  ! disable A_PV
  o_asym(3) = 0

  print *, "done."

end subroutine disable_distributions

subroutine finalise_distributions

  implicit none

  integer :: ndiv_sig

  print*, "Printing histograms..."
  open(unit = 10, file = 'Output/'//output_file, status = "replace", action = "write")

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

  ! plot distributions in all asymmetries
  if((o_sigp == 1) .and. (o_sigm == 1))then
    do jasy=1,n_asymmetries
      if(o_asym(jasy) == 0)then
        continue
      else
        ! snorm(jasy)=0.d0
        sfxsigptot(jasy)=0d0
        do j=1,ndiv_sigp
          fxsigptot(jasy,j)=0.d0
          do i=1,it
            fxsigp(jasy,j,i)=fxsigp(jasy,j,i)*sigma/cnorm(i)/sigpw
            fxsigptot(jasy,j)=fxsigptot(jasy,j)+fxsigp(jasy,j,i)
          end do
          sfxsigptot(jasy)=sfxsigptot(jasy)+fxsigptot(jasy,j)*sigpw
        end do
        sfxsigmtot(jasy)=0d0
        do j=1,ndiv_sigm
          do i=1,it
            fxsigm(jasy,j,i)=fxsigm(jasy,j,i)*sigma/cnorm(i)/sigmw
            fxsigmtot(jasy,j)=fxsigmtot(jasy,j)+fxsigm(jasy,j,i)
          end do
          sfxsigmtot(jasy)=sfxsigmtot(jasy)+fxsigmtot(jasy,j)*sigmw
        end do
        write(10,*)'ASYMMETRY'
        if(jasy == 1)then
          write(10,*)'ALL'
          write(10,*)'A_{ll}'
        else if(jasy == 2)then
          write(10,*)'AL'
          write(10,*)'A_{l}'
        else if(jasy == 3)then
          write(10,*)'APV'
          write(10,*)'A_{pv}'
        else if(jasy == 4)then
          write(10,*)'AFB'
          write(10,*)'A_{fb}'
        else if(jasy == 5)then
          write(10,*)'AFBSTAR'
          write(10,*)'A_{fb^{*}}'
        else if(jasy == 6)then
          write(10,*)'AFBSTAR_reco'
          write(10,*)'A_{fb^{*}}(reco)'
        else if(jasy == 7)then
          write(10,*)'AtRFB'
          write(10,*)'a^{t}_{rfb}'
        else if(jasy == 8)then
          write(10,*)'AttbRFB'
          write(10,*)"a^{b\bar{b}}_{rfb}"
        else if(jasy == 9)then
          write(10,*)'ARFB'
          write(10,*)"A_{rfb}"
         else if(jasy == 10)then
          write(10,*)'ARFB_reco'
          write(10,*)"A_{rfb}(reco)"
        else if(jasy == 11)then
          write(10,*)'A_l'
          write(10,*)'A_{l^+}'
        else if(jasy == 12)then
          write(10,*)'AlFB'
          write(10,*)'A^{l^{+}}_{FB}'
        end if
        write(10,*)'M_{tt}'
        ndiv_sig=(ndiv_sigm+ndiv_sigp)/2
        do i=1,ndiv_sig
          if(fxsigptot(jasy,i)+fxsigmtot(jasy,i) == 0.d0)then
            write(10,*)(xsigm(i)+xsigp(i))/2.d0,0.d0,0.d0,0.d0
            !           snorm(jasy)=snorm(jasy)+0.d0
          else
            write(10,*)(xsigm(i)+xsigp(i))/2.d0, &
            (fxsigptot(jasy,i)-fxsigmtot(jasy,i))/ &
            (fxsigptot(jasy,i)+fxsigmtot(jasy,i)), &
            fxsigptot(jasy,i),fxsigmtot(jasy,i)
            !             snorm(jasy)=snorm(jasy)+
            !    &               (fxsigptot(jasy,i)-fxsigmtot(jasy,i))/
            !    &               (fxsigptot(jasy,i)+fxsigmtot(jasy,i))
            !    &               *fxmtttot(i)*mttw/sigma
          end if
        end do
        asym_int(jasy)=(sfxsigptot(jasy)-sfxsigmtot(jasy))/(sfxsigptot(jasy)+sfxsigmtot(jasy))
        write(10,*)'END'
        !         write(10,*)'(total asymmetry:',asym_int(jasy),')'
        !         write(10,*)'(integrated asymmetry:',snorm(jasy),' )'
      end if
    end do
  end if
  write(10,*) 'CLOSE'
  close(10)
  print *, "...complete."
end subroutine finalise_distributions

end module distributions
