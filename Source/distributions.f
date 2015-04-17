module distributions

  use configuration, only: print_all_distributions, print_2d_distributions, include_transverse, include_asymmetries, o_asym, &
                           n_fb_asymmetries, n_asymmetries, final_state, initial_state, n_final, tops_decay, include_errors
  use mathematics, only: pi
  use kinematics, only: sigma
  use integration, only: it
  use class_Histogram

  implicit none

  public

  integer :: o_etmiss
  integer :: o_pt356
  integer :: o_eta356
  integer :: o_phi356
  integer :: o_pt478
  integer :: o_eta478
  integer :: o_phi478
  integer :: o_mtt
  integer :: o_mtt_reco
  integer :: o_m478
  integer :: o_m356_reco
  integer :: o_beta
  integer :: o_cost
  integer :: o_et
  integer :: o_delta_y
  integer :: o_fl
  integer :: o_cosfl
  integer :: o_dphi
  integer :: o_cost5
  integer :: o_cost7
  integer :: o_ct7ct5

  type(histogram) :: histetmiss
  type(histogram) :: histpt356
  type(histogram) :: histeta356
  type(histogram) :: histphi356
  type(histogram) :: histpt478
  type(histogram) :: histeta478
  type(histogram) :: histphi478
  type(histogram) :: histmtt
  type(histogram) :: histmtt_reco
  type(histogram) :: histm478
  type(histogram) :: histm356_reco
  type(histogram) :: histbeta
  type(histogram) :: histcost
  type(histogram) :: histet
  type(histogram) :: histdelta_y
  type(histogram) :: histfl
  type(histogram) :: histcosfl
  type(histogram) :: histdphi
  type(histogram) :: histcost5
  type(histogram) :: histcost7
  type(histogram) :: histct7ct5

  ! distributions in pts of asymptotic particles
  real :: ptmax(8),ptmin(8),ptw(8)
  real :: xpt(8,500),fxpt(8,500,20),fxpttot(8,500)
  real :: sumw2pt(8,500,20),sumw2pttot(8,500)
  integer :: o_pt(8)
  integer :: ndiv_pt(8)

  ! distributions in etas of asymptotic particles
  real :: etamax(8),etamin(8),etaw(8)
  real :: xeta(8,500),fxeta(8,500,20),fxetatot(8,500)
  real :: sumw2eta(8,500,20),sumw2etatot(8,500)
  integer :: o_eta(8)
  integer :: ndiv_eta(8)

  ! distributions in phis of asymptotic particles
  real :: phimax(8),phimin(8),phiw(8)
  real :: xphi(8,500),fxphi(8,500,20),fxphitot(8,500)
  real :: sumw2phi(8,500,20),sumw2phitot(8,500)
  integer :: o_phi(8)
  integer :: ndiv_phi(8)

  ! distributions in ycol of asymptotic particles
  real :: ycolmax(8),ycolmin(8),ycolw(8)
  real :: xycol(8,500),fxycol(8,500,20),fxycoltot(8,500)
  real :: sumw2ycol(8,500,20),sumw2ycoltot(8,500)
  integer :: o_ycol(8)
  integer :: ndiv_ycol(8)

  ! distributions in transverse variables
  integer :: ntrans
  parameter (ntrans=10)
  real :: transmax(ntrans),transmin(ntrans),transw(ntrans)
  real :: xtrans(ntrans,500),fxtrans(ntrans,500,20) ,fxtranstot(ntrans,500)
  real :: sumw2trans(ntrans,500,20) ,sumw2transtot(ntrans,500)
  integer :: o_tran(ntrans)
  integer :: ndiv_trans(ntrans)
  real :: sfxtranstot(ntrans),sfxtransdptot(ntrans)

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

  ! trans 2d dist
  real :: xtransdp(ntrans,500,500), fxtransdp(ntrans,500,500,20),fxtransdptot(ntrans,500,500)
  integer :: include_transversedp(ntrans)

  public :: initialise_distributions
  public :: generate_bins
  public :: print_distributions
  public :: check_distributions

  
!   real :: cnorm(20)
  real :: atot(n_asymmetries),atoterr(n_asymmetries)
  integer, private :: i, j, k, ip, iasy, jasy, itrans
  real, private :: sfxpttot(8), sfxetatot(8), sfxphitot(8), sfxycoltot(8)
  real, private :: sfxsigptot(n_asymmetries), sfxsigmtot(n_asymmetries), asym_int(n_asymmetries)

contains

  subroutine initialise_distributions

    implicit none 

    print*, "Initialising distributions"
  
    ! sets flags, binning range and divisions.

    do i=1,n_final
      o_pt(i)=print_all_distributions
      ptmax(i)=7000.d0/(1+tops_decay*6)
      ptmin(i)=0.d0
      ndiv_pt(i)=70
      ! eta distributions
      o_eta(i)=print_all_distributions
      etamax(i)=+10
      etamin(i)=-10
      ndiv_eta(i)=50
      ! phi distributions
      o_phi(i)=print_all_distributions
      phimax(i)=+pi
      phimin(i)=-pi
      ndiv_phi(i)=50
      ! ycol distributions
      o_ycol(i)=print_all_distributions
      ycolmax(i)=+4.d0
      ycolmin(i)=-4.d0
      ndiv_ycol(i)=100
    end do
    o_etmiss = print_all_distributions
    o_pt356 = print_all_distributions
    o_eta356 = print_all_distributions
    o_phi356 = print_all_distributions
    o_pt478 = print_all_distributions
    o_eta478 = print_all_distributions
    o_phi478 = print_all_distributions
    o_m478 = print_all_distributions
    o_m356_reco = print_all_distributions
    o_mtt = 1
    o_mtt_reco  =  print_all_distributions
    o_beta = print_all_distributions
    o_cost = print_all_distributions
    o_et = print_all_distributions
    o_delta_y = print_all_distributions
    o_fl = include_asymmetries
    o_cosfl = include_asymmetries
    o_dphi = include_asymmetries
    o_cost5 = include_asymmetries
    o_cost7 = include_asymmetries
    o_ct7ct5 = include_asymmetries    

    histetmiss = histogram("ETmiss","d#sigma-/dE_{tt}--[pb/GeV]", 'M_{tt}--[GeV]', 0.d0, 1000.d0, 140)
    histmtt = histogram("Mtt","d#sigma-/dM_{tt}--[pb/GeV]", 'M_{tt}--[GeV]', 0.d0, 14000.d0, 140)
    histmtt_reco = histogram("Mtt_reco","d#sigma-/dM_{tt}--[pb/GeV]", 'M_{tt}--[GeV]', 0.d0, 14000.d0, 140)
    histbeta = histogram("beta","d#sigma-/d#beta_{tt}--[pb/GeV]", 'beta_{t}', 0.d0, 1000.d0, 100)
    histpt356 = histogram("pt356","d#sigma-/p^{t}_{T}--[pb]", "p^{t}_{T}", 0.d0, 1000.d0, 100)
    histeta356 = histogram("eta356","d#sigma-/-#eta_{t}-[pb]", "#eta_{t}", 0.d0, 1000.d0, 100)
    histphi356 = histogram("phi356","d#sigma-/-#phi_{t}-[pb]", "#phi_{t}", -10.d0, 10.d0, 50)
    histpt478 = histogram("pt478","d#sigma-/p^{#bar{t}}_{T}--[pb]", "p^{#bar{t}}_{T}", 0.d0, 1000.d0, 100)
    histeta478 = histogram("eta478","d#sigma-/-#eta_{#bar{t}}-[pb]", "#eta_{#bar{t}}", -10, 10, 100)
    histphi478 = histogram("phi478","d#sigma-/#phi_{#bar{t}}--[pb]", "#phi_{#bar{t}}", -pi, pi, 100)
    histm478 = histogram("m478","d#sigma-/--[pb]", "", 0.d0, 1000.d0, 100)
    histm356_reco = histogram("m356_reco","d#sigma-/--[pb]", "", 0.d0, 1000.d0, 100)
    histbeta = histogram("beta","d#sigma-/--[pb]", "", 0.d0, 1000.d0, 100)
    histcost = histogram("cost","d#sigma-/--[pb]", "", -1.d0, 1000.d0, 100)
    histet = histogram("Et","d#sigma-/--[pb]", "", 0.d0, 1000.d0, 100)
    histdelta_y = histogram("delta_y","d#sigma-/--[pb]", "", -4, 4, 100)
    histfl = histogram("phil","d#sigma-/--[pb]", "", 0.d0, 1000.d0, 100)
    histcosfl = histogram("cosphil","d#sigma-/--[pb]", "", -1, 1, 100)
    histdphi = histogram("dphi","d#sigma-/--[pb]", "", 0.d0, 2*pi, 100)
    histcost5 = histogram("cost5","d#sigma-/--[pb]", "", -1.d0, 1.d0, 100)
    histcost7 = histogram("cost7","d#sigma-/--[pb]", "", -1.d0, 1.d0, 100)
    histct7ct5 = histogram("ct7ct5","d#sigma-/--[pb]", "", -1.d0, 1.d0, 100)

    if (o_etmiss == 1) call histetmiss%initialise()
    if (o_pt356 == 1) call histpt356%initialise()
    if (o_eta356 == 1) call histeta356%initialise()
    if (o_phi356 == 1) call histphi356%initialise()
    if (o_pt478 == 1) call histpt478%initialise()
    if (o_eta478 == 1) call histeta478%initialise()
    if (o_phi478 == 1) call histphi478%initialise()
    if (o_m478 == 1) call histm478%initialise()
    if (o_m356_reco == 1) call histm356_reco%initialise()
    if (o_beta == 1) call histbeta%initialise()
    if (o_mtt == 1) call histmtt%initialise
    if (o_mtt_reco == 1) call histmtt_reco%initialise()
    if (o_beta == 1) call histbeta%initialise()
    if (o_cost == 1) call histcost%initialise()
    if (o_et == 1) call histet%initialise()
    if (o_delta_y == 1) call histdelta_y%initialise()
    if (o_fl == 1) call histfl%initialise()
    if (o_cosfl == 1) call histcosfl%initialise()
    if (o_dphi == 1) call histdphi%initialise()
    if (o_cost5 == 1) call histcost5%initialise()
    if (o_cost7 == 1) call histcost7%initialise()
    if (o_ct7ct5 == 1) call histct7ct5%initialise()

    ! transverse variables
    do itrans=1,ntrans
      if(final_state == 0 )then
        o_tran(itrans)=0
      else
        o_tran(itrans)=include_transverse
      end if
    end do
    ! invarient mass of the visible decay products of the tt pair
    transmax(1)=4000
    transmin(1)=0.d0
    ndiv_trans(1)=40
    ! sum of tranvserse energy
    transmax(2)=4000
    transmin(2)=0.d0
    ndiv_trans(2)=40
    ! transverse mass 1
    transmax(3)=4000
    transmin(3)=0.d0
    ndiv_trans(3)=40
    ! transverse mass 2
    transmax(4)=4000
    transmin(4)=0.d0
    ndiv_trans(4)=40
    ! transverse mass 3
    transmax(5)=4000
    transmin(5)=0.d0
    ndiv_trans(5)=40
    !  lepton transverse mass
    transmax(6)=500
    transmin(6)=0.d0
    ndiv_trans(6)=40
    ! contransverse mass 1
    transmax(7)=4000
    transmin(7)=0.d0
    ndiv_trans(7)=40
    ! contransverse mass 2
    transmax(8)=4000
    transmin(8)=0.d0
    ndiv_trans(8)=40
    ! contransverse mass 3
    transmax(9)=4000
    transmin(9)=0.d0
    ndiv_trans(9)=40
    ! lepton contransverse mass
    transmax(10)=500
    transmin(10)=0.d0
    ndiv_trans(10)=50

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
    ! dphi2d
    if((o_dphi == 1) .and. (o_mtt == 1))then
      o_dphi2d=print_2d_distributions
    else
      o_dphi2d=0
    end if
    ! dtransph
    do itrans=1, ntrans
      if((o_dphi == 1) .and. (o_tran(itrans) == 1))then
        include_transversedp(itrans)=print_2d_distributions
      else
        include_transversedp(itrans)=0
      end if
    end do
    ! asymmetries
    do iasy=1,n_asymmetries
      o_asym(iasy)=include_asymmetries
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

    ! disable A_PV
    o_asym(3) = 0

    if (final_state == 0) then
      ! disable 2to6 variables 
      do ip = 5, 8
        o_pt(i) = 0
        o_eta(i) = 0
        o_phi(i) = 0
      end do
      do itrans = 1, ntrans
        o_tran(itrans) = 0
      end do
      o_tran = 0
      o_pt356 = 0
      o_eta356 = 0
      o_phi356 = 0
      o_pt478  = 0
      o_eta478 = 0
      o_phi478 = 0
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
      o_m356_reco = 0
      o_m478 = 0
    end if

    if (final_state > 0) then
      ! disable non 2to6 variables
      o_asym(1) = 0
      o_asym(2) = 0
      o_asym(3) = 0
    end if

    if (final_state == 1) then 
      ! disable non-useful variables in dileptonic
      do i = 4, 10
        o_asym(i) = 0
      end do
      o_m356_reco = 0
      o_m478 = 0
    end if

    if (final_state == 2) then
      ! disable non-useful variables in semi-hadronic
      do itrans = 1, ntrans
        o_tran(itrans) = 0
        include_transversedp(itrans) = 0
      end do
      o_cost7 = 0
      o_ct7ct5 = 0
      o_dphi = 0
      o_dphi2d = 0
      o_ct7ct52d = 0
      o_cost72d = 0
      o_cost52d = 0
      o_etmiss = 0
      o_asym(11) = 0
    end if

    if (final_state == 3) then
      ! disable non-useful variables in fully hadronic
      do itrans = 1, ntrans
        o_tran(itrans) = 0
        include_transversedp(itrans) = 0
      end do
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
      o_asym(6) = 0
      o_asym(10) = 0
      o_asym(11) = 0
      o_asym(12) = 0
      o_m356_reco = 0
      o_m478 = 0
    end if

  end subroutine initialise_distributions

  subroutine generate_bins
    ! (finds bin width, finds midpoints.)
    implicit none

    print*, "Generating bins"

    do ip=3,n_final
      if(o_pt(ip) == 1)then
        ptw(ip)=(ptmax(ip)-ptmin(ip))/ndiv_pt(ip)
        do j=1,ndiv_pt(ip)
          xpt(ip,j)=ptmin(ip)+ptw(ip)*(j-1)+ptw(ip)/2.d0
        end do
      end if
      if(o_eta(ip) == 1)then
        etaw(ip)=(etamax(ip)-etamin(ip))/ndiv_eta(ip)
        do j=1,ndiv_eta(ip)
          xeta(ip,j)=etamin(ip)+etaw(ip)*(j-1)+etaw(ip)/2.d0
        end do
      end if
      if(o_phi(ip) == 1)then
        phiw(ip)=(phimax(ip)-phimin(ip))/ndiv_phi(ip)
        do j=1,ndiv_phi(ip)
          xphi(ip,j)=phimin(ip)+phiw(ip)*(j-1)+phiw(ip)/2.d0
        end do
      end if
      if(o_ycol(ip) == 1)then
        ycolw(ip)=(ycolmax(ip)-ycolmin(ip))/ndiv_ycol(ip)
        do j=1,ndiv_ycol(ip)
          xycol(ip,j)=ycolmin(ip)+ycolw(ip)*(j-1)+ycolw(ip)/2.d0
        end do
      end if
    end do

    do itrans=1,ntrans
      if(o_tran(itrans) == 1)then
        transw(itrans)=(transmax(itrans)-transmin(itrans)) &
        /ndiv_trans(itrans)
        do i=1,ndiv_trans(itrans)
          xtrans(itrans,i)=transmin(itrans)+transw(itrans)*(i-1) &
          +transw(itrans)/2.d0
        end do
      end if
    end do

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

  subroutine print_distributions

    implicit none

    integer :: ndiv_sig

    print*, "Printing histograms"

    write(10,*) '-----------------------------------------------------'
    write(10,*)'HISTOGRAMS'
    do ip=3,8
      ! plot distributions in pt
      if(o_pt(ip) == 1)then
        sfxpttot(ip)=0d0
        do j=1,ndiv_pt(ip)
          fxpttot(ip,j)=0.d0
          do i=1,it
            fxpt(ip,j,i)=fxpt(ip,j,i)*sigma/cnorm(i)/ptw(ip)
            fxpttot(ip,j)=fxpttot(ip,j)+fxpt(ip,j,i)
            if (include_errors == 1) then 
            sumw2pt(ip,j,i)=sumw2pt(ip,j,i)*sigma/cnorm(i)/ptw(ip)*sigma/cnorm(i)/ptw(ip)
            sumw2pttot(ip,j)=sumw2pttot(ip,j)+sumw2pt(ip,j,i)
            else
              sumw2pttot(ip,j)=0
            end if
          end do
          sfxpttot(ip)=sfxpttot(ip)+fxpttot(ip,j)*ptw(ip)
        end do
        write(10,*)'DISTRIBUTION'
        write(10,'(a,i1)')'pt',ip
        write(10,'(a,i1,a)')'d#sigma-/dp_{t}(',ip,')--[pb/GeV]'
        write(10,'(a,i1,a)')'p_{t}(',ip,')--[GeV]'
        do i=1,ndiv_pt(ip)
          write(10,*)xpt(ip,i),fxpttot(ip,i)
        end do
        write(10,*)'END'
      end if
      ! plot distributions in eta
      if(o_eta(ip) == 1)then
        sfxetatot(ip)=0d0
        do j=1,ndiv_eta(ip)
          fxetatot(ip,j)=0.d0
          do i=1,it
            fxeta(ip,j,i)=fxeta(ip,j,i)*sigma/cnorm(i)/etaw(ip)
            fxetatot(ip,j)=fxetatot(ip,j)+fxeta(ip,j,i)
            if (include_errors == 1) then 
            sumw2eta(ip,j,i)=sumw2eta(ip,j,i)*sigma/cnorm(i)/etaw(ip)*sigma/cnorm(i)/etaw(ip)
            sumw2etatot(ip,j)=sumw2etatot(ip,j)+sumw2eta(ip,j,i)
            else
              sumw2etatot(ip,j)=0
            end if
          end do
          sfxetatot(ip)=sfxetatot(ip)+fxetatot(ip,j)*etaw(ip)
        end do
        write(10,*)'DISTRIBUTION'
        write(10,'(a,i1)')'eta',ip
        write(10,'(a,i1,a)')'d#sigma-/d#eta(',ip,')--[pb]'
        write(10,'(a,i1,a)')'#eta(',ip,')'
        do i=1,ndiv_eta(ip)
          write(10,*)xeta(ip,i),fxetatot(ip,i)
        end do
        write(10,*)'END'
      end if
      ! plot distributions in phi
      if(o_phi(ip) == 1)then
        sfxphitot(ip)=0d0
        do j=1,ndiv_phi(ip)
          fxphitot(ip,j)=0.d0
          do i=1,it
            fxphi(ip,j,i)=fxphi(ip,j,i)*sigma/cnorm(i)/phiw(ip)
            fxphitot(ip,j)=fxphitot(ip,j)+fxphi(ip,j,i)
            if (include_errors == 1) then 
            sumw2phi(ip,j,i)=sumw2phi(ip,j,i)*sigma/cnorm(i)/phiw(ip)*sigma/cnorm(i)/phiw(ip)
            sumw2phitot(ip,j)=sumw2phitot(ip,j)+sumw2phi(ip,j,i)
            else
              sumw2phitot(ip,j)=0
            end if
          end do
          sfxphitot(ip)=sfxphitot(ip)+fxphitot(ip,j)*phiw(ip)
        end do
        write(10,*)'DISTRIBUTION'
        write(10,'(a,i1)')'phi',ip
        write(10,'(a,i1,a)')'d#sigma-/d#phi(',ip,')--[pb/rad]'
        write(10,'(a,i1,a)')'#phi(',ip,')--[rad]'
        do i=1,ndiv_phi(ip)
          write(10,*)xphi(ip,i),fxphitot(ip,i)
        end do
        write(10,*)'END'
      end if
      ! plot distributions in ycol
      if(o_ycol(ip) == 1)then
        sfxycoltot(ip)=0d0
        do j=1,ndiv_ycol(ip)
          fxycoltot(ip,j)=0.d0
          do i=1,it
            fxycol(ip,j,i)=fxycol(ip,j,i)*sigma/cnorm(i)/ycolw(ip)
            fxycoltot(ip,j)=fxycoltot(ip,j)+fxycol(ip,j,i)
            if (include_errors == 1) then 
            sumw2ycol(ip,j,i)=sumw2ycol(ip,j,i)*sigma/cnorm(i)/ycolw(ip)*sigma/cnorm(i)/ycolw(ip)
            sumw2ycoltot(ip,j)=sumw2ycoltot(ip,j)+sumw2ycol(ip,j,i)
            else
              sumw2ycoltot(ip,j)=0
            end if
          end do
          sfxycoltot(ip)=sfxycoltot(ip)+fxycoltot(ip,j)*ycolw(ip)
        end do
        write(10,*)'DISTRIBUTION'
        write(10,'(a,i1)')'y',ip
        write(10,'(a,i1,a)')'d#sigma-/dy(',ip,')--[pb]'
        write(10,'(a,i1,a)')'y(',ip,')'
        do i=1,ndiv_ycol(ip)
          write(10,*)xycol(ip,i),fxycoltot(ip,i)
        end do
        write(10,*)'END'
      end if
    end do
       
    if(o_etmiss == 1) call histetmiss % finalise()
    if(o_pt356 == 1) call histpt356 % finalise()
    if(o_eta356 == 1) call histeta356 % finalise()
    if(o_phi356 == 1) call histphi356 % finalise()
    if(o_pt478 == 1) call histpt478 % finalise()
    if(o_eta478 == 1) call histeta478 % finalise()
    if(o_phi478 == 1) call histphi478 % finalise()
    if(o_mtt == 1) call histmtt % finalise()
    if(o_mtt_reco == 1) call histmtt_reco % finalise()
    if(o_m478 == 1) call histm478 % finalise()
    if(o_m356_reco == 1) call histm356_reco % finalise()
    if(o_beta == 1) call histbeta % finalise()
    if(o_cost == 1) call histcost % finalise()
    if(o_et == 1) call histet % finalise()
    if(o_delta_y == 1) call histdelta_y % finalise()
    if(o_fl == 1) call histfl%finalise()
    if(o_cosfl == 1) call histcosfl%finalise()
    if(o_dphi == 1) call histdphi%finalise()
    if(o_cost5 == 1) call histcost5%finalise()
    if(o_cost7 == 1) call histcost7%finalise()
    if(o_ct7ct5 == 1) call histct7ct5%finalise()

!     ! plot distributions in all transverse variables
!     do itrans=1,ntrans
!       if(o_tran(itrans) == 1)then
!         sfxtranstot(itrans)=0d0
!         do j=1,ndiv_trans(itrans)
!           fxtranstot(itrans,j)=0.d0
!           do i=1,it
!             fxtrans(itrans,j,i)=fxtrans(itrans,j,i) &
!             *sigma/cnorm(i)/transw(itrans)
!             fxtranstot(itrans,j)=fxtranstot(itrans,j) &
!             +fxtrans(itrans,j,i)
!             if (include_errors == 1) then 
!               sumw2trans(itrans,j,i)=sumw2trans(itrans,j,i) &
!               *sigma/cnorm(i)/transw(itrans)*sigma/cnorm(i)/transw(itrans)
!               sumw2transtot(itrans,j)=sumw2transtot(itrans,j) &
!               +sumw2trans(itrans,j,i)
!             else
!               sumw2transtot(itrans,j)=0
!             end if
!           end do
!           sfxtranstot(itrans)=sfxtranstot(itrans)+ &
!           fxtranstot(itrans,j)*transw(itrans)
!         end do
!         write(10,*)'DISTRIBUTION'
!         if (itrans == 1)then
!           write(10,*)'Mvis'
!           write(10,*)'d#sigma-/dM_{vis}--[pb/GeV]'
!           write(10,*)'M_{vis}--[GeV]'
!         else if (itrans == 2)then
!           write(10,*)'HT'
!           write(10,*)'d#sigma-/dH_{T}--[pb/GeV]'
!           write(10,*)'H_{T}--[GeV]'
!         else if (itrans == 3)then
!           write(10,*)'M_T1'
!           write(10,*)'d#sigma-/dM_{T1}--[pb/GeV]'
!           write(10,*)'M_{T1}--[GeV]'
!         else if (itrans == 4)then
!           write(10,*)'M_T2'
!           write(10,*)'d#sigma-/dM_{T2}--[pb/GeV]'
!           write(10,*)'M_{T2}--[GeV]'
!         else if (itrans == 5)then
!           write(10,*)'M_T3'
!           write(10,*)'d#sigma-/dM_{T3}--[pb/GeV]'
!           write(10,*)'M_{T3}--[GeV]'
!         else if (itrans == 6)then
!           write(10,*)'MlT'
!           write(10,*)'d#sigma-/dM^{l}_{T}--[pb/GeV]'
!           write(10,*)'M^{l}_{T}--[GeV]'
!         else if (itrans == 7)then
!           write(10,*)'M_CT1'
!           write(10,*)'d#sigma-/dM_{CT1}--[pb/GeV]'
!           write(10,*)'M_{CT1}--[GeV]'
!         else if (itrans == 8)then
!           write(10,*)'M_CT2'
!           write(10,*)'d#sigma-/dM_{CCT2}--[pb/GeV]'
!           write(10,*)'M_{CT2}--[GeV]'
!         else if (itrans == 9)then
!           write(10,*)'M_CT3'
!           write(10,*)'d#sigma-/dM_{CT3}--[pb/GeV]'
!           write(10,*)'M_{CT3}--[GeV]'
!         else if (itrans == 10)then
!           write(10,*)'MlCT'
!           write(10,*)'d#sigma-/dM^{l}_{CT}--[pb/GeV]'
!           write(10,*)'m^{l}_{CT}--[GeV]'
!         else
!           continue
!         end if
!         do i=1,ndiv_trans(itrans)
!           write(10,*)xtrans(itrans,i),fxtranstot(itrans,i)
!         end do
!         write(10,*)'END'
!       end if
!     end do

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

  end subroutine print_distributions

  subroutine  check_distributions
    real :: diff_max = 1e-12
    integer :: n_error = 0

    print*, "Checking histograms"

    do ip=3,n_final
      if(o_pt(ip) == 1)then
        if(abs(sigma-sfxpttot(ip))>diff_max)then
          write(10,*)'pt',ip,' error:',sfxpttot(ip)
          n_error=n_error+1
        end if
      end if
      if(o_eta(ip) == 1)then
        if(abs(sigma-sfxetatot(ip))>diff_max)then
          write(10,*)'eta',ip,' error:',sfxetatot(ip)
          n_error=n_error+1
        end if
      end if
      if(o_phi(ip) == 1)then
        if(abs(sigma-sfxphitot(ip))>diff_max)then
          write(10,*)'phi',ip,' error:',sfxphitot(ip)
          n_error=n_error+1
        end if
      end if
    end do

    do iasy=1,n_asymmetries
      if(o_asym(iasy) == 0)then
        continue
      else
        if(abs(atot(iasy)-asym_int(iasy))>diff_max)then
          if(iasy /= 7)then ! fuck AttbRFB
            write(10,*)'a error:',iasy,asym_int(iasy)
            n_error=n_error+1
          end if
        end if
      end if
    end do
    if(n_error > 0)write(10,*)'integration errors:',n_error
    write(10,*)'CLOSE'
    close(10)

  end subroutine check_distributions

end module distributions
