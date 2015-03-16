module distributions

  use configuration, only: o_distros, o_dist2d, o_trans, o_asyms, o_asym, nspat, nasym, ifinal_state, initial_state, nfinal
  use kinematics, only: pi, sigma
  use integration, only: it

  ! distributions in pts of asymptotic particles
  real :: ptmax(8),ptmin(8),ptw(8)
  real :: xpt(8,500),fxpt(8,500,20),fxpttot(8,500)
  integer :: o_pt(8)
  integer :: ndiv_pt(8)

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

  ! distribution in etmiss
  real :: etmissmax,etmissmin,etmissw
  real :: xetmiss(500),fxetmiss(500,20),fxetmisstot(500)
  integer :: o_etmiss
  integer :: ndiv_etmiss

  ! distribution in pt of the top
  real :: pt356max,pt356min,pt356w
  real :: xpt356(500),fxpt356(500,20),fxpt356tot(500)
  integer :: o_pt356
  integer :: ndiv_pt356

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

  ! distribution in pt of the anti-top
  real :: pt478max,pt478min,pt478w
  real :: xpt478(500),fxpt478(500,20),fxpt478tot(500)
  integer :: o_pt478
  integer :: ndiv_pt478

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
  real :: rmttmax,rmttmin,rmttw
  real :: xrmtt(500),fxrmtt(500,20),fxrmtttot(500)
  integer :: o_rmtt
  integer :: ndiv_rmtt

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
  real :: etmax,etmin,etw
  real :: xet(500),fxet(500,20),fxettot(500)
  integer :: o_et
  integer :: ndiv_et

  ! distribution in delta_y
  real :: delta_ymax,delta_ymin,delta_yw
  real :: xdelta_y(500),fxdelta_y(500,20), fxdelta_ytot(500)
  integer :: o_delta_y
  integer :: ndiv_delta_y

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

  public :: initialise_distributions
  public :: generate_bins
  public :: print_distributions
  public :: check_distributions

  real :: sfxpttot(8), sfxetatot(8), sfxphitot(8), sfxycoltot(8)
  real :: sfxsigptot(nasym), sfxsigmtot(nasym), asym_int(nasym)
  real :: cnorm(20)
  real :: atot(nasym),atoterr(nasym)


  contains

    subroutine initialise_distributions
    
      ! sets flags, binning range and divisions.

      do i=1,nfinal
        o_pt(i)=o_distros
        ptmax(i)=7000.d0/(1+initial_state*6)
        ptmin(i)=0.d0
        ndiv_pt(i)=70
        ! eta distributions
        o_eta(i)=o_distros
        etamax(i)=+10
        etamin(i)=-10
        ndiv_eta(i)=50
        ! phi distributions
        o_phi(i)=o_distros
        phimax(i)=+pi
        phimin(i)=-pi
        ndiv_phi(i)=50
        ! ycol distributions
        o_ycol(i)=o_distros
        ycolmax(i)=+4.d0
        ycolmin(i)=-4.d0
        ndiv_ycol(i)=100
      end do
      ! missing transverse momentum
      o_etmiss=o_distros
      etmissmax=7000.d0/(1+initial_state*6)
      etmissmin=0.d0
      ndiv_etmiss=70
      ! top transverse momentum
      o_pt356=o_distros
      pt356max=7000.d0/(1+initial_state*6)
      pt356min=0.d0
      ndiv_pt356=70
      ! 2to6 top pseudorapidity
      o_eta356=o_distros
      eta356max=+10
      eta356min=-10
      ndiv_eta356=50
      ! 2to6 top pseudorapidity
      o_phi356=o_distros
      phi356max=+pi
      phi356min=-pi
      ndiv_phi356=50
      ! anti-top transverse momentum
      o_pt478=o_distros
      pt478max=7000.d0/(1+initial_state*6)
      pt478min=0.d0
      ndiv_pt478=70
      ! 2to6 anti-top pseudorapidity
      o_eta478=o_distros
      eta478max=+10
      eta478min=-10
      ndiv_eta478=50
      ! 2to6 top pseudorapidity
      o_phi478=o_distros
      phi478max=+pi
      phi478min=-pi
      ndiv_phi478=50
      ! invarient mass of tt pair (always on)
      o_rmtt=1
      rmttmax=14000.d0/(1+initial_state*6)
      rmttmin=0.d0
      ndiv_rmtt=140
      ! boost of parton com
      o_beta=o_distros
      betamax=1000.d0
      betamin=0.d0
      ndiv_beta=100
      ! costheta
      o_cost=o_distros
      costmax=+1.d0
      costmin=-1.d0
      ndiv_cost=50
      ! top energy
      o_et=o_distros
      etmax=7000.d0/(1+initial_state*6)
      etmin=0.d0
      ndiv_et=70
      ! delta_y
      o_delta_y=o_distros
      delta_ymax=4.d0
      delta_ymin=-4.d0
      ndiv_delta_y=100
      ! transverse variables
      do itrans=1,ntrans
        if(ifinal_state == 0)then
          o_tran(itrans)=0
        else
          o_tran(itrans)=o_trans
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
      ! phi_l
      o_fl=o_asyms
      flmax=+2*pi
      flmin=0
      ndiv_fl=100
      ! cosphi_l
      o_cosfl=o_asyms
      cosflmax=+1.d0
      cosflmin=-1.d0
      ndiv_cosfl=100
      ! delta phi
      o_dphi=o_asyms
      dphimax=+pi
      dphimin=0
      ndiv_dphi=10
      ! cost5
      o_cost5=o_asyms
      cost5max=+1
      cost5min=-1
      ndiv_cost5=10
      ! cost7
      o_cost7=o_asyms
      cost7max=+1
      cost7min=-1
      ndiv_cost7=10
      !  ct7ct5
      o_ct7ct5=o_asyms
      ct7ct5max=+1
      ct7ct5min=-1
      ndiv_ct7ct5=10
      ! sigp
      o_sigp=o_asyms
      sigpmax=rmttmax
      sigpmin=rmttmin
      ndiv_sigp=ndiv_rmtt/5
      ! sigm
      o_sigm=o_asyms
      sigmmax=rmttmax
      sigmmin=rmttmin
      ndiv_sigm=ndiv_rmtt/5
      ! dphi2d
      if((o_dphi == 1) .and. (o_rmtt == 1))then
        o_dphi2d=o_dist2d
      else
        o_dphi2d=0
      end if
      ! dtransph
      do itrans=1, ntrans
        if((o_dphi == 1) .and. (o_tran(itrans) == 1))then
          o_transdp(itrans)=o_dist2d
        else
          o_transdp(itrans)=0
        end if
      end do
      ! asymmetries
      do iasy=1,nasym
        o_asym(iasy)=o_asyms
      end do

      if (ifinal_state == 0)then
        do ip=5,8
          o_pt(i)   = 0
          o_eta(i)  = 0
          o_phi(i)  = 0
        end do
        do itrans=1,ntrans
          o_tran(itrans)=0
        end do
        o_tran  = 0
        o_pt356  = 0
        o_eta356 = 0
        o_phi356 = 0
        o_pt478  = 0
        o_eta478 = 0
        o_phi478 = 0
        o_etmiss = 0
        o_ht     = 0
        o_fl     = 0
        o_dphi   = 0
        o_cosfl  = 0
        o_cost7  = 0
        o_cost5  = 0
        o_ct7ct5 = 0
        o_dphi2d = 0
        o_asym(9) = 0    ! turn off a_l
      else if (ifinal_state > 0)then
        o_asym(1) = 0   ! turn off a_ll
        o_asym(2) = 0   ! turn off a_l
        o_asym(3) = 0   ! turn off a_pv
      end if

    end subroutine initialise_distributions

    subroutine generate_bins
      ! (finds bin width, finds midpoints.)

      do ip=3,nfinal
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

      if(o_etmiss == 1)then
        etmissw=(etmissmax-etmissmin)/ndiv_etmiss
        do i=1,ndiv_etmiss
          xetmiss(i)=etmissmin+etmissw*(i-1)+etmissw/2.d0
        end do
      end if

      if(o_pt356 == 1)then
        pt356w=(pt356max-pt356min)/ndiv_pt356
        do i=1,ndiv_pt356
          xpt356(i)=pt356min+pt356w*(i-1)+pt356w/2.d0
        end do
      end if

      if(o_eta356 == 1)then
        eta356w=(eta356max-eta356min)/ndiv_eta356
        do i=1,ndiv_eta356
          xeta356(i)=eta356min+eta356w*(i-1)+eta356w/2.d0
        end do
      end if

      if(o_phi356 == 1)then
        phi356w=(phi356max-phi356min)/ndiv_phi356
        do i=1,ndiv_phi356
          xphi356(i)=phi356min+phi356w*(i-1)+phi356w/2.d0
        end do
      end if

      if(o_pt478 == 1)then
        pt478w=(pt478max-pt478min)/ndiv_pt478
        do i=1,ndiv_pt478
          xpt478(i)=pt478min+pt478w*(i-1)+pt478w/2.d0
        end do
      end if

      if(o_eta478 == 1)then
        eta478w=(eta478max-eta478min)/ndiv_eta478
        do i=1,ndiv_eta478
          xeta478(i)=eta478min+eta478w*(i-1)+eta478w/2.d0
        end do
      end if

      if(o_phi478 == 1)then
        phi478w=(phi478max-phi478min)/ndiv_phi478
        do i=1,ndiv_phi478
          xphi478(i)=phi478min+phi478w*(i-1)+phi478w/2.d0
        end do
      end if

      if(o_rmtt == 1)then
        rmttw=(rmttmax-rmttmin)/ndiv_rmtt
        do i=1,ndiv_rmtt
          xrmtt(i)=rmttmin+rmttw*(i-1)+rmttw/2.d0
        end do
      end if

      if(o_beta == 1)then
        betaw=(betamax-betamin)/ndiv_beta
        do i=1,ndiv_beta
          xbeta(i)=betamin+betaw*(i-1)+betaw/2.d0
        end do
      end if

      if(o_cost == 1)then
        costw=(costmax-costmin)/ndiv_cost
        do i=1,ndiv_cost
          xcost(i)=costmin+costw*(i-1)+costw/2.d0
        end do
      end if

      if(o_et == 1)then
        etw=(etmax-etmin)/ndiv_et
        do i=1,ndiv_et
          xet(i)=etmin+etw*(i-1)+etw/2.d0
        end do
      end if

      if(o_delta_y == 1)then
        delta_yw=(delta_ymax-delta_ymin)/ndiv_delta_y
        do i=1,ndiv_delta_y
          xdelta_y(i)=delta_ymin+delta_yw*(i-1)+delta_yw/2.d0
        end do
      end if

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

      if(o_fl == 1)then
        flw=(flmax-flmin)/ndiv_fl
        do i=1,ndiv_fl
          xfl(i)=flmin+flw*(i-1)+flw/2.d0
        end do
      end if

      if(o_cosfl == 1)then
        cosflw=(cosflmax-cosflmin)/ndiv_cosfl
        do i=1,ndiv_cosfl
          xcosfl(i)=cosflmin+cosflw*(i-1)+cosflw/2.d0
        end do
      end if

      if(o_dphi == 1)then
        dphiw=(dphimax-dphimin)/ndiv_dphi
        do i=1,ndiv_dphi
          xdphi(i)=dphimin+dphiw*(i-1)+dphiw/2.d0
        end do
      end if

      if(o_cost5 == 1)then
        cost5w=(cost5max-cost5min)/ndiv_cost5
        do i=1,ndiv_cost5
          xcost5(i)=cost5min+cost5w*(i-1)+cost5w/2.d0
        end do
      end if

      if(o_cost7 == 1)then
        cost7w=(cost7max-cost7min)/ndiv_cost7
        do i=1,ndiv_cost7
          xcost7(i)=cost7min+cost7w*(i-1)+cost7w/2.d0
        end do
      end if

      if(o_ct7ct5 == 1)then
        ct7ct5w=(ct7ct5max-ct7ct5min)/ndiv_ct7ct5
        do i=1,ndiv_ct7ct5
          xct7ct5(i)=ct7ct5min+ct7ct5w*(i-1)+ct7ct5w/2.d0
        end do
      end if

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

      integer :: ip

      write(*,*)'sigma',sigma

      write(*,*)''
      write(*,*)'HISTOGRAMS'
      do ip=3,8
        ! plot distributions in pt
        if(o_pt(ip) == 1)then
          sfxpttot(ip)=0d0
          do j=1,ndiv_pt(ip)
            fxpttot(ip,j)=0.d0
            do i=1,it
              fxpt(ip,j,i)=fxpt(ip,j,i)*sigma/cnorm(i)/ptw(ip)
              fxpttot(ip,j)=fxpttot(ip,j)+fxpt(ip,j,i)
            end do
            sfxpttot(ip)=sfxpttot(ip)+fxpttot(ip,j)*ptw(ip)
          end do
          write(*,*)'DISTRIBUTION'
          write(*,'(a,i1)')'pt',ip
          write(*,'(a,i1,a)')'d#sigma-/dp_{t}(',ip,')--[pb/gev]'
          write(*,'(a,i1,a)')'p_{t}(',ip,')--[gev]'
          do i=1,ndiv_pt(ip)
            write(*,*)xpt(ip,i),fxpttot(ip,i)
          end do
          write(*,*)'END'
        end if
        ! plot distributions in eta
        if(o_eta(ip) == 1)then
          sfxetatot(ip)=0d0
          do j=1,ndiv_eta(ip)
            fxetatot(ip,j)=0.d0
            do i=1,it
              fxeta(ip,j,i)=fxeta(ip,j,i)*sigma/cnorm(i)/etaw(ip)
              fxetatot(ip,j)=fxetatot(ip,j)+fxeta(ip,j,i)
            end do
            sfxetatot(ip)=sfxetatot(ip)+fxetatot(ip,j)*etaw(ip)
          end do
          write(*,*)'DISTRIBUTION'
          write(*,'(a,i1)')'eta',ip
          write(*,'(a,i1,a)')'d#sigma-/d#eta(',ip,')--[pb]'
          write(*,'(a,i1,a)')'#eta(',ip,')'
          do i=1,ndiv_eta(ip)
            write(*,*)xeta(ip,i),fxetatot(ip,i)
          end do
          write(*,*)'END'
        end if
        ! plot distributions in phi
        if(o_phi(ip) == 1)then
          sfxphitot(ip)=0d0
          do j=1,ndiv_phi(ip)
            fxphitot(ip,j)=0.d0
            do i=1,it
              fxphi(ip,j,i)=fxphi(ip,j,i)*sigma/cnorm(i)/phiw(ip)
              fxphitot(ip,j)=fxphitot(ip,j)+fxphi(ip,j,i)
            end do
            sfxphitot(ip)=sfxphitot(ip)+fxphitot(ip,j)*phiw(ip)
          end do
          write(*,*)'DISTRIBUTION'
          write(*,'(a,i1)')'phi',ip
          write(*,'(a,i1,a)')'d#sigma-/d#phi(',ip,')--[pb/rad]'
          write(*,'(a,i1,a)')'#phi(',ip,')--[rad]'
          do i=1,ndiv_phi(ip)
            write(*,*)xphi(ip,i),fxphitot(ip,i)
          end do
          write(*,*)'END'
        end if
        ! plot distributions in ycol
        if(o_ycol(ip) == 1)then
          sfxycoltot(ip)=0d0
          do j=1,ndiv_ycol(ip)
            fxycoltot(ip,j)=0.d0
            do i=1,it
              fxycol(ip,j,i)=fxycol(ip,j,i)*sigma/cnorm(i)/ycolw(ip)
              fxycoltot(ip,j)=fxycoltot(ip,j)+fxycol(ip,j,i)
            end do
            sfxycoltot(ip)=sfxycoltot(ip)+fxycoltot(ip,j)*ycolw(ip)
          end do
          write(*,*)'DISTRIBUTION'
          write(*,'(a,i1)')'y',ip
          write(*,'(a,i1,a)')'d#sigma-/dy(',ip,')--[pb]'
          write(*,'(a,i1,a)')'y(',ip,')'
          do i=1,ndiv_ycol(ip)
            write(*,*)xycol(ip,i),fxycoltot(ip,i)
          end do
          write(*,*)'END'
        end if
      end do
         
      ! plot distribution in etmiss
      if(o_etmiss == 1)then
        sfxetmisstot=0d0
        do j=1,ndiv_etmiss
          fxetmisstot(j)=0.d0
          do i=1,it
            fxetmiss(j,i)=fxetmiss(j,i)*sigma/cnorm(i)/etmissw
            fxetmisstot(j)=fxetmisstot(j)+fxetmiss(j,i)
          end do
          sfxetmisstot=sfxetmisstot+fxetmisstot(j)*etmissw
        end do
        write(*,*)'DISTRIBUTION'
        write(*,*)'etmiss'
        write(*,*)'d#sigma-/dp_{tmiss}--[pb/gev]'
        write(*,*)'p_{t}(miss)--[gev]'
        do i=1,ndiv_etmiss
          write(*,*)xetmiss(i),fxetmisstot(i)
        end do
        write(*,*)'END'
      end if

      ! plot distribution in pt356
      if(o_pt356 == 1)then
        sfxpt356tot=0d0
        do j=1,ndiv_pt356
          fxpt356tot(j)=0.d0
          do i=1,it
            fxpt356(j,i)=fxpt356(j,i)*sigma/cnorm(i)/pt356w
            fxpt356tot(j)=fxpt356tot(j)+fxpt356(j,i)
          end do
          sfxpt356tot=sfxpt356tot+fxpt356tot(j)*pt356w
        end do
        write(*,*)'DISTRIBUTION'
        write(*,*)'pt356'
        write(*,*)'d#sigma-/dp_{t}--[pb/gev]'
        write(*,*)'p_{t}(t)--[gev]'
        do i=1,ndiv_pt356
          write(*,*)xpt356(i),fxpt356tot(i)
        end do
        write(*,*)'END'
      end if
      ! plot distribution in eta356
      if(o_eta356 == 1)then
        sfxeta356tot=0d0
        do j=1,ndiv_eta356
          fxeta356tot(j)=0.d0
          do i=1,it
            fxeta356(j,i)=fxeta356(j,i)*sigma/cnorm(i)/eta356w
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
      ! plot distribution in phi356
      if(o_phi356 == 1)then
        sfxphi356tot=0d0
        do j=1,ndiv_phi356
          fxphi356tot(j)=0.d0
          do i=1,it
            fxphi356(j,i)=fxphi356(j,i)*sigma/cnorm(i)/phi356w
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
      ! plot distribution in pt478
      if(o_pt478 == 1)then
        sfxpt478tot=0d0
        do j=1,ndiv_pt478
          fxpt478tot(j)=0.d0
          do i=1,it
            fxpt478(j,i)=fxpt478(j,i)*sigma/cnorm(i)/pt478w
            fxpt478tot(j)=fxpt478tot(j)+fxpt478(j,i)
          end do
          sfxpt478tot=sfxpt478tot+fxpt478tot(j)*pt478w
        end do
        write(*,*)'DISTRIBUTION'
        write(*,*)'pt478'
        write(*,*)'d#sigma-/dp_{t}--[pb/gev]'
        write(*,*)'p_{t}(#bar{t})--[gev]'
        do i=1,ndiv_pt478
          write(*,*)xpt478(i),fxpt478tot(i)
        end do
        write(*,*)'END'
      end if
      ! plot distribution in eta478
      if(o_eta478 == 1)then
        sfxeta478tot=0d0
        do j=1,ndiv_eta478
          fxeta478tot(j)=0.d0
          do i=1,it
            fxeta478(j,i)=fxeta478(j,i)*sigma/cnorm(i)/eta478w
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
      ! plot distribution in phi478
      if(o_phi478 == 1)then
        sfxphi478tot=0d0
        do j=1,ndiv_phi478
          fxphi478tot(j)=0.d0
          do i=1,it
            fxphi478(j,i)=fxphi478(j,i)*sigma/cnorm(i)/phi478w
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
      ! plot distribution in mtt
      if(o_rmtt == 1)then
        sfxrmtttot=0d0
        do j=1,ndiv_rmtt
          fxrmtttot(j)=0.d0
          do i=1,it
            fxrmtt(j,i)=fxrmtt(j,i)*sigma/cnorm(i)/rmttw
            fxrmtttot(j)=fxrmtttot(j)+fxrmtt(j,i)
          end do
          sfxrmtttot=sfxrmtttot+fxrmtttot(j)*rmttw
        end do
        write(*,*)'DISTRIBUTION'
        write(*,*)'Mtt'
        write(*,*)'d#sigma-/dm_{tt}--[pb/gev]'
        write(*,*)'m_{tt}--[gev]'
        do i=1,ndiv_rmtt
          write(*,*)xrmtt(i),fxrmtttot(i)
        end do
        write(*,*)'END'
      end if
      ! plot distribution in beta.
      if(o_beta == 1)then
        sfxbetatot=0d0
        do j=1,ndiv_beta
          fxbetatot(j)=0.d0
          do i=1,it
            fxbeta(j,i)=fxbeta(j,i)*sigma/cnorm(i)/betaw
            fxbetatot(j)=fxbetatot(j)+fxbeta(j,i)
          end do
          sfxbetatot=sfxbetatot+fxbetatot(j)*betaw
        end do
        write(*,*)'DISTRIBUTION'
        write(*,*)'beta'
        write(*,*)'d#sigma-/d#beta_{t}--[pb]'
        write(*,*)'#beta_{t}'
        do i=1,ndiv_beta
          write(*,*)xbeta(i),fxbetatot(i)
        end do
        write(*,*)'END'
      end if
      ! plot distribution in cost
      if(o_cost == 1)then
        sfxcosttot=0d0
        do j=1,ndiv_cost
          fxcosttot(j)=0.d0
          do i=1,it
            fxcost(j,i)=fxcost(j,i)*sigma/cnorm(i)/costw
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
      ! plot distribution in et
      if(o_et == 1)then

        sfxettot=0d0
        do j=1,ndiv_et
          fxettot(j)=0.d0
          do i=1,it
            fxet(j,i)=fxet(j,i)*sigma/cnorm(i)/etw
            fxettot(j)=fxettot(j)+fxet(j,i)
          end do
          sfxettot=sfxettot+fxettot(j)*etw
        end do
        write(*,*)'DISTRIBUTION'
        write(*,*)'et'
        write(*,*)'d#sigma-/de_{t}--[pb/gev]'
        write(*,*)'e_{t}--[gev]'
        do i=1,ndiv_et
          write(*,*)xet(i),fxettot(i)
        end do
        write(*,*)'END'
      end if
      ! plot distribution in delta y
      if(o_delta_y == 1)then
        sfxdelta_ytot=0d0
        do j=1,ndiv_delta_y
          fxdelta_ytot(j)=0.d0
          do i=1,it
            fxdelta_y(j,i)=fxdelta_y(j,i)*sigma/cnorm(i)/delta_yw
            fxdelta_ytot(j)=fxdelta_ytot(j)+fxdelta_y(j,i)
          end do
          sfxdelta_ytot=sfxdelta_ytot+fxdelta_ytot(j)*delta_yw
        end do
        write(*,*)'DISTRIBUTION'
        write(*,*)'delta_y'
        write(*,*)'d#sigma-/d#delta-y--[pb]'
        write(*,*)'#delta-y'
        do i=1,ndiv_delta_y
          write(*,*)xdelta_y(i),fxdelta_ytot(i)
        end do
        write(*,*)'END'
      end if
      ! plot distributions in all transverse variables
      do itrans=1,ntrans
        if(o_tran(itrans) == 1)then
          sfxtranstot(itrans)=0d0
          do j=1,ndiv_trans(itrans)
            fxtranstot(itrans,j)=0.d0
            do i=1,it
              fxtrans(itrans,j,i)=fxtrans(itrans,j,i) &
              *sigma/cnorm(i)/transw(itrans)
              fxtranstot(itrans,j)=fxtranstot(itrans,j) &
              +fxtrans(itrans,j,i)
            end do
            sfxtranstot(itrans)=sfxtranstot(itrans)+ &
            fxtranstot(itrans,j)*transw(itrans)
          end do
          write(*,*)'DISTRIBUTION'
          if (itrans == 1)then
            write(*,*)'mvis'
            write(*,*)'d#sigma-/dm_{vis}--[pb/gev]'
            write(*,*)'m_{vis}--[gev]'
          else if (itrans == 2)then
            write(*,*)'ht'
            write(*,*)'d#sigma-/dh_{t}--[pb/gev]'
            write(*,*)'h_{t}--[gev]'
          else if (itrans == 3)then
            write(*,*)'m_t1'
            write(*,*)'d#sigma-/dm_{t1}--[pb/gev]'
            write(*,*)'m_{t}--[gev]'
          else if (itrans == 4)then
            write(*,*)'m_t2'
            write(*,*)'d#sigma-/dm_{t2}--[pb/gev]'
            write(*,*)'m_{t}--[gev]'
          else if (itrans == 5)then
            write(*,*)'m_t3'
            write(*,*)'d#sigma-/dm_{t3}--[pb/gev]'
            write(*,*)'m_{t}--[gev]'
          else if (itrans == 6)then
            write(*,*)'mlt'
            write(*,*)'d#sigma-/dm^{l}_{t}--[pb/gev]'
            write(*,*)'m_{t}--[gev]'
          else if (itrans == 7)then
            write(*,*)'m_ct1'
            write(*,*)'d#sigma-/dm_{t1}--[pb/gev]'
            write(*,*)'m_{t}--[gev]'
          else if (itrans == 8)then
            write(*,*)'m_ct2'
            write(*,*)'d#sigma-/dm_{t2}--[pb/gev]'
            write(*,*)'m_{t}--[gev]'
          else if (itrans == 9)then
            write(*,*)'m_ct3'
            write(*,*)'d#sigma-/dm_{t3}--[pb/gev]'
            write(*,*)'m_{t}--[gev]'
          else if (itrans == 10)then
            write(*,*)'mlct'
            write(*,*)'d#sigma-/dm^{l}_{ct}--[pb/gev]'
            write(*,*)'m^{l}_{ct}--[gev]'
          else
            continue
          end if
          do i=1,ndiv_trans(itrans)
            write(*,*)xtrans(itrans,i),fxtranstot(itrans,i)
          end do
          write(*,*)'END'
        end if
      end do
      ! plot distribution in fl
      if(o_fl == 1)then
        sfxfltot=0d0
        do j=1,ndiv_fl
          fxfltot(j)=0.d0
          do i=1,it
            fxfl(j,i)=fxfl(j,i)*sigma/cnorm(i)/flw
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
      ! plot distribution in cosfl
      if(o_cosfl == 1)then
        sfxcosfltot=0d0
        do j=1,ndiv_cosfl
          fxcosfltot(j)=0.d0
          do i=1,it
            fxcosfl(j,i)=fxcosfl(j,i)*sigma/cnorm(i)/cosflw
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
      ! plot distribution in delta phi
      if(o_dphi == 1)then
        sfxdphitot=0d0
        do j=1,ndiv_dphi
          fxdphitot(j)=0.d0
          do i=1,it
            fxdphi(j,i)=fxdphi(j,i)*sigma/cnorm(i)/dphiw
            fxdphitot(j)=fxdphitot(j)+fxdphi(j,i)
          end do
          sfxdphitot=sfxdphitot+fxdphitot(j)*dphiw
        end do
        write(*,*)'DISTRIBUTION'
        write(*,*)'dphi'
        write(*,*)'d#sigma-/d#delta#phi--[pb]'
        write(*,*)'#delta#phi'
        do i=1,ndiv_dphi
          write(*,*)xdphi(i),fxdphitot(i)
        end do
        write(*,*)'END'
      end if
      ! plot 2d-distribution in delta phi
      if(o_dphi2d == 1)then
        sfxdphi2dtot=0d0
        do i=1,ndiv_dphi
          do j=1,ndiv_rmtt
            fxdphi2dtot(i,j)=0.d0
            do k=1,it
              fxdphi2d(i,j,k)=fxdphi2d(i,j,k)*sigma/cnorm(k)/dphiw/rmttw
              fxdphi2dtot(i,j)=fxdphi2dtot(i,j)+fxdphi2d(i,j,k)
            end do
            sfxdphitot=sfxdphitot+fxdphitot(j)*dphiw
          end do
        end do
        write(*,*)'2D-DISTRIBUTION'
        write(*,*)'dphi2d'
        write(*,*)'d^2#sigma-/d#delta#phi-dm_{tt}--[pb/gev]'
        write(*,*)'#delta#phi--[rad]'
        write(*,*) ndiv_dphi
        write(*,*) dphimin
        write(*,*) dphimax
        write(*,*)'#m_{tt}--[gev]'
        write(*,*) ndiv_rmtt
        write(*,*) rmttmin
        write(*,*) rmttmax
        do i=1,ndiv_dphi
          do j=1,ndiv_rmtt
            write(*,*)xdphi(i),xrmtt(j),fxdphi2dtot(i,j)
          end do
        end do
        write(*,*)'END'
      end if
      ! plot 2d-distributions in delta_phi and all transverse variables
      do itrans=1,ntrans
        if(o_transdp(itrans) == 1)then
          sfxtransdptot(itrans)=0d0
          do i=1,ndiv_dphi
            do j=1,ndiv_trans(itrans)
              fxtransdptot(itrans,i,j)=0.d0
              do k=1,it
                fxtransdp(itrans,i,j,k)=fxtransdp(itrans,i,j,k) &
                *sigma/cnorm(k)/transw(itrans)/dphiw
                fxtransdptot(itrans,i,j)=fxtransdptot(itrans,i,j) &
                +fxtransdp(itrans,i,j,k)
              end do
              sfxtransdptot(itrans)=sfxtransdptot(itrans)+ &
              fxtransdptot(itrans,i,j)*transw(itrans)
            end do
          end do
          write(*,*)'2D-DISTRIBUTION'
          if (itrans == 1)then
            write(*,*)'dphimvis'
            write(*,*)'d#sigma-/dm_{vis}--[pb/gev]'
          else if (itrans == 2)then
            write(*,*)'ht'
            write(*,*)'d#sigma-/dh_{t}--[pb/gev]'
          else if (itrans == 3)then
            write(*,*)'dphim_t1'
            write(*,*)'d#sigma-/dm_{t1}--[pb/gev]'
          else if (itrans == 4)then
            write(*,*)'dphim_t2'
            write(*,*)'d#sigma-/dm_{t2}--[pb/gev]'
          else if (itrans == 5)then
            write(*,*)'dphim_t3'
            write(*,*)'d#sigma-/dm_{t3}--[pb/gev]'
          else if (itrans == 6)then
            write(*,*)'dphimlt'
            write(*,*)'d#sigma-/dm^{l}_{t}--[pb/gev]'
          else if (itrans == 7)then
            write(*,*)'dphim_ct1'
            write(*,*)'d#sigma-/dm_{t1}--[pb/gev]'
          else if (itrans == 8)then
            write(*,*)'dphim_ct2'
            write(*,*)'d#sigma-/dm_{t2}--[pb/gev]'
          else if (itrans == 9)then
            write(*,*)'dphim_ct3'
            write(*,*)'d#sigma-/dm_{t3}--[pb/gev]'
          else if (itrans == 10)then
            write(*,*)'dphimlct'
            write(*,*)'d#sigma-/dm^{l}_{ct}--[pb/gev]'
          else
            continue
          end if
          write(*,*)'#delta#phi'
          write(*,*) ndiv_dphi
          write(*,*) dphimin
          write(*,*) dphimax
          if (itrans == 1)then
            write(*,*)'m_{vis}--[gev]'
          else if (itrans == 2)then
            write(*,*)'h_{t}--[gev]'
          else if (itrans == 3)then
            write(*,*)'m_{t}--[gev]'
          else if (itrans == 4)then
            write(*,*)'m_{t}--[gev]'
          else if (itrans == 5)then
            write(*,*)'m_{t}--[gev]'
          else if (itrans == 6)then
            write(*,*)'m_{t}--[gev]'
          else if (itrans == 7)then
            write(*,*)'m_{t}--[gev]'
          else if (itrans == 8)then
            write(*,*)'m_{t}--[gev]'
          else if (itrans == 9)then
            write(*,*)'m_{t}--[gev]'
          else if (itrans == 10)then
            write(*,*)'m^{l}_{ct}--[gev]'
          else
            continue
          end if
          write(*,*)ndiv_trans(itrans)
          write(*,*)transmin(itrans)
          write(*,*)transmax(itrans)
          do i=1,ndiv_dphi
            do j=1,ndiv_trans(itrans)
              write(*,*)xdphi(i),xtrans(itrans,j) &
              ,fxtransdptot(itrans,i,j)
            end do
          end do
          write(*,*)'END'
        end if
      end do
      ! plot distribution in cost5
      if(o_cost5 == 1)then
        sfxcost5tot=0d0
        do j=1,ndiv_cost5
          fxcost5tot(j)=0.d0
          do i=1,it
            fxcost5(j,i)=fxcost5(j,i)*sigma/cnorm(i)/cost5w
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
      ! plot distribution in cost7
      if(o_cost7 == 1)then
        sfxcost7tot=0d0
        do j=1,ndiv_cost7
          fxcost7tot(j)=0.d0
          do i=1,it
            fxcost7(j,i)=fxcost7(j,i)*sigma/cnorm(i)/cost7w
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
      ! plot distribution in ct7ct5
      if(o_ct7ct5 == 1)then
        sfxct7ct5tot=0d0
        do j=1,ndiv_ct7ct5
          fxct7ct5tot(j)=0.d0
          do i=1,it
            fxct7ct5(j,i)=fxct7ct5(j,i)*sigma/cnorm(i)/ct7ct5w
            fxct7ct5tot(j)=fxct7ct5tot(j)+fxct7ct5(j,i)
          end do
          sfxct7ct5tot=sfxct7ct5tot+fxct7ct5tot(j)*ct7ct5w
        end do
        write(*,*)'DISTRIBUTION'
        write(*,*)'ct7ct5'
        write(*,*) &
        'd^{2}#sigma-/d(cos#theta^{*}_{+}cos#theta^{*}_{-})--[pb]'
        write(*,*)'cos#theta_{+}cos#theta_{-}'
        do i=1,ndiv_ct7ct5
          write(*,*)xct7ct5(i),fxct7ct5tot(i)
        end do
        write(*,*)'END'
      end if
      ! plot distributions in all asymmetries
      if((o_sigp == 1) .and. (o_sigm == 1))then
        do jasy=1,nasym
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
            write(*,*)'ASYMMETRY'
            if(jasy == 1)then
              write(*,*)'ALL'
              write(*,*)'a_{ll}'
            else if(jasy == 2)then
              write(*,*)'AL'
              write(*,*)'a_{l}'
            else if(jasy == 3)then
              write(*,*)'APV'
              write(*,*)'a_{pv}'
            else if(jasy == 4)then
              write(*,*)'AFB'
              write(*,*)'a_{fb}'
            else if(jasy == 5)then
              write(*,*)'AFBSTAR'
              write(*,*)'a_{fb^{*}}'
            else if(jasy == 6)then
              write(*,*)'AtRFB'
              write(*,*)'a^{t}_{rfb}'
            else if(jasy == 7)then
              write(*,*)'AttbRFB'
              write(*,*)"a^{b\bar{b}}_{rfb}"
            else if(jasy == 8)then
              write(*,*)'ARFB'
              write(*,*)"a_{rfb}"
            else if(jasy == 9)then
              write(*,*)'A_l'
              write(*,*)'a_{l}'
            end if
            write(*,*)'m_{tt}'
            ndiv_sig=(ndiv_sigm+ndiv_sigp)/2
            do i=1,ndiv_sig
              if(fxsigptot(jasy,i)+fxsigmtot(jasy,i) == 0.d0)then
                write(*,*)(xsigm(i)+xsigp(i))/2.d0,0.d0,0.d0,0.d0
                !           snorm(jasy)=snorm(jasy)+0.d0
              else
                write(*,*)(xsigm(i)+xsigp(i))/2.d0, &
                (fxsigptot(jasy,i)-fxsigmtot(jasy,i))/ &
                (fxsigptot(jasy,i)+fxsigmtot(jasy,i)), &
                fxsigptot(jasy,i),fxsigmtot(jasy,i)
                !             snorm(jasy)=snorm(jasy)+
                !    &               (fxsigptot(jasy,i)-fxsigmtot(jasy,i))/
                !    &               (fxsigptot(jasy,i)+fxsigmtot(jasy,i))
                !    &               *fxrmtttot(i)*rmttw/sigma
              end if
            end do
            write(*,*)'test'
            asym_int(jasy)=(sfxsigptot(jasy)-sfxsigmtot(jasy))/(sfxsigptot(jasy)+sfxsigmtot(jasy))
            write(*,*)'END'
            !         write(*,*)'(total asymmetry:',asym_int(jasy),')'
            !         write(*,*)'(integrated asymmetry:',snorm(jasy),' )'
          end if
        end do
      end if
    end subroutine print_distributions



    subroutine  check_distributions
      real :: diff_max=1e-12
      n_error=0
      do ip=3,nfinal
        if(o_pt(ip) == 1)then
          if(abs(cross-sfxpttot(ip))>diff_max)then
            write(*,*)'pt',ip,' error:',sfxpttot(ip)
            n_error=n_error+1
          end if
        end if
        if(o_eta(ip) == 1)then
          if(abs(cross-sfxetatot(ip))>diff_max)then
            write(*,*)'eta',ip,' error:',sfxetatot(ip)
            n_error=n_error+1
          end if
        end if
        if(o_phi(ip) == 1)then
          if(abs(cross-sfxphitot(ip))>diff_max)then
            write(*,*)'phi',ip,' error:',sfxphitot(ip)
            n_error=n_error+1
          end if
        end if
      end do
      if(o_etmiss == 1)then
        if(abs(cross-sfxetmisstot)>diff_max)then
          write(*,*)'etmiss error:',sfxetmisstot
          n_error=n_error+1
        end if
      end if
      if(o_pt356 == 1)then
        if(abs(cross-sfxpt356tot)>diff_max)then
          write(*,*)'pt356 error:',sfxpt356tot
          n_error=n_error+1
        end if
      end if
      if(o_eta356 == 1)then
        if(abs(cross-sfxeta356tot)>diff_max)then
          write(*,*)'eta356 error:',sfxeta356tot
          n_error=n_error+1
        end if
      end if
      if(o_phi356 == 1)then
        if(abs(cross-sfxphi356tot)>diff_max)then
          write(*,*)'phi356 error:',sfxphi356tot
          n_error=n_error+1
        end if
      end if
      if(o_pt478 == 1)then
        if(abs(cross-sfxpt478tot)>diff_max)then
          write(*,*)'pt478 error:',sfxpt478tot
          n_error=n_error+1
        end if
      end if
      if(o_eta478 == 1)then
        if(abs(cross-sfxeta478tot)>diff_max)then
          write(*,*)'eta478 error:',sfxeta478tot
          n_error=n_error+1
        end if
      end if
      if(o_phi478 == 1)then
        if(abs(cross-sfxphi478tot)>diff_max)then
          write(*,*)'phi478 error:',sfxphi478tot
          n_error=n_error+1
        end if
      end if
      if(o_rmtt == 1)then
        if(abs(cross-sfxrmtttot)>diff_max)then
          write(*,*)'rmtt error:',sfxrmtttot
          n_error=n_error+1
        end if
      end if

      if(o_beta == 1)then
        if(abs(cross-sfxbetatot)>diff_max)then
          write(*,*)'beta error:',sfxbetatot
          n_error=n_error+1
        end if
      end if

      if(o_cost == 1)then
        if(abs(cross-sfxcosttot)>diff_max)then
          write(*,*)'cost error:',sfxcosttot
          n_error=n_error+1
        end if
      end if
      if(o_et == 1)then
        if(abs(cross-sfxettot)>diff_max)then
          write(*,*)'et error:',sfxettot
          n_error=n_error+1
        end if
      end if
      if(o_fl == 1)then
        if(abs(cross-sfxfltot)>diff_max)then
          write(*,*)'fl error:',sfxfltot
          n_error=n_error+1
        end if
      end if
      if(o_cosfl == 1)then
        if(abs(cross-sfxcosfltot)>diff_max)then
          write(*,*)'cosfl error:',sfxcosfltot
          n_error=n_error+1
        end if
      end if
      if(o_dphi == 1)then
        if(abs(cross-sfxdphitot)>diff_max)then
          write(*,*)'dphi error:',sfxdphitot
          n_error=n_error+1
        end if
      end if 
      do iasy=1,nasym
        if(o_asym(iasy) == 0)then
          continue
        else
          if(abs(atot(iasy)-asym_int(iasy))>diff_max)then
            write(*,*)'a error:',iasy,asym_int(iasy)
            n_error=n_error+1
          end if
        end if
      end do
      write(*,*)'integration errors:',n_error
      write(*,*)'CLOSE'

    end subroutine check_distributions

end module distributions
