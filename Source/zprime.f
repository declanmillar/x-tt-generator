program zprime

  ! calculates the cross section and generates distributions for
  ! pp -> tt,
  ! pp -> tt -> bw^+bbarw^- -> bbbare^+nue^-nubarc
  ! (future:) pp -> bw^+bbarw^- -> b bbar e^+ nu qqbar'
  ! uses adapted madgraph functions.
  ! uses cteq6 and mrs99 pdf subroutines.
  ! authors: declan millar, stefano moretti.

  use configuration
  use kinematics
  use distributions

  implicit real (a-h,p-z)
  implicit integer (l-o)

  external differential_cross_section

  common/ew/a_em,s2w
  common/fermions/ fmass,     fwidth
  real :: fmass(12), fwidth(12)
  common/vmass1/rm_w,gamma_w,rm_z,gamma_z
  common/vmass2/rm_a,gamma_a,rm_h,gamma_h
  common/zp/rmzp(5),gamzp(5)
  common/zpparam/paramzp(5)
  common/coupzpva/gp(5),gv_d(5),ga_d(5),gv_u(5),ga_u(5)
  common/coupzp/gzpd(2,5),gzpu(2,5)
  common/qcd/rlambdaqcd4,nloops

  ! local variables
  ! flag for zp width specification
  real :: o_width(5)
  ! polarised/hemispherised cross sections
  real :: cnorm(20)
  ! imension snorm(6)  !,ave(4)
  real :: poltot(-1:1,-1:1),polchi(-1:1,-1:1)
  real :: spattot(nspat,-1:1),spatchi(nspat,-1:1)
  real :: sfxpttot(8),sfxetatot(8),sfxphitot(8),sfxycoltot(8)

  real :: sfxsigptot(nasym),sfxsigmtot(nasym)
  real :: asym_int(nasym)
  real :: atot(nasym),atoterr(nasym)
  ! test particular matrix element
  real :: testp1(0:3),testp2(0:3),testp3(0:3),testp4(0:3)

  ! local constant
  integer :: today(3), now(3)
  ! branching ratio for t->bev=bmv=btv (with qcd corrections?)
  real :: brtbln=0.10779733d0
  ! branching ratio for t->bev=bmv=btv=1/9 (tree level)
  ! real brtbln/0.11111111d0/
  ! branching ratio for t->beq=bmq=btq=6/9 (tree level)
  ! real :: brtbeq=0.66666666d0

  integer i, j, k

  call read_config

  call modify_config

  call initialise_distributions

  s=collider_energy*collider_energy

  ! qcdl4 is qcd lambda4 (to match pdf fits).
  ! (pdfs are intrinsically linked to the value of lamda_qcd; alpha_qcd)
  if(o_structure == 1)qcdl4=0.326d0
  if(o_structure == 2)qcdl4=0.326d0
  if(o_structure == 3)qcdl4=0.326d0
  if(o_structure == 4)qcdl4=0.215d0
  if(o_structure == 5)qcdl4=0.300d0
  if(o_structure == 6)qcdl4=0.300d0
  if(o_structure == 7)qcdl4=0.300d0
  if(o_structure == 8)qcdl4=0.229d0
  if(o_structure == 9)qcdl4=0.383d0
  rlambdaqcd4=qcdl4
  ! initialise cteq grids.
  if(o_structure <= 4)then
    icteq=o_structure
    call setctq6(icteq)
  end if

  ! use appropriately evolved alphas.
  if(o_structure == 1)nloops=2
  if(o_structure == 2)nloops=2
  if(o_structure == 3)nloops=1
  if(o_structure == 4)nloops=1
  if(o_structure == 5)nloops=1
  if(o_structure == 6)nloops=1
  if(o_structure == 7)nloops=1
  if(o_structure == 8)nloops=1
  if(o_structure == 9)nloops=1

  ! initialise madgraph - masses and coupling constants of particles
  call initialise_madgraph(o_nwa,model)

  ! vegas parameters
  ! real ::s of integration
  if(ifinal_state == 0)then
    ndim=3
  else if(ifinal_state > 0)then
    ndim=15
  end if
  ! if nprn<0 no print-out
  nprn=0
  if(ifinal_state == 0)then
    rm3=fmass(11)
    rm4=rm3
    rm5=0.d0
    rm6=0.d0
    rm7=0.d0
    rm8=0.d0
    ! integrates on:

    ! x(3) = (x1-tau)/(1-tau),
    ! x(2) = (ecm-rm3-rm4)/(ecm_max-rm3-rm4),
    ! x(1) = cos(theta3_cm)

    ! limits:
    do i=3,2,-1
      xl(i)=0.d0
      xu(i)=1.d0
    end do
    do i=1,1
      xl(i)=-1.d0
      xu(i)=1.d0
    end do

  else if(ifinal_state >= 1)then
    rm3=fmass(12)
    rm4=rm3
    rm5=0.d0
    rm6=0.d0
    rm7=0.d0
    rm8=0.d0
    ! integrates on:
         
    ! x(15) = (x1-tau)/(1-tau),
    ! x(14) = (ecm-rm3-rm4-rm5-rm6-rm7-rm8)
    !        /(ecm_max-rm3-rm4-rm5-rm6-rm7-rm8),
    ! x(13) = (xx356-xx356min)/(xx356max-xx356min),
    ! where xx356 = arctg((rm356**2-rm3**2)/rm3/gamt),
    ! x(12) = (xx478-xx478min)/(xx478max-xx478min),
    ! where xx478 = arctg((rm478**2-rm3**2)/rm3/gamt),
    ! x(11) = (xx56-xx56min)/(xx56max-xx56min),
    ! where xx56 = arctg((rm56**2-rm_w**2)/rm_w/gamw),
    ! x(10) = (xx78-xx78min)/(xx78max-xx78min),
    ! where xx78 = arctg((rm78**2-rm_w**2)/rm_w/gamw),
    ! x(9) = cos(theta_cm_356) = -cos(theta_cm_478)
    ! x(8) = cos(theta56_cm_356),
    ! x(7) = cos(theta78_cm_478),
    ! x(6) = cos(theta5_cm_56),
    ! x(5) = cos(theta7_cm_78),
    ! x(4) = fi56_cm_356,
    ! x(3) = fi78_cm_478,
    ! x(2) = fi5_cm_56,
    ! x(1) = fi8_cm_78;

    ! limits:
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

  ! generate bins
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

  ! output information before integration
  write(*,*)'====================================================='
  call idate(today)     ! today(1)=day, (2)=month, (3)=year
  call itime(now)       ! now(1)=hour, (2)=minute, (3)=second
  write(*,*)'date ',today(3),today(2),today(1)
  write(*,*)'time ',now(1),now(2),now(3)
  write(*,*)'-----------------------------------------------------'
  write(*,*)'process'
  if(initial_state == 0)then
    if(ifinal_state == 0) &
    write(*,*)'pp #rightarrow t#bar{t}', &
    ' #times br(t#rightarrow bl#nu)^{2}'
    if(ifinal_state == 1) &
    write(*,*)'pp #rightarrow t#bar{t}', &
    '#rightarrow b#bar{b} w^{+}w^{-}', &
    '#rightarrow b#bar{b} l^{+}l^{-} #nu#bar{#nu}'
    if(ifinal_state == 2) &
    write(*,*)'pp #rightarrow t#bar{t}', &
    '#rightarrow b#bar{b} w^{+}w^{-}', &
    '#rightarrow b#bar{b} q#bar{q} l #nu'
    if(ifinal_state == 3) &
    write(*,*)'pp #rightarrow t#bar{t}', &
    '#rightarrow b#bar{b} w^{+}w^{-}', &
    "#rightarrow b#bar{b} q#bar{q}q'#bar{q}'"
  else if(initial_state == 1)then
    if(ifinal_state == 0) &
    write(*,*)'p#bar{p} #rightarrow t#bar{t}', &
    ' #times br(t#rightarrow bl#nu)^{2}'
    if(ifinal_state == 1) &
    write(*,*)'p#bar{p} #rightarrow t#bar{t}', &
    '#rightarrow b#bar{b} w^{+}w^{-}', &
    '#rightarrow b#bar{b} l^{+}l^{-} #nu#bar{#nu}'
    if(ifinal_state == 2) &
    write(*,*)'p#bar{p} #rightarrow t#bar{t}', &
    '#rightarrow b#bar{b} w^{+}w^{-}', &
    '#rightarrow b#bar{b} q#bar{q} l #nu'
    if(ifinal_state == 3) &
    write(*,*)'p#bar{p} #rightarrow t#bar{t}', &
    '#rightarrow b#bar{b} w^{+}w^{-}', &
    "#rightarrow b#bar{b} q#bar{q}q'#bar{q}'"
  end if
  write(*,*)'-----------------------------------------------------'
  write(*,*)'notes'
  write(*,*)'units: gev'
  write(*,*)'quarks: all massless except t, b.'
  if(o_structure == 1)write(*,*)'pdfs: cteq6m.'
  if(o_structure == 2)write(*,*)'pdfs: cteq6d.'
  if(o_structure == 3)write(*,*)'pdfs: cteq6l.'
  if(o_structure == 4)write(*,*)'pdfs: cteq6l1.'
  if(o_structure == 5)write(*,*)'pdfs: mrs99 (cor01).'
  if(o_structure == 6)write(*,*)'pdfs: mrs99 (cor02).'
  if(o_structure == 7)write(*,*)'pdfs: mrs99 (cor03).'
  if(o_structure == 8)write(*,*)'pdfs: mrs99 (cor04).'
  if(o_structure == 9)write(*,*)'pdfs: mrs99 (cor05).'
  if((ifinal_state >= 1) .and. (o_nwa == 0))write(*,*)'tops: off-shell.'
  if((ifinal_state >= 1) .and. (o_nwa == 1))write(*,*)'tops: nwa.'
  write(*,*)'bsm model: ',model
  if(include_qcd == 1)write(*,*)'qcd: on '
  if(include_qcd == 0)write(*,*)'qcd: off'
  if(include_ew == 1) write(*,*)'ew:  on '
  if(include_ew == 0) write(*,*)'ew:  off'
  if(include_bsm == 1)write(*,*)'bsm: on '
  if(include_bsm == 0)write(*,*)'bsm: off'
  if(interference == 0)write(*,*)'interference: none'
  if(interference == 1)write(*,*)'interference: sm'
  if(interference == 2)write(*,*)'interference: full'
  if(interference == 3)write(*,*)'interference: no square terms.'
  if(o_m_eq_1 == 1)write(*,*)'phase space only'
  if(o_symx1x2 == 1)write(*,*)'symmetrical integration over x1<->x2'
  if(o_symcost == 1)write(*,*)'symmetrical integration over cost'
  write(*,*)'iseed: ',iseed
  write(*,*)'-----------------------------------------------------'
  write(*,*)'parameters'
  write(*,*)'#sqrt{s}              ',collider_energy
  write(*,*)'at |y| <              ',abs(ytmax)
  write(*,*)'loops a_s evaluated at',nloops
  write(*,*)'a_{s}(m_{z})          ',alfas(rm_z,rlambdaqcd4,nloops)
  write(*,*)'#lambda_{qcd}(4)      ',qcdl4
  write(*,*)'m_{b}                 ',fmass(12)
  write(*,*)'#gamma_{b}            ',fwidth(12)
  write(*,*)'m_{t}                 ',fmass(11)
  write(*,*)'#gamma_{t}            ',fwidth(11)
  write(*,*)'m_{z}                 ',rm_z
  write(*,*)'#gamma_{z}            ',gamma_z
  write(*,*)'m_{w}                 ',rm_w
  write(*,*)'#gamma_{w}            ',gamma_w
  write(*,*)'m_{h}                 ',rm_h
  write(*,*)'#gamma_{h}            ',gamma_h
  write(*,*)'-----------------------------------------------------'
  write(*,*)'zprime parameters'
  do i=1,5
    if(rmzp(i) > 0)then
      write(*,*)'z#prime               ',i
      write(*,*)'m_{z#prime}           ',rmzp(i)
      write(*,*)'o_width:              ',o_width(i)
      write(*,*)'#gamma_{z#prime}      ',gamzp(i)
      write(*,*)'g_{p}                 ',gp(i)
      write(*,*)'g_{v}^{u}             ',gv_u(i)
      write(*,*)'g_{a}^{u}             ',ga_u(i)
      write(*,*)'g_{v}^{d}             ',gv_d(i)
      write(*,*)'g_{a}^{d}             ',ga_d(i)
      write(*,*)
    end if
  end do
  write(*,*)'-----------------------------------------------------'
  write(*,*)'cuts'
  write(*,*)'-----------------------------------------------------'


  ! reset counter
  npoints=0

  ! reset 
  if(ifinal_state == 0)then
    do i=1,20
      resl(i)=0.d0
      standdevl(i)=0.d0
      cnorm(i)=0.d0
      do iphel=-1,+1,2
        do jphel=-1,+1,2
          xsec_polar(i,iphel,jphel)=0.d0
          error_polar(i,iphel,jphel)=0.d0
        end do
      end do
      do ispat=1,nspat
        do iasy=-1,+1,2
          xsec_fb(ispat,i,iasy)=0.d0
          error_fb(ispat,i,iasy)=0.d0
        end do
      end do
    end do
  end if

  ! integrate 
  write(*,*)'starting integration'
  it=0  
  call vegas(ndim,differential_cross_section,avgi,sd,chi2a)

  if(ifinal_state == 0 .and. o_br == 1)then
    ! multiply by branching ratios
    avgi=avgi*brtbln*brtbln
    sd=sd*brtbln*brtbln
  end if

  ! collect total cross-section
  cross=avgi
  error=sd

  ! print integrated cross section
  write(*,*)''
  write(*,*)'integrated cross section'
  if(cross == 0d0)then
    write(*,*)'sigma = 0  ! check permitted gauge sectors.'
    stop
  else
    write(*,*)'sigma (pb)','error (same units)'
    write(*,*)cross,error
    write(*,*)'(using ',npoints,' points)'
  end if
  ! re-weight distributions for different iterations
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

  ! total asymmetries
  ! collect polarised cross sections.
  if(o_asyms == 1)then
    if(ifinal_state == 0)then
      do iphel=-1,+1,2
        do jphel=-1,+1,2
          do i=1,it
            xsec_polar(i,iphel,jphel)=xsec_polar(i,iphel,jphel) &
            *avgi/cnorm(i)
            error_polar(i,iphel,jphel)=xsec_polar(i,iphel,jphel) &
            *sd/cnorm(i)
          end do
          poltot(iphel,jphel)=0.d0
          polchi(iphel,jphel)=0.d0
          do i=1,it
            poltot(iphel,jphel)=poltot(iphel,jphel) &
            +xsec_polar(i,iphel,jphel)
            polchi(iphel,jphel)=polchi(iphel,jphel) &
            +error_polar(i,iphel,jphel)
          end do
          polchi(iphel,jphel)=polchi(iphel,jphel) &
          /poltot(iphel,jphel)
          !        polchi(iphel,jphel)=
          ! & sqrt(abs(polchi(iphel,jphel)
          ! &         -poltot(iphel,jphel)**2*dfloat(ncall)))
          ! & /dfloat(ncall)
        end do
      end do
    end if

    ! collect unpolarised spatial asymmetry
    do ispat=1,nspat
      if(o_asym(ispat+3) == 0)then
        continue
      else
        do iab=-1,+1,2
          do i=1,it
            xsec_fb(ispat,i,iab)=xsec_fb(ispat,i,iab) &
            *avgi/cnorm(i)
            error_fb(ispat,i,iab)=xsec_fb(ispat,i,iab) &
            *sd/cnorm(i)
          end do
          spattot(ispat,iab)=0.d0
          spatchi(ispat,iab)=0.d0
          do i=1,it     ! add up each iteration
            spattot(ispat,iab)=spattot(ispat,iab) &
            +xsec_fb(ispat,i,iab)
            spatchi(ispat,iab)=spatchi(ispat,iab) &
            +error_fb(ispat,i,iab)
          end do
          spatchi(ispat,iab)=spatchi(ispat,iab) &
          /spattot(ispat,iab)
          !         spatchi(iasy)=
          !    & sqrt(abs(spatchi(iasy)
          !    &         -spattot(iasy)**2*dfloat(ncall)))
          !    & /dfloat(ncall)
        end do
      end if
    end do

    ! define asymmetries
    if(ifinal_state == 0)then
      ! all
      atot(1)= &
      +(poltot(+1,+1)-poltot(+1,-1) &
      -poltot(-1,+1)+poltot(-1,-1)) &
      /cross
      atoterr(1)= &
      +(polchi(+1,+1)+polchi(+1,-1) &
      +polchi(-1,+1)+polchi(-1,-1)) &
      /4.d0*atot(1)
      ! al
      atot(2)= &
      +(poltot(-1,-1)-poltot(+1,-1) &
      +poltot(-1,+1)-poltot(+1,+1)) &
      /cross
      atoterr(2)= &
      +(polchi(-1,-1)+polchi(+1,-1) &
      +polchi(-1,+1)+polchi(+1,+1)) &
      /4.d0*atot(2)
      ! apv
      atot(3)= &
      +(poltot(-1,-1)-poltot(+1,+1)) &
      /cross/2.d0
      atoterr(3)= &
      +(polchi(-1,-1)+polchi(+1,+1)) &
      /2.d0*atot(3)
    end if

    do iasy=4,nasym
      ispat=iasy-3
      if(o_asym(iasy) > 0)then
        atot(iasy)= &
        +(spattot(ispat,+1)-spattot(ispat,-1)) &
        /cross
        atoterr(iasy)= &
        +sd/avgi*atot(iasy)
      end if
    end do

    ! print asymmetries
    write(*,*)'total asymmetries'
    if(ifinal_state == 0)then
      write(*,*)'all:                  uncertainty (same units):'
      write(*,*)atot(1),atoterr(1)
      write(*,*)'al:                   uncertainty (same units):'
      write(*,*)atot(2),atoterr(2)
      write(*,*)'apv:                  uncertainty (same units):'
      write(*,*)atot(3),atoterr(3)
      write(*,*)'afb:                 uncertainty (same units):'
      write(*,*)atot(4),atoterr(4)
      write(*,*)'afb*:                   uncertainty (same units):'
      write(*,*)atot(5),atoterr(5)
      write(*,*)'atrfb:                  uncertainty (same units):'
      write(*,*)atot(6),atoterr(6)
      write(*,*)"attbrfb/a:              uncertainty (same units):"
      write(*,*)atot(7),atoterr(7)
      write(*,*)"arfb/a':              uncertainty (same units):"
      write(*,*)atot(8),atoterr(8)
    else if(ifinal_state > 0)then
      write(*,*)'a_l:                  uncertainty (same units):'
      write(*,*)atot(9),atoterr(9)
    end if
  end if

  ! plot distributions

  write(*,*)''
  write(*,*)'histograms'
  do ip=3,8
    ! plot distributions in pt
    if(o_pt(ip) == 1)then
      sfxpttot(ip)=0d0
      do j=1,ndiv_pt(ip)
        fxpttot(ip,j)=0.d0
        do i=1,it
          fxpt(ip,j,i)=fxpt(ip,j,i)*avgi/cnorm(i)/ptw(ip)
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
          fxeta(ip,j,i)=fxeta(ip,j,i)*avgi/cnorm(i)/etaw(ip)
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
          fxphi(ip,j,i)=fxphi(ip,j,i)*avgi/cnorm(i)/phiw(ip)
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
          fxycol(ip,j,i)=fxycol(ip,j,i)*avgi/cnorm(i)/ycolw(ip)
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
        fxetmiss(j,i)=fxetmiss(j,i)*avgi/cnorm(i)/etmissw
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
        fxpt356(j,i)=fxpt356(j,i)*avgi/cnorm(i)/pt356w
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
  ! plot distribution in phi356
  if(o_phi356 == 1)then
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
  ! plot distribution in pt478
  if(o_pt478 == 1)then
    sfxpt478tot=0d0
    do j=1,ndiv_pt478
      fxpt478tot(j)=0.d0
      do i=1,it
        fxpt478(j,i)=fxpt478(j,i)*avgi/cnorm(i)/pt478w
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
  ! plot distribution in phi478
  if(o_phi478 == 1)then
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
  ! plot distribution in mtt
  if(o_rmtt == 1)then
    sfxrmtttot=0d0
    do j=1,ndiv_rmtt
      fxrmtttot(j)=0.d0
      do i=1,it
        fxrmtt(j,i)=fxrmtt(j,i)*avgi/cnorm(i)/rmttw
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
        fxbeta(j,i)=fxbeta(j,i)*avgi/cnorm(i)/betaw
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
  ! plot distribution in et
  if(o_et == 1)then

    sfxettot=0d0
    do j=1,ndiv_et
      fxettot(j)=0.d0
      do i=1,it
        fxet(j,i)=fxet(j,i)*avgi/cnorm(i)/etw
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
        fxdelta_y(j,i)=fxdelta_y(j,i)*avgi/cnorm(i)/delta_yw
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
          *avgi/cnorm(i)/transw(itrans)
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
  ! plot distribution in cosfl
  if(o_cosfl == 1)then
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
  ! plot distribution in delta phi
  if(o_dphi == 1)then
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
          fxdphi2d(i,j,k)=fxdphi2d(i,j,k)*avgi/cnorm(k)/dphiw/rmttw
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
            *avgi/cnorm(k)/transw(itrans)/dphiw
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
  ! plot distribution in cost7
  if(o_cost7 == 1)then
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
  ! plot distribution in ct7ct5
  if(o_ct7ct5 == 1)then
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
            !    &               *fxrmtttot(i)*rmttw/avgi
          end if
        end do
        asym_int(jasy)=(sfxsigptot(jasy)-sfxsigmtot(jasy))/ &
        (sfxsigptot(jasy)+sfxsigmtot(jasy))
        write(*,*)'END'
        !         write(*,*)'(total asymmetry:',asym_int(jasy),')'
        !         write(*,*)'(integrated asymmetry:',snorm(jasy),' )'
      end if
    end do
  end if

  ! check distributions
  diff_max=1e-12
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
    if(abs(cross-sfxbetatot)>diff_max*10)then
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
  write(*,*)'close'
  stop
end program zprime
