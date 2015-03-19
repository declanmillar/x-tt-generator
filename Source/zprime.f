program zprime

  ! calculates the cross section and generates distributions for
  ! pp -> tt,
  ! pp -> tt -> bw^+bbarw^- -> bbbare^+nue^-nubarc
  ! (future:) pp -> bw^+bbarw^- -> b bbar e^+ nu qqbar'
  ! uses adapted madgraph functions.
  ! uses cteq6 and mrs99 pdf subroutines.
  ! authors: declan millar, stefano moretti.

  use configuration
  use quantum_field_theory
  use kinematics
  use distributions
  use integration

  implicit none

  external differential_cross_section

  real :: poltot(-1:1, -1:1), polchi(-1:1, -1:1)
  real :: spattot(n_fb_asymmetries, -1:1), spatchi(n_fb_asymmetries, -1:1)

  real :: avgi, chi2a, cross, error, sd, stantot, alfas, qcdl4
  integer :: iab, ndimensions, iasy, icteq, iphel,ispat,jphel

  ! branching ratio for t->bev=bmv=btv (with qcd corrections?)
  real :: brtbln = 0.10779733d0

  ! branching ratio for t->bev=bmv=btv=1/9 (tree level)
  ! real brtbln/0.11111111d0/

  ! branching ratio for t->beq=bmq=btq=6/9 (tree level)
  ! real :: brtbeq=0.66666666d0

  integer :: i, j, k
  integer :: today(3), now(3)

  call read_config

  call modify_config

  call initialise_distributions

  s = collider_energy*collider_energy

  ! (pdfs are intrinsically linked to the value of lamda_qcd; alpha_qcd)
  if (structure_function == 1) qcdl4 = 0.326d0
  if (structure_function == 2) qcdl4 = 0.326d0
  if (structure_function == 3) qcdl4 = 0.326d0
  if (structure_function == 4) qcdl4 = 0.215d0
  if (structure_function == 5) qcdl4 = 0.300d0
  if (structure_function == 6) qcdl4 = 0.300d0
  if (structure_function == 7) qcdl4 = 0.300d0
  if (structure_function == 8) qcdl4 = 0.229d0
  if (structure_function == 9) qcdl4 = 0.383d0

  ! qcdl4 is qcd lambda4 (to match pdf fits).
  rlambdaqcd4 = qcdl4

  ! initialise cteq grids.
  if (structure_function <= 4) then
    icteq = structure_function
    call setctq6(icteq)
  end if

  ! use appropriately evolved alphas.
  if (structure_function == 1) nloops = 2
  if (structure_function == 2) nloops = 2
  if (structure_function == 3) nloops = 1
  if (structure_function == 4) nloops = 1
  if (structure_function == 5) nloops = 1
  if (structure_function == 6) nloops = 1
  if (structure_function == 7) nloops = 1
  if (structure_function == 8) nloops = 1
  if (structure_function == 9) nloops = 1

  ! initialise madgraph - masses and coupling constants of particles
  call initialise_standard_model
  call initialise_zprimes

  ! dimensions
  if (final_state == 0) then
    ndimensions = 3
  else if (final_state > 0) then
    ndimensions = 15
  end if

  ! if nprn<0 no print-out
  nprn = 0

  if (final_state == 0) then
    m3 = fmass(11)
    m4 = m3
    m5 = 0.d0
    m6 = 0.d0
    m7 = 0.d0
    m8 = 0.d0

    ! integrates on:
    ! x(3) = (x1 - tau)/(1 - tau),
    ! x(2) = (ecm - m3 - m4)/(ecm_max - m3 - m4),
    ! x(1) = cos(theta3_cm)

    ! limits:
    do i = 3, 2, -1
      xl(i) = 0.d0
      xu(i) = 1.d0
    end do
    do i = 1, 1
      xl(i) = -1.d0
      xu(i) = 1.d0
    end do

  else if (final_state >= 1) then
    m3 = fmass(12)
    m4 = m3
    m5 = 0.d0
    m6 = 0.d0
    m7 = 0.d0
    m8 = 0.d0

    ! integrates on:         
    ! x(15) = (x1 - tau)/(1 - tau),
    ! x(14) = (ecm - m3 - m4 - m5 - m6 - m7 - m8)
    !        /(ecm_max - m3 - m4 - m5 - m6 - m7 - m8),
    ! x(13) = (xx356 - xx356min)/(xx356max - xx356min),
    ! where xx356 = arctg((m356**2 - m3**2)/m3/gamt),
    ! x(12) = (xx478 - xx478min)/(xx478max - xx478min),
    ! where xx478 = arctg((m478**2 - m3**2)/m3/gamt),
    ! x(11) = (xx56 - xx56min)/(xx56max - xx56min),
    ! where xx56 = arctg((m56**2 - rm_w**2)/rm_w/gamw),
    ! x(10) = (xx78 - xx78min)/(xx78max - xx78min),
    ! where xx78 = arctg((m78**2 - rm_w**2)/rm_w/gamw),
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
    do i = 15, 14, -1
      xl(i) = 0.d0
      xu(i) = 1.d0
    end do
    do i = 13, 10, -1
      xl(i) = 0.d0
      xu(i) = 1.d0
    end do
    do i = 9, 5, -1
      xl(i) = -1.d0
      xu(i) = 1.d0
    end do
    do i = 4, 1, -1
      xl(i) = 0.d0
      xu(i) = 2.d0*pi
    end do
  end if

  call generate_bins

  ! output information before integration
  print *, '====================================================='
  call idate(today)     ! today(1) = day, (2) = month, (3) = year
  call itime(now)       ! now(1) = hour, (2) = minute, (3) = second
  print *, 'date ', today(3), today(2), today(1)
  print *, 'time ', now(1), now(2), now(3)
  print *, '-----------------------------------------------------'
  print *, 'process'
  if (initial_state == 0) then
    if (final_state == 0) then
      print *, 'pp #rightarrow t#bar{t}', &
               ' #times br(t#rightarrow bl#nu)^{2}'
    else if (final_state == 1) then
      print *, 'pp #rightarrow t#bar{t}', &
               '#rightarrow b#bar{b} w^{+}w^{-}', &
               '#rightarrow b#bar{b} l^{+}l^{-} #nu#bar{#nu}'
    else if (final_state == 2) then
      print *, 'pp #rightarrow t#bar{t}', &
               '#rightarrow b#bar{b} w^{+}w^{-}', &
               '#rightarrow b#bar{b} q#bar{q} l #nu'
    else if (final_state == 3) then
      print *, 'pp #rightarrow t#bar{t}', &
               '#rightarrow b#bar{b} w^{+}w^{-}', &
               "#rightarrow b#bar{b} q#bar{q}q'#bar{q}'"
    end if
  else if (initial_state == 1) then
    if (final_state == 0) then
      print *, 'p#bar{p} #rightarrow t#bar{t}', &
               ' #times br(t#rightarrow bl#nu)^{2}'
    else if (final_state == 1) then
      print *, 'p#bar{p} #rightarrow t#bar{t}', &
               '#rightarrow b#bar{b} w^{+}w^{-}', &
               '#rightarrow b#bar{b} l^{+}l^{-} #nu#bar{#nu}'
    else if (final_state == 2) then
      print *, 'p#bar{p} #rightarrow t#bar{t}', &
               '#rightarrow b#bar{b} w^{+}w^{-}', &
               '#rightarrow b#bar{b} q#bar{q} l #nu'
    else if (final_state == 3) then
      print *, 'p#bar{p} #rightarrow t#bar{t}', &
               '#rightarrow b#bar{b} w^{+}w^{-}', &
               "#rightarrow b#bar{b} q#bar{q}q'#bar{q}'"
    end if
  end if
  print *, '-----------------------------------------------------'
  print *, 'notes'
  print *, 'units: gev'
  print *, 'quarks: all massless except t, b.'
  if (structure_function == 1) print *, 'pdfs: cteq6m.'
  if (structure_function == 2) print *, 'pdfs: cteq6d.'
  if (structure_function == 3) print *, 'pdfs: cteq6l.'
  if (structure_function == 4) print *, 'pdfs: cteq6l1.'
  if (structure_function == 5) print *, 'pdfs: mrs99 (cor01).'
  if (structure_function == 6) print *, 'pdfs: mrs99 (cor02).'
  if (structure_function == 7) print *, 'pdfs: mrs99 (cor03).'
  if (structure_function == 8) print *, 'pdfs: mrs99 (cor04).'
  if (structure_function == 9) print *, 'pdfs: mrs99 (cor05).'
  if ((final_state >= 1) .and. (use_nwa == 0)) print *, 'tops: off-shell.'
  if ((final_state >= 1) .and. (use_nwa == 1)) print *, 'tops: nwa.'
  print *, 'bsm model_name: ', model_name
  if (include_qcd == 1) print *, 'qcd: on '
  if (include_qcd == 0) print *, 'qcd: off'
  if (include_ew == 1) print *, 'ew:  on '
  if (include_ew == 0) print *, 'ew:  off'
  if (include_bsm == 1) print *, 'bsm: on '
  if (include_bsm == 0) print *, 'bsm: off'
  if (interference == 0) print *, 'interference: none'
  if (interference == 1) print *, 'interference: sm'
  if (interference == 2) print *, 'interference: full'
  if (interference == 3) print *, 'interference: no square terms.'
  if (phase_space_only == 1) print *, 'phase space only'
  if (symmetrise_x1x2 == 1) print *, 'symmetrical integration over x1<->x2'
  if (symmetrise_costheta_t == 1) print *, 'symmetrical integration over cost'
  print *, 'seed: ', seed
  print *, '-----------------------------------------------------'
  print *, 'parameters'
  print *, '#sqrt{s}              ', collider_energy
  print *, 'at |y| <              ', abs(ytmax)
  print *, 'loops a_s evaluated at', nloops
  print *, 'a_{s}(m_{z})          ', alfas(rm_z, rlambdaqcd4, nloops)
  print *, '#lambda_{qcd}(4)      ', qcdl4
  print *, 'm_{b}                 ', fmass(12)
  print *, '#gamma_{b}            ', fwidth(12)
  print *, 'm_{t}                 ', fmass(11)
  print *, '#gamma_{t}            ', fwidth(11)
  print *, 'm_{z}                 ', rm_z
  print *, '#gamma_{z}            ', gamma_z
  print *, 'm_{w}                 ', rm_w
  print *, '#gamma_{w}            ', gamma_w
  print *, 'm_{h}                 ', rm_h
  print *, '#gamma_{h}            ', gamma_h
  print *, '-----------------------------------------------------'
  print *, 'zprime parameters'
  do i = 1, 5
    if (rmzp(i) > 0) then
      print *, 'z#prime               ', i
      print *, 'm_{z#prime}           ', rmzp(i)
      print *, '#gamma_{z#prime}      ', gamzp(i)
      print *, 'g_{p}                 ', gp(i)
      print *, 'g_{v}^{u}             ', gv_u(i)
      print *, 'g_{a}^{u}             ', ga_u(i)
      print *, 'g_{v}^{d}             ', gv_d(i)
      print *, 'g_{a}^{d}             ', ga_d(i)
      print *,
    end if
  end do
  print *, '-----------------------------------------------------'
  print *, 'cuts'
  print *, '-----------------------------------------------------'

  ! reset counter
  npoints = 0

  ! reset 
  if (final_state == 0) then
    do i = 1, 20
      resl(i) = 0.d0
      standdevl(i) = 0.d0
      cnorm(i) = 0.d0
      do iphel = -1, +1, 2
        do jphel = -1, +1, 2
          xsec_polar(i, iphel, jphel) = 0.d0
          error_polar(i, iphel, jphel) = 0.d0
        end do
      end do
      do ispat = 1, n_fb_asymmetries
        do iasy = -1, +1, 2
          xsec_fb(ispat, i, iasy) = 0.d0
          error_fb(ispat, i, iasy) = 0.d0
        end do
      end do
    end do
  end if

  ! integrate 
  print *, 'starting integration'
  it = 0 
  call vegas(ndimensions, differential_cross_section, avgi, sd, chi2a)

  if (final_state == 0 .and. use_branching_ratio == 1) then
    ! multiply by branching ratios
    avgi = avgi*brtbln*brtbln
    sd = sd*brtbln*brtbln
  end if

  ! collect total cross-section
  cross = avgi
  sigma = avgi
  error = sd

  ! print integrated cross section
  print *, ''
  print *, 'integrated cross section'
  if (cross == 0.d0) then
    print *, 'sigma = 0  ! check permitted gauge sectors.'
    stop
  else
    print *, 'sigma (pb)', 'error (same units)'
    print *, cross, error
    print *, '(using ', npoints, ' points)'
  end if

  ! re-weight distributions for different iterations
  stantot = 0.d0
  do i = 1, it
    stantot = stantot + 1.d0/standdevl(i)/standdevl(i)
  end do
  do i = 1, it
    standdevl(i) = standdevl(i)*standdevl(i)*stantot
  end do
  do i = 1, it
    cnorm(i) = resl(i)*standdevl(i)
  end do

  ! total asymmetries
  ! collect polarised cross sections.
  if (include_asymmetries == 1) then
    if (final_state == 0) then
      do iphel = -1, +1, 2
        do jphel = -1, +1, 2
          do i = 1, it
            xsec_polar(i, iphel, jphel) = xsec_polar(i, iphel, jphel) &
                                          *avgi/cnorm(i)
            error_polar(i, iphel, jphel) = xsec_polar(i, iphel, jphel) &
                                           *sd/cnorm(i)
          end do
          poltot(iphel, jphel) = 0.d0
          polchi(iphel, jphel) = 0.d0
          do i = 1, it
            poltot(iphel, jphel) = poltot(iphel, jphel) &
                                   + xsec_polar(i, iphel, jphel)
            polchi(iphel, jphel) = polchi(iphel, jphel) &
                                   + error_polar(i, iphel, jphel)
          end do
          polchi(iphel, jphel) = polchi(iphel, jphel) &
          /poltot(iphel, jphel)
          !        polchi(iphel,jphel)=
          ! & sqrt(abs(polchi(iphel,jphel)
          ! &         -poltot(iphel,jphel)**2*dfloat(ncall)))
          ! & /dfloat(ncall)
        end do
      end do
    end if

    ! collect unpolarised spatial asymmetry
    do ispat = 1, n_fb_asymmetries
      if (o_asym(ispat + 3) == 0) then
        continue
      else
        do iab = -1,+1, 2
          do i = 1, it
            xsec_fb(ispat, i, iab) = xsec_fb(ispat, i, iab) &
                                     *avgi/cnorm(i)
            error_fb(ispat, i, iab) = xsec_fb(ispat, i, iab) &
                                      *sd/cnorm(i)
          end do
          spattot(ispat, iab) = 0.d0
          spatchi(ispat, iab) = 0.d0
          do i = 1, it
            spattot(ispat, iab) = spattot(ispat, iab) &
                                  + xsec_fb(ispat, i, iab)
            spatchi(ispat, iab) = spatchi(ispat, iab) &
                                  + error_fb(ispat, i, iab)
          end do
          spatchi(ispat, iab) = spatchi(ispat, iab) &
                                /spattot(ispat, iab)
          !         spatchi(iasy)=
          !    & sqrt(abs(spatchi(iasy)
          !    &         -spattot(iasy)**2*dfloat(ncall)))
          !    & /dfloat(ncall)
        end do
      end if
    end do

    ! define asymmetries
    if (final_state == 0) then
      ! all
      atot(1) = (poltot(+1, +1) - poltot(+1, -1) &
                 - poltot(-1, +1) + poltot(-1, -1)) &
                /cross
      atoterr(1) = (polchi(+1, +1) + polchi(+1, -1) &
                    +polchi(-1, +1) + polchi(-1, -1)) &
                   /4.d0*atot(1)
      ! al
      atot(2) = (poltot(-1, -1) - poltot(+1, -1) &
                 + poltot(-1, +1) - poltot(+1, +1)) &
                /cross
      atoterr(2) = (polchi(-1, -1) + polchi(+1, -1) &
                    +polchi(-1, +1) + polchi(+1, +1)) &
                   /4.d0*atot(2)
      ! apv
      atot(3) = (poltot(-1, -1) - poltot(+1, +1)) &
                /cross/2.d0
      atoterr(3) = (polchi(-1, -1) + polchi(+1, +1)) &
                   /2.d0*atot(3)
    end if

    do iasy = 4, n_asymmetries
      ispat = iasy - 3
      if (o_asym(iasy) > 0) then
        atot(iasy) = (spattot(ispat, +1) - spattot(ispat, -1))/cross
        atoterr(iasy) = sd/avgi*atot(iasy)
      end if
    end do

    ! print asymmetries
    print *, 'total asymmetries'
    if (final_state == 0) then
      print *, 'ALL:                    uncertainty (same units):'
      print *, atot(1), atoterr(1)
      print *, 'AL:                     uncertainty (same units):'
      print *, atot(2), atoterr(2)
      print *, 'APV:                    uncertainty (same units):'
      print *, atot(3), atoterr(3)
      print *, 'AFB:                    uncertainty (same units):'
      print *, atot(4), atoterr(4)
      print *, 'AFB*:                   uncertainty (same units):'
      print *, atot(5), atoterr(5)
      print *, 'AtRFB:                  uncertainty (same units):'
      print *, atot(6), atoterr(6)
      print *, "AttbRFB:                uncertainty (same units):"
      print *, atot(7), atoterr(7)
      print *, "ARFB/:                  uncertainty (same units):"
      print *, atot(8), atoterr(8)
    else if (final_state > 0) then
      print *, 'A_l:                    uncertainty (same units):'
      print *, atot(9), atoterr(9)
    end if
  end if

  call print_distributions
  call check_distributions

  stop
end program zprime
