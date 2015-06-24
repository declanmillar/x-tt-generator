program zprime

  ! Calculates the cross section, asymmetries 
  ! and generates events or distributions for

  ! pp -> tt,
  ! pp -> tt -> bbWW -> bbllnunu.

  ! Uses adapted MadGraph functions.
  ! Uses HELAS subroutines.
  ! Uses CTEQ6 and MRS99 PDF subroutines.

  ! Authors: Declan Millar, Stefano Moretti.

  use mathematics, only: pi
  use configuration
  use modelling
  use scattering
  use kinematics
  use integration

  implicit none

  external dsigma

  real :: sigma_ee, sigma_emu, sigma_eq, sigma_qq
  real :: error_sigma_ee, error_sigma_emu, error_sigma_eq, error_sigma_qq
  real :: sigma_pol_tot(-1:1,-1:1), error_pol_tot(-1:1,-1:1)
  real :: sigma_fb_tot(-1:1,n_fb_asymmetries), error_fb_tot(-1:1,n_fb_asymmetries)
  real :: all, error_all, al, error_al, apv, error_apv
  real :: afb(n_fb_asymmetries), error_afb(n_fb_asymmetries)
  real :: chi2_sigma, error_sigma, stantot
  real :: alfas, qcdl4
  ! branching ratio for t->benu=bmuv=btaumu (1/9 with QCD corrections)
  real, parameter :: brtbln = 0.10779733d0
  ! branching ratio for t->bqq is leftover after 3 l generations
  real, parameter :: brtbqq = 1 - brtbln*3
  integer :: ndimensions, iab, iasy, ifb, icteq, lam3, lam4
  integer :: i, j, k
  integer :: today(3), now(3)
  double precision :: start_time, finish_time

  call cpu_time(start_time)

  call read_config

  call modify_config

  print*, "Output Ntuple will be written to ", ntuple_file
  call rootinit(ntuple_file)

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

  ! Match Lambda QCD for PDF fits.
  lambdaqcd4 = qcdl4

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

  if (final_state == 0) then
    ! multiply by branching ratios
    fac_ee = brtbln*brtbln
    fac_emu = 2*brtbln*brtbln
    fac_eq = brtbln*brtbqq
    fac_qq = brtbqq*brtbqq
  else if (final_state > 0) then
    ! scale dilepton to other classifications
    fac_ee = 1
    fac_emu = 2
    fac_eq = 12
    fac_qq = 36  
  end if

  if (final_state == 0) then
    m3 = tmass
    m4 = tmass
    m5 = 0.d0
    m6 = 0.d0
    m7 = 0.d0
    m8 = 0.d0
  else if (final_state > 0) then
    m3 = bmass
    m4 = bmass
    m5 = emass
    m6 = nuemass
    m7 = emass
    m8 = nuemass
  end if

  ! VEGAS parameters
  if (use_rambo == 0) then
    if (final_state == 0) then
      ! integrates on:
      ! x(3) = (x1 - tau)/(1 - tau),
      ! x(2) = (ecm - m3 - m4)/(ecm_max - m3 - m4),
      ! x(1) = cos(theta3_cm)
      ndimensions = 3

      ! limits:
      do i = 3, 2, -1
        xl(i) = 0.d0
        xu(i) = 1.d0
      end do
      do i = 1, 1
        xl(i) = -1.d0
        xu(i) = 1.d0
      end do

    else if (final_state > 0) then
      ! integrates on:         
      ! x(15) = (x1 - tau)/(1 - tau),
      ! x(14) = (ecm - m3 - m4 - m5 - m6 - m7 - m8)
      !        /(ecm_max - m3 - m4 - m5 - m6 - m7 - m8),
      ! x(13) = (xx356 - xx356min)/(xx356max - xx356min),
      ! where xx356 = arctg((m356**2 - m3**2)/m3/gamt),
      ! or x(13) = (m356 - m356min)/(m356max - m356min),
      ! where m356min = m3 + m5 + m6, m356max = ecm_max - m4 - m7 - m8
      ! x(12) = (xx478 - xx478min)/(xx478max - xx478min),
      ! where xx478 = arctg((m478**2 - m3**2)/m3/gamt),
      ! or x(12) = (m478 - m478min)/(m478max - m478min),
      ! where m478min = m4 + m7 + m8, m478max = ecm_max - m356
      ! x(11) = (xx56 - xx56min)/(xx56max - xx56min),
      ! where xx56 = arctg((m56**2 - rm_w**2)/rm_w/gamw),
      ! or x(11) = (m56 - m56min)/(m56max - m56min),
      ! where m56min = m5 + m6, m56max = m356 - m3
      ! x(10) = (xx78 - xx78min)/(xx78max - xx78min),
      ! where xx78 = arctg((m78**2 - rm_w**2)/rm_w/gamw),
      ! or x(10) = (m78 - m78min)/(m78max - m78min),
      ! where m78min = m7 + m8, m78max = m478 - m4
      ! x(9) = cos(theta_cm_356) = -cos(theta_cm_478),
      ! x(8) = cos(theta56_cm_356),
      ! x(7) = cos(theta78_cm_478),
      ! x(6) = cos(theta5_cm_56),
      ! x(5) = cos(theta7_cm_78),
      ! x(4) = fi56_cm_356,
      ! x(3) = fi78_cm_478,
      ! x(2) = fi5_cm_56,
      ! x(1) = fi8_cm_78;
      ndimensions = 15

      ! set integration limits:
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
  else if (use_rambo == 1) then
    ! integrates on:
    ! x(2) = (x1 - tau)/(1 - tau),
    ! x(1) = (ecm - m3 - m4 - m5 - m6 - m7 - m8)/
    !        (ecm_max - m3 - m4 - m5 - m6 - m7 - m8)
    ndimensions = 2
    do i = 2, 1, -1
      xl(i) = 0.d0
      xu(i) = 1.d0
    end do
  end if

  ! if nprn<0 no print-out
  nprn = 0   

  ! output information before integration
  print*, '====================================================='
  call idate(today)     ! today(1) = day, (2) = month, (3) = year
  call itime(now)       ! now(1) = hour, (2) = minute, (3) = second
  print*, 'date ', today(3), today(2), today(1)
  print*, 'time ', now(1), now(2), now(3)
  print*, '-----------------------------------------------------'
  print*, 'process'
  if (initial_state == 0) then
    if (final_state == 0) then
      print*, 'pp #rightarrow t#bar{t}', &
               ' #times br(t#rightarrow bl#nu)^{2}'
    else if (final_state == 1) then
      print*, 'pp #rightarrow t#bar{t}', &
               '#rightarrow b#bar{b} w^{+}w^{-}', &
               '#rightarrow b#bar{b} l^{+}l^{-} #nu#bar{#nu}'
    else if (final_state == 2) then
      print*, 'pp #rightarrow t#bar{t}', &
               '#rightarrow b#bar{b} w^{+}w^{-}', &
               '#rightarrow b#bar{b} q#bar{q} l #nu'
    else if (final_state == 3) then
      print*, 'pp #rightarrow t#bar{t}', &
               '#rightarrow b#bar{b} w^{+}w^{-}', &
               "#rightarrow b#bar{b} q#bar{q}q'#bar{q}'"
    end if
  else if (initial_state == 1) then
    if (final_state == 0) then
      print*, 'p#bar{p} #rightarrow t#bar{t}', &
               ' #times br(t#rightarrow bl#nu)^{2}'
    else if (final_state == 1) then
      print*, 'p#bar{p} #rightarrow t#bar{t}', &
               '#rightarrow b#bar{b} w^{+}w^{-}', &
               '#rightarrow b#bar{b} l^{+}l^{-} #nu#bar{#nu}'
    else if (final_state == 2) then
      print*, 'p#bar{p} #rightarrow t#bar{t}', &
               '#rightarrow b#bar{b} w^{+}w^{-}', &
               '#rightarrow b#bar{b} q#bar{q} l #nu'
    else if (final_state == 3) then
      print*, 'p#bar{p} #rightarrow t#bar{t}', &
               '#rightarrow b#bar{b} w^{+}w^{-}', &
               "#rightarrow b#bar{b} q#bar{q}q'#bar{q}'"
    end if
  end if
  print*, '-----------------------------------------------------'
  print*, 'notes'
  print*, 'units: gev'
  print*, 'quarks: all massless except t, b.'
  if (structure_function == 1) print*, 'PDFs: CTEQ6m.'
  if (structure_function == 2) print*, 'PDFs: CTEQ6d.'
  if (structure_function == 3) print*, 'PDFs: CTEQ6l.'
  if (structure_function == 4) print*, 'PDFs: CTEQ6l1.'
  if (structure_function == 5) print*, 'PDFs: MRS99 (cor01).'
  if (structure_function == 6) print*, 'PDFs: MRS99 (cor02).'
  if (structure_function == 7) print*, 'PDFs: MRS99 (cor03).'
  if (structure_function == 8) print*, 'PDFs: MRS99 (cor04).'
  if (structure_function == 9) print*, 'PDFs: MRS99 (cor05).'
  if ((final_state >= 1) .and. (use_nwa == 0)) print*, 'tops: off-shell.'
  if ((final_state >= 1) .and. (use_nwa == 1)) print*, 'tops: nwa.'
  print*, 'bsm model_name: ', model_name
  if (include_qcd == 1) print*, 'QCD: on '
  if (include_qcd == 0) print*, 'QCD: off'
  if (include_ew == 1) print*, 'EW:  on '
  if (include_ew == 0) print*, 'EW:  off'
  if (include_bsm == 1) print*, 'BSM: on '
  if (include_bsm == 0) print*, 'BSM: off'
  if (include_gg == 1) print*, 'gg: on '
  if (include_gg == 0) print*, 'gg: off'
  if (include_qq == 1) print*, 'qq: on '
  if (include_qq == 0) print*, 'qq: off'
  if (interference == 0) print*, 'interference: none'
  if (interference == 1) print*, 'interference: sm'
  if (interference == 2) print*, 'interference: full'
  if (interference == 3) print*, 'interference: no square terms.'
  if (phase_space_only == 1) print*, 'phase space only'
  if (symmetrise_x1x2 == 1) print*, 'symmetrical integration over x1<->x2'
  if (symmetrise_costheta_t == 1) print*, 'symmetrical integration over cost'
  if (include_errors == 1) print*, 'Distribution errors included.'
  if (include_errors == 0) print*, 'Distribution errors excluded.'
  if (use_rambo == 0) print*, 'PS is MANUAL'
  if (use_rambo == 1) print*, 'PS is RAMBO'
  if (map_phase_space == 0) print*, "Phase space is unmapped."
  if (map_phase_space == 1) print*, "Phase space is mapped."
  print*, 'seed: ', seed
  print*, '-----------------------------------------------------'
  print*, 'parameters'
  print*, '#sqrt{s}              ', collider_energy
  print*, 'at |y| <              ', 100
  print*, 'loops a_s evaluated at', nloops
  print*, 'a_{s}(m_{z})          ', alfas(rm_z, lambdaqcd4, nloops)
  print*, '#lambda_{qcd}(4)      ', qcdl4
  print*, 'm_{b}                 ', fmass(12)
  print*, '#gamma_{b}            ', fwidth(12)
  print*, 'm_{t}                 ', fmass(11)
  print*, '#gamma_{t}            ', fwidth(11)
  print*, 'm_{z}                 ', rm_z
  print*, '#gamma_{z}            ', gamma_z
  print*, 'm_{w}                 ', rm_w
  print*, '#gamma_{w}            ', gamma_w
  print*, 'm_{h}                 ', rm_h
  print*, '#gamma_{h}            ', gamma_h
  print*, '-----------------------------------------------------'
  print*, 'zprime parameters'
  do i = 1, 5
    if (mass_zp(i) > 0) then
      print*, 'z#prime               ', i
      print*, 'm_{z#prime}           ', mass_zp(i)
      print*, '#Gamma_{z#prime}      ', gamzp(i)
      print*, 'Gamma/m               ', gamzp(i)/mass_zp(i)             
      print*, 'g_{L}^{u}             ', gZpu(1,i)
      print*, 'g_{R}^{u}             ', gZpu(2,i)
      print*, 'g_{L}^{d}             ', gZpd(1,i)
      print*, 'g_{R}^{d}             ', gZpd(2,i)
      print*, 'g_{L}^{l}             ', gZpl(1,i)
      print*, 'g_{R}^{l}             ', gZpl(2,i)
      print*, 'g_{L}^{n}             ', gZpn(1,i)
      print*, 'g_{R}^{n}             ', gZpn(2,i)
      print*, 'g_{L}^{t}             ', gZpt(1,i)
      print*, 'g_{R}^{t}             ', gZpt(2,i)
      print*, 'g_{L}^{b}             ', gZpb(1,i)
      print*, 'g_{R}^{b}             ', gZpb(2,i)
      print*, 'g_{L}^{l3}            ', gZpl3(1,i)
      print*, 'g_{R}^{l3}            ', gZpl3(2,i)
      print*, 'g_{L}^{n3}            ', gZpn3(1,i)
      print*, 'g_{R}^{n3}            ', gZpn3(2,i)
    end if
  end do

  ! reset counter
  npoints = 0

  ! reset 
  if (final_state == 0) then
    do i = 1, 20
      resl(i) = 0.d0
      standdevl(i) = 0.d0
      cnorm(i) = 0.d0
      do lam3 = -1, +1, 2
        do lam4 = -1, +1, 2
          sigma_pol(lam3, lam4, i) = 0.d0
          error_pol(lam3, lam4, i) = 0.d0
        end do
      end do
      do ifb = 1, n_fb_asymmetries
        do iab = -1, +1, 2
          sigma_fb(ifb, iab, i) = 0.d0
          error_fb(ifb, iab, i) = 0.d0
        end do
      end do
    end do
  end if

  ! integrate 
  print*, 'Starting integration... '
  it = 0 
  call vegas(ndimensions, dsigma, sigma, error_sigma, chi2_sigma)

  print*, "...complete."  

  ! convert results to different tt classifications
  sigma_ee = sigma*fac_ee
  error_sigma_ee = error_sigma*fac_ee
  sigma_emu = sigma*fac_emu
  error_sigma_emu = error_sigma*fac_emu
  sigma_eq = sigma*fac_eq
  error_sigma_eq = error_sigma*fac_eq
  sigma_qq = sigma*fac_qq
  error_sigma_qq = error_sigma*fac_qq

  if (sigma == 0.d0) then
    print*, "Error: sigma = 0. Are any gauge sectors active?"
    stop
  else
    print*, "Using ", npoints, " points:"
    if (final_state == 0) then
      print*, "sigma_tt (pb)", "Uncertainty (pb)"
      print*, sigma, error_sigma
    end if
    print*, "sigma_ee (pb)", "Uncertainty (pb)"
    print*, sigma_ee, error_sigma_ee
    print*, "sigma_emu (pb)", "Uncertainty (pb)"
    print*, sigma_emu, error_sigma_emu
    print*, "sigma_eq (pb)", "Uncertainty (pb)"
    print*, sigma_eq, error_sigma_eq
    print*, "sigma_qq (pb)", "Uncertainty (pb)"
    print*, sigma_qq, error_sigma_qq
  end if

  print*, "Calculating factor to re-weight for different iterations."
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
  print*, "...complete."

  ! write sigma and cnorms to a txt file
  i = len(ntuple_file)
  do while(ntuple_file(i:i) == '')
    i = i-1
  end do
  open(unit = 11, file = ntuple_file(1:i)//".txt", status = "replace", action = "write")
  write(11,*) sigma
  do i = 1, it
  	write(11,*) cnorm(i)
  end do
  close(11)

  print*, "Collating polar cross sections..."
  if (final_state == 0) then
    do lam3 = -1, +1, 2
      do lam4 = -1, +1, 2
        sigma_pol_tot(lam3,lam4) = 0.d0
        error_pol_tot(lam3,lam4) = 0.d0
        do i = 1, it
          sigma_pol(lam3,lam4,i) = sigma_pol(lam3,lam4,i)*sigma/cnorm(i)
          error_pol(lam3,lam4,i) = sigma_pol(lam3,lam4,i)*error_sigma/cnorm(i)
          sigma_pol_tot(lam3,lam4) = sigma_pol_tot(lam3,lam4) + sigma_pol(lam3,lam4,i)
          error_pol_tot(lam3,lam4) = sigma_pol_tot(lam3,lam4) + error_pol(lam3,lam4,i)
        end do
        error_pol_tot(lam3,lam4) = error_pol_tot(lam3,lam4)/sigma_pol_tot(lam3,lam4)
      end do
    end do
  end if
  print*, "...complete."

  print*, "Calculating polar asymmetries..."
  if (final_state == 0) then
    all = (sigma_pol_tot(+1, +1) - sigma_pol_tot(+1, -1) &
         - sigma_pol_tot(-1, +1) + sigma_pol_tot(-1, -1))/sigma
    error_all = (sigma_pol_tot(+1, +1) + sigma_pol_tot(+1, -1) &
               + sigma_pol_tot(-1, +1) + sigma_pol_tot(-1, -1))/4.d0*all

    al = (sigma_pol_tot(-1, -1) - sigma_pol_tot(+1, -1) &
        + sigma_pol_tot(-1, +1) - sigma_pol_tot(+1, +1))/sigma
    error_al = (sigma_pol_tot(-1, -1) + sigma_pol_tot(+1, -1) &
              + sigma_pol_tot(-1, +1) + sigma_pol_tot(+1, +1))/4.d0*al

    apv = (sigma_pol_tot(-1, -1) - sigma_pol_tot(+1, +1))/sigma/2.d0
    error_apv = (sigma_pol_tot(-1, -1) + sigma_pol_tot(+1, +1))/2.d0*apv
  end if
  print*, "...complete."

  print*, "Printing total asymmetries..."
  print*, "total asymmetries"
  if (final_state == 0) then
    print*, "ALL:                    uncertainty:"
    print*, all, error_all
    print*, "AL:                     uncertainty:"
    print*, al, error_al
    print*, "APV:                    uncertainty:"
    print*, apv, error_apv
  end if

  print*, "...complete."

  call rootclose
  call cpu_time(finish_time)
  print '(" Program completed in ",f10.6," seconds.")', finish_time - start_time
  stop
end program zprime
