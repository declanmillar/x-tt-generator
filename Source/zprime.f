program zprime

  ! Calculates the cross section, asymmetries, and generates events for

  ! pp -> ff,
  ! pp -> tt -> bbWW -> bbllvv.

  ! Uses :
  !  - adapted MadGraph functions
  !  - HELAS subroutines
  !  - VEGAS Monte Carlo integration
  !  - CTEQ6 and MRS99 PDF subroutines
  !  - RootTuple for filling ntuples

  ! Authors: Declan Millar, Stefano Moretti.

  use mathematics, only: pi
  use configuration
  use modelling
  use scattering
  use integration

  implicit none

  real :: sigma_pol_tot(-1:1,-1:1), error_pol_tot(-1:1,-1:1)
  real :: all, error_all, al, error_al, apv, error_apv
  real :: chi2_sigma, error_sigma, stantot
  real :: alfas
  integer :: ndimensions, lam3, lam4, i, j, k, today(3), now(3)
  integer, parameter :: log = 27
  double precision :: start_time, finish_time

  call cpu_time(start_time)
  call read_config
  call modify_config
  call rootinit(ntuple_file)

  i = len(log_file)
  do while(log_file(i:i) == '')
    i = i - 1
  end do
  open(unit = log, file = log_file(1:i), status = "replace", action = "write")

  print*, "Ntuple: ", ntuple_file
  print*, "Log: ", log_file

  s = collider_energy*collider_energy

  ! pdfs are intrinsically linked to the value of lamda_qcd; alpha_qcd
  if (structure_function == 1) lambdaqcd4 = 0.326d0
  if (structure_function == 2) lambdaqcd4 = 0.326d0
  if (structure_function == 3) lambdaqcd4 = 0.326d0
  if (structure_function == 4) lambdaqcd4 = 0.215d0
  if (structure_function == 5) lambdaqcd4 = 0.300d0
  if (structure_function == 6) lambdaqcd4 = 0.300d0
  if (structure_function == 7) lambdaqcd4 = 0.300d0
  if (structure_function == 8) lambdaqcd4 = 0.229d0
  if (structure_function == 9) lambdaqcd4 = 0.383d0

  ! initialise cteq grids.
  if (structure_function <= 4) call setctq6(structure_function)

  ! use appropriately evolved alphas.
  if (structure_function <= 2) then
    nloops = 2
  else
    nloops = 1
  end if

  ! initialise madgraph - masses and coupling constants of particles
  call initialise_model

  if (final_state <= 0) then
    m3 = fmass(ffinal)
    m4 = fmass(ffinal)
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
    if (final_state <= 0) then
      ! integrates on:
      ! x(3) = (x1 - tau)/(1 - tau),
      ! x(2) = (ecm - ecm_min)/(ecm_max - ecm_min),
      ! x(1) = cos(theta3_cm).
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
      ! x(14) = (ecm - ecm_min)/(ecm_max - ecm_min),
      ! x(13) = (xx356 - xx356min)/(xx356max - xx356min),
      !   where xx356 = arctg((m356**2 - m3**2)/m3/gamt),
      !   or x(13) = (m356 - m356min)/(m356max - m356min),
      !   where m356min = m3 + m5 + m6, m356max = ecm_max - m4 - m7 - m8
      ! x(12) = (xx478 - xx478min)/(xx478max - xx478min),
      !   where xx478 = arctg((m478**2 - m3**2)/m3/gamt),
      !   or x(12) = (m478 - m478min)/(m478max - m478min),
      !   where m478min = m4 + m7 + m8, m478max = ecm_max - m356
      ! x(11) = (xx56 - xx56min)/(xx56max - xx56min),
      !   where xx56 = arctg((m56**2 - rm_w**2)/rm_w/gamw),
      !   or x(11) = (m56 - m56min)/(m56max - m56min),
      !   where m56min = m5 + m6, m56max = m356 - m3
      ! x(10) = (xx78 - xx78min)/(xx78max - xx78min),
      !   where xx78 = arctg((m78**2 - rm_w**2)/rm_w/gamw),
      !   or x(10) = (m78 - m78min)/(m78max - m78min),
      !   where m78min = m7 + m8, m78max = m478 - m4
      ! x(9) = cos(theta_cm_356) = -cos(theta_cm_478),
      ! x(8) = cos(theta56_cm_356),
      ! x(7) = cos(theta78_cm_478),
      ! x(6) = cos(theta5_cm_56),
      ! x(5) = cos(theta7_cm_78),
      ! x(4) = fi56_cm_356,
      ! x(3) = fi78_cm_478,
      ! x(2) = fi5_cm_56,
      ! x(1) = fi8_cm_78.
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
    ! x(1) = (ecm - ecm_min)/(ecm_max - ecm_min)
    ndimensions = 2
    do i = 2, 1, -1
      xl(i) = 0.d0
      xu(i) = 1.d0
    end do
  end if


  if (initial_state == 0) then
    if (final_state == -1) then
      write(log,*) "Process:p p -> l+ l-"
    else if (final_state == 0) then
      write(log,*) 'Process:p p -> t t~'
    else if (final_state == 1) then
      write(log,*) 'Process:p p -> t t~ -> b b W+ W- -> b b l+ l- nu nu~'
    end if
  else if (initial_state == 1) then
    if (final_state == -1) then
      write(log,*) 'pp~ -> l+l-'
    else if (final_state == 0) then
      write(log,*) 'pp~ -> tt~'
    else if (final_state == 1) then
      write(log,*) 'Process:p p~ -> t t~ -> b b W+ W- -> b b l+ l- nu nu~'
    end if
  end if
  if (structure_function == 1) write(log,*) 'PDFs:CTEQ6m'
  if (structure_function == 2) write(log,*) 'PDFs:CTEQ6d'
  if (structure_function == 3) write(log,*) 'PDFs:CTEQ6l'
  if (structure_function == 4) write(log,*) 'PDFs:CTEQ6l1'
  if (structure_function == 5) write(log,*) 'PDFs:MRS99 (cor01)'
  if (structure_function == 6) write(log,*) 'PDFs:MRS99 (cor02)'
  if (structure_function == 7) write(log,*) 'PDFs:MRS99 (cor03)'
  if (structure_function == 8) write(log,*) 'PDFs:MRS99 (cor04)'
  if (structure_function == 9) write(log,*) 'PDFs:MRS99 (cor05)'
  if ((final_state >= 1) .and. (use_nwa == 1)) write(log,*) 'NWA:ON'
  if ((final_state >= 1) .and. (use_nwa == 0)) write(log,*) 'NWA:OFF'
  write(log,*) 'Model:', model_name
  if (include_qcd == 1) write(log,*) 'QCD:ON '
  if (include_qcd == 0) write(log,*) 'QCD:OFF'
  if (include_qfd == 1) write(log,*) 'QFD:ON '
  if (include_qfd == 0) write(log,*) 'QFD:OFF'
  if (include_bsm == 1) write(log,*) 'BSM:ON '
  if (include_bsm == 0) write(log,*) 'BSM:OFF'
  if (include_gg == 1) write(log,*) 'gg:ON '
  if (include_gg == 0) write(log,*) 'gg:OFF'
  if (include_qq == 1) write(log,*) 'qq:ON '
  if (include_qq == 0) write(log,*) 'qq:OFF'
  if (interference == 0) write(log,*) "Interference:none"
  if (interference == 1) write(log,*) "Interference:(gamma + Z) + (Z')"
  if (interference == 2) write(log,*) "Interference:(gamma + Z + Z')"
  if (interference == 3) write(log,*) "Interference:(gamma + Z + Z') - (gamma) - (Z)"
  if (interference == 4) write(log,*) "Interference:(gamma + Z + Z') - (gamma) - (Z) - (Z')"
  if (symmetrise_x1x2 == 1) write(log,*) 'Symmetrising integration:x1<->x2'
  if (symmetrise_costheta_t == 1) write(log,*) 'symmetrising integration:costhetat'
  if (symmetrise_costheta_5 == 1) write(log,*) 'symmetrising integration:costheta5'
  if (symmetrise_costheta_7 == 1) write(log,*) 'symmetrising integration:costheta7'
  if (use_rambo == 1) write(log,*) 'RAMBO:ON'
  if (map_phase_space == 0) write(log,*) "Phase space mapping:ON"
  write(log,*) 'Seed:', seed
  write(log,*) 'Collider energy:', collider_energy
  if (ecm_low > 0) write(log,*) "E_CM low         :", ecm_low
  if (ecm_up > 0) write(log,*) "E_CM up             ", ecm_up
  write(log,*) 'Loops:', nloops
  write(log,*) 'alpha_s(m_Z):', alfas(rm_z, lambdaqcd4, nloops)
  write(log,*) 'lambdaQCD^4:', lambdaqcd4
  write(log,*) 'm_b:', fmass(12)
  write(log,*) 'Gamma_b:', fwidth(12)
  write(log,*) 'm_t:', fmass(11)
  write(log,*) 'Gamma_t:', fwidth(11)
  write(log,*) 'm_Z:', rm_z
  write(log,*) 'Gamma_Z:', gamma_z
  write(log,*) 'm_W:', rm_w
  write(log,*) 'Gamma_W:', gamma_w
  write(log,*) 'm_h:', rm_h
  write(log,*) 'Gamma_h:', gamma_h
  if (include_bsm == 1) then
    do i = 1, 5
      if (mass_zp(i) > 0) then
        write(log,*) "Z' no.:", i
        write(log,*) "m_Z':", mass_zp(i), "[GeV]"
        write(log,*) "Gamma_Z':", gamzp(i), "[GeV]"
        write(log,*) "Gamma_Z'/m_Z':", gamzp(i)/mass_zp(i)
        write(log,*) "gLu:", gZpu(1,i)
        write(log,*) "gRu:", gZpu(2,i)
        write(log,*) "gLd:", gZpd(1,i)
        write(log,*) "gRd:", gZpd(2,i)
        write(log,*) "gLl:", gZpl(1,i)
        write(log,*) "gRl:", gZpl(2,i)
        write(log,*) "gLn:", gZpn(1,i)
        write(log,*) "gRn:", gZpn(2,i)
        write(log,*) "gLt:", gZpt(1,i)
        write(log,*) "gRt:", gZpt(2,i)
        write(log,*) "gLb:", gZpb(1,i)
        write(log,*) "gRb:", gZpb(2,i)
        write(log,*) "gLl3:", gZpl3(1,i)
        write(log,*) "gRl3:", gZpl3(2,i)
        write(log,*) "gLn3:", gZpn3(1,i)
        write(log,*) "gRn3:", gZpn3(2,i)
      end if
    end do
    print*, '------'
  end if

  ! reset counter
  npoints = 0

  ! reset
  if (final_state <= 0) then
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
    end do
  end if

  ! integrate
  call debug("Starting integration... ")
  it = 0
  ! if nprn < 0 no print-out
  nprn = 0
  call vegas(ndimensions, dsigma, sigma, error_sigma, chi2_sigma)

  call debug("...complete.")

  call debug("Calculating factor to re-weight for different iterations.")
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
  call debug("...complete.")

  do i = 1, it
  	write(log,*) "Iteration weighting:", i, ":", cnorm(i)
  end do

  if (sigma == 0.d0) then
    write(log,*) "Error = sigma = 0"
    stop
  else
    write(log,*) "Cross section:", sigma, ":", error_sigma, ":[pb]"
  end if

  if (final_state == 0) then
    call debug("Collating polar cross sections...")
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
    call debug("...complete.")

    call debug("Calculating polar asymmetries...")
    all = (sigma_pol_tot(+1, +1) - sigma_pol_tot(+1, -1) - sigma_pol_tot(-1, +1) + sigma_pol_tot(-1, -1))/sigma
    error_all = (sigma_pol_tot(+1, +1) + sigma_pol_tot(+1, -1) + sigma_pol_tot(-1, +1) + sigma_pol_tot(-1, -1))/4.d0*all

    al = (sigma_pol_tot(-1, -1) - sigma_pol_tot(+1, -1) + sigma_pol_tot(-1, +1) - sigma_pol_tot(+1, +1))/sigma
    error_al = (sigma_pol_tot(-1, -1) + sigma_pol_tot(+1, -1) + sigma_pol_tot(-1, +1) + sigma_pol_tot(+1, +1))/4.d0*al

    apv = (sigma_pol_tot(-1, -1) - sigma_pol_tot(+1, +1))/sigma/2.d0
    error_apv = (sigma_pol_tot(-1, -1) + sigma_pol_tot(+1, +1))/2.d0*apv

    write(log,*) "ALL:", all, ":", error_all
    write(log,*) "AL:", al, ":", error_al
    write(log,*) "APV:", apv, ":", error_apv
    call debug("...complete.")

  end if
  write(log,*) "VEGAS points:", npoints
  call rootclose
  write(log,*) 'Author:Declan Millar'
  call idate(today)     ! today(1):day, (2):month, (3):year
  call itime(now)       ! now(1):hour, (2):minute, (3):second
  write(log,*) 'Date:', today(3), today(2), today(1)
  write(log,*) 'Time:', now(1), now(2), now(3)
  call cpu_time(finish_time)
  write(log,*), "Run-time:", finish_time - start_time
  close(log)
  stop
end program zprime
