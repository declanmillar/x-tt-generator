program zprime

    ! calculates the cross section, asymmetries, and generates events for
    !   p p -> f f~
    !   p p -> t t~ -> b b~ W+ W- -> b b~ l+ l- vl vl~
    !
    ! requirements
    !   HELAS subroutines
    !   VEGAS Monte Carlo integration
    !   CT10, CTEQ6 and MRS99 PDF subroutines
    !   RootTuple for filling n-tuples
    !
    ! authors
    !   Declan Millar <declan.millar@cern.ch>
    !   Stefano Moretti <s.moretti@soton.ac.uk>

    use kinds
    use configuration
    use modelling
    use scattering
    use lhef
    use exceptions
    use tao_random_numbers
    use vamp

    implicit none

    real(kind=default) :: sigma_pol_tot(-1:1, -1:1), error_pol_tot(-1:1, -1:1)
    real(kind=default) :: all, error_all, al, error_al, apv, error_apv
    real(kind=default) :: chi2_sigma, error_sigma, stantot
    real(kind=default) :: alfas
    integer :: ndimensions, lam3, lam4, i, j, today(3), now(3)
    double precision :: start_time, finish_time, runtime
    character(40) tablefile
    integer :: idbm1, idbm2, pdfg(2), pdfs(2)
    real(kind=default) :: ebm(2)

    logical :: use_vamp
    real(kind=default), dimension(15) :: x
    type(exception) :: exc
    type(tao_random_state) :: rng
    type(vamp_grid) :: grid
    real(kind=default), dimension(2,15) :: domain
    integer, dimension(2,2) :: calls

    print*, "starting program zprime ..."

    call cpu_time(start_time)
    call read_config
    call modify_config
    call print_config

    print*, "log:     ", trim(log_file)
    open(unit = log, file = log_file, status = "replace", action = "write")

    if (ntuple_out) then
        print*, "n-tuple: ", trim(ntuple_file)
        print*, "initiating n-tuple ..."
        call rootinit(ntuple_file)
    end if
    if (lhef_out) then
        print*, "lhe:    ", trim(lhe_file)
        print*, "opening lhe file ..."
        call lhe_open(lhe_file)
    end if

    s = sqrts * sqrts

    ! pdfs are intrinsically linked to the value of lamda_qcd; alpha_qcd
    print*, "setting lambdaqcd based on PDF set ..."
    if (ipdf ==  1) lambdaqcd4 = 0.326d0
    if (ipdf ==  2) lambdaqcd4 = 0.326d0
    if (ipdf ==  3) lambdaqcd4 = 0.326d0
    if (ipdf ==  4) lambdaqcd4 = 0.215d0
    if (ipdf ==  5) lambdaqcd4 = 0.300d0
    if (ipdf ==  6) lambdaqcd4 = 0.300d0
    if (ipdf ==  7) lambdaqcd4 = 0.300d0
    if (ipdf ==  8) lambdaqcd4 = 0.229d0
    if (ipdf ==  9) lambdaqcd4 = 0.383d0
    if (ipdf == 10) lambdaqcd4 = 0.326d0
    if (ipdf == 11) lambdaqcd4 = 0.215d0

    ! initialise cteq grids
    if (ipdf <= 4) call setctq6(ipdf)
    if (ipdf == 10) tablefile = "ct14ln.pds"
    if (ipdf == 11) tablefile = "ct14ll.pds"
    if (ipdf > 9) call setct14(tablefile)

    print*, "using appropriately evolved alphas ..."
    if (ipdf <= 2) then
        nloops = 2
    else if (ipdf == 10) then
        nloops = 2
    else
        nloops = 1
    end if

    write(log, *) 'loops: ', nloops
    write(log, *) 'lambda_QCD^4: ', lambdaqcd4
    write(log, *) 'alpha_s(m_Z): ', alfas(zmass, lambdaqcd4, nloops)

    print*, "initialising masses and coupling constants ..."
    call initialise_model

    print*,"setting external particle masses ..."
    if (final_state <= 0) then
        m3 = fmass(ffinal)
        m4 = fmass(ffinal)
        m5 = 0.d0
        m6 = 0.d0
        m7 = 0.d0
        m8 = 0.d0
    else if (final_state == 1) then
        m3 = fmass(12)
        m4 = fmass(12)
        m5 = 0.d0
        m6 = 0.d0
        m7 = 0.d0
        m8 = 0.d0
    end if

    ! use_vamp = .true.
    ! if (use_vamp) then
    print*, "setting VAMP integration limits ..."
    if (use_rambo) then
        ndimensions = 2

        do i = 2, 1, -1
            domain(1, i) = 0.d0
            domain(2, i) = 1.d0
        end do
    else
        if (final_state <= 0) then
            ndimensions = 3

            do i = 3, 2, -1
                domain(1, i) = 0.d0
                domain(2, i) = 1.d0
            end do
            do i = 1, 1
                domain(1, i) = -1.d0
                domain(2, i) = 1.d0
            end do

        else if (final_state > 0) then
            ndimensions = 15

            do i = 15, 14, -1
                domain(1, i) = 0.d0
                domain(2, i) = 1.d0
            end do
            do i = 13, 10, -1
                domain(1, i) = 0.d0
                domain(2, i) = 1.d0
            end do
            do i = 9, 5, -1
                domain(1, i) = -1.d0
                domain(2, i) = 1.d0
            end do
            do i = 4, 1, -1
                domain(1, i) = 0.d0
                domain(2, i) = 2.d0 * pi
            end do
        end if
    end if
    ! else
    !     print*, "setting VEGAS integration limits ..."
    !     if (use_rambo) then
    !         ! integrates on
    !         !     x(2) = (x1 - tau) / (1 - tau),
    !         !     x(1) = (ecm - ecm_min) / (ecm_max - ecm_min)
    !         ndimensions = 2
    !         do i = 2, 1, -1
    !             xl(i) = 0.d0
    !             xu(i) = 1.d0
    !         end do
    !     else
    !         if (final_state <= 0) then
    !             ! integrates on
    !             !     x(3) = (x1 - tau) / (1 - tau)
    !             !     x(2) = (ecm - ecm_min) / (ecm_max - ecm_min)
    !             !     x(1) = cos(theta3_cm)
    !             ndimensions = 3

    !             ! limits:
    !             do i = 3, 2, -1
    !                 xl(i) = 0.d0
    !                 xu(i) = 1.d0
    !             end do
    !             do i = 1, 1
    !                 xl(i) = -1.d0
    !                 xu(i) = 1.d0
    !             end do

    !         else if (final_state > 0) then
    !             ! integrates on
    !             !     x(15) = (x1 - tau) / (1 - tau)
    !             !     x(14) = (ecm - ecm_min) / (ecm_max - ecm_min),
    !             !     x(13) = (xx356 - xx356min) / (xx356max - xx356min),
    !             !       where xx356 = arctg((m356**2 - m3**2) / m3 / gamt),
    !             !       or x(13) = (m356 - m356min) / (m356max - m356min),
    !             !       where m356min = m3 + m5 + m6, m356max = ecm_max - m4 - m7 - m8
    !             !     x(12) = (xx478 - xx478min) / (xx478max - xx478min)
    !             !       where xx478 = arctg((m478**2 - m3**2) / m3 / gamt),
    !             !       or x(12) = (m478 - m478min) / (m478max - m478min),
    !             !       where m478min = m4 + m7 + m8, m478max = ecm_max - m356
    !             !     x(11) = (xx56 - xx56min) / (xx56max - xx56min),
    !             !       where xx56 = arctg((m56**2 - rm_w**2) / rm_w / gamw),
    !             !       or x(11) = (m56 - m56min) / (m56max - m56min),
    !             !       where m56min = m5 + m6, m56max = m356 - m3
    !             !     x(10) = (xx78 - xx78min) / (xx78max - xx78min),
    !             !       where xx78 = arctg((m78**2 - rm_w**2) / rm_w / gamw),
    !             !       or x(10) = (m78 - m78min) / (m78max - m78min),
    !             !       where m78min = m7 + m8, m78max = m478 - m4
    !             !     x(9) = cos(theta_cm_356) <--> -cos(theta_cm_478),
    !             !     x(8) = cos(theta56_cm_356),
    !             !     x(7) = cos(theta78_cm_478),
    !             !     x(6) = cos(theta5_cm_56),
    !             !     x(5) = cos(theta7_cm_78),
    !             !     x(4) = fi56_cm_356,
    !             !     x(3) = fi78_cm_478,
    !             !     x(2) = fi5_cm_56,
    !             !     x(1) = fi8_cm_78
    !             ndimensions = 15

    !             ! set integration limits:
    !             do i = 15, 14, -1
    !                 xl(i) = 0.d0
    !                 xu(i) = 1.d0
    !             end do
    !             do i = 13, 10, -1
    !                 xl(i) = 0.d0
    !                 xu(i) = 1.d0
    !             end do
    !             do i = 9, 5, -1
    !                 xl(i) = -1.d0
    !                 xu(i) = 1.d0
    !             end do
    !             do i = 4, 1, -1
    !                 xl(i) = 0.d0
    !                 xu(i) = 2.d0 * pi
    !             end do
    !         end if
    !     end if
    ! end if

    ! reset counter
    ! npoints = 0

    ! reset
    ! if (final_state <= 0) then
    !     do i = 1, 20
    !         resl(i) = 0.d0
    !         standdevl(i) = 0.d0
    !         cnorm(i) = 0.d0
    !         do lam3 = -1, +1, 2
    !             do lam4 = -1, +1, 2
    !                 sigma_pol(lam3, lam4, i) = 0.d0
    !                 error_pol(lam3, lam4, i) = 0.d0
    !             end do
    !         end do
    !     end do
    ! end if

    ! if (use_vamp) then
    event_mode = .false.
    print*, "integrating using VAMP ..."
    call tao_random_create (rng, seed = 0)

    print*, "creating VAMP grid..."
    call clear_exception (exc)
    call vamp_create_grid (grid, domain, num_calls = ncall / 10, exc = exc) 
    call handle_exception (exc)

    print*, "initial sampling of VAMP grid with ...", ncall / 10, " points"
    call clear_exception (exc)
    call vamp_sample_grid (rng, grid, dsigma, NO_DATA, itmx + 1, sigma, error_sigma, chi2_sigma, exc = exc)
    call handle_exception (exc)
    print *, "preliminary integral = ", sigma, "+/-", error_sigma, " (chi^2 = ", chi2_sigma, ")"

    print*, "discarding first integral ..."
    call clear_exception (exc)
    call vamp_discard_integral (grid, num_calls = ncall, exc = exc)
    call handle_exception (exc)

    print*, "full sampling of VAMP grid with ...", ncall, " points"
    call clear_exception (exc)
    call vamp_sample_grid (rng, grid, dsigma, NO_DATA, itmx - 1, sigma, error_sigma, chi2_sigma, exc = exc)
    call handle_exception (exc)
    print *, "integral = ", sigma, "+/-", error_sigma, " (chi^2 = ", chi2_sigma, ")"

    ! alternative method
    ! calls(:,1) = (/ 6, ncall / 10 /) 
    ! calls(:,2) = (/ 4, ncall /) 
    ! call clear_exception(exc)
    ! call vamp_integrate(rng, domain, dsigma, calls, sigma, error_sigma, chi2_sigma, exc = exc)
    ! call handle_exception(exc)
    ! else
    !     print*, "integrating using VEGAS ..."
    !     call vegas(ndimensions, dsigma, sigma, error_sigma, chi2_sigma)
    ! end if

    if (sigma == 0.d0) then
        print*, "ERROR: Cross section = 0. Stopping."
        stop
    end if

    event_mode = .true.
    print*, "generating...", ncall, " events"
    call clear_exception (exc)
    call vamp_next_event_single(x, rng, grid, dsigma, NO_DATA)!, NO_DATA, weight, channel, weights, grids, exc)
    call handle_exception (exc)


    ! print*, "calculating iteration weightings ..."
    ! stantot = 0.d0
    ! do i = 1, it
    !     stantot = stantot + 1.d0 / standdevl(i) / standdevl(i)
    ! end do
    ! do i = 1, it
    !     standdevl(i) = standdevl(i) * standdevl(i) * stantot
    ! end do
    ! do i = 1, it
    !     cnorm(i) = resl(i) * standdevl(i)
    ! end do

    ! do i = 1, it
    !     write(log, *) "iteration weighting:", i, ":", cnorm(i)
    ! end do

    ! if (final_state == 0) then
    !     do lam3 = -1, 1, 2
    !         do lam4 = -1, 1, 2
    !           sigma_pol_tot(lam3, lam4) = 0.d0
    !           error_pol_tot(lam3, lam4) = 0.d0
    !           do i = 1,  it
    !               sigma_pol(lam3, lam4, i) = sigma_pol(lam3, lam4, i) * sigma / cnorm(i)
    !               error_pol(lam3, lam4, i) = sigma_pol(lam3, lam4, i) * error_sigma / cnorm(i)
    !               sigma_pol_tot(lam3, lam4) = sigma_pol_tot(lam3, lam4) + sigma_pol(lam3, lam4, i)
    !               error_pol_tot(lam3, lam4) = sigma_pol_tot(lam3, lam4) + error_pol(lam3, lam4, i)
    !           end do
    !           if (sigma_pol_tot(lam3, lam4) == 0) then
    !               error_pol_tot(lam3, lam4) = 0
    !           else
    !               error_pol_tot(lam3, lam4) = error_pol_tot(lam3, lam4) / sigma_pol_tot(lam3, lam4)
    !           end if
    !         end do
    !     end do

    !     all = (sigma_pol_tot(+1, +1) - sigma_pol_tot(+1, -1) - sigma_pol_tot(-1, +1) + sigma_pol_tot(-1, -1)) / sigma
    !     error_all = (sigma_pol_tot(+1, +1) + sigma_pol_tot(+1, -1) + sigma_pol_tot(-1, +1) + sigma_pol_tot(-1, -1)) / 4.d0 * all

    !     al = (sigma_pol_tot(-1, -1) - sigma_pol_tot(+1, -1) + sigma_pol_tot(-1, +1) - sigma_pol_tot(+1, +1)) / sigma
    !     error_al = (sigma_pol_tot(-1, -1) + sigma_pol_tot(+1, -1) + sigma_pol_tot(-1, +1) + sigma_pol_tot(+1, +1)) / 4.d0 * al

    !     apv = (sigma_pol_tot(-1, -1) - sigma_pol_tot(+1, +1)) / sigma / 2.d0
    !     error_apv = (sigma_pol_tot(-1, -1) + sigma_pol_tot(+1, +1)) / 2.d0 * apv

    !     write(log, *) "ALL:", all, ":", error_all
    !     write(log, *) "AL:", al, ":", error_al
    !     write(log, *) "APV:", apv, ":", error_apv
    ! end if
    ! write(log, *) "VEGAS points: ", npoints

    if (ntuple_out) then
        print*, "closing n-tuple ..."
        call rootclose
    end if

    call idate(today)     ! today(1):day, (2):month, (3):year
    call itime(now)       ! now(1):hour, (2):minute, (3):second
    call cpu_time(finish_time)
    runtime = finish_time - start_time
    write(log, *) "runtime: ", runtime

    if (lhef_out) then
        print*, "calculating lhe beam info ..."
        idbm1 = 2212
        if (initial_state == 0) then
            idbm2 = 2212
        else
            idbm2 = -2212
        end if
        do i = 1, 2
            ebm(i) = sqrts / 2
            pdfg(i) = 1
            pdfs(i) = 1
        end do

        print*, "printing lhe footer ..."
        call lhe_footer()
        print*, "printing lhe header ..."
        call lhe_header(today(3), today(2), today(1), now(1), now(2), now(3), runtime)
        print*, "printing lhe beam info ..."
        call lhe_beam(idbm1, idbm2, ebm(1), ebm(2), pdfg(1), pdfg(2), pdfs(2), pdfs(2), 2)
        print*, "printing lhe process info ..."
        call lhe_process(sigma, error_sigma, 1.d0, 81)
        print*, "closing lhe file ..."
        call lhe_close
    end if

    write(log, *) 'author: Declan Millar'
    write(log, *) 'date: ', today(3), today(2), today(1)
    write(log, *) 'time: ', now(1), now(2), now(3)
    write(log, *) "runtime: ", runtime
    close(log)
    print*, "zprime program completed in", runtime, "seconds"
    stop
end program zprime
