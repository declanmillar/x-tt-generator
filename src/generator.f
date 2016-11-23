program generator

    ! calculates the cross section, asymmetries, and generates events for
    !   p p -> f f~
    !   p p -> t t~ -> b b~ W+ W- -> b b~ l+ l- vl vl~
    !
    ! uses
    !   HELAS for probability amplitudes
    !   VAMP for Monte Carlo integration and event generation
    !   CT10, CTEQ6 and MRS99 for PDFs
    !   RootTuple for filling n-tuples
    !
    ! authors
    !   Declan Millar <declan.millar@cern.ch>
    !   Stefano Moretti <s.moretti@soton.ac.uk>

    use kinds
    use progress
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
    integer :: lam3, lam4, i, j

    integer :: today(3), now(3)
    real(kind = default) :: start_time, finish_time, runtime
    real(kind = default) :: event_start, event_finish, event_time

    character(40) tablefile
    integer :: idbm(2), pdfg(2), pdfs(2)
    real(kind=default) :: ebm(2)

    ! VAMP
    integer :: ndimensions
    real(kind=default), dimension(15) :: x
    type(exception) :: exc
    type(tao_random_state) :: rng
    type(vamp_grid) :: grid
    real(kind=default), dimension(2,15) :: domain
    real(kind=default) event

    ! ---

    call cpu_time(start_time)
    call read_config
    call modify_config
    call print_config

    s = sqrts * sqrts

    ! pdfs are intrinsically linked to the value of lamda_qcd; alpha_qcd
    print*, "initialisation: setting lambdaqcd based on PDF set ..."
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

    print*, "initialisation: using appropriately evolved alphas ..."
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

    print*, "initialisation: masses and coupling constants ..."
    call initialise_model

    print*,"initialisation: setting external particle masses ..."
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

    print*, "integration: setting VAMP limits ..."
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

    ! cut = .false.
    record_events = .false.
    print*, "integration: integrating using VAMP ..."
    call tao_random_create (rng, seed = 0)

    print*, "integration: creating VAMP grid..."
    call vamp_create_grid (grid, domain, num_calls = ncall / 10, exc = exc) 

    print*, "integration: initial sampling of VAMP grid with ", ncall / 10, " points and ", itmx + 1, "iterations ..."
    call vamp_sample_grid (rng, grid, dsigma, NO_DATA, itmx + 1, sigma, error_sigma, chi2_sigma, exc = exc)
    print *, "integration: preliminary integral = ", sigma, "+/-", error_sigma, " (chi^2 = ", chi2_sigma, ")"

    print*, "integration: discarding preliminary integral ..."
    call vamp_discard_integral (grid, num_calls = ncall, exc = exc)

    print*, "integration: full sampling of VAMP grid with ...", ncall, " points and ", itmx - 1, "iterations ..."
    call vamp_sample_grid (rng, grid, dsigma, NO_DATA, itmx - 1, sigma, error_sigma, chi2_sigma, exc = exc)
    print *, "integration: integral = ", sigma, "+/-", error_sigma, " (chi^2 = ", chi2_sigma, ")"

    if (sigma == 0.d0) then
        print*, "integration: integration: integral = 0. Stopping."
        stop
    else
        write(log,*) "cross section: ", sigma
        write(log,*) "uncertainty: ", error_sigma
        write(log,*) "chi^2: ", chi2_sigma
    end if


    if (ntuple_out) then
        print*, "n-tuple: ", trim(ntuple_file)
        print*, "initiating n-tuple ..."
        call rootinit(ntuple_file)
    end if

    if (lhef_out) then
        print*, "lhe: calculating beam info ..."
        idbm(1) = 2212
        if (initial_state == 0) then
            idbm(2) = 2212
        else
            idbm(2) = -2212
        end if
        do i = 1, 2
            ebm(i) = sqrts / 2
            pdfg(i) = 1
            pdfs(i) = 1
        end do

        call lhe_open(lhe_file)
        call lhe_header()
        call lhe_beam(idbm(1), idbm(2), ebm(1), ebm(2), pdfg(1), pdfg(2), pdfs(2), pdfs(2), 2)
        call lhe_process(sigma, error_sigma, 1.d0, 81)
    end if

    ! print*, "vamp: warming up grid ..."
    ! call vamp_warmup_grid(rng, grid, dsigma, NO_DATA, itmx, exc = exc)

    symmetrise = .false.
    ! cut = .true.

    print*, "vamp: generating single event ..."
    call cpu_time(event_start)
    call vamp_next_event_single(x, rng, grid, dsigma, NO_DATA, exc = exc)
    call cpu_time(event_finish)
    event_time = event_finish - event_start
    print*, "event time: ", event_time, " [secs]"

    print*, "vamp: generating", nevents, " events ...  (estimate:", nevents * event_time, "[secs])"
    do i = 1, nevents
        call vamp_next_event_single(x, rng, grid, dsigma, NO_DATA, exc = exc)
        record_events = .true.
        event = dsigma(x, NO_DATA)
        record_events = .false. 
        call ProgressPercentage(i, nevents, 50)
    end do

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
    !
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
    !
    ! do i = 1, it
    !     write(log, *) "iteration weighting:", i, ":", cnorm(i)
    ! end do
    !
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
    !
    !     all = (sigma_pol_tot(+1, +1) - sigma_pol_tot(+1, -1) - sigma_pol_tot(-1, +1) + sigma_pol_tot(-1, -1)) / sigma
    !     error_all = (sigma_pol_tot(+1, +1) + sigma_pol_tot(+1, -1) + sigma_pol_tot(-1, +1) + sigma_pol_tot(-1, -1)) / 4.d0 * all
    !
    !     al = (sigma_pol_tot(-1, -1) - sigma_pol_tot(+1, -1) + sigma_pol_tot(-1, +1) - sigma_pol_tot(+1, +1)) / sigma
    !     error_al = (sigma_pol_tot(-1, -1) + sigma_pol_tot(+1, -1) + sigma_pol_tot(-1, +1) + sigma_pol_tot(+1, +1)) / 4.d0 * al
    !
    !     apv = (sigma_pol_tot(-1, -1) - sigma_pol_tot(+1, +1)) / sigma / 2.d0
    !     error_apv = (sigma_pol_tot(-1, -1) + sigma_pol_tot(+1, +1)) / 2.d0 * apv
    !
    !     write(log, *) "ALL:", all, ":", error_all
    !     write(log, *) "AL:", al, ":", error_al
    !     write(log, *) "APV:", apv, ":", error_apv
    ! end if
    ! write(log, *) "VEGAS points: ", npoints

    call idate(today)
    call itime(now)
    call cpu_time(finish_time)
    runtime = finish_time - start_time
    write(log, *) "runtime: ", runtime

    if (ntuple_out) then
        print*, "n-tuple: closing ..."
        call rootclose
    end if

    if (lhef_out) then
        print*, "lhe: printing footer ..."
        call lhe_footer()
        print*, "lhe: closing file ..."
        call lhe_close
    end if

    write(log, *) 'author: Declan Millar'
    write(log, *) 'date: ', today(3), today(2), today(1)
    write(log, *) 'time: ', now(1), now(2), now(3)
    write(log, *) "runtime: ", runtime
    close(log)
    print*, "runtime: ", runtime, "seconds"
    stop
end program generator
