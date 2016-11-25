program generator

    ! calculates the cross section, asymmetries, and generates events for
    !   p p -> f f~
    !   p p -> t t~ -> b b~ W+ W- -> b b~ l+ l- vl vl~
    !
    ! uses
    !   HELAS for probability amplitudes
    !   VAMP for Monte Carlo integration and event generation
    !   CTEQ6, CT10 for PDFs
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

    real(kind = default) :: sigma, error, chi2
    integer :: i, j, k

    ! time keeping
    integer, dimension(8) :: now
    real(kind = default) :: start_time, finish_time, runtime
    real(kind = default) :: event_start, event_finish, event_time

    ! lhe
    integer :: idbm(2), pdfg(2), pdfs(2)
    real(kind=default) :: ebm(2)

    ! VAMP
    real(kind = default) event
    integer :: ndimensions
    real(kind = default), dimension(15) :: x
    real(kind = default), dimension(2,15) :: domain
    type(exception) :: exc
    type(tao_random_state) :: rng
    type(vamp_grid) :: grid
    type(vamp_data_t) :: data
    integer, dimension(2,2) :: calls 


    ! real(kind = default) :: sigma_pol_tot(-1:1, -1:1), error_pol_tot(-1:1, -1:1), stantot
    ! real(kind = default) :: all, error_all, al, error_al, apv, error_apv

    ! ---

    print*, precision(1.0)

    call cpu_time(start_time)
    call read_config
    call modify_config
    call print_config
    call initialise_pdfs
    call initialise_model
    call initialise_masses
    call initialise_s
    call set_energy_limits


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
            do i = 15, 10, -1
            ! do i = 15, 14, -1
                domain(1, i) = 0.d0
                domain(2, i) = 1.d0
            ! end do
            ! do i = 13, 10, -1
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

    record_events = .false.
    print*, "integration: integrating using VAMP ..."
    ! call tao_random_create (rng, seed = seed)
    call tao_random_create (rng, seed = 0)

    print*, "integration: creating VAMP grid with", ncall / 10, "calls ..."
    call clear_exception(exc)
    call vamp_create_grid(grid, domain, num_calls = ncall / 10, exc = exc) 
    call handle_exception(exc)

    ! print*, "integration: warm up VAMP grid with", itmx + 1, "iterations ..."
    ! call clear_exception(exc)
    ! call vamp_warmup_grid(rng, grid, dsigma, NO_DATA, itmx + 1, exc = exc)
    ! call handle_exception(exc)

    print*, "integration: initial sampling of VAMP grid with", itmx + 1, "iterations ..."
    call clear_exception(exc)
    call vamp_sample_grid(rng, grid, dsigma, NO_DATA, itmx + 1, exc = exc)
    call handle_exception(exc)

    print*, "integration: discarding preliminary integral with", ncall, "calls ..."
    call clear_exception(exc)
    call vamp_discard_integral(grid, num_calls = ncall, exc = exc)
    call handle_exception(exc)

    print*, "integration: full sampling of VAMP grid with ", itmx - 1, "iterations ..."
    call clear_exception(exc)
    call vamp_sample_grid(rng, grid, dsigma, NO_DATA, itmx - 1, sigma, error, chi2, exc = exc)
    call handle_exception(exc)
    print *, "integration: integral = ", sigma, "+/-", error, " (chi^2 = ", chi2, ")"  

    print*, "integration: sampling grid 0 ..."
    call clear_exception(exc)
    call vamp_sample_grid0(rng, grid, dsigma, NO_DATA, exc = exc)
    call handle_exception(exc)

    if (sigma == 0.d0) then
        print*, "integration: integration: integral = 0. Stopping."
        stop
    else
        write(log,*) "cross section: ", sigma
        write(log,*) "uncertainty: ", error
        write(log,*) "chi^2: ", chi2
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
        call lhe_beam(idbm(1), idbm(2), ebm(1), ebm(2), pdfg(1), pdfg(2), pdfs(1), pdfs(2), 2)
        call lhe_process(sigma, error, 1.d0, 81)
    end if

    symmetrise = .false.

    print*, "vamp: generating", nevents, " events ..."
    if (.not. batch) call set_total(nevents)
    do i = 1, nevents
        ! call clear_exception (exc)
        call vamp_next_event(x, rng, grid, dsigma, NO_DATA, exc = exc)
        ! call handle_exception(exc)
        ! record_events = .true.
        event = dsigma(x, NO_DATA)
        ! record_events = .false. 
        if (.not. batch) call progress_percentage(i)
    end do

    ! reset
    ! if (final_state <= 0) then
    !     do i = 1, 20
    !         resl(i) = 0.d0
    !         standdevl(i) = 0.d0
    !         cnorm(i) = 0.d0
    !         do j = -1, +1, 2
    !             do k = -1, +1, 2
    !                 sigma_pol(j, k, i) = 0.d0
    !                 error_pol(j, k, i) = 0.d0
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
    !     do j = -1, 1, 2
    !         do k = -1, 1, 2
    !           sigma_pol_tot(j, k) = 0.d0
    !           error_pol_tot(j, k) = 0.d0
    !           do i = 1,  it
    !               sigma_pol(j, k, i) = sigma_pol(j, k, i) * sigma / cnorm(i)
    !               error_pol(j, k, i) = sigma_pol(j, k, i) * error / cnorm(i)
    !               sigma_pol_tot(j, k) = sigma_pol_tot(j, k) + sigma_pol(j, k, i)
    !               error_pol_tot(j, k) = sigma_pol_tot(j, k) + error_pol(j, k, i)
    !           end do
    !           if (sigma_pol_tot(j, k) == 0) then
    !               error_pol_tot(j, k) = 0
    !           else
    !               error_pol_tot(j, k) = error_pol_tot(j, k) / sigma_pol_tot(j, k)
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

    call date_and_time(VALUES = now)
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
    write(log, *) 'date: ', now(1), now(2), now(3)
    write(log, *) 'time: ', now(5), now(6), now(7)
    write(log, *) "runtime: ", runtime
    close(log)
    print*, "runtime: ", runtime, "seconds"
    stop
end program generator
