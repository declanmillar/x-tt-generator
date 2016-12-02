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
    use mpi90
    use exceptions
    use tao_random_numbers
    use vampi

    implicit none

    real(kind = default) :: sigma, error, chi2
    integer :: i, j, k

    ! time keeping
    integer, dimension(8) :: now
    integer :: start_time, end_time, runtime
    integer :: integrate_start, integrate_end
    integer :: event_start, event_end

    ! lhe
    integer :: idbm(2), pdfg(2), pdfs(2)
    real(kind=default) :: ebm(2)

    ! VAMP
    real(kind = default) event
    integer :: ndimensions
    real(kind = default), dimension(:), allocatable :: x
    real(kind = default), dimension(:,:), allocatable :: domain
    type(exception) :: exc
    type(tao_random_state) :: rng
    type(vamp_grid) :: grid
    type(vamp_data_t) :: data
    integer, dimension(2,2) :: calls
    type(vamp_history), dimension(10) :: history
    real(kind=default), parameter :: acceptable = 4
    integer :: failures

    real(kind = default) :: all, error_all, al, error_al, apv, error_apv

    ! ---

    call date_and_time(values = now)
    write(log, *) 'date: ', now(1), now(2), now(3)
    write(log, *) 'time: ', now(5), now(6), now(7)

    call read_config
    call modify_config
    call print_config
    call initialise_pdfs
    call initialise_model
    call initialise_masses
    call initialise_s
    call set_energy_limits


    if (use_rambo) then
        ndimensions = 2
    else
        if (final_state <= 0) then
            ndimensions = 3
        else if (final_state > 0) then
            ndimensions = 15
        end if
    end if

    print*, "integration: allocating domain with", ndimensions, "dimensions ..."
    allocate(x(ndimensions))
    allocate(domain(2, ndimensions))

    print*, "integration: setting VAMP limits ..."
    if (use_rambo) then
        do i = 2, 1, -1
            domain(1, i) = 0.d0
            domain(2, i) = 1.d0
        end do
    else
        if (final_state <= 0) then
            do i = 3, 2, -1
                domain(1, i) = 0.d0
                domain(2, i) = 1.d0
            end do
            do i = 1, 1
                domain(1, i) = -1.d0
                domain(2, i) = 1.d0
            end do

        else if (final_state > 0) then
            do i = 15, 10, -1
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
    print*, "integration: initialising MPI ..."
    call mpi90_init

    print*, "integration: integrating using VAMP ..."
    call tao_random_create (rng, 0)
    call system_clock (integrate_start) 
    call tao_random_seed (rng, integrate_start)
    call vamp_create_history (history)

    print*, "integration: creating VAMP grid with", ncall / 10, "calls ..."
    call clear_exception(exc)
    call vamp_create_grid(grid, domain, num_calls = ncall / 10, exc = exc) 
    call handle_exception(exc)

    print*, "integration: initial sampling of VAMP grid with", itmx + 1, "iterations ..."
    call clear_exception(exc)
    call vamp_sample_grid(rng, grid, dsigma, itmx + 1, exc = exc, history = history)
    call handle_exception(exc)

    print*, "integration: discarding preliminary integral with", ncall, "calls ..."
    call clear_exception(exc)
    call vamp_discard_integral(grid, num_calls = ncall, exc = exc)
    call handle_exception(exc)

    print*, "integration: full sampling of VAMP grid with ", itmx - 1, "iterations ..."
    call clear_exception(exc)
    call vamp_sample_grid(rng, grid, dsigma, itmx - 1, sigma, error, chi2, exc = exc, history = history(itmx + 2:))
    call handle_exception(exc)

    ! print*, "integration: refining grid ..."
    ! call clear_exception(exc)
    ! call vamp_sample_grid0(rng, grid, dsigma, NO_DATA, exc = exc, history = history(2 * itmx + 1:))
    ! call handle_exception(exc)

    print*, "integration: printing history ..."
    call vamp_print_history(history, "history")
    call vamp_delete_history(history)

    print *, "integration: integral = ", sigma, "+/-", error, " (chi^2 = ", chi2, ")"
    call system_clock(integrate_end)
    print *, "integration: time = ", integrate_end - integrate_start

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

    if (final_state <= 0) then
        print*, "initialisation: reset polarised arrays ..."
        do j = -1, +1, 2
            do k = -1, +1, 2
                sigma_pol(j, k) = 0.d0
                error_pol(j, k) = 0.d0
            end do
        end do
    end if

    symmetrise = .false.

    print*, "vamp: generating", nevents, " events ..."
    call system_clock(event_start)
    if (.not. batch) call set_total(nevents)
    do i = 1, nevents
        call clear_exception(exc)
        call vamp_next_event(x, rng, grid, dsigma, exc = exc)
        call handle_exception(exc)
        record_events = .true.
        event = dsigma(x, NO_DATA)
        record_events = .false. 
        if (.not. batch) call progress_percentage(i)
    end do
    call system_clock(event_end)
    print *, "integration: time = ", event_end - event_start

    if (final_state <= 0) then
        print *, "finalisation: calculating asymmetries for polarized final state"
        do j = -1, 1, 2
            do k = -1, 1, 2
                sigma_pol(j, k) = sigma_pol(j, k) * sigma
                error_pol(j, k) = sigma_pol(j, k) * error
            end do
        end do
    
        all = (sigma_pol(1, 1) - sigma_pol(1, -1) - sigma_pol(-1, 1) + sigma_pol(-1, -1)) / sigma
        error_all = (sigma_pol(1, 1) + sigma_pol(1, -1) + sigma_pol(-1, 1) + sigma_pol(-1, -1)) / 4.d0 * all
    
        al = (sigma_pol(-1, -1) - sigma_pol(1, -1) + sigma_pol(-1, 1) - sigma_pol(1, 1)) / sigma
        error_al = (sigma_pol(-1, -1) + sigma_pol(1, -1) + sigma_pol(-1, 1) + sigma_pol(1, 1)) / 4.d0 * al
    
        apv = (sigma_pol(-1, -1) - sigma_pol(1, 1)) / sigma / 2.d0
        error_apv = (sigma_pol(-1, -1) + sigma_pol(1, 1)) / 2.d0 * apv
    
        write(log, *) "ALL:", all, ":", error_all
        write(log, *) "AL:", al, ":", error_al
        write(log, *) "APV:", apv, ":", error_apv
    end if


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

    call system_clock(end_time)
    runtime = end_time - start_time
    write(log, *) "runtime: ", runtime
    write(log, *) 'author: Declan Millar'
    close(log)
    print*, "runtime: ", runtime, "seconds"
    stop
end program generator
