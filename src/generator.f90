program generator

    ! calculates the cross section, asymmetries, and generates events for
    !   p p -> f f~
    !   p p -> t t~ -> b b~ W+ W- -> b b~ l+ l- vl vl~
    !
    ! uses
    !   HELAS for probability amplitudes
    !   VAMP for Monte Carlo integration and event generation
    !   MRS99, CTEQ6, CT10 for PDFs
    !   RootTuple for filling n-tuples
    !
    ! authors
    !   Declan Millar <declan.millar@cern.ch>
    !   Stefano Moretti 

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
    real(kind = default) :: all, error_all, al, error_al, apv, error_apv
    integer :: i, j

    ! time keeping
    integer, dimension(8) :: now
    real(kind = default) :: start_time, end_time, runtime
    real(kind = default) :: integrate_start, integrate_end
    real(kind = default) :: event_start, event_end
    integer :: ticks, num_proc, proc_id

    ! lhe
    integer :: idbm(2), pdfg(2), pdfs(2)
    real(kind=default) :: ebm(2)

    ! VAMP
    real(kind = default) event
    integer :: ndimensions
    real(kind = default), dimension(:), allocatable :: x
    real(kind = default), dimension(:, :), allocatable :: domain
    type(exception) :: exc
    type(tao_random_state) :: rng
    type(vamp_grid) :: grid
    type(vamp_data_t) :: data
    integer, dimension(2, 2) :: calls
    type(vamp_history), dimension(:), allocatable :: history

    call cpu_time(start_time) 
    call date_and_time(values = now)
    print*, 'date: ', now(1), now(2), now(3)
    print*, 'time: ', now(5), now(6), now(7)

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
            ndimensions = 16
        end if
    end if

    print*, "integration: allocating domain with", ndimensions, "dimensions ..."
    allocate(x(ndimensions))
    allocate(domain(2, ndimensions))
    allocate(history(2 * itmx))

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
            do i = 16, 10, -1
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
    call mpi90_init()
    call mpi90_size(num_proc)
    call mpi90_rank(proc_id)

    print*, "integration: integrating using VAMP ..."
    call system_clock(ticks)
    call tao_random_create(rng, 0)
    call tao_random_seed(rng, ticks + proc_id)
    call vamp_create_history(history)

    calls(:, 1) = (/6, ncall / 10 /)
    calls(:, 2) = (/4, ncall /)

    call cpu_time(integrate_start) 

    print*, "integration: creating VAMP grid with", calls(2, 1), "calls ..."
    call vamp_create_grid(grid, domain, num_calls = calls(2, 1)) 

    print*, "integration: initial sampling of VAMP grid with", calls(1, 1), "iterations ..."
    call clear_exception(exc)
    call vamp_sample_grid(rng, grid, dsigma, calls(1, 1), history = history, exc = exc)
    call handle_exception(exc)

    print*, "integration: discarding preliminary integral with", calls(2, 2), "calls ..."
    call vamp_discard_integral(grid, num_calls = calls(2, 2))

    ! print*, "integration: warm up VAMP grid with ", calls(1, 2), "iterations ..."
    ! call clear_exception(exc)
    ! call vamp_warmup_grid(rng, grid, dsigma, calls(1, 2), history = history(calls(1, 1) + 1:), exc = exc)
    ! call handle_exception(exc)

    print*, "integration: full sampling of VAMP grid with ", calls(1, 2), "iterations ..."
    call clear_exception(exc)
    call vamp_sample_grid(rng, grid, dsigma, calls(1, 2), sigma, error, chi2, history = history(calls(1, 1) + 1:), exc = exc)
    call handle_exception(exc)

    print*, "integration: refining grid ..."
    call clear_exception(exc)
    call vamp_sample_grid0(rng, grid, dsigma, no_data, exc = exc)
    call handle_exception(exc)

    print *, "integration: integral = ", sigma, "+/-", error, " (chi^2 = ", chi2, ")"
    call cpu_time(integrate_end)
    print *, "integration: time = ", (integrate_end - integrate_start) / 60, "[mins]"
    if (sigma <= 0) stop

    print*, "integration: printing history ..."
    call vamp_print_history(history, "history")
    call vamp_delete_history(history)

    print*, "process: calculating beam info ..."
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

    if (ntuple_out) then
        print*, "n-tuple: ", trim(ntuple_file)
        print*, "initiating n-tuple ..."
        call rootinit(ntuple_file)
        call rootaddprocessdouble(idbm(1), "idbm1")
        call rootaddprocessdouble(idbm(2), "idbm2")
        call rootaddprocessdouble(ebm(1), "ebm1")
        call rootaddprocessdouble(ebm(2), "ebm2")
        call rootaddprocessdouble(pdfg(1), "pdfg1")
        call rootaddprocessdouble(pdfg(2), "pdfg2")
        call rootaddprocessdouble(pdfs(1), "pdfs1")
        call rootaddprocessdouble(pdfs(2), "pdfs2")
        call rootaddprocessdouble(sigma, "cross_section")
        call rootaddprocessdouble(error, "cross_section_uncertainty")
        call rootaddprocessdouble(chi2, "vamp_chi2")
    end if

    if (lhef_out) then
        call lhe_open(lhe_file)
        call lhe_header()
        call lhe_beam(idbm(1), idbm(2), ebm(1), ebm(2), pdfg(1), pdfg(2), pdfs(1), pdfs(2), 2)
        call lhe_process(sigma, error, 1.d0, 81)
    end if

    if (final_state <= 0) then
        print*, "initialisation: reset polarised arrays ..."
        do i = -1, +1, 2
            do j = -1, +1, 2
                sigma_pol(i, j) = 0.d0
                error_pol(i, j) = 0.d0
            end do
        end do
    end if

    ! symmetrise = .false.

    print*, "vamp: generating", nevents, " events ..."
    call cpu_time(event_start)
    if (.not. batch) call set_total(nevents)
    do i = 1, nevents
        call clear_exception(exc)
        call vamp_next_event(x, rng, grid, dsigma, exc = exc)
        call handle_exception(exc)
        record_events = .true.
        ! symmetrise = .true.
        event = dsigma(x, no_data)
        print*, "event", i, "weight = ", event 
        record_events = .false.
        ! symmetrise = .false.
        if (.not. batch) call progress_percentage(i)
    end do
    call cpu_time(event_end)
    print *, "event generation: time = ", (event_end - event_start) / 60, "[mins]"

    call vamp_delete_grid(grid)

    if (final_state <= 0) then
        print *, "finalisation: calculating asymmetries for polarized final state"
        do i = -1, 1, 2
            do j = -1, 1, 2
                sigma_pol(i, j) = sigma_pol(i, j) * sigma
                error_pol(i, j) = sigma_pol(i, j) * error
            end do
        end do
    
        all = (sigma_pol(1, 1) - sigma_pol(1, -1) - sigma_pol(-1, 1) + sigma_pol(-1, -1)) / sigma
        error_all = (sigma_pol(1, 1) + sigma_pol(1, -1) + sigma_pol(-1, 1) + sigma_pol(-1, -1)) / 4.d0 * all
    
        al = (sigma_pol(-1, -1) - sigma_pol(1, -1) + sigma_pol(-1, 1) - sigma_pol(1, 1)) / sigma
        error_al = (sigma_pol(-1, -1) + sigma_pol(1, -1) + sigma_pol(-1, 1) + sigma_pol(1, 1)) / 4.d0 * al
    
        apv = (sigma_pol(-1, -1) - sigma_pol(1, 1)) / sigma / 2.d0
        error_apv = (sigma_pol(-1, -1) + sigma_pol(1, 1)) / 2.d0 * apv
    
        print*, "ALL:", all, ":", error_all
        print*, "AL:", al, ":", error_al
        print*, "APV:", apv, ":", error_apv
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

    call cpu_time(end_time)
    runtime = end_time - start_time
    print*, "runtime: ", runtime / 60, "[mins]"
    stop
end program generator
