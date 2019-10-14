program generator

    ! calculates the cross section, asymmetries, and generates events for
    !   p p -> f f~
    !   p p -> t t~ -> b b~ W+ W- -> b b~ l+ l- vl vl~
    !
    ! uses
    !   HELAS for probability amplitudes
    !   VAMP for Monte Carlo integration and event generation
    !   MRS99, CTEQ6, CT14 PDFs
    !
    ! authors
    !   Declan Millar
    !   Stefano Moretti

    use kinds
    use runtime
    use progress
    use configuration
    use modelling
    use scattering
    use lhef
    ! use mpi90
    use exceptions
    use tao_random_numbers
    use vamp

    implicit none

    integer :: i, j
    real(kind = default) :: input_xsec(2), xsec, xerr, chi2, weight, all, error_all, al, error_al, apv, error_apv

    ! time keeping
    integer :: ticks, now(8)
    real(kind = default) :: time0, time1

    ! lhe
    integer :: idbm(2), pdfg(2), pdfs(2)
    real(kind = default) :: ebm(2)

    ! VAMP
    integer :: ndimensions, nweighted, calls(2, 3)
    real(kind = default), allocatable :: x(:), domain(:, :), weights(:)
    type(exception) :: exc
    type(tao_random_state) :: rng
    type(vamp_grid) :: grid
    type(vamp_grids) :: grids
    type(vamp_history), allocatable :: history(:), histories(:, :)

    call cpu_time(time0)
    call date_and_time(values = now)
    write(*, "(a8, i4, a1, i2.2, a1, i2.2, a1, i2.2, a1, i2.2, a1, i2.2)") &
        "time = ", now(1), "-", now(2), "-", now(3), " ", now(5), ":", now(6), ":", now(7)

    call read_config
    call modify_config
    call print_config
    call initialise_pdfs
    call initialise_model
    call initialise_masses
    call initialise_s
    call set_energy_limits

    if (verbose) print*, "generator: creating random number seed ..."
    call system_clock(ticks)
    call tao_random_create(rng, 0)
    call tao_random_seed(rng, ticks)

    if (verbose) print*, "generator: setting dimensions ..."
    if (use_rambo) then
        ndimensions = 2
    else
        if (final_state < 1) then
            ndimensions = 3
        else
            ndimensions = 16
        end if
    end if

    if (verbose) print*, "generator: allocating x with", ndimensions, "dimensions ..."
    allocate(x(ndimensions))

    xsec = 1.0
    xerr = 0.0

    if (new_grid) then
        if (verbose) print*, "generator: integrating using VAMP ..."
        record_events = .false.

        if (verbose) print*, "generator: allocating domain with", ndimensions, "dimensions ..."
        allocate(weights(3))
        allocate(domain(2, ndimensions))
        allocate(history(2 * itmx - 1))
        allocate(histories(3 * itmx, size(weights)))

        if (verbose) print*, "generator: setting VAMP limits ..."
        if (use_rambo) then
            do i = 2, 1, -1
                domain(1, i) = 0.d0
                domain(2, i) = 1.d0
            end do
        else
            if (final_state < 1) then
                do i = 3, 2, -1
                    domain(1, i) = 0.d0
                    domain(2, i) = 1.d0
                end do
                do i = 1, 1
                    domain(1, i) = -1.d0
                    domain(2, i) = 1.d0
                end do
            else
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
                    domain(2, i) = twopi
                end do
            end if
        end if

        calls(:, 1) = (/ itmx, ncall / 10 /)
        calls(:, 2) = (/ itmx, ncall / 10 /)
        calls(:, 3) = (/ itmx, ncall /)

        call cpu_time(time1)

        call vamp_create_history(history)

        ! preliminary sampling
        call vamp_create_grid(grid, domain, num_calls = calls(2, 1))
        call clear_exception(exc)
        call vamp_sample_grid(rng, grid, event, no_data, calls(1, 1), &
                                integral = xsec, std_dev = xerr, avg_chi2 = chi2, &
                                exc = exc, history = history)
        call handle_exception(exc)
        call vamp_print_history(history, "prelim")

        if (xsec <= 0) stop

        ! full sampling
        call vamp_discard_integral(grid, num_calls = calls(2, 3))
        call clear_exception(exc)
        call vamp_sample_grid(rng, grid, event, no_data, calls(1, 3) - 1, &
                                integral = xsec, std_dev = xerr, avg_chi2 = chi2, &
                                exc = exc, history = history(calls(1, 1) + 1:))
        call handle_exception(exc)
        call vamp_print_history(history(calls(1, 1) + 1:), "full")
        call vamp_delete_history(history)
        call write_cross_section(xsec, xerr)

        ! final refinement
        call clear_exception(exc)
        call vamp_sample_grid0(rng, grid, event, no_data, exc = exc)
        call handle_exception(exc)
        call vamp_write_grid(grid, grid_file)
    else ! read grid file
        print*, "input VAMP grid = ", trim(grid_file)
        call vamp_read_grid(grid, grid_file)
        if (unweighted) then
            input_xsec = read_cross_section()
            xsec = input_xsec(1)
            xerr = input_xsec(2)
        end if
    end if

    print *, "xsec = ", xsec, "+/-", xerr
    if (xsec <= 0) stop

    if (nevents > 0) then
        if (verbose) print*, "process: calculating beam info ..."

        idbm(1) = 2212
        if (ppbar == 0) then
            idbm(2) = 2212
        else
            idbm(2) = -2212
        end if
        do i = 1, 2
            ebm(i) = sqrts / 2
            pdfg(i) = pdf_group
            pdfs(i) = pdf_set
        end do

        if (.not. unweighted) then
            xsec = 0.0
            xerr = 0.0
            sigma_pol = 0.d0
            error_pol = 0.d0
            nweighted = nevents

            print*, "weighted events: ", nweighted
            if (.not. batch) call set_total(nweighted)
            do i = 1,  nweighted
                if (.not. unweighted) record_events = .true.
                call clear_exception (exc)
                call vamp_next_event_single(x, rng, grid, event, no_data)
                call handle_exception (exc)
                ! if (.not. unweighted) call rootaddevent(weight) // TODO might need to write weighted events to LHEF in future
                xsec = xsec + weight
                xerr = xerr + weight * weight
                if (.not. batch) call progress_bar(i)
            end do

            xsec = xsec / nweighted
            xerr = xerr / nweighted / nweighted

            xerr = sqrt(xerr)

            print *, "xsec(e+e-)", xsec, "+/-", xerr

            if (final_state < 1) then
                if (verbose) print *, "finalisation: calculating asymmetries for polarized final state"
                do i = -1, 1, 2
                    do j = -1, 1, 2
                        sigma_pol(i, j) = sigma_pol(i, j) * xsec
                        error_pol(i, j) = sigma_pol(i, j) * xerr
                    end do
                end do

                all = (sigma_pol(1, 1) - sigma_pol(1, -1) - sigma_pol(-1, 1) + sigma_pol(-1, -1)) / xsec
                error_all = (sigma_pol(1, 1) + sigma_pol(1, -1) + sigma_pol(-1, 1) + sigma_pol(-1, -1)) / 4.d0 * all

                al = (sigma_pol(-1, -1) - sigma_pol(1, -1) + sigma_pol(-1, 1) - sigma_pol(1, 1)) / xsec
                error_al = (sigma_pol(-1, -1) + sigma_pol(1, -1) + sigma_pol(-1, 1) + sigma_pol(1, 1)) / 4.d0 * al

                apv = (sigma_pol(-1, -1) - sigma_pol(1, 1)) / xsec / 2.d0
                error_apv = (sigma_pol(-1, -1) + sigma_pol(1, 1)) / 2.d0 * apv

                print*, "ALL = ", all, "+/-", error_all
                print*, "AL = ", al, "+/-", error_al
                print*, "APV = ", apv, "+/-", error_apv
            end if
        end if

        print*, "lhef = ", trim(lhe_file)
        if (verbose) print*, "initiating lhef file ..."
        call lhe_open(lhe_file)
        call lhe_beam(idbm(1), idbm(2), ebm(1), ebm(2), pdfg(1), pdfg(2), pdfs(1), pdfs(2), idw)
        call lhe_process(xsec, xerr, 1.d0, final_state)

        if (unweighted) then
            print*, "unweighted events = ", nevents
            call cpu_time(time1)
            if (.not. batch) call set_total(nevents)
            do i = 1, nevents
                call clear_exception(exc)
                call vamp_next_event_single(x, rng, grid, event)
                call handle_exception(exc)
                record_events = .true.
                weight = event(x, no_data)
                record_events = .false.
                if (.not. batch) call progress_bar(i)
            end do
        end if

        call lhe_close
    end if

    call vamp_delete_grid(grid)

    call cpu_time(time1)
    call print_runtime(time1 - time0)
    stop
end program generator
