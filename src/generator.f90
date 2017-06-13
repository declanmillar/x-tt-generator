program generator

    ! calculates the cross section, asymmetries, and generates events for
    !   p p -> f f~
    !   p p -> t t~ -> b b~ W+ W- -> b b~ l+ l- vl vl~
    !
    ! uses
    !   HELAS for probability amplitudes
    !   VAMP for Monte Carlo integration and event generation
    !   MRS99, CTEQ6, CT14 PDFs
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

    integer :: i, j
    real(kind = default) :: sigma, error, chi2, weight, all, error_all, al, error_al, apv, error_apv

    ! time keeping
    integer :: ticks, proc_id, now(8)
    real(kind = default) :: time0, time1, time2, time3

    ! lhe
    integer :: idbm(2), pdfg(2), pdfs(2)
    real(kind = default) :: ebm(2)

    ! VAMP
    integer :: ndimensions, nweighted, calls(2, 3)
    real(kind = default), allocatable :: x(:), domain(:,:), weights(:)
    type(exception) :: exc
    type(tao_random_state) :: rng
    type(vamp_grid) :: grid
    type(vamp_grids) :: grids
    type(vamp_history), allocatable :: history(:), histories(:,:)

    real(kind = default) :: cross_section(2)

    call cpu_time(time0)
    call date_and_time(values = now)
    write(*,"(a8,     i4,      a1,  i2.2,   a1,  i2.2,   a1,  i2.2,   a1,  i2.2,   a1,  i2.2)") &
            "time = ", now(1), "-", now(2), "-", now(3), " ", now(5), ":", now(6), ":", now(7)

    call read_config
    call modify_config
    call print_config
    call initialise_pdfs
    call initialise_model
    call initialise_masses
    call initialise_s
    call set_energy_limits


    if (verbose) print*, "generator: initialising MPI ..."
    call mpi90_init()
    call mpi90_rank(proc_id)

    if (verbose) print*, "generator: creating random number seed ..."
    call system_clock(ticks)
    call tao_random_create(rng, 0)
    call tao_random_seed(rng, ticks + proc_id)

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

        calls(:, 1) = (/itmx, ncall / 10 /)
        calls(:, 2) = (/itmx, ncall / 10 /)
        calls(:, 3) = (/itmx, ncall /)

        call cpu_time(time1)

        if (multichannel) then

            call vamp_create_history(history)
            call vamp_create_history(histories)

            print*, "VAMP sampling points: ", calls(2, 1)

            weights = 1
            call vamp_create_grids(grids, domain, calls(2, 1), weights)

            print*, "generator: initial sampling of VAMP grid with", calls(1, 1), "iterations ..."
            call clear_exception(exc)
            call vamp_sample_grids(rng, grids, dsigma, calls(1, 1), sigma, error, chi2, exc = exc, &
                                   history = history, histories = histories)
            call clear_exception(exc)
            call vamp_print_history (history, "multi")
            call vamp_print_history (histories, "multi")

            print *, "generator: integral = ", sigma, "+/-", error, " [pb]"
            if (sigma <= 0) stop

            print*, "generator: discarding integral and re-sampling grid with ", calls(2, 2), "calls ..."
            call vamp_discard_integrals(grids, calls(2, 2))

            print*, "generator: refining weights for VAMP grid with ", calls(1, 2), "iterations ..."
            do i = 1, calls(1, 2)
                call clear_exception(exc)
                call vamp_sample_grids(rng, grids, dsigma, 1, sigma, error, chi2, exc = exc, &
                                       history = history(calls(1, 1) + i:), &
                                       histories = histories(calls(1, 1) + i:, :))
                call handle_exception(exc)
                call clear_exception(exc)
                call vamp_refine_weights(grids)
                call handle_exception(exc)
            end do

            print *, "generator: integral = ", sigma, "+/-", error, " (chi^2 = ", chi2, ")"
            if (sigma <= 0) stop

            print*, "generator: discarding integral and re-sampling grid with ", calls(2, 3), "calls ..."
            call vamp_discard_integrals(grids, calls(2, 3))

            print*, "generator: warming up grid with ", calls(1, 3), "iterations ..."
            call clear_exception(exc)
            ! call vamp_warmup_grids(rng, grids, dsigma, calls(1, 3), &
            !                        history = history(calls(1, 1) + calls(1, 2) + 1:), &
            !                        histories = histories(calls(1, 1) + calls(1, 2) + 1:, :))
            call vamp_sample_grids(rng, grids, dsigma, calls(1, 3), sigma, error, chi2, &
                                   history = history(calls(1, 1) + calls(1, 2):), &
                                   histories = histories(calls(1, 1) + calls(1, 2):, :))
            call clear_exception(exc)

            print *, "generator: integral = ", sigma, "+/-", error, " (chi^2 = ", chi2, ")"
            if (sigma <= 0) stop

            if (verbose) print*, "generator: printing history ..."
            call vamp_print_history(history, "history")
            call vamp_print_history(histories, "histories")
            call vamp_delete_history(history(itmx:))
            call vamp_delete_history(histories(itmx:,:))

            if (verbose) print*, "saving vamp grid to ", grid_file
            call vamp_write_grids(grids, grid_file)

        else
            call vamp_create_history(history)

            call cpu_time(time1)
            print*, "prelim points = ", calls(2, 1)
            call vamp_create_grid(grid, domain, num_calls = calls(2, 1))

            print*, "prelim iterations = ", calls(1, 1)
            call clear_exception(exc)
            call vamp_sample_grid(rng, grid, dsigma, calls(1, 1), sigma, error, chi2, exc = exc, history = history)
            call handle_exception(exc)
            print *, "sigma = ", sigma, "+/-", error
            if (sigma <= 0) stop

            call cpu_time(time2)
            print *, "prelim time = ", (time2 - time1) / 60, "[mins]"

            call cpu_time(time1)
            print*, "full points = ", calls(2, 3)
            call vamp_discard_integral(grid, num_calls = calls(2, 3))

            print*, "full iterations = ", calls(1, 3)
            call clear_exception(exc)
            call vamp_sample_grid(rng, grid, dsigma, calls(1, 3) - 1, sigma, error, chi2, exc = exc, &
                                  history = history(calls(1, 1) + 1:))
            call handle_exception(exc)
            call clear_exception(exc)
            call vamp_sample_grid0(rng, grid, dsigma, no_data, exc = exc)
            call handle_exception(exc)
            print *, "sigma = ", sigma, "+/-", error
            if (sigma <= 0) stop

            if (verbose) print*, "generator: printing history ..."
            call vamp_print_history(history, "history")
            call vamp_delete_history(history)

            if (verbose) print*, "generator: saving vamp grid to ", grid_file
            call vamp_write_grid(grid, grid_file)
        end if
        call cpu_time(time2)
        print *, "generator: time = ", (time2 - time1) / 60, "[mins]"
        print *, "generator: complete"

    else
        if (multichannel) then
            print*, "input VAMP grids = ", trim(grid_file)
            call vamp_read_grids(grids, grid_file)
        else
            print*, "input VAMP grid = ", trim(grid_file)
            call vamp_read_grid(grid, grid_file)
        end if
        cross_section = read_cross_section("/scratch/dam1g09/zprime/qq-tt-bbllvv_SM_13TeV_CT14LL.txt")
        sigma = cross_section(1)
        error = cross_section(2)
    end if

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

        if (ntuple_out) then
            print*, "n-tuple = ", trim(ntuple_file)
            if (verbose) print*, "initiating n-tuple ..."
            call rootinit(ntuple_file)
        end if

        sigma = 0.0
        error = 0.0
        sigma_pol = 0.d0
        error_pol = 0.d0

        if (.not. unweighted) then
        !     nweighted = ncall
        ! else
            nweighted = nevents
        ! end if

          print*, "weighted events: ", nweighted
          if (.not. batch) call set_total(nweighted)
          do i = 1,  nweighted
              if (.not. unweighted) record_events = .true.
              call clear_exception (exc)
              if (multichannel) then
                 call vamp_next_event(x, rng, grids, dsigma, phi, weight, exc = exc)
              else
                 call vamp_next_event(x, rng, grid, dsigma, weight, exc = exc)
              end if
              call handle_exception (exc)
              if (.not. unweighted) call rootaddevent(weight)
              sigma = sigma + weight
              error = error + weight * weight
              if (.not. batch) call progress_bar(i)
          end do

          print *, "SIGMA(e+e-) = ", sigma, "+/-", sqrt(error)

          sigma = sigma / nweighted
          error = error / nweighted / nweighted

          print *, "sigma(e+e-) = ", sigma, "+/-", sqrt(error)

          ! dilepton full!!!
          sigma = sigma
          error = error

          print *, "sigma(l+l-) = ", sigma, "+/-", sqrt(error)

          if (final_state < 1) then
              if (verbose) print *, "finalisation: calculating asymmetries for polarized final state"
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

              print*, "ALL = ", all, "+/-", error_all
              print*, "AL = ", al, "+/-", error_al
              print*, "APV = ", apv, "+/-", error_apv
          end if
        end if

        if (ntuple_out) then
            call rootaddprocessdouble(idbm(1), "idbm1")
            call rootaddprocessdouble(idbm(2), "idbm2")
            call rootaddprocessdouble(ebm(1), "ebm1")
            call rootaddprocessdouble(ebm(2), "ebm2")
            call rootaddprocessdouble(pdfg(1), "pdfg1")
            call rootaddprocessdouble(pdfg(2), "pdfg2")
            call rootaddprocessdouble(pdfs(1), "pdfs1")
            call rootaddprocessdouble(pdfs(2), "pdfs2")
            call rootaddprocessdouble(sigma, "cross_section")
            call rootaddprocessdouble(sqrt(error), "cross_section_uncertainty")
        end if

        if (lhef_out) then
            print*, "lhef = ", trim(lhe_file)
            if (verbose) print*, "initiating lhef file ..."
            call lhe_open(lhe_file)
            call lhe_beam(idbm(1), idbm(2), ebm(1), ebm(2), pdfg(1), pdfg(2), pdfs(1), pdfs(2), idw)
            call lhe_process(sigma, sqrt(error), 1.d0, 9999)
        end if

        if (unweighted) then
            print*, "unweighted events = ", nevents
            call cpu_time(time1)
            if (.not. batch) call set_total(nevents)
            do i = 1, nevents
                call clear_exception(exc)
                if (multichannel) then
                    call vamp_next_event(x, rng, grids, dsigma, phi, exc = exc)
                else
                    call vamp_next_event(x, rng, grid, dsigma, exc = exc)
                end if
                call handle_exception(exc)
                record_events = .true.
                weight = dsigma(x, no_data)
                if (ntuple_out) call rootaddevent(1.d0)
                record_events = .false.
                if (.not. batch) call progress_bar(i)
            end do
            call cpu_time(time2)
            print *, "event generation time = ", (time2 - time1) / 60, "[mins]"
        end if

        if (ntuple_out) then
            if (verbose) print*, "n-tuple: closing ..."
            call rootclose
        end if

        if (lhef_out) call lhe_close

        if (verbose) print*, "generator: updating vamp grid at ", grid_file
        if (multichannel) then
            call vamp_write_grids(grids, grid_file)
        else
            call vamp_write_grid(grid, grid_file)
        end if
    else
        if (verbose) print*, "skipping event generation"
    end if

    if (.not. multichannel) then
        if (verbose) print *, "deleting grid ..."
        call vamp_delete_grid(grid)
    end if

    call cpu_time(time3)
    print*, "runtime = ", (time3 - time0) / 60, "[mins]"
    stop
end program generator
