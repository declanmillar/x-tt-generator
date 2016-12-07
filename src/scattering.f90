module scattering

    use kinds
    use configuration
    use modelling
    use tt_bbeevv

    implicit none

    public :: dsigma, initialise_masses, initialise_s, initialise_pdfs, set_energy_limits

    real(kind=default) :: sigma_pol(-1:1, -1:1), error_pol(-1:1, -1:1)
    real(kind=default), private :: m3, m4, m5, m6, m7, m8
    real(kind=default), private :: s, ecm_max, ecm_min
    real(kind=default), parameter, public :: unit_conv = 0.38937966d9 ! GeV^{-2} to nb (pb?)
    logical, public :: record_events
    real(kind=default), private :: qq

    ! temporary top mass and width
    real(kind=default) :: mt, gamt

contains


subroutine initialise_pdfs

    character(40) tablefile
    real(kind=default) :: alfas

    print*, "initialisation: setting lambda_QCD based on PDF set ..."
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
    print*, 'lambda_QCD^4: ', lambdaqcd4

    print*, "initialisation: opening pdf tables ..."
    if (ipdf <=  4) call setctq6(ipdf)
    if (ipdf == 10) tablefile = "ct14ln.pds"
    if (ipdf == 11) tablefile = "ct14ll.pds"
    if (ipdf >= 10) call setct14(tablefile)

    print*, "initialisation: using appropriately evolved alphas ..."
    if (ipdf <= 2 .or. ipdf == 10) then
        nloops = 2
    else
        nloops = 1
    end if
    print*, 'loops: ', nloops
    print*, 'alpha_s(m_Z): ', alfas(zmass, lambdaqcd4, nloops)

    ! scale for the pdfs
    if (final_state >= 0) then
        print*, "initialisation: scale = 2 * m_top"
        qq = 2.d0 * mt
    else
        qq = 0.d0
        print*, "initialisation: scale = ecm"
    end if

end subroutine initialise_pdfs

subroutine initialise_masses

    print*, "initialisation: setting external particle masses ..."
    if (final_state <= 0) then
        m3 = fmass(ffinal)
        m4 = fmass(ffinal)
    else if (final_state == 1) then
        m3 = fmass(12)
        m4 = fmass(12)
    end if
    m5 = 0.d0
    m6 = 0.d0
    m7 = 0.d0
    m8 = 0.d0

    ! store top parameters
    mt = fmass(ffinal)
    gamt = fwidth(ffinal)

end subroutine initialise_masses

subroutine set_energy_limits
    print*, "initialisation: setting energy limits ..."
    ecm_min = m3 + m4 + m5 + m6 + m7 + m8
    ecm_max = sqrts
    if (ecm_low > ecm_min) then
        ecm_min = ecm_low
    end if
    if (ecm_up > 0) then
        ecm_max = ecm_up
    end if
    print*, "initialisation: ecm_min = ", ecm_min
    print*, "initialisation: ecm_max = ", ecm_max
end subroutine set_energy_limits


subroutine initialise_s

    print*, "initialisation: setting s ..."

    s = sqrts * sqrts

    print*, "initialisation: s = ", s

end subroutine initialise_s



function dsigma(x, data, weights, channel, grids)

    ! computes the differential cross section for:
    !     p p -> f f~
    !     p p -> t t~ -> b b~ W+ W- -> b b~ l+ l- vl vl~

    ! requirements:
    !     HELAS subroutines
    !     VAMP Monte Carlo integration
    !     CT10, CTEQ6 and MRS99 PDF subroutines
    !     RootTuple for filling n-tuples

    ! authors:
    !     Declan Millar <declan.millar@cern.ch>
    !     Stefano Moretti

    use configuration
    use lhef
    use vamp, only: vamp_data_t, vamp_grid

    implicit none

    real(kind=default) :: dsigma
    real(kind=default), dimension(:), intent(in) :: x
    class(vamp_data_t), intent(in) :: data
    real(kind=default), dimension(:), intent(in), optional :: weights
    integer, intent(in), optional :: channel
    type(vamp_grid), dimension(:), intent(in), optional :: grids

    real(kind=default) ddsigma

    ! alphas
    real(kind=default) :: alfas, gs4, gs2, a_s

    ! square matrix elements
    real(kind=default) :: qfduu1, qfduu2, qfddd1, qfddd2, qcdqq, qcdgg
    real(kind=default) :: sgg_tt, sqq_tt, sqq_ff
    real(kind=default) :: sqq_tt_bbeevv_qcd, sqq_tt_bbeevv
    real(kind=default) :: sgg_bbemuvevm, sqq_bbemuvevm, suu_bbemuvevm, sdd_bbemuvevm
    real(kind=default) :: sgg_bbtatavtvt, sqq_bbtatavtvt, suu_bbtatavtvt, sdd_bbtatavtvt

    ! pdfs
    real(kind=default) :: ctq6pdf, ct14pdf
    real(kind=default) :: fx1(13), fx2(13), x1, x2, xx1, xx2, x1x2(2,2)
    real(kind=default) :: d1, u1, str1, chm1, btm1, glu1
    real(kind=default) :: d2, u2, str2, chm2, btm2, glu2
    real(kind=default) :: dbar1, ubar1, dbar2, ubar2
    real(kind=default) :: dsea1, usea1, dsea2, usea2
    integer :: ix, ixmax

    ! energies
    real(kind=default) :: shat, tau, ecm

    real(kind=default) :: qcm, pcm, qcm2
    real(kind=default) :: pq5, pq52, pq56, pq7, pq78
    real(kind=default) :: rl356, rl478, rl56, rl78, rpl356, rpl478
    real(kind=default) :: rps, rps356, rps478

    ! transfer invariant masses
    real(kind=default) :: q, q2, rq5, rq52, rq56, rq562, rq7, rq72, rq78, rq782

    ! boost
    real(kind=default) :: gamma, v

    ! angles
    real(kind=default) :: phi
    real(kind=default) :: st, st5, st56, st7, st78
    real(kind=default) :: ct, ct5, ct56, ct7, ct78
    real(kind=default) :: sf5, sf56, sf7, sf78
    real(kind=default) :: cf5, cf56, cf7, cf78

    ! rambo
    real(kind=default) :: rmass(100), prambo(4,100), wgtr
    integer :: jps

    ! for cuts
    real(kind=default) :: arg, eps, rpl, eta, pt

    ! arctan
    real(kind=default) :: at356, at356max, at356min, at56, at56max, at56min
    real(kind=default) :: at478, at478max, at478min, at78, at78max, at78min

    ! 4-momenta
    real(kind=default) :: p(4,8), p356(4), p478(4), q56(4), q78(4), p56(4), p78(4), q5(4), q7(4)
    real(kind=default) :: pcol(4,8), pcol56(4), pcol78(4), pcol356(4), pcol478(4)
    real(kind=default) :: p1(0:3), p2(0:3), p3(0:3), p4(0:3), p5(0:3), p6(0:3), p7(0:3), p8(0:3)

    ! invariant masses
    real(kind=default) :: m356, m356max, m356min, m478, m478max, m478min
    real(kind=default) :: m356_2, m478_2, m56, m56_2, m56max, m56min, m78, m78_2, m78max, m78min

    ! polarised
    real(kind=default) :: qcdpolqq(-1:1, -1:1), qcdpolgg(-1:1, -1:1)
    real(kind=default) :: qfdpoluu1(-1:1, -1:1), qfdpoldd1(-1:1, -1:1), qfdpolbb1(-1:1, -1:1)
    real(kind=default) :: qfdpoluu2(-1:1, -1:1), qfdpoldd2(-1:1, -1:1), qfdpolbb2(-1:1, -1:1)
    real(kind=default) :: weight_pol(-1:1, -1:1), dsigmapol(-1:1, -1:1)

    real(kind=default) :: random
    integer :: seed

    ! temporary iterators
    integer :: i, j

    ! ---

    if (verbose) print*, "dsigma: begin"

    call random_number(random)

    if (use_rambo) then
        ecm = x(1) * (ecm_max - ecm_min) + ecm_min
    else 
        ecm = x(2 + 12 * tops_decay) * (ecm_max - ecm_min) + ecm_min
    end if
    shat = ecm * ecm
    tau = shat / s

    if (qq == 0.d0) qq = ecm

    ! x1 and x2 of the partons
    if (use_rambo) then
        xx1 = x(2) * (1.d0 - tau) + tau
    else
        xx1 = x(3 + 12 * tops_decay) * (1.d0 - tau) + tau
    end if
    xx2 = tau / xx1
    x1x2(1, 1) = xx1
    x1x2(1, 2) = xx2
    x1x2(2, 1) = xx2
    x1x2(2, 2) = xx1

    ! symmetrise phase space with x1 <-> x2
    if (symmetrise) then
        ixmax = 2
    else
        ixmax = 1
    end if

    dsigma = 0.d0
    do ix = 1, ixmax
        ddsigma = 0.d0
        x1 = x1x2(ix, 1)
        x2 = x1x2(ix, 2)

        ! initialisation
        ddsigma = 0.d0
        do i = 1, 100
            rmass(i) = 0.d0
            do j = 1, 4
                prambo(j,i) = 0.d0
            end do
        end do
        do i = 1, 4
            do j = 1, 8
                p(i,j) = 0.d0
                pcol(i,j) = 0.d0
            end do
        end do

        if (verbose) print*, "constructing hadronic structure functions ..."
        if (ipdf <= 4) then
            if ((x1 <= 1.d-6) .or. (x1 >= 1.d0)) return
            if ((x2 <= 1.d-6) .or. (x2 >= 1.d0)) return
            if ((qq <= 1.3d0) .or. (qq >= 1.d4)) return

            ! *x for compatibility with MRS which return xf(x)
            u1 = x1 * ctq6pdf(1, x1, qq)
            d1 = x1 * ctq6pdf(2, x1, qq)
            ubar1 = x1 * ctq6pdf(-1, x1, qq)
            dbar1 = x1 * ctq6pdf(-2, x1, qq)
            str1 = x1 * ctq6pdf(3, x1, qq)
            chm1 = x1 * ctq6pdf(4, x1, qq)
            btm1 = x1 * ctq6pdf(5, x1, qq)
            glu1 = x1 * ctq6pdf(0, x1, qq)
            u2 = x2 * ctq6pdf(1, x2, qq)
            d2 = x2 * ctq6pdf(2, x2, qq)
            ubar2 = x2 * ctq6pdf(-1, x2, qq)
            dbar2 = x2 * ctq6pdf(-2, x2, qq)
            str2 = x2 * ctq6pdf(3, x2, qq)
            chm2 = x2 * ctq6pdf(4, x2, qq)
            btm2 = x2 * ctq6pdf(5, x2, qq)
            glu2 = x2 * ctq6pdf(0, x2, qq)

        else if (ipdf > 4 .and. ipdf < 10) then
            if ((x1 <= 1.d-5) .or. (x1 >= 1.d0)) return
            if ((x2 <= 1.d-5) .or. (x2 >= 1.d0)) return
            if ((qq * qq <= 1.25d0) .or. (qq * qq >= 1.d7)) return

            call mrs99(x1, qq, ipdf-4, u1, d1, usea1, dsea1, str1, chm1, btm1, glu1)
            call mrs99(x2, qq, ipdf-4, u2, d2, usea2, dsea2, str2, chm2, btm2, glu2)

            u1 = u1 + usea1
            d1 = d1 + dsea1
            u2 = u2 + usea2
            d2 = d2 + dsea2
            ubar1 = usea1
            dbar1 = dsea1
            ubar2 = usea2
            dbar2 = dsea2
        else if (ipdf > 9) then
            if ((x1 <= 1.d-9) .or. (x1 >= 1.d0)) then
                if (verbose) print*, "invalid x1"
                return
            end if
            if ((x2 <= 1.d-9) .or. (x2 >= 1.d0)) then
                if (verbose) print*, "invalid x2"
                return
            end if
            if ((qq <= 1.3d0) .or. (qq >= 1.d5)) then
                if (verbose) print*, "invalid qq"
                return
            end if

            ! *x for compatibility with MRS which return xf(x)
            u1 = x1 * ct14pdf(1, x1, qq)
            d1 = x1 * ct14pdf(2, x1, qq)
            ubar1 = x1 * ct14pdf(-1, x1, qq)
            dbar1 = x1 * ct14pdf(-2, x1, qq)
            str1 = x1 * ct14pdf(3, x1, qq)
            chm1 = x1 * ct14pdf(4, x1, qq)
            btm1 = x1 * ct14pdf(5, x1, qq)
            glu1 = x1 * ct14pdf(0, x1, qq)
            u2 = x2 * ct14pdf(1, x2, qq)
            d2 = x2 * ct14pdf(2, x2, qq)
            ubar2 = x2 * ct14pdf(-1, x2, qq)
            dbar2 = x2 * ct14pdf(-2, x2, qq)
            str2 = x2 * ct14pdf(3, x2, qq)
            chm2 = x2 * ct14pdf(4, x2, qq)
            btm2 = x2 * ct14pdf(5, x2, qq)
            glu2 = x2 * ct14pdf(0, x2, qq)
        end if

        ! initialise pdfs
        fx1(1) = d1
        fx1(2) = u1
        fx1(3) = str1
        fx1(4) = chm1
        fx1(5) = btm1
        fx1(6) = 0.d0
        fx1(7) = dbar1
        fx1(8) = ubar1
        fx1(9) = fx1(3)
        fx1(10) = fx1(4)
        fx1(11) = fx1(5)
        fx1(12) = fx1(6)
        fx1(13) = glu1
        fx2(1) = d2 * (1 - initial_state) + dbar2 * initial_state
        fx2(2) = u2 * (1 - initial_state) + ubar2 * initial_state
        fx2(3) = str2
        fx2(4) = chm2
        fx2(5) = btm2
        fx2(6) = 0.d0
        fx2(7) = d2 * initial_state + dbar2 * (1 - initial_state)
        fx2(8) = u2 * initial_state + ubar2 * (1 - initial_state)
        fx2(9) = fx2(3)
        fx2(10) = fx2(4)
        fx2(11) = fx2(5)
        fx2(12) = fx2(6)
        fx2(13) = glu2
        do i = 1, 13
            fx1(i) = fx1(i) / x1
            fx2(i) = fx2(i) / x2
        end do

        pcm = ecm / 2.d0
        p(4,1) = pcm
        p(3,1) = pcm
        p(2,1) = 0.d0
        p(1,1) = 0.d0
        p(4,2) = pcm
        p(3,2) = -pcm
        p(2,2) = 0.d0
        p(1,2) = 0.d0

        if (use_rambo) then
            call system_clock(seed) 
        end if

        if (final_state < 1) then
            if (use_rambo) then
                rmass(1) = m3
                rmass(2) = m4
                jps = 2
                call rambo(seed, jps, ecm, rmass, prambo, wgtr)
                do i = 3, jps + 2
                    do j = 1, 4
                        p(j, i) = prambo(j, i - 2)
                    end do
                end do
            else
                ! Calculate 2->2 final state momenta in the parton CoM frame manually

                phi = 2.d0 * pi * random
                ct = x(1)
                st = sqrt(1.d0 - ct * ct)

                ! magnitude of 3 momentum for products in general two body decay
                qcm2 = ((ecm * ecm - m3 * m3 - m4 * m4)**2 - (2.d0 * m3 * m4)**2) / (4.d0 * ecm * ecm)
                if (qcm2 < 0.d0) then
                    if (verbose) print*, "invalid qcm2"
                    return
                end if
                qcm = sqrt(qcm2)

                p(4,3) = sqrt(qcm2 + m3 * m3)
                p(3,3) = qcm * ct
                p(2,3) = qcm * st * cos(phi)
                p(1,3) = qcm * st * sin(phi)
                p(4,4) = sqrt(qcm2 + m4 * m4)
                p(3,4) = -qcm * ct
                p(2,4) = -qcm * st * cos(phi)
                p(1,4) = -qcm * st * sin(phi)
            end if

        else if (final_state > 0) then
            if (use_rambo) then
                rmass(1) = m3
                rmass(2) = m4
                rmass(3) = m5
                rmass(4) = m6
                rmass(5) = m7
                rmass(6) = m8
                jps = 6
                call rambo(seed, jps, ecm, rmass, prambo, wgtr)
                do i = 3, jps + 2
                    do j = 1, 4
                       p(j,i) = prambo(j, i - 2)
                    end do
                end do
            else
                phi = 2.d0 * pi * random

                m356min = m3 + m5 + m6
                m356max = ecm - m4 - m7 - m8
                if (map_phase_space) then
                    at356min = atan((m356min * m356min - mt * mt) / mt / gamt)
                    at356max = atan((m356max * m356max - mt * mt) / mt / gamt)
                    at356 = x(13) * (at356max - at356min) + at356min
                    rl356 = tan(at356) * mt * gamt
                    m356_2 = mt * mt + rl356
                    if (m356_2 < 0.d0) then
                        if (verbose) print*, "invalid"
                        return
                    end if
                    m356 = sqrt(m356_2)
                else
                    m356 = x(13) * (m356max - m356min) + m356min
                end if

                m478min = m4 + m7 + m8
                m478max = ecm - m356
                if (map_phase_space) then
                    at478min = atan((m478min * m478min - mt * mt) / mt / gamt)
                    at478max = atan((m478max * m478max - mt * mt) / mt / gamt)
                    at478 = x(12) * (at478max - at478min) + at478min
                    rl478 = tan(at478) * mt * gamt
                    m478_2 = mt * mt + rl478
                    if (m478_2 < 0.d0) then
                        if (verbose) print*, "invalid"
                        return
                    end if
                    m478 = sqrt(m478_2)
                else 
                    m478 = x(12) * (m478max - m478min) + m478min
                end if

                m56min = m5 + m6
                m56max = m356 - m3
                if (map_phase_space) then
                    at56min = atan((m56min * m56min - wmass * wmass) / wmass / wwidth)
                    at56max = atan((m56max * m56max - wmass * wmass) / wmass / wwidth)
                    at56 = x(11) * (at56max - at56min) + at56min
                    rl56 = tan(at56) * wmass * wwidth
                    m56_2 = wmass * wmass + rl56
                    if (m56_2 < 0.d0) then
                        if (verbose) print*, "invalid"
                        return
                    end if
                    m56 = sqrt(m56_2)
                else
                    m56 = x(11) * (m56max - m56min) + m56min
                end if

                m78min = m7 + m8
                m78max = m478 - m4
                if (map_phase_space) then
                    at78min = atan((m78min*m78min - wmass * wmass) / wmass / wwidth)
                    at78max = atan((m78max * m78max - wmass * wmass) / wmass / wwidth)
                    at78 = x(10) * (at78max - at78min) + at78min
                    rl78 = tan(at78) * wmass * wwidth
                    m78_2 = wmass * wmass + rl78
                    if (m78_2 < 0.d0) then
                        if (verbose) print*, "invalid"
                        return
                    end if
                    m78 = sqrt(m78_2)
                else
                    m78 = x(10)*(m78max - m78min) + m78min
                end if

                ! assign angles
                ct = x(9)
                st = sqrt(abs(1.d0 - ct * ct))
                ct56 = x(8)
                st56 = sqrt(1.d0 - ct56 * ct56)
                ct78 = x(7)
                st78 = sqrt(1.d0 - ct78 * ct78)
                ct5 = x(6)
                st5 = sqrt(1.d0 - ct5 * ct5)
                ct7 = x(5)
                st7 = sqrt(1.d0 - ct7 * ct7)
                cf56 = cos(x(4))
                sf56 = sin(x(4))
                cf78 = cos(x(3))
                sf78 = sin(x(3))
                cf5 = cos(x(2))
                sf5 = sin(x(2))
                cf7 = cos(x(1))
                sf7 = sin(x(1))

                ! two body decay of s-channel mediating boson
                q2 = ((ecm * ecm - m356 * m356 - m478 * m478)**2 - (2.d0 * m356 * m478)**2) / (4.d0 * ecm * ecm)
                if (q2 < 0.d0) then
                    if (verbose) print*, "invalid"
                    return
                end if
                q = sqrt(q2)

                p356(3) = q * ct
                p356(2) = q * st * cos(phi)
                p356(1) = q * st * sin(phi)
                p356(4) = sqrt(q2 + m356 * m356)

                do i = 1, 3
                    p478(i) = -p356(i)
                end do
                p478(4) = sqrt(q2 + m478 * m478)

                ! two body decay of the top
                rq562 = ((m356 * m356 - m3 * m3 - m56 * m56)**2 - (2.d0 * m3 * m56)**2) / (4.d0 * m356 * m356)
                if (rq562 < 0.d0) then
                    if (verbose) print*, "invalid"
                    return
                end if
                rq56 = sqrt(rq562)
                q56(3) = rq56 * st56 * cf56
                q56(2) = rq56 * st56 * sf56
                q56(1) = rq56 * ct56
                q56(4) = sqrt(rq562 + m56 * m56)
                pq56 = 0.d0
                do i = 1,3
                    pq56 = pq56 + p356(i) * q56(i)
                end do
                p56(4) = (p356(4) * q56(4) + pq56) / m356
                p(4,3) = p356(4) - p56(4)
                do i = 1,3
                    p56(i) = q56(i) + p356(i) * (p56(4) + q56(4)) / (p356(4) + m356)
                    p(i,3) = p356(i) - p56(i)
                end do

                ! two body decay of the anti-top
                rq782 = ((m478 * m478 - m4 * m4 - m78 * m78)**2 - (2.d0 * m4 * m78)**2) / (4.d0 * m478 * m478)
                if (rq782 < 0.d0) then
                    if (verbose) print*, "invalid"
                    return
                end if
                rq78 = sqrt(rq782)
                q78(3) = rq78 * st78 * cf78
                q78(2) = rq78 * st78 * sf78
                q78(1) = rq78 * ct78
                q78(4) = sqrt(rq782 + m78 * m78)
                pq78 = 0.d0
                do i = 1, 3
                  pq78 = pq78 + p478(i) * q78(i)
                end do
                p78(4) = (p478(4) * q78(4) + pq78) / m478
                p(4,4) = p478(4) - p78(4)
                do i = 1, 3
                    p78(i) = q78(i) + p478(i) * (p78(4) + q78(4)) / (p478(4) + m478)
                    p(i,4) = p478(i) - p78(i)
                end do

                ! two body decay of the W+
                rq52 = ((m56 * m56 - m5 * m5 - m6 * m6)**2 - (2.d0 * m5 * m6)**2) / (4.d0 * m56 * m56)
                if (rq52 < 0.d0) then
                    if (verbose) print*, "invalid"
                    return
                end if
                rq5 = sqrt(rq52)
                q5(3) = rq5 * st5 * cf5
                q5(2) = rq5 * st5 * sf5
                q5(1) = rq5 * ct5
                q5(4) = sqrt(rq52 + m5 * m5)
                pq5 = 0.d0
                do i = 1, 3
                    pq5 = pq5 + p56(i) * q5(i)
                end do
                p(4,5) = (p56(4) * q5(4) + pq5) / m56
                p(4,6) = p56(4) - p(4,5)
                do i = 1,3
                    p(i,5) = q5(i) + p56(i) * (p(4,5) + q5(4)) / (p56(4) + m56)
                    p(i,6) = p56(i) - p(i,5)
                end do

                ! two body decay of the W-
                rq72 = ((m78 * m78 - m7 * m7 - m8 * m8)**2 - (2.d0 * m7 * m8)**2) / (4.d0 * m78 * m78)
                if (rq72 < 0.d0) then
                    if (verbose) print*, "invalid"
                    return
                end if
                rq7 = sqrt(rq72)
                q7(3) = rq7 * st7 * cf7
                q7(2) = rq7 * st7 * sf7
                q7(1) = rq7 * ct7
                q7(4) = sqrt(rq72 + m7 * m7)
                pq7 = 0.d0
                do i = 1, 3
                    pq7 = pq7 + p78(i) * q7(i)
                end do
                p(4,7) = (p78(4) * q7(4) + pq7) / m78
                p(4,8) = p78(4) - p(4,7)
                do i = 1, 3
                    p(i,7) = q7(i) + p78(i) * (p(4,7) + q7(4)) / (p78(4) + m78)
                    p(i,8) = p78(i) - p(i,7)
                end do
            end if
        end if

        if (verbose) print*, "kinematics: boost initial and final state momenta to the collider frame ..."
        v = (x1 - x2) / (x1 + x2)
        gamma = (x1 + x2) / 2.d0 / sqrt(x1 * x2)
        do i = 1, nfinal
            pcol(4,i) = gamma * (p(4,i) + v * p(3,i))
            pcol(3,i) = gamma * (p(3,i) + v * p(4,i))
            pcol(2,i) = p(2,i)
            pcol(1,i) = p(1,i)
        end do

        if (cut) then
            if (verbose) print*, "kinematics: applying fiducial cuts ..."
            do i = 3, nfinal
                pt = sqrt(pcol(1,i) * pcol(1,i) + pcol(2,i) * pcol(2,i))
                if (pt < 25) return
                eta = atanh(pcol(3,i) / sqrt(pcol(1,i) * pcol(1,i) + pcol(2,i) * pcol(2,i) + pcol(3,i) * pcol(3,i)))
                if (abs(eta) > 2.5) return
            end do
        end if

        if (verbose) print*, "kinematics: copying parton CoM 4-momenta for MadGraph ..."
        p1(0) = p(4, 1)
        p2(0) = p(4, 2)
        p3(0) = p(4, 3)
        p4(0) = p(4, 4)
        p5(0) = p(4, 5)
        p6(0) = p(4, 6)
        p7(0) = p(4, 7)
        p8(0) = p(4, 8)
        do i = 1, 3
            p1(i) = p(i, 1)
            p2(i) = p(i, 2)
            p3(i) = p(i, 3)
            p4(i) = p(i, 4)
            p5(i) = p(i, 5)
            p6(i) = p(i, 6)
            p7(i) = p(i, 7)
            p8(i) = p(i, 8)
        end do

        if (verbose) then
            print*, "kinematics: checking ..."
            print*, "p1 = ", p1
            print*, "p2 = ", p2
            print*, "p3 = ", p3
            print*, "p4 = ", p4
            print*, "p5 = ", p5
            print*, "p6 = ", p6
            print*, "p7 = ", p7
            print*, "p8 = ", p8
            print*, "Delta E  = ", p1(0) + p2(0) - p3(0) - p4(0) - p5(0) - p6(0) - p7(0) - p8(0)
            print*, "Delta Px = ", p1(1) + p2(1) - p3(1) - p4(1) - p5(1) - p6(1) - p7(1) - p8(1)
            print*, "Delta Py = ", p1(2) + p2(2) - p3(2) - p4(2) - p5(2) - p6(2) - p7(2) - p8(2)
            print*, "Delta Pz = ", p1(3) + p2(3) - p3(3) - p4(3) - p5(3) - p6(3) - p7(3) - p8(3)
            print*, "m1 = ", sqrt(abs(p1(0) * p1(0) - p1(1) * p1(1) - p1(2) * p1(2) - p1(3) * p1(3)))
            print*, "m2 = ", sqrt(abs(p2(0) * p2(0) - p2(1) * p2(1) - p2(2) * p2(2) - p2(3) * p2(3)))
            print*, "m3 = ", sqrt(abs(p3(0) * p3(0) - p3(1) * p3(1) - p3(2) * p3(2) - p3(3) * p3(3)))
            print*, "m4 = ", sqrt(abs(p4(0) * p4(0) - p4(1) * p4(1) - p4(2) * p4(2) - p4(3) * p4(3)))
            print*, "m5 = ", sqrt(abs(p5(0) * p5(0) - p5(1) * p5(1) - p5(2) * p5(2) - p5(3) * p5(3)))
            print*, "m6 = ", sqrt(abs(p6(0) * p6(0) - p6(1) * p6(1) - p6(2) * p6(2) - p6(3) * p6(3)))
            print*, "m7 = ", sqrt(abs(p7(0) * p7(0) - p7(1) * p7(1) - p7(2) * p7(2) - p7(3) * p7(3)))
            print*, "m8 = ", sqrt(abs(p8(0) * p8(0) - p8(1) * p8(1) - p8(2) * p8(2) - p8(3) * p8(3)))
        end if

        if (.not. include_gg .and. .not. include_qq .and. .not. include_uu .and. .not. include_dd) then
            if (verbose) print*, "matrix elements: setting |M|^2 = 1 and skipping ..."
            if (ix == 1) then
                ddsigma = 0.5 / x1
            else if (ix == 2) then
                ddsigma = 0.5 / x2
            end if
            if (final_state <= 0) then
                do i = -1, 1, 2
                    do j = -1, 1, 2
                      if (ix == 1) then
                          dsigmapol(i,j) = 0.5 / x1 / (2 * ddsigma)
                      else if (ix == 2) then
                          dsigmapol(i,j) = 0.5 / x2 / (2 * ddsigma)
                      end if
                    end do
                end do
            end if
            go to 666
        end if

        if (verbose) print*,  "matrix elements: calculate QCD coupling ..."
        gs2 = 4.d0 * pi * alfas(qq, lambdaqcd4, nloops)
        gs4 = gs2 * gs2

        if (verbose) print*, "matrix elements: initialising ..."
        qcdqq = 0.d0
        qcdgg = 0.d0
        qfduu1 = 0.d0
        qfddd1 = 0.d0
        qfduu2 = 0.d0
        qfddd2 = 0.d0
        do i = -1, 1
            do j = -1, 1
                qcdpolqq(i, j) = 0.d0
                qcdpolgg(i, j) = 0.d0
                qfdpoluu1(i, j) = 0.d0
                qfdpoldd1(i, j) = 0.d0
                qfdpoluu2(i, j) = 0.d0
                qfdpoldd2(i, j) = 0.d0
                dsigmapol(i,j) = 0.d0
                weight_pol(i, j) = 0.d0
            end do
        end do

        if (final_state <= 0) then
            if (verbose) print*, "matrix elements: calculating 2 -> 2 ..."
            do i = -1, 1, 2
                do j = -1, 1, 2
                    if (include_gg) then
                        qcdpolgg(i,j) = sgg_tt(p1, p2, p3, p4, i, j) * gs4
                    end if
                    if (include_qq) then
                        qcdpolqq(i,j) = sqq_tt(3, p1, p2, p3, p4, i, j) * gs4
                    end if
                    if (include_uu) then
                        qfdpoluu1(i, j) = sqq_ff(3, ffinal, p1, p2, p3, p4, i, j)
                        qfdpoluu2(i, j) = sqq_ff(3, ffinal, p2, p1, p3, p4, i, j)
                    end if
                    if (include_dd) then
                        qfdpoldd1(i, j) = sqq_ff(4, ffinal, p1, p2, p3, p4, i, j)
                        qfdpoldd2(i, j) = sqq_ff(4, ffinal, p2, p1, p3, p4, i, j)
                    end if
                end do
            end do

        else if (final_state == 1) then
            if (verbose) print*, "matrix elements: calculating 2 -> 6 ..."
            if (include_gg) then
                qcdgg = sgg_tt_bbeevv(p1, p2, p3, p4, p5, p7, p6, p8)
            end if
            if (include_qq) then
                qcdqq = sqq_tt_bbeevv_qcd(3, p1, p2, p3, p4, p5, p7, p6, p8)
            end if
            if (include_uu) then
                qfduu1 = sqq_tt_bbeevv(3, 11, p1, p2, p3, p4, p5, p7, p6, p8)
                qfduu2 = sqq_tt_bbeevv(3, 11, p2, p1, p3, p4, p5, p7, p6, p8)
            end if
            if (include_dd) then
                qfddd1 = sqq_tt_bbeevv(4, 11, p1, p2, p3, p4, p5, p7, p6, p8)
                qfddd2 = sqq_tt_bbeevv(4, 11, p2, p1, p3, p4, p5, p7, p6, p8)
            end if
        else
            stop "error: invalid final state"
        end if

        if (verbose) print*, "matrix elements: multiply qcd |m|^2 by g_s^4 ..."
        ! madgraph gs is set to one due to its scale dependence
        qcdqq = qcdqq * gs4
        qcdgg = qcdgg * gs4

        ddsigma = 0.d0
        if (final_state <= 0) then
            if (verbose) print*, "scattering: summing over 2->2 |m|^2 with pdfs of initial partons ..."
            do i = -1, 1, 2
                do j = -1, 1, 2
                    dsigmapol(i,j) = fx1(13) * fx2(13) * qcdpolgg(i,j) &
                                   + fx1( 1) * fx2( 7) * (qcdpolqq(i,j) + qfdpoldd1(i,j)) &
                                   + fx1( 2) * fx2( 8) * (qcdpolqq(i,j) + qfdpoluu1(i,j)) &
                                   + fx1( 3) * fx2( 9) * (qcdpolqq(i,j) + qfdpoldd1(i,j)) &
                                   + fx1( 4) * fx2(10) * (qcdpolqq(i,j) + qfdpoluu1(i,j)) &
                                   + fx1( 5) * fx2(11) * (qcdpolqq(i,j) + qfdpoldd1(i,j)) &
                                   + fx1( 7) * fx2( 1) * (qcdpolqq(i,j) + qfdpoldd2(i,j)) &
                                   + fx1( 8) * fx2( 2) * (qcdpolqq(i,j) + qfdpoluu2(i,j)) &
                                   + fx1( 9) * fx2( 3) * (qcdpolqq(i,j) + qfdpoldd2(i,j)) &
                                   + fx1(10) * fx2( 4) * (qcdpolqq(i,j) + qfdpoluu2(i,j)) &
                                   + fx1(11) * fx2( 5) * (qcdpolqq(i,j) + qfdpoldd2(i,j))
                    if (ix == 1) then
                        dsigmapol(i, j) = dsigmapol(i, j) / x1
                    else if (ix == 2) then
                        dsigmapol(i, j) = dsigmapol(i, j) / x2
                    end if
                    ddsigma = ddsigma + dsigmapol(i,j)
                end do
            end do
        else if (final_state > 0) then
            if (verbose) print*, "scattering: summing over 2->6 |m|^2 with PDFs of initial partons ..."
            ddsigma = fx1(13) * fx2(13) * qcdgg &
                    + fx1( 1) * fx2( 7) * (qcdqq + qfddd1) &
                    + fx1( 2) * fx2( 8) * (qcdqq + qfduu1) &
                    + fx1( 3) * fx2( 9) * (qcdqq + qfddd1) &
                    + fx1( 4) * fx2(10) * (qcdqq + qfduu1) &
                    + fx1( 5) * fx2(11) * (qcdqq + qfddd1) &
                    + fx1( 7) * fx2( 1) * (qcdqq + qfddd2) &
                    + fx1( 8) * fx2( 2) * (qcdqq + qfduu2) &
                    + fx1( 9) * fx2( 3) * (qcdqq + qfddd2) &
                    + fx1(10) * fx2( 4) * (qcdqq + qfduu2) &
                    + fx1(11) * fx2( 5) * (qcdqq + qfddd2)
            if (ix == 1) then
                ddsigma = ddsigma / x1
            else if (ix  ==  2) then
                ddsigma = ddsigma / x2
            end if
        end if

        if (ddsigma == 0.d0) return

        if (final_state < 1) then
            do i = -1, 1, 2
                do j = -1, 1, 2
                    dsigmapol(i, j) = dsigmapol(i, j) / ddsigma
                end do
            end do
        end if

        666 continue

        if (verbose) print*, "scattering: multiplying by Jacobian for dx1 dx2 -> dx(2) dx(3) ..."
        ddsigma = ddsigma * (1.d0 - tau) * 2.d0 * ecm / s * (ecm_max - ecm_min)

        if (verbose) print*, "scattering: applying unit conversion ..."
        ddsigma = ddsigma * unit_conv

        if (verbose) print*, "scattering: multiply by phase space, azimuthal integration and flux factors ..."
        if (final_state <= 0) then
            if (use_rambo) then
                ddsigma = ddsigma * wgtr
            else
                ddsigma = ddsigma * qcm / (2.d0 * pcm) * 2.d0 ** (4 - 3 * 2) * 2.d0 * pi
            end if
            ddsigma = ddsigma / 2.d0 / ecm / ecm * (2.d0 * pi) ** (4 - 3 * 2)

        else if (final_state > 0) then
            if (use_rambo) then
                ddsigma = ddsigma * wgtr
            else
                ddsigma = ddsigma * q * rq56 * rq78 * rq5 * rq7/ecm * 256.d0 * 2.d0**(4 - 3 * (6)) * 2.d0 * pi
                if (map_phase_space) then
                    ddsigma = ddsigma * ((m356 * m356 - mt * mt)**2 + mt * mt * gamt * gamt) &
                                      * (at356max - at356min) / m356 / mt / gamt / 2.d0
                    ddsigma = ddsigma * ((m478 * m478 - mt * mt)**2 + mt * mt * gamt * gamt) &
                                      * (at478max - at478min) / m478 / mt / gamt / 2.d0
                    ddsigma = ddsigma * ((m56 * m56 - wmass * wmass)**2 + wmass * wmass * wwidth * wwidth) &
                                      * (at56max - at56min) / m56 / wmass / wwidth / 2.d0
                    ddsigma = ddsigma * ((m78 * m78 - wmass * wmass)**2 + wmass * wmass * wwidth * wwidth) &
                                      * (at78max - at78min) / m78 / wmass / wwidth / 2.d0
                    ddsigma = ddsigma * gamt / gamma_t * gamt / gamma_t ! nwa
                else
                    ddsigma = ddsigma * (m356max - m356min)
                    ddsigma = ddsigma * (m478max - m478min)
                    ddsigma = ddsigma * (m56max - m56min)
                    ddsigma = ddsigma * (m78max - m78min)
                end if
            end if

            ddsigma = ddsigma / 2.d0 / ecm / ecm * (2.d0 * pi)**(4 - 3 * (6))
        end if

        if (symmetrise) ddsigma = ddsigma / 2

        if (record_events) then
            if (ntuple_out) then
                if (verbose) print*, "scattering: write final particle collider frame momenta to n-tuple..."

                if (final_state == -1) then
                    call rootaddparticle(11,  pcol(1, 3), pcol(2, 3), pcol(3, 3), pcol(4, 3))
                    call rootaddparticle(-11, pcol(1, 4), pcol(2, 4), pcol(3, 4), pcol(4, 4))

                else if (final_state == 0) then
                    call rootaddparticle(6,   pcol(1, 3), pcol(2, 3), pcol(3, 3), pcol(4, 3))
                    call rootaddparticle(-6,  pcol(1, 4), pcol(2, 4), pcol(3, 4), pcol(4, 4))

                else if (final_state == 1) then
                    call rootaddparticle(5,   pcol(1, 3), pcol(2, 3), pcol(3, 3), pcol(4, 3))
                    call rootaddparticle(-5,  pcol(1, 4), pcol(2, 4), pcol(3, 4), pcol(4, 4))
                    call rootaddparticle(-11, pcol(1, 5), pcol(2, 5), pcol(3, 5), pcol(4, 5))
                    call rootaddparticle(12,  pcol(1, 6), pcol(2, 6), pcol(3, 6), pcol(4, 6))
                    call rootaddparticle(11,  pcol(1, 7), pcol(2, 7), pcol(3, 7), pcol(4, 7))
                    call rootaddparticle(-12, pcol(1, 8), pcol(2, 8), pcol(3, 8), pcol(4, 8))
                end if

                if (final_state < 1) then

                    ! Compute polarised event weightings
                    do i = -1, 1, 2
                        do j = -1, 1, 2
                            sigma_pol(i, j) = sigma_pol(i, j) + ddsigma * dsigmapol(i, j)
                            weight_pol(i, j) = ddsigma * dsigmapol(i, j)
                            error_pol(i, j) = error_pol(i, j) + sigma_pol(i, j)**2
                        end do
                    end do

                    call rootadddouble(weight_pol(-1, -1), "weightLL")
                    call rootadddouble(weight_pol(-1,  1), "weightLR")
                    call rootadddouble(weight_pol( 1, -1), "weightRL")
                    call rootadddouble(weight_pol( 1,  1), "weightRR")
                end if

                call rootaddevent(1.d0)
            end if

            if (lhef_out) then
                if (verbose) print*, "scattering: writing final particle collider frame momenta ..."
                pcol56 = pcol(1:4, 5) + pcol(1:4, 6)
                pcol78 = pcol(1:4, 7) + pcol(1:4, 8)
                pcol356 = pcol56 + pcol(1:4, 3)
                pcol478 = pcol78 + pcol(1:4, 4)

                a_s = alfas(zmass, lambdaqcd4, nloops)
                 
                if (final_state == 1) then
                    if (include_gg) then

                        if (verbose) print*, "lhe: g g -> t t~ -> b W+ b~ W- -> b b~ e+ ve e- ve~"

                        call lhe_add_event(12, 11, 1.d0, qq, a_em, a_s)

                        call lhe_add_particle( 21, -1,  0,  0, 101, 102, pcol(1:4,1), 0.d0, 9.d0) ! 1:  g
                        call lhe_add_particle( 21, -1,  0,  0, 103, 101, pcol(1:4,2), 0.d0, 9.d0) ! 2:  g
                        call lhe_add_particle(  6,  2,  1,  2, 103,   0, pcol356    , 0.d0, 9.d0) ! 3:  t
                        call lhe_add_particle( -6,  2,  1,  2,   0, 102, pcol478    , 0.d0, 9.d0) ! 4:  t~
                        call lhe_add_particle( 24,  2,  3,  3,   0,   0, pcol56     , 0.d0, 9.d0) ! 5:  W+
                        call lhe_add_particle(  5,  1,  3,  3, 103,   0, pcol(1:4,3), 0.d0, 9.d0) ! 6:  b
                        call lhe_add_particle(-24,  2,  4,  4,   0,   0, pcol78     , 0.d0, 9.d0) ! 7:  W-
                        call lhe_add_particle( -5,  1,  4,  4,   0, 102, pcol(1:4,4), 0.d0, 9.d0) ! 8:  b~
                        call lhe_add_particle(-11,  1,  6,  6,   0,   0, pcol(1:4,5), 0.d0, 9.d0) ! 9:  e+
                        call lhe_add_particle( 12,  1,  6,  6,   0,   0, pcol(1:4,6), 0.d0, 9.d0) ! 10: ve
                        call lhe_add_particle( 11,  1, 10, 10,   0,   0, pcol(1:4,7), 0.d0, 9.d0) ! 11: e-
                        call lhe_add_particle(-12,  1, 10, 10,   0,   0, pcol(1:4,8), 0.d0, 9.d0) ! 12: ve~

                    else if (include_qq) then

                        if (verbose) print*, "lhe:  q q~ -> t t~ -> b  W+ b~ W- -> b b~ e+ ve e- ve~"

                        call lhe_add_event(12, 12, 1.d0, qq, a_em, a_s)

                        call lhe_add_particle(  2, -1,  0,  0, 101,   0, pcol(1:4,1), 0.d0, 9.d0) ! 1:  u
                        call lhe_add_particle( -2, -1,  0,  0,   0, 102, pcol(1:4,2), 0.d0, 9.d0) ! 2:  u~
                        call lhe_add_particle(  6,  2,  1,  2, 101,   0, pcol356    , 0.d0, 9.d0) ! 3:  t
                        call lhe_add_particle( -6,  2,  1,  2,   0, 102, pcol478    , 0.d0, 9.d0) ! 4:  t~
                        call lhe_add_particle( 24,  2,  3,  0,   0,   0, pcol56     , 0.d0, 9.d0) ! 5:  W+
                        call lhe_add_particle(  5,  1,  3,  0, 101,   0, pcol(1:4,3), 0.d0, 9.d0) ! 6:  b
                        call lhe_add_particle(-24,  2,  4,  0,   0,   0, pcol78     , 0.d0, 9.d0) ! 7:  W-
                        call lhe_add_particle( -5,  1,  4,  0,   0, 102, pcol(1:4,4), 0.d0, 9.d0) ! 8:  b~
                        call lhe_add_particle(-11,  1,  5,  0,   0,   0, pcol(1:4,5), 0.d0, 9.d0) ! 9:  e+
                        call lhe_add_particle( 12,  1,  5,  0,   0,   0, pcol(1:4,6), 0.d0, 9.d0) ! 10: ve
                        call lhe_add_particle( 11,  1,  7,  0,   0,   0, pcol(1:4,7), 0.d0, 9.d0) ! 11: e-
                        call lhe_add_particle(-12,  1,  7,  0,   0,   0, pcol(1:4,8), 0.d0, 9.d0) ! 12: ve~

                    else if (include_uu) then

                        if (verbose) print*, "lhe: u u~ -> t t~ -> b b~ W+ W- -> b b~ e+ ve e- ve~"

                        call lhe_add_event(12, 13, 1.d0, qq, a_em, a_s)

                        call lhe_add_particle(  2, -1,  0,  0, 101,   0, pcol(1:4,1), 0.d0, 9.d0) ! 1:  u
                        call lhe_add_particle( -2, -1,  0,  0,   0, 101, pcol(1:4,2), 0.d0, 9.d0) ! 2:  u~
                        call lhe_add_particle(  6,  2,  1,  2, 102,   0, pcol356    , 0.d0, 9.d0) ! 3:  t
                        call lhe_add_particle( -6,  2,  1,  2,   0, 102, pcol478    , 0.d0, 9.d0) ! 4:  t~
                        call lhe_add_particle( 24,  2,  3,  0,   0,   0, pcol56     , 0.d0, 9.d0) ! 5:  W+
                        call lhe_add_particle(  5,  1,  3,  0, 102,   0, pcol(1:4,3), 0.d0, 9.d0) ! 6:  b
                        call lhe_add_particle(-24,  2,  4,  0,   0,   0, pcol78     , 0.d0, 9.d0) ! 7:  W-
                        call lhe_add_particle( -5,  1,  4,  0,   0, 102, pcol(1:4,4), 0.d0, 9.d0) ! 8:  b~
                        call lhe_add_particle(-11,  1,  5,  0,   0,   0, pcol(1:4,5), 0.d0, 9.d0) ! 9:  e+
                        call lhe_add_particle( 12,  1,  5,  0,   0,   0, pcol(1:4,6), 0.d0, 9.d0) ! 10: ve
                        call lhe_add_particle( 11,  1,  7,  0,   0,   0, pcol(1:4,7), 0.d0, 9.d0) ! 11: e-
                        call lhe_add_particle(-12,  1,  7,  0,   0,   0, pcol(1:4,8), 0.d0, 9.d0) ! 12: ve~

                    else if (include_dd) then

                        if (verbose) print*, "lhe: d d~ -> t t~ -> b b~ W+ W- -> b b~ e+ ve e- ve~"

                        call lhe_add_event(12, 14, 1.d0, qq, a_em, a_s)

                        call lhe_add_particle(  1, -1,  0,  0, 101,   0, pcol(1:4,1), 0.d0, 9.d0) ! 1:  d
                        call lhe_add_particle( -1, -1,  0,  0,   0, 101, pcol(1:4,2), 0.d0, 9.d0) ! 2:  d~
                        call lhe_add_particle(  6,  2,  1,  2, 102,   0, pcol356    , 0.d0, 9.d0) ! 3:  t
                        call lhe_add_particle( -6,  2,  1,  2,   0, 102, pcol478    , 0.d0, 9.d0) ! 4:  t~
                        call lhe_add_particle( 24,  2,  3,  0,   0,   0, pcol56     , 0.d0, 9.d0) ! 5:  W+
                        call lhe_add_particle(  5,  1,  3,  0, 102,   0, pcol(1:4,3), 0.d0, 9.d0) ! 6:  b
                        call lhe_add_particle(-24,  2,  4,  0,   0,   0, pcol78     , 0.d0, 9.d0) ! 7:  W-
                        call lhe_add_particle( -5,  1,  4,  0,   0, 102, pcol(1:4,4), 0.d0, 9.d0) ! 8:  b~
                        call lhe_add_particle(-11,  1,  5,  0,   0,   0, pcol(1:4,5), 0.d0, 9.d0) ! 9:  e+
                        call lhe_add_particle( 12,  1,  5,  0,   0,   0, pcol(1:4,6), 0.d0, 9.d0) ! 10: ve
                        call lhe_add_particle( 11,  1,  7,  0,   0,   0, pcol(1:4,7), 0.d0, 9.d0) ! 11: e-
                        call lhe_add_particle(-12,  1,  7,  0,   0,   0, pcol(1:4,8), 0.d0, 9.d0) ! 12: ve~

                    end if
                else if (final_state == 2) then
                    if (include_gg) then

                        if (verbose) print*, "lhe: g g -> t t~ -> b W+ b~ W- -> b b~ e+ ve u~ d"

                        call lhe_add_event(12, 21, 1.d0, qq, a_em, a_s)

                        call lhe_add_particle( 21, -1,  0,  0, 101, 102, pcol(1:4,1), 0.d0, 9.d0) ! 1:  g
                        call lhe_add_particle( 21, -1,  0,  0, 103, 101, pcol(1:4,2), 0.d0, 9.d0) ! 2:  g
                        call lhe_add_particle(  6,  2,  1,  2, 103,   0, pcol356    , 0.d0, 9.d0) ! 3:  t
                        call lhe_add_particle( -6,  2,  1,  2,   0, 102, pcol478    , 0.d0, 9.d0) ! 4:  t~
                        call lhe_add_particle( 24,  2,  3,  0,   0,   0, pcol56     , 0.d0, 9.d0) ! 5:  W+
                        call lhe_add_particle(  5,  1,  3,  0, 103,   0, pcol(1:4,3), 0.d0, 9.d0) ! 6:  b
                        call lhe_add_particle(-24,  2,  4,  0,   0,   0, pcol78     , 0.d0, 9.d0) ! 7:  W-
                        call lhe_add_particle( -5,  1,  4,  0,   0, 102, pcol(1:4,4), 0.d0, 9.d0) ! 8:  b~
                        call lhe_add_particle(-11,  1,  5,  0,   0,   0, pcol(1:4,5), 0.d0, 9.d0) ! 9:  e+
                        call lhe_add_particle( 12,  1,  5,  0,   0,   0, pcol(1:4,6), 0.d0, 9.d0) ! 10: ve
                        call lhe_add_particle(  1,  1,  7,  0, 104,   0, pcol(1:4,7), 0.d0, 9.d0) ! 11: u~
                        call lhe_add_particle( -2,  1,  7,  0,   0, 104, pcol(1:4,8), 0.d0, 9.d0) ! 12: d

                    else if (include_qq) then

                        if (verbose) print*, "lhe:  q q~ -> t t~ -> b  W+ b~ W- -> b b~ e+ ve u~ d"

                        call lhe_add_event(12, 22, 1.d0, qq, a_em, a_s)

                        call lhe_add_particle(  2, -1,  0,  0, 101,   0, pcol(1:4,1), 0.d0, 9.d0) ! 1:  u
                        call lhe_add_particle( -2, -1,  0,  0,   0, 102, pcol(1:4,2), 0.d0, 9.d0) ! 2:  u~
                        call lhe_add_particle(  6,  2,  1,  2, 101,   0, pcol356    , 0.d0, 9.d0) ! 3:  t
                        call lhe_add_particle( -6,  2,  1,  2,   0, 102, pcol478    , 0.d0, 9.d0) ! 4:  t~
                        call lhe_add_particle( 24,  2,  3,  0,   0,   0, pcol56     , 0.d0, 9.d0) ! 5:  W+
                        call lhe_add_particle(  5,  1,  3,  0, 101,   0, pcol(1:4,3), 0.d0, 9.d0) ! 6:  b
                        call lhe_add_particle(-24,  2,  4,  0,   0,   0, pcol78     , 0.d0, 9.d0) ! 7:  W-
                        call lhe_add_particle( -5,  1,  4,  0,   0, 102, pcol(1:4,4), 0.d0, 9.d0) ! 8:  b~
                        call lhe_add_particle(-11,  1,  5,  0,   0,   0, pcol(1:4,5), 0.d0, 9.d0) ! 9:  e+
                        call lhe_add_particle( 12,  1,  5,  0,   0,   0, pcol(1:4,6), 0.d0, 9.d0) ! 10: ve
                        call lhe_add_particle( 11,  1,  7,  0, 103,   0, pcol(1:4,7), 0.d0, 9.d0) ! 11: u~
                        call lhe_add_particle(-12,  1,  7,  0,   0, 103, pcol(1:4,8), 0.d0, 9.d0) ! 12: d

                    else if (include_uu) then

                        if (verbose) print*, "lhe: u u~ -> t t~ -> b b~ W+ W- -> b b~ e+ ve u~ d"

                        call lhe_add_event(12, 23, 1.d0, qq, a_em, a_s)

                        call lhe_add_particle(  2, -1,  0,  0, 101,   0, pcol(1:4,1), 0.d0, 9.d0) ! 1:  u
                        call lhe_add_particle( -2, -1,  0,  0,   0, 101, pcol(1:4,2), 0.d0, 9.d0) ! 2:  u~
                        call lhe_add_particle(  6,  2,  1,  2, 102,   0, pcol356    , 0.d0, 9.d0) ! 3:  t
                        call lhe_add_particle( -6,  2,  1,  2,   0, 102, pcol478    , 0.d0, 9.d0) ! 4:  t~
                        call lhe_add_particle( 24,  2,  3,  0,   0,   0, pcol56     , 0.d0, 9.d0) ! 5:  W+
                        call lhe_add_particle(  5,  1,  3,  0, 102,   0, pcol(1:4,3), 0.d0, 9.d0) ! 6:  b
                        call lhe_add_particle(-24,  2,  4,  0,   0,   0, pcol78     , 0.d0, 9.d0) ! 7:  W-
                        call lhe_add_particle( -5,  1,  4,  0,   0, 102, pcol(1:4,4), 0.d0, 9.d0) ! 8:  b~
                        call lhe_add_particle( -1,  1,  5,  0,   0,   0, pcol(1:4,5), 0.d0, 9.d0) ! 9:  e+
                        call lhe_add_particle(  2,  1,  5,  0,   0,   0, pcol(1:4,6), 0.d0, 9.d0) ! 10: ve
                        call lhe_add_particle(  1,  1,  7,  0, 103,   0, pcol(1:4,7), 0.d0, 9.d0) ! 11: u~
                        call lhe_add_particle( -2,  1,  7,  0,   0, 103, pcol(1:4,8), 0.d0, 9.d0) ! 12: d

                    else if (include_dd) then

                        if (verbose) print*, "lhe: d d~ -> t t~ -> b b~ W+ W- -> b b~ e+ ve u~ d"

                        call lhe_add_event(12, 24, 1.d0, qq, a_em, a_s)

                        call lhe_add_particle(  1, -1,  0,  0, 101,   0, pcol(1:4,1), 0.d0, 9.d0) ! 1:  d
                        call lhe_add_particle( -1, -1,  0,  0,   0, 101, pcol(1:4,2), 0.d0, 9.d0) ! 2:  d~
                        call lhe_add_particle(  6,  2,  1,  2, 102,   0, pcol356    , 0.d0, 9.d0) ! 3:  t
                        call lhe_add_particle( -6,  2,  1,  2,   0, 102, pcol478    , 0.d0, 9.d0) ! 4:  t~
                        call lhe_add_particle( 24,  2,  3,  0,   0,   0, pcol56     , 0.d0, 9.d0) ! 5:  W+
                        call lhe_add_particle(  5,  1,  3,  0, 102,   0, pcol(1:4,3), 0.d0, 9.d0) ! 6:  b
                        call lhe_add_particle(-24,  2,  4,  0,   0,   0, pcol78     , 0.d0, 9.d0) ! 7:  W-
                        call lhe_add_particle( -5,  1,  4,  0,   0, 102, pcol(1:4,4), 0.d0, 9.d0) ! 8:  b~
                        call lhe_add_particle( -1,  1,  5,  0,   0,   0, pcol(1:4,5), 0.d0, 9.d0) ! 9:  e+
                        call lhe_add_particle(  2,  1,  5,  0,   0,   0, pcol(1:4,6), 0.d0, 9.d0) ! 10: ve
                        call lhe_add_particle(  1,  1,  7,  0, 103,   0, pcol(1:4,7), 0.d0, 9.d0) ! 11: u~
                        call lhe_add_particle( -2,  1,  7,  0,   0, 103, pcol(1:4,8), 0.d0, 9.d0) ! 12: d

                    end if
                end if

                call lhe_end_event
            end if
        end if
        if (verbose) print*, "ddsigma: ", ddsigma
        dsigma = dsigma + ddsigma
    end do
    if (verbose) print*, "dsigma: ", dsigma
    return
end function dsigma

end module scattering
