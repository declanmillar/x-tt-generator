module scattering

    use kinds
    use configuration
    use modelling
    use tt_bbeevv
    use vampi

    implicit none

    public :: dsigma, initialise_masses, initialise_s, initialise_pdfs, set_energy_limits, phi
    logical, public :: record_events
    real(kind=default), public :: sigma_pol(-1:1, -1:1), error_pol(-1:1, -1:1)
    real(kind=default), private :: m3, m4, m5, m6, m7, m8, mt, gamt
    real(kind=default), private :: s, ecm_max, ecm_min, scale, a_s

contains

pure function phi (xi, channel) result (x)
    real(kind=default), dimension(:), intent(in) :: xi
    integer, intent(in) :: channel
    real(kind=default), dimension(size(xi)) :: x
    x = xi
end function phi

subroutine initialise_pdfs

    character(40) tablefile
    real(kind=default) :: alfas

    if (verbose) print*, "scattering: opening pdf tables ..."
    if (pdf <=  4) call setctq6(pdf)
    if (pdf == 10) tablefile = "ct14ln.pds"
    if (pdf == 11) tablefile = "ct14ll.pds"
    if (pdf >= 10) call setct14(tablefile)

    if (verbose) print*, "scattering: setting lambda_QCD based on PDF set ..."
    if (pdf ==  1) lambdaqcd4 = 0.326d0
    if (pdf ==  2) lambdaqcd4 = 0.326d0
    if (pdf ==  3) lambdaqcd4 = 0.326d0
    if (pdf ==  4) lambdaqcd4 = 0.215d0
    if (pdf ==  5) lambdaqcd4 = 0.300d0
    if (pdf ==  6) lambdaqcd4 = 0.300d0
    if (pdf ==  7) lambdaqcd4 = 0.300d0
    if (pdf ==  8) lambdaqcd4 = 0.229d0
    if (pdf ==  9) lambdaqcd4 = 0.383d0
    if (pdf == 10) lambdaqcd4 = 0.326d0
    if (pdf == 11) lambdaqcd4 = 0.215d0
    write(*,"(a16,f5.3)") "Lambda_QCD^4 = ", lambdaqcd4

    if (verbose) print*, "scattering: setting n_loops ..."
    if (pdf <= 2 .or. pdf == 10) then
        nloops = 2
    else
        nloops = 1
    end if
    write(*,"(a16,i1)")   "loops = ", nloops

    if (verbose) print*, "scattering: calculating  alpha_s(zmass) ..."
    a_s = alfas(zmass, lambdaqcd4, nloops)
    write(*,"(a16,f5.3)") "alpha_s(m_Z) = ", a_s

    ! scale for the pdfs
    if (final_state >= 0) then
        write(*,"(a25)") "Q = 2 * m_top"
        scale = 2.d0 * mt
    else
        write(*,"(a25)") "Q = Ecm"
        scale = 0.d0
    end if

end subroutine initialise_pdfs

subroutine initialise_masses

    if (verbose) print*, "scattering: setting external particle masses ..."
    if (final_state < 1) then
        m3 = fmass(ffinal)
        m4 = fmass(ffinal)
    else
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
    if (verbose) print*, "scattering: setting energy limits ..."
    ecm_min = m3 + m4 + m5 + m6 + m7 + m8
    ecm_max = sqrts
    if (ecm_low > ecm_min) then
        ecm_min = ecm_low
    end if
    if (ecm_up > 0) then
        ecm_max = ecm_up
    end if
    if (verbose) print*, "scattering: ecm_min = ", ecm_min
    if (verbose) print*, "scattering: ecm_max = ", ecm_max
end subroutine set_energy_limits


subroutine initialise_s

    if (verbose) print*, "scattering: setting s ..."

    s = sqrts * sqrts

    if (verbose) print*, "scattering: s = ", s

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
    real(kind=default), parameter :: unit_conv = 0.38937966d9 ! GeV^{-2} to nb (pb?)

    ! alphas
    real(kind=default) :: alfas, gs4, gs2, a_s

    ! square matrix elements
    real(kind=default) :: suu1, suu2, sdd1, sdd2, sqq, sgg
    real(kind=default) :: sgg_tt, sqq_tt, sqq_ff
    real(kind=default) :: sgg_bbemuvevm, sqq_bbemuvevm, suu_bbemuvevm, sdd_bbemuvevm
    real(kind=default) :: sgg_bbtatavtvt, sqq_bbtatavtvt, suu_bbtatavtvt, sdd_bbtatavtvt

    ! pdfs
    real(kind=default) :: ctq6pdf, ct14pdf
    real(kind=default) :: fx1(13), fx2(13), x1, x2, xx1, xx2, x1x2(2,2)
    real(kind=default) :: d1, u1, str1, chm1, btm1, glu1
    real(kind=default) :: d2, u2, str2, chm2, btm2, glu2
    real(kind=default) :: dbar1, ubar1, dbar2, ubar2
    real(kind=default) :: dsea1, usea1, dsea2, usea2

    ! energies
    real(kind=default) :: shat, tau, ecm

    real(kind=default) :: qcm, pcm, qcm2
    real(kind=default) :: pq5, pq52, pq56, pq7, pq78
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
    real(kind=default) :: eta, pt

    ! arctan
    real(kind=default) :: at356, at356max, at356min, at56, at56max, at56min, at34min, at34
    real(kind=default) :: at478, at478max, at478min, at78, at78max, at78min, at34max

    ! 4-momenta
    real(kind=default) :: p(4,8), p356(4), p478(4), q56(4), q78(4), p56(4), p78(4), q5(4), q7(4)
    real(kind=default) :: pcol(4,8), pcol56(4), pcol78(4), pcol356(4), pcol478(4)
    real(kind=default) :: p1(0:3), p2(0:3), p3(0:3), p4(0:3), p5(0:3), p6(0:3), p7(0:3), p8(0:3)

    ! invariant masses
    real(kind=default) :: m356, m356max, m356min, m478, m478max, m478min
    real(kind=default) :: m356_2, m478_2, m56, m56_2, m56max, m56min, m78, m78_2, m78max, m78min

    ! polarised
    real(kind=default) :: spolqq(-1:1, -1:1), spolgg(-1:1, -1:1)
    real(kind=default) :: spoluu1(-1:1, -1:1), spoldd1(-1:1, -1:1), spolbb1(-1:1, -1:1)
    real(kind=default) :: spoluu2(-1:1, -1:1), spoldd2(-1:1, -1:1), spolbb2(-1:1, -1:1)
    real(kind=default) :: weight_pol(-1:1, -1:1), dsigmapol(-1:1, -1:1)

    integer :: seed
    real(kind=default) :: random

    ! temporary iterators
    integer :: i, j

    ! ---

    if (verbose) print*, "dsigma: begin"

    dsigma = 0

    if (present(channel) .and. (include_dd .or. include_uu) .and. (channel == 2 .or. channel == 3)) then
        if (channel == 2) then
            if (verbose) print*, "dsigma: flattening Breit Wigner" 
            at34min = atan((ecm_min * ecm_min - zmass * zmass) / (zmass * zwidth))
            at34max = atan((ecm_max * ecm_max - zmass * zmass) / (zmass * zwidth))
            if (use_rambo) then
                at34 = x(1) * (at34max - at34min) + at34min
            else 
                at34 = x(2 + 12 * tops_decay) * (at34max - at34min) + at34min
            end if
            shat = zmass * zmass + tan(at34) * zmass * zwidth
            if (shat < 0.d0) then
                if (verbose) print*, "invalid"
                return
            end if
            ecm = sqrt(shat)
        else if (channel == 3) then
            if (verbose) print*, "dsigma: flattening Breit Wigner" 
            at34min = atan((ecm_min * ecm_min - xmass(1) * xmass(1)) / (xmass(1) * xwidth(1)))
            at34max = atan((ecm_max * ecm_max - xmass(1) * xmass(1)) / (xmass(1) * xwidth(1)))
            if (use_rambo) then
                at34 = x(1) * (at34max - at34min) + at34min
            else 
                at34 = x(2 + 12 * tops_decay) * (at34max - at34min) + at34min
            end if
            shat = xmass(1) * xmass(1) + tan(at34) * xmass(1) * xwidth(1)
            if (shat < 0.d0) then
                if (verbose) print*, "invalid"
                return
            end if
            ecm = sqrt(shat)
        end if
    else
        if (verbose) print*, "dsigma: not flattening Breit Wigner" 
        if (use_rambo) then
            ecm = x(1) * (ecm_max - ecm_min) + ecm_min
        else 
            ecm = x(2 + 12 * tops_decay) * (ecm_max - ecm_min) + ecm_min
        end if
        shat = ecm * ecm
    end if

    if (verbose) print*, "dsigma: ecm calculated" 
    
    tau = shat / s

    if (scale == 0.d0) scale = ecm

    ! x1 and x2 of the partons
    if (use_rambo) then
        x1 = x(2) * (1.d0 - tau) + tau
    else
        x1 = x(3 + 12 * tops_decay) * (1.d0 - tau) + tau
    end if
    x2 = tau / x1

    ! initialisation
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
    if (pdf <= 4) then
        if ((x1 <= 1.d-6) .or. (x1 >= 1.d0)) then
            if (verbose) print*, "invalid x1"
            return
        end if
        if ((x2 <= 1.d-6) .or. (x2 >= 1.d0)) then
            if (verbose) print*, "invalid x2"
            return
        end if
        if ((scale <= 1.3d0) .or. (scale >= 1.d4)) then
            if (verbose) print*, "invalid scale"
            return
        end if

        ! *x for compatibility with MRS which return xf(x)
        u1 = x1 * ctq6pdf(1, x1, scale)
        d1 = x1 * ctq6pdf(2, x1, scale)
        ubar1 = x1 * ctq6pdf(-1, x1, scale)
        dbar1 = x1 * ctq6pdf(-2, x1, scale)
        str1 = x1 * ctq6pdf(3, x1, scale)
        chm1 = x1 * ctq6pdf(4, x1, scale)
        btm1 = x1 * ctq6pdf(5, x1, scale)
        glu1 = x1 * ctq6pdf(0, x1, scale)
        u2 = x2 * ctq6pdf(1, x2, scale)
        d2 = x2 * ctq6pdf(2, x2, scale)
        ubar2 = x2 * ctq6pdf(-1, x2, scale)
        dbar2 = x2 * ctq6pdf(-2, x2, scale)
        str2 = x2 * ctq6pdf(3, x2, scale)
        chm2 = x2 * ctq6pdf(4, x2, scale)
        btm2 = x2 * ctq6pdf(5, x2, scale)
        glu2 = x2 * ctq6pdf(0, x2, scale)

    else if (pdf > 4 .and. pdf < 10) then
        if ((x1 <= 1.d-5) .or. (x1 >= 1.d0)) return
        if ((x2 <= 1.d-5) .or. (x2 >= 1.d0)) return
        if ((scale * scale <= 1.25d0) .or. (scale * scale >= 1.d7)) return

        call mrs99(x1, scale, pdf-4, u1, d1, usea1, dsea1, str1, chm1, btm1, glu1)
        call mrs99(x2, scale, pdf-4, u2, d2, usea2, dsea2, str2, chm2, btm2, glu2)

        u1 = u1 + usea1
        d1 = d1 + dsea1
        u2 = u2 + usea2
        d2 = d2 + dsea2
        ubar1 = usea1
        dbar1 = dsea1
        ubar2 = usea2
        dbar2 = dsea2
    else if (pdf > 9) then
        if ((x1 <= 1.d-9) .or. (x1 >= 1.d0)) then
            if (verbose) print*, "invalid x1"
            return
        end if
        if ((x2 <= 1.d-9) .or. (x2 >= 1.d0)) then
            if (verbose) print*, "invalid x2"
            return
        end if
        if ((scale <= 1.3d0) .or. (scale >= 1.d5)) then
            if (verbose) print*, "invalid scale"
            return
        end if

        ! *x for compatibility with MRS which return xf(x)
        u1 = x1 * ct14pdf(1, x1, scale)
        d1 = x1 * ct14pdf(2, x1, scale)
        ubar1 = x1 * ct14pdf(-1, x1, scale)
        dbar1 = x1 * ct14pdf(-2, x1, scale)
        str1 = x1 * ct14pdf(3, x1, scale)
        chm1 = x1 * ct14pdf(4, x1, scale)
        btm1 = x1 * ct14pdf(5, x1, scale)
        glu1 = x1 * ct14pdf(0, x1, scale)
        u2 = x2 * ct14pdf(1, x2, scale)
        d2 = x2 * ct14pdf(2, x2, scale)
        ubar2 = x2 * ct14pdf(-1, x2, scale)
        dbar2 = x2 * ct14pdf(-2, x2, scale)
        str2 = x2 * ct14pdf(3, x2, scale)
        chm2 = x2 * ct14pdf(4, x2, scale)
        btm2 = x2 * ct14pdf(5, x2, scale)
        glu2 = x2 * ct14pdf(0, x2, scale)
    end if

    if (verbose) print*, "initialising pdfs..."
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

    if (verbose) print*, "setting incoming 4-momenta ..."
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
        if (verbose) print*, "setting RAMBO random number seed ..."
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
            call random_number(random)
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
            phi = 2.d0 * pi * x(16)

            m356min = m3 + m5 + m6
            m356max = ecm - m4 - m7 - m8
            if (map_phase_space) then
                at356min = atan((m356min * m356min - mt * mt) / (mt * gamt))
                at356max = atan((m356max * m356max - mt * mt) / (mt * gamt))
                at356 = x(13) * (at356max - at356min) + at356min
                m356_2 = mt * mt + tan(at356) * mt * gamt
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
                at478min = atan((m478min * m478min - mt * mt) / (mt * gamt))
                at478max = atan((m478max * m478max - mt * mt) / (mt * gamt))
                at478 = x(12) * (at478max - at478min) + at478min
                m478_2 = mt * mt + tan(at478) * mt * gamt
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
                at56min = atan((m56min * m56min - wmass * wmass) / (wmass * wwidth))
                at56max = atan((m56max * m56max - wmass * wmass) / (wmass * wwidth))
                at56 = x(11) * (at56max - at56min) + at56min
                m56_2 = wmass * wmass + tan(at56) * wmass * wwidth
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
                at78min = atan((m78min * m78min - wmass * wmass) / (wmass * wwidth))
                at78max = atan((m78max * m78max - wmass * wmass) / (wmass * wwidth))
                at78 = x(10) * (at78max - at78min) + at78min
                m78_2 = wmass * wmass + tan(at78) * wmass * wwidth
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
        if (verbose) print*, "kinematics: applying detector cuts ..."
        do i = 3, nfinal
            pt = sqrt(pcol(1,i) * pcol(1,i) + pcol(2,i) * pcol(2,i))
            if (pt < 25.0) return
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
        dsigma = 1.d0 / x1
        if (final_state < 1) then
            do i = -1, 1, 2
                do j = -1, 1, 2
                  dsigmapol(i,j) = 1.d0 / (x1 * dsigma)
                end do
            end do
        end if
        go to 666
    end if

    if (verbose) print*,  "matrix elements: calculate QCD coupling ..."
    gs2 = 4.d0 * pi * alfas(scale, lambdaqcd4, nloops)
    gs4 = gs2 * gs2

    if (verbose) print*, "matrix elements: initialising ..."
    sqq = 0.d0
    sgg = 0.d0
    suu1 = 0.d0
    sdd1 = 0.d0
    suu2 = 0.d0
    sdd2 = 0.d0
    do i = -1, 1
        do j = -1, 1
            spolqq(i, j) = 0.d0
            spolgg(i, j) = 0.d0
            spoluu1(i, j) = 0.d0
            spoldd1(i, j) = 0.d0
            spoluu2(i, j) = 0.d0
            spoldd2(i, j) = 0.d0
            dsigmapol(i, j) = 0.d0
            weight_pol(i, j) = 0.d0
        end do
    end do

    if (final_state < 1) then
        if (verbose) print*, "matrix elements: calculating 2 -> 2 ..."
        do i = -1, 1, 2
            do j = -1, 1, 2
                if (include_gg) then
                    spolgg(i,j) = sgg_tt(p1, p2, p3, p4, i, j) * gs4
                end if
                if (include_qq) then
                    spolqq(i,j) = sqq_tt(3, p1, p2, p3, p4, i, j) * gs4
                end if
                if (include_uu) then
                    spoluu1(i, j) = sqq_ff(3, ffinal, p1, p2, p3, p4, i, j)
                    spoluu2(i, j) = sqq_ff(3, ffinal, p2, p1, p3, p4, i, j)
                end if
                if (include_dd) then
                    spoldd1(i, j) = sqq_ff(4, ffinal, p1, p2, p3, p4, i, j)
                    spoldd2(i, j) = sqq_ff(4, ffinal, p2, p1, p3, p4, i, j)
                end if
            end do
        end do

    else
        if (verbose) print*, "matrix elements: calculating 2 -> 6 ..."
        if (include_gg) then
            sgg = sgg_tt_bbeevv(p1, p2, p3, p4, p5, p7, p6, p8, channel)
        end if
        if (include_qq) then
            sqq = sqq_tt_bbeevv(3, p1, p2, p3, p4, p5, p7, p6, p8)
        end if
        if (include_uu) then
            suu1 = sqq_tt_bbeevv_ew(3, 11, p1, p2, p3, p4, p5, p7, p6, p8, channel)
            suu2 = sqq_tt_bbeevv_ew(3, 11, p2, p1, p3, p4, p5, p7, p6, p8, channel)
        end if
        if (include_dd) then
            sdd1 = sqq_tt_bbeevv_ew(4, 11, p1, p2, p3, p4, p5, p7, p6, p8, channel)
            sdd2 = sqq_tt_bbeevv_ew(4, 11, p2, p1, p3, p4, p5, p7, p6, p8, channel)
        end if
    end if

    if (verbose) print*, "matrix elements: multiply qcd |m|^2 by g_s^4 ..."
    sqq = sqq * gs4
    sgg = sgg * gs4

    dsigma = 0.d0
    if (final_state < 1) then
        if (verbose) print*, "scattering: summing over 2->2 |m|^2 with pdfs of initial partons ..."
        do i = -1, 1, 2
            do j = -1, 1, 2
                dsigmapol(i,j) = fx1(13) * fx2(13) *  spolgg(i, j) &
                               + fx1( 1) * fx2( 7) * (spolqq(i, j) + spoldd1(i, j)) &
                               + fx1( 2) * fx2( 8) * (spolqq(i, j) + spoluu1(i, j)) &
                               + fx1( 3) * fx2( 9) * (spolqq(i, j) + spoldd1(i, j)) &
                               + fx1( 4) * fx2(10) * (spolqq(i, j) + spoluu1(i, j)) &
                               + fx1( 5) * fx2(11) * (spolqq(i, j) + spoldd1(i, j)) &
                               + fx1( 7) * fx2( 1) * (spolqq(i, j) + spoldd2(i, j)) &
                               + fx1( 8) * fx2( 2) * (spolqq(i, j) + spoluu2(i, j)) &
                               + fx1( 9) * fx2( 3) * (spolqq(i, j) + spoldd2(i, j)) &
                               + fx1(10) * fx2( 4) * (spolqq(i, j) + spoluu2(i, j)) &
                               + fx1(11) * fx2( 5) * (spolqq(i, j) + spoldd2(i, j))
                dsigmapol(i, j) = dsigmapol(i, j) / x1
                dsigma = dsigma + dsigmapol(i, j)
            end do
        end do
    else if (final_state > 0) then
        if (verbose) print*, "scattering: summing over 2->6 |m|^2 with PDFs of initial partons ..."
        dsigma = fx1(13) * fx2(13) *  sgg &
               + fx1( 1) * fx2( 7) * (sqq + sdd1) &
               + fx1( 2) * fx2( 8) * (sqq + suu1) &
               + fx1( 3) * fx2( 9) * (sqq + sdd1) &
               + fx1( 4) * fx2(10) * (sqq + suu1) &
               + fx1( 5) * fx2(11) * (sqq + sdd1) &
               + fx1( 7) * fx2( 1) * (sqq + sdd2) &
               + fx1( 8) * fx2( 2) * (sqq + suu2) &
               + fx1( 9) * fx2( 3) * (sqq + sdd2) &
               + fx1(10) * fx2( 4) * (sqq + suu2) &
               + fx1(11) * fx2( 5) * (sqq + sdd2)
        dsigma = dsigma / x1
    end if

    if (dsigma == 0.d0) return

    if (final_state < 1) then
        do i = -1, 1, 2
            do j = -1, 1, 2
                dsigmapol(i, j) = dsigmapol(i, j) / dsigma
            end do
        end do
    end if 

    666 continue

    if (verbose) print*, "scattering: multiplying by Jacobian for dx1 dx2 -> dx(2) dx(3) ..."
    dsigma = dsigma * (1.d0 - tau) * 2.d0 * ecm / s * (ecm_max - ecm_min)

    if (verbose) print*, "scattering: applying unit conversion ..."
    dsigma = dsigma * unit_conv

    if (verbose) print*, "scattering: multiply by phase space, azimuthal integration and flux factors ..."
    if (final_state < 1) then
        if (use_rambo) then
            dsigma = dsigma * wgtr
        else
            dsigma = dsigma * qcm / (2.d0 * pcm) * 2.d0 ** (4 - 3 * 2) * 2.d0 * pi
        end if
        dsigma = dsigma / (2.d0 * ecm * ecm) * (2.d0 * pi) ** (4 - 3 * 2)

    else if (final_state > 0) then
        if (use_rambo) then
            dsigma = dsigma * wgtr
        else
            dsigma = dsigma * q * rq56 * rq78 * rq5 * rq7 / ecm * 256.d0 * 2.d0 ** (4 - 3 * 6) * 2.d0 * pi
            if (map_phase_space) then
                if (present(channel) .and. (include_dd .or. include_uu)) then

                    if (channel == 2) then
                        dsigma = dsigma * ((ecm*ecm - zmass*zmass)**2 + zmass*zmass*zwidth*zwidth) &
                                        * (at34max - at34min) / (ecm * zmass * zwidth * 2.d0)
                    else if (channel == 3) then
                        dsigma = dsigma * ((ecm*ecm - xmass(1)*xmass(1))**2 + xmass(1)*xmass(1)*xwidth(1)*xwidth(1)) &
                                        * (at34max - at34min) / (ecm * xmass(1) * xwidth(1) * 2.d0)
                    end if
                end if
                dsigma = dsigma * ((m356 * m356 - mt * mt)**2 + mt * mt * gamt * gamt) &
                                * (at356max - at356min) / (m356 * mt * gamt * 2.d0)
                dsigma = dsigma * ((m478 * m478 - mt * mt)**2 + mt * mt * gamt * gamt) &
                                * (at478max - at478min) / (m478 * mt * gamt * 2.d0)
                dsigma = dsigma * ((m56 * m56 - wmass * wmass)**2 + wmass * wmass * wwidth * wwidth) &
                                * (at56max - at56min) / (m56 * wmass * wwidth * 2.d0)
                dsigma = dsigma * ((m78 * m78 - wmass * wmass)**2 + wmass * wmass * wwidth * wwidth) &
                                * (at78max - at78min) / (m78 * wmass * wwidth * 2.d0)
                dsigma = dsigma * gamt * gamt / (gamma_t * gamma_t) ! nwa
            else
                dsigma = dsigma * (m356max - m356min)
                dsigma = dsigma * (m478max - m478min)
                dsigma = dsigma * (m56max - m56min)
                dsigma = dsigma * (m78max - m78min)
            end if
        end if

        dsigma = dsigma / (2.d0 * ecm * ecm) * (2.d0 * pi)**(4 - 3 * (6))
    end if

    if (verbose) print*, "dsigma: ", dsigma

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
                        sigma_pol(i, j) = sigma_pol(i, j) + dsigma * dsigmapol(i, j)
                        weight_pol(i, j) = dsigma * dsigmapol(i, j)
                        error_pol(i, j) = error_pol(i, j) + sigma_pol(i, j) * sigma_pol(i, j)
                    end do
                end do

                call rootadddouble(weight_pol(-1, -1), "weightLL")
                call rootadddouble(weight_pol(-1,  1), "weightLR")
                call rootadddouble(weight_pol( 1, -1), "weightRL")
                call rootadddouble(weight_pol( 1,  1), "weightRR")
            end if
        end if

        if (lhef_out) then
            if (verbose) print*, "scattering: writing final particle collider frame momenta ..."
            pcol56 = pcol(1:4, 5) + pcol(1:4, 6)
            pcol78 = pcol(1:4, 7) + pcol(1:4, 8)
            pcol356 = pcol56 + pcol(1:4, 3)
            pcol478 = pcol78 + pcol(1:4, 4)
             
            do i = 11, 15, 2
                do j = 11, 15, 2
                    call lhe_add_event(12, 9999, 1.d0, scale, a_em, a_s)
                    if (include_gg) then
                        call lhe_add_particle( 21, -1,  0,  0, 101, 102, pcol(1:4,1)) ! 1:  g
                        call lhe_add_particle( 21, -1,  0,  0, 102, 103, pcol(1:4,2)) ! 2:  g
                        call lhe_add_particle(  6,  2,  1,  2, 101,   0, pcol356    ) ! 3:  t
                        call lhe_add_particle( -6,  2,  1,  2,   0, 103, pcol478    ) ! 4:  t~
                        call lhe_add_particle( 24,  2,  3,  0,   0,   0, pcol56     ) ! 5:  W+
                        call lhe_add_particle(  5,  1,  3,  0, 101,   0, pcol(1:4,3)) ! 6:  b
                        call lhe_add_particle(-24,  2,  4,  0,   0,   0, pcol78     ) ! 7:  W-
                        call lhe_add_particle( -5,  1,  4,  0,   0, 103, pcol(1:4,4)) ! 8:  b~
                    end if

                    if (include_qq) then
                        call lhe_add_particle(  1, -1,  0,  0, 101,   0, pcol(1:4,1)) ! 1:  d
                        call lhe_add_particle( -1, -1,  0,  0,   0, 102, pcol(1:4,2)) ! 2:  d~
                        call lhe_add_particle(  6,  2,  1,  2, 101,   0, pcol356    ) ! 3:  t
                        call lhe_add_particle( -6,  2,  1,  2,   0, 102, pcol478    ) ! 4:  t~
                        call lhe_add_particle( 24,  2,  3,  0,   0,   0, pcol56     ) ! 5:  W+
                        call lhe_add_particle(  5,  1,  3,  0, 101,   0, pcol(1:4,3)) ! 6:  b
                        call lhe_add_particle(-24,  2,  4,  0,   0,   0, pcol78     ) ! 7:  W-
                        call lhe_add_particle( -5,  1,  4,  0,   0, 102, pcol(1:4,4)) ! 8:  b~
                    end if

                    if (include_dd) then
                        call lhe_add_particle(  1, -1,  0,  0, 101,   0, pcol(1:4,1)) ! 1:  d
                        call lhe_add_particle( -1, -1,  0,  0,   0, 101, pcol(1:4,2)) ! 2:  d~
                    end if

                    if (include_uu) then
                        call lhe_add_particle(  2, -1,  0,  0, 101,   0, pcol(1:4,1)) ! 1:  u
                        call lhe_add_particle( -2, -1,  0,  0,   0, 101, pcol(1:4,2)) ! 2:  u~
                    end if

                    if (include_dd .or. include_uu) then
                        call lhe_add_particle(  6,  2,  1,  2, 102,   0, pcol356    ) ! 3:  t
                        call lhe_add_particle( -6,  2,  1,  2,   0, 102, pcol478    ) ! 4:  t~
                        call lhe_add_particle( 24,  2,  3,  0,   0,   0, pcol56     ) ! 5:  W+
                        call lhe_add_particle(  5,  1,  3,  0, 102,   0, pcol(1:4,3)) ! 6:  b
                        call lhe_add_particle(-24,  2,  4,  0,   0,   0, pcol78     ) ! 7:  W-
                        call lhe_add_particle( -5,  1,  4,  0,   0, 102, pcol(1:4,4)) ! 8:  b~
                    end if

                    call lhe_add_particle(-i,       1,  5,  0,   0,   0, pcol(1:4,5)) ! 9:  e+,  mu,  ta+
                    call lhe_add_particle( i + 1,   1,  5,  0,   0,   0, pcol(1:4,6)) ! 10: ve,  vm,  vt
                    call lhe_add_particle( j,       1,  7,  0,   0,   0, pcol(1:4,7)) ! 11: e-,  mu-, ta-
                    call lhe_add_particle(-j - 1,   1,  7,  0,   0,   0, pcol(1:4,8)) ! 12: ve~, vm~, vt~

                    call lhe_end_event
                end do
            end do
        end if
    end if
    return
end function dsigma

end module scattering
