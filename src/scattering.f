module scattering

  implicit none

  real :: sigma, sigma_pol(-1:1, -1:1, 20), error_pol(-1:1, -1:1, 20)
  real :: m3, m4, m5, m6, m7, m8, s
  real, parameter :: unit_conv = 0.38937966d9 ! GeV^{-2} to nb (pb?)
  public :: dsigma

contains

function dsigma(x, wgt)

  ! computes the differential cross section for
  !   p p -> f f~
  !   p p -> t t~ -> b b~ W+ W- -> b b~ l+ l- vl vl~

  ! requirements
  !   HELAS subroutines
  !   VEGAS Monte Carlo integration
  !   CT10, CTEQ6 and MRS99 PDF subroutines
  !   RootTuple for filling n-tuples

  ! authors
  !   Declan Millar <declan.millar@cern.ch>
  !   Stefano Moretti

  use configuration
  use modelling
  use integration
  use lhef

  implicit none

  ! differential cross section
  real :: dsigma, ddsigma, x(100), wgt, weight

  ! alphas
  real :: alfas, gs4, gs2, a_s

  ! square matrix elements
  real :: qfduu1, qfduu2, qfddd1, qfddd2, qcdqq, qcdgg
  real :: sgg_tt, sqq_tt, sqq_ff
  real :: sgg_tt_bbeevv, sqq_tt_bbeevv_qcd, sqq_tt_bbeevv
  real :: sgg_bbemuvevm, sqq_bbemuvevm, suu_bbemuvevm, sdd_bbemuvevm
  real :: sgg_bbtatavtvt, sqq_bbtatavtvt, suu_bbtatavtvt, sdd_bbtatavtvt

  ! pdfs
  real :: ctq6pdf, ct14pdf, qq
  real :: fx1(13), fx2(13), x1, x2, xx1, xx2, x1x2(2,2)
  real :: d1, u1, str1, chm1, btm1, glu1
  real :: d2, u2, str2, chm2, btm2, glu2
  real :: dbar1, ubar1, dbar2, ubar2
  real :: dsea1, usea1, dsea2, usea2
  integer :: ix, ixmax

  ! energies
  real :: shat, tau, ecm, ecm_max, ecm_min

  real :: qcm, pcm, qcm2
  real :: pq5, pq52, pq56, pq7, pq78
  real :: rl356, rl478, rl56, rl78, rpl356, rpl478
  real :: rps, rps356, rps478

  ! transfer invariant masses
  real :: q, q2, rq5, rq52, rq56, rq562, rq7, rq72, rq78, rq782

  ! boost
  real :: gamma, v

  ! angles
  real :: phi
  real :: st, st5, st56, st7, st78
  real :: ct, ct5, ct56, ct7, ct78
  real :: sf5, sf56, sf7, sf78
  real :: cf5, cf56, cf7, cf78

  ! temporary top mass and width
  real :: mt, gamt

  ! rambo
  real :: rmass(100), prambo(4,100), wgtr
  integer :: jps

  ! for cuts
  real :: eta, pt

  ! arctan
  real :: at356, at356max, at356min, at56, at56max, at56min
  real :: at478, at478max, at478min, at78, at78max, at78min

  ! 4-momenta
  real :: p(4,8), p356(4), p478(4), q56(4), q78(4), p56(4), p78(4), q5(4), q7(4)
  real :: pcol(4,8), pcol56(4), pcol78(4), pcol356(4), pcol478(4)
  real :: p1(0:3), p2(0:3), p3(0:3), p4(0:3), p5(0:3), p6(0:3), p7(0:3), p8(0:3)

  ! invariant masses
  real :: m356, m356max, m356min, m478, m478max, m478min
  real :: m356_2, m478_2, m56, m56_2, m56max, m56min, m78, m78_2, m78max, m78min

  ! polarised
  real :: qcdpolqq(-1:1, -1:1), qcdpolgg(-1:1, -1:1)
  real :: qfdpoluu1(-1:1, -1:1), qfdpoldd1(-1:1, -1:1), qfdpolbb1(-1:1, -1:1)
  real :: qfdpoluu2(-1:1, -1:1), qfdpoldd2(-1:1, -1:1), qfdpolbb2(-1:1, -1:1)
  real :: weight_pol(-1:1, -1:1, 20), pfx(-1:1, -1:1), pfxtot

  ! internal random number seed
  integer, parameter :: jseed = 987654321

  ! temporary iterators
  integer :: i, j, k

  ! ---

  if (verbose) print*, "--- begin dsigma ---"

  ! store top parameters
  mt = fmass(ffinal)
  gamt = fwidth(ffinal)

  ! limits
  if (ecm_up == 0.d0) then
    ecm_max = sqrts
  else
    ecm_max = ecm_up
  end if
  if (ecm_low == 0.d0) then
    ecm_min = m3 + m4 + m5 + m6 + m7 + m8
  else
    ecm_min = ecm_low
  end if

  if (use_rambo) then
    ecm = x(1) * (ecm_max - ecm_min) + ecm_min
  else 
    ecm = x(2 + 12 * tops_decay) * (ecm_max - ecm_min) + ecm_min
  end if
  shat = ecm * ecm
  tau = shat / s

  ! scale for the pdfs
  if (final_state < 0) then
    qq = ecm
  else
    qq = 2.d0 * mt
  end if
  if (qq == 0.d0) then
    qq = ecm
  end if

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
  dsigma = 0.d0
  if (symmetrise) then
    ixmax = 2
  else
    ixmax = 1
  end if

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

    ! construct hadronic structure functions
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
      if ((x1 <= 1.d-9) .or. (x1 >= 1.d0)) return
      if ((x2 <= 1.d-9) .or. (x2 >= 1.d0)) return
      if ((qq <= 1.3d0) .or. (qq >= 1.d5)) return

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

        phi = 2.d0 * pi * ran(jseed)
        ct = x(1)
        st = sqrt(1.d0 - ct * ct)

        ! magnitude of 3 momentum for products in general two body decay
        qcm2 = ((ecm * ecm - m3 * m3 - m4 * m4)**2 - (2.d0 * m3 * m4)**2) / (4.d0 * ecm * ecm)
        if (qcm2 < 0.d0) return
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
        phi = 2.d0 * pi * ran(jseed)

        m356min = m3 + m5 + m6
        m356max = ecm - m4 - m7 - m8
        if (map_phase_space) then
          at356min = atan((m356min * m356min - mt * mt) / mt / gamt)
          at356max = atan((m356max * m356max - mt * mt) / mt / gamt)
          at356 = x(13) * (at356max - at356min) + at356min
          rl356 = tan(at356) * mt * gamt
          m356_2 = mt * mt + rl356
          if (m356_2 < 0.d0) return
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
          if (m478_2 < 0.d0) return
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
          if (m56_2 < 0.d0) return
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
          if (m78_2 < 0.d0) return
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
        if (q2 < 0.d0) return
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
        if (rq562 < 0.d0) return
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
        if (rq782 < 0.d0) return
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
        if (rq52 < 0.d0) return
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
        if (rq72 < 0.d0) return
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

    ! velocity of tt system in collider frame
    v = (x1 - x2) / (x1 + x2)

    ! gamma factor
    gamma = (x1 + x2) / 2.d0 / sqrt(x1 * x2)

    ! boost initial and final state momenta to the collider frame
    do i = 1, n_final
      pcol(4,i) = gamma * (p(4,i) + v * p(3,i))
      pcol(3,i) = gamma * (p(3,i) + v * p(4,i))
      pcol(2,i) = p(2,i)
      pcol(1,i) = p(1,i)
    end do

    ! fiducial cuts
    if (cut) then
      do i = 3, n_final

        pt = sqrt(pcol(1,i) * pcol(1,i) + pcol(2,i) * pcol(2,i))
        if (pt < 25) then
          dsigma = 0.d0
          return
        end if

        eta = atanh(pcol(3,i) / sqrt(pcol(1,i) * pcol(1,i) + pcol(2,i) * pcol(2,i) + pcol(3,i) * pcol(3,i)))
        if (abs(eta) > 2.5) then
          dsigma = 0.d0
          return
        end if
      end do
    end if

    ! parton CoM 4-momenta
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
      print*, "Check 4-momenta:"
      print*, "p1 = ", p1
      print*, "p2 = ", p2
      print*, "p3 = ", p3
      print*, "p4 = ", p4
      print*, "p5 = ", p5
      print*, "p6 = ", p6
      print*, "p7 = ", p7
      print*, "p8 = ", p8
      print*, ""
      print*, "Check conservation of 4-momentum:"
      print*, "Delta E  = ", p1(0) + p2(0) - p3(0) - p4(0) - p5(0) - p6(0) - p7(0) - p8(0)
      print*, "Delta Px = ", p1(1) + p2(1) - p3(1) - p4(1) - p5(1) - p6(1) - p7(1) - p8(1)
      print*, "Delta Py = ", p1(2) + p2(2) - p3(2) - p4(2) - p5(2) - p6(2) - p7(2) - p8(2)
      print*, "Delta Pz = ", p1(3) + p2(3) - p3(3) - p4(3) - p5(3) - p6(3) - p7(3) - p8(3)
      print*, ""
      print*, "Check invariant masses:"
      print*, "m1 = ", sqrt(abs(p1(0) * p1(0) - p1(1) * p1(1) - p1(2) * p1(2) - p1(3) * p1(3)))
      print*, "m2 = ", sqrt(abs(p2(0) * p2(0) - p2(1) * p2(1) - p2(2) * p2(2) - p2(3) * p2(3)))
      print*, "m3 = ", sqrt(abs(p3(0) * p3(0) - p3(1) * p3(1) - p3(2) * p3(2) - p3(3) * p3(3)))
      print*, "m4 = ", sqrt(abs(p4(0) * p4(0) - p4(1) * p4(1) - p4(2) * p4(2) - p4(3) * p4(3)))
      print*, "m5 = ", sqrt(abs(p5(0) * p5(0) - p5(1) * p5(1) - p5(2) * p5(2) - p5(3) * p5(3)))
      print*, "m6 = ", sqrt(abs(p6(0) * p6(0) - p6(1) * p6(1) - p6(2) * p6(2) - p6(3) * p6(3)))
      print*, "m7 = ", sqrt(abs(p7(0) * p7(0) - p7(1) * p7(1) - p7(2) * p7(2) - p7(3) * p7(3)))
      print*, "m8 = ", sqrt(abs(p8(0) * p8(0) - p8(1) * p8(1) - p8(2) * p8(2) - p8(3) * p8(3)))
      print*, ""
    end if

    if (.not. include_gg .and. .not. include_qq .and. .not. include_uu .and. .not. include_dd) then
      ! set |M| = 1 and skip matrix element calculation
      if (ix == 1) then
        pfxtot = 0.5 / x1
      else if (ix == 2) then
        pfxtot = 0.5 / x2
      end if
      if (final_state <= 0) then
        do i = -1, 1, 2
          do j = -1, 1, 2
            if (ix == 1) then
              pfx(i,j) = 0.5 / x1 / (pfxtot + pfxtot)
            else if (ix == 2) then
              pfx(i,j) = 0.5 / x2 / (pfxtot + pfxtot)
            end if
          end do
        end do
      end if
      go to 666
    end if

    ! Calculate QCD coupling
    gs2 = 4.d0 * pi * alfas(qq, lambdaqcd4, nloops)
    gs4 = gs2 * gs2

    ! initilise
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
        pfx(i,j) = 0.d0
        do k = 1, 20
          weight_pol(i, j, k) = 0.d0
        end do
      end do
    end do

    if (final_state <= 0) then
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
      stop "ERROR! Invalid final state."
    end if

    ! multiple qcd |m|^2 by g_s^4 (madgraph gs is set to one due to scale dependence)
    qcdqq = qcdqq * gs4
    qcdgg = qcdgg * gs4

    pfxtot = 0.d0
    if (final_state <= 0) then
      ! Summing over 2to2 |m|^2 with pdfs of all initial partons
      do i = -1, 1, 2
        do j = -1, 1, 2
                   pfx(i,j) = qcdpolgg(i,j)   * fx1(13) * fx2(13) &
           + (qcdpolqq(i,j) + qfdpoldd1(i,j)) * fx1( 1) * fx2( 7) &
           + (qcdpolqq(i,j) + qfdpoluu1(i,j)) * fx1( 2) * fx2( 8) &
           + (qcdpolqq(i,j) + qfdpoldd1(i,j)) * fx1( 3) * fx2( 9) &
           + (qcdpolqq(i,j) + qfdpoluu1(i,j)) * fx1( 4) * fx2(10) &
           + (qcdpolqq(i,j) + qfdpoldd1(i,j)) * fx1( 5) * fx2(11) &
           + (qcdpolqq(i,j) + qfdpoldd2(i,j)) * fx1( 7) * fx2( 1) &
           + (qcdpolqq(i,j) + qfdpoluu2(i,j)) * fx1( 8) * fx2( 2) &
           + (qcdpolqq(i,j) + qfdpoldd2(i,j)) * fx1( 9) * fx2( 3) &
           + (qcdpolqq(i,j) + qfdpoluu2(i,j)) * fx1(10) * fx2( 4) &
           + (qcdpolqq(i,j) + qfdpoldd2(i,j)) * fx1(11) * fx2( 5)
          if (ix == 1) then
            pfx(i, j) = pfx(i, j) / x1
          else if (ix == 2) then
            pfx(i, j) = pfx(i, j) / x2
          end if
          pfxtot = pfxtot + pfx(i,j)
        end do
      end do
    else if (final_state > 0) then
      ! sum over 2->6 |m|^2 with PDFs of all initial partons"
      pfxtot = fx1( 1) * fx2( 7) * (qcdqq + qfddd1) &
             + fx1( 2) * fx2( 8) * (qcdqq + qfduu1) &
             + fx1( 3) * fx2( 9) * (qcdqq + qfddd1) &
             + fx1( 4) * fx2(10) * (qcdqq + qfduu1) &
             + fx1( 5) * fx2(11) * (qcdqq + qfddd1) &
             + fx1( 7) * fx2( 1) * (qcdqq + qfddd2) &
             + fx1( 8) * fx2( 2) * (qcdqq + qfduu2) &
             + fx1( 9) * fx2( 3) * (qcdqq + qfddd2) &
             + fx1(10) * fx2( 4) * (qcdqq + qfduu2) &
             + fx1(11) * fx2( 5) * (qcdqq + qfddd2) &
             + fx1(13) * fx2(13) * qcdgg
      if (ix == 1) then
        pfxtot = (pfxtot)/x1
      else if (ix  ==  2) then
        pfxtot = (pfxtot)/x2
      end if
    end if

    if (pfxtot == 0.d0) return

    if (final_state < 1) then
      ! weight for distributions
      do i = -1, 1, 2
        do j = -1, 1, 2
          pfx(i,j) = pfx(i, j) / pfxtot
        end do
      end do
    end if

    666 continue

    ! multiply by jacobian from dx1 dx2 -> dx(2) dx(3)
    pfxtot = pfxtot * (1.d0 - tau) * 2.d0 * ecm / s * (ecm_max - ecm_min)

    ddsigma = pfxtot

    ! apply unit converstion
    ddsigma = ddsigma * unit_conv

    ! multiply by phase space factor, azimuthal integration and flux factor
    if (final_state <= 0) then
      if (use_rambo) then
        ddsigma = ddsigma * wgtr
      else
        ! ddsigma = ddsigma * qcm / (2.d0 * pcm) * 2.d0**(4 - 3 * (2)) * 2.d0 * pi
        ddsigma = ddsigma * qcm * pi / ecm / 2
      end if
      ! ddsigma = ddsigma / 2.d0 / ecm / ecm * (2.d0 * pi)**(4 - 3 * (2))
      ddsigma = ddsigma / ecm / ecm / pi / pi / 8

    else if (final_state > 0) then
      if (use_rambo) then
        ddsigma = ddsigma * wgtr
      else
        ! ddsigma = ddsigma * q * rq56 * rq78 * rq5 * rq7/ecm * 256.d0 * 2.d0**(4 - 3 * (6)) * 2.d0 * pi
        ddsigma = ddsigma * q * rq56 * rq78 * rq5 * rq7 * pi / ecm / 32
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

      ! ddsigma = ddsigma / 2.d0 / ecm / ecm * (2.d0 * pi)**(4 - 3 * (6))
      ddsigma = ddsigma / ecm / ecm / pi**14 / 32768
    end if

    ddsigma = ddsigma
    if (symmetrise) ddsigma = ddsigma/ 2

    weight = ddsigma * wgt

    ! Compute polarised event weightings
    if (final_state < 1) then
      do i = -1, 1, 2
        do j = -1, 1, 2
          sigma_pol(i, j, it) = sigma_pol(i, j, it) + ddsigma * wgt * pfx(i, j)
          weight_pol(i, j, it) = ddsigma * wgt * pfx(i, j)
          error_pol(i, j, it) = error_pol(i, j, it) + sigma_pol(i, j, it)**2
        end do
      end do
    end if

    if (ntuple_out) then
      call rootaddint(it, "iteration")

      ! write final particle collider frame momenta to n-tuple
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
        call rootadddouble(weight_pol(-1, -1, it), "weightLL")
        call rootadddouble(weight_pol(-1,  1, it), "weightLR")
        call rootadddouble(weight_pol( 1, -1, it), "weightRL")
        call rootadddouble(weight_pol( 1,  1, it), "weightRR")
      end if

      call rootaddevent(weight)
    end if

    if (lhef_out) then
      pcol56 = pcol(1:4, 5) + pcol(1:4, 6)
      pcol78 = pcol(1:4, 7) + pcol(1:4, 8)
      pcol356 = pcol56 + pcol(1:4, 3)
      pcol478 = pcol78 + pcol(1:4, 4)

      a_s = alfas(zmass, lambdaqcd4, nloops)

      if (include_gg) then
        ! g g -> t t~ -> b W+ b~ W- -> b b~ e+ ve u~ d
        call lhe_add_event(12, 81, weight, qq, a_em, a_s)
        call lhe_add_particle( 21, -1,  0,  0, 101, 102, pcol(1:4,1), 0.0, 9.0) ! g
        call lhe_add_particle( 21, -1,  0,  0, 103, 101, pcol(1:4,2), 0.0, 9.0) ! g
        call lhe_add_particle(  6,  2,  1,  2, 103,   0, pcol356    , 0.0, 9.0) ! t
        call lhe_add_particle( -6,  2,  1,  2,   0, 102, pcol478    , 0.0, 9.0) ! t~
        call lhe_add_particle( 24,  2,  3,  3,   0,   0, pcol56     , 0.0, 9.0) ! W+
        call lhe_add_particle(  5,  1,  3,  3, 103,   0, pcol(1:4,3), 0.0, 9.0) ! b
        call lhe_add_particle(-24,  2,  4,  4,   0,   0, pcol78     , 0.0, 9.0) ! W-
        call lhe_add_particle( -5,  1,  4,  4,   0, 102, pcol(1:4,4), 0.0, 9.0) ! b~
        call lhe_add_particle(-11,  1,  6,  6,   0,   0, pcol(1:4,5), 0.0, 9.0) ! e+
        call lhe_add_particle( 12,  1,  6,  6,   0,   0, pcol(1:4,6), 0.0, 9.0) ! ve
        call lhe_add_particle(  1,  1, 10, 10, 104,   0, pcol(1:4,7), 0.0, 9.0) ! u~
        call lhe_add_particle( -2,  1, 10, 10,   0, 104, pcol(1:4,8), 0.0, 9.0) ! d
      else if (include_qq) then
        ! q q~ -> t t~ -> b  W+ b~ W- -> b b~ e+ ve u~ d
        call lhe_add_event(12, 82, weight, qq, a_em, a_s)
        call lhe_add_particle(  2, -1,  0,  0, 101,   0, pcol(1:4,1), 0.0, 9.0) ! u
        call lhe_add_particle( -2, -1,  0,  0,   0, 102, pcol(1:4,2), 0.0, 9.0) ! u~
        call lhe_add_particle(  6,  2,  1,  2, 101,   0, pcol356    , 0.0, 9.0) ! t
        call lhe_add_particle( -6,  2,  1,  2,   0, 102, pcol478    , 0.0, 9.0) ! t~
        call lhe_add_particle( 24,  2,  3,  0,   0,   0, pcol56     , 0.0, 9.0) ! W+
        call lhe_add_particle(  5,  1,  3,  0, 101,   0, pcol(1:4,3), 0.0, 9.0) ! b
        call lhe_add_particle(-24,  2,  4,  0,   0,   0, pcol78     , 0.0, 9.0) ! W-
        call lhe_add_particle( -5,  1,  4,  0,   0, 102, pcol(1:4,4), 0.0, 9.0) ! b~
        call lhe_add_particle(-11,  1,  5,  0,   0,   0, pcol(1:4,5), 0.0, 9.0) ! e+
        call lhe_add_particle( 12,  1,  5,  0,   0,   0, pcol(1:4,6), 0.0, 9.0) ! ve
        call lhe_add_particle( 11,  1,  7,  0,   0,   0, pcol(1:4,7), 0.0, 9.0) ! e-
        call lhe_add_particle(-12,  1,  7,  0,   0,   0, pcol(1:4,8), 0.0, 9.0) ! ve~
      else if (include_uu) then
        ! u u~ -> t t~ -> b b~ W+ W- -> b b~ e+ ve u~ d
        call lhe_add_event(12, 83, weight, qq, a_em, a_s)
        call lhe_add_particle(  2, -1,  0,  0, 101,   0, pcol(1:4,1), 0.0, 9.0) ! u
        call lhe_add_particle( -2, -1,  0,  0,   0, 101, pcol(1:4,2), 0.0, 9.0) ! u~
        call lhe_add_particle(  6,  2,  1,  2, 102,   0, pcol356    , 0.0, 9.0) ! t
        call lhe_add_particle( -6,  2,  1,  2,   0, 102, pcol478    , 0.0, 9.0) ! t~
        call lhe_add_particle( 24,  2,  3,  3,   0,   0, pcol56     , 0.0, 9.0) ! W+
        call lhe_add_particle(  5,  1,  3,  3, 102,   0, pcol(1:4,3), 0.0, 9.0) ! b
        call lhe_add_particle(-24,  2,  4,  4,   0,   0, pcol78     , 0.0, 9.0) ! W-
        call lhe_add_particle( -5,  1,  4,  4,   0, 102, pcol(1:4,4), 0.0, 9.0) ! b~
        call lhe_add_particle( -1,  1,  6,  6,   0,   0, pcol(1:4,5), 0.0, 9.0) ! e+
        call lhe_add_particle(  2,  1,  6,  6,   0,   0, pcol(1:4,6), 0.0, 9.0) ! ve
        call lhe_add_particle(  1,  1, 10, 10, 103,   0, pcol(1:4,7), 0.0, 9.0) ! u~
        call lhe_add_particle( -2,  1, 10, 10,   0, 103, pcol(1:4,8), 0.0, 9.0) ! d
      else if (include_dd) then
        ! d d~ -> t t~ -> b b~ W+ W- -> b b~ e+ ve u~ d
        call lhe_add_event(12, 84, weight, qq, a_em, a_s)
        call lhe_add_particle(  1, -1,  0,  0, 101,   0, pcol(1:4,1), 0.0, 9.0) ! d
        call lhe_add_particle( -1, -1,  0,  0,   0, 101, pcol(1:4,2), 0.0, 9.0) ! d~
        call lhe_add_particle(  6,  2,  1,  2, 102,   0, pcol356    , 0.0, 9.0) ! t
        call lhe_add_particle( -6,  2,  1,  2,   0, 102, pcol478    , 0.0, 9.0) ! t~
        call lhe_add_particle( 24,  2,  3,  3,   0,   0, pcol56     , 0.0, 9.0) ! W+
        call lhe_add_particle(  5,  1,  3,  3, 102,   0, pcol(1:4,3), 0.0, 9.0) ! b
        call lhe_add_particle(-24,  2,  4,  4,   0,   0, pcol78     , 0.0, 9.0) ! W-
        call lhe_add_particle( -5,  1,  4,  4,   0, 102, pcol(1:4,4), 0.0, 9.0) ! b~
        call lhe_add_particle( -1,  1,  6,  6,   0,   0, pcol(1:4,5), 0.0, 9.0) ! e+
        call lhe_add_particle(  2,  1,  6,  6,   0,   0, pcol(1:4,6), 0.0, 9.0) ! ve
        call lhe_add_particle(  1,  1, 10, 10, 103,   0, pcol(1:4,7), 0.0, 9.0) ! u~
        call lhe_add_particle( -2,  1, 10, 10,   0, 103, pcol(1:4,8), 0.0, 9.0) ! d
      end if

      call lhe_end_event
    end if

    npoints = npoints + 1
    if (verbose) print*, "Event", npoints, "complete."
    dsigma = dsigma + ddsigma
  end do
  if (verbose) print*, "--- end dsigma ---"
  return
end function dsigma

end module scattering
