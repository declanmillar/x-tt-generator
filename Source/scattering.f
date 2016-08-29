module scattering

  implicit none

  real :: sigma, sigma_pol(-1:1,-1:1,20), error_pol(-1:1,-1:1,20)
  real :: m3, m4, m5, m6, m7, m8, s

  real, parameter :: unit_conv = 0.38937966d9 ! GeV^{-2} to nb (pb?)
  public :: dsigma

contains

function dsigma(x, wgt)

  ! Computes the fully differential cross section for

  ! pp -> ff,
  ! pp -> t t~ -> b b~ W+ W- -> b b~ l+ l- vl vl~.

  ! Uses adapted MadGraph functions.
  ! Uses HELAS subroutines.
  ! Uses VEGAS Monte Carlo integration.
  ! Uses CTEQ6 and MRS99 PDF subroutines.
  ! Uses RootTuple for filling ntuples.

  ! Authors: Declan Millar, Stefano Moretti.

  use configuration
  use modelling
  use integration
  use lhef

  implicit none

  ! differential cross section
  real :: dsigma

  ! vegas arguments
  real :: x(100), wgt

  ! external functions
  real :: alfas
  real :: sgg_tt, sqq_tt, sqq_ff
  real :: sgg_tt_bbeevv, sqq_tt_bbeevv_qcd, sqq_tt_bbeevv
  real :: sgg_bbemuvevm, sqq_bbemuvevm, suu_bbemuvevm, sdd_bbemuvevm
  real :: sgg_bbtatavtvt, sqq_bbtatavtvt, suu_bbtatavtvt, sdd_bbtatavtvt
  real :: ctq6pdf, ct14pdf

  real :: ecm, ecm_max, ecm_min
  real :: hist
  real :: gs, gs2
  real :: gcol, qcm, pcm, qcm2
  real :: pt356, pt478, phi356, phi478, ycol356, ycol478
  real :: phi
  real :: pq5, pq52, pq56, pq7, pq78
  real :: qq
  real :: resall
  real :: rl356, rl478, rl56, rl78, rpl356, rpl478
  real :: rps356, rps, rps478
  real :: rq, rq2, rq5, rq52, rq56, rq562, rq7, rq72, rq78, rq782
  real :: sf5, sf56, sf7, sf78
  real :: shat
  real :: st, st5, st56, st7, st78
  real :: tau
  real :: vcol
  real :: a_s
  real :: arg356, arg478
  real :: beta
  real :: cf5, cf56, cf7, cf78
  real :: ct, ct5, ct56, ct7, ct78

  ! parton momentum fraction
  real :: x1, x2, xx, xx1, xx2

  ! structure functions
  real :: d1, d2, dbar1, dbar2, u1, u2, ubar1, ubar2, str1, str2, &
  chm1, chm2, btm1, btm2, glu1, glu2, ggd, dsea1, usea1, usea2, dsea2

  ! temporary dsigma
  real :: ddsigma

  ! temporary top mass and width
  real :: rmt, gamt

  ! rambo
  real :: xmass(100), prambo(4,100), wgtr

  ! cuts
  real :: eta, arg, pt, rpl

  ! arctan
  real :: xx356max, xx356min, xx478max, xx478min, xx56max, xx56min, xx78max, xx78min

  ! iterators
  integer :: i, j, k, ix, nbin, ibin, jbin, imode, lam3, lam4, jps

  ! phase space vectors.
  real :: q356(4), q478(4), q56(4), q78(4), p56(4), p78(4), q5(4), q7(4), q(4,8)
  real :: qcol(4,8), qcol56(4), qcol78(4), qcol356(4), qcol478(4)

  ! 4-momenta
  real :: p1(0:3), p2(0:3), p3(0:3), p4(0:3), p5(0:3), p6(0:3), p7(0:3), p8(0:3)

  ! check phase space
  real :: delta_e, delta_x, delta_y, delta_z

  ! invariant masses
  real :: mass1, mass2, mass3, mass4, mass5, mass6, mass7, mass8
  real :: m356, m356_2, m356max, m356min, m478, m478_2, m478max, m478min
  real :: m56, m56_2, m56max, m56min, m78, m78_2, m78max, m78min

  ! polarised square matrix elements q-qbar
  real :: qcdpolqq(-1:1, -1:1), qcdpolgg(-1:1, -1:1)!, qcdpolbb(-1:1, -1:1)
  real :: qfdpoluu1(-1:1, -1:1), qfdpoldd1(-1:1, -1:1), qfdpolbb1(-1:1, -1:1)
  real :: pfx(-1:1, -1:1)
  real :: pfxtot

  ! polarised square matrix elements qbar-q
  real :: qfdpoluu2(-1:1, -1:1), qfdpoldd2(-1:1, -1:1), qfdpolbb2(-1:1, -1:1)

  ! square matrix elements
  real :: qfduu1 ,qfduu2, qfddd1, qfddd2, qcdqq, qcdgg, qcdbb

  ! weight per polarisation
  real :: weight(-1:1, -1:1, 20)

  ! pdfs
  real :: fx1(13), fx2(13), x1x2(2,2)

  ! internal random number seed
  integer, parameter :: jseed  = 987654321

  double precision :: time1, time2, time3, time4

  ! store top parameters
  rmt = fmass(ffinal)
  gamt = fwidth(ffinal)

  ! limits
  if (ecm_up == 0.d0) then
    ecm_max = collider_energy
  else
    ecm_max = ecm_up
  end if
  if (ecm_low == 0.d0) then
    ecm_min = m3 + m4 + m5 + m6 + m7 + m8
  else
    ecm_min = ecm_low
  end if

  ecm = x((2 + 12*tops_decay)*(1 - use_rambo) + use_rambo)*(ecm_max - ecm_min) + ecm_min
  shat = ecm*ecm
  tau = shat/s

  ! scale for the pdfs
  if (final_state < 0) then
    qq = ecm
  else
    qq = 2.d0*rmt
  end if
  if (qq == 0.d0) then
    ! qq = 0! Setting to Ecm.
    qq = ecm
  end if

  ! x1 and x2 of the partons
  xx1 = x((3 + 12*tops_decay)*(1 - use_rambo) + 2*use_rambo)*(1.d0 - tau) + tau
  xx2 = tau/xx1
  x1x2(1, 1) = xx1
  x1x2(1, 2) = xx2
  x1x2(2, 1) = xx2
  x1x2(2, 2) = xx1

  ! Symmetrise phase space with x1 <-> x2
  dsigma = 0.d0
  do ix = 1, symmetrise + 1
    ddsigma = 0.d0
    x1 = x1x2(ix, 1)
    x2 = x1x2(ix, 2)

    ! initialisation
    ddsigma = 0.d0
    do i = 1, 100
      xmass(i) = 0.d0
      do j = 1, 4
        prambo(j,i) = 0.d0
      end do
    end do
    do i = 1, 4
      do j = 1, 8
        q(i,j) = 0.d0
        qcol(i,j) = 0.d0
      end do
    end do

    ! Construct hadronic structure functions
    if (structure_function <= 4) then
      if ((x1 <= 1.d-6) .or. (x1 >= 1.d0)) return
      if ((x2 <= 1.d-6) .or. (x2 >= 1.d0)) return
      if ((qq <= 1.3d0) .or. (qq >= 1.d4)) return

      ! *x for compatibility with MRS which return xf(x)
      u1 = x1*ctq6pdf(1, x1, qq)
      d1 = x1*ctq6pdf(2, x1, qq)
      ubar1 = x1*ctq6pdf(-1, x1, qq)
      dbar1 = x1*ctq6pdf(-2, x1, qq)
      str1 = x1*ctq6pdf(3, x1, qq)
      chm1 = x1*ctq6pdf(4, x1, qq)
      btm1 = x1*ctq6pdf(5, x1, qq)
      glu1 = x1*ctq6pdf(0, x1, qq)
      u2 = x2*ctq6pdf(1, x2, qq)
      d2 = x2*ctq6pdf(2, x2, qq)
      ubar2 = x2*ctq6pdf(-1, x2, qq)
      dbar2 = x2*ctq6pdf(-2, x2, qq)
      str2 = x2*ctq6pdf(3, x2, qq)
      chm2 = x2*ctq6pdf(4, x2, qq)
      btm2 = x2*ctq6pdf(5, x2, qq)
      glu2 = x2*ctq6pdf(0, x2, qq)

    else if (structure_function > 4 .and. structure_function < 10) then
      imode = structure_function - 4
      if ((x1 <= 1.d-5) .or. (x1 >= 1.d0)) return
      if ((x2 <= 1.d-5) .or. (x2 >= 1.d0)) return
      if ((qq**2 <= 1.25d0) .or. (qq**2 >= 1.d7)) return

      call mrs99(x1, qq, imode, u1, d1, usea1, dsea1, str1, chm1, btm1, glu1)
      call mrs99(x2, qq, imode, u2, d2, usea2, dsea2, str2, chm2, btm2, glu2)

      u1 = u1 + usea1
      d1 = d1 + dsea1
      u2 = u2 + usea2
      d2 = d2 + dsea2
      ubar1 = usea1
      dbar1 = dsea1
      ubar2 = usea2
      dbar2 = dsea2
    else if (structure_function > 9) then
      if ((x1 <= 1.d-9) .or. (x1 >= 1.d0)) return
      if ((x2 <= 1.d-9) .or. (x2 >= 1.d0)) return
      if ((qq <= 1.3d0) .or. (qq >= 1.d5)) return

      ! *x for compatibility with MRS which return xf(x)
      u1 = x1*ct14pdf(1, x1, qq)
      d1 = x1*ct14pdf(2, x1, qq)
      ubar1 = x1*ct14pdf(-1, x1, qq)
      dbar1 = x1*ct14pdf(-2, x1, qq)
      str1 = x1*ct14pdf(3, x1, qq)
      chm1 = x1*ct14pdf(4, x1, qq)
      btm1 = x1*ct14pdf(5, x1, qq)
      glu1 = x1*ct14pdf(0, x1, qq)
      u2 = x2*ct14pdf(1, x2, qq)
      d2 = x2*ct14pdf(2, x2, qq)
      ubar2 = x2*ct14pdf(-1, x2, qq)
      dbar2 = x2*ct14pdf(-2, x2, qq)
      str2 = x2*ct14pdf(3, x2, qq)
      chm2 = x2*ct14pdf(4, x2, qq)
      btm2 = x2*ct14pdf(5, x2, qq)
      glu2 = x2*ct14pdf(0, x2, qq)
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
    fx2(1) = d2*(1 - initial_state) + dbar2*initial_state
    fx2(2) = u2*(1 - initial_state) + ubar2*initial_state
    fx2(3) = str2
    fx2(4) = chm2
    fx2(5) = btm2
    fx2(6) = 0.d0
    fx2(7) = d2*initial_state + dbar2*(1 - initial_state)
    fx2(8) = u2*initial_state + ubar2*(1 - initial_state)
    fx2(9) = fx2(3)
    fx2(10) = fx2(4)
    fx2(11) = fx2(5)
    fx2(12) = fx2(6)
    fx2(13) = glu2
    do i = 1, 13
      fx1(i) = fx1(i)/x1
      fx2(i) = fx2(i)/x2
    end do

    pcm = ecm/2.d0
    q(4,1) = pcm
    q(3,1) = pcm
    q(2,1) = 0.d0
    q(1,1) = 0.d0
    q(4,2) = pcm
    q(3,2) = -pcm
    q(2,2) = 0.d0
    q(1,2) = 0.d0

    if (final_state < 1) then
      if (use_rambo == 0) then
        ! Calculate 2->2 final state momenta in the parton CoM frame manually

        phi = 2.d0*pi*ran(jseed)
        ct = x(1)
        st = sqrt(1.d0 - ct*ct)

        ! magnitude of 3 momentum for products in general two body decay
        qcm2 = ((ecm*ecm - m3*m3 - m4*m4)**2 - (2.d0*m3*m4)**2)/(4.d0*ecm*ecm)
        if (qcm2 < 0.d0) return
        qcm = sqrt(qcm2)

        q(4,3) = sqrt(qcm2 + m3*m3)
        q(3,3) = qcm*ct
        q(2,3) = qcm*st*cos(phi)
        q(1,3) = qcm*st*sin(phi)
        q(4,4) = sqrt(qcm2 + m4*m4)
        q(3,4) = -qcm*ct
        q(2,4) = -qcm*st*cos(phi)
        q(1,4) = -qcm*st*sin(phi)
      else if (use_rambo == 1) then
        ! calculate 2->2 final state momenta in the parton CoM frame using RAMBO
        xmass(1) = m3
        xmass(2) = m4
        jps = 2
        call rambo(seed, jps, ecm, xmass, prambo, wgtr)
        do i = 3, jps + 2
          do j = 1, 4
            q(j, i) = prambo(j, i - 2)
          end do
        end do
      end if

    else if (final_state > 0) then
      if (use_rambo == 0) then
        ! calculate 2->6 final state momenta in the parton CoM frame manually
        phi = 2.d0*pi*ran(jseed)

        m356min = m3 + m5 + m6
        m356max = ecm - m4 - m7 - m8
        if (map_phase_space == 0) then
          m356 = x(13)*(m356max - m356min) + m356min
        else
          ! flatten the integrand around the top propagator
          xx356min = atan(((m356min)**2 - rmt**2)/rmt/gamt)
          xx356max = atan(((m356max)**2 - rmt**2)/rmt/gamt)
          xx = x(13)*(xx356max - xx356min) + xx356min
          rl356 = tan(xx)*rmt*gamt
          m356_2 = (rmt**2 + rl356)
          if (m356_2 < 0.d0) return
          m356 = sqrt(m356_2)
        end if

        m478min = m4 + m7 + m8
        m478max = ecm - m356
        if (map_phase_space == 0) then
          m478 = x(12)*(m478max - m478min) + m478min
        else
          ! flatten the integrand around the anti-top propagator
          xx478min = atan(((m478min)**2 - rmt**2)/rmt/gamt)
          xx478max = atan(((m478max)**2 - rmt**2)/rmt/gamt)
          xx = x(12)*(xx478max - xx478min) + xx478min
          rl478 = tan(xx)*rmt*gamt
          m478_2 = (rmt**2 + rl478)
          if (m478_2 < 0.d0) return
          m478 = sqrt(m478_2)
        end if

        m56min = m5 + m6
        m56max = m356 - m3
        if (map_phase_space == 0) then
          m56 = x(11)*(m56max - m56min) + m56min
        else
          ! flatten the integrand around the W+ propagator
          xx56min = atan(((m56min)**2 - wmass**2)/wmass/wwidth)
          xx56max = atan(((m56max)**2 - wmass**2)/wmass/wwidth)
          xx = x(11)*(xx56max - xx56min) + xx56min
          rl56 = tan(xx)*wmass*wwidth
          m56_2 = (wmass**2 + rl56)
          if (m56_2 < 0.d0) return
          m56 = sqrt(m56_2)
        end if

        m78min = m7 + m8
        m78max = m478 - m4
        if (map_phase_space == 0) then
          m78 = x(10)*(m78max - m78min) + m78min
        else
          ! flatten the integrand around the W- propagator
          xx78min = atan(((m78min)**2 - wmass**2)/wmass/wwidth)
          xx78max = atan(((m78max)**2 - wmass**2)/wmass/wwidth)
          xx = x(10)*(xx78max - xx78min) + xx78min
          rl78 = tan(xx)*wmass*wwidth
          m78_2 = (wmass**2 + rl78)
          if (m78_2 < 0.d0) return
          m78 = sqrt(m78_2)
        end if

        ! assign angles
        ct = x(9)
        st = sqrt(abs(1.d0 - ct*ct))
        ct56 = x(8)
        st56 = sqrt(1.d0 - ct56*ct56)
        ct78 = x(7)
        st78 = sqrt(1.d0 - ct78*ct78)
        ct5 = x(6)
        st5 = sqrt(1.d0 - ct5*ct5)
        ct7 = x(5)
        st7 = sqrt(1.d0 - ct7*ct7)
        cf56 = cos(x(4))
        sf56 = sin(x(4))
        cf78 = cos(x(3))
        sf78 = sin(x(3))
        cf5 = cos(x(2))
        sf5 = sin(x(2))
        cf7 = cos(x(1))
        sf7 = sin(x(1))

        ! two body decay of s-channel mediating boson
        rq2 = ((ecm*ecm - m356*m356 - m478*m478)**2 - (2.d0*m356*m478)**2)/(4.d0*ecm*ecm)
        if (rq2 < 0.d0) return
        rq = sqrt(rq2)

        q356(3) = rq*ct
        q356(2) = rq*st*cos(phi)
        q356(1) = rq*st*sin(phi)
        q356(4) = sqrt(rq2 + m356*m356)

        do i = 1, 3
          q478(i) =  - q356(i)
        end do
        q478(4) = sqrt(rq2 + m478*m478)

        ! two body decay of the top
        rq562 = ((m356*m356 - m3*m3 - m56*m56)**2 - (2.d0*m3*m56)**2)/(4.d0*m356*m356)
        if (rq562 < 0.d0) return
        rq56 = sqrt(rq562)
        q56(3) = rq56*st56*cf56
        q56(2) = rq56*st56*sf56
        q56(1) = rq56*ct56
        q56(4) = sqrt(rq562 + m56*m56)
        pq56 = 0.d0
        do i = 1,3
          pq56 = pq56 + q356(i)*q56(i)
        end do
        p56(4) = (q356(4)*q56(4) + pq56)/m356
        q(4,3) = q356(4) - p56(4)
        do i = 1,3
          p56(i) = q56(i) + q356(i)*(p56(4) + q56(4))/(q356(4) + m356)
          q(i,3) = q356(i) - p56(i)
        end do

        ! two body decay of the anti-top
        rq782 = ((m478*m478 - m4*m4 - m78*m78)**2 - (2.d0*m4*m78)**2)/(4.d0*m478*m478)
        if (rq782 < 0.d0) return
        rq78 = sqrt(rq782)
        q78(3) = rq78*st78*cf78
        q78(2) = rq78*st78*sf78
        q78(1) = rq78*ct78
        q78(4) = sqrt(rq782 + m78*m78)
        pq78 = 0.d0
        do i = 1, 3
          pq78 = pq78 + q478(i)*q78(i)
        end do
        p78(4) = (q478(4)*q78(4) + pq78)/m478
        q(4,4) = q478(4) - p78(4)
        do i = 1, 3
          p78(i) = q78(i) + q478(i)*(p78(4) + q78(4))/(q478(4) + m478)
          q(i,4) = q478(i) - p78(i)
        end do

        ! two body decay of the W+
        rq52 = ((m56*m56 - m5*m5 - m6*m6)**2 - (2.d0*m5*m6)**2)/(4.d0*m56*m56)
        if (rq52 < 0.d0) return
        rq5 = sqrt(rq52)
        q5(3) = rq5*st5*cf5
        q5(2) = rq5*st5*sf5
        q5(1) = rq5*ct5
        q5(4) = sqrt(rq52 + m5*m5)
        pq5 = 0.d0
        do i = 1, 3
          pq5 = pq5 + p56(i)*q5(i)
        end do
        q(4,5) = (p56(4)*q5(4) + pq5)/m56
        q(4,6) = p56(4) - q(4,5)
        do i = 1,3
          q(i,5) = q5(i) + p56(i)*(q(4,5) + q5(4))/(p56(4) + m56)
          q(i,6) = p56(i) - q(i,5)
        end do

        ! two body decay of the W-
        rq72 = ((m78*m78 - m7*m7 - m8*m8)**2 - (2.d0*m7*m8)**2)/(4.d0*m78*m78)
        if (rq72 < 0.d0) return
        rq7 = sqrt(rq72)
        q7(3) = rq7*st7*cf7
        q7(2) = rq7*st7*sf7
        q7(1) = rq7*ct7
        q7(4) = sqrt(rq72 + m7*m7)
        pq7 = 0.d0
        do i = 1, 3
          pq7 = pq7 + p78(i)*q7(i)
        end do
        q(4,7) = (p78(4)*q7(4) + pq7)/m78
        q(4,8) = p78(4) - q(4,7)
        do i = 1, 3
          q(i,7) = q7(i) + p78(i)*(q(4,7) + q7(4))/(p78(4) + m78)
          q(i,8) = p78(i) - q(i,7)
        end do
      else if (use_rambo == 1) then
        ! Calculate 2->6 final state momenta in the parton CoM using RAMBO
        xmass(1) = m3
        xmass(2) = m4
        xmass(3) = m5
        xmass(4) = m6
        xmass(5) = m7
        xmass(6) = m8
        jps = 6
        call rambo(seed, jps, ecm, xmass, prambo, wgtr)
        do i = 3, jps + 2
          do j = 1, 4
            q(j,i) = prambo(j,i-2)
          end do
        end do
      end if
    end if

    ! velocity of tt system in collider frame
    vcol = (x1 - x2)/(x1 + x2)

    ! gamma factor
    gcol = (x1 + x2)/2.d0/sqrt(x1*x2)

    ! boost initial and final state momenta to the collider frame
    do i = 1, n_final
      qcol(4, i) = gcol*(q(4, i) + vcol*q(3, i))
      qcol(3, i) = gcol*(q(3, i) + vcol*q(4, i))
      qcol(2, i) = q(2, i)
      qcol(1, i) = q(1, i)
    end do

    ! fiducial cuts
    if (cut == 1) then
      do i = 3, n_final

        pt = sqrt(q(1,i)*q(1,i) + q(2,i)*q(2,i))
        if (pt < 25) then
          dsigma = 0.d0
          return
        end if

        eta = atanh(q(3,i) / sqrt(q(1,i)*q(1,i) + q(2,i)*q(2,i) + q(3,i)*q(3,i)))
        if (abs(eta) > 2.5) then
          dsigma = 0.d0
          return
        end if
      end do
    end if

    ! parton CoM 4-momenta
    p1(0) = q(4,1)
    p2(0) = q(4,2)
    p3(0) = q(4,3)
    p4(0) = q(4,4)
    p5(0) = q(4,5)
    p6(0) = q(4,6)
    p7(0) = q(4,7)
    p8(0) = q(4,8)
    do i = 1, 3
      p1(i) = q(i,1)
      p2(i) = q(i,2)
      p3(i) = q(i,3)
      p4(i) = q(i,4)
      p5(i) = q(i,5)
      p6(i) = q(i,6)
      p7(i) = q(i,7)
      p8(i) = q(i,8)
    end do

    if (verbose == 1) then
      if (final_state <= 0) then
        print*, "2->2 kinematics"
        print*, "p1  = ", p1
        print*, "p2  = ", p2
        print*, "p3  = ", p3
        print*, "p4  = ", p4
        delta_e = p1(0) + p2(0) - p3(0) - p4(0)
        delta_x = p1(1) + p2(1) - p3(1) - p4(1)
        delta_y = p1(2) + p2(2) - p3(2) - p4(2)
        delta_z = p1(3) + p2(3) - p3(3) - p4(3)
        print*, "delta_E  = ", delta_e
        print*, "delta_Px = ", delta_x
        print*, "delta_Py = ", delta_y
        print*, "delta_Pz = ", delta_z
        mass1 = sqrt(abs(p1(0)**2 - p1(1)**2 - p1(2)**2 - p1(3)**2))
        mass2 = sqrt(abs(p2(0)**2 - p2(1)**2 - p2(2)**2 - p2(3)**2))
        mass3 = sqrt(abs(p3(0)**2 - p3(1)**2 - p3(2)**2 - p3(3)**2))
        mass4 = sqrt(abs(p4(0)**2 - p4(1)**2 - p4(2)**2 - p4(3)**2))
        print*, "m1 = ", mass1
        print*, "m2 = ", mass2
        print*, "m3 = ", mass3, m3
        print*, "m4 = ", mass4, m4

      else if (final_state > 0) then
        print*, "2->6 kinematics:"
        print*, "p1  = ", p1
        print*, "p2  = ", p2
        print*, "p3  = ", p3
        print*, "p4  = ", p4
        print*, "p5  = ", p5
        print*, "p6  = ", p6
        print*, "p7  = ", p7
        print*, "p8  = ", p8
        ! check conservation of 4-momentum.
        delta_e = p1(0) + p2(0) - p3(0) - p4(0) - p5(0) - p6(0) - p7(0) - p8(0)
        delta_x = p1(1) + p2(1) - p3(1) - p4(1) - p5(1) - p6(1) - p7(1) - p8(1)
        delta_y = p1(2) + p2(2) - p3(2) - p4(2) - p5(2) - p6(2) - p7(2) - p8(2)
        delta_z = p1(3) + p2(3) - p3(3) - p4(3) - p5(3) - p6(3) - p7(3) - p8(3)
        print*, "delta_E  = ", delta_e
        print*, "delta_Px = ", delta_x
        print*, "delta_Py = ", delta_y
        print*, "delta_Pz = ", delta_z
        ! check invarient mass.
        mass1 = sqrt(abs(p1(0)**2 - p1(1)**2 - p1(2)**2 - p1(3)**2))
        mass2 = sqrt(abs(p2(0)**2 - p2(1)**2 - p2(2)**2 - p2(3)**2))
        mass3 = sqrt(abs(p3(0)**2 - p3(1)**2 - p3(2)**2 - p3(3)**2))
        mass4 = sqrt(abs(p4(0)**2 - p4(1)**2 - p4(2)**2 - p4(3)**2))
        mass5 = sqrt(abs(p5(0)**2 - p5(1)**2 - p5(2)**2 - p5(3)**2))
        mass6 = sqrt(abs(p6(0)**2 - p6(1)**2 - p6(2)**2 - p6(3)**2))
        mass7 = sqrt(abs(p7(0)**2 - p7(1)**2 - p7(2)**2 - p7(3)**2))
        mass8 = sqrt(abs(p8(0)**2 - p8(1)**2 - p8(2)**2 - p8(3)**2))
        print*, "m1 = ", mass1
        print*, "m2 = ", mass2
        print*, "m3 = ", mass3
        print*, "m4 = ", mass4
        print*, "m5 = ", mass5
        print*, "m6 = ", mass6
        print*, "m7 = ", mass7
        print*, "m8 = ", mass8
      end if
    end if

    if (phase_space_only == 1) then
      ! Set |M|=1 and skip matrix element calculation
      if (ix == 1) then
        pfxtot = 0.5/x1
      else if (ix == 2) then
        pfxtot = 0.5/x2
      end if
      if (final_state <= 0) then
        do lam3 = -1, 1, 2
          do lam4 = -1, 1, 2
            if (ix == 1) then
              pfx(lam3,lam4) = 0.5/x1/(pfxtot + pfxtot)
            else if (ix == 2) then
              pfx(lam3,lam4) = 0.5/x2/(pfxtot + pfxtot)
            end if
          end do
        end do
      end if
      go to 666
    end if

    ! Calculate QCD coupling
    a_s = alfas(qq, lambdaqcd4, nloops)
    gs2 = 4.d0*pi*a_s
    gs = sqrt(gs2)

    ! initilise
    qcdqq = 0.d0
    qcdgg = 0.d0
    qfduu1 = 0.d0
    qfddd1 = 0.d0
    qfduu2 = 0.d0
    qfddd2 = 0.d0
    do lam3 = -1, 1
      do lam4 = -1, 1
        qcdpolqq(lam3,lam4) = 0.d0
        qcdpolgg(lam3,lam4) = 0.d0
        qfdpoluu1(lam3,lam4) = 0.d0
        qfdpoldd1(lam3,lam4) = 0.d0
        qfdpoluu2(lam3,lam4) = 0.d0
        qfdpoldd2(lam3,lam4) = 0.d0
        pfx(lam3,lam4) = 0.d0
        do i = 1, 20
          weight(lam3,lam4,i) = 0.d0
        end do
      end do
    end do

    resall = 0
    if (final_state <= 0) then
      do lam3 = -1, 1, 2
        do lam4 = -1, 1, 2
          if (include_gg == 1) then
            qcdpolgg(lam3,lam4) = sgg_tt(p1, p2, p3, p4, lam3, lam4)*gs**4
          end if
          if (include_qq == 1) then
            qcdpolqq(lam3,lam4) = sqq_tt(3, p1, p2, p3, p4, lam3, lam4)*gs**4
          end if
          if (include_uu == 1) then
            qfdpoluu1(lam3, lam4) = sqq_ff(3, ffinal, p1, p2, p3, p4, lam3, lam4)
            qfdpoluu2(lam3, lam4) = sqq_ff(3, ffinal, p2, p1, p3, p4, lam3, lam4)
          end if
          if (include_dd == 1) then
            qfdpoldd1(lam3, lam4) = sqq_ff(4, ffinal, p1, p2, p3, p4, lam3, lam4)
            qfdpoldd2(lam3, lam4) = sqq_ff(4, ffinal, p2, p1, p3, p4, lam3, lam4)
          end if
          resall = resall &
          + qcdpolgg(lam3,lam4) + qcdpolqq(lam3,lam4) &
          + qfdpoluu1(lam3,lam4) + qfdpoluu2(lam3,lam4) &
          + qfdpoldd1(lam3,lam4) + qfdpoldd2(lam3,lam4)
        end do
      end do

    else if (final_state == 1) then
      if (include_gg == 1) then
        qcdgg = sgg_tt_bbeevv(p1, p2, p3, p4, p5, p7, p6, p8)
      end if
      if (include_qq == 1) then
        qcdqq = sqq_tt_bbeevv_qcd(3, p1, p2, p3, p4, p5, p7, p6, p8)
      end if
      if (include_uu == 1) then
        qfduu1 = sqq_tt_bbeevv(3, 11, p1, p2, p3, p4, p5, p7, p6, p8)
        qfduu2 = sqq_tt_bbeevv(3, 11, p2, p1, p3, p4, p5, p7, p6, p8)
      end if
      if (include_dd == 1) then
        qfddd1 = sqq_tt_bbeevv(4, 11, p1, p2, p3, p4, p5, p7, p6, p8)
        qfddd2 = sqq_tt_bbeevv(4, 11, p2, p1, p3, p4, p5, p7, p6, p8)
      end if
      resall = qcdqq + qcdgg + qfduu1 + qfddd1 + qfduu2 + qfddd2

    ! else if (final_state == 2) then
    !   if (include_gg == 1) then
    !       qcdgg = sgg_bbtatavtvt(p1, p2, p3, p4, p5, p7, p6, p8)
    !   end if
    !   if (include_qq == 1) then
    !       qcdqq = sqq_bbtatavtvt(p1, p2, p3, p4, p5, p7, p6, p8)
    !   end if
    !   if (include_uu == 1) then
    !       qfduu1 = suu_bbtatavtvt(p1, p2, p3, p4, p5, p7, p6, p8)
    !       qfduu2 = suu_bbtatavtvt(p2, p1, p3, p4, p5, p7, p6, p8)
    !   end if
    !   if (include_dd == 1) then
    !       qfddd1 = sdd_bbtatavtvt(p1, p2, p3, p4, p5, p7, p6, p8)
    !       qfddd2 = sdd_bbtatavtvt(p2, p1, p3, p4, p5, p7, p6, p8)
    !   end if
    !   resall = qcdqq + qcdgg + qfduu1 + qfddd1 + qfduu2 + qfddd2

    else if (final_state == 12) then
      if (include_gg == 1) then
          qcdgg = sgg_bbemuvevm(p1, p2, p3, p4, p5, p7, p6, p8)
      end if
      if (include_qq == 1) then
          qcdqq = sqq_bbemuvevm(p1, p2, p3, p4, p5, p7, p6, p8)
      end if
      if (include_uu == 1) then
          qfduu1 = suu_bbemuvevm(p1, p2, p3, p4, p5, p7, p6, p8)
          qfduu2 = suu_bbemuvevm(p2, p1, p3, p4, p5, p7, p6, p8)
      end if
      if (include_dd == 1) then
          qfddd1 = sdd_bbemuvevm(p1, p2, p3, p4, p5, p7, p6, p8)
          qfddd2 = sdd_bbemuvevm(p2, p1, p3, p4, p5, p7, p6, p8)
      end if
      resall = qcdqq + qcdgg + qfduu1 + qfddd1 + qfduu2 + qfddd2
    else
      stop "ERROR! Invalid final state."
    end if

    ! if all |M|^2 contributions zero skip
    if (resall == 0.d0) return

    ! multiple qcd |m|^2 by g_s^4 (madgraph gs is set to one due to scale dependence.)
    qcdqq = qcdqq*gs**4
    qcdgg = qcdgg*gs**4

    pfxtot = 0.d0
    if (final_state <= 0) then
      ! Summing over 2to2 |m|^2 with pdfs of all initial partons
      do lam3 = -1, 1, 2
        do lam4 = -1, 1, 2
                    pfx(lam3,lam4) = qcdpolgg(lam3,lam4) *fx1(13)*fx2(13) &
           + (qcdpolqq(lam3,lam4) + qfdpoldd1(lam3,lam4))*fx1( 1)*fx2( 7) &
           + (qcdpolqq(lam3,lam4) + qfdpoluu1(lam3,lam4))*fx1( 2)*fx2( 8) &
           + (qcdpolqq(lam3,lam4) + qfdpoldd1(lam3,lam4))*fx1( 3)*fx2( 9) &
           + (qcdpolqq(lam3,lam4) + qfdpoluu1(lam3,lam4))*fx1( 4)*fx2(10) &
           + (qcdpolqq(lam3,lam4) + qfdpoldd1(lam3,lam4))*fx1( 5)*fx2(11) &
           + (qcdpolqq(lam3,lam4) + qfdpoldd2(lam3,lam4))*fx1( 7)*fx2( 1) &
           + (qcdpolqq(lam3,lam4) + qfdpoluu2(lam3,lam4))*fx1( 8)*fx2( 2) &
           + (qcdpolqq(lam3,lam4) + qfdpoldd2(lam3,lam4))*fx1( 9)*fx2( 3) &
           + (qcdpolqq(lam3,lam4) + qfdpoluu2(lam3,lam4))*fx1(10)*fx2( 4) &
           + (qcdpolqq(lam3,lam4) + qfdpoldd2(lam3,lam4))*fx1(11)*fx2( 5)
          if (ix == 1) then
            pfx(lam3,lam4) = pfx(lam3,lam4)/x1
          else if (ix == 2) then
            pfx(lam3,lam4) = pfx(lam3,lam4)/x2
          end if
          pfxtot = pfxtot + pfx(lam3,lam4)
        end do
      end do
    else if (final_state > 0) then
      ! sum over 2->6 |m|^2 with PDFs of all initial partons"
      pfxtot = fx1( 1)*fx2( 7)*(qcdqq + qfddd1) &
             + fx1( 2)*fx2( 8)*(qcdqq + qfduu1) &
             + fx1( 3)*fx2( 9)*(qcdqq + qfddd1) &
             + fx1( 4)*fx2(10)*(qcdqq + qfduu1) &
             + fx1( 5)*fx2(11)*(qcdqq + qfddd1) &
             + fx1( 7)*fx2( 1)*(qcdqq + qfddd2) &
             + fx1( 8)*fx2( 2)*(qcdqq + qfduu2) &
             + fx1( 9)*fx2( 3)*(qcdqq + qfddd2) &
             + fx1(10)*fx2( 4)*(qcdqq + qfduu2) &
             + fx1(11)*fx2( 5)*(qcdqq + qfddd2) &
             + fx1(13)*fx2(13)*qcdgg
      if (ix == 1) then
        pfxtot = (pfxtot)/x1
      else if (ix  ==  2) then
        pfxtot = (pfxtot)/x2
      end if
    end if

    if (pfxtot == 0.d0) return

    if (final_state < 1) then
      ! weight for distributions
      do lam3 = -1, 1, 2
        do lam4 = -1, 1, 2
          pfx(lam3,lam4) = pfx(lam3,lam4)/(pfxtot)
        end do
      end do
    end if

    666 continue

    ! multiply by jacobian from dx1 dx2 -> dx(2) dx(3)
    pfxtot = pfxtot*(1.d0 - tau)*2.d0*ecm/s*(ecm_max - ecm_min)

    ddsigma = pfxtot

    ! apply unit converstion
    ddsigma = ddsigma*unit_conv

    ! multiply by phase space volume and flux factor and azimuthal integration
    if (final_state <= 0) then
      if (use_rambo == 0) then
        ! 2-body phase space factor + azimuthal integration
        ddsigma = ddsigma*qcm/(2.d0*pcm)*2.d0**(4 - 3*(2))*2.d0*pi
      else if (use_rambo == 1) then
        ddsigma = ddsigma*wgtr
      end if

      ! flux factor
      ddsigma = ddsigma/2.d0/ecm/ecm*(2.d0*pi)**(4 - 3*(2))

    else if (final_state > 0) then
      if (use_rambo == 0) then

        ! phase space factor
        ddsigma = ddsigma*rq*rq56*rq78*rq5*rq7/ecm*256.d0*2.d0**(4 - 3*(6))*2.d0*pi

        if (map_phase_space == 1) then
          ddsigma = ddsigma*((m356*m356 - rmt*rmt)**2 + rmt**2*gamt**2)*(xx356max - xx356min)/(2.d0*m356)/rmt/gamt
          ddsigma = ddsigma*((m478*m478 - rmt*rmt)**2 + rmt**2*gamt**2)*(xx478max - xx478min)/(2.d0*m478)/rmt/gamt
          ddsigma = ddsigma*((m56*m56 - wmass*wmass)**2 + wmass**2*wwidth**2)*(xx56max - xx56min)/(2.d0*m56)/wmass/wwidth
          ddsigma = ddsigma*((m78*m78 - wmass*wmass)**2 + wmass**2*wwidth**2)*(xx78max - xx78min)/(2.d0*m78)/wmass/wwidth
          ! nwa
          ddsigma = ddsigma*gamt/gamma_t*gamt/gamma_t
        else
          ddsigma = ddsigma*(m356max - m356min)
          ddsigma = ddsigma*(m478max - m478min)
          ddsigma = ddsigma*(m56max - m56min)
          ddsigma = ddsigma*(m78max - m78min)
        end if
      else if (use_rambo == 1) then
        ddsigma = ddsigma*wgtr
      end if

      ! flux factor
      ddsigma = ddsigma/2.d0/ecm/ecm*(2.d0*pi)**(4 - 3*(6))
    end if

    ddsigma = ddsigma/real(symmetrise + 1)

    ! bin
    hist = ddsigma*wgt

    ! Compute polarised event weightings
    if (final_state < 1) then
      do lam3 = -1, 1, 2
        do lam4 = -1, 1, 2
          sigma_pol(lam3,lam4,it) = sigma_pol(lam3,lam4,it) + ddsigma*wgt*pfx(lam3,lam4)
          ! print*, pfx(lam3,lam4), sigma_pol(lam3,lam4,it)
          weight(lam3,lam4,it) = ddsigma*wgt*pfx(lam3,lam4)
          error_pol(lam3,lam4,it) = error_pol(lam3,lam4,it) + sigma_pol(lam3,lam4,it)**2
        end do
      end do
    end if

    if (ntuple_out == 1) then
      call rootaddint(it, "iteration")

      ! Write final particle collider frame momenta to Ntuple
      if (final_state == -1) then
        call rootaddparticle(11,  qcol(1,3), qcol(2,3), qcol(3,3), qcol(4,3))
        call rootaddparticle(-11, qcol(1,4), qcol(2,4), qcol(3,4), qcol(4,4))
      else if (final_state == 0) then
        call rootaddparticle(6,  qcol(1,3), qcol(2,3), qcol(3,3), qcol(4,3))
        call rootaddparticle(-6, qcol(1,4), qcol(2,4), qcol(3,4), qcol(4,4))
      else if (final_state == 1) then
        call rootaddparticle(5,   qcol(1,3), qcol(2,3), qcol(3,3), qcol(4,3))
        call rootaddparticle(-5,  qcol(1,4), qcol(2,4), qcol(3,4), qcol(4,4))
        call rootaddparticle(-11, qcol(1,5), qcol(2,5), qcol(3,5), qcol(4,5))
        call rootaddparticle(12,  qcol(1,6), qcol(2,6), qcol(3,6), qcol(4,6))
        call rootaddparticle(11,  qcol(1,7), qcol(2,7), qcol(3,7), qcol(4,7))
        call rootaddparticle(-12, qcol(1,8), qcol(2,8), qcol(3,8), qcol(4,8))
      end if

      if (final_state < 1) then
        call rootadddouble(weight(-1,-1,it), "weightLL")
        call rootadddouble(weight(-1, 1,it), "weightLR")
        call rootadddouble(weight( 1,-1,it), "weightRL")
        call rootadddouble(weight( 1, 1,it), "weightRR")
      end if

      call rootaddevent(hist)
    end if

    if (lhef_out == 1) then
      qcol56 = qcol(1:4,5) + qcol(1:4,6)
      qcol78 = qcol(1:4,7) + qcol(1:4,8)
      qcol356 = qcol56 + qcol(1:4,3)
      qcol478 = qcol78 + qcol(1:4,4)

      if (include_gg == 1) then
        call lhe_add_event(12, 1, hist, qq, 0.0078125, alfas(zmass, lambdaqcd4, nloops))

        ! g g -> b b~ W+ W- -> b b~ e+ ve u~ d
        call lhe_add_particle(-1,  21,  0,  0, 501, 502, qcol(1:4,1), 9.0, 9.0) ! 1
        call lhe_add_particle(-1,  21,  0,  0, 503, 501, qcol(1:4,2), 9.0, 9.0) ! 2
        call lhe_add_particle( 2,   6,  1,  2, 503,   0, qcol356    , 9.0, 9.0) ! 3
        call lhe_add_particle( 2,  -6,  1,  2,   0, 502, qcol478    , 9.0, 9.0) ! 4
        call lhe_add_particle( 1,   5,  3,  3, 503,   0, qcol(1:4,3), 9.0, 9.0) ! 5
        call lhe_add_particle( 2,  24,  3,  3,   0,   0, qcol56     , 9.0, 9.0) ! 6
        call lhe_add_particle( 1,  -1,  6,  6,   0,   0, qcol(1:4,5), 9.0, 9.0) ! 7
        call lhe_add_particle( 1,   2,  6,  6,   0,   0, qcol(1:4,6), 9.0, 9.0) ! 8
        call lhe_add_particle( 1,  -5,  4,  4,   0, 502, qcol(1:4,4), 9.0, 9.0) ! 9
        call lhe_add_particle( 2, -24,  4,  4,   0,   0, qcol78     , 9.0, 9.0) ! 10
        call lhe_add_particle( 1,   4, 10, 10, 504,   0, qcol(1:4,7), 9.0, 9.0) ! 11
        call lhe_add_particle( 1,  -3, 10, 10,   0, 504, qcol(1:4,8), 9.0, 9.0) ! 12

        ! g g -> b b~ W+ W- -> b b~ e+ ve c~ s
        ! call lhe_add_particle(-1,  21,  0,  0, 501, 502) ! 1
        ! call lhe_add_particle(-1,  21,  0,  0, 503, 501) ! 2
        ! call lhe_add_particle( 2,  -6,  1,  2,   0, 502) ! 3
        ! call lhe_add_particle( 2,   6,  1,  2, 503,   0) ! 4
        ! call lhe_add_particle( 1,  -5,  3,  3,   0, 502) ! 5
        ! call lhe_add_particle( 2, -24,  3,  3,   0,   0) ! 6
        ! call lhe_add_particle( 1,  -1,  6,  6,   0,   0) ! 7
        ! call lhe_add_particle( 1,   2,  6,  6,   0,   0) ! 8
        ! call lhe_add_particle( 1,  -5,  4,  4, 503,   0) ! 9
        ! call lhe_add_particle( 2,  24,  4,  4,   0,   0) ! 10
        ! call lhe_add_particle( 1,  -7, 10, 10, 504,   0) ! 11
        ! call lhe_add_particle( 1,   8, 10, 10,   0, 504) ! 12

        ! g g -> b b~ W+ W- -> b b~ d s~ e- ve~
        ! call lhe_add_particle(-1,  21,  0,  0, 501, 502) ! 1
        ! call lhe_add_particle(-1,  21,  0,  0, 503, 501) ! 2
        ! call lhe_add_particle( 2,  -6,  1,  2,   0, 502) ! 3
        ! call lhe_add_particle( 2,   6,  1,  2, 503,   0) ! 4
        ! call lhe_add_particle( 1,  -5,  3,  3,   0, 502) ! 5
        ! call lhe_add_particle( 2, -24,  3,  3,   0,   0) ! 6
        ! call lhe_add_particle( 1,  -3,  6,  6, 504,   0) ! 7
        ! call lhe_add_particle( 1,   2,  6,  6,   0, 504) ! 8
        ! call lhe_add_particle( 1,  -5,  4,  4, 503,   0) ! 9
        ! call lhe_add_particle( 2,  24,  4,  4,   0,   0) ! 10
        ! call lhe_add_particle( 1,  -3, 10, 10,   0,   0) ! 11
        ! call lhe_add_particle( 1,   4, 10, 10,   0,   0) ! 12

        ! g g -> b b~ W+ W- -> b b~ d s~ e- ve~
        ! call lhe_add_particle(-1,  21,  0,  0, 501, 502) ! 1
        ! call lhe_add_particle(-1,  21,  0,  0, 503, 501) ! 2
        ! call lhe_add_particle( 2,  -6,  1,  2,   0, 502) ! 3
        ! call lhe_add_particle( 2,   6,  1,  2, 503,   0) ! 4
        ! call lhe_add_particle( 1,  -5,  3,  3,   0, 502) ! 5
        ! call lhe_add_particle( 2, -24,  3,  3,   0,   0) ! 6
        ! call lhe_add_particle( 1,  -1,  6,  6, 504,   0) ! 7
        ! call lhe_add_particle( 1,   2,  6,  6,   0, 504) ! 8
        ! call lhe_add_particle( 1,  -5,  4,  4, 503,   0) ! 9
        ! call lhe_add_particle( 2,  24,  4,  4,   0,   0) ! 10
        ! call lhe_add_particle( 1,  -7, 10, 10, 505,   0) ! 11
        ! call lhe_add_particle( 1,   8, 10, 10,   0, 505) ! 12

        ! g g -> b b~ W+ W- -> b b~ e+ ve u~ d
        ! call lhe_add_particle(-1,  21,  0,  0, 501, 502) ! 1
        ! call lhe_add_particle(-1,  21,  0,  0, 503, 501) ! 2
        ! call lhe_add_particle( 2,  -6,  1,  2,   0, 502) ! 3
        ! call lhe_add_particle( 2,   6,  1,  2, 503,   0) ! 4
        ! call lhe_add_particle( 1,  -5,  3,  3,   0, 502) ! 5
        ! call lhe_add_particle( 2, -24,  3,  3,   0,   0) ! 6
        ! call lhe_add_particle( 1,  -1,  6,  6,   0,   0) ! 7
        ! call lhe_add_particle( 1,   2,  6,  6,   0,   0) ! 8
        ! call lhe_add_particle( 1,  -5,  4,  4, 503,   0) ! 9
        ! call lhe_add_particle( 2,  24,  4,  4,   0,   0) ! 10
        ! call lhe_add_particle( 1,   4, 10, 10, 504,   0) ! 11
        ! call lhe_add_particle( 1,  -3, 10, 10,   0, 504) ! 12

        ! g g -> b b~ W+ W- -> b b~ e+ ve c~ s
        ! call lhe_add_particle(-1,  21,  0,  0, 501, 502) ! 1
        ! call lhe_add_particle(-1,  21,  0,  0, 503, 501) ! 2
        ! call lhe_add_particle( 2,  -6,  1,  2,   0, 502) ! 3
        ! call lhe_add_particle( 2,   6,  1,  2, 503,   0) ! 4
        ! call lhe_add_particle( 1,  -5,  3,  3,   0, 502) ! 5
        ! call lhe_add_particle( 2, -24,  3,  3,   0,   0) ! 6
        ! call lhe_add_particle( 1,  -1,  6,  6,   0,   0) ! 7
        ! call lhe_add_particle( 1,   2,  6,  6,   0,   0) ! 8
        ! call lhe_add_particle( 1,  -5,  4,  4, 503,   0) ! 9
        ! call lhe_add_particle( 2,  24,  4,  4,   0,   0) ! 10
        ! call lhe_add_particle( 1,  -7, 10, 10, 504,   0) ! 11
        ! call lhe_add_particle( 1,   8, 10, 10,   0, 504) ! 12

        ! g g -> b b~ W+ W- -> b b~ d s~ e- ve~
        ! call lhe_add_particle(-1,  21,  0,  0, 501, 502) ! 1
        ! call lhe_add_particle(-1,  21,  0,  0, 503, 501) ! 2
        ! call lhe_add_particle( 2,  -6,  1,  2,   0, 502) ! 3
        ! call lhe_add_particle( 2,   6,  1,  2, 503,   0) ! 4
        ! call lhe_add_particle( 1,  -5,  3,  3,   0, 502) ! 5
        ! call lhe_add_particle( 2, -24,  3,  3,   0,   0) ! 6
        ! call lhe_add_particle( 1,  -3,  6,  6, 504,   0) ! 7
        ! call lhe_add_particle( 1,   2,  6,  6,   0, 504) ! 8
        ! call lhe_add_particle( 1,  -5,  4,  4, 503,   0) ! 9
        ! call lhe_add_particle( 2,  24,  4,  4,   0,   0) ! 10
        ! call lhe_add_particle( 1,  -3, 10, 10,   0,   0) ! 11
        ! call lhe_add_particle( 1,   4, 10, 10,   0,   0) ! 12

        ! g g -> b b~ W+ W- -> b b~ d s~ e- ve~
        ! call lhe_add_particle(-1,  21,  0,  0, 501, 502) ! 1
        ! call lhe_add_particle(-1,  21,  0,  0, 503, 501) ! 2
        ! call lhe_add_particle( 2,  -6,  1,  2,   0, 502) ! 3
        ! call lhe_add_particle( 2,   6,  1,  2, 503,   0) ! 4
        ! call lhe_add_particle( 1,  -5,  3,  3,   0, 502) ! 5
        ! call lhe_add_particle( 2, -24,  3,  3,   0,   0) ! 6
        ! call lhe_add_particle( 1,  -1,  6,  6, 504,   0) ! 7
        ! call lhe_add_particle( 1,   2,  6,  6,   0, 504) ! 8
        ! call lhe_add_particle( 1,  -5,  4,  4, 503,   0) ! 9
        ! call lhe_add_particle( 2,  24,  4,  4,   0,   0) ! 10
        ! call lhe_add_particle( 1,  -7, 10, 10, 505,   0) ! 11
        ! call lhe_add_particle( 1,   8, 10, 10,   0, 505) ! 12
      end if

      if (include_qq == 1) then
        ! q q -> b b~ W+ W- -> b b~ e+ ve u~ d
        call lhe_add_particle(-1,   3,  0,  0, 501,   0, qcol(1:4,1), 9.0, 9.0) ! 1
        call lhe_add_particle(-1,  -3,  0,  0,   0, 502, qcol(1:4,2), 9.0, 9.0) ! 2
        call lhe_add_particle( 2,   6,  1,  2, 501,   0, qcol356    , 9.0, 9.0) ! 3
        call lhe_add_particle( 2,  -6,  1,  2,   0, 502, qcol478    , 9.0, 9.0) ! 4
        call lhe_add_particle( 1,   5,  3,  3, 501,   0, qcol(1:4,3), 9.0, 9.0) ! 5
        call lhe_add_particle( 2,  24,  3,  3,   0,   0, qcol56     , 9.0, 9.0) ! 6
        call lhe_add_particle( 1,  -1,  6,  6,   0,   0, qcol(1:4,5), 9.0, 9.0) ! 7
        call lhe_add_particle( 1,   2,  6,  6,   0,   0, qcol(1:4,6), 9.0, 9.0) ! 8
        call lhe_add_particle( 1,  -5,  4,  4,   0, 502, qcol(1:4,4), 9.0, 9.0) ! 9
        call lhe_add_particle( 2, -24,  4,  4,   0,   0, qcol78     , 9.0, 9.0) ! 10
        call lhe_add_particle( 1,   4, 10, 10, 504,   0, qcol(1:4,7), 9.0, 9.0) ! 11
        call lhe_add_particle( 1,  -3, 10, 10,   0, 504, qcol(1:4,8), 9.0, 9.0) ! 12
      end if

      if (include_uu == 1) then
        ! u u -> b b~ W+ W- -> b b~ e+ ve u~ d
        call lhe_add_particle(-1,   3,  0,  0, 501,   0, qcol(1:4,1), 9.0, 9.0) ! 1
        call lhe_add_particle(-1,  -3,  0,  0,   0, 501, qcol(1:4,2), 9.0, 9.0) ! 2
        call lhe_add_particle( 2,   6,  1,  2, 502,   0, qcol356    , 9.0, 9.0) ! 3
        call lhe_add_particle( 2,  -6,  1,  2,   0, 502, qcol478    , 9.0, 9.0) ! 4
        call lhe_add_particle( 1,   5,  3,  3, 502,   0, qcol(1:4,3), 9.0, 9.0) ! 5
        call lhe_add_particle( 2,  24,  3,  3,   0,   0, qcol56     , 9.0, 9.0) ! 6
        call lhe_add_particle( 1,  -1,  6,  6,   0,   0, qcol(1:4,5), 9.0, 9.0) ! 7
        call lhe_add_particle( 1,   2,  6,  6,   0,   0, qcol(1:4,6), 9.0, 9.0) ! 8
        call lhe_add_particle( 1,  -5,  4,  4,   0, 502, qcol(1:4,4), 9.0, 9.0) ! 9
        call lhe_add_particle( 2, -24,  4,  4,   0,   0, qcol78     , 9.0, 9.0) ! 10
        call lhe_add_particle( 1,   4, 10, 10, 503,   0, qcol(1:4,7), 9.0, 9.0) ! 11
        call lhe_add_particle( 1,  -3, 10, 10,   0, 503, qcol(1:4,8), 9.0, 9.0) ! 12
      end if

      if (include_dd == 1) then
        ! d d -> b b~ W+ W- -> b b~ e+ ve u~ d
        call lhe_add_particle(-1,   4,  0,  0, 501,   0, qcol(1:4,1), 9.0, 9.0) ! 1
        call lhe_add_particle(-1,  -4,  0,  0,   0, 501, qcol(1:4,2), 9.0, 9.0) ! 2
        call lhe_add_particle( 2,   6,  1,  2, 502,   0, qcol356    , 9.0, 9.0) ! 3
        call lhe_add_particle( 2,  -6,  1,  2,   0, 502, qcol478    , 9.0, 9.0) ! 4
        call lhe_add_particle( 1,   5,  3,  3, 502,   0, qcol(1:4,3), 9.0, 9.0) ! 5
        call lhe_add_particle( 2,  24,  3,  3,   0,   0, qcol56     , 9.0, 9.0) ! 6
        call lhe_add_particle( 1,  -1,  6,  6,   0,   0, qcol(1:4,5), 9.0, 9.0) ! 7
        call lhe_add_particle( 1,   2,  6,  6,   0,   0, qcol(1:4,6), 9.0, 9.0) ! 8
        call lhe_add_particle( 1,  -5,  4,  4,   0, 502, qcol(1:4,4), 9.0, 9.0) ! 9
        call lhe_add_particle( 2, -24,  4,  4,   0,   0, qcol78     , 9.0, 9.0) ! 10
        call lhe_add_particle( 1,   4, 10, 10, 503,   0, qcol(1:4,7), 9.0, 9.0) ! 11
        call lhe_add_particle( 1,  -3, 10, 10,   0, 503, qcol(1:4,8), 9.0, 9.0) ! 12
      end if

      ! g g -> b b~ W+ W-
      ! call lhe_add_particle(-1,  21, 0, 0, 501, 502)
      ! call lhe_add_particle(-1,  21, 0, 0, 503, 501)
      ! call lhe_add_particle( 2,  -6, 1, 2,   0, 502)
      ! call lhe_add_particle( 2,   6, 1, 2, 503,   0)
      ! call lhe_add_particle( 1,  -5, 3, 3,   0, 502)
      ! call lhe_add_particle( 1, -24, 3, 3,   0,   0)
      ! call lhe_add_particle( 1,   5, 4, 4, 503,   0)
      ! call lhe_add_particle( 1,  24, 4, 4,   0,   0)
    end if

    npoints = npoints + 1
    if (verbose == 1) print*, "Event", npoints, "complete."
    dsigma = dsigma + ddsigma
  end do
  return
end function dsigma

end module scattering
