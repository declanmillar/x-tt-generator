function dsigma(x,wgt)

  ! Computes the fully differential cross section for
  ! * pp -> tt,
  ! * pp -> tt -> bw^+bbarw^- -> bbbare^+nue^-nubar.
  ! Calculates matrix element using helas subroutines; scales by PDFs, phase
  ! space volume factor and populates root ntuple with final state collider momenta.
  ! authors: declan millar, stefano moretti

  use mathematics, only: pi
  use configuration
  use modelling
  use scattering
  use kinematics
  use integration
  
  implicit none

  real :: x(100), wgt

  ! external functions
  real :: dsigma, alfas, sqqff_qcd, sggff_qcd, sqqff_ewp, sqqbbffff_qcd, sggbbffff_qcd, sqqbbffff_ewp, ctq6pdf

  real :: ecm, ecm_max, ecm_min, pcm, qcm2
  real :: hist, hist1, hist2
  real :: gs, gs2
  real :: gcol, qcm
  real :: pt356, pt478, phi356, phi478, ycol356, ycol478
  real :: phit
  real :: pq5, pq52, pq56, pq7, pq78
  real :: qq
  real :: qqd1, qqd2
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

  ! temporary dsigmas
  real :: ffxn, fffxn, ffffxn, fffffxn

  ! temporary top mass and width
  real :: rmt, gamt

  ! rambo
  real :: xmass(100), prambo(4,100), wgtr

  ! arctan
  real :: xx356max, xx356min, xx478max, xx478min, xx56max, xx56min, xx78max, xx78min

  ! iterators
  integer :: i, j, k, jx, ix, nbin, ibin, jbin, imode, lam3, lam4, jps, i5, i7

  ! phase space vectors.
  real :: q356(4), q478(4), q56(4), q78(4), p56(4), p78(4), q5(4), q7(4), q(4,8), qcol(4,8)

  ! 4-momenta
  real :: p1(0:3), p2(0:3), p3(0:3), p4(0:3), p5(0:3), p6(0:3), p7(0:3), p8(0:3)

  ! check phase space
  real :: delta_e, delta_x, delta_y, delta_z

  ! invarient masses
  real :: mass1, mass2, mass3, mass4, mass5, mass6, mass7, mass8
  real :: m356, m356_2, m356max, m356min, m478, m478_2, m478max, m478min
  real :: m56, m56_2, m56max, m56min, m78, m78_2, m78max, m78min  

  ! polarised square matrix elements q-qbar
  real :: qcdpolqq(-1:1, -1:1), qcdpolbb(-1:1, -1:1), qcdpolgg(-1:1, -1:1)
  real :: ewzpoluu1(-1:1, -1:1), ewzpoldd1(-1:1, -1:1), ewzpolbb1(-1:1, -1:1)
  real :: pfx(-1:1, -1:1)
  real :: pfxtot

  ! polarised square matrix elements qbar-q
  real :: ewzpoluu2(-1:1, -1:1), ewzpoldd2(-1:1, -1:1), ewzpolbb2(-1:1, -1:1)

  ! square matrix elements
  real :: ewzuu1 ,ewzuu2, ewzdd1, ewzdd2, ewzbb1, ewzbb2, qcdqq, qcdgg, qcdbb

  ! weight per polarisation
  real :: weight(-1:1, -1:1, 20)
  real :: weightLL = -999, weightLR = -999, weightRL = -999, weightRR = -999

  ! pdfs
  real :: fx1(13), fx2(13), x1x2(2,2)

  ! internal random number seed
  integer, parameter :: jseed  = 987654321

  ! --- Method ---

  call debug("Start of fxn")

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

  ! x1 and x2 of the partons
  xx1 = x((3 + 12*tops_decay)*(1 - use_rambo) + 2*use_rambo) * (1.d0 - tau) + tau
  xx2 = tau/xx1
  x1x2(1, 1) = xx1
  x1x2(1, 2) = xx2
  x1x2(2, 1) = xx2
  x1x2(2, 2) = xx1

  ! loop over x1 and x2  
  dsigma = 0.d0
  do ix = 1, ixmax
    ffxn = 0.d0
    x1 = x1x2(ix, 1)
    x2 = x1x2(ix, 2)

    do jx = 1, jxmax ! loop over costheta_cm
      ffffxn = 0 
      do i5 = 1, i5max ! loop over costheta5
        fffffxn = 0
        do i7 = 1, i7max ! loop over costheta7
          if (verbose == 1) print*, "Generating event ", npoints + 1, ", x it:", ix, ", c it: ", jx, ", i5: ", i5, ",  i7: ", i7

          ! initialisation
          fffxn = 0.d0
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

          ! scale for the pdfs
          if (final_state < 0) then
            qq = rm_z
          else
            qq = 2.d0*rmt
          end if
          if (qq == 0.d0) then
            call debug("qq = 0! Setting to Z mass.")
            qq = rm_z
          end if


          call debug("Constructing hadronic structure functions...")
          if (structure_function <= 4) then
            if ((x1 <= 1.d-6) .or. (x1 >= 1.d0)) then
              fffxn = 0.d0
              call debug("x1 out of range. Setting fxn = 0 and Skipping.")
              go to 999
            end if
            if ((x2 <= 1.d-6) .or. (x2 >= 1.d0)) then
              fffxn = 0.d0
              call debug("x2 out of range. Setting fxn = 0 and Skipping.")
              go to 999
            end if
            if ((qq <= 1.3d0) .or. (qq >= 1.d4)) then
              fffxn = 0.d0
              call debug("qq out of range. Setting fxn = 0 and Skipping.")
              go to 999
            end if

            ! cteq pdfs are multiplied by x for compatibility with mrs, which return xf(x).

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

          else if (structure_function == 5) then
            imode = 1
            if ((x1 <= 1.d-5) .or. (x1 >= 1.d0)) then
              fffxn = 0.d0
              call debug("fxn = 0. Skipping.")
              go to 999
            end if
            if ((x2 <= 1.d-5) .or. (x2 >= 1.d0)) then
              fffxn = 0.d0
              call debug("fxn = 0. Skipping.")
              go to 999
            end if
            if ((qq**2 <= 1.25d0) .or. (qq**2 >= 1.d7)) then
              fffxn = 0.d0
              call debug("fxn = 0. Skipping.")
              go to 999
            end if
            call mrs99(x1, qq, imode, u1, d1, usea1, dsea1, str1, chm1, btm1, glu1)
            call mrs99(x2, qq, imode, u2, d2, usea2, dsea2, str2, chm2, btm2, glu2)
          else if (structure_function == 6) then
            imode = 2
            if ((x1 <= 1.d-5) .or. (x1 >= 1.d0)) then
              fffxn = 0.d0
              call debug("fxn = 0. Skipping.")
              go to 999
            end if
            if ((x2 <= 1.d-5) .or. (x2 >= 1.d0)) then
              fffxn = 0.d0
              call debug("fxn = 0. Skipping.")
              go to 999
            end if
            if ((qq**2 <= 1.25d0) .or. (qq**2 >= 1.d7)) then
              fffxn = 0.d0
              call debug("fxn = 0. Skipping.")
              go to 999
            end if
            call mrs99(x1, qq, imode, u1, d1, usea1, dsea1, str1, chm1, btm1, glu1)
            call mrs99(x2, qq, imode, u2, d2, usea2, dsea2, str2, chm2, btm2, glu2)
          else if (structure_function == 7) then
            imode = 3
            if ((x1 <= 1.d-5) .or. (x1 >= 1.d0)) then
              fffxn = 0.d0
              call debug("fxn = 0. Skipping.")
              go to 999
            end if
            if ((x2 <= 1.d-5) .or. (x2 >= 1.d0)) then
              fffxn = 0.d0
              call debug("fxn = 0. Skipping.")
              go to 999
            end if
            if ((qq**2 <= 1.25d0) .or. (qq**2 >= 1.d7)) then
              fffxn = 0.d0
              call debug("fxn = 0. Skipping.")
              go to 999
            end if
            call mrs99(x1, qq, imode, u1, d1, usea1, dsea1, str1, chm1, btm1, glu1)
            call mrs99(x2, qq, imode, u2, d2, usea2, dsea2, str2, chm2, btm2, glu2)
          else if (structure_function == 8) then
            imode = 4
            if ((x1 <= 1.d-5) .or. (x1 >= 1.d0)) then
              fffxn = 0.d0
              call debug("fxn = 0. Skipping.")
              go to 999
            end if
            if ((x2 <= 1.d-5) .or. (x2 >= 1.d0)) then
              fffxn = 0.d0
              call debug("fxn = 0. Skipping.")
              go to 999
            end if
            if ((qq**2 <= 1.25d0) .or. (qq**2 >= 1.d7)) then
              fffxn = 0.d0
              call debug("fxn = 0. Skipping.")
              go to 999
            end if
            call mrs99(x1, qq, imode, u1, d1, usea1, dsea1, str1, chm1, btm1, glu1)
            call mrs99(x2, qq, imode, u2, d2, usea2, dsea2, str2, chm2, btm2, glu2)
          else if (structure_function == 9) then
            imode = 5
            if ((x1 <= 1.d-5) .or. (x1 >= 1.d0)) then
              fffxn = 0.d0
              call debug("fxn = 0. Skipping.")
              go to 999
            end if
            if ((x2 <= 1.d-5) .or. (x2 >= 1.d0)) then
              fffxn = 0.d0
              call debug("fxn = 0. Skipping.")
              go to 999
            end if
            if ((qq**2 <= 1.25d0) .or. (qq**2 >= 1.d7)) then
              fffxn = 0.d0
              call debug("fxn = 0. Skipping.")
              go to 999
            end if
            call mrs99(x1, qq, imode, u1, d1, usea1, dsea1, str1, chm1, btm1, glu1)
            call mrs99(x2, qq, imode, u2, d2, usea2, dsea2, str2, chm2, btm2, glu2)
          end if
          call debug("...complete.")
          if (structure_function > 4) then
            u1 = u1 + usea1
            d1 = d1 + dsea1
            u2 = u2 + usea2
            d2 = d2 + dsea2
            ubar1 = usea1
            dbar1 = dsea1
            ubar2 = usea2
            dbar2 = dsea2
          end if

          call debug("Constructing PDFs...")

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
          do i = 1, 13
            fx1(i) = fx1(i)/x1
          end do
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
            fx2(i) = fx2(i)/x2
          end do
          call debug("...complete.")

          call debug("Creating initial (massless) parton momenta.")
          pcm = ecm/2.d0
          q(4,1) = pcm
          q(3,1) = pcm
          q(2,1) = 0.d0
          q(1,1) = 0.d0
          q(4,2) = pcm
          q(3,2) = -pcm
          q(2,2) = 0.d0
          q(1,2) = 0.d0
          call debug("...complete.")

          if (final_state <= 0) then
            if (use_rambo == 0) then
              call debug("Calculating 2to2 final state momenta in the parton CoM frame manually...")
              ! give vegas assigned values
              phit = 2.d0*pi*ran(jseed)
              if (jx == 1) then
                ct = x(1)
              else if (jx == 2) then
                ct = -x(1)
              else
                print*, "Error: invalid jx."
              end if
              st = sqrt(1.d0 - ct*ct)

              ! magnitude of 3 momentum for products in general two body decay
              qcm2 = ((ecm*ecm - m3*m3 - m4*m4)**2 - (2.d0*m3*m4)**2)/(4.d0*ecm*ecm)
              if (qcm2 < 0.d0) then
                fffxn = 0.d0
                call debug("fxn = 0. Skipping.")
                go to 999
              else
                qcm = sqrt(qcm2)
              endif

              q(4,3) = sqrt(qcm2 + m3*m3)
              q(3,3) = qcm*ct
              q(2,3) = qcm*st*cos(phit)
              q(1,3) = qcm*st*sin(phit)
              q(4,4) = sqrt(qcm2 + m4*m4)
              q(3,4) = -qcm*ct
              q(2,4) = -qcm*st*cos(phit)
              q(1,4) = -qcm*st*sin(phit)
            else if (use_rambo == 1) then
              call debug("Calculating 2to2 final state momenta in the parton CoM frame using RAMBO...")
              xmass(1) = m3
              xmass(2) = m4
              jps = 2
              call rambo(seed,jps,ecm,xmass,prambo,wgtR)
              do i = 3, jps + 2
                do j = 1, 4
                  q(j,i) = prambo(j, i-2)
                end do
              end do
            end if
            ! initialise unrequired array elements
            do i = 1, 4
              q(i,5) = 0.d0
              q(i,6) = 0.d0
              q(i,7) = 0.d0
              q(i,8) = 0.d0
            end do
            call debug("...complete.")

          else if (final_state > 0) then
            if (use_rambo == 0) then
              call debug("Calculating 2to6 final state momenta in the parton CoM frame manually...")
              phit = 2.d0*pi*ran(jseed)

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
                if (m356_2 < 0.d0) then
                  fffxn = 0.d0
                  call debug("fxn = 0. Skipping.")
                  go to 999
                else
                  m356 = sqrt(m356_2)
                endif
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
                if (m478_2 < 0.d0) then
                  fffxn = 0.d0
                  call debug("fxn = 0. Skipping.")
                  go to 999
                else
                  m478 = sqrt(m478_2)
                endif
              end if

              m56min = m5 + m6
              m56max = m356 - m3
              if (map_phase_space == 0) then
                m56 = x(11)*(m56max - m56min) + m56min
              else
                ! flatten the integrand around the W+ propagator
                xx56min = atan(((m56min)**2 - rm_w**2)/rm_w/gamma_w)
                xx56max = atan(((m56max)**2 - rm_w**2)/rm_w/gamma_w)
                xx = x(11)*(xx56max - xx56min) + xx56min
                rl56 = tan(xx)*rm_w*gamma_w
                m56_2 = (rm_w**2 + rl56)
                if (m56_2 < 0.d0) then
                  fffxn = 0.d0
                  call debug("fxn = 0. Skipping.")
                  go to 999
                else
                  m56 = sqrt(m56_2)
                endif
              end if

              m78min = m7 + m8
              m78max = m478 - m4
              if (map_phase_space == 0) then
                m78 = x(10)*(m78max - m78min) + m78min
              else
                ! flatten the integrand around the W- propagator
                xx78min = atan(((m78min)**2 - rm_w**2)/rm_w/gamma_w)
                xx78max = atan(((m78max)**2 - rm_w**2)/rm_w/gamma_w)
                xx = x(10)*(xx78max - xx78min) + xx78min
                rl78 = tan(xx)*rm_w*gamma_w
                m78_2 = (rm_w**2 + rl78)
                if (m78_2 < 0.d0) then
                  fffxn = 0.d0
                  call debug("fxn = 0. Skipping.")
                  go to 999
                else
                  m78 = sqrt(m78_2)
                endif
              end if

              if (jx == 1) ct = x(9)
              if (jx == 2) ct = -x(9)

              ! assign angles
              st = sqrt(abs(1.d0 - ct*ct))
              ct56 = x(8)
              st56 = sqrt(1.d0 - ct56*ct56)
              ct78 = x(7)
              st78 = sqrt(1.d0 - ct78*ct78)
              if (i5 == 1) ct5 = x(6)
              if (i5 == 2) ct5 = -x(6)
              st5 = sqrt(1.d0 - ct5*ct5)
              if (i7 == 1) ct7 = x(5)
              if (i7 == 2) ct7 = -x(5)
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
              if (rq2 < 0.d0) then
                fffxn = 0.d0
                call debug("fxn = 0. Skipping.")
                go to 999
              else
                rq = sqrt(rq2)
              endif

              q356(3) = rq*ct
              q356(2) = rq*st*cos(phit)
              q356(1) = rq*st*sin(phit)
              q356(4) = sqrt(rq2 + m356*m356)

              do i = 1, 3
                q478(i) =  - q356(i)
              end do
              q478(4) = sqrt(rq2 + m478*m478)

              ! two body decay of the top
              rq562 = ((m356*m356 - m3*m3 - m56*m56)**2 - (2.d0*m3*m56)**2)/(4.d0*m356*m356)
              if (rq562 < 0.d0) then
                fffxn = 0.d0
                call debug("fxn = 0. Skipping.")
                go to 999
              else
                rq56 = sqrt(rq562)
              endif
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
              if (rq782 < 0.d0) then
                fffxn = 0.d0
                call debug("fxn = 0. Skipping.")
                go to 999
              else
                rq78 = sqrt(rq782)
              endif
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
              if (rq52 < 0.d0) then
                fffxn = 0.d0
                call debug("fxn = 0. Skipping.")
                go to 999
              else
                rq5 = sqrt(rq52)
              endif
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
              if (rq72 < 0.d0) then
                fffxn = 0.d0
                call debug("fxn = 0. Skipping.")
                go to 999
              else
                rq7 = sqrt(rq72)
              endif
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
              call debug("Calculating 2to6 final state momenta in the parton CoM frame using RAMBO...")
              xmass(1) = m3
              xmass(2) = m4
              xmass(3) = m5 
              xmass(4) = m6
              xmass(5) = m7
              xmass(6) = m8
              jps = 6
              call rambo(seed,jps,ecm,xmass,prambo,wgtr)
              do i = 3, jps + 2
                do j = 1, 4
                  q(j,i) = prambo(j,i-2)
                end do
              end do
            end if

              call debug("...complete.")
            end if

            call debug("Boosting parton CoM momenta to collider frame...")
            ! velocity of ttbar system in collider frame
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
            call debug("...complete.")

            call debug("Assigning particle 4-momenta...")

            ! parton CoM 4 -momenta
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

          call debug("..done.")

          if (verbose == 1) then
            if (final_state <= 0) then
              print*, "Checking 2to2 kinematics..."
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
              print*, "Checking 2to6 kinematics..."
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
            print*, "...complete."
          end if

          if (phase_space_only == 1) then
            call debug("Setting |M|=1 and skipping matrix element calculation...")
            if (ix == 1) then
              pfxtot = 0.5/x1
            else if (ix == 2) then
              pfxtot = 0.5/x2
            end if
            if (final_state <= 0) then
              do lam3 = -1,1,2
                do lam4 = -1,1,2
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
        
          call debug("Calculating QCD coupling...")
          a_s = alfas(qq,lambdaqcd4,nloops)
          gs2 = 4.d0*pi*a_s
          gs = sqrt(gs2)
          call debug("...complete.")

          ! initilise
          qcdqq = 0.d0
          qcdbb = 0.d0
          qcdgg = 0.d0
          ewzuu1 = 0.d0
          ewzdd1 = 0.d0
          ewzbb1 = 0.d0
          ewzuu2 = 0.d0
          ewzdd2 = 0.d0
          ewzbb2 = 0.d0
          do lam3 = -1, 1
            do lam4 = -1, 1
              qcdpolqq(lam3,lam4) = 0.d0
              qcdpolbb(lam3,lam4) = 0.d0
              qcdpolgg(lam3,lam4) = 0.d0
              ewzpoluu1(lam3,lam4) = 0.d0
              ewzpoldd1(lam3,lam4) = 0.d0
              ewzpolbb1(lam3,lam4) = 0.d0
              ewzpoluu2(lam3,lam4) = 0.d0
              ewzpoldd2(lam3,lam4) = 0.d0
              ewzpolbb2(lam3,lam4) = 0.d0
              pfx(lam3,lam4) = 0.d0
              do i = 1, 20
                weight(lam3,lam4,i) = 0.d0
              end do
            end do
          end do

          resall = 0
          if (final_state <= 0) then
            call debug("Computing 2to2 square matrix elements...")
            if (include_qcd == 1) then
              call debug("Computing QCD matrix elements...")
              do lam3 = -1, 1, 2
                do lam4 = -1, 1, 2
                  if (include_gg == 1) then 
                    qcdpolgg(lam3,lam4) = sggff_qcd(p1,p2,p3,p4,lam3,lam4)*gs**4
                  end if
                  if (include_qq == 1) then
                    qcdpolqq(lam3,lam4) = sqqff_qcd(3 ,p1,p2,p3,p4,lam3,lam4)*gs**4
                    qcdpolbb(lam3,lam4) = sqqff_qcd(12,p1,p2,p3,p4,lam3,lam4)*gs**4
                  end if
                  resall = resall + qcdpolgg(lam3,lam4) + qcdpolqq(lam3,lam4) + qcdpolbb(lam3,lam4)
                end do
              end do
              call debug("...complete.")
            end if
            if ((include_qfd == 1) .or. (include_bsm == 1)) then
              call debug("Computing EW+Z' matrix elements...")
              do lam3 = -1,1,2
                do lam4 = -1,1,2
                  if (include_qq == 1) then
                    ewzpoluu1(lam3,lam4) = sqqff_ewp( 3,ffinal,p1,p2,p3,p4,lam3,lam4)
                    ewzpoluu2(lam3,lam4) = sqqff_ewp( 3,ffinal,p2,p1,p3,p4,lam3,lam4)
                    ewzpoldd1(lam3,lam4) = sqqff_ewp( 4,ffinal,p1,p2,p3,p4,lam3,lam4)
                    ewzpoldd2(lam3,lam4) = sqqff_ewp( 4,ffinal,p2,p1,p3,p4,lam3,lam4)
                    ewzpolbb1(lam3,lam4) = sqqff_ewp(12,ffinal,p1,p2,p3,p4,lam3,lam4)
                    ewzpolbb2(lam3,lam4) = sqqff_ewp(12,ffinal,p2,p1,p3,p4,lam3,lam4)
                  end if
                  resall = resall &
                 + ewzpoluu1(lam3,lam4) + ewzpoluu2(lam3,lam4) &
                 + ewzpoldd1(lam3,lam4) + ewzpoldd2(lam3,lam4) &
                 + ewzpolbb1(lam3,lam4) + ewzpolbb2(lam3,lam4)
                end do
              end do
              call debug("...complete.")
            end if

          else if (final_state > 0) then
            call debug("Computing 2to6 square matrix elements...")
            ! (Do not change the deliberate order of p6 and p7.)
            if (include_qcd == 1) then
              call debug("Computing QCD matrix elements...")
              if (include_gg == 1) then 
                qcdgg = sggbbffff_qcd(p1, p2, p3, p4, p5, p7, p6, p8)
              end if
              if (include_qq == 1) then 
                qcdqq = sqqbbffff_qcd(3 , p1, p2, p3, p4, p5, p7, p6, p8)
                qcdbb = sqqbbffff_qcd(12, p1, p2, p3, p4, p5, p7, p6, p8)
              end if              
              call debug("...complete.")
            end if
            if ((include_qfd == 1) .or. (include_bsm == 1)) then
              call debug("Computing EW+Z' matrix elements...")
              if (include_qq == 1) then 
                ewzuu1 = sqqbbffff_ewp( 3,11, p1, p2, p3, p4, p5, p7, p6, p8)
                ewzuu2 = sqqbbffff_ewp( 3,11, p2, p1, p3, p4, p5, p7, p6, p8)
                ewzdd1 = sqqbbffff_ewp( 4,11, p1, p2, p3, p4, p5, p7, p6, p8)
                ewzdd2 = sqqbbffff_ewp( 4,11, p2, p1, p3, p4, p5, p7, p6, p8)
                ewzbb1 = sqqbbffff_ewp(12,11, p1, p2, p3, p4, p5, p7, p6, p8)
                ewzbb2 = sqqbbffff_ewp(12,11, p2, p1, p3, p4, p5, p7, p6, p8)
              end if
              call debug("...complete.")
            end if
            resall = qcdqq + qcdgg + qcdbb + ewzuu1 + ewzdd1 + ewzbb1 &
                   + ewzuu2 + ewzdd2 + ewzbb2
          end if

          if (resall == 0.d0) then
            fffxn = 0.d0
            call debug("No |M^2| contributions. Setting fxn = 0 and skipping.")
            go to 999
          end if

          ! multiple qcd |m|^2 by g_s^4 (madgraph gs is set to one due to scale dependence.)
          qcdqq = qcdqq*gs**4
          qcdgg = qcdgg*gs**4
      
          pfxtot = 0.d0
          if (final_state <= 0) then
            call debug("Summing over 2to2 |m|^2 with pdfs of all initial partons..." )
            do lam3 = -1, 1, 2
              do lam4 = -1, 1, 2
                pfx(lam3,lam4) = qcdpolgg(lam3,lam4)*fx1(13)*fx2(13) &
                + (qcdpolqq(lam3,lam4) + ewzpoldd1(lam3,lam4))*fx1( 1)*fx2( 7) &
                + (qcdpolqq(lam3,lam4) + ewzpoluu1(lam3,lam4))*fx1( 2)*fx2( 8) &
                + (qcdpolqq(lam3,lam4) + ewzpoldd1(lam3,lam4))*fx1( 3)*fx2( 9) &
                + (qcdpolqq(lam3,lam4) + ewzpoluu1(lam3,lam4))*fx1( 4)*fx2(10) &
                + (qcdpolbb(lam3,lam4) + ewzpolbb1(lam3,lam4))*fx1( 5)*fx2(11) &
                + (qcdpolqq(lam3,lam4) + ewzpoldd2(lam3,lam4))*fx1( 7)*fx2( 1) &
                + (qcdpolqq(lam3,lam4) + ewzpoluu2(lam3,lam4))*fx1( 8)*fx2( 2) &
                + (qcdpolqq(lam3,lam4) + ewzpoldd2(lam3,lam4))*fx1( 9)*fx2( 3) &
                + (qcdpolqq(lam3,lam4) + ewzpoluu2(lam3,lam4))*fx1(10)*fx2( 4) &
                + (qcdpolbb(lam3,lam4) + ewzpolbb2(lam3,lam4))*fx1(11)*fx2( 5)
                if (ix == 1) then
                  pfx(lam3,lam4) = pfx(lam3,lam4)/x1
                else if (ix == 2) then
                  pfx(lam3,lam4) = pfx(lam3,lam4)/x2
                end if
                pfxtot = pfxtot + pfx(lam3,lam4)
              end do
            end do
            call debug("...complete.")
          else if (final_state > 0) then
            call debug("Summing over 2to6 |m|^2 with PDFs of all initial partons..." )
            pfxtot = fx1( 1)*fx2( 7)*(qcdqq + ewzdd1) &
                   + fx1( 2)*fx2( 8)*(qcdqq + ewzuu1) &
                   + fx1( 3)*fx2( 9)*(qcdqq + ewzdd1) &
                   + fx1( 4)*fx2(10)*(qcdqq + ewzuu1) &
                   + fx1( 5)*fx2(11)*(qcdbb + ewzbb1) &
                   + fx1( 7)*fx2( 1)*(qcdqq + ewzdd2) &
                   + fx1( 8)*fx2( 2)*(qcdqq + ewzuu2) &
                   + fx1( 9)*fx2( 3)*(qcdqq + ewzdd2) &
                   + fx1(10)*fx2( 4)*(qcdqq + ewzuu2) &
                   + fx1(11)*fx2( 5)*(qcdbb + ewzbb2) &
                   + fx1(13)*fx2(13)*qcdgg
            if (ix == 1) then
              pfxtot = (pfxtot)/x1
            else if (ix  ==  2) then
              pfxtot = (pfxtot)/x2
            end if
            call debug("...complete." )
          end if

          if (pfxtot == 0.d0) then
              fffxn = 0.d0
              call debug("fxn = 0. Skipping.")
            go to 999
          end if

          if (final_state <= 0) then
            ! weight for distributions
            do lam3 = -1, 1, 2
              do lam4 = -1, 1, 2
                pfx(lam3,lam4) = pfx(lam3,lam4)/(pfxtot)
              end do
            end do
          end if

          666 continue
          if (phase_space_only == 1) call debug("...complete.")

          call debug("Multiplying by jacobian from dx1 dx2 -> dx(2) dx(3)...")
          pfxtot = pfxtot*(1.d0 - tau)*2.d0*ecm/s &
          *(ecm_max - m3 - m4 - m5 - m6 - m7 - m8)
          call debug("...complete.")

          fffxn = pfxtot

          call debug("Applying unit converstion")
          fffxn = fffxn*unit_conv
          call debug("...complete.")

          call debug("Multiplying by phase space volume and flux factor and azimuthal integration...")
          if (final_state <= 0) then
            if (use_rambo == 0) then
              ! 2-body phase space factor + azimuthal integration
              fffxn = fffxn*qcm/(2.d0*pcm)*2.d0**(4 - 3*(2))*2.d0*pi
            else if (use_rambo == 1) then
              fffxn = fffxn*wgtr
            end if

            ! flux factor
            fffxn = fffxn/2.d0/ecm/ecm*(2.d0*pi)**(4 - 3*(2))

          else if (final_state > 0) then
            if (use_rambo == 0) then

              ! phase space factor
              fffxn = fffxn*rq*rq56*rq78*rq5*rq7/ecm*256.d0*2.d0**(4 - 3*(6))*2.d0*pi

              if (map_phase_space == 1) then
                fffxn = fffxn*((m356*m356 - rmt*rmt)**2 + rmt**2*gamt**2)*(xx356max - xx356min)/(2.d0*m356)/rmt/gamt        
                fffxn = fffxn*((m478*m478 - rmt*rmt)**2 + rmt**2*gamt**2)*(xx478max - xx478min)/(2.d0*m478)/rmt/gamt
                fffxn = fffxn*((m56*m56 - rm_w*rm_w)**2 + rm_w**2*gamma_w**2)*(xx56max - xx56min)/(2.d0*m56)/rm_w/gamma_w
                fffxn = fffxn*((m78*m78 - rm_w*rm_w)**2 + rm_w**2*gamma_w**2)*(xx78max - xx78min)/(2.d0*m78)/rm_w/gamma_w
                ! nwa
                fffxn = fffxn*gamt/gamma_t*gamt/gamma_t
              else
                fffxn = fffxn*(m356max - m356min)   
                fffxn = fffxn*(m478max - m478min)
                fffxn = fffxn*(m56max - m56min)
                fffxn = fffxn*(m78max - m78min)
              end if
            else if (use_rambo == 1) then
              fffxn = fffxn*wgtr
            end if
          
            ! flux factor
            fffxn = fffxn/2.d0/ecm/ecm*(2.d0*pi)**(4 - 3*(6))
          end if

          fffxn = fffxn/real(ixmax)/real(jxmax)/real(i5max)/real(i7max)

          ! binning
          hist = fffxn*wgt
          
          call debug("...complete.")

          if (final_state <= 0) then
            call debug("Computing polarised event weightings...")
            do lam3 = -1, +1, 2
              do lam4 = -1, +1, 2
                sigma_pol(lam3,lam4,it) = sigma_pol(lam3,lam4,it) + fffxn*wgt*pfx(lam3,lam4)
                weight(lam3,lam4,it) = fffxn*wgt*pfx(lam3,lam4)
                error_pol(lam3,lam4,it) = error_pol(lam3,lam4,it) + sigma_pol(lam3,lam4,it)**2
              end do
            end do
            weightLL = weight(-1,-1,it)
            weightLR = weight(-1, 1,it)
            weightRL = weight( 1,-1,it)
            weightRR = weight( 1, 1,it)
            call debug("...complete.")
          end if

          call debug("Writing final particle collider frame momenta to Ntuple...")
          if (final_state == -1) then
            call rootaddparticle(11,qcol(1,3),qcol(2,3),qcol(3,3),qcol(4,3))
            call rootaddparticle(-11,qcol(1,4),qcol(2,4),qcol(3,4),qcol(4,4))
          else if (final_state == 0) then
            call rootaddparticle(6,qcol(1,3),qcol(2,3),qcol(3,3),qcol(4,3))
            call rootaddparticle(-6,qcol(1,4),qcol(2,4),qcol(3,4),qcol(4,4))
          else if (final_state == 1) then
            call rootaddparticle(5,qcol(1,3),qcol(2,3),qcol(3,3),qcol(4,3))
            call rootaddparticle(-5,qcol(1,4),qcol(2,4),qcol(3,4),qcol(4,4))
            call rootaddparticle(-11,qcol(1,5),qcol(2,5),qcol(3,5),qcol(4,5))
            call rootaddparticle(12,qcol(1,6),qcol(2,6),qcol(3,6),qcol(4,6))
            call rootaddparticle(11,qcol(1,7),qcol(2,7),qcol(3,7),qcol(4,7))
            call rootaddparticle(-12,qcol(1,8),qcol(2,8),qcol(3,8),qcol(4,8))
          end if
             
          call rootadddouble(weightLL, "weightLL")
          call rootadddouble(weightLR, "weightLR")
          call rootadddouble(weightRL, "weightRL")
          call rootadddouble(weightRR, "weightRR")
          
          ! convert results to different tt classifications
          call rootadddouble(hist*fac_ee,"weight_ee")
          call rootadddouble(hist*fac_emu,"weight_emu")
          call rootadddouble(hist*fac_eq,"weight_eq")
          call rootadddouble(hist*fac_qq,"weight_qq")
          call rootaddint(it,"iteration")
          call rootaddevent(hist)

          call debug("...complete.")

          ! stats
          npoints = npoints + 1
          if (verbose == 1 ) print*, "...event", npoints, "complete."
          999 continue
          fffffxn = fffffxn + fffxn
        end do
        ffffxn = ffffxn + fffffxn
      end do
      ffxn = ffxn + ffffxn
    end do
    dsigma = dsigma + ffxn  
  end do
  return
end function dsigma
