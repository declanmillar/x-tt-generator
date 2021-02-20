module modelling

  use vamp_kinds
  use helas
  use configuration

  implicit none

  real(kind=default), public :: lambdaqcd4

  ! SM couplings
  real(kind=default) :: gw, gwwa, gwwZ
  real(kind=default) :: gal(2), gau(2), gad(2), gwf(2)
  real(kind=default) :: gZn(2), gZl(2), gZu(2), gZd(2), g1(2)
  real(kind=default) :: gwwh, gZZh, ghhh, gwwhh, gZZhh, ghhhh
  complex(kind=complex) :: gchf(2, 12)
  real(kind=default) :: gg(2), g

  ! fermion masses and widths
  real(kind=default) :: fmass(12), fwidth(12), qmass(6), lmass(6)

  real(kind=default) :: tmass2, wmass2

  ! quark masses
  real(kind=default), parameter :: umass = 0.d0, cmass = 0.d0, tmass = 172.5d0
  real(kind=default), parameter :: dmass = 0.d0, smass = 0.d0, bmass = 4.18d0

  ! leptons masses
  real(kind=default), parameter :: emass = 0.d0, mumass = 0.d0, taumass = 1.78d0
  real(kind=default), parameter :: nuemass = 0.d0, numumass = 0.d0, nutaumass = 0.d0

  ! quark widths
  real(kind=default), parameter :: uwidth = 0.d0, cwidth = 0.d0, twidth = 1.3d0
  real(kind=default), parameter :: dwidth = 0.d0, swidth = 0.d0, bwidth = 0.d0

  ! lepton widths
  real(kind=default), parameter :: ewidth = 0.d0, muwidth = 0.d0, tauwidth = 0.d0
  real(kind=default), parameter :: nuewidth = 0.d0, numuwidth = 0.d0, nutauwidth = 0.d0

  ! SM boson masses
  real(kind=default), parameter :: wmass = 80.4d0, zmass = 91.19d0, wwidth = 2.08d0, zwidth = 2.5d0
  real(kind=default), parameter :: amass = 0.d0, awidth = 0.d0, hmass = 125.d0, hwidth = 3.1278d-3

  ! Other SM parameters
  real(kind=default), parameter :: a_em = 0.0078125d0, s2w = 0.2320d0, vev = 246.d0
  real(kind=default) :: cotw

  ! zprime parameters
  real(kind=default) :: xmass(5), xwidth(5)
  real(kind=default) :: xparam(5)
  real(kind=default) :: gp(5), gV_d(5), gA_d(5), gV_u(5), gA_u(5), ga_l(5), gv_l(5), gv_v(5), ga_v(5)
  real(kind=default) :: gxd(2, 5), gxu(2, 5), gxl(2, 5), gxn(2, 5), gxb(2, 5), gxt(2, 5), gxl3(2, 5), gxn3(2, 5)
  integer :: manual_width(5)

  ! model parameters
  integer :: model_type
  real(kind=default) :: paramx, paramsin2phi

  ! methods
  public :: initialise_model
  private :: initialise_standard_model
  private :: initialise_zprimes
  private :: reset_zprimes
  private :: convert_zprime_couplings
  private :: width_zprimes
  private :: print_model

contains

  subroutine initialise_model

    if (verbose) print *, "modelling: initialising model ..."

    call initialise_standard_model
    call initialise_zprimes
    call print_model

  end subroutine initialise_model

  subroutine initialise_standard_model

    integer i

    if (verbose) print *, "modelling: initialising standard model ..."

    fmass(1) = emass
    fmass(2) = nuemass
    fmass(3) = umass
    fmass(4) = dmass
    fmass(5) = mumass
    fmass(6) = numumass
    fmass(7) = cmass
    fmass(8) = smass
    fmass(9) = taumass
    fmass(10) = nutaumass
    fmass(11) = tmass
    fmass(12) = bmass

    qmass(1) = umass
    qmass(2) = dmass
    qmass(3) = cmass
    qmass(4) = smass
    qmass(5) = tmass
    qmass(6) = bmass

    lmass(1) = emass
    lmass(2) = nuemass
    lmass(3) = mumass
    lmass(4) = numumass
    lmass(5) = taumass
    lmass(6) = nutaumass

    fwidth(1) = ewidth
    fwidth(2) = nuewidth
    fwidth(3) = uwidth
    fwidth(4) = dwidth
    fwidth(5) = muwidth
    fwidth(6) = numuwidth
    fwidth(7) = cwidth
    fwidth(8) = swidth
    fwidth(9) = tauwidth
    fwidth(10) = nutauwidth

    if (nwa) then
      fwidth(11) = 1.d-5
    else
      fwidth(11) = twidth
    end if

    fwidth(12) = bwidth

    tmass2 = tmass*tmass
    wmass2 = wmass*wmass

    ! call SM HELAS couplings
    call coup1x(s2w, gw, gwwa, gwwZ)
    call coup2x(s2w, gal, gau, gad, gwf, gZn, gZl, gZu, gZd, g1)
    call coup3x(s2w, zmass, hmass, gwwh, gZZh, ghhh, gwwhh, gZZhh, ghhhh)
    do i = 1, 12
      call coup4x(s2w, zmass, fmass(i), gchf(1, i))
    end do

    ! QCD couplings (set to one, multiply correctly in zprime)
    g = 1.d0
    gg(1) = -g
    gg(2) = -g

  end subroutine initialise_standard_model

  subroutine initialise_zprimes

    integer :: i
    integer, parameter :: mdl = 30

    call reset_zprimes

    if (verbose) print *, "modelling: initialising Z' parameters ..."

    ! read model file
    open (unit=mdl, file='models/'//trim(model_name)//'.mdl', status='old')
    read (mdl, *) model_type
    if (model_type == 0) then
      read (mdl, *) xmass
      read (mdl, *) xwidth
      read (mdl, *) gp
      read (mdl, *) xparam
      read (mdl, *) gV_u
      read (mdl, *) gA_u
      read (mdl, *) gV_d
      read (mdl, *) gA_d
      read (mdl, *) gV_l
      read (mdl, *) gA_l
      read (mdl, *) gV_v
      read (mdl, *) gA_v
    else if (model_type == 1) then
      read (mdl, *) paramx
      read (mdl, *) paramsin2phi
    else
      print *, "error: invalid model type"
    end if
    close (mdl)

    if (model_type == 0) then
      do i = 1, 5
        if ((xmass(i) > 0.d0) .and. (xwidth(i) < 0.d0)) then
          manual_width(i) = 0
        else
          manual_width(i) = 1
        end if
      end do

      ! convert from VA to LR couplings
      call convert_zprime_couplings

    else if (model_type == 1) then
      call initialise_non_universal
    end if

    ! Calculate benchmark widths
    call width_zprimes

    return

  end subroutine initialise_zprimes

  subroutine reset_zprimes

    integer :: i, j

    if (verbose) print *, "modelling: resetting Z' parameters"

    do i = 1, 5
      xmass(i) = 0
      xwidth(i) = 0
      gV_u(i) = 0
      gV_d(i) = 0
      gV_l(i) = 0
      gV_v(i) = 0
      gA_u(i) = 0
      gA_d(i) = 0
      gA_l(i) = 0
      gA_v(i) = 0
      do j = 1, 2
        gxd(j, i) = 0
        gxu(j, i) = 0
        gxl(j, i) = 0
        gxn(j, i) = 0
        gxb(j, i) = 0
        gxt(j, i) = 0
        gxl3(j, i) = 0
        gxn3(j, i) = 0
      end do
    end do

  end subroutine reset_zprimes

  subroutine initialise_non_universal

    integer :: i, j
    real(kind=default) :: x, sin2phi, e, st, ct, sp, cp, m0

    if (verbose) print *, "modelling: initialising non-universal model ..."

    call reset_zprimes

    x = paramx
    sin2phi = paramsin2phi

    e = sqrt(4.d0*pi*a_em)
    st = sqrt(s2w)
    ct = sqrt(1.d0 - s2w)
    sp = sqrt(sin2phi)
    cp = sqrt(1.d0 - sin2phi)
    m0 = e*vev/2/st
    xmass(1) = sqrt(m0*m0*(x/sp/sp/cp/cp + sp*sp/cp/cp))

    ! leptons
    gxl(1, 1) = e/2/st*(sp/cp + sp*sp*sp*cp/(x*ct*ct)*(1 - 2*st*st))
    gxn(1, 1) = e/2/st*(-sp/cp - sp*sp*sp*cp/(x*ct*ct))
    gxl3(1, 1) = gxl(1, 1) - e/(2*st*sp*cp)
    gxn3(1, 1) = gxn(1, 1) + e/(2*st*sp*cp)

    ! quarks
    gxu(1, 1) = e/(2*st)*(-sp/cp - sp*sp*sp*cp/(x*ct*ct)*(1 - 4*st*st/3))
    gxd(1, 1) = e/(2*st)*(sp/cp + sp*sp*sp*cp/(x*ct*ct)*(1 - 2*st*st/3))
    gxt(1, 1) = gxu(1, 1) + e/(2*st*sp*cp)
    gxb(1, 1) = gxd(1, 1) - e/(2*st*sp*cp)

    ! right handed couplings
    gxl(2, 1) = e/2*st*2*st*st*sp*sp*sp*cp/(x*ct*ct)
    gxn(2, 1) = 0
    gxl3(2, 1) = gxl(2, 1)
    gxn3(2, 1) = gxn(2, 1)

    gxu(2, 1) = gxl(2, 1)*2/3
    gxd(2, 1) = gxl(2, 1)*1/3
    gxt(2, 1) = gxu(2, 1)
    gxb(2, 1) = gxd(2, 1)

    manual_width(1) = 0
    do i = 2, 5
      manual_width(i) = 1
    end do

    return

  end subroutine initialise_non_universal

  subroutine convert_zprime_couplings

    ! input: vector and axial Z' couplings to up and down quarks
    ! output: left and right chiral couplings to up and down quarks

    integer i

    if (verbose) print *, "modelling: converting couplings ..."
    do i = 1, 5
      gxd(1, i) = gp(i)*(gv_d(i) + ga_d(i))/2.d0
      gxd(2, i) = gp(i)*(gv_d(i) - ga_d(i))/2.d0
      gxu(1, i) = gp(i)*(gv_u(i) + ga_u(i))/2.d0
      gxu(2, i) = gp(i)*(gv_u(i) - ga_u(i))/2.d0
      gxl(1, i) = gp(i)*(gv_l(i) + ga_l(i))/2.d0
      gxl(2, i) = gp(i)*(gv_l(i) - ga_l(i))/2.d0
      gxn(1, i) = gp(i)*(gv_v(i) + ga_v(i))/2.d0
      gxn(2, i) = gp(i)*(gv_v(i) - ga_v(i))/2.d0
      gxb(1, i) = gxd(1, i)
      gxb(2, i) = gxd(2, i)
      gxt(1, i) = gxu(1, i)
      gxt(2, i) = gxu(2, i)
      gxl3(1, i) = gxl(1, i)
      gxl3(2, i) = gxl(2, i)
      gxn3(1, i) = gxn(1, i)
      gxn3(2, i) = gxn(2, i)
    end do

    return
  end subroutine convert_zprime_couplings

  subroutine width_zprimes

    implicit none

    integer :: i, n
    real(kind=default) :: mq, ml, mx
    real(kind=default) :: gv, ga
    real(kind=default) :: width, widthqq, widthll, widthqq_tmp, widthll_tmp
    real(kind=default) :: widthzh, widthww, fzh
    real(kind=default) :: pi
    real(kind=default) :: a_s, alfas
    real(kind=default) :: e, st, ct, cotw, gz
    real(kind=default) :: stmix, ctmix
    real(kind=default) :: ez, pz

    if (verbose) print *, "modelling: calculating Z' -> ff widths ..."

    ! couplings.
    pi = dacos(-1.d0)

    e = sqrt(4.d0*pi*a_em)

    st = sqrt(s2w)
    ct = sqrt(1.d0 - s2w)

    stmix = sqrt(s2mix)
    ctmix = sqrt(1.d0 - s2mix)

    cotw = ct/st
    gz = e/(ct*st)

    do n = 1, 5
      width = 0.d0
      mx = xmass(n)
      if (manual_width(n) == 0) then
        a_s = alfas(mx, lambdaqcd4, nloops)
        ! quarks
        widthqq = 0.d0
        do i = 1, 6
          widthqq_tmp = 0.d0
          mq = qmass(i)
          if (xmass(n) > 2.d0*mq) then

            if (i == 1 .or. i == 3) then
              gv = gxu(1, n) + gxu(2, n)
              ga = gxu(1, n) - gxu(2, n)
            else if (i == 2 .or. i == 4) then
              gv = gxd(1, n) + gxd(2, n)
              ga = gxd(1, n) - gxd(2, n)
            else if (i == 5) then
              gv = gxt(1, n) + gxt(2, n)
              ga = gxt(1, n) - gxt(2, n)
            else if (i == 6) then
              gv = gxb(1, n) + gxb(2, n)
              ga = gxb(1, n) - gxb(2, n)
            end if

            ! with QCD kfactor
            widthqq_tmp = 3.d0/48.d0/pi*mx*sqrt(1.d0 - 4.d0*mq*mq/mx/mx) &
                          *(gv*gv*(1.d0 + 2.d0*mq*mq/mx/mx) + ga*ga*(1.d0 - 4.d0*mq*mq/mx/mx)) &
                          *(1.d0 + 1.045d0*a_s/pi)

            widthqq = widthqq + widthqq_tmp

          end if
        end do

        ! leptons
        widthll = 0.d0
        do i = 1, 6
          widthll_tmp = 0.d0
          ml = lmass(i)

          if (mx > 2.d0*ml) then

            if (i == 1 .or. i == 3) then
              gv = gxl(1, n) + gxl(2, n)
              ga = gxl(1, n) - gxl(2, n)
            else if (i == 2 .or. i == 4) then
              gv = gxn(1, n) + gxn(2, n)
              ga = gxn(1, n) - gxn(2, n)
            else if (i == 5) then
              gv = gxl3(1, n) + gxl3(2, n)
              ga = gxl3(1, n) - gxl3(2, n)
            else if (i == 6) then
              gv = gxn3(1, n) + gxn3(2, n)
              ga = gxn3(1, n) - gxn3(2, n)
            end if

            widthll_tmp = 1.d0/48.d0/pi*mx*sqrt(1.d0 - 4.d0*mq*mq/mx/mx) &
                          *(gv*gv*(1.d0 + 2.d0*mq*mq/mx/mx) + ga*ga*(1.d0 - 4.d0*mq*mq/mx/mx))

            widthll = widthll + widthll_tmp
          end if
        end do

        width = widthqq + widthll

        print *, " Gamma_Z'(", n, ")->qq = ", widthqq, " [GeV]"
        print *, " Gamma_Z'(", n, ")->ll = ", widthll, " [GeV]"
        print *, " Gamma_Z'(", n, ")->ff = ", width, " [GeV]"

        if (z_mixing == 1) then
          if (verbose) print *, "modelling: calculating Z' mixing ..."
          widthww = 1/(48.d0*pi)*e*e*cotw*cotw*stmix*mx*sqrt(1 - 4*wmass*wmass/mx/mx) &
                    *(0.25*(mx/wmass)**4 + 4*mx*mx/wmass/wmass - 17 - 12*wmass*wmass/mx/mx)

          fzh = -gz*mx/zmass*ctmix*stmix/(ctmix*ctmix + mx*mx/zmass/zmass*stmix*stmix)**(3/2)

          ez = (mx*mx + zmass*zmass - hmass*hmass)/(2*mx)

          pz = sqrt(ez*ez + zmass*zmass)

          widthzh = 1/(24.d0*pi)*fzh*fzh*pz*(ez*ez/zmass/zmass + 2)

          width = width + widthww + widthzh

          print *, "Gamma_Z'(", n, ") -> WW = ", widthww, " [GeV]"
          print *, "Gamma_Z'(", n, ") -> hZ = ", widthzh, " [GeV]"
        end if

        xwidth(n) = width
      end if
    end do
    if (verbose) print *, "modelling: calculated Z' -> ff widths"

    return
  end subroutine width_zprimes

  subroutine print_model

    integer :: i

    if (verbose) print *, "modelling: printing model parameters ..."

    print *, "m_b = ", fmass(12)
    print *, "Gamma_b = ", fwidth(12)
    print *, "m_t = ", fmass(11)
    print *, "Gamma_t = ", fwidth(11)
    print *, "m_Z = ", zmass
    print *, "Gamma_Z = ", zwidth
    print *, "m_W = ", wmass
    print *, "Gamma_W = ", wwidth
    print *, "m_h = ", hmass
    print *, "Gamma_h = ", hwidth

    if (include_x) then
      do i = 1, 5
        if (xmass(i) > 0) then
          print *, "m_Z'(", i, ") = ", xmass(i)
          print *, "Gamma_Z'(", i, ") = ", xwidth(i)
          print *, "gL_uZ'(", i, ") = ", gxu(1, i)
          print *, "gR_uZ'(", i, ") = ", gxu(2, i)
          print *, "gL_dZ'(", i, ") = ", gxd(1, i)
          print *, "gR_dZ'(", i, ") = ", gxd(2, i)
          print *, "gL_lZ'(", i, ") = ", gxl(1, i)
          print *, "gR_lZ'(", i, ") = ", gxl(2, i)
          print *, "gL_nZ'(", i, ") = ", gxn(1, i)
          print *, "gR_nZ'(", i, ") = ", gxn(2, i)
          print *, "gL_tZ'(", i, ") = ", gxt(1, i)
          print *, "gR_tZ'(", i, ") = ", gxt(2, i)
          print *, "gL_bZ'(", i, ") = ", gxb(1, i)
          print *, "gR_bZ'(", i, ") = ", gxb(2, i)
          print *, "gL_taZ'(", i, ") = ", gxl3(1, i)
          print *, "gR_taZ'(", i, ") = ", gxl3(2, i)
          print *, "gL_vtZ'(", i, ") = ", gxn3(1, i)
          print *, "gR_vtZ'(", i, ") = ", gxn3(2, i)
        end if
      end do
    end if
  end subroutine print_model

end module modelling
