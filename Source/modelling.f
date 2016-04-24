module modelling

  use configuration

  implicit none

  ! SM couplings
  real :: gw, gwwa, gwwZ
  real :: gal(2), gau(2), gad(2), gwf(2)
  real :: gZn(2), gZl(2), gZu(2), gZd(2), g1(2)
  real :: gwwh, gZZh, ghhh, gwwhh, gZZhh, ghhhh
  complex*16 :: gchf(2, 12)
  real :: gg(2), g

  ! fermion masses and widths
  real :: fmass(12), fwidth(12), qmass(6), lmass(6)

  ! quark masses
  real, parameter :: umass = 0.d0, cmass = 0.d0, tmass = 173.d0
  real, parameter :: dmass = 0.d0, smass = 0.d0, bmass = 4.18d0

  ! leptons masses
  real, parameter :: emass = 0.d0, mumass = 0.d0, taumass = 1.78d0
  real, parameter :: nuemass = 0d0, numumass = 0d0, nutaumass = 0d0

  ! quark widths
  real, parameter :: uwidth = 0.d0,  cwidth = 0.d0, twidth = 1.55d0
  real, parameter :: dwidth = 0.d0,  swidth = 0.d0, bwidth = 0.d0

  ! lepton widths
  real, parameter :: ewidth = 0.d0,   muwidth = 0.d0, tauwidth = 0.d0
  real, parameter :: nuewidth = 0.d0, numuwidth = 0.d0, nutauwidth = 0.d0

  ! SM boson masses
  real, parameter :: wmass = 80.23d0, zmass = 91.19d0, wwidth = 2.08d0, zwidth = 2.5d0
  real, parameter :: amass = 0d0, awidth = 0d0, hmass = 125.d0, hwidth = 0.31278d-2

  real :: Gamma_t = twidth

  ! Other SM parameters
  real, parameter :: a_em = 0.0078125, s2w = 0.2320d0, vev = 246.d0
  real :: cotw

  ! zprime parameters
  real :: mass_zp(5), gamZp(5)
  real :: paramZp(5)
  real :: gp(5), gV_d(5), gA_d(5), gV_u(5), gA_u(5), ga_l(5), gv_l(5), gv_nu(5), ga_nu(5)
  real :: gZpd(2,5),gZpu(2,5),gZpl(2,5),gZpn(2,5),gZpb(2,5), gZpt(2,5), gZpl3(2,5), gZpn3(2,5)
  integer :: manual_width(5)

  ! model parameters
  integer :: model_type
  real :: xparam, sin2phiparam

  ! methods
  public :: initialise_model
  private :: initialise_standard_model
  private :: initialise_zprimes
  private :: reset_zprimes
  private :: convert_zprime_couplings
  private :: width_zprimes
  private :: width_zprime_ssm
  private :: print_model

contains

subroutine initialise_model

  call initialise_standard_model
  call initialise_zprimes
  call print_model

end subroutine initialise_model

subroutine initialise_standard_model

  integer i

  fmass(1)  = emass
  fmass(2)  = nuemass
  fmass(3)  = umass
  fmass(4)  = dmass
  fmass(5)  = mumass
  fmass(6)  = numumass
  fmass(7)  = cmass
  fmass(8)  = smass
  fmass(9)  = taumass
  fmass(10) = nutaumass
  fmass(11) = tmass
  fmass(12) = bmass

  qmass(1)  = umass
  qmass(2)  = dmass
  qmass(3)  = cmass
  qmass(4)  = smass
  qmass(5)  = tmass
  qmass(6)  = bmass

  lmass(1)  = emass
  lmass(2)  = nuemass
  lmass(3)  = mumass
  lmass(4)  = numumass
  lmass(5)  = taumass
  lmass(6)  = nutaumass

  fwidth(1)  = ewidth
  fwidth(2)  = nuewidth
  fwidth(3)  = uwidth
  fwidth(4)  = dwidth
  fwidth(5)  = muwidth
  fwidth(6)  = numuwidth
  fwidth(7)  = cwidth
  fwidth(8)  = swidth
  fwidth(9)  = tauwidth
  fwidth(10) = nutauwidth

  if (use_NWA == 1) then
    fwidth(11) = 1.d-5
  else
    fwidth(11) = twidth
  end if

  fwidth(12) = bwidth

  ! call SM HELAS couplings
  call coup1x(s2w, gw, gwwa, gwwZ)
  call coup2x(s2w, gal, gau, gad, gwf, gZn, gZl, gZu, gZd, g1)
  call coup3x(s2w, zmass, hmass, gwwh, gZZh, ghhh, gwwhh, gZZhh, ghhhh)
  do i = 1, 12
    call coup4x(s2w, zmass, fmass(i), gchf(1,i))
  enddo

  ! QCD couplings (set to one, multiply correctly in zprime)
  g = 1.d0
  gg(1) = -g
  gg(2) = -g

end subroutine initialise_standard_model

subroutine initialise_zprimes

  integer imodel_name, i

  call reset_zprimes

  ! Extract model_name filename (Remove white space.)
  imodel_name = len(model_name)
  do while(model_name(imodel_name:imodel_name) == '')
    imodel_name = imodel_name - 1
  end do

  ! read model file
  open(unit = 42, file = 'Models/'//model_name(1:imodel_name)//'.mdl', status = 'old')
  read(42,*) model_type
  if (model_type == 0) then
    read(42,*) mass_zp
    read(42,*) gamZp
    read(42,*) gp
    read(42,*) paramZp
    read(42,*) gV_u
    read(42,*) gA_u
    read(42,*) gV_d
    read(42,*) gA_d
    read(42,*) gV_l
    read(42,*) gA_l
    read(42,*) gV_nu
    read(42,*) gA_nu
  else if (model_type == 1) then
    read(42,*) xparam
    read(42,*) sin2phiparam
  else
    print*, "Error: invalid model type! Must be 0-1."
  end if
  close(42)

  if (model_type == 0) then
    ! If gamZp is negative, the function widthZp is used.
    do i = 1, 5
      if ((mass_zp(i) > 0.d0) .and. (gamzp(i) < 0.d0)) then
        manual_width(i) = 0
      else
        manual_width(i) = 1
      end if
    enddo

    ! convert from VA to LR couplings
    call convert_zprime_couplings

  else if (model_type == 1) then
    call initialise_non_universal
  end if

  ! Calculate benchmark widths
  call width_zprimes

  ! igw=0 ! don't include w width effects
  ! call topwid(fmass(11),wmass,fmass(12),wwidth,igw,fwidth(11))
  ! call printconstants
  return

end subroutine initialise_zprimes

subroutine reset_zprimes
  integer :: i, j

  do i = 1, 5
    mass_zp(i) = 0
    gamZp(i) = 0
    gV_u(i) = 0
    gV_d(i) = 0
    gV_l(i) = 0
    gV_nu(i) = 0
    gA_u(i) = 0
    gA_d(i) = 0
    gA_l(i) = 0
    gA_nu(i) = 0
    do j = 1, 2
      gZpd(j,i) = 0
      gZpu(j,i) = 0
      gZpl(j,i) = 0
      gZpn(j,i) = 0
      gZpb(j,i) = 0
      gZpt(j,i) = 0
      gZpl3(j,i) = 0
      gZpn3(j,i) = 0
    end do
  end do

end subroutine reset_zprimes

subroutine initialise_non_universal

  integer :: i, j
  real :: x, sin2phi, e, st, ct, sp, cp, m0

  call reset_zprimes

!   x = 56
!   sin2phi = 0.1
  x = xparam
  sin2phi = sin2phiparam

  e = sqrt(4.d0*pi*a_em)
  st = sqrt(s2w)
  ct = sqrt(1.d0 - s2w)
  sp = sqrt(sin2phi)
  cp = sqrt(1.d0 - sin2phi)
  m0 = e*vev/2/st
  mass_zp(1) = sqrt(m0*m0*(x/sp/sp/cp/cp + sp*sp/cp/cp))

  ! leptons
  gZpl(1,1) = e/2/st*(sp/cp + sp*sp*sp*cp/(x*ct*ct)*(1 - 2*st*st))
  gZpn(1,1) = e/2/st*(-sp/cp - sp*sp*sp*cp/(x*ct*ct))
  gZpl3(1,1) = gZpl(1,1) - e/2/st/sp/cp
  gZpn3(1,1) = gZpn(1,1) + e/2/st/sp/cp

  ! quarks
  gZpu(1,1) = e/2/st*(-sp/cp - sp*sp*sp*cp/(x*ct*ct)*(1 - 4*st*st/3))
  gZpd(1,1) = e/2/st*(sp/cp + sp*sp*sp*cp/(x*ct*ct)*(1 - 2*st*st/3))
  gZpt(1,1) = gZpu(1,1) + e/2/st/sp/cp
  gZpb(1,1) = gZpd(1,1) - e/2/st/sp/cp

  ! right handed couplings
  gZpl(2,1) = e/2*st*2*st*st*sp*sp*sp*cp/x/ct/ct
  gZpn(2,1) = 0
  gZpl3(2,1) = gZpl(2,1)
  gZpn3(2,1) = gZpn(2,1)

  gZpu(2,1) = gZpl(2,1)*2/3
  gZpd(2,1) = gZpl(2,1)*1/3
  gZpt(2,1) = gzpu(2,1)
  gZpb(2,1) = gZpd(2,1)

  manual_width(1) = 0
  do i = 2, 5
    manual_width(i) = 1
  enddo

  return

end subroutine initialise_non_universal

subroutine convert_zprime_couplings

  ! input: vector and axial Zp couplings to up and down quarks
  ! output: left and right chiral couplings to up and down quarks

  integer i


  do i = 1, 5
      gZpd(1,i) = gp(i)*(gv_d(i) + ga_d(i))/2.d0
      gZpd(2,i) = gp(i)*(gv_d(i) - ga_d(i))/2.d0
      gZpu(1,i) = gp(i)*(gv_u(i) + ga_u(i))/2.d0
      gZpu(2,i) = gp(i)*(gv_u(i) - ga_u(i))/2.d0
      gZpl(1,i) = gp(i)*(gv_l(i) + ga_l(i))/2.d0
      gZpl(2,i) = gp(i)*(gv_l(i) - ga_l(i))/2.d0
      gZpn(1,i) = gp(i)*(gv_nu(i) + ga_nu(i))/2.d0
      gZpn(2,i) = gp(i)*(gv_nu(i) - ga_nu(i))/2.d0
      gZpb(1,i) = gZpd(1,i)
      gZpb(2,i) = gZpd(2,i)
      gZpt(1,i) = gZpu(1,i)
      gZpt(2,i) = gZpu(2,i)
      gZpl3(1,i) = gZpl(1,i)
      gZpl3(2,i) = gZpl(2,i)
      gZpn3(1,i) = gZpn(1,i)
      gZpn3(2,i) = gZpn(2,i)
  enddo


  return
end subroutine convert_zprime_couplings

subroutine width_zprimes

  ! calculates Z' width contributions from decay to fermions

  implicit none

  integer :: i, n
  real :: mq, ml, mzp
  real :: gv, ga
  real :: width, widthqq, widthll, widthqq_tmp, widthll_tmp
  real :: widthzh, widthww, fzh
  real :: pi
  real :: a_s, alfas
  real :: e, st, ct, cotw, gz
  real :: stmix, ctmix
  real :: ez, pz

  e = sqrt(4.d0*pi*a_em)

  st = sqrt(s2w)
  ct = sqrt(1.d0 - s2w)

  stmix = sqrt(s2mix)
  ctmix = sqrt(1.d0 - s2mix)

  cotw = ct/st
  gz = e/ct/st


  ! couplings.
  pi = dacos(-1.d0)

  do n = 1, 5
    width = 0.d0
    mzp = mass_zp(n)
    if (manual_width(n) == 0) then
      a_s = alfas(mzp, lambdaQCD4, nloops)
      ! quarks
      widthqq = 0.d0
      do i = 1, 6
        widthqq_tmp = 0.d0
        mq = qmass(i)
        if (mass_zp(n) > 2.d0*mq) then

          if (i == 1 .or. i == 3) then
            gv = gzpu(1,n) + gzpu(2,n)
            ga = gzpu(1,n) - gzpu(2,n)
          else if (i == 2 .or. i == 4) then
            gv = gzpd(1,n) + gzpd(2,n)
            ga = gzpd(1,n) - gzpd(2,n)
          else if (i == 5) then
            gv = gzpt(1,n) + gzpt(2,n)
            ga = gzpt(1,n) - gzpt(2,n)
          else if (i == 6) then
            gv = gzpb(1,n) + gzpb(2,n)
            ga = gzpb(1,n) - gzpb(2,n)
          end if

          ! with QCD kfactor
          widthqq_tmp = 3.d0/48.d0/pi*mzp &
                      *sqrt(1.d0 - 4.d0*mq**2/mzp**2) &
                      *(gv**2*(1.d0 + 2.d0*mq**2/mzp**2) &
                      + ga**2*(1.d0 - 4.d0*mq**2/mzp**2)) &
                      *(1.d0 + 1.045d0*a_s/pi)

          widthqq = widthqq + widthqq_tmp

        end if
      end do

      ! leptons
      widthll = 0.d0
      do i = 1, 6
        widthll_tmp = 0.d0
        ml = lmass(i)

        if (mzp > 2.d0*ml) then

          if (i == 1 .or. i == 3) then
            gv = gzpl(1,n) + gzpl(2,n)
            ga = gzpl(1,n) - gzpl(2,n)
          else if (i == 2 .or. i == 4) then
            gv = gzpn(1,n) + gzpn(2,n)
            ga = gzpn(1,n) - gzpn(2,n)
          else if (i == 5) then
            gv = gzpl3(1,n) + gzpl3(2,n)
            ga = gzpl3(1,n) - gzpl3(2,n)
          else if (i == 6) then
            gv = gzpn3(1,n) + gzpn3(2,n)
            ga = gzpn3(1,n) - gzpn3(2,n)
          end if

          widthll_tmp = 1.d0/48.d0/pi*mzp &
                        *sqrt(1.d0 - 4.d0*mq**2/mzp**2) &
                        *(gv**2*(1.d0 + 2.d0*mq**2/mzp**2) &
                        + ga**2*(1.d0 - 4.d0*mq**2/mzp**2))

          widthll = widthll + widthll_tmp
        end if
      end do

      width = widthqq + widthll

      print*, 'Gamma(Zp(', n, ')->ff)=', width,' [GeV]'
      print*, 'Gamma(Zp(', n, ')->qq)=', widthqq,' [GeV]'
      print*, 'Gamma(Zp(', n, ')->ll)=', widthll,' [GeV]'

      if (z_mixing == 1) then
        print*, "I am being run."
        widthww = 1/(48.d0*pi)*e*e*cotw*cotw*stmix*mzp*sqrt(1 - 4*wmass*wmass/mzp/mzp) &
                  *(0.25*(mzp/wmass)**4 + 4*mzp*mzp/wmass/wmass - 17 - 12*wmass*wmass/mzp/mzp)

        fzh = -gz*mzp/zmass*ctmix*stmix/(ctmix*ctmix + mzp*mzp/zmass/zmass*stmix*stmix)**(3/2)

        ez = (mzp*mzp + zmass*zmass - hmass*hmass)/(2*mzp)

        pz = sqrt(ez*ez + zmass*zmass)

        widthzh = 1/(24.d0*pi)*fzh*fzh*pz*(ez*ez/zmass/zmass + 2)

        width = width + widthww + widthzh

        print*, 'Gamma(Zp(', n, ')->qq)=', widthqq,' [GeV]'
        print*, 'Gamma(Zp(', n, ')->qq)=', widthww,' [GeV]'
      end if

      gamZp(n) = width
    end if
  end do
  return
end subroutine width_zprimes



function width_zprime_ssm(rm_Zp)

  ! Calculates the width of a Z' in the SSM.

  implicit none

  ! implicit to explicit variable dump
  real :: width_zprime_ssm
  real :: rm_zp
  real :: ctw
  real :: e
  real :: eq
  real :: gweak
  real :: gamt
  real :: gf
  integer :: i
  real :: pi
  real :: mq, rmt
  real :: t3q
  real :: temp, temp1, temp2
  real :: alfas, a_s


  rmt=fmass(11)
  gamt=fwidth(11)

  ! couplings.
  pi=dacos(-1.d0)
  ctw=sqrt(1.d0-s2w)
  e=sqrt(4.d0*pi*a_em)
  gweak=e/sqrt(s2w)
  a_s=alfas(rm_Zp,lambdaqcd4,nloops)
  GF=1.16639D-5
  ! renormalise e.
  e=sqrt(s2w*8.d0*zmass*zmass*ctw*ctw*GF/sqrt(2.d0))
  ! Zp width.
  width_zprime_ssm=0.d0
  do i=1,6
    mq=0.d0
    if(i == 6)mq=rmt
    if((i == 2) .OR. (i == 4) .OR. (i == 6)) then

    ! Quarks
    ! u-quark.
      t3q=+1.d0/2.d0
      eq=+2.d0/3.d0
    else if((i == 1) .OR. (i == 3) .OR. (i == 5)) then
    ! d-quark.
      t3q=-1.d0/2.d0
      eq=-1.d0/3.d0
    end if
    if(rm_Zp <= 2.d0*mq)goto 123
    width_zprime_ssm=width_zprime_ssm+3.d0*zmass**2*rm_Zp*GF/24.d0/pi/sqrt(2.d0) &
    *sqrt(1.d0-4.d0*mq**2/rm_Zp**2) &
    *((2.d0*t3q)**2 &
    *(1.d0-4.d0*mq**2/rm_Zp**2) &
    +(2.d0*t3q-4.d0*eq*s2w)**2 &
    *(1.d0+2.d0*mq**2/rm_Zp**2)) &
    *(1.d0+1.045d0*a_s/pi)
  !        cr=-eq*s2w
  !        cl=t3q-eq*s2w
  !        gv=cl+cr
  !        ga=cl-cr
  !        width_zprime_ssm=width_zprime_ssm+3.d0*e**2/12.d0/pi*rm_Zp/16.d0/s2w/ctw**2
  !     &               *sqrt(1.d0-4.d0*mq**2/rm_Zp**2)
  !     &               *(4.d0*gv**2*(1.d0+2.d0*mq**2/rm_Zp**2)
  !     &                +4.d0*ga**2*(1.d0-4.d0*mq**2/rm_Zp**2))
  !     &               *(1.d0+1.045d0*a_s/pi)
    123 continue
  end do

  !       print *,'Z'' width due to quarks: ',width_zprime_ssm,' [GeV]'

  ! Leptons
  temp=0.d0
  temp1=0.d0
  temp2=0.d0

  do i=1,6
    mq=0.d0
    if((i == 2) .OR. (i == 4) .OR. (i == 6)) then
    ! neutrino.
      t3q=+1.d0/2.d0
      eq=0.d0
    else if((i == 1) .OR. (i == 3) .OR. (i == 5)) then
    ! lepton.
      t3q=-1.d0/2.d0
      eq=-1.d0
    end if
    width_zprime_ssm=width_zprime_ssm+1.d0*zmass**2*rm_Zp*GF/24.d0/pi/sqrt(2.d0) &
    *sqrt(1.d0-4.d0*mq**2/rm_Zp**2) &
    *((2.d0*t3q)**2 &
    *(1.d0-4.d0*mq**2/rm_Zp**2) &
    +(2.d0*t3q-4.d0*eq*s2w)**2 &
    *(1.d0+2.d0*mq**2/rm_Zp**2))
  !        cr=-eq*s2w
  !        cl=t3q-eq*s2w
  !        gv=cl+cr
  !        ga=cl-cr
  !        width_zprime_ssm=width_zprime_ssm+1.d0*e**2/12.d0/pi*rm_Zp/16.d0/s2w/ctw**2
  !     &               *sqrt(1.d0-4.d0*mq**2/rm_Zp**2)
  !     &               *(4.d0*gv**2*(1.d0+2.d0*mq**2/rm_Zp**2)
  !     &                +4.d0*ga**2*(1.d0-4.d0*mq**2/rm_Zp**2))

    temp=temp+1.d0*zmass**2*rm_Zp*GF/24.d0/pi/sqrt(2.d0) &
    *sqrt(1.d0-4.d0*mq**2/rm_Zp**2) &
    *((2.d0*t3q)**2 &
    *(1.d0-4.d0*mq**2/rm_Zp**2) &
    +(2.d0*t3q-4.d0*eq*s2w)**2 &
    *(1.d0+2.d0*mq**2/rm_Zp**2))
    if((i == 2) .OR. (i == 4) .OR. (i == 6)) &
    temp1=temp1+1.d0*zmass**2*rm_Zp*GF/24.d0/pi/sqrt(2.d0) &
    *sqrt(1.d0-4.d0*mq**2/rm_Zp**2) &
    *((2.d0*t3q)**2 &
    *(1.d0-4.d0*mq**2/rm_Zp**2) &
    +(2.d0*t3q-4.d0*eq*s2w)**2 &
    *(1.d0+2.d0*mq**2/rm_Zp**2))
    if((i == 1) .OR. (i == 3) .OR. (i == 5)) &
    temp2=temp2+1.d0*zmass**2*rm_Zp*GF/24.d0/pi/sqrt(2.d0) &
    *sqrt(1.d0-4.d0*mq**2/rm_Zp**2) &
    *((2.d0*t3q)**2 &
    *(1.d0-4.d0*mq**2/rm_Zp**2) &
    +(2.d0*t3q-4.d0*eq*s2w)**2 &
    *(1.d0+2.d0*mq**2/rm_Zp**2))

  end do

!         print *,'Z'' width due to quarks+leptons:',width_zprime_ssm,' [GeV]'
  !       print *,'(so that due to leptons are:',temp,' [GeV])'
  !       print *,'(of which due to e/mu/tau:',temp2,' [GeV])'
  !       print *,'(of which due to their neutrinos are:',temp1,' [GeV])'
  return
end function width_zprime_ssm

subroutine print_model

  integer :: i

  write(log,*) 'm_b:', fmass(12)
  write(log,*) 'Gamma_b:', fwidth(12)
  write(log,*) 'm_t:', fmass(11)
  write(log,*) 'Gamma_t:', fwidth(11)
  write(log,*) 'm_Z:', zmass
  write(log,*) 'Gamma_Z:', zwidth
  write(log,*) 'm_W:', wmass
  write(log,*) 'Gamma_W:', wwidth
  write(log,*) 'm_h:', hmass
  write(log,*) 'Gamma_h:', hwidth
  if (include_x == 1) then
    do i = 1, 5
      if (mass_zp(i) > 0) then
        write(log,*) "Z' no.:", i
        write(log,*) "m_Z':", mass_zp(i), "[GeV]"
        write(log,*) "Gamma_Z':", gamzp(i), "[GeV]"
        write(log,*) "Gamma_Z'/m_Z':", gamzp(i)/mass_zp(i)
        write(log,*) "gLu:", gZpu(1,i)
        write(log,*) "gRu:", gZpu(2,i)
        write(log,*) "gLd:", gZpd(1,i)
        write(log,*) "gRd:", gZpd(2,i)
        write(log,*) "gLl:", gZpl(1,i)
        write(log,*) "gRl:", gZpl(2,i)
        write(log,*) "gLn:", gZpn(1,i)
        write(log,*) "gRn:", gZpn(2,i)
        write(log,*) "gLt:", gZpt(1,i)
        write(log,*) "gRt:", gZpt(2,i)
        write(log,*) "gLb:", gZpb(1,i)
        write(log,*) "gRb:", gZpb(2,i)
        write(log,*) "gLl3:", gZpl3(1,i)
        write(log,*) "gRl3:", gZpl3(2,i)
        write(log,*) "gLn3:", gZpn3(1,i)
        write(log,*) "gRn3:", gZpn3(2,i)
      end if
    end do
    print*, '------'
  end if
end subroutine print_model

end module modelling
