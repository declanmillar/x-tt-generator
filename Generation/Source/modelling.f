module modelling

  use configuration, only: use_nwa, model_name, verbose, nloops, lambdaqcd4

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
  real, parameter :: umass = 0.d0, cmass = 0.d0, tmass = 175.d0
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
  real, parameter :: rm_w = 80.23d0,rm_z = 91.19d0,gamma_w = 2.08d0,gamma_z = 2.5d0
  real, parameter :: rm_a = 0d0, gamma_a = 0d0, rm_h = 125.d0, gamma_h = 0.31278d-2

  real :: Gamma_t = twidth

! Other SM parameters
  real, parameter :: a_em = 0.0078125, s2w = 0.2320d0

! zprime parameters
  real :: mass_zp(5),gamZp(5)
  real :: paramZp(5)
  real :: gp(5),gV_d(5),gA_d(5),gV_u(5),gA_u(5), ga_l(5), gv_l(5), gv_nu(5), ga_nu(5)
  real :: gZpd(2,5),gZpu(2,5),gZpl(2,5),gZpn(2,5)

  public :: initialise_standard_model
  public :: initialise_zprimes
  public :: width_zprime_ssm
  public :: convert_zprime_couplings
  public :: width_zprime_benchmark


contains

subroutine initialise_standard_model

  integer i

	print*, "Initialising standard model..."

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
  call coup3x(s2w, rm_Z, rm_h, gwwh, gZZh, ghhh, gwwhh, gZZhh, ghhhh)
  do i = 1, 12
    call coup4x(s2w, rm_Z, fmass(i), gchf(1,i))
  enddo

  ! QCD couplings (set to one, multiply correctly in zprime)
  g = 1.d0
  gg(1) = -g
  gg(2) = -g

!   print*, lambdaqcd4, nloops

  print*, "...done"

end subroutine initialise_standard_model

subroutine initialise_zprimes 

  integer manual_width(5), imodel_name, i

  print*, "Initialising zprimes..."

  ! Extract model_name filename (Remove white space.)
  imodel_name = len(model_name)
  do while(model_name(imodel_name:imodel_name) == '')
    imodel_name = imodel_name - 1
  end do

  ! read model file
  open(unit = 42, file = 'Models/'//model_name(1:imodel_name)//'.mdl', status = 'old')
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
  close(42)
  print*, "Reading of model file complete."


  ! Check whether width has been specified
  ! (If gamZp is negative, the function widthZp is used instead.)
  do i = 1, 5
    if ((mass_zp(i) > 0.d0) .and. (gamzp(i) < 0.d0)) then
      manual_width(i) = 0
    else
      manual_width(i) = 1
    end if
  enddo

  ! Calculate sequential Zp widths
    do i = 1, 5
      if (manual_width(i) == 0) gamZp(i) = width_zprime_ssm(mass_zp(i))
    end do

  ! convert from VA to LR couplings
  call convert_zprime_couplings

  ! Calculate benchmark widths
  call width_zprime_benchmark

  ! igw=0 ! don't include w width effects
  ! call topwid(fmass(11),wmass,fmass(12),wwidth,igw,fwidth(11))
  ! call printconstants
  return

  print*, "...done."
end subroutine initialise_zprimes

subroutine convert_zprime_couplings

  ! input: vector and axial Zp couplings to up and down quarks
  ! output: left and right chiral couplings to up and down quarks

  integer i

	print*, "Converting zprime couplings from AV to LR..."

  do i = 1, 5
      gZpd(1,i) = gp(i)*(gv_d(i) + ga_d(i))/2.d0
      gZpd(2,i) = gp(i)*(gv_d(i) - ga_d(i))/2.d0
      gZpu(1,i) = gp(i)*(gv_u(i) + ga_u(i))/2.d0
      gZpu(2,i) = gp(i)*(gv_u(i) - ga_u(i))/2.d0
      gZpl(1,i) = gp(i)*(gv_l(i) + ga_l(i))/2.d0
      gZpl(2,i) = gp(i)*(gv_l(i) - ga_l(i))/2.d0
      gZpn(1,i) = gp(i)*(gv_nu(i) + ga_nu(i))/2.d0
      gZpn(2,i) = gp(i)*(gv_nu(i) - ga_nu(i))/2.d0
  enddo

  print*, "...done."
   
  return
end subroutine convert_zprime_couplings

subroutine width_zprime_benchmark

  ! NOTE 23 Jan 2014 found a factor two wrong in the convention for gv and ga
  ! in the E6 models coming from the factor 2 in the feynman rule for Elenas paper
  ! All widths have been rescaled by 1/4
  ! calculates Zp width contributions from light fermions

  implicit none

  integer :: i, n
  real :: mq, ml, mzp
  real :: gv, ga
  real :: widthqq, widthll, widthqq_tmp, widthll_tmp
  real :: pi
  real :: a_s, alfas

  ! couplings.
  pi = dacos(-1.d0)

  do n = 1, 5
    widthqq = 0.d0
    mzp = mass_zp(n)
    if (mzp > 0.d0) then
      a_s = alfas(mzp, lambdaQCD4, nloops)
      ! quarks
      do i = 1, 6
        widthqq_tmp = 0.d0
        mq = qmass(i)
        if (mass_zp(n) <= 2.d0*mq) goto 123

        if (i == 1 .or. i == 3 .or. i == 5) then
          gv = gzpu(1,n) + gzpu(2,n)
          ga = gzpu(1,n) - gzpu(2,n)
        else if (i == 2 .or. i == 4 .or. i == 6) then
          gv = gzpd(1,n) + gzpd(2,n)
          ga = gzpd(1,n) - gzpd(2,n)
        end if        

        ! with QCD kfactor
        widthqq_tmp = 3.d0/48.d0/pi*mzp &
                    *sqrt(1.d0 - 4.d0*mq**2/mzp**2) &
                    *(gv**2*(1.d0 + 2.d0*mq**2/mzp**2) &
                    + ga**2*(1.d0 - 4.d0*mq**2/mzp**2)) &
                    *(1.d0 + 1.045d0*a_s/pi)

        widthqq = widthqq + widthqq_tmp
              
        123 continue
      end do

      ! leptons
      do i = 1, 6
        widthll_tmp = 0.d0
        ml = lmass(i)

        if (mzp <= 2.d0*ml) goto 456

        if (i == 1 .or. i == 3 .or. i == 5) then
          gv = gzpl(1,n) + gzpl(2,n)
          ga = gzpl(1,n) - gzpl(2,n)
        else if (i == 2 .or. i == 4 .or. i == 6) then
          gv = gzpn(1,n) + gzpn(2,n)
          ga = gzpn(1,n) - gzpn(2,n)
        end if        
          
        widthll_tmp = 1.d0/48.d0/pi*mzp &
                      *sqrt(1.d0 - 4.d0*mq**2/mzp**2) &
                      *(gv**2*(1.d0 + 2.d0*mq**2/mzp**2) &
                      +ga**2*(1.d0 - 4.d0*mq**2/mzp**2))

        widthll = widthll + widthll_tmp
        456 continue
      end do

      gamZp(n) = widthqq + widthll

      print*, 'ZPRIME WIDTHS'
      print*, 'Gamma(Zp(', n, ')->ff)=', gamZp(n),' [GeV]'
      print*, 'Gamma(Zp(', n, ')->ll)=', widthll,' [GeV]'
      print*, 'Gamma(Zp(', n, ')->qq)=', widthqq,' [GeV]'
    end if
  end do
  return
end subroutine width_zprime_benchmark

function width_zprime_ssm(rm_Zp)

  ! Calculates the width of the Zp in the SSM.
  ! Authors: stefano moretti and declan millar <d.millar@soton.ac.uk>

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
  real :: mq
  real :: rmt
  real :: t3q
  real :: temp
  real :: temp1
  real :: temp2
  real :: alfas
  real :: a_s

  print*, "Calculating zprime widths..."

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
  e=sqrt(s2w*8.d0*rm_Z*rm_Z*ctw*ctw*GF/sqrt(2.d0))
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
    width_zprime_ssm=width_zprime_ssm+3.d0*rm_Z**2*rm_Zp*GF/24.d0/pi/sqrt(2.d0) &
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
    width_zprime_ssm=width_zprime_ssm+1.d0*rm_Z**2*rm_Zp*GF/24.d0/pi/sqrt(2.d0) &
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

    temp=temp+1.d0*rm_Z**2*rm_Zp*GF/24.d0/pi/sqrt(2.d0) &
    *sqrt(1.d0-4.d0*mq**2/rm_Zp**2) &
    *((2.d0*t3q)**2 &
    *(1.d0-4.d0*mq**2/rm_Zp**2) &
    +(2.d0*t3q-4.d0*eq*s2w)**2 &
    *(1.d0+2.d0*mq**2/rm_Zp**2))
    if((i == 2) .OR. (i == 4) .OR. (i == 6)) &
    temp1=temp1+1.d0*rm_Z**2*rm_Zp*GF/24.d0/pi/sqrt(2.d0) &
    *sqrt(1.d0-4.d0*mq**2/rm_Zp**2) &
    *((2.d0*t3q)**2 &
    *(1.d0-4.d0*mq**2/rm_Zp**2) &
    +(2.d0*t3q-4.d0*eq*s2w)**2 &
    *(1.d0+2.d0*mq**2/rm_Zp**2))
    if((i == 1) .OR. (i == 3) .OR. (i == 5)) &
    temp2=temp2+1.d0*rm_Z**2*rm_Zp*GF/24.d0/pi/sqrt(2.d0) &
    *sqrt(1.d0-4.d0*mq**2/rm_Zp**2) &
    *((2.d0*t3q)**2 &
    *(1.d0-4.d0*mq**2/rm_Zp**2) &
    +(2.d0*t3q-4.d0*eq*s2w)**2 &
    *(1.d0+2.d0*mq**2/rm_Zp**2))

  end do

        print *,'Z'' width due to quarks+leptons:',width_zprime_ssm,' [GeV]'
  !       print *,'(so that due to leptons are:',temp,' [GeV])'
  !       print *,'(of which due to e/mu/tau:',temp2,' [GeV])'
  !       print *,'(of which due to their neutrinos are:',temp1,' [GeV])'
  print*, "...done."
  return
end function width_zprime_ssm

end module modelling
