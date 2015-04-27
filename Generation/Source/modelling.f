module modelling

  use configuration, only: use_nwa, model_name, verbose, nloops, lambdaqcd4

  implicit none

  ! SM couplings
  real :: gw, gwwa, gwwZ
  real :: gal(2),gau(2),gad(2),gwf(2)
  real :: gZn(2),gZl(2),gZu(2),gZd(2),g1(2)
  real :: gwwh,gZZh,ghhh,gwwhh,gZZhh,ghhhh
  complex*16 :: gchf(2,12)
  real :: gg(2), g

  ! fermion masses and widths
  real :: fmass(12), fwidth(12)

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
  real, parameter :: eWidth = 0.d0,   muwidth = 0.d0, tauwidth = 0.d0
  real, parameter :: nuewidth = 0.d0, numuwidth = 0.d0, nutauwidth = 0.d0

  ! SM boson masses
  real, parameter :: rm_w = 80.23d0,rm_z = 91.19d0,gamma_w = 2.08d0,gamma_z = 2.5d0
  real, parameter :: rm_a = 0d0, gamma_a = 0d0, rm_h = 125.d0, gamma_h = 0.31278d-2

  real :: Gamma_t = twidth

! Other SM parameters
  real, parameter :: a_em = 0.0078125, s2w = 0.2320d0

! Zprime parameters
  real :: rmZp(5),gamZp(5)
  real :: paramZp(5)
  real :: gp(5),gV_d(5),gA_d(5),gV_u(5),gA_u(5)
  real :: gZpd(2,5),gZpu(2,5)

  public :: initialise_standard_model
  public :: initialise_zprimes
  public :: widthZp
  public :: coupZpx


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

  if(use_NWA == 1)then
    fwidth(11)=1.d-5
  else
    fwidth(11)=twidth
  end if

  fwidth(12) = bwidth

  ! call SM HELAS couplings
  call coup1x(s2w,gw,gwwa,gwwZ)
  call coup2x(s2w,gal,gau,gad,gwf,gZn,gZl,gZu,gZd,g1)
  call coup3x(s2w,rm_Z,rm_h,gwwh,gZZh,ghhh,gwwhh,gZZhh,ghhhh)
  do i=1,12
    call coup4x(s2w,rm_Z,fmass(i),gchf(1,i))
  enddo

  ! QCD couplings (set to one, multiply correctly in zprime)
  g = 1.d0
  gg(1)=-g
  gg(2)=-g

  print*, lambdaqcd4, nloops

  print*, "...done"

end subroutine initialise_standard_model

subroutine initialise_zprimes 

  integer o_width(5), imodel_name, i

  print*, "Initialising standard model..."

  ! Extract model_name filename (Remove white space.)
  imodel_name = len(model_name)
  do while(model_name(imodel_name:imodel_name) == '')
    imodel_name = imodel_name-1
  end do

  ! read model file
  open(unit=42,file='Models/'//model_name(1:imodel_name)//'.mdl',status='old')
  read(42,*) rmZp
  read(42,*) gamZp
  read(42,*) gp
  read(42,*) paramZp
  read(42,*) gV_u
  read(42,*) gA_u
  read(42,*) gV_d
  read(42,*) gA_d
  close(42)

  ! Check whether width has been specified
  ! (If gamZp is zero, the function widthZp is used instead.)
  do i=1,5
    if ((gamZp(i) == 0d0) .AND. (rmZp(i) > 0d0)) then
      o_width(i) = 0
    else
      o_width(i) = 1
    end if
  enddo

  ! Calculate sequential Zp widths
  do i=1,5
    if (o_width(i) == 0) gamZp(i)= &
    widthZp(rmZp(i))
  end do

  ! convert from VA to LR couplings
  call coupZpx

  ! igw=0 ! don't include w width effects
  ! call topwid(fmass(11),wmass,fmass(12),wwidth,igw,fwidth(11))
  ! call printconstants
  return

  print*, "...done."
end subroutine initialise_zprimes

subroutine coupZpx

  ! input: vector and axial Zp couplings to up and down quarks
  ! output: left and right chiral couplings to up and down quarks

  integer i

	print*, "Converting Zprime couplings from AV to LR..."

  do i=1,5
      gZpd(1,i) = gp(i)*(gV_d(i)+gA_d(i))/2.d0
      gZpd(2,i) = gp(i)*(gV_d(i)-gA_d(i))/2.d0
      gZpu(1,i) = gp(i)*(gV_u(i)+gA_u(i))/2.d0
      gZpu(2,i) = gp(i)*(gV_u(i)-gA_u(i))/2.d0
  enddo

  print*, "...done."
   
  return
end subroutine coupZpx

function widthZp(rm_Zp)

  ! Calculates the width of the Zp in the SSM.
  ! Authors: stefano moretti and declan millar <d.millar@soton.ac.uk>

  implicit none

  ! implicit to explicit variable dump
  real :: widthzp
  real :: rm_zp
  real :: ctw
  real :: e
  real :: eq
  real :: gweak
  real :: gamt
  real :: gf
  integer :: i
  real :: pi
  real :: rmq
  real :: rmt
  real :: t3q
  real :: temp
  real :: temp1
  real :: temp2
  real :: alfas
  real :: a_s

  print*, "Calculating Zprime widths..."

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
  widthZp=0.d0
  do i=1,6
    rmq=0.d0
    if(i == 6)rmq=rmt
    if((i == 2) .OR. (i == 4) .OR. (i == 6))then

    ! Quarks
    ! u-quark.
      t3q=+1.d0/2.d0
      eq=+2.d0/3.d0
    else if((i == 1) .OR. (i == 3) .OR. (i == 5))then
    ! d-quark.
      t3q=-1.d0/2.d0
      eq=-1.d0/3.d0
    end if
    if(rm_Zp <= 2.d0*rmq)goto 123
    widthZp=widthZp+3.d0*rm_Z**2*rm_Zp*GF/24.d0/pi/sqrt(2.d0) &
    *sqrt(1.d0-4.d0*rmq**2/rm_Zp**2) &
    *((2.d0*t3q)**2 &
    *(1.d0-4.d0*rmq**2/rm_Zp**2) &
    +(2.d0*t3q-4.d0*eq*s2w)**2 &
    *(1.d0+2.d0*rmq**2/rm_Zp**2)) &
    *(1.d0+1.045d0*a_s/pi)
  !        cr=-eq*s2w
  !        cl=t3q-eq*s2w
  !        gv=cl+cr
  !        ga=cl-cr
  !        widthZp=widthZp+3.d0*e**2/12.d0/pi*rm_Zp/16.d0/s2w/ctw**2
  !     &               *sqrt(1.d0-4.d0*rmq**2/rm_Zp**2)
  !     &               *(4.d0*gv**2*(1.d0+2.d0*rmq**2/rm_Zp**2)
  !     &                +4.d0*ga**2*(1.d0-4.d0*rmq**2/rm_Zp**2))
  !     &               *(1.d0+1.045d0*a_s/pi)
    123 continue
  end do

  !       print *,'Z'' width due to quarks: ',widthZp,' [GeV]'

  ! Leptons
  temp=0.d0
  temp1=0.d0
  temp2=0.d0

  do i=1,6
    rmq=0.d0
    if((i == 2) .OR. (i == 4) .OR. (i == 6))then
    ! neutrino.
      t3q=+1.d0/2.d0
      eq=0.d0
    else if((i == 1) .OR. (i == 3) .OR. (i == 5))then
    ! lepton.
      t3q=-1.d0/2.d0
      eq=-1.d0
    end if
    widthZp=widthZp+1.d0*rm_Z**2*rm_Zp*GF/24.d0/pi/sqrt(2.d0) &
    *sqrt(1.d0-4.d0*rmq**2/rm_Zp**2) &
    *((2.d0*t3q)**2 &
    *(1.d0-4.d0*rmq**2/rm_Zp**2) &
    +(2.d0*t3q-4.d0*eq*s2w)**2 &
    *(1.d0+2.d0*rmq**2/rm_Zp**2))
  !        cr=-eq*s2w
  !        cl=t3q-eq*s2w
  !        gv=cl+cr
  !        ga=cl-cr
  !        widthZp=widthZp+1.d0*e**2/12.d0/pi*rm_Zp/16.d0/s2w/ctw**2
  !     &               *sqrt(1.d0-4.d0*rmq**2/rm_Zp**2)
  !     &               *(4.d0*gv**2*(1.d0+2.d0*rmq**2/rm_Zp**2)
  !     &                +4.d0*ga**2*(1.d0-4.d0*rmq**2/rm_Zp**2))

    temp=temp+1.d0*rm_Z**2*rm_Zp*GF/24.d0/pi/sqrt(2.d0) &
    *sqrt(1.d0-4.d0*rmq**2/rm_Zp**2) &
    *((2.d0*t3q)**2 &
    *(1.d0-4.d0*rmq**2/rm_Zp**2) &
    +(2.d0*t3q-4.d0*eq*s2w)**2 &
    *(1.d0+2.d0*rmq**2/rm_Zp**2))
    if((i == 2) .OR. (i == 4) .OR. (i == 6)) &
    temp1=temp1+1.d0*rm_Z**2*rm_Zp*GF/24.d0/pi/sqrt(2.d0) &
    *sqrt(1.d0-4.d0*rmq**2/rm_Zp**2) &
    *((2.d0*t3q)**2 &
    *(1.d0-4.d0*rmq**2/rm_Zp**2) &
    +(2.d0*t3q-4.d0*eq*s2w)**2 &
    *(1.d0+2.d0*rmq**2/rm_Zp**2))
    if((i == 1) .OR. (i == 3) .OR. (i == 5)) &
    temp2=temp2+1.d0*rm_Z**2*rm_Zp*GF/24.d0/pi/sqrt(2.d0) &
    *sqrt(1.d0-4.d0*rmq**2/rm_Zp**2) &
    *((2.d0*t3q)**2 &
    *(1.d0-4.d0*rmq**2/rm_Zp**2) &
    +(2.d0*t3q-4.d0*eq*s2w)**2 &
    *(1.d0+2.d0*rmq**2/rm_Zp**2))

  end do

  !       print *,'Z'' width due to quarks+leptons:',widthZp,' [GeV]'
  !       print *,'(so that due to leptons are:',temp,' [GeV])'
  !       print *,'(of which due to e/mu/tau:',temp2,' [GeV])'
  !       print *,'(of which due to their neutrinos are:',temp1,' [GeV])'
  print*, "...done."
  return
end function widthZp

end module modelling
