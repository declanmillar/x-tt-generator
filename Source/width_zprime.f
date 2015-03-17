function widthZp(rmW,rmZ,rmZp,a_em,s2w,rlambdaQCD4,nloop)

  ! Calculates the width of the Zp in the SSM.
  ! Authors: stefano moretti and declan millar <d.millar@soton.ac.uk>

  implicit none

  common/fermions/ fmass,     fwidth
  real :: fmass(12), fwidth(12)


  ! implicit to explicit variable dump
  real :: ME2
  real :: widthzp
  real :: rmw
  real :: rmz
  real :: rmzp
  real :: a_em
  real :: s2w
  real :: rlambdaqcd4
  integer :: nloop
  real :: a_s
  real :: ctw
  real :: e
  real :: eq
  real :: g
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

  rmt=fmass(11)
  gamt=fwidth(11)
  ! couplings.
  pi=dacos(-1.d0)
  ctw=sqrt(1.d0-s2w)
  e=sqrt(4.d0*pi*a_em)
  g=e/sqrt(s2w)
  a_s=alfas(rmZp,rlambdaQCD4,nloop,4)
  GF=1.16639D-5
  ! renormalise e.
  e=sqrt(s2w*8.d0*rmZ*rmZ*ctw*ctw*GF/sqrt(2.d0))
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
    if(rmZp <= 2.d0*rmq)goto 123
    widthZp=widthZp+3.d0*rmZ**2*rmZp*GF/24.d0/pi/sqrt(2.d0) &
    *sqrt(1.d0-4.d0*rmq**2/rmZp**2) &
    *((2.d0*t3q)**2 &
    *(1.d0-4.d0*rmq**2/rmZp**2) &
    +(2.d0*t3q-4.d0*eq*s2w)**2 &
    *(1.d0+2.d0*rmq**2/rmZp**2)) &
    *(1.d0+1.045d0*a_s/pi)
  !        cr=-eq*s2w
  !        cl=t3q-eq*s2w
  !        gv=cl+cr
  !        ga=cl-cr
  !        widthZp=widthZp+3.d0*e**2/12.d0/pi*rmZp/16.d0/s2w/ctw**2
  !     &               *sqrt(1.d0-4.d0*rmq**2/rmZp**2)
  !     &               *(4.d0*gv**2*(1.d0+2.d0*rmq**2/rmZp**2)
  !     &                +4.d0*ga**2*(1.d0-4.d0*rmq**2/rmZp**2))
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
    widthZp=widthZp+1.d0*rmZ**2*rmZp*GF/24.d0/pi/sqrt(2.d0) &
    *sqrt(1.d0-4.d0*rmq**2/rmZp**2) &
    *((2.d0*t3q)**2 &
    *(1.d0-4.d0*rmq**2/rmZp**2) &
    +(2.d0*t3q-4.d0*eq*s2w)**2 &
    *(1.d0+2.d0*rmq**2/rmZp**2))
  !        cr=-eq*s2w
  !        cl=t3q-eq*s2w
  !        gv=cl+cr
  !        ga=cl-cr
  !        widthZp=widthZp+1.d0*e**2/12.d0/pi*rmZp/16.d0/s2w/ctw**2
  !     &               *sqrt(1.d0-4.d0*rmq**2/rmZp**2)
  !     &               *(4.d0*gv**2*(1.d0+2.d0*rmq**2/rmZp**2)
  !     &                +4.d0*ga**2*(1.d0-4.d0*rmq**2/rmZp**2))

    temp=temp+1.d0*rmZ**2*rmZp*GF/24.d0/pi/sqrt(2.d0) &
    *sqrt(1.d0-4.d0*rmq**2/rmZp**2) &
    *((2.d0*t3q)**2 &
    *(1.d0-4.d0*rmq**2/rmZp**2) &
    +(2.d0*t3q-4.d0*eq*s2w)**2 &
    *(1.d0+2.d0*rmq**2/rmZp**2))
    if((i == 2) .OR. (i == 4) .OR. (i == 6)) &
    temp1=temp1+1.d0*rmZ**2*rmZp*GF/24.d0/pi/sqrt(2.d0) &
    *sqrt(1.d0-4.d0*rmq**2/rmZp**2) &
    *((2.d0*t3q)**2 &
    *(1.d0-4.d0*rmq**2/rmZp**2) &
    +(2.d0*t3q-4.d0*eq*s2w)**2 &
    *(1.d0+2.d0*rmq**2/rmZp**2))
    if((i == 1) .OR. (i == 3) .OR. (i == 5)) &
    temp2=temp2+1.d0*rmZ**2*rmZp*GF/24.d0/pi/sqrt(2.d0) &
    *sqrt(1.d0-4.d0*rmq**2/rmZp**2) &
    *((2.d0*t3q)**2 &
    *(1.d0-4.d0*rmq**2/rmZp**2) &
    +(2.d0*t3q-4.d0*eq*s2w)**2 &
    *(1.d0+2.d0*rmq**2/rmZp**2))

  end do

!       print *,'Z'' width due to quarks+leptons:',widthZp,' [GeV]'
!       print *,'(so that due to leptons are:',temp,' [GeV])'
!       print *,'(of which due to e/mu/tau:',temp2,' [GeV])'
!       print *,'(of which due to their neutrinos are:',temp1,' [GeV])'

  return
end function widthZp
