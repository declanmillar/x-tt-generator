subroutine initialise_madgraph(o_NWA,model_name)
  ! sets up masses and coupling constants of particles

  use modelling

  implicit none

  ! arguments
  integer :: o_NWA
  character(50) :: model_name

  real ::  widthZp
  external widthZp

! other local variables
  integer :: i, o_width(5), imodel_name

! enter global fermion masses
  fmass(1) = emass
  fmass(2) = nuemass
  fmass(3) = umass
  fmass(4) = dmass
  fmass(5) = mumass
  fmass(6) = numumass
  fmass(7) = cmass
  fmass(8) = smass
  fmass(9) = taumass
  fmass(10)= nutaumass
  fmass(11)= tmass
  fmass(12)= bmass

  fwidth(1) = ewidth
  fwidth(2) = nuewidth
  fwidth(3) = uwidth
  fwidth(4) = dwidth
  fwidth(5) = muwidth
  fwidth(6) = numuwidth
  fwidth(7) = cwidth
  fwidth(8) = swidth
  fwidth(9) = tauwidth
  fwidth(10)= nutauwidth
  if(o_NWA == 1)then
    fwidth(11)=1.d-5
  else
    fwidth(11)=twidth
  end if
  fwidth(12)= bwidth

  ! call SM HELAS couplings
  call coup1x(s2w,gw,gwwa,gwwZ)
  call coup2x(s2w,gal,gau,gad,gwf,gZn,gZl,gZu,gZd,g1)
  call coup3x(s2w,rm_Z,rm_h,gwwh,gZZh,ghhh,gwwhh,gZZhh,ghhhh)
  do i=1,12
    call coup4x(s2w,rm_Z,fmass(i),gchf(1,i))
  enddo

  ! QCD couplings
  g = 1.d0
  gg(1)=-g
  gg(2)=-g

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
    widthZp(rm_W,rm_Z,rmZp(i),a_em,s2w,rlambdaQCD4,nloops)
  end do

  ! convert from VA to LR couplings
  call coupZpx

  ! igw=0 ! don't include w width effects
  ! call topwid(fmass(11),wmass,fmass(12),wwidth,igw,fwidth(11))
  ! call printconstants
  return
end subroutine initialise_madgraph

! subroutine printconstants

!   implicit none

! !     local

! !      integer i

!   write(*,'(a)') 'MADGRAPH Parameters'
!   write(*,'(a)') 'Boson masses and Widths:'
!   write(*,10) 'W',wmass,'W',wwidth,'Z',Zmass,'Z',Zwidth
!   write(*,10) 'A',Amass,'A',awidth,'H',hmass,'H',hwidth
!   write(*,'(a)') 'Quark masses and Widths:'
!   write(*,10) 't',fmass(11),'t',fwidth(11), &
!   'b',fmass(12),'b',fwidth(12)


!   10 format(a,1x,5Hmass=,f6.2,1H:,2x,a,1x,6Hwidth=,f6.2,1H:,2x, &
!   a,1x,5Hmass=,f6.2,1H:,2x,a,1x,6Hwidth=,f6.2)

! end subroutine printconstants

! subroutine topwid(rmt,rmw,rmb,rgw,igw,rgt)
! !*************************************************************************
! !     the total weak decay width of the top quark, including
! !     the effects of bottom mass and, if igw=1,  a finite w width.
! !     from james stirling 6-10-94
! !*************************************************************************
!   implicit complex*16(a-h,o-Z)
!   real :: rmt,rmb,rmw,xw,xb,rgw,rgt

!   pi=4.*atan(1.d0)
!   gf=1.16637d-05
!   gw=sqrt(rmw**2*gf/dsqrt(2.d0))
! !                            flavour & colour
!   xb=rmb/rmt
!   xw=rmw/rmt
!   if(igw == 1) goto 10
!   if(rmt <= (rmw+rmb)) then
!     write(6,*)'WARNING: mt < mw + mb !!!!'
!     stop
!   endif
!   rgt = gf*rmt**3/8d0/pi/dsqrt(2d0) &
!   * dsqrt( (1d0-(xw+xb)**2)*(1d0-(xw-xb)**2) ) &
!   * ( (1d0-xb**2)**2 + (1d0+xb**2)*xw**2 - 2d0*xw**4 )
!   return
!   10 continue
!   rm=xb**2
!   om=1.+rm-dcmplx(rmw**2,rmw*rgw)/rmt**2
!   y1=om+sqrt(om*om-4.*rm)
!   y0=om-sqrt(om*om-4.*rm)
!   Z1=2.
!   Z0=2.*sqrt(rm)

!   d0=(-y0**8+3.*y0**7*rm+3.*y0**7-8.*y0**6*rm-12.*y0**5*rm** &
!   &  2-12.*y0**5*rm+96.*y0**4*rm**2-48.*y0**3*rm**3-48.*y0**3* &
!   rm**2-128.*y0**2*rm**3+192.*y0*rm**4+192.*y0*rm**3-256.* &
!   rm**4)/(24.*y0**4*(y1-y0))
!   d1=(-y1**8+3.*y1**7*rm+3.*y1**7-8.*y1**6*rm-12.*y1**5*rm** &
!   &  2-12.*y1**5*rm+96.*y1**4*rm**2-48.*y1**3*rm**3-48.*y1**3* &
!   rm**2-128.*y1**2*rm**3+192.*y1*rm**4+192.*y1*rm**3-256.* &
!   rm**4)/(24.*y1**4*(y1-y0))
!   a4=(32.*rm**4*(y1-y0))/(3.*y1*y0*(y1-y0))
!   a3=(8.*rm**3*(-3.*y1**2*y0*rm-3.*y1**2*y0+4.*y1**2*rm+3.* &
!   y1*y0**2*rm+3.*y1*y0**2-4.*y0**2*rm))/(3.*y1**2*y0**2*(y1 &
!   -y0))
!   a2=(8.*rm**3*(2.*y1**3*y0**2-3.*y1**3*y0*rm-3.*y1**3*y0+4. &
!   *y1**3*rm-2.*y1**2*y0**3+3.*y1*y0**3*rm+3.*y1*y0**3-4.*y0 &
!   **3*rm))/(3.*y1**3*y0**3*(y1-y0))
!   a1=(2.*rm**2*(3.*y1**4*y0**3*rm+3.*y1**4*y0**3+8.*y1**4*y0 &
!   **2*rm-12.*y1**4*y0*rm**2-12.*y1**4*y0*rm+16.*y1**4*rm**2 &
!   -3.*y1**3*y0**4*rm-3.*y1**3*y0**4-8.*y1**2*y0**4*rm+12.* &
!   y1*y0**4*rm**2+12.*y1*y0**4*rm-16.*y0**4*rm**2))/(3.*y1** &
!   &  4*y0**4*(y1-y0))
!   b0=(y1**3-3.*y1**2*rm-3.*y1**2+8.*y1*rm-y0**3+3.*y0**2*rm+ &
!   &  3.*y0**2-8.*y0*rm)/(24.*(y1-y0))
!   b1=(y1+y0-3.*rm-3.)/24.
!   b2=1./24.

!   rint=d0*log((Z1-y0)/(Z0-y0)) &
!   -d1*log((y1-Z1)/(y1-Z0)) &
!   -a4/3.*(1./Z1**3-1./Z0**3) &
!   -a3/2.*(1./Z1**2-1./Z0**2) &
!   -a2   *(1./Z1   -1./Z0   ) &
!   +a1*log(Z1/Z0) &
!   +b0   *(Z1   -Z0   ) &
!   +b1/2.*(Z1**2-Z0**2) &
!   +b2/3.*(Z1**3-Z0**3)

!   gw4=gw**4

! ! total width includes flavour & colour factors
!   rgt=rmt**3/(rmw*rgw)*gw4/(8.*pi**3)*aimag(rint)
!   rgt=9.*rgt
!   return
! end subroutine topwid
