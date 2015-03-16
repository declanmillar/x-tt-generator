! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! introduced a new common for calculating more precise distributions
! common/reslocal/resl(20),standdevl(20)


! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------
!  start vegas section

block data
  implicit double precision(a-h,o-z)
!   makes default parameter assignments for vegas
  common/bveg1/ncall,itmx,nprn,ndev,xl(100),xu(100),acc
  common/bveg2/it,ndo,si,swgt,schi,xi(50,100)
  common/bveg3/alph,ndmx,mds
  data ncall/1000/,itmx/10/,nprn/0/,acc/1.d-2/, &
  xl/100*0.d0/,xu/100*1.d0/, &
  alph/1.5d0/,ndmx/50/,mds/1/,ndev/6/, &
  ndo/1/,xi/5000*1.d0/,it/0/,si,swgt,schi/3*0.d0/
end


subroutine vegas(ndim,fxn,avgi,sd,chi2a)

!     !  subroutine performs ndim-dimensional monte carlo integ'n
!      - by g.p. lepage    sept 1976/(rev)aug 1979
!      - algorithm described in j comp phys 27,192(1978)

  implicit double precision(a-h,o-z)
  common/bveg1/ncall,itmx,nprn,ndev,xl(100),xu(100),acc
  common/bveg2/it,ndo,si,swgt,schi,xi(50,100)
  common/bveg3/alph,ndmx,mds
  common/bveg4/calls,ti,tsi
  common/reslocal/resl(20),standdevl(20)
  dimension d(50,100),di(50,100),xin(50),r(50),dx(100),ia(100), &
  kg(100),dt(100),x(100)
  dimension rand(100)
  data one/1.d0/

  ndo=1
  do 1 j=1,ndim
    xi(1,j)=one
  1 END DO

  entry vegas1(ndim,fxn,avgi,sd,chi2a)
!          initializes cummulative variables, but not grid
  it=0
  si=0.d0
  swgt=si
  schi=si

  entry vegas2(ndim,fxn,avgi,sd,chi2a)
!         - no initialization
  nd=ndmx
  ng=1
  if(mds == 0) go to 2
  ng=(ncall/2.d0)**(1.d0/ndim)
  mds=-1
  if((2*ng-ndmx) < 0) go to 2
  mds=1
  npg=ng/ndmx+1
  nd=ng/npg
  ng=npg*nd
  2 k=ng**ndim
  npg=ncall/k
  if(npg < 2) npg=2
  calls=npg*k
  dxg=one/ng
  dv2g=(calls*dxg**ndim)**2/npg/npg/(npg-one)
  xnd=nd
  ndm=nd-1
  dxg=dxg*xnd
  xjac=one/calls
  do 3 j=1,ndim
    dx(j)=xu(j)-xl(j)
    xjac=xjac*dx(j)
  3 END DO

!   rebin, preserving bin density
  if(nd == ndo) go to 8
  rc=ndo/xnd
  do 7 j=1,ndim
    k=0
    xn=0.d0
    dr=xn
    i=k
    4 k=k+1
    dr=dr+one
    xo=xn
    xn=xi(k,j)
    5 if(rc > dr) go to 4
    i=i+1
    dr=dr-rc
    xin(i)=xn-(xn-xo)*dr
    if(i < ndm) go to 5
    do 6 i=1,ndm
      xi(i,j)=xin(i)
    6 END DO
    xi(nd,j)=one
  7 END DO
  ndo=nd

  8 if(nprn >= 0) write(ndev,200) ndim,calls,it,itmx,acc,nprn, &
  alph,mds,nd,(xl(j),xu(j),j=1,ndim)

  entry vegas3(ndim,fxn,avgi,sd,chi2a)
!         - main integration loop
  9 it=it+1
  ti=0.d0
  tsi=ti
  do 10 j=1,ndim
    kg(j)=1
    do 10 i=1,nd
      d(i,j)=ti
      di(i,j)=ti
  10 END DO

  11 fb=0.d0
  f2b=fb
  k=0
  12 k=k+1
  call randa(ndim,rand)
  wgt=xjac
  do 15 j=1,ndim
    xn=(kg(j)-rand(j))*dxg+one
    ia(j)=xn
    if(ia(j) > 1) go to 13
    xo=xi(ia(j),j)
    rc=(xn-ia(j))*xo
    go to 14
    13 xo=xi(ia(j),j)-xi(ia(j)-1,j)
    rc=xi(ia(j)-1,j)+(xn-ia(j))*xo
    14 x(j)=xl(j)+rc*dx(j)
    wgt=wgt*xo*xnd
  15 END DO

  f=wgt
  f=f*fxn(x,wgt)
  f2=f*f
  fb=fb+f
  f2b=f2b+f2
  do 16 j=1,ndim
    di(ia(j),j)=di(ia(j),j)+f
    if(mds >= 0) d(ia(j),j)=d(ia(j),j)+f2
  16 END DO
  if(k < npg) go to 12

  f2b=sqrt(f2b*npg)
  f2b=(f2b-fb)*(f2b+fb)
  ti=ti+fb
  tsi=tsi+f2b
  if(mds >= 0) go to 18
  do 17 j=1,ndim
    d(ia(j),j)=d(ia(j),j)+f2b
  17 END DO
  18 k=ndim
  19 kg(k)=mod(kg(k),ng)+1
  if(kg(k) /= 1) go to 11
  k=k-1
  if(k > 0) go to 19

!   compute final results for this iteration
  tsi=tsi*dv2g
  ti2=ti*ti
! ccccccccccccccc modifica 1: permette integrale nullo  cccccccccccccccc`

  if(tsi == 0.) then
    wgt=0.d0
    si=0.d0
    swgt=0.d0
    schi=0.d0
    avgi=0.d0
    chi2a=0.d0
    sd=0.d0
    if(nprn < 0) then
      continue
    else
      tsi=sqrt(tsi)
      write(ndev,201) it,ti,tsi,avgi,sd,chi2a
      if(nprn == 0) then
        continue
      else
        do j=1,ndim
          write(ndev,202) j,(xi(i,j),di(i,j),i=1,nd,nprn)
        end do
      end if
    end if
    return
  end if

! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc`

! protect against accidental single precision negative values
  xxx=tsi/ti**2
  if(xxx > -1.d-4) f2b=abs(f2b)

  wgt=one/tsi
  si=si+ti*wgt
  swgt=swgt+wgt
  schi=schi+ti2*wgt
  avgi=si/swgt
  chi2a=(schi-si*avgi)/(it-.9999d0)
  sd=sqrt(one/swgt)

  if(nprn < 0) go to 21
  tsi=sqrt(tsi)
  write(ndev,201) it,ti,tsi,avgi,sd,chi2a
  resl(it)=ti
  standdevl(it)=tsi
  if(nprn == 0) go to 21
  do 20 j=1,ndim
    write(ndev,202) j,(xi(i,j),di(i,j),i=1,nd,nprn)
  20 END DO

!   refine grid
  21 do 23 j=1,ndim
    xo=d(1,j)
    xn=d(2,j)
    d(1,j)=(xo+xn)/2.d0
    dt(j)=d(1,j)
    do 22 i=2,ndm
      d(i,j)=xo+xn
      xo=xn
      xn=d(i+1,j)
      d(i,j)=(d(i,j)+xn)/3.d0
      dt(j)=dt(j)+d(i,j)
    22 END DO
    d(nd,j)=(xn+xo)/2.d0
    dt(j)=dt(j)+d(nd,j)
  23 END DO

  do 28 j=1,ndim
    rc=0.d0
    do 24 i=1,nd
      r(i)=0.d0
      if(d(i,j) <= 0.d0) go to 24
      xo=dt(j)/d(i,j)
      r(i)=((xo-one)/xo/log(xo))**alph
      rc=rc+r(i)
    24 END DO
    rc=rc/xnd
    k=0
    xn=0.d0
    dr=xn
    i=k
    25 k=k+1
    dr=dr+r(k)
    xo=xn
    xn=xi(k,j)
    26 if(rc > dr) go to 25
    i=i+1
    dr=dr-rc
    xin(i)=xn-(xn-xo)*dr/r(k)
    if(i < ndm) go to 26
    do 27 i=1,ndm
      xi(i,j)=xin(i)
    27 END DO
    xi(nd,j)=one
  28 END DO

  if(it < itmx .AND. acc*abs(avgi) < sd) go to 9
  200 format(/35h input parameters for vegas:  ndim=,i3,8h  ncall=,f8.0 &
  /28x,5h  it=,i5,7h  itmx=,i5/28x,6h  acc=,g9.3 &
  /28x,7h  nprn=,i3,7h  alph=,f5.2/28x,6h  mds=,i3,6h   nd=,i4 &
  /28x,10h  (xl,xu)=,(t40,2h( ,g12.6,3h , ,g12.6,2h )))
  201 format(///21h integration by vegas//14h iteration no.,i3, &
  &   14h:   integral =,g14.8/21x,10hstd dev  =,g10.4/ &
  &   34h accumulated results:   integral =,g14.8/ &
  &   24x,10hstd dev  =,g10.4/24x,17hchi**2 per it'n =,g10.4)
  202 format(/15h data for axis ,i2/25h    x       delta i       , &
  &   24h   x       delta i       ,18h   x       delta i &
  /(1h ,f7.6,1x,g11.4,5x,f7.6,1x,g11.4,5x,f7.6,1x,g11.4))
  return
end subroutine vegas

!***************************************************************
! load vegas data if desired
subroutine load_vegas(ndim,name)
  implicit double precision(a-h,o-z)
  common/bveg2/it,ndo,si,swgt,schi,xi(50,100)
  character(30) :: name
!      open(unit=13,file=name,status='old',shared)
  read(13,210) it,ndo,si,swgt,schi
  do 190 j=1,ndim
    read(13,*)  (xi(i,j),i=1,50)
  190 END DO
  210 format(2i8,3z16)
  close(unit=13)
  return
end subroutine load_vegas

!*******************************************************************
! store vegas data for possible later use
subroutine store_vegas(ndim,name)
  implicit double precision(a-h,o-z)
  common/bveg2/it,ndo,si,swgt,schi,xi(50,100)
  character(30) :: name
  open(unit=12,file=name,status='new')
  write(12,210) it,ndo,si,swgt,schi
  do 190 j=1,ndim
    write(12,*)  (xi(i,j),i=1,50)
  190 END DO
  210 format(2i8,3z16)
  close(unit=12)
  return
end subroutine store_vegas

subroutine randa(n,rand)

!   subroutine generates uniformly distributed random no's x(i),i=1,n
  implicit double precision(a-h,o-z)
  common/rndm/iseed
  dimension rand(n)
  do 1 i=1,n
    rand(i)=ran2(iseed)
  1 END DO
  return
end subroutine randa

!--------- end the vegas section -------------------------------------
