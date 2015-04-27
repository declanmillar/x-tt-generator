function alfas(mu,lam,nloop)

! EWNG version from EERAD

! two loop strong coupling constant at scale rq
! alfas in terms of lambda of four flavors for one or two loops.
! matching achieved using renorm group eqn. approximately
! above mu=mb,mu=mt

implicit none

integer :: nloop
real :: mu,lam,mc,mb,mt
real :: b4,b4p,b5,b5p,b6,b6p,one,two
real :: alfas,atinv,abinv,atinv1,abinv1,asinv,xqc,xqb,xqt,xb,xt
parameter(one=1.d0,two=2.d0)
parameter(b4=1.326291192d0,b5=1.220187897d0,b6=1.114084601d0)
parameter(b4p=0.490197225d0,b5p=0.401347248d0,b6p=0.295573466d0)
! b=(33.d0-2.d0*nf)/6.d0/pi
! bp=(153.d0-19.d0*nf)/(2.d0*pi*(33.d0-2.d0*nf))
parameter(mc=1.5d0,mb=4.5d0,mt=175.d0)

if (mu < mc) then
    write(6,*) 'unimplemented, mu too low',mu
    alfas=0.d0
    return
endif

xb=log(mb/lam)
abinv=b4*xb

if (mu < mb) then
    xqc=log(mu/lam)
    asinv=b4*xqc
    alfas=one/asinv
elseif (mu < mt) then
    xqb=log(mu/mb)
    asinv=abinv+b5*xqb
    alfas=one/asinv
else
    xqt=log(mu/mt)
    xt=log(mt/mb)
    atinv=abinv+b5*xt
    asinv=atinv+b6*xqt
    alfas=one/asinv
endif
       
if (nloop == 1) return

abinv1=abinv/(one-b4p*log(two*xb)/abinv)
if (mu < mb) then
    asinv=asinv/(one-b4p*log(two*xqc)/asinv)
    alfas=one/asinv
elseif (mu < mt) then
    asinv=abinv1+b5*xqb+b5p*log((asinv+b5p)/(abinv+b5p))
    alfas=one/asinv
else
    atinv1=abinv1+b5*xt+b5p*log((b5p+atinv)/(b5p+abinv))
    asinv=atinv1+b6*xqt+b6p*log((b6p+asinv)/(atinv+b6p))
    alfas=one/asinv
endif
return
end function alfas
