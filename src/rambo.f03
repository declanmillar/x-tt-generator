    subroutine rambo(iseed,n,et,xm,p,wt)
    use kinds
!------------------------------------------------------

!                       rambo

!    ra(ndom)  m(omenta)  b(eautifully)  o(rganized)

!    a democratic multi-particle phase space generator
!    authors:  s.d. ellis,  r. kleiss,  w.j. stirling
!    this is version 1.0 -  written by r. kleiss

!    n  = number of particles (>1, in this version <101)
!    et = total centre-of-mass energy
!    xm = particle masses ( dim=100 )
!    p  = particle momenta ( dim=(4,100) )
!    wt = weight of the event

!------------------------------------------------------
    implicit real(kind=default)(a-h,o-z)
    dimension xm(100),p(4,100),q(4,100),z(100),r(4), &
    b(3),p2(100),xm2(100),e(100),v(100),iwarn(5)
    data acc/1.d-14/,itmax/6/,ibegin/0/,iwarn/5*0/
    save twopi,po2log,z

! initialization step: factorials for the phase space weight
    if(ibegin /= 0) goto 103
    ibegin=1
    twopi=8.d0*datan(1.d0)
    po2log=dlog(twopi/4.d0)
    z(2)=po2log
    do 101 k=3,100
        z(k)=z(k-1)+po2log-2.d0*dlog(real(k-2,kind=default))
    101 END DO
    do 102 k=3,100
        z(k)=(z(k)-dlog(real(k-1,kind=default)))
    102 END DO

! check on the number of particles
    103 if(n > 1 .AND. n < 101) goto 104
    print 1001,n
    stop

! check whether total energy is sufficient; count nonzero masses
    104 xmt=0.d0
    nm=0
    do 105 i=1,n
        if(xm(i) /= 0.d0) nm=nm+1
        xmt=xmt+dabs(xm(i))
    105 END DO
    if(xmt <= et) goto 201
    print 1002,xmt,et
    stop

! the parameter values are now accepted

! generate n massless momenta in infinite phase space
    201 do 202 i=1,n
        c=2.d0*ran2(iseed)-1.d0
        s=dsqrt(dabs(1.d0-c*c))
        f=twopi*ran2(iseed)
        999 continue
        arg=ran2(iseed)*ran2(iseed)
        if(arg <= 0.d0)go to 999
        q(4,i)=-dlog(arg)
        q(3,i)=q(4,i)*c
        q(2,i)=q(4,i)*s*dcos(f)
        q(1,i)=q(4,i)*s*dsin(f)
    202 END DO

! calculate the parameters of the conformal transformation
    do 203 i=1,4
        r(i)=0.d0
    203 END DO
    do i=1,n
        do k=1,4
            r(k)=r(k)+q(k,i)
        end do
    end do
    rmas=dsqrt(dabs(r(4)**2-r(3)**2-r(2)**2-r(1)**2))
    if(rmas == 0.d0)goto 201
    do 205 k=1,3
        b(k)=-r(k)/rmas
    205 END DO
    g=r(4)/rmas
    a=1.d0/(1.d0+g)
    x=et/rmas

! transform the q's conformally into the p's
    do 207 i=1,n
        bq=b(1)*q(1,i)+b(2)*q(2,i)+b(3)*q(3,i)
        do 206 k=1,3
            p(k,i)=x*(q(k,i)+b(k)*(q(4,i)+a*bq))
        206 END DO
        p(4,i)=x*(g*q(4,i)+bq)
    207 END DO

! calculate weight and possible warnings
    wt=po2log
    if(n /= 2) wt=(2.d0*n-4.d0)*dlog(et)+z(n)
    if(wt >= -180.d0) goto 208
    if(iwarn(1) <= 5) print 1004,wt
    iwarn(1)=iwarn(1)+1
    208 if(wt <= 174.d0) goto 209
    if(iwarn(2) <= 5) print 1005,wt
    iwarn(2)=iwarn(2)+1

! return for weighted massless momenta
    209 if(nm /= 0) goto 210
    wt=dexp(wt)
    return

! massive particles: rescale the momenta by a factor x
    210 xmax=dsqrt(1.d0-(xmt/et)**2)
    do 301 i=1,n
        xm2(i)=xm(i)**2
        p2(i)=p(4,i)**2
    301 END DO
    iter=0
    x=xmax
    accu=et*acc
    302 f0=-et
    g0=0.d0
    x2=x*x
    do 303 i=1,n
        e(i)=dsqrt(xm2(i)+x2*p2(i))
        f0=f0+e(i)
        g0=g0+p2(i)/e(i)
    303 END DO
    if(dabs(f0) <= accu) goto 305
    iter=iter+1
    if(iter <= itmax) goto 304
    print 1006,itmax
    goto 305
    304 x=x-f0/(x*g0)
    goto 302
    305 do 307 i=1,n
        v(i)=x*p(4,i)
        do 306 k=1,3
            p(k,i)=x*p(k,i)
        306 END DO
        p(4,i)=e(i)
    307 END DO

! calculate the mass-effect weight factor
    wt2=1.d0
    wt3=0.d0
    do 308 i=1,n
        wt2=wt2*v(i)/e(i)
        wt3=wt3+v(i)**2/e(i)
    308 END DO
    wtm=(2.d0*n-3.d0)*dlog(x)+dlog(dabs(wt2/wt3*et))

! return for  weighted massive momenta
    wt=wt+wtm
    if(wt >= -180.d0) goto 309
    if(iwarn(3) <= 5) print 1004,wt
    iwarn(3)=iwarn(3)+1
    309 if(wt <= 174.d0) goto 310
    if(iwarn(4) <= 5) print 1005,wt
    iwarn(4)=iwarn(4)+1
    310 wt=dexp(wt)
    return

    1001 format(' rambo fails: # of particles =',i5,' is not allowed')
    1002 format(' rambo fails: total mass =',d15.6,' is not', &
    ' smaller than total energy =',d15.6)
    1004 format(' rambo warns: weight = exp(',f20.9,') may underflow')
    1005 format(' rambo warns: weight = exp(',f20.9,') may  overflow')
    1006 format(' rambo warns:',i3,' iterations did not give the', &
    ' desired accuracy =',d15.6)
    end subroutine rambo

    function random(seed)
    use kinds
!     -----------------
! ref.: k. park and k.w. miller, comm. of the acm 31 (1988) p.1192
! use seed = 1 as first value.

    implicit integer(a-z)
    double precision :: minv,random
    save
    parameter(m=2147483647,a=16807,q=127773,r=2836)
    parameter(minv=0.46566128752458d-09)
    hi = seed/q
    lo = mod(seed,q)
    seed = a*lo - r*hi
    if(seed <= 0) seed = seed + m
    random = seed*minv
    end function random

    function rn(dummy)
    use kinds
    integer :: dummy,seed
    common/mag/ seed
    real(kind=default) :: rn,random
    rn=random(seed)
    return
    end function rn
