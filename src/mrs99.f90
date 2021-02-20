    subroutine mrs99(x,q,mode,upv,dnv,usea,dsea,str,chm,bot,glu)
    use vamp_kinds
!****************************************************************C
!       C
!     This is a package for the new **corrected** MRST parton    C
!     distributions. The format is similar to the previous       C
!     (1998) MRST series.                                        C
!       C
!     NOTE: 7 new sets are added here, corresponding to shifting C
!     the small x HERA data up and down by 2.5%, and by varying  C
!     the charm and strange distributions, and by forcing a      C
!     larger d/u ratio at large x.                               C
!       C
!     As before, x times the parton distribution is returned,    C
!     q is the scale in GeV, MSbar factorization is assumed,     C
!     and Lambda(MSbar,nf=4) is given below for each set.        C
!       C
!     NAMING SCHEME:                                             C
!                       C
!  mode  set    comment             L(4)/MeV  a_s(M_Z)  grid#1   C
!  ----  ---    -------             --------  -------   ------   C
!       C
!  1     COR01  central gluon, a_s    300      0.1175   0.00537  C
!  2     COR02  higher gluon          300      0.1175   0.00497  C
!  3     COR03  lower gluon           300      0.1175   0.00398  C
!  4     COR04  lower a_s             229      0.1125   0.00585  C
!  5     COR05  higher a_s            383      0.1225   0.00384  C
!  6     COR06  quarks up             303.3    0.1178   0.00497  C
!  7     COR07  quarks down           290.3    0.1171   0.00593  C
!  8     COR08  strange up            300      0.1175   0.00524  C
!  9     COR09  strange down          300      0.1175   0.00524  C
!  10    C0R10  charm up              300      0.1175   0.00525  C
!  11    COR11  charm down            300      0.1175   0.00524  C
!  12    COR12  larger d/u            300      0.1175   0.00515  C
!                       C
!      The corresponding grid files are called cor01.dat etc.    C
!               C
!      The reference is:                                         C
!      A.D. Martin, R.G. Roberts, W.J. Stirling, R.S Thorne      C
!      Univ. Durham preprint DTP/99/64, hep-ph/9907231 (1999)    C
!                                                                C
!      Comments to : W.J.Stirling@durham.ac.uk                   C
!                                                                C
!       C
!****************************************************************C
    implicit real(a-h,o-z)
    data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
    q2=q*q
    if(q2 < qsqmin .OR. q2 > qsqmax) print 99
    if(x < xmin .OR. x > xmax)       print 98
    if(mode == 1) then
        call mrs991(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu)
    elseif(mode == 2) then
        call mrs992(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu)
    elseif(mode == 3) then
        call mrs993(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu)
    elseif(mode == 4) then
        call mrs994(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu)
    elseif(mode == 5) then
        call mrs995(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu)
    elseif(mode == 6) then
        call mrs996(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu)
    elseif(mode == 7) then
        call mrs997(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu)
    elseif(mode == 8) then
        call mrs998(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu)
    elseif(mode == 9) then
        call mrs999(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu)
    elseif(mode == 10) then
        call mrs9910(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu)
    elseif(mode == 11) then
        call mrs9911(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu)
    elseif(mode == 12) then
        call mrs9912(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu)
    endif
    99 format('  WARNING:  Q^2 VALUE IS OUT OF RANGE   ')
    98 format('  WARNING:   X  VALUE IS OUT OF RANGE   ')
    return
    end subroutine mrs99

    subroutine mrs991(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
    use vamp_kinds
    implicit real(a-h,o-z)
    parameter(nx=49,nq=37,ntenth=23,np=8)
    real(kind=default) :: f(np,nx,nq+1),qq(nq),xx(nx),g(np),n0(np)
    data xx/1d-5,2d-5,4d-5,6d-5,8d-5, &
    &             1d-4,2d-4,4d-4,6d-4,8d-4, &
    &             1d-3,2d-3,4d-3,6d-3,8d-3, &
    &             1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2, &
    .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0, &
    .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0, &
    .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0, &
    .8d0,.9d0,1d0/
    data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1, &
    &         1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2, &
    &         1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4, &
    &         1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6, &
    &         1.8d6,3.2d6,5.6d6,1d7/
    data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
    data n0/3,4,5,9,9,9,9,9/
    data init/0/
    save
    xsave=x
    q2save=qsq
    if(init /= 0) goto 10
    open(unit=1,file='cor01.dat',status='old')
    do n=1,nx-1
        do m=1,nq
            read(1,50)f(1,n,m),f(2,n,m),f(3,n,m),f(4,n,m), &
            f(5,n,m),f(7,n,m),f(6,n,m),f(8,n,m)
        ! notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
            do 25 i=1,np
                f(i,n,m)=f(i,n,m)/(1d0-xx(n))**n0(i)
            25 END DO
        END DO
    end do
    do j=1,ntenth-1
        xx(j)=log10(xx(j)/xx(ntenth))+xx(ntenth)
        do i=1,8
            if(i == 5 .OR. i == 7) continue
            do 30 k=1,nq
                f(i,j,k)=log10(f(i,j,k)/f(i,ntenth,k))+f(i,ntenth,k)
            30 END DO
        END DO
    end do
    50 format(8f10.5)
    do i=1,np
        do m=1,nq
            f(i,nx,m)=0d0
        END DO
    end do
    init=1
    10 continue
    if(x < xmin) x=xmin
    if(x > xmax) x=xmax
    if(qsq < qsqmin)      qsq=qsqmin
    if(qsq > qsqmax)      qsq=qsqmax
    xxx=x
    if(x < xx(ntenth)) xxx=log10(x/xx(ntenth))+xx(ntenth)
    n=0
    70 n=n+1
    if(xxx > xx(n+1)) goto 70
    a=(xxx-xx(n))/(xx(n+1)-xx(n))
    m=0
    80 m=m+1
    if(qsq > qq(m+1)) goto 80
    b=(qsq-qq(m))/(qq(m+1)-qq(m))
    do 60 i=1,np
        g(i)= (1d0-a)*(1d0-b)*f(i,n,m)   + (1d0-a)*b*f(i,n,m+1) &
        +       a*(1d0-b)*f(i,n+1,m) +       a*b*f(i,n+1,m+1)
        if(n >= ntenth) goto 65
        if(i == 5 .OR. i == 7) goto 65
        fac=(1d0-b)*f(i,ntenth,m)+b*f(i,ntenth,m+1)
        g(i)=fac*10d0**(g(i)-fac)
        65 continue
        g(i)=g(i)*(1d0-x)**n0(i)
    60 END DO
    upv=g(1)
    dnv=g(2)
    usea=g(4)
    dsea=g(8)
    str=g(6)
    chm=g(5)
    glu=g(3)
    bot=g(7)
    x=xsave
    qsq=q2save
    return
    end subroutine mrs991

    subroutine mrs992(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
    use vamp_kinds
    implicit real(a-h,o-z)
    parameter(nx=49,nq=37,ntenth=23,np=8)
    real(kind=default) :: f(np,nx,nq+1),qq(nq),xx(nx),g(np),n0(np)
    data xx/1d-5,2d-5,4d-5,6d-5,8d-5, &
    &             1d-4,2d-4,4d-4,6d-4,8d-4, &
    &             1d-3,2d-3,4d-3,6d-3,8d-3, &
    &             1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2, &
    .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0, &
    .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0, &
    .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0, &
    .8d0,.9d0,1d0/
    data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1, &
    &         1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2, &
    &         1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4, &
    &         1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6, &
    &         1.8d6,3.2d6,5.6d6,1d7/
    data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
    data n0/3,4,5,9,9,9,9,9/
    data init/0/
    save
    xsave=x
    q2save=qsq
    if(init /= 0) goto 10
    open(unit=1,file='cor02.dat',status='old')
    do n=1,nx-1
        do m=1,nq
            read(1,50)f(1,n,m),f(2,n,m),f(3,n,m),f(4,n,m), &
            f(5,n,m),f(7,n,m),f(6,n,m),f(8,n,m)
        ! notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
            do 25 i=1,np
                f(i,n,m)=f(i,n,m)/(1d0-xx(n))**n0(i)
            25 END DO
        END DO
    end do
    do j=1,ntenth-1
        xx(j)=log10(xx(j)/xx(ntenth))+xx(ntenth)
        do i=1,8
            if(i == 5 .OR. i == 7) continue
            do 30 k=1,nq
                f(i,j,k)=log10(f(i,j,k)/f(i,ntenth,k))+f(i,ntenth,k)
            30 END DO
        END DO
    end do
    50 format(8f10.5)
    do i=1,np
        do m=1,nq
            f(i,nx,m)=0d0
        END DO
    end do
    init=1
    10 continue
    if(x < xmin) x=xmin
    if(x > xmax) x=xmax
    if(qsq < qsqmin)      qsq=qsqmin
    if(qsq > qsqmax)      qsq=qsqmax
    xxx=x
    if(x < xx(ntenth)) xxx=log10(x/xx(ntenth))+xx(ntenth)
    n=0
    70 n=n+1
    if(xxx > xx(n+1)) goto 70
    a=(xxx-xx(n))/(xx(n+1)-xx(n))
    m=0
    80 m=m+1
    if(qsq > qq(m+1)) goto 80
    b=(qsq-qq(m))/(qq(m+1)-qq(m))
    do 60 i=1,np
        g(i)= (1d0-a)*(1d0-b)*f(i,n,m)   + (1d0-a)*b*f(i,n,m+1) &
        +       a*(1d0-b)*f(i,n+1,m) +       a*b*f(i,n+1,m+1)
        if(n >= ntenth) goto 65
        if(i == 5 .OR. i == 7) goto 65
        fac=(1d0-b)*f(i,ntenth,m)+b*f(i,ntenth,m+1)
        g(i)=fac*10d0**(g(i)-fac)
        65 continue
        g(i)=g(i)*(1d0-x)**n0(i)
    60 END DO
    upv=g(1)
    dnv=g(2)
    usea=g(4)
    dsea=g(8)
    str=g(6)
    chm=g(5)
    glu=g(3)
    bot=g(7)
    x=xsave
    qsq=q2save
    return
    end subroutine mrs992

    subroutine mrs993(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
    use vamp_kinds
    implicit real(a-h,o-z)
    parameter(nx=49,nq=37,ntenth=23,np=8)
    real(kind=default) :: f(np,nx,nq+1),qq(nq),xx(nx),g(np),n0(np)
    data xx/1d-5,2d-5,4d-5,6d-5,8d-5, &
    &             1d-4,2d-4,4d-4,6d-4,8d-4, &
    &             1d-3,2d-3,4d-3,6d-3,8d-3, &
    &             1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2, &
    .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0, &
    .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0, &
    .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0, &
    .8d0,.9d0,1d0/
    data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1, &
    &         1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2, &
    &         1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4, &
    &         1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6, &
    &         1.8d6,3.2d6,5.6d6,1d7/
    data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
    data n0/3,4,5,9,9,9,9,9/
    data init/0/
    save
    xsave=x
    q2save=qsq
    if(init /= 0) goto 10
    open(unit=1,file='cor03.dat',status='old')
    do n=1,nx-1
        do m=1,nq
            read(1,50)f(1,n,m),f(2,n,m),f(3,n,m),f(4,n,m), &
            f(5,n,m),f(7,n,m),f(6,n,m),f(8,n,m)
        ! notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
            do 25 i=1,np
                f(i,n,m)=f(i,n,m)/(1d0-xx(n))**n0(i)
            25 END DO
        END DO
    end do
    do j=1,ntenth-1
        xx(j)=log10(xx(j)/xx(ntenth))+xx(ntenth)
        do i=1,8
            if(i == 5 .OR. i == 7) continue
            do 30 k=1,nq
                f(i,j,k)=log10(f(i,j,k)/f(i,ntenth,k))+f(i,ntenth,k)
            30 END DO
        END DO
    end do
    50 format(8f10.5)
    do i=1,np
        do m=1,nq
            f(i,nx,m)=0d0
        END DO
    end do
    init=1
    10 continue
    if(x < xmin) x=xmin
    if(x > xmax) x=xmax
    if(qsq < qsqmin)      qsq=qsqmin
    if(qsq > qsqmax)      qsq=qsqmax
    xxx=x
    if(x < xx(ntenth)) xxx=log10(x/xx(ntenth))+xx(ntenth)
    n=0
    70 n=n+1
    if(xxx > xx(n+1)) goto 70
    a=(xxx-xx(n))/(xx(n+1)-xx(n))
    m=0
    80 m=m+1
    if(qsq > qq(m+1)) goto 80
    b=(qsq-qq(m))/(qq(m+1)-qq(m))
    do 60 i=1,np
        g(i)= (1d0-a)*(1d0-b)*f(i,n,m)   + (1d0-a)*b*f(i,n,m+1) &
        +       a*(1d0-b)*f(i,n+1,m) +       a*b*f(i,n+1,m+1)
        if(n >= ntenth) goto 65
        if(i == 5 .OR. i == 7) goto 65
        fac=(1d0-b)*f(i,ntenth,m)+b*f(i,ntenth,m+1)
        g(i)=fac*10d0**(g(i)-fac)
        65 continue
        g(i)=g(i)*(1d0-x)**n0(i)
    60 END DO
    upv=g(1)
    dnv=g(2)
    usea=g(4)
    dsea=g(8)
    str=g(6)
    chm=g(5)
    glu=g(3)
    bot=g(7)
    x=xsave
    qsq=q2save
    return
    end subroutine mrs993


    subroutine mrs994(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
    use vamp_kinds
    implicit real(a-h,o-z)
    parameter(nx=49,nq=37,ntenth=23,np=8)
    real(kind=default) :: f(np,nx,nq+1),qq(nq),xx(nx),g(np),n0(np)
    data xx/1d-5,2d-5,4d-5,6d-5,8d-5, &
    &             1d-4,2d-4,4d-4,6d-4,8d-4, &
    &             1d-3,2d-3,4d-3,6d-3,8d-3, &
    &             1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2, &
    .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0, &
    .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0, &
    .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0, &
    .8d0,.9d0,1d0/
    data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1, &
    &         1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2, &
    &         1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4, &
    &         1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6, &
    &         1.8d6,3.2d6,5.6d6,1d7/
    data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
    data n0/3,4,5,9,9,9,9,9/
    data init/0/
    save
    xsave=x
    q2save=qsq
    if(init /= 0) goto 10
    open(unit=1,file='cor04.dat',status='old')
    do n=1,nx-1
        do m=1,nq
            read(1,50)f(1,n,m),f(2,n,m),f(3,n,m),f(4,n,m), &
            f(5,n,m),f(7,n,m),f(6,n,m),f(8,n,m)
        ! notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
            do 25 i=1,np
                f(i,n,m)=f(i,n,m)/(1d0-xx(n))**n0(i)
            25 END DO
        END DO
    end do
    do j=1,ntenth-1
        xx(j)=log10(xx(j)/xx(ntenth))+xx(ntenth)
        do i=1,8
            if(i == 5 .OR. i == 7) continue
            do 30 k=1,nq
                f(i,j,k)=log10(f(i,j,k)/f(i,ntenth,k))+f(i,ntenth,k)
            30 END DO
        END DO
    end do
    50 format(8f10.5)
    do i=1,np
        do m=1,nq
            f(i,nx,m)=0d0
        END DO
    end do
    init=1
    10 continue
    if(x < xmin) x=xmin
    if(x > xmax) x=xmax
    if(qsq < qsqmin)      qsq=qsqmin
    if(qsq > qsqmax)      qsq=qsqmax
    xxx=x
    if(x < xx(ntenth)) xxx=log10(x/xx(ntenth))+xx(ntenth)
    n=0
    70 n=n+1
    if(xxx > xx(n+1)) goto 70
    a=(xxx-xx(n))/(xx(n+1)-xx(n))
    m=0
    80 m=m+1
    if(qsq > qq(m+1)) goto 80
    b=(qsq-qq(m))/(qq(m+1)-qq(m))
    do 60 i=1,np
        g(i)= (1d0-a)*(1d0-b)*f(i,n,m)   + (1d0-a)*b*f(i,n,m+1) &
        +       a*(1d0-b)*f(i,n+1,m) +       a*b*f(i,n+1,m+1)
        if(n >= ntenth) goto 65
        if(i == 5 .OR. i == 7) goto 65
        fac=(1d0-b)*f(i,ntenth,m)+b*f(i,ntenth,m+1)
        g(i)=fac*10d0**(g(i)-fac)
        65 continue
        g(i)=g(i)*(1d0-x)**n0(i)
    60 END DO
    upv=g(1)
    dnv=g(2)
    usea=g(4)
    dsea=g(8)
    str=g(6)
    chm=g(5)
    glu=g(3)
    bot=g(7)
    x=xsave
    qsq=q2save
    return
    end subroutine mrs994

    subroutine mrs995(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
    use vamp_kinds
    implicit real(a-h,o-z)
    parameter(nx=49,nq=37,ntenth=23,np=8)
    real(kind=default) :: f(np,nx,nq+1),qq(nq),xx(nx),g(np),n0(np)
    data xx/1d-5,2d-5,4d-5,6d-5,8d-5, &
    &             1d-4,2d-4,4d-4,6d-4,8d-4, &
    &             1d-3,2d-3,4d-3,6d-3,8d-3, &
    &             1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2, &
    .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0, &
    .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0, &
    .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0, &
    .8d0,.9d0,1d0/
    data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1, &
    &         1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2, &
    &         1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4, &
    &         1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6, &
    &         1.8d6,3.2d6,5.6d6,1d7/
    data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
    data n0/3,4,5,9,9,9,9,9/
    data init/0/
    save
    xsave=x
    q2save=qsq
    if(init /= 0) goto 10
    open(unit=1,file='cor05.dat',status='old')
    do n=1,nx-1
        do m=1,nq
            read(1,50)f(1,n,m),f(2,n,m),f(3,n,m),f(4,n,m), &
            f(5,n,m),f(7,n,m),f(6,n,m),f(8,n,m)
        ! notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
            do 25 i=1,np
                f(i,n,m)=f(i,n,m)/(1d0-xx(n))**n0(i)
            25 END DO
        END DO
    end do
    do j=1,ntenth-1
        xx(j)=log10(xx(j)/xx(ntenth))+xx(ntenth)
        do i=1,8
            if(i == 5 .OR. i == 7) continue
            do 30 k=1,nq
                f(i,j,k)=log10(f(i,j,k)/f(i,ntenth,k))+f(i,ntenth,k)
            30 END DO
        END DO
    end do
    50 format(8f10.5)
    do i=1,np
        do m=1,nq
            f(i,nx,m)=0d0
        END DO
    end do
    init=1
    10 continue
    if(x < xmin) x=xmin
    if(x > xmax) x=xmax
    if(qsq < qsqmin)      qsq=qsqmin
    if(qsq > qsqmax)      qsq=qsqmax
    xxx=x
    if(x < xx(ntenth)) xxx=log10(x/xx(ntenth))+xx(ntenth)
    n=0
    70 n=n+1
    if(xxx > xx(n+1)) goto 70
    a=(xxx-xx(n))/(xx(n+1)-xx(n))
    m=0
    80 m=m+1
    if(qsq > qq(m+1)) goto 80
    b=(qsq-qq(m))/(qq(m+1)-qq(m))
    do 60 i=1,np
        g(i)= (1d0-a)*(1d0-b)*f(i,n,m)   + (1d0-a)*b*f(i,n,m+1) &
        +       a*(1d0-b)*f(i,n+1,m) +       a*b*f(i,n+1,m+1)
        if(n >= ntenth) goto 65
        if(i == 5 .OR. i == 7) goto 65
        fac=(1d0-b)*f(i,ntenth,m)+b*f(i,ntenth,m+1)
        g(i)=fac*10d0**(g(i)-fac)
        65 continue
        g(i)=g(i)*(1d0-x)**n0(i)
    60 END DO
    upv=g(1)
    dnv=g(2)
    usea=g(4)
    dsea=g(8)
    str=g(6)
    chm=g(5)
    glu=g(3)
    bot=g(7)
    x=xsave
    qsq=q2save
    return
    end subroutine mrs995


    subroutine mrs996(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
    use vamp_kinds
    implicit real(a-h,o-z)
    parameter(nx=49,nq=37,ntenth=23,np=8)
    real(kind=default) :: f(np,nx,nq+1),qq(nq),xx(nx),g(np),n0(np)
    data xx/1d-5,2d-5,4d-5,6d-5,8d-5, &
    &             1d-4,2d-4,4d-4,6d-4,8d-4, &
    &             1d-3,2d-3,4d-3,6d-3,8d-3, &
    &             1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2, &
    .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0, &
    .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0, &
    .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0, &
    .8d0,.9d0,1d0/
    data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1, &
    &         1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2, &
    &         1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4, &
    &         1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6, &
    &         1.8d6,3.2d6,5.6d6,1d7/
    data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
    data n0/3,4,5,9,9,9,9,9/
    data init/0/
    save
    xsave=x
    q2save=qsq
    if(init /= 0) goto 10
    open(unit=1,file='cor06.dat',status='old')
    do n=1,nx-1
        do m=1,nq
            read(1,50)f(1,n,m),f(2,n,m),f(3,n,m),f(4,n,m), &
            f(5,n,m),f(7,n,m),f(6,n,m),f(8,n,m)
        ! notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
            do 25 i=1,np
                f(i,n,m)=f(i,n,m)/(1d0-xx(n))**n0(i)
            25 END DO
        END DO
    end do
    do j=1,ntenth-1
        xx(j)=log10(xx(j)/xx(ntenth))+xx(ntenth)
        do i=1,8
            if(i == 5 .OR. i == 7) continue
            do 30 k=1,nq
                f(i,j,k)=log10(f(i,j,k)/f(i,ntenth,k))+f(i,ntenth,k)
            30 END DO
        END DO
    end do
    50 format(8f10.5)
    do i=1,np
        do m=1,nq
            f(i,nx,m)=0d0
        END DO
    end do
    init=1
    10 continue
    if(x < xmin) x=xmin
    if(x > xmax) x=xmax
    if(qsq < qsqmin)      qsq=qsqmin
    if(qsq > qsqmax)      qsq=qsqmax
    xxx=x
    if(x < xx(ntenth)) xxx=log10(x/xx(ntenth))+xx(ntenth)
    n=0
    70 n=n+1
    if(xxx > xx(n+1)) goto 70
    a=(xxx-xx(n))/(xx(n+1)-xx(n))
    m=0
    80 m=m+1
    if(qsq > qq(m+1)) goto 80
    b=(qsq-qq(m))/(qq(m+1)-qq(m))
    do 60 i=1,np
        g(i)= (1d0-a)*(1d0-b)*f(i,n,m)   + (1d0-a)*b*f(i,n,m+1) &
        +       a*(1d0-b)*f(i,n+1,m) +       a*b*f(i,n+1,m+1)
        if(n >= ntenth) goto 65
        if(i == 5 .OR. i == 7) goto 65
        fac=(1d0-b)*f(i,ntenth,m)+b*f(i,ntenth,m+1)
        g(i)=fac*10d0**(g(i)-fac)
        65 continue
        g(i)=g(i)*(1d0-x)**n0(i)
    60 END DO
    upv=g(1)
    dnv=g(2)
    usea=g(4)
    dsea=g(8)
    str=g(6)
    chm=g(5)
    glu=g(3)
    bot=g(7)
    x=xsave
    qsq=q2save
    return
    end subroutine mrs996

    subroutine mrs997(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
    use vamp_kinds
    implicit real(a-h,o-z)
    parameter(nx=49,nq=37,ntenth=23,np=8)
    real(kind=default) :: f(np,nx,nq+1),qq(nq),xx(nx),g(np),n0(np)
    data xx/1d-5,2d-5,4d-5,6d-5,8d-5, &
    &             1d-4,2d-4,4d-4,6d-4,8d-4, &
    &             1d-3,2d-3,4d-3,6d-3,8d-3, &
    &             1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2, &
    .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0, &
    .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0, &
    .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0, &
    .8d0,.9d0,1d0/
    data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1, &
    &         1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2, &
    &         1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4, &
    &         1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6, &
    &         1.8d6,3.2d6,5.6d6,1d7/
    data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
    data n0/3,4,5,9,9,9,9,9/
    data init/0/
    save
    xsave=x
    q2save=qsq
    if(init /= 0) goto 10
    open(unit=1,file='cor07.dat',status='old')
    do n=1,nx-1
        do m=1,nq
            read(1,50)f(1,n,m),f(2,n,m),f(3,n,m),f(4,n,m), &
            f(5,n,m),f(7,n,m),f(6,n,m),f(8,n,m)
        ! notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
            do 25 i=1,np
                f(i,n,m)=f(i,n,m)/(1d0-xx(n))**n0(i)
            25 END DO
        END DO
    end do
    do j=1,ntenth-1
        xx(j)=log10(xx(j)/xx(ntenth))+xx(ntenth)
        do i=1,8
            if(i == 5 .OR. i == 7) continue
            do 30 k=1,nq
                f(i,j,k)=log10(f(i,j,k)/f(i,ntenth,k))+f(i,ntenth,k)
            30 END DO
        END DO
    end do
    50 format(8f10.5)
    do i=1,np
        do m=1,nq
            f(i,nx,m)=0d0
        END DO
    end do
    init=1
    10 continue
    if(x < xmin) x=xmin
    if(x > xmax) x=xmax
    if(qsq < qsqmin)      qsq=qsqmin
    if(qsq > qsqmax)      qsq=qsqmax
    xxx=x
    if(x < xx(ntenth)) xxx=log10(x/xx(ntenth))+xx(ntenth)
    n=0
    70 n=n+1
    if(xxx > xx(n+1)) goto 70
    a=(xxx-xx(n))/(xx(n+1)-xx(n))
    m=0
    80 m=m+1
    if(qsq > qq(m+1)) goto 80
    b=(qsq-qq(m))/(qq(m+1)-qq(m))
    do 60 i=1,np
        g(i)= (1d0-a)*(1d0-b)*f(i,n,m)   + (1d0-a)*b*f(i,n,m+1) &
        +       a*(1d0-b)*f(i,n+1,m) +       a*b*f(i,n+1,m+1)
        if(n >= ntenth) goto 65
        if(i == 5 .OR. i == 7) goto 65
        fac=(1d0-b)*f(i,ntenth,m)+b*f(i,ntenth,m+1)
        g(i)=fac*10d0**(g(i)-fac)
        65 continue
        g(i)=g(i)*(1d0-x)**n0(i)
    60 END DO
    upv=g(1)
    dnv=g(2)
    usea=g(4)
    dsea=g(8)
    str=g(6)
    chm=g(5)
    glu=g(3)
    bot=g(7)
    x=xsave
    qsq=q2save
    return
    end subroutine mrs997



    subroutine mrs998(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
    use vamp_kinds
    implicit real(a-h,o-z)
    parameter(nx=49,nq=37,ntenth=23,np=8)
    real(kind=default) :: f(np,nx,nq+1),qq(nq),xx(nx),g(np),n0(np)
    data xx/1d-5,2d-5,4d-5,6d-5,8d-5, &
    &             1d-4,2d-4,4d-4,6d-4,8d-4, &
    &             1d-3,2d-3,4d-3,6d-3,8d-3, &
    &             1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2, &
    .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0, &
    .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0, &
    .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0, &
    .8d0,.9d0,1d0/
    data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1, &
    &         1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2, &
    &         1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4, &
    &         1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6, &
    &         1.8d6,3.2d6,5.6d6,1d7/
    data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
    data n0/3,4,5,9,9,9,9,9/
    data init/0/
    save
    xsave=x
    q2save=qsq
    if(init /= 0) goto 10
    open(unit=1,file='cor08.dat',status='old')
    do n=1,nx-1
        do m=1,nq
            read(1,50)f(1,n,m),f(2,n,m),f(3,n,m),f(4,n,m), &
            f(5,n,m),f(7,n,m),f(6,n,m),f(8,n,m)
        ! notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
            do 25 i=1,np
                f(i,n,m)=f(i,n,m)/(1d0-xx(n))**n0(i)
            25 END DO
        END DO
    end do
    do j=1,ntenth-1
        xx(j)=log10(xx(j)/xx(ntenth))+xx(ntenth)
        do i=1,8
            if(i == 5 .OR. i == 7) continue
            do 30 k=1,nq
                f(i,j,k)=log10(f(i,j,k)/f(i,ntenth,k))+f(i,ntenth,k)
            30 END DO
        END DO
    end do
    50 format(8f10.5)
    do i=1,np
        do m=1,nq
            f(i,nx,m)=0d0
        END DO
    end do
    init=1
    10 continue
    if(x < xmin) x=xmin
    if(x > xmax) x=xmax
    if(qsq < qsqmin)      qsq=qsqmin
    if(qsq > qsqmax)      qsq=qsqmax
    xxx=x
    if(x < xx(ntenth)) xxx=log10(x/xx(ntenth))+xx(ntenth)
    n=0
    70 n=n+1
    if(xxx > xx(n+1)) goto 70
    a=(xxx-xx(n))/(xx(n+1)-xx(n))
    m=0
    80 m=m+1
    if(qsq > qq(m+1)) goto 80
    b=(qsq-qq(m))/(qq(m+1)-qq(m))
    do 60 i=1,np
        g(i)= (1d0-a)*(1d0-b)*f(i,n,m)   + (1d0-a)*b*f(i,n,m+1) &
        +       a*(1d0-b)*f(i,n+1,m) +       a*b*f(i,n+1,m+1)
        if(n >= ntenth) goto 65
        if(i == 5 .OR. i == 7) goto 65
        fac=(1d0-b)*f(i,ntenth,m)+b*f(i,ntenth,m+1)
        g(i)=fac*10d0**(g(i)-fac)
        65 continue
        g(i)=g(i)*(1d0-x)**n0(i)
    60 END DO
    upv=g(1)
    dnv=g(2)
    usea=g(4)
    dsea=g(8)
    str=g(6)
    chm=g(5)
    glu=g(3)
    bot=g(7)
    x=xsave
    qsq=q2save
    return
    end subroutine mrs998

    subroutine mrs999(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
    use vamp_kinds
    implicit real(a-h,o-z)
    parameter(nx=49,nq=37,ntenth=23,np=8)
    real(kind=default) :: f(np,nx,nq+1),qq(nq),xx(nx),g(np),n0(np)
    data xx/1d-5,2d-5,4d-5,6d-5,8d-5, &
    &             1d-4,2d-4,4d-4,6d-4,8d-4, &
    &             1d-3,2d-3,4d-3,6d-3,8d-3, &
    &             1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2, &
    .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0, &
    .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0, &
    .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0, &
    .8d0,.9d0,1d0/
    data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1, &
    &         1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2, &
    &         1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4, &
    &         1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6, &
    &         1.8d6,3.2d6,5.6d6,1d7/
    data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
    data n0/3,4,5,9,9,9,9,9/
    data init/0/
    save
    xsave=x
    q2save=qsq
    if(init /= 0) goto 10
    open(unit=1,file='cor09.dat',status='old')
    do n=1,nx-1
        do m=1,nq
            read(1,50)f(1,n,m),f(2,n,m),f(3,n,m),f(4,n,m), &
            f(5,n,m),f(7,n,m),f(6,n,m),f(8,n,m)
        ! notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
            do 25 i=1,np
                f(i,n,m)=f(i,n,m)/(1d0-xx(n))**n0(i)
            25 END DO
        END DO
    end do
    do j=1,ntenth-1
        xx(j)=log10(xx(j)/xx(ntenth))+xx(ntenth)
        do i=1,8
            if(i == 5 .OR. i == 7) continue
            do 30 k=1,nq
                f(i,j,k)=log10(f(i,j,k)/f(i,ntenth,k))+f(i,ntenth,k)
            30 END DO
        END DO
    end do
    50 format(8f10.5)
    do i=1,np
        do m=1,nq
            f(i,nx,m)=0d0
        END DO
    end do
    init=1
    10 continue
    if(x < xmin) x=xmin
    if(x > xmax) x=xmax
    if(qsq < qsqmin)      qsq=qsqmin
    if(qsq > qsqmax)      qsq=qsqmax
    xxx=x
    if(x < xx(ntenth)) xxx=log10(x/xx(ntenth))+xx(ntenth)
    n=0
    70 n=n+1
    if(xxx > xx(n+1)) goto 70
    a=(xxx-xx(n))/(xx(n+1)-xx(n))
    m=0
    80 m=m+1
    if(qsq > qq(m+1)) goto 80
    b=(qsq-qq(m))/(qq(m+1)-qq(m))
    do 60 i=1,np
        g(i)= (1d0-a)*(1d0-b)*f(i,n,m)   + (1d0-a)*b*f(i,n,m+1) &
        +       a*(1d0-b)*f(i,n+1,m) +       a*b*f(i,n+1,m+1)
        if(n >= ntenth) goto 65
        if(i == 5 .OR. i == 7) goto 65
        fac=(1d0-b)*f(i,ntenth,m)+b*f(i,ntenth,m+1)
        g(i)=fac*10d0**(g(i)-fac)
        65 continue
        g(i)=g(i)*(1d0-x)**n0(i)
    60 END DO
    upv=g(1)
    dnv=g(2)
    usea=g(4)
    dsea=g(8)
    str=g(6)
    chm=g(5)
    glu=g(3)
    bot=g(7)
    x=xsave
    qsq=q2save
    return
    end subroutine mrs999



    subroutine mrs9910(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
    use vamp_kinds
    implicit real(a-h,o-z)
    parameter(nx=49,nq=37,ntenth=23,np=8)
    real(kind=default) :: f(np,nx,nq+1),qq(nq),xx(nx),g(np),n0(np)
    data xx/1d-5,2d-5,4d-5,6d-5,8d-5, &
    &             1d-4,2d-4,4d-4,6d-4,8d-4, &
    &             1d-3,2d-3,4d-3,6d-3,8d-3, &
    &             1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2, &
    .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0, &
    .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0, &
    .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0, &
    .8d0,.9d0,1d0/
    data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1, &
    &         1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2, &
    &         1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4, &
    &         1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6, &
    &         1.8d6,3.2d6,5.6d6,1d7/
    data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
    data n0/3,4,5,9,9,9,9,9/
    data init/0/
    save
    xsave=x
    q2save=qsq
    if(init /= 0) goto 10
    open(unit=1,file='cor10.dat',status='old')
    do n=1,nx-1
        do m=1,nq
            read(1,50)f(1,n,m),f(2,n,m),f(3,n,m),f(4,n,m), &
            f(5,n,m),f(7,n,m),f(6,n,m),f(8,n,m)
        ! notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
            do 25 i=1,np
                f(i,n,m)=f(i,n,m)/(1d0-xx(n))**n0(i)
            25 END DO
        END DO
    end do
    do j=1,ntenth-1
        xx(j)=log10(xx(j)/xx(ntenth))+xx(ntenth)
        do i=1,8
            if(i == 5 .OR. i == 7) continue
            do 30 k=1,nq
                f(i,j,k)=log10(f(i,j,k)/f(i,ntenth,k))+f(i,ntenth,k)
            30 END DO
        END DO
    end do
    50 format(8f10.5)
    do i=1,np
        do m=1,nq
            f(i,nx,m)=0d0
        END DO
    end do
    init=1
    10 continue
    if(x < xmin) x=xmin
    if(x > xmax) x=xmax
    if(qsq < qsqmin)      qsq=qsqmin
    if(qsq > qsqmax)      qsq=qsqmax
    xxx=x
    if(x < xx(ntenth)) xxx=log10(x/xx(ntenth))+xx(ntenth)
    n=0
    70 n=n+1
    if(xxx > xx(n+1)) goto 70
    a=(xxx-xx(n))/(xx(n+1)-xx(n))
    m=0
    80 m=m+1
    if(qsq > qq(m+1)) goto 80
    b=(qsq-qq(m))/(qq(m+1)-qq(m))
    do 60 i=1,np
        g(i)= (1d0-a)*(1d0-b)*f(i,n,m)   + (1d0-a)*b*f(i,n,m+1) &
        +       a*(1d0-b)*f(i,n+1,m) +       a*b*f(i,n+1,m+1)
        if(n >= ntenth) goto 65
        if(i == 5 .OR. i == 7) goto 65
        fac=(1d0-b)*f(i,ntenth,m)+b*f(i,ntenth,m+1)
        g(i)=fac*10d0**(g(i)-fac)
        65 continue
        g(i)=g(i)*(1d0-x)**n0(i)
    60 END DO
    upv=g(1)
    dnv=g(2)
    usea=g(4)
    dsea=g(8)
    str=g(6)
    chm=g(5)
    glu=g(3)
    bot=g(7)
    x=xsave
    qsq=q2save
    return
    end subroutine mrs9910

    subroutine mrs9911(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
    use vamp_kinds
    implicit real(a-h,o-z)
    parameter(nx=49,nq=37,ntenth=23,np=8)
    real(kind=default) :: f(np,nx,nq+1),qq(nq),xx(nx),g(np),n0(np)
    data xx/1d-5,2d-5,4d-5,6d-5,8d-5, &
    &             1d-4,2d-4,4d-4,6d-4,8d-4, &
    &             1d-3,2d-3,4d-3,6d-3,8d-3, &
    &             1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2, &
    .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0, &
    .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0, &
    .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0, &
    .8d0,.9d0,1d0/
    data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1, &
    &         1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2, &
    &         1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4, &
    &         1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6, &
    &         1.8d6,3.2d6,5.6d6,1d7/
    data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
    data n0/3,4,5,9,9,9,9,9/
    data init/0/
    save
    xsave=x
    q2save=qsq
    if(init /= 0) goto 10
    open(unit=1,file='cor11.dat',status='old')
    do n=1,nx-1
        do m=1,nq
            read(1,50)f(1,n,m),f(2,n,m),f(3,n,m),f(4,n,m), &
            f(5,n,m),f(7,n,m),f(6,n,m),f(8,n,m)
        ! notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
            do 25 i=1,np
                f(i,n,m)=f(i,n,m)/(1d0-xx(n))**n0(i)
            25 END DO
        END DO
    end do
    do j=1,ntenth-1
        xx(j)=log10(xx(j)/xx(ntenth))+xx(ntenth)
        do i=1,8
            if(i == 5 .OR. i == 7) continue
            do 30 k=1,nq
                f(i,j,k)=log10(f(i,j,k)/f(i,ntenth,k))+f(i,ntenth,k)
            30 END DO
        END DO
    end do
    50 format(8f10.5)
    do i=1,np
        do m=1,nq
            f(i,nx,m)=0d0
        END DO
    end do
    init=1
    10 continue
    if(x < xmin) x=xmin
    if(x > xmax) x=xmax
    if(qsq < qsqmin)      qsq=qsqmin
    if(qsq > qsqmax)      qsq=qsqmax
    xxx=x
    if(x < xx(ntenth)) xxx=log10(x/xx(ntenth))+xx(ntenth)
    n=0
    70 n=n+1
    if(xxx > xx(n+1)) goto 70
    a=(xxx-xx(n))/(xx(n+1)-xx(n))
    m=0
    80 m=m+1
    if(qsq > qq(m+1)) goto 80
    b=(qsq-qq(m))/(qq(m+1)-qq(m))
    do 60 i=1,np
        g(i)= (1d0-a)*(1d0-b)*f(i,n,m)   + (1d0-a)*b*f(i,n,m+1) &
        +       a*(1d0-b)*f(i,n+1,m) +       a*b*f(i,n+1,m+1)
        if(n >= ntenth) goto 65
        if(i == 5 .OR. i == 7) goto 65
        fac=(1d0-b)*f(i,ntenth,m)+b*f(i,ntenth,m+1)
        g(i)=fac*10d0**(g(i)-fac)
        65 continue
        g(i)=g(i)*(1d0-x)**n0(i)
    60 END DO
    upv=g(1)
    dnv=g(2)
    usea=g(4)
    dsea=g(8)
    str=g(6)
    chm=g(5)
    glu=g(3)
    bot=g(7)
    x=xsave
    qsq=q2save
    return
    end subroutine mrs9911


    subroutine mrs9912(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu)
    use vamp_kinds
    implicit real(a-h,o-z)
    parameter(nx=49,nq=37,ntenth=23,np=8)
    real(kind=default) :: f(np,nx,nq+1),qq(nq),xx(nx),g(np),n0(np)
    data xx/1d-5,2d-5,4d-5,6d-5,8d-5, &
    &             1d-4,2d-4,4d-4,6d-4,8d-4, &
    &             1d-3,2d-3,4d-3,6d-3,8d-3, &
    &             1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2, &
    .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0, &
    .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0, &
    .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0, &
    .8d0,.9d0,1d0/
    data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1, &
    &         1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2, &
    &         1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4, &
    &         1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6, &
    &         1.8d6,3.2d6,5.6d6,1d7/
    data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
    data n0/3,4,5,9,9,9,9,9/
    data init/0/
    save
    xsave=x
    q2save=qsq
    if(init /= 0) goto 10
    open(unit=1,file='cor12.dat',status='old')
    do n=1,nx-1
        do m=1,nq
            read(1,50)f(1,n,m),f(2,n,m),f(3,n,m),f(4,n,m), &
            f(5,n,m),f(7,n,m),f(6,n,m),f(8,n,m)
        ! notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
            do 25 i=1,np
                f(i,n,m)=f(i,n,m)/(1d0-xx(n))**n0(i)
            25 END DO
        END DO
    end do
    do j=1,ntenth-1
        xx(j)=log10(xx(j)/xx(ntenth))+xx(ntenth)
        do i=1,8
            if(i == 5 .OR. i == 7) continue
            do 30 k=1,nq
                f(i,j,k)=log10(f(i,j,k)/f(i,ntenth,k))+f(i,ntenth,k)
            30 END DO
        END DO
    end do
    50 format(8f10.5)
    do i=1,np
        do m=1,nq
            f(i,nx,m)=0d0
        END DO
    end do
    init=1
    10 continue
    if(x < xmin) x=xmin
    if(x > xmax) x=xmax
    if(qsq < qsqmin)      qsq=qsqmin
    if(qsq > qsqmax)      qsq=qsqmax
    xxx=x
    if(x < xx(ntenth)) xxx=log10(x/xx(ntenth))+xx(ntenth)
    n=0
    70 n=n+1
    if(xxx > xx(n+1)) goto 70
    a=(xxx-xx(n))/(xx(n+1)-xx(n))
    m=0
    80 m=m+1
    if(qsq > qq(m+1)) goto 80
    b=(qsq-qq(m))/(qq(m+1)-qq(m))
    do 60 i=1,np
        g(i)= (1d0-a)*(1d0-b)*f(i,n,m)   + (1d0-a)*b*f(i,n,m+1) &
        +       a*(1d0-b)*f(i,n+1,m) +       a*b*f(i,n+1,m+1)
        if(n >= ntenth) goto 65
        if(i == 5 .OR. i == 7) goto 65
        fac=(1d0-b)*f(i,ntenth,m)+b*f(i,ntenth,m+1)
        g(i)=fac*10d0**(g(i)-fac)
        65 continue
        g(i)=g(i)*(1d0-x)**n0(i)
    60 END DO
    upv=g(1)
    dnv=g(2)
    usea=g(4)
    dsea=g(8)
    str=g(6)
    chm=g(5)
    glu=g(3)
    bot=g(7)
    x=xsave
    qsq=q2save
    return
    end subroutine mrs9912
