c ====================================================================== 
      subroutine initialise_madgraph(o_NWA,model)
c sets up masses and coupling constants of particles
c ----------------------------------------------------------------------
      implicit none

c arguments
      integer o_NWA
      character*50 model

c SM couplings
      real*8        gw, gwwa, gwwZ
      common/coup1/ gw, gwwa, gwwZ
      real*8        gal(2),gau(2),gad(2),gwf(2)
      common/coup2a/gal,   gau,   gad,   gwf
      real*8        gZn(2),gZl(2),gZu(2),gZd(2),g1(2)
      common/coup2b/gZn,   gZl,   gZu,   gZd,   g1
      real*8        gwwh,gZZh,ghhh,gwwhh,gZZhh,ghhhh
      common/coup3/ gwwh,gZZh,ghhh,gwwhh,gZZhh,ghhhh
      complex*16    gchf(2,12)
      common/coup4/ gchf
      real*8          gg(2), g
      common/coupqcd/ gg,    g

c SM masses and widths
      real*8           fmass(12), fwidth(12)
      common/fermions/ fmass,     fwidth
      real*8        rm_W,Gamma_W,rm_Z,Gamma_Z
      common/vmass1/rm_W,Gamma_W,rm_Z,Gamma_Z
      real*8        rm_A,Gamma_a,rm_h,Gamma_h
      common/vmass2/rm_A,Gamma_a,rm_h,Gamma_h
      real*8        Gamma_t
      common/Gammas/Gamma_t
c Other SM parameters
      real*8    a_EM,s2w
      common/EW/a_EM,s2w
      real*8     rlambdaQCD4
      integer                nloops
      common/QCD/rlambdaQCD4,nloops


c Zprime parameters
      real*8    rmZp(5),gamZp(5)
      common/Zp/rmZp   ,gamZp
      real*8         paramZp(5)
      common/Zpparam/paramZp
      real*8          gp(5),gV_d(5),gA_d(5),gV_u(5),gA_u(5)
      common/coupZpVA/gp   ,gV_d   ,gA_d   ,gV_u   ,gA_u
      real*8        gZpd(2,5),gZpu(2,5)
      common/coupZp/gZpd     ,gZpu
c   constants
      real*8     alpha_EM             ,sin2theta_weinberg
      parameter (alpha_EM=0.0078125   ,sin2theta_weinberg=0.2320d0)
      ! integer igw

      real*8   widthZp
      external widthZp

c   quark masses
      real*8     umass,       cmass,            tmass
      parameter (umass=0.d0,  cmass=0.d0,       tmass=175.d0)
      real*8     dmass,       smass,            bmass
      parameter (dmass=0.d0,  smass=0.d0,       bmass=4.18d0)

c   leptons masses
      real*8     emass,        mumass,          taumass
      parameter (emass=0.d0,   mumass=0.d0,     taumass=1.78d0)
      real*8     nuemass,      numumass,        nutaumass
      parameter (nuemass=0d0,  numumass=0d0,    nutaumass=0d0)

c   quark widths
      real*8     uwidth,       cwidth,          twidth
      parameter (uwidth=0.d0,  cwidth=0.d0,     twidth=1.55d0)
      real*8     dwidth,       swidth,          bwidth
      parameter (dwidth=0.d0,  swidth=0.d0,     bwidth=0.d0)

c   lepton widths
      real*8     eWidth,        muwidth,        tauwidth
      parameter (eWidth=0.d0,   muwidth=0.d0,   tauwidth=0.d0)
      real*8     nueWidth,      numuwidth,      nutauwidth
      parameter (nueWidth=0.d0, numuwidth=0.d0, nutauwidth=0.d0)

c   SM boson masses
      real*8     Wmass,      Zmass,      Wwidth,     Zwidth
      parameter (Wmass=80.23d0,Zmass=91.19d0,Wwidth=2.08d0,Zwidth=2.5d0)
      real*8     Amass,     Awidth,     hmass,        hwidth
      parameter (Amass=0d0, Awidth=0d0, hmass=125.d0, hwidth=0.31278d-2)

c other local variables
      integer i,imodel,o_width(5)
c ----------------------------------------------------------------------
c enter global fermion masses  
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
      if(o_NWA.eq.1)then 
        fwidth(11)=1.d-5
      else     
        fwidth(11)=twidth
      end if
      fwidth(12)= bwidth

!    Save EW parameters to common blocks
      a_em=alpha_EM
      s2w=sin2theta_weinberg      

!     enter global boson masses and widths
      rm_W=Zmass*dsqrt(1.d0-s2w)
      rm_Z=Zmass
      rm_A=Amass
      rm_h=hmass
      Gamma_t=twidth
      Gamma_W=wwidth
      Gamma_Z=Zwidth
      Gamma_A=awidth
      Gamma_h=hwidth

!   call SM HELAS couplings
      call coup1x(s2w,gw,gwwa,gwwZ)
      call coup2x(s2w,gal,gau,gad,gwf,gZn,gZl,gZu,gZd,g1)
      call coup3x(s2w,Zmass,hmass,gwwh,gZZh,ghhh,gwwhh,gZZhh,ghhhh)
      do i=1,12
         call coup4x(s2w,Zmass,fmass(i),gchf(1,i))
      enddo

!   QCD couplings
      g = 1.d0
      gg(1)=-g
      gg(2)=-g

      ! Extract model filename (Remove white space.)
      imodel = len(model)
      do while(model(imodel:imodel).eq.'') 
        imodel = imodel-1
      end do

      ! Read model file
      open(unit=42,file='Models/'//model(1:imodel)//'.mdl',status='old')
      read(42,*) rmZp
      read(42,*) gamZp      
      read(42,*) gp 
      read(42,*) paramZp
      read(42,*) gV_u
      read(42,*) gA_u
      read(42,*) gV_d
      read(42,*) gA_d

  !   Check whether width has been specified
  !   (If gamZp is zero, the function widthZp is used instead.)
      do i=1,5
        if ((gamZp(i).eq.0d0).and.(rmZp(i).gt.0d0)) then
          o_width(i) = 0
        else
          o_width(i) = 1
        end if
      enddo

  ! Calculate sequential Zp widths
      do i=1,5
        if (o_width(i).eq.0) gamZp(i)=
     &            widthZp(rm_W,rm_Z,rmZp(i),a_em,s2w,rlambdaQCD4,nloops)
      end do

      ! convert from VA to LR couplings
      call coupZpx

      ! igw=0 ! don't include w width effects
      ! call topwid(fmass(11),wmass,fmass(12),wwidth,igw,fwidth(11))
      ! call printconstants
      return
      end
c ====================================================================== 

!       subroutine printconstants
! c*************************************************************************
! c     prints out all masses, widths, and couplings in common blocks
! c*************************************************************************
!       implicit none
! c
! c     local
! c
! c      integer i
! c
! c     global
! c
!       real*8          gw, gwwa, gwwZ
!       common /coup1/ gw, gwwa, gwwZ
!       real*8         gal(2),gau(2),gad(2),gwf(2)
!       common /coup2a/gal,   gau,   gad,   gwf
!       real*8         gZn(2),gZl(2),gZu(2),gZd(2),g1(2)
!       common /coup2b/gZn,   gZl,   gZu,   gZd,   g1
!       real*8         gwwh,gZZh,ghhh,gwwhh,gZZhh,ghhhh
!       common /coup3/ gwwh,gZZh,ghhh,gwwhh,gZZhh,ghhhh
!       complex*16     gchf(2,12)
!       common /coup4/ gchf
!       real*8         wmass,wwidth,Zmass,Zwidth
!       common /vmass1/wmass,wwidth,Zmass,Zwidth
!       real*8         Amass,awidth,hmass,hwidth
!       common /vmass2/Amass,awidth,hmass,hwidth
!       real*8            fmass(12), fwidth(12)
!       common /fermions/ fmass,     fwidth
!       real*8           gg(2), g
!       common /coupqcd/ gg,    g

! c --------------------------------- method -------------------------------------
!       write(*,'(a)') 'MADGRAPH Parameters'
!       write(*,'(a)') 'Boson masses and Widths:'
!       write(*,10) 'W',wmass,'W',wwidth,'Z',Zmass,'Z',Zwidth
!       write(*,10) 'A',Amass,'A',awidth,'H',hmass,'H',hwidth
!       write(*,'(a)') 'Quark masses and Widths:'
!       write(*,10) 't',fmass(11),'t',fwidth(11),
!      &            'b',fmass(12),'b',fwidth(12)


!  10   format(a,1x,5Hmass=,f6.2,1H:,2x,a,1x,6Hwidth=,f6.2,1H:,2x,
!      &       a,1x,5Hmass=,f6.2,1H:,2x,a,1x,6Hwidth=,f6.2)

!       end
! c ================================ end =========================================

!       subroutine topwid(rmt,rmw,rmb,rgw,igw,rgt)
! c*************************************************************************
! c     the total weak decay width of the top quark, including
! c     the effects of bottom mass and, if igw=1,  a finite w width.
! c     from james stirling 6-10-94
! c*************************************************************************
!       implicit complex*16(a-h,o-Z)
!       real*8 rmt,rmb,rmw,xw,xb,rgw,rgt
! *
!       pi=4.*datan(1.d0)
!       gf=1.16637d-05
!       gw=cdsqrt(rmw**2*gf/dsqrt(2.d0))
! *                            flavour & colour
!       xb=rmb/rmt
!       xw=rmw/rmt
!       if(igw.eq.1) goto 10
!       if(rmt.le.(rmw+rmb)) then
!           write(6,*)'WARNING: mt < mw + mb !!!!'
!           stop
!           endif
!       rgt = gf*rmt**3/8d0/pi/dsqrt(2d0) 
!      .   * dsqrt( (1d0-(xw+xb)**2)*(1d0-(xw-xb)**2) )
!      .   * ( (1d0-xb**2)**2 + (1d0+xb**2)*xw**2 - 2d0*xw**4 )
!       return
!   10  continue
!       rm=xb**2
!       om=1.+rm-dcmplx(rmw**2,rmw*rgw)/rmt**2
!       y1=om+cdsqrt(om*om-4.*rm)
!       y0=om-cdsqrt(om*om-4.*rm)
!       Z1=2.
!       Z0=2.*cdsqrt(rm)
! *
!       d0=(-y0**8+3.*y0**7*rm+3.*y0**7-8.*y0**6*rm-12.*y0**5*rm**
!      . 2-12.*y0**5*rm+96.*y0**4*rm**2-48.*y0**3*rm**3-48.*y0**3*
!      . rm**2-128.*y0**2*rm**3+192.*y0*rm**4+192.*y0*rm**3-256.*
!      . rm**4)/(24.*y0**4*(y1-y0))
!       d1=(-y1**8+3.*y1**7*rm+3.*y1**7-8.*y1**6*rm-12.*y1**5*rm**
!      . 2-12.*y1**5*rm+96.*y1**4*rm**2-48.*y1**3*rm**3-48.*y1**3*
!      . rm**2-128.*y1**2*rm**3+192.*y1*rm**4+192.*y1*rm**3-256.*
!      . rm**4)/(24.*y1**4*(y1-y0))
!       a4=(32.*rm**4*(y1-y0))/(3.*y1*y0*(y1-y0))
!       a3=(8.*rm**3*(-3.*y1**2*y0*rm-3.*y1**2*y0+4.*y1**2*rm+3.*
!      . y1*y0**2*rm+3.*y1*y0**2-4.*y0**2*rm))/(3.*y1**2*y0**2*(y1
!      . -y0))
!       a2=(8.*rm**3*(2.*y1**3*y0**2-3.*y1**3*y0*rm-3.*y1**3*y0+4.
!      . *y1**3*rm-2.*y1**2*y0**3+3.*y1*y0**3*rm+3.*y1*y0**3-4.*y0
!      . **3*rm))/(3.*y1**3*y0**3*(y1-y0))
!       a1=(2.*rm**2*(3.*y1**4*y0**3*rm+3.*y1**4*y0**3+8.*y1**4*y0
!      . **2*rm-12.*y1**4*y0*rm**2-12.*y1**4*y0*rm+16.*y1**4*rm**2
!      . -3.*y1**3*y0**4*rm-3.*y1**3*y0**4-8.*y1**2*y0**4*rm+12.*
!      . y1*y0**4*rm**2+12.*y1*y0**4*rm-16.*y0**4*rm**2))/(3.*y1**
!      . 4*y0**4*(y1-y0))
!       b0=(y1**3-3.*y1**2*rm-3.*y1**2+8.*y1*rm-y0**3+3.*y0**2*rm+
!      . 3.*y0**2-8.*y0*rm)/(24.*(y1-y0))
!       b1=(y1+y0-3.*rm-3.)/24.
!       b2=1./24.
! *
!       rint=d0*cdlog((Z1-y0)/(Z0-y0))
!      .    -d1*cdlog((y1-Z1)/(y1-Z0))
!      .    -a4/3.*(1./Z1**3-1./Z0**3)
!      .    -a3/2.*(1./Z1**2-1./Z0**2)
!      .    -a2   *(1./Z1   -1./Z0   )
!      .    +a1*cdlog(Z1/Z0)
!      .    +b0   *(Z1   -Z0   )
!      .    +b1/2.*(Z1**2-Z0**2)
!      .    +b2/3.*(Z1**3-Z0**3)
! *
!       gw4=gw**4
! *
! * total width includes flavour & colour factors
!       rgt=rmt**3/(rmw*rgw)*gw4/(8.*pi**3)*dimag(rint)
!       rgt=9.*rgt
!       return
!       end
