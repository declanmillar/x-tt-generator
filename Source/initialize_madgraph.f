c ============================== subroutine ====================================      
      subroutine initialise_madgraph(rmt,gamt)
c     sets up masses and coupling constants of particles
c ----------------------------- declarations -----------------------------------
      implicit none
c arguments
      real*8 rmt,gamt

c constants
      real*8     sw2
      parameter (sw2 = .2320d0) ! sin(theta_weinberg)
      integer igw

c define masses and widths of fermions

      real*8     tmass,      bmass,    cmass,    smass,    umass
      parameter (tmass=175.d0,bmass=4.18d0,cmass=0.d0,
     &                        smass=0.d0,umass=0.d0)
      real*8     twidth,    bwidth,    cwidth,    swidth,    uwidth
      parameter (twidth=0.d0,bwidth=0.d0,cwidth=0.d0,
     &                         swidth=0.d0,uwidth=0.d0)
      real*8     dmass,     emass,    mumass,    taumass
      parameter (dmass=0.d0,emass=0.d0,mumass=0.d0,
     &                                 taumass=1.78d0)
      real*8     dwidth,    ewidth,    muwidth,    tauwidth
      parameter (dwidth=0d0,ewidth=0d0,muwidth=0d0,tauwidth=0d0)

c define masses and widths of the SM bosons
      real*8     wmass,      zmass,      wwidth,     zwidth
      parameter (wmass=80.23d0, zmass=91.19d0, 
     &           wwidth=2.08d0, zwidth=2.50d0)
      real*8     amass,     awidth,     hmass,        hwidth
      parameter (amass=0d0, awidth=0d0, hmass=125.d0, hwidth=0.31278d-2)

c local variables
      integer i

c global variables
      real*8         gw, gwwa, gwwz
      common /coup1/ gw, gwwa, gwwz
      real*8         gal(2),gau(2),gad(2),gwf(2)
      common /coup2a/gal,   gau,   gad,   gwf
      real*8         gzn(2),gzl(2),gzu(2),gzd(2),g1(2)
      common /coup2b/gzn,   gzl,   gzu,   gzd,   g1
      real*8        gwwh,gzzh,ghhh,gwwhh,gzzhh,ghhhh
      common /coup3/ gwwh,gzzh,ghhh,gwwhh,gzzhh,ghhhh
      complex*16     gchf(2,12)
      common /coup4/ gchf
      real*8         Wmass1,Wwidth1,zmass1,zwidth1
      common /vmass1/wmass1,wwidth1,zmass1,zwidth1
      real*8         amass1,awidth1,hmass1,hwidth1
      common /vmass2/amass1,awidth1,hmass1,hwidth1
      real*8            fmass(12), fwidth(12)
      common /fermions/ fmass,     fwidth
      real*8           gg(2), g
      common /coupqcd/ gg,    g

c ---------------------------------- method ----------------------------
c enter fermion masses  
      fmass(1) = emass
      fmass(2) = 0d0
      fmass(3) = umass
      fmass(4) = dmass
      fmass(5) = mumass
      fmass(6) = 0d0
      fmass(7) = cmass
      fmass(8) = smass
      fmass(9) = taumass
      fmass(10)= 0d0
      fmass(11)= rmt
      fmass(12)= bmass

      fwidth(1) = ewidth
      fwidth(2) = 0d0
      fwidth(3) = uwidth
      fwidth(4) = dwidth
      fwidth(5) = muwidth
      fwidth(6) = 0d0
      fwidth(7) = cwidth
      fwidth(8) = swidth
      fwidth(9) = tauwidth
      fwidth(10)= 0d0
      fwidth(11)= gamt
      fwidth(12)= bwidth

!      wmass1=wmass??!!!!!
      wmass1=zmass*dsqrt(1.d0-sw2)
      zmass1=zmass
      amass1=amass
      hmass1=hmass
      wwidth1=wwidth
      zwidth1=zwidth
      awidth1=awidth
      hwidth1=hwidth
      call coup1x(sw2,gw,gwwa,gwwz)
      call coup2x(sw2,gal,gau,gad,gwf,gzn,gzl,gzu,gzd,g1)
      call coup3x(sw2,zmass,hmass,gwwh,gzzh,ghhh,gwwhh,gzzhh,ghhhh)
      do i=1,12
         call coup4x(sw2,zmass,fmass(i),gchf(1,i))
      enddo

c qcd couplings
      g = 1d0
      gg(1)=-g
      gg(2)=-g
      igw=0 ! don't include w width effects

!      call topwid(fmass(11),wmass,fmass(12),wwidth,igw,fwidth(11))
!      call printconstants
      return
      end
c ============================= end ============================================ 

      subroutine printconstants
c*************************************************************************
c     prints out all masses, widths, and couplings in common blocks
c*************************************************************************
      implicit none
c
c     local
c
c      integer i
c
c     global
c
      real*8          gw, gwwa, gwwz
      common /coup1/ gw, gwwa, gwwz
      real*8         gal(2),gau(2),gad(2),gwf(2)
      common /coup2a/gal,   gau,   gad,   gwf
      real*8         gzn(2),gzl(2),gzu(2),gzd(2),g1(2)
      common /coup2b/gzn,   gzl,   gzu,   gzd,   g1
      real*8         gwwh,gzzh,ghhh,gwwhh,gzzhh,ghhhh
      common /coup3/ gwwh,gzzh,ghhh,gwwhh,gzzhh,ghhhh
      complex*16     gchf(2,12)
      common /coup4/ gchf
      real*8         wmass,wwidth,zmass,zwidth
      common /vmass1/wmass,wwidth,zmass,zwidth
      real*8         amass,awidth,hmass,hwidth
      common /vmass2/amass,awidth,hmass,hwidth
      real*8            fmass(12), fwidth(12)
      common /fermions/ fmass,     fwidth
      real*8           gg(2), g
      common /coupqcd/ gg,    g

c --------------------------------- method -------------------------------------
      write(*,'(a)') 'MADGRAPH Parameters'
      write(*,'(a)') 'Boson masses and Widths:'
      write(*,10) 'W',wmass,'W',wwidth,'Z',zmass,'Z',zwidth
      write(*,10) 'A',amass,'A',awidth,'H',hmass,'H',hwidth
      write(*,'(a)') 'Quark masses and Widths:'
      write(*,10) 't',fmass(11),'t',fwidth(11),
     &            'b',fmass(12),'b',fwidth(12)


 10   format(a,1x,5Hmass=,f6.2,1H:,2x,a,1x,6Hwidth=,f6.2,1H:,2x,
     &       a,1x,5Hmass=,f6.2,1H:,2x,a,1x,6Hwidth=,f6.2)

      end
c ================================ end =========================================

      subroutine topwid(rmt,rmw,rmb,rgw,igw,rgt)
c*************************************************************************
c     the total weak decay width of the top quark, including
c     the effects of bottom mass and, if igw=1,  a finite w width.
c     from james stirling 6-10-94
c*************************************************************************
      implicit complex*16(a-h,o-z)
      real*8 rmt,rmb,rmw,xw,xb,rgw,rgt
*
      pi=4.*datan(1.d0)
      gf=1.16637d-05
      gw=cdsqrt(rmw**2*gf/dsqrt(2.d0))
*                            flavour & colour
      xb=rmb/rmt
      xw=rmw/rmt
      if(igw.eq.1) goto 10
      if(rmt.le.(rmw+rmb)) then
          write(6,*)'WARNING: mt < mw + mb !!!!'
          stop
          endif
      rgt = gf*rmt**3/8d0/pi/dsqrt(2d0) 
     .   * dsqrt( (1d0-(xw+xb)**2)*(1d0-(xw-xb)**2) )
     .   * ( (1d0-xb**2)**2 + (1d0+xb**2)*xw**2 - 2d0*xw**4 )
      return
  10  continue
      rm=xb**2
      om=1.+rm-dcmplx(rmw**2,rmw*rgw)/rmt**2
      y1=om+cdsqrt(om*om-4.*rm)
      y0=om-cdsqrt(om*om-4.*rm)
      z1=2.
      z0=2.*cdsqrt(rm)
*
      d0=(-y0**8+3.*y0**7*rm+3.*y0**7-8.*y0**6*rm-12.*y0**5*rm**
     . 2-12.*y0**5*rm+96.*y0**4*rm**2-48.*y0**3*rm**3-48.*y0**3*
     . rm**2-128.*y0**2*rm**3+192.*y0*rm**4+192.*y0*rm**3-256.*
     . rm**4)/(24.*y0**4*(y1-y0))
      d1=(-y1**8+3.*y1**7*rm+3.*y1**7-8.*y1**6*rm-12.*y1**5*rm**
     . 2-12.*y1**5*rm+96.*y1**4*rm**2-48.*y1**3*rm**3-48.*y1**3*
     . rm**2-128.*y1**2*rm**3+192.*y1*rm**4+192.*y1*rm**3-256.*
     . rm**4)/(24.*y1**4*(y1-y0))
      a4=(32.*rm**4*(y1-y0))/(3.*y1*y0*(y1-y0))
      a3=(8.*rm**3*(-3.*y1**2*y0*rm-3.*y1**2*y0+4.*y1**2*rm+3.*
     . y1*y0**2*rm+3.*y1*y0**2-4.*y0**2*rm))/(3.*y1**2*y0**2*(y1
     . -y0))
      a2=(8.*rm**3*(2.*y1**3*y0**2-3.*y1**3*y0*rm-3.*y1**3*y0+4.
     . *y1**3*rm-2.*y1**2*y0**3+3.*y1*y0**3*rm+3.*y1*y0**3-4.*y0
     . **3*rm))/(3.*y1**3*y0**3*(y1-y0))
      a1=(2.*rm**2*(3.*y1**4*y0**3*rm+3.*y1**4*y0**3+8.*y1**4*y0
     . **2*rm-12.*y1**4*y0*rm**2-12.*y1**4*y0*rm+16.*y1**4*rm**2
     . -3.*y1**3*y0**4*rm-3.*y1**3*y0**4-8.*y1**2*y0**4*rm+12.*
     . y1*y0**4*rm**2+12.*y1*y0**4*rm-16.*y0**4*rm**2))/(3.*y1**
     . 4*y0**4*(y1-y0))
      b0=(y1**3-3.*y1**2*rm-3.*y1**2+8.*y1*rm-y0**3+3.*y0**2*rm+
     . 3.*y0**2-8.*y0*rm)/(24.*(y1-y0))
      b1=(y1+y0-3.*rm-3.)/24.
      b2=1./24.
*
      rint=d0*cdlog((z1-y0)/(z0-y0))
     .    -d1*cdlog((y1-z1)/(y1-z0))
     .    -a4/3.*(1./z1**3-1./z0**3)
     .    -a3/2.*(1./z1**2-1./z0**2)
     .    -a2   *(1./z1   -1./z0   )
     .    +a1*cdlog(z1/z0)
     .    +b0   *(z1   -z0   )
     .    +b1/2.*(z1**2-z0**2)
     .    +b2/3.*(z1**3-z0**3)
*
      gw4=gw**4
*
* total width includes flavour & colour factors
      rgt=rmt**3/(rmw*rgw)*gw4/(8.*pi**3)*dimag(rint)
      rgt=9.*rgt
      return
      end
