ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c e' stato introdotto un nuovo common per facilitare un calcolo piu' preciso
c delle distribuzioni
c common/reslocal/resl(20),standdevl(20)
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c-----------------------------------------------------------------------
c  start vegas section
c
      block data
      implicit double precision(a-h,o-z)
c   makes default parameter assignments for vegas
      common/bveg1/ncall,itmx,nprn,ndev,xl(100),xu(100),acc
      common/bveg2/it,ndo,si,swgt,schi,xi(50,100)
      common/bveg3/alph,ndmx,mds
      data ncall/1000/,itmx/10/,nprn/0/,acc/1.d-2/,
     1     xl/100*0.d0/,xu/100*1.d0/,
     3     alph/1.5d0/,ndmx/50/,mds/1/,ndev/6/,
     4     ndo/1/,xi/5000*1.d0/,it/0/,si,swgt,schi/3*0.d0/
      end


      subroutine vegas(ndim,fxn,avgi,sd,chi2a)
c
c     !  subroutine performs ndim-dimensional monte carlo integ'n
c      - by g.p. lepage    sept 1976/(rev)aug 1979
c      - algorithm described in j comp phys 27,192(1978)
c
      implicit double precision(a-h,o-z)
      common/bveg1/ncall,itmx,nprn,ndev,xl(100),xu(100),acc
      common/bveg2/it,ndo,si,swgt,schi,xi(50,100)
      common/bveg3/alph,ndmx,mds
      common/bveg4/calls,ti,tsi
      common/reslocal/resl(20),standdevl(20)
      dimension d(50,100),di(50,100),xin(50),r(50),dx(100),ia(100),
     1          kg(100),dt(100),x(100)
      dimension rand(100)
      data one/1.d0/
c
      ndo=1
      do 1 j=1,ndim
 1    xi(1,j)=one
c
      entry vegas1(ndim,fxn,avgi,sd,chi2a)
c          initializes cummulative variables, but not grid
      it=0
      si=0.d0
      swgt=si
      schi=si
c
      entry vegas2(ndim,fxn,avgi,sd,chi2a)
c         - no initialization
      nd=ndmx
      ng=1
      if(mds.eq.0) go to 2
      ng=(ncall/2.d0)**(1.d0/ndim)
      mds=-1
      if((2*ng-ndmx).lt.0) go to 2
      mds=1
      npg=ng/ndmx+1
      nd=ng/npg
      ng=npg*nd
 2    k=ng**ndim
      npg=ncall/k
      if(npg.lt.2) npg=2
      calls=npg*k
      dxg=one/ng
      dv2g=(calls*dxg**ndim)**2/npg/npg/(npg-one)
      xnd=nd
      ndm=nd-1
      dxg=dxg*xnd
      xjac=one/calls
      do 3 j=1,ndim
      dx(j)=xu(j)-xl(j)
 3    xjac=xjac*dx(j)
c
c   rebin, preserving bin density
      if(nd.eq.ndo) go to 8
      rc=ndo/xnd
      do 7 j=1,ndim
      k=0
      xn=0.d0
      dr=xn
      i=k
 4    k=k+1
      dr=dr+one
      xo=xn
      xn=xi(k,j)
 5    if(rc.gt.dr) go to 4
      i=i+1
      dr=dr-rc
      xin(i)=xn-(xn-xo)*dr
      if(i.lt.ndm) go to 5
      do 6 i=1,ndm
 6    xi(i,j)=xin(i)
 7    xi(nd,j)=one
      ndo=nd
c
 8     if(nprn.ge.0) write(ndev,200) ndim,calls,it,itmx,acc,nprn,
     1                    alph,mds,nd,(xl(j),xu(j),j=1,ndim)
c
      entry vegas3(ndim,fxn,avgi,sd,chi2a)
c         - main integration loop
 9    it=it+1
      ti=0.d0
      tsi=ti
      do 10 j=1,ndim
      kg(j)=1
      do 10 i=1,nd
      d(i,j)=ti
 10   di(i,j)=ti
c
 11   fb=0.d0
      f2b=fb
      k=0
 12   k=k+1
      call randa(ndim,rand)
      wgt=xjac
      do 15 j=1,ndim
      xn=(kg(j)-rand(j))*dxg+one
      ia(j)=xn
      if(ia(j).gt.1) go to 13
      xo=xi(ia(j),j)
      rc=(xn-ia(j))*xo
      go to 14
 13   xo=xi(ia(j),j)-xi(ia(j)-1,j)
      rc=xi(ia(j)-1,j)+(xn-ia(j))*xo
 14   x(j)=xl(j)+rc*dx(j)
 15   wgt=wgt*xo*xnd
c
      f=wgt
      f=f*fxn(x,wgt)
      f2=f*f
      fb=fb+f
      f2b=f2b+f2
      do 16 j=1,ndim
      di(ia(j),j)=di(ia(j),j)+f
 16   if(mds.ge.0) d(ia(j),j)=d(ia(j),j)+f2
      if(k.lt.npg) go to 12
c
      f2b=sqrt(f2b*npg)
      f2b=(f2b-fb)*(f2b+fb)
      ti=ti+fb
      tsi=tsi+f2b
      if(mds.ge.0) go to 18
      do 17 j=1,ndim
 17   d(ia(j),j)=d(ia(j),j)+f2b
 18   k=ndim
 19   kg(k)=mod(kg(k),ng)+1
      if(kg(k).ne.1) go to 11
      k=k-1
      if(k.gt.0) go to 19
c
c   compute final results for this iteration
      tsi=tsi*dv2g
      ti2=ti*ti
ccccccccccccccccc modifica 1: permette integrale nullo  cccccccccccccccc`
c
      if(tsi.eq.0.) then
        wgt=0.d0
        si=0.d0
        swgt=0.d0
        schi=0.d0
        avgi=0.d0
        chi2a=0.d0
        sd=0.d0
        if(nprn.lt.0) then
          continue
        else
          tsi=sqrt(tsi) 
          write(ndev,201) it,ti,tsi,avgi,sd,chi2a       
          if(nprn.eq.0) then
            continue
          else
            do j=1,ndim
              write(ndev,202) j,(xi(i,j),di(i,j),i=1,nd,nprn)       
            end do
          end if
        end if
        return
      end if
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc`
c
c protect against accidental single precision negative values
      xxx=tsi/ti**2
      if(xxx.gt.-1.d-4) f2b=abs(f2b)

      wgt=one/tsi
      si=si+ti*wgt
      swgt=swgt+wgt
      schi=schi+ti2*wgt
      avgi=si/swgt
      chi2a=(schi-si*avgi)/(it-.9999d0)
      sd=sqrt(one/swgt)
c
      if(nprn.lt.0) go to 21
      tsi=sqrt(tsi)
      write(ndev,201) it,ti,tsi,avgi,sd,chi2a
      resl(it)=ti
      standdevl(it)=tsi
      if(nprn.eq.0) go to 21
      do 20 j=1,ndim
 20   write(ndev,202) j,(xi(i,j),di(i,j),i=1,nd,nprn)
c
c   refine grid
 21   do 23 j=1,ndim
      xo=d(1,j)
      xn=d(2,j)
      d(1,j)=(xo+xn)/2.d0
      dt(j)=d(1,j)
      do 22 i=2,ndm
      d(i,j)=xo+xn
      xo=xn
      xn=d(i+1,j)
      d(i,j)=(d(i,j)+xn)/3.d0
 22   dt(j)=dt(j)+d(i,j)
      d(nd,j)=(xn+xo)/2.d0
 23   dt(j)=dt(j)+d(nd,j)
c
      do 28 j=1,ndim
      rc=0.d0
      do 24 i=1,nd
      r(i)=0.d0
      if(d(i,j).le.0.d0) go to 24
      xo=dt(j)/d(i,j)
      r(i)=((xo-one)/xo/log(xo))**alph
 24   rc=rc+r(i)
      rc=rc/xnd
      k=0
      xn=0.d0
      dr=xn
      i=k
 25   k=k+1
      dr=dr+r(k)
      xo=xn
      xn=xi(k,j)
 26   if(rc.gt.dr) go to 25
      i=i+1
      dr=dr-rc
      xin(i)=xn-(xn-xo)*dr/r(k)
      if(i.lt.ndm) go to 26
      do 27 i=1,ndm
 27   xi(i,j)=xin(i)
 28   xi(nd,j)=one
c
      if(it.lt.itmx.and.acc*abs(avgi).lt.sd) go to 9
 200  format(/35h input parameters for vegas:  ndim=,i3,8h  ncall=,f8.0
     1  /28x,5h  it=,i5,7h  itmx=,i5/28x,6h  acc=,g9.3
     2  /28x,7h  nprn=,i3,7h  alph=,f5.2/28x,6h  mds=,i3,6h   nd=,i4
     3  /28x,10h  (xl,xu)=,(t40,2h( ,g12.6,3h , ,g12.6,2h )))
 201  format(///21h integration by vegas//14h iteration no.,i3,
     1  14h:   integral =,g14.8/21x,10hstd dev  =,g10.4/
     2  34h accumulated results:   integral =,g14.8/
     3  24x,10hstd dev  =,g10.4/24x,17hchi**2 per it'n =,g10.4)
 202  format(/15h data for axis ,i2/25h    x       delta i       ,
     1  24h   x       delta i       ,18h   x       delta i
     2  /(1h ,f7.6,1x,g11.4,5x,f7.6,1x,g11.4,5x,f7.6,1x,g11.4))
      return
      end

c***************************************************************
c load vegas data if desired
      subroutine load_vegas(ndim,name)
      implicit double precision(a-h,o-z)
      common/bveg2/it,ndo,si,swgt,schi,xi(50,100)
      character*30 name
c      open(unit=13,file=name,status='old',shared)
      read(13,210) it,ndo,si,swgt,schi
      do 190 j=1,ndim
 190  read(13,*)  (xi(i,j),i=1,50)
 210  format(2i8,3z16)
      close(unit=13)
      return
      end

c*******************************************************************
c store vegas data for possible later use
      subroutine store_vegas(ndim,name)
      implicit double precision(a-h,o-z)
      common/bveg2/it,ndo,si,swgt,schi,xi(50,100)
      character*30 name
      open(unit=12,file=name,status='new')
      write(12,210) it,ndo,si,swgt,schi
      do 190 j=1,ndim
 190  write(12,*)  (xi(i,j),i=1,50)
 210  format(2i8,3z16)
      close(unit=12)
      return
      end

      subroutine randa(n,rand)
c
c   subroutine generates uniformly distributed random no's x(i),i=1,n
      implicit double precision(a-h,o-z)
      common/rndm/iseed
      dimension rand(n)
      do 1 i=1,n
 1    rand(i)=ran2(iseed)
      return
      end

c--------- end the vegas section -------------------------------------
