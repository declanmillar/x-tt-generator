c ======================================================================
      real*8 function sqqff_EWp(iq,jf,p1,p2,p3,p4,lam3,lam4)

c Returns amplitude squared summed/avg over colors and helicities
c for the point in phase space p1 ,p2 ,p3 ,p4, lam3, lam4 for
c   q q  -> t t (via A ,Z ,{Z'})
c note iq (initial quark) has now been added to the arguments.
c 1 = down, 2 = up, 3 = strange, 4 = charm, 5 = bottom, 6 = top
c ---------------------------------------------------------------------- 
      implicit none

c arguments
      integer iq,jf ! incoming quark type (up/down)
      integer lam3,lam4 ! ttbar helicities
      real*8 p1(0:3),p2(0:3),p3(0:3),p4(0:3)      

c constants
      integer nexternal ! number of external legs   
      integer ncomb ! number of helicity combinations                   
      parameter ( nexternal=4, ncomb= 16 )
 
c local variables 
      integer nhel(nexternal,ncomb),ntry
      real*8 t
      real*8 qqff_EWp
      integer ihel
      logical goodhel(ncomb)
      data goodhel/ncomb*.false./
      data ntry/0/
c   All possible helicity combinations
      data (nhel(ihel,  1),ihel=1,4) / -1, -1, -1, -1/
      data (nhel(ihel,  2),ihel=1,4) / -1, -1, -1,  1/
      data (nhel(ihel,  3),ihel=1,4) / -1, -1,  1, -1/
      data (nhel(ihel,  4),ihel=1,4) / -1, -1,  1,  1/
      data (nhel(ihel,  5),ihel=1,4) / -1,  1, -1, -1/
      data (nhel(ihel,  6),ihel=1,4) / -1,  1, -1,  1/
      data (nhel(ihel,  7),ihel=1,4) / -1,  1,  1, -1/
      data (nhel(ihel,  8),ihel=1,4) / -1,  1,  1,  1/
      data (nhel(ihel,  9),ihel=1,4) /  1, -1, -1, -1/
      data (nhel(ihel, 10),ihel=1,4) /  1, -1, -1,  1/
      data (nhel(ihel, 11),ihel=1,4) /  1, -1,  1, -1/
      data (nhel(ihel, 12),ihel=1,4) /  1, -1,  1,  1/
      data (nhel(ihel, 13),ihel=1,4) /  1,  1, -1, -1/
      data (nhel(ihel, 14),ihel=1,4) /  1,  1, -1,  1/
      data (nhel(ihel, 15),ihel=1,4) /  1,  1,  1, -1/
      data (nhel(ihel, 16),ihel=1,4) /  1,  1,  1,  1/
c ----------------------------------------------------------------------
      sqqff_EWp = 0d0
      ntry=ntry+1
      do ihel=1,ncomb
         !if (goodhel(ihel) .or. ntry .lt. 10) then
             t=qqff_EWp(iq,jf,p1, p2, p3, p4,lam3,lam4,nhel(1,ihel))
             sqqff_EWp = sqqff_EWp + t
!              if (t .gt. 0d0 .and. .not. goodhel(ihel)) then
!                  goodhel(ihel)=.true.
!                  ! write(*,*) ihel!,t
!              endif
        !endif
      enddo
      sqqff_EWp = sqqff_EWp /  4d0
!       write(*,*)sqqff_EWp
      end
c ======================================================================

c ======================================================================
      real*8 function qqff_EWp(iq,jf,p1,p2,p3,p4,lam3,lam4,nhel)
c returns amplitude squared summed/avg over colors
c for the point in phase space p1,p2,p3,p4
c and helicity nhel(1),nhel(2) for process : 
c   qqb -> ttb (via A,Z,{Z'})    
c ----------------------------------------------------------------------
      implicit none

c Local constants
      integer    ngraphs ,nexternal
      parameter( ngraphs=7 ,nexternal=4 )

c Arguments
      integer iq,jf,lam3,lam4
      real*8 p1(0:3),p2(0:3),p3(0:3),p4(0:3) ! momenta
      integer nhel(nexternal) ! n_hel

c Global variables
      real*8         gW, gWWA, gWWZ
      common /coup1/ gW, gWWA, gWWZ
      real*8         gAl(2),gAu(2),gAd(2),gWf(2)
      common /coup2a/gAl,   gAu,   gAd,   gWf
      real*8         gZn(2),gZl(2),gZu(2),gZd(2),g1(2)
      common /coup2b/gZn,   gZl,   gZu,   gZd,   g1
      real*8         gWWh,gZZh,ghhh,gWWhh,gZZhh,ghhhh
      common /coup3/ gWWh,gZZh,ghhh,gWWhh,gZZhh,ghhhh
      complex*16     gh(2,12)
      common /coup4/ gh
      real*8         Wmass,Wwidth,Zmass,Zwidth
      common /vmass1/Wmass,Wwidth,Zmass,Zwidth
      real*8         Amass,Awidth,hmass,hwidth
      common /vmass2/Amass,Awidth,hmass,hwidth
      real*8            fmass(12), fwidth(12)
      common /fermions/ fmass,     fwidth
c Zprime parameters
      real*8    rmZp(5),gamZp(5)
      common/Zp/rmZp   ,gamZp
      real*8         paramZp(5)
      common/Zpparam/paramZp
      real*8          gp(5),gV_d(5),gA_d(5),gV_u(5),gA_u(5)
      common/coupZpVA/gp   ,gV_d   ,gA_d   ,gV_u   ,gA_u
      real*8        gZpd(2,5),gZpu(2,5)
      common/coupZp/gZpd     ,gZpu
      integer     npoints
      common/stat/npoints
      integer       o_QCD,o_EW,o_BSM
      common/igauge/o_QCD,o_EW,o_BSM
      integer             o_int
      common/interference/o_int


c Local variables
      integer i,j
      ! parameter ( jf=11 ) ! final state tops
      complex*16 amp_tmp
      complex*16 amp( ngraphs )
      complex*16 w1(6) ,w2(6) ,w3(6) ,w4(6) ! external      
      complex*16 w5(6) ,w6(6) ,w7(6) ! interal
      real*8 gAq(2),gAf(2)
      real*8 gZq(2),gZf(2)
      real*8 gZpq(2,5),gZpf(2,5)
      real*8 gZpq_tmp(2),gZpf_tmp(2) ! necessary to pass 2d arrays

c ----------------------------------------------------------------------
c select only final state spins from shell script
      if((nhel(3).eq.lam3).and.(nhel(4).eq.lam4))then
        continue
      else
        qqff_EWp = 0.d0 
        return
      end if

c up/down type couplings 
      if((iq.eq. 3).or.(iq.eq. 7).or.(iq.eq.11))then
        do i=1,2
          gAq(i)=gAu(i)
          gZq(i)=gZu(i)
          do j=1,5
            gZpq(i,j)=gZpu(i,j)
          end do
        enddo
      else if((iq.eq. 4).or.(iq.eq. 8).or.(iq.eq.12))then
        do i=1,2
          gAq(i)=gAd(i)
          gZq(i)=gZd(i)
          do j=1,5
            gZpq(i,j)=gZpd(i,j)
          end do
        enddo
      else
        write(*,*)'Incorrect quark ID number.'
      end if

      if((jf.eq. 1).or.(jf.eq. 3).or.(jf.eq. 5).or.
     &   (jf.eq. 7).or.(jf.eq. 9).or.(jf.eq.11))then
        do i=1,2
          gAf(i)=gAu(i)
          gZf(i)=gZu(i)
          do j=1,5
            gZpf(i,j)=gZpu(i,j)
          end do
        enddo
      else if((jf.eq. 2).or.(jf.eq. 4).or.(jf.eq. 6).or.
     &        (jf.eq. 8).or.(jf.eq.10).or.(jf.eq.12))then
        do i=1,2
          gAf(i)=gAd(i)
          gZf(i)=gZd(i)
          do j=1,5
            gZpf(i,j)=gZpd(i,j)
          end do
        enddo
      else
        write(*,*)'Incorrect fermion ID number.'
      end if

c initialise amplitudes
      do i=1,ngraphs
        amp(i)=0d0
      enddo      

c wavefunctions
      call ixxxxx( p1 ,fmass(iq) ,nhel(1) , 1   ,w1 )                       
      call oxxxxx( p2 ,fmass(iq) ,nhel(2) ,-1   ,w2 )                       
      call oxxxxx( p3 ,fmass(jf) ,nhel(3) , 1   ,w3 )                       
      call ixxxxx( p4 ,fmass(jf) ,nhel(4) ,-1   ,w4 )

      if (o_EW.eq.1)then
c A diagram               
        call jioxxx( w1  ,w2  ,gAq ,Amass ,Awidth ,w5 )
        call iovxxx( w4  ,w3  ,w5  ,gAf   ,amp(1) )

c Z diagram
        call jioxxx( w1 ,w2 ,gZq ,Zmass ,Zwidth ,w6 )
        call iovxxx( w4 ,w3 ,w6  ,gZf ,amp(2) )
      else
        continue
      end if
c Z' diagrams
      if (o_BSM.eq.1)then
        do i =1,5
          if (rmZp(i).gt.0) then
            do j=1,2
              gZpq_tmp(j)=gZpq(j,i)
              gZpf_tmp(j)=gZpf(j,i)
            end do      
            call jioxxx( w1 ,w2 ,gZpq_tmp ,rmZp(i) ,gamZp(i) ,w7 )
            call iovxxx( w4 ,w3 ,w7   ,gZpf_tmp ,amp(2+i) )
          else
            continue
          end if 
        end do
      else
        continue  
      end if

c total M*M for given helicity combination
      qqff_EWp = 0.d0 
      amp_tmp = (0.d0,0.d0)
      if (o_int.eq.0)then ! no interference
        do i=1,ngraphs
          qqff_EWp = qqff_EWp+amp(i)*conjg(amp(i))
        end do
      else if (o_int.eq.1)then ! SM interference
        do i = 1, 2       
          amp_tmp = amp_tmp + amp(i)
        end do
        qqff_EWp =qqff_EWp+amp_tmp*conjg(amp_tmp)
        do i=3,ngraphs
          qqff_EWp = qqff_EWp+amp(i)*conjg(amp(i))
        end do
      else if (o_int.eq.2)then ! full interference
        do i = 1, ngraphs       
          amp_tmp = amp_tmp + amp(i)
        end do
        qqff_EWp =qqff_EWp+amp_tmp*conjg(amp_tmp)
      else if (o_int.eq.3)then ! interference only
        do i = 1, ngraphs       
          amp_tmp = amp_tmp + amp(i)
        end do
        qqff_EWp =qqff_EWp+amp_tmp*conjg(amp_tmp)
        do i=3,ngraphs
          qqff_EWp = qqff_EWp-amp(i)*conjg(amp(i))
        end do
      else
        write(*,*)'Error: interference flag not set.'
        stop
      end if

c print individual amplitudes
!       if (npoints.lt.10)then
!         do i=1,ngraphs
!           write(*,*)'M: ' ,i ,amp(i)
!         enddo
!       end if

      end
c ======================================================================