function sqqff_ewp(iq,jf,p1,p2,p3,p4,lam3,lam4)

! Returns amplitude squared summed/avg over colors and helicities
! for the point in phase space p1 ,p2 ,p3 ,p4, lam3, lam4 for
!   q q  -> t t (via A ,Z ,{Z'})
! note iq (initial quark) has now been added to the arguments.
! 1 = down, 2 = up, 3 = strange, 4 = charm, 5 = bottom, 6 = top

  implicit none

! arguments
  integer :: iq,jf ! incoming quark type (up/down)
  integer :: lam3,lam4 ! ttbar helicities
  real :: p1(0:3),p2(0:3),p3(0:3),p4(0:3)
  real :: sqqff_EWp

! constants
  integer :: nexternal ! number of external legs
  integer :: ncomb ! number of helicity combinations
  parameter ( nexternal=4, ncomb= 16 )
   
! local variables
  integer :: nhel(nexternal,ncomb),ntry
  real :: t
  real :: qqff_EWp
  integer :: ihel
  logical :: goodhel(ncomb)
  data goodhel/ncomb* .FALSE. /
  data ntry/0/
!   All possible helicity combinations
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

  sqqff_EWp = 0d0
  ntry=ntry+1
  do ihel=1,ncomb
  ! f (goodhel(ihel) .or. ntry .lt. 10) then
    t=qqff_EWp(iq,jf,p1, p2, p3, p4,lam3,lam4,nhel(1,ihel))
    sqqff_EWp = sqqff_EWp + t
  !              if (t .gt. 0d0 .and. .not. goodhel(ihel)) then
  !                  goodhel(ihel)=.true.
  !                  ! write(*,*) ihel!,t
  !              endif
  ! ndif
  enddo
  sqqff_EWp = sqqff_EWp /  4d0
!       write(*,*)sqqff_EWp
end function sqqff_ewp



function qqff_ewp(iq,jf,p1,p2,p3,p4,lam3,lam4,nhel)
! returns amplitude squared summed/avg over colors
! for the point in phase space p1,p2,p3,p4
! and helicity nhel(1),nhel(2) for process :
!   qqb -> ttb (via A,Z,{Z'})

  implicit none

  real :: qqff_EWp

! Local constants
  integer ::    ngraphs ,nexternal
  parameter( ngraphs=7 ,nexternal=4 )

! Arguments
  integer :: iq,jf,lam3,lam4
  real :: p1(0:3),p2(0:3),p3(0:3),p4(0:3) ! momenta
  integer :: nhel(nexternal) ! n_hel

! Global variables
  real ::         gW, gWWA, gWWZ
  common /coup1/ gW, gWWA, gWWZ
  real ::         gAl(2),gAu(2),gAd(2),gWf(2)
  common /coup2a/gAl,   gAu,   gAd,   gWf
  real ::         gZn(2),gZl(2),gZu(2),gZd(2),g1(2)
  common /coup2b/gZn,   gZl,   gZu,   gZd,   g1
  real ::         gWWh,gZZh,ghhh,gWWhh,gZZhh,ghhhh
  common /coup3/ gWWh,gZZh,ghhh,gWWhh,gZZhh,ghhhh
  complex*16     gh(2,12)
  common /coup4/ gh
  real ::         Wmass,Wwidth,Zmass,Zwidth
  common /vmass1/Wmass,Wwidth,Zmass,Zwidth
  real ::         Amass,Awidth,hmass,hwidth
  common /vmass2/Amass,Awidth,hmass,hwidth
  real ::            fmass(12), fwidth(12)
  common /fermions/ fmass,     fwidth
! Zprime parameters
  real ::    rmZp(5),gamZp(5)
  common/Zp/rmZp   ,gamZp
  real ::         paramZp(5)
  common/Zpparam/paramZp
  real ::          gp(5),gV_d(5),gA_d(5),gV_u(5),gA_u(5)
  common/coupZpVA/gp   ,gV_d   ,gA_d   ,gV_u   ,gA_u
  real ::        gZpd(2,5),gZpu(2,5)
  common/coupZp/gZpd     ,gZpu
  integer ::     npoints
  common/stat/npoints
  integer ::       o_QCD,o_EW,o_BSM
  common/igauge/o_QCD,o_EW,o_BSM
  integer ::             o_int
  common/interference/o_int


! Local variables
  integer :: i,j
! parameter ( jf=11 ) ! final state tops
  complex*16 amp_tmp
  complex*16 amp( ngraphs )
  complex*16 w1(6) ,w2(6) ,w3(6) ,w4(6) ! external
  complex*16 w5(6) ,w6(6) ,w7(6) ! interal
  real :: gAq(2),gAf(2)
  real :: gZq(2),gZf(2)
  real :: gZpq(2,5),gZpf(2,5)
  real :: gZpq_tmp(2),gZpf_tmp(2) ! necessary to pass 2d arrays


! select only final state spins from shell script
  if((nhel(3) == lam3) .AND. (nhel(4) == lam4))then
    continue
  else
    qqff_EWp = 0.d0
    return
  end if

! up/down type couplings
  if((iq == 3) .OR. (iq == 7) .OR. (iq == 11))then
    do i=1,2
      gAq(i)=gAu(i)
      gZq(i)=gZu(i)
      do j=1,5
        gZpq(i,j)=gZpu(i,j)
      end do
    enddo
  else if((iq == 4) .OR. (iq == 8) .OR. (iq == 12))then
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

  if((jf == 1) .OR. (jf == 3) .OR. (jf == 5) .OR. &
  (jf == 7) .OR. (jf == 9) .OR. (jf == 11))then
    do i=1,2
      gAf(i)=gAu(i)
      gZf(i)=gZu(i)
      do j=1,5
        gZpf(i,j)=gZpu(i,j)
      end do
    enddo
  else if((jf == 2) .OR. (jf == 4) .OR. (jf == 6) .OR. &
    (jf == 8) .OR. (jf == 10) .OR. (jf == 12))then
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

! initialise amplitudes
  do i=1,ngraphs
    amp(i)=0d0
  enddo

! wavefunctions
  call ixxxxx( p1 ,fmass(iq) ,nhel(1) , 1   ,w1 )
  call oxxxxx( p2 ,fmass(iq) ,nhel(2) ,-1   ,w2 )
  call oxxxxx( p3 ,fmass(jf) ,nhel(3) , 1   ,w3 )
  call ixxxxx( p4 ,fmass(jf) ,nhel(4) ,-1   ,w4 )

  if (o_EW == 1)then
  ! A diagram
    call jioxxx( w1  ,w2  ,gAq ,Amass ,Awidth ,w5 )
    call iovxxx( w4  ,w3  ,w5  ,gAf   ,amp(1) )
  ! Z diagram
    call jioxxx( w1 ,w2 ,gZq ,Zmass ,Zwidth ,w6 )
    call iovxxx( w4 ,w3 ,w6  ,gZf ,amp(2) )
  else
    continue
  end if
! Z' diagrams
  if (o_BSM == 1)then
    do i =1,5
      if (rmZp(i) > 0) then
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

! total M*M for given helicity combination
  qqff_EWp = 0.d0
  amp_tmp = (0.d0,0.d0)
  if (o_int == 0)then ! no interference
    do i=1,ngraphs
      qqff_EWp = qqff_EWp+amp(i)*conjg(amp(i))
    end do
  else if (o_int == 1)then ! SM interference
    do i = 1, 2
      amp_tmp = amp_tmp + amp(i)
    end do
    qqff_EWp =qqff_EWp+amp_tmp*conjg(amp_tmp)
    do i=3,ngraphs
      qqff_EWp = qqff_EWp+amp(i)*conjg(amp(i))
    end do
  else if (o_int == 2)then ! full interference
    do i = 1, ngraphs
      amp_tmp = amp_tmp + amp(i)
    end do
    qqff_EWp =qqff_EWp+amp_tmp*conjg(amp_tmp)
  else if (o_int == 3)then ! interference only
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

! print individual amplitudes
!       if (npoints.lt.10)then
!         do i=1,ngraphs
!           write(*,*)'M: ' ,i ,amp(i)
!         enddo
!       end if

end function qqff_ewp