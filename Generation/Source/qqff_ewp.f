function sqqff_ewp(iq,jf,p1,p2,p3,p4,lam3,lam4)

  ! Returns amplitude squared summed/avg over colors and helicities
  ! for the point in phase space p1 ,p2 ,p3 ,p4, lam3, lam4 for
  !   q q  -> t t (via A ,Z ,{Z'})
  ! note iq (initial quark) has now been added to the arguments.
  ! 1 = down, 2 = up, 3 = strange, 4 = charm, 5 = bottom, 6 = top

  implicit none

  ! arguments
  integer :: iq,jf  ! incoming quark type (up/down)
  integer :: lam3,lam4  ! ttbar helicities
  real :: p1(0:3),p2(0:3),p3(0:3),p4(0:3)
  real :: sqqff_EWp

  ! constants
  integer :: nexternal  ! number of external legs
  integer :: ncomb  ! number of helicity combinations
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
    !                   ! write(*,*) ihel!,t
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

  use Configuration, only: include_EW, include_BSM, interference
  use modelling

  implicit none

  real :: qqff_EWp

  ! Local constants
  integer ::    ngraphs ,nexternal
  parameter( ngraphs=7 ,nexternal=4 )

  ! Arguments
  integer :: iq,jf,lam3,lam4
  real :: p1(0:3),p2(0:3),p3(0:3),p4(0:3)   ! momenta
  integer :: nhel(nexternal)  ! n_hel

  ! Local variables
  integer :: i,j
  ! parameter ( jf=11 )   ! final state tops
  complex*16 amp_tmp, amp_tmp2
  complex*16 amp( ngraphs )
  complex*16 w1(6) ,w2(6) ,w3(6) ,w4(6)   ! external
  complex*16 w5(6) ,w6(6) ,w7(6)  ! interal
  real :: gAq(2),gAf(2)
  real :: gZq(2),gZf(2)
  real :: gZpq(2,5),gZpf(2,5)
  real :: gZpq_tmp(2),gZpf_tmp(2)   ! necessary to pass 2d arrays


  ! select only final state spins from shell script
  if ((nhel(3) == lam3) .and. (nhel(4) == lam4)) then
    continue
  else
    qqff_EWp = 0.d0
    return
  end if

  if ((iq == 3) .or. (iq == 7)) then
    do i = 1,2
      gAq(i) = gAu(i)
      gZq(i) = gZu(i)
      do j = 1,5
        gZpq(i,j) = gZpu(i,j)
      end do
    enddo
  else if ((iq == 4) .or. (iq == 8)) then
    do i = 1,2
      gAq(i) = gAd(i)
      gZq(i) = gZd(i)
      do j = 1,5
        gZpq(i,j) = gZpd(i,j)
      end do
    enddo
  else if (iq == 11) then
    do i = 1,2
      gAq(i) = gAu(i)
      gZq(i) = gZu(i)
      do j = 1,5
        gZpq(i,j) = gZpt(i,j)
      end do
    enddo
  else if (iq == 12) then
    do i = 1,2
      gAq(i) = gAd(i)
      gZq(i) = gZd(i)
      do j = 1,5
        gZpq(i,j) = gZpb(i,j)
      end do
    enddo
  end if

  if ((jf == 3) .or. (jf == 7)) then
    do i = 1,2
      gAf(i) = gAu(i)
      gZf(i) = gZu(i)
      do j = 1,5
        gZpf(i,j) = gZpu(i,j)
      end do
    enddo
  else if ((jf == 4) .or. (jf == 8)) then
    do i = 1,2
      gAf(i) = gAd(i)
      gZf(i) = gZd(i)
      do j = 1,5
        gZpf(i,j) = gZpd(i,j)
      end do
    enddo
  else if (jf == 11) then
    do i = 1,2
      gAf(i) = gAd(i)
      gZf(i) = gZd(i)
      do j = 1,5
        gZpf(i,j) = gZpt(i,j)
      end do
    enddo
  else if (jf == 12) then
    do i = 1,2
      gAf(i) = gAd(i)
      gZf(i) = gZd(i)
      do j = 1,5
        gZpf(i,j) = gZpb(i,j)
      end do
    enddo
  end if

  if ((jf == 1) .or. (jf == 5)) then
    do i = 1,2
      gAf(i) = gAl(i)
      gZf(i) = gZl(i)
      do j = 1,5
        gZpf(i,j) = gZpl(i,j)
      end do
    enddo
  else if ((jf == 2) .or. (jf == 6)) then
    do i = 1,2
      gAf(i) = 0
      gZf(i) = gZn(i)
      do j = 1,5
        gZpf(i,j) = gZpn(i,j)
      end do
    enddo
  else if (jf == 9) then
    do i = 1,2
      gAf(i) = gAl(i)
      gZf(i) = gZl(i)
      do j = 1,5
        gZpf(i,j) = gZpl3(i,j)
      end do
    enddo
  else if (jf == 10) then
    do i = 1,2
      gAf(i) = 0
      gZf(i) = gZn(i)
      do j = 1,5
        gZpf(i,j) = gZpn3(i,j)
      end do
    enddo
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

  if (include_EW == 1) then
    ! A diagram
    call jioxxx( w1  ,w2  ,gAq ,rm_a ,gamma_A ,w5 )
    call iovxxx( w4  ,w3  ,w5  ,gAf   ,amp(1) )
    ! Z diagram
    call jioxxx( w1 ,w2 ,gZq ,rm_Z ,gamma_Z ,w6 )
    call iovxxx( w4 ,w3 ,w6  ,gZf ,amp(2) )
  else
    continue
  end if
  ! Z' diagrams
  if (include_BSM == 1) then
    do i =1,5
      if (mass_zp(i) > 0) then
        do j=1,2
          gZpq_tmp(j)=gZpq(j,i)
          gZpf_tmp(j)=gZpf(j,i)
        end do
        call jioxxx( w1 ,w2 ,gZpq_tmp ,mass_zp(i) ,gamZp(i) ,w7 )
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
  amp_tmp = (0.d0, 0.d0)
  amp_tmp2 = (0.d0, 0.d0)

  if (interference == 0) then   ! no interference
    do i = 1, ngraphs
      qqff_EWp = qqff_EWp+amp(i)*conjg(amp(i))
    end do

  else if (interference == 1) then  ! SM interference
    do i = 1, 2
      amp_tmp = amp_tmp + amp(i)
    end do
    qqff_EWp = qqff_EWp + amp_tmp*conjg(amp_tmp)
    do i = 3, ngraphs
      qqff_EWp = qqff_EWp + amp(i)*conjg(amp(i))
    end do

  else if (interference == 2) then  ! full interference
    do i = 1, ngraphs
      amp_tmp = amp_tmp + amp(i)
    end do
    qqff_EWp =qqff_EWp + amp_tmp*conjg(amp_tmp)

  else if (interference == 3) then  ! Z' and Z',SM interference only
    do i = 1, ngraphs
      amp_tmp = amp_tmp + amp(i)
    end do
    do i = 1, 2 
      amp_tmp2 = amp_tmp2 + amp(i)
    end do
    qqff_EWp = qqff_EWp + amp_tmp*conjg(amp_tmp) - amp_tmp2*conjg(amp_tmp2)

  else if (interference == 4) then  ! Z',SM interference only
    do i = 1, ngraphs
      amp_tmp = amp_tmp + amp(i)
    end do
    do i = 1, 2 
      amp_tmp2 = amp_tmp2 + amp(i)
    end do
    qqff_EWp = qqff_EWp + amp_tmp*conjg(amp_tmp) - amp_tmp2*conjg(amp_tmp2)
    do i = 3, ngraphs
      qqff_EWp = qqff_EWp - amp(i)*conjg(amp(i))
    end do
  else
    write(*,*)'Error: interference flag not set.'
    stop
  end if

  ! print individual amplitudes
  !       if (npoints.lt.10) then
  !         do i=1,ngraphs
  !           write(*,*)'M: ' ,i ,amp(i)
  !         enddo
  !       end if

end function qqff_ewp
