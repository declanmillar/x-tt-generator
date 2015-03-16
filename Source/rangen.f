function ran2(iseed)!declan: changed ran to ran2 for gfortran compatability
  implicit none
  integer :: iseed,nevhep,nrn(2)
  double precision :: ran2,ran1,dummy,hwrset,hwrget
  data nevhep/0/
  nevhep=nevhep+1
  nrn(1)=iseed
  nrn(2)=abs(iseed-111111111)
  if(nevhep == 1)dummy = hwrset(nrn)
  ran2=ran1(0)
  return
end function ran2

function ran1(i)

!     main random number generator
!     uses method of l'ecuyer, (via f.james, comp phys comm 60(1990)329)

  implicit none
  double precision :: ran1,hwrset,hwrget
  integer :: i,iseed(2),k,iz,jseed(2)
  save iseed
  data iseed/12345,67890/
  k=iseed(1)/53668
  iseed(1)=40014*(iseed(1)-k*53668)-k*12211
  if (iseed(1) < 0) iseed(1)=iseed(1)+2147483563
  k=iseed(2)/52774
  iseed(2)=40692*(iseed(2)-k*52774)-k*3791
  if (iseed(2) < 0) iseed(2)=iseed(2)+2147483399
  iz=iseed(1)-iseed(2)
  if (iz < 1) iz=iz+2147483562
  ran1=dble(iz)*4.656613001013252d-10
!--->                (4.656613001013252d-10 = 1.d0/2147483589)
  return
!-----------------------------------------------------------------------
  entry hwrset(jseed)
!-----------------------------------------------------------------------
  iseed(1)=jseed(1)
  iseed(2)=jseed(2)
  999 return
!-----------------------------------------------------------------------
  entry hwrget(jseed)
!-----------------------------------------------------------------------
  jseed(1)=iseed(1)
  jseed(2)=iseed(2)
  return
end function ran1
