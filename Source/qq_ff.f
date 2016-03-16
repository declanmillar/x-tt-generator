function sqq_ff(iq, jf, p1, p2, p3, p4, lam3, lam4)

  ! function generated by madgraph
  ! returns amplitude squared summed/avg over colors and helicities
  ! for the point in phase space p1, p2, p3, p4, lam3, lam4
  ! for process: q q~  -> A, Z, Z' -> t t~

  implicit none

  ! functions
  real :: sqq_ff
  real :: qq_ff

  ! arguments
  integer :: iq, jf
  real :: p1(0:3), p2(0:3), p3(0:3), p4(0:3)
  integer :: lam3, lam4

  ! constants
  integer, parameter :: nexternal = 4, ncomb = 16

  ! local variables
  integer :: nhel(nexternal,ncomb),ntry
  real :: t
  integer :: ihel
  logical :: goodhel(ncomb)
  data goodhel /ncomb*.false./
  data ntry /0/

  ! store possible helicity combinations
  data (nhel(ihel,  1), ihel = 1, 4) /-1, -1, -1, -1/
  data (nhel(ihel,  2), ihel = 1, 4) /-1, -1, -1,  1/
  data (nhel(ihel,  3), ihel = 1, 4) /-1, -1,  1, -1/
  data (nhel(ihel,  4), ihel = 1, 4) /-1, -1,  1,  1/
  data (nhel(ihel,  5), ihel = 1, 4) /-1,  1, -1, -1/
  data (nhel(ihel,  6), ihel = 1, 4) /-1,  1, -1,  1/
  data (nhel(ihel,  7), ihel = 1, 4) /-1,  1,  1, -1/
  data (nhel(ihel,  8), ihel = 1, 4) /-1,  1,  1,  1/
  data (nhel(ihel,  9), ihel = 1, 4) / 1, -1, -1, -1/
  data (nhel(ihel, 10), ihel = 1, 4) / 1, -1, -1,  1/
  data (nhel(ihel, 11), ihel = 1, 4) / 1, -1,  1, -1/
  data (nhel(ihel, 12), ihel = 1, 4) / 1, -1,  1,  1/
  data (nhel(ihel, 13), ihel = 1, 4) / 1,  1, -1, -1/
  data (nhel(ihel, 14), ihel = 1, 4) / 1,  1, -1,  1/
  data (nhel(ihel, 15), ihel = 1, 4) / 1,  1,  1, -1/
  data (nhel(ihel, 16), ihel = 1, 4) / 1,  1,  1,  1/

  sqq_ff = 0d0
  ntry = ntry + 1
  do ihel = 1, ncomb
    t = qq_ff(iq, jf, p1, p2, p3, p4, lam3, lam4, nhel(1, ihel))
    sqq_ff = sqq_ff + t
  enddo
  sqq_ff = sqq_ff/4d0
end function sqq_ff



function qqff_ewp(iq, jf, p1, p2, p3, p4, lam3, lam4, nhel)

  ! function generated by madgraph
  ! returns amplitude squared summed/avg over colors and helicities
  ! for the point in phase space p1, p2, p3, p4, lam3, lam4
  ! and helicity nhel(1),nhel(2)
  ! for process: q q~  -> A, Z, Z' -> t t~

  use configuration, only: include_a, include_z, include_x, interference
  use modelling

  implicit none

  ! functions
  real :: qq_ff

  ! local constants
  integer, parameter :: ngraphs = 7 ,nexternal = 4

  ! arguments
  integer :: iq, jf, lam3,lam4
  real :: p1(0:3), p2(0:3), p3(0:3), p4(0:3)
  integer :: nhel(nexternal)

  ! local variables
  integer :: i, j
  complex*16 amp_tmp, amp_tmp2
  complex*16 amp(ngraphs)
  complex*16 w1(6), w2(6), w3(6), w4(6)
  complex*16 w5(6), w6(6), w7(6)
  real :: gAq(2), gAf(2)
  real :: gZq(2), gZf(2)
  real :: gZpq(2,5), gZpf(2,5)
  real :: gZpq_tmp(2), gZpf_tmp(2)

  ! select only specified helicities
  if ((nhel(3) == lam3) .and. (nhel(4) == lam4)) then
    continue
  else
    qq_ff = 0.d0
    return
  end if

  ! initial state

  ! uu or cc
  if ((iq == 3) .or. (iq == 7)) then
    do i = 1, 2
      gAq(i) = gAu(i)
      gZq(i) = gZu(i)
      do j = 1, 5
        gZpq(i,j) = gZpu(i,j)
      end do
    enddo

  ! dd or ss
  else if ((iq == 4) .or. (iq == 8)) then
    do i = 1, 2
      gAq(i) = gAd(i)
      gZq(i) = gZd(i)
      do j = 1, 5
        gZpq(i,j) = gZpd(i,j)
      end do
    enddo

  ! tt
  else if (iq == 11) then
    do i = 1, 2
      gAq(i) = gAu(i)
      gZq(i) = gZu(i)
      do j = 1, 5
        gZpq(i,j) = gZpt(i,j)
      end do
    enddo

  ! bb
  else if (iq == 12) then
    do i = 1, 2
      gAq(i) = gAd(i)
      gZq(i) = gZd(i)
      do j = 1, 5
        gZpq(i,j) = gZpb(i,j)
      end do
    enddo
  end if

  ! uu or cc
  if ((jf == 3) .or. (jf == 7)) then
    do i = 1, 2
      gAf(i) = gAu(i)
      gZf(i) = gZu(i)
      do j = 1, 5
        gZpf(i,j) = gZpu(i,j)
      end do
    enddo

  ! dd or ss
  else if ((jf == 4) .or. (jf == 8)) then
    do i = 1, 2
      gAf(i) = gAd(i)
      gZf(i) = gZd(i)
      do j = 1, 5
        gZpf(i,j) = gZpd(i,j)
      end do
    enddo

  ! tt
  else if (jf == 11) then
    do i = 1, 2
      gAf(i) = gAd(i)
      gZf(i) = gZd(i)
      do j = 1, 5
        gZpf(i,j) = gZpt(i,j)
      end do
    enddo

  ! bb
  else if (jf == 12) then
    do i = 1, 2
      gAf(i) = gAd(i)
      gZf(i) = gZd(i)
      do j = 1, 5
        gZpf(i,j) = gZpb(i,j)
      end do
    enddo
  end if

  ! ee or mumu
  if ((jf == 1) .or. (jf == 5)) then
    do i = 1, 2
      gAf(i) = gAl(i)
      gZf(i) = gZl(i)
      do j = 1, 5
        gZpf(i,j) = gZpl(i,j)
      end do
    enddo

  ! veve or vmvm
  else if ((jf == 2) .or. (jf == 6)) then
    do i = 1, 2
      gAf(i) = 0
      gZf(i) = gZn(i)
      do j = 1, 5
        gZpf(i,j) = gZpn(i,j)
      end do
    enddo

  ! tata
  else if (jf == 9) then
    do i = 1, 2
      gAf(i) = gAl(i)
      gZf(i) = gZl(i)
      do j = 1, 5
        gZpf(i,j) = gZpl3(i,j)
      end do
    enddo

  ! vtvt
  else if (jf == 10) then
    do i = 1, 2
      gAf(i) = 0
      gZf(i) = gZn(i)
      do j = 1, 5
        gZpf(i,j) = gZpn3(i,j)
      end do
    enddo
  end if

  ! initialise amplitudes
  do i = 1, ngraphs
    amp(i) = 0.d0
  enddo

  ! wavefunctions
  call ixxxxx(p1, fmass(iq), nhel(1),  1, w1)
  call oxxxxx(p2, fmass(iq), nhel(2), -1, w2)
  call oxxxxx(p3, fmass(jf), nhel(3),  1, w3)
  call ixxxxx(p4, fmass(jf), nhel(4), -1, w4)

  if (include_a == 1) then
    call jioxxx(w1, w2, gAq, amass, awidth, w5)
    call iovxxx(w4, w3, w5, gAf, amp(1))
  end if
  if (include_z == 1) then
    call jioxxx(w1, w2, gZq, zmass, zwidth, w6)
    call iovxxx(w4, w3, w6, gZf, amp(2))
  end if

  ! Z' diagrams
  if (include_x == 1) then
    do i = 1, 5
      if (mass_zp(i) > 0) then
        do j = 1, 2
          gZpq_tmp(j) = gZpq(j,i)
          gZpf_tmp(j) = gZpf(j,i)
        end do
        call jioxxx(w1, w2, gZpq_tmp, mass_zp(i), gamZp(i), w7)
        call iovxxx(w4, w3, w7, gZpf_tmp, amp(2+i))
      end if
    end do
  end if

  qq_ff = 0.d0
  amp_tmp = (0.d0, 0.d0)
  amp_tmp2 = (0.d0, 0.d0)

  if (interference == 0) then
    ! |A|^2 + |Z|^2 + |Z'|^2 + ...
    do i = 1, ngraphs
      qq_ff = qq_ff + amp(i)*conjg(amp(i))
    end do

  else if (interference == 1) then
    ! |A + Z|^2 + |Z'|^2 + ...
    do i = 1, 2
      amp_tmp = amp_tmp + amp(i)
    end do
    qq_ff = qq_ff + amp_tmp*conjg(amp_tmp)
    do i = 3, ngraphs
      qq_ff = qq_ff + amp(i)*conjg(amp(i))
    end do

  else if (interference == 2) then
    ! |A + Z + Z' + ...|^2
    do i = 1, ngraphs
      amp_tmp = amp_tmp + amp(i)
    end do
    qq_ff = qq_ff + amp_tmp*conjg(amp_tmp)

  else if (interference == 3) then
    ! |A + Z + Z' + ...|^2 - |A + Z|^2
    do i = 1, ngraphs
      amp_tmp = amp_tmp + amp(i)
    end do
    do i = 1, 2
      amp_tmp2 = amp_tmp2 + amp(i)
    end do
    qq_ff = qq_ff + amp_tmp*conjg(amp_tmp) - amp_tmp2*conjg(amp_tmp2)

  else if (interference == 4) then
    ! |A + Z + Z' + ...|^2 - |A|^2 + |Z|^2 + |Z'|^2 + ...
    do i = 1, ngraphs
      amp_tmp = amp_tmp + amp(i)
    end do
    do i = 1, 2
      amp_tmp2 = amp_tmp2 + amp(i)
    end do
    qq_ff = qq_ff + amp_tmp*conjg(amp_tmp) - amp_tmp2*conjg(amp_tmp2)
    do i = 3, ngraphs
      qq_ff = qq_ff - amp(i)*conjg(amp(i))
    end do
  end if

end function qqff_ewp
