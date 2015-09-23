module kinematics

  implicit none

  real :: m3, m4, m5, m6, m7, m8
  real :: s
  real, parameter :: unit_conv = 0.38937966d9 ! GeV^{-2} to nb

end module kinematics

! function longitudinal_neutrino_momentum(p_l, pT_nu)

!   use mathematics, only: solve_quadratic
!   use modelling, only: rm_W

!   implicit none

!   ! finds the longitudinal neutrino momentum for semi-hadronic decay
!   ! assuming all particles are massless

!   real :: longitudinal_neutrino_momentum
!   real :: p_l(0:3), pT_nu(1:2), p_l0
!   complex, dimension(2) :: pz_nu
!   real :: a, b, c, k
!   integer :: i

! !   ! recalculate lepton energy in zero mass approximation
! !   p_l0 = sqrt(p_l(1)*p_l(1) + p_l(2)*p_l(2) + p_l(3)*p_l(3))

! !   ! check this matches the 4-vector energy
! !   if ( abs(p_l0 - p_l(0)) > 1.d-12) print *, "p_l0 doesn't match"

!   k = rm_W*rm_W/2 + p_l(1)*pT_nu(1) + p_l(2)*pT_nu(2)

!   a = p_l(1)*p_l(1) + p_l(2)*p_l(2)

!   b = -2*k*p_l(3)

!   c = (pT_nu(1)*pT_nu(1) + pT_nu(2)*pT_nu(2))*p_l(0)*p_l(0) - k*k

!   pz_nu = solve_quadratic(a, b, c)

!   ! select single solution
!   if (aimag(pz_nu(1)) == 0 .and. aimag(pz_nu(1)) == 0) then
!     ! two real solutions - pick smallest one
!     if (abs(real(pz_nu(1))) < abs(real(pz_nu(2)))) then
!       ! solution 1 < than solution 2
!       longitudinal_neutrino_momentum = real(pz_nu(1))
!     else if (abs(real(pz_nu(1))) > abs(real(pz_nu(2)))) then 
!       ! solution 1 > than solution 2
!       longitudinal_neutrino_momentum = real(pz_nu(2))
!     else 
!       ! solutions are equal pick 2
!       longitudinal_neutrino_momentum = real(pz_nu(2))
!     end if
!   else
!     ! no real solutions - take the real part of 1
!     longitudinal_neutrino_momentum = real(pz_nu(1))
!   end if

!   return

! end function longitudinal_neutrino_momentum

! function mass(p)

!   ! calculates the invarient length of a 4-vector

!     implicit none

!   real :: mass, mass2
!   real :: p(0:3)
!   integer :: i
!   mass2 = p(0)*p(0) 
!   do i = 1, 3
!     mass2 = mass2 - p(i)*p(i)
!   end do
!   mass = sqrt(mass2)
!   return

! end function mass