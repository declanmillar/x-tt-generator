module kinematics

  use configuration, only: n_fb_asymmetries

  implicit none

  real :: m3,m4,m5,m6,m7,m8
  real :: s
  real :: sigma
  integer, parameter :: pi = 3.14159265358979323846d0
  real :: unit_conv ! GeV^{-2} to nb
  parameter (unit_conv=0.38937966d9)
  real :: Xsec_polar(20,-1:1,-1:1),error_polar(20,-1:1,-1:1)
  real :: Xsec_FB(n_fb_asymmetries,20,-1:1),error_FB(n_fb_asymmetries,20,-1:1)
  integer :: npoints

  public :: reconstruct_neutrino

contains

  subroutine reconstruct_neutrino

    ! finds the longitudinal neutrino momentum for semi-hadronic decay

    

  end subroutine reconstruct_neutrino

end module kinematics
