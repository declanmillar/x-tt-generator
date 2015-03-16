module kinematics
  use configuration, only: nspat
  real :: rm3,rm4,rm5,rm6,rm7,rm8,s
  real:: pi
  real :: sigma
  parameter (pi=3.14159265358979323846d0)
  real :: unit_conv ! GeV^{-2} to nb
  parameter (unit_conv=0.38937966d9)
  real :: Xsec_polar(20,-1:1,-1:1),error_polar(20,-1:1,-1:1)
  real :: Xsec_FB(nspat,20,-1:1),error_FB(nspat,20,-1:1)
  integer :: npoints
end module Kinematics