module scattering

  implicit none

  integer, parameter :: n_fb_asymmetries = 9
  real :: fac_ee, fac_emu, fac_eq, fac_qq
  real :: sigma
  real :: sigma_pol(-1:1,-1:1,20), error_pol(-1:1,-1:1,20)
  real :: sigma_fb(-1:1,9,20), error_fb(-1:1,9,20)

end module scattering