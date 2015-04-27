module scattering

  use configuration, only: n_fb_asymmetries

  implicit none

  real :: sigma
  real :: sigma_pol(-1:1,-1:1,20), error_pol(-1:1,-1:1,20)
  real :: sigma_fb(-1:1,n_fb_asymmetries,20), error_fb(-1:1,n_fb_asymmetries,20)
  real :: fac_ee, fac_emu, fac_eq, fac_qq

end module scattering