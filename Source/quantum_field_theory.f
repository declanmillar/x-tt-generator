module quantum_field_theory

  implicit none 

  ! SM couplings
  real :: gw, gwwa, gwwZ
  real :: gal(2),gau(2),gad(2),gwf(2)
  real :: gZn(2),gZl(2),gZu(2),gZd(2),g1(2)
  real :: gwwh,gZZh,ghhh,gwwhh,gZZhh,ghhhh
  complex*16 :: gchf(2,12)
  real :: gg(2), g

! SM masses and widths
  real :: fmass(12), fwidth(12)
  real :: rm_W,Gamma_W,rm_Z,Gamma_Z
  real :: rm_A,Gamma_a,rm_h,Gamma_h
  real :: Gamma_t

! Other SM parameters
  real :: a_EM,s2w
  real :: rlambdaQCD4
  integer :: nloops

! Zprime parameters
  real :: rmZp(5),gamZp(5)
  real :: paramZp(5)
  real :: gp(5),gV_d(5),gA_d(5),gV_u(5),gA_u(5)
  real :: gZpd(2,5),gZpu(2,5)

end module quantum_field_theory