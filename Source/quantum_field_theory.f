module quantum_field_theory

  implicit none 

  ! SM couplings
  real :: gw, gwwa, gwwZ
  real :: gal(2),gau(2),gad(2),gwf(2)
  real :: gZn(2),gZl(2),gZu(2),gZd(2),g1(2)
  real :: gwwh,gZZh,ghhh,gwwhh,gZZhh,ghhhh
  complex*16 :: gchf(2,12)
  real :: gg(2), g

  ! fermion masses and widths
  real :: fmass(12), fwidth(12)

  ! quark masses
  real, parameter :: umass=0.d0,  cmass=0.d0,       tmass=175.d0
  real, parameter :: dmass=0.d0,  smass=0.d0,       bmass=4.18d0

  ! leptons masses
  real, parameter :: emass=0.d0,   mumass=0.d0,     taumass=1.78d0
  real, parameter :: nuemass=0d0,  numumass=0d0,    nutaumass=0d0

  ! quark widths
  real, parameter :: uwidth=0.d0,  cwidth=0.d0,     twidth=1.55d0
  real, parameter :: dwidth=0.d0,  swidth=0.d0,     bwidth=0.d0

  ! lepton widths
  real, parameter :: eWidth=0.d0,   muwidth=0.d0,   tauwidth=0.d0
  real, parameter :: nueWidth=0.d0, numuwidth=0.d0, nutauwidth=0.d0

  ! SM boson masses
  real, parameter :: rm_w=80.23d0,rm_z=91.19d0,gamma_w=2.08d0,gamma_z=2.5d0
  real, parameter :: rm_A=0d0, gamma_a=0d0, rm_h=125.d0, gamma_h=0.31278d-2

  real :: Gamma_t

! Other SM parameters
  real, parameter :: a_EM = 0.0078125, s2w = 0.2320d0
  real :: rlambdaQCD4
  integer :: nloops

! Zprime parameters
  real :: rmZp(5),gamZp(5)
  real :: paramZp(5)
  real :: gp(5),gV_d(5),gA_d(5),gV_u(5),gA_u(5)
  real :: gZpd(2,5),gZpu(2,5)

contains



end module quantum_field_theory