module mathematics

  real, parameter :: pi = 3.14159265358979323846d0

  public :: solve_quadratic

contains

  function solve_quadratic(a, b, c)

    ! solves quadratic for both roots 
    ! returns both as complex values in a complex vector x(2)

    implicit none

    complex, dimension(2) :: solve_quadratic
    real :: a, b, c
    complex :: term1
    complex :: term2
    complex :: discriminator

    term1 = -b/(2*a)
!     print *, term1
    discriminator = b*b - 4*a*c
!     print *, discriminator
    term2 = sqrt(discriminator)/(2*a)
!     print *, term2    

    solve_quadratic(1) = term1+term2
    solve_quadratic(2) = term1-term2
    
    return
    
  end function solve_quadratic
end module mathematics