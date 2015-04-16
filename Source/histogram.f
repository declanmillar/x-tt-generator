module class_Histogram
  implicit none
  private

  type, public :: Histogram
     real :: radius
   contains
     procedure :: bin => histogram_bin
     procedure :: print => histogram_print
  end type Histogram
contains
  function histogram_area(this) result(area)
    class(Histogram), intent(in) :: this
    real :: area
    area = pi * this%radius**2
  end function histogram_area

  subroutine histogram_print(this)
    write(10,*)'DISTRIBUTION'
    write(10,*) name
!     write(10,*) ytitle
!     write(10,*) xtitle
!     do i = 1, ndiv
!       write(10,*) x(i), fxtot(i)
!     end do
    write(10,*)'END'
  end subroutine histogram_print
end module class_Histogram


program histogram_test
  use class_Histogram
  implicit none

  type(Histogram) :: hist 
  hist = Histogram()      
  call hist%print         
end program histogram_test