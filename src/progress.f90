module progress

    use kinds

    implicit none

    ! public :: progress_percentage
    public :: progress_bar
    integer, private :: n
    integer, private, parameter :: w = 50

contains

subroutine set_total(m)
    integer, intent(in) :: m
    integer :: i
    n = m

    write(*,"(a2)", advance = "no") "["
    do i = 0, w - 1
        write(*,"(a1)", advance = "no") "-"
    end do
    write(*,"(a3, i3, a2)", advance = "no") "] (", 0, "%)"
    write(*,"(1a1)", advance = "no") char(13)
end subroutine set_total

! subroutine progress_percentage(x)
!     integer :: x, c
!     real(kind=default) :: ratio

!     if ((x .ne. n) .and. (mod(x, (n / 100 + 1)) .ne. 0)) return

!     ratio = x / real(n)
!     c = ratio * w

!     write(*,"(1a1, a11, i3, a1)", advance = "no") char(13), "progress: ", int(ratio * 100), "%"
!     if (x == n) then
!         print*, ""
!     else 
!         write(*,"(1a1)", advance = "no") char(13)
!     end if
! end subroutine progress_percentage

subroutine progress_bar(x)
    integer :: x, c, i
    real(kind=default) :: ratio

    if ((x .ne. n) .and. (mod(x, (n / 100 + 1)) .ne. 0)) return

    ratio = x / real(n)
    c = ratio * w

    write(*,"(a2)", advance = "no") "["
    do i = 0, c - 1
        write(*,"(a1)", advance = "no") "#"
    end do
    do i = c, w -1
        write(*,"(a1)", advance = "no") "-"
    end do
    write(*,"(a3, i3, a2)", advance = "no") "] (", int(ratio * 100), "%)"
    if (x == n) then
        print*, ""
    else
        write(*,"(1a1)", advance = "no") char(13)
    end if
end subroutine progress_bar

end module progress