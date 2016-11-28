module progress

    use kinds

    implicit none

    public :: progress_percentage
    public :: progress_bar
    integer, private :: n
    integer, private, parameter :: w = 50

contains

subroutine set_total(m)
    integer, intent(in) :: m
    n = m
end subroutine set_total

subroutine progress_percentage(x)
    integer :: x, c
    real(kind=default) :: ratio

    if ((x .ne. n) .and. (mod(x, (n / 100 + 1)) .ne. 0)) return

    ratio = x / real(n)
    c = ratio * w

    print*, "progress: ", int(ratio * 100), "%";
end subroutine progress_percentage

subroutine progress_bar(x)
    integer :: x
    ! if ((x .ne. n) .and. (mod(x, (n / 100 + 1)) .ne. 0)) return

    ! float ratio = x / (float) n;
    ! unsigned int c = ratio * w;

    ! std::cout << "progress: " << std::setw(3) << (int)(ratio * 100) << '%' << '[';
    ! for (unsigned int i = 0; i < c; i++) std::cout << '=';
    ! for (unsigned int i = c; i < w; i++) std::cout << ' ';
    ! if (x == n) std::cout << '\n' << std::flush;
    ! else std::cout << ']' << '\r' << std::flush;ma
end subroutine progress_bar

end module progress