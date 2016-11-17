module progress

    use kinds

    implicit none

    public :: ProgressPercentage

contains

subroutine ProgressPercentage(x, n, w)
    integer :: x, n, w, c
    real(kind=default) :: ratio

    if ((x .ne. n) .and. (mod(x, (n / 100 + 1)) .ne. 0)) return

    ratio = x / real(n)
    c = ratio * w

    print*, "progress: ", int(ratio * 100), "%";
end subroutine ProgressPercentage

end module progress