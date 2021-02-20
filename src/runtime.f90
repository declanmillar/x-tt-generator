module runtime

    use vamp_kinds

    implicit none

    public :: print_runtime

contains

subroutine print_runtime( seconds )

    real( kind = default ) :: seconds

    write( *, "( a10, i2.2, a1, i2.2, a1, i2.2 )" ) &
        "runtime = ", floor( seconds / 3600 ), ":", floor( modulo( seconds / 60 , 60.d0 ) ), ":", floor( modulo( seconds, 60.d0 ) )

end subroutine print_runtime

end module runtime
