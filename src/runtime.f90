module runtime

    use kinds

    implicit none

    public :: print_runtime

contains

subroutine print_runtime( secs )

    real( kind = default) :: secs

    write( *, "( a10, i2.2, a1, i2.2, a1, i2.2 )" ) &
        "runtime = ", floor( secs / 3600 ), ":", floor( modulo( secs / 60 , 60.d0 ) ), ":", floor( modulo( secs, 60.d0 ) )

end subroutine print_runtime

end module runtime
