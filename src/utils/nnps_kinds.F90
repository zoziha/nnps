!> Float point kind
module nnps_kinds

    implicit none

    private
    public :: rk

#ifdef REAL64
    integer, parameter :: rk = kind(0.0d0)
#else
    integer, parameter :: rk = kind(0.0)
#endif

end module nnps_kinds
