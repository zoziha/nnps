!> 浮点型工作精度
!> Float point kind
module nnps_kinds

    implicit none

    !> 工作精度, 默认为单精度, working precision, default is single precision
#ifdef REAL64
    integer, parameter :: wp = kind(0.0d0)
#else
    integer, parameter :: wp = kind(0.0)
#endif

end module nnps_kinds
