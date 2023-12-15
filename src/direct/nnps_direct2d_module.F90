!> 2维直接搜索
!> 2D direct neighbor search
module nnps_direct2d_module

    use nnps_kinds, only: wp
    use nnps_vector, only: vector
    use nnps_math, only: distance2d
    use, intrinsic :: iso_fortran_env, only: error_unit
#ifndef SERIAL
    use omp_lib, only: omp_get_thread_num, omp_get_max_threads
#endif
    implicit none

    private
    public :: nnps_direct2d

    !> 二维直接搜索类型
    !> nnps_direct2d
    type nnps_direct2d
        real(wp), pointer :: loc(:, :)                          !! particle 2d coordinate
        type(vector), allocatable, private :: threads_pairs(:)  !! thread local pairs
    contains
        procedure :: init, query
    end type nnps_direct2d

contains

    !> 初始化二维直接搜索
    !> initialize
    subroutine init(self, loc, n)
        class(nnps_direct2d), intent(inout), target :: self
        real(wp), intent(in), target :: loc(:, :)               !! 粒子坐标 (地址)
        integer, intent(in) :: n                                !! 粒子数量, n = size(loc, 2)
        integer :: thread_num

        self%loc => loc
#ifdef SERIAL
        thread_num = 0
#else
        thread_num = omp_get_max_threads() - 1
#endif
        allocate (self%threads_pairs(0:thread_num))
        call self%threads_pairs(:)%init(2, n)

    end subroutine init

    !> 查询粒子对
    !> query
    subroutine query(self, radius, pairs, rdxs, n)
        class(nnps_direct2d), intent(inout), target :: self
        real(wp), intent(in) :: radius                      !! query radius
        integer, pointer :: pairs(:)                        !! particle pairs
        integer, intent(in) :: n                            !! 粒子数量
        real(wp), dimension(:), pointer :: rdxs             !! particle pair distance, r, dx(2)
        integer :: i, j
        real(wp) :: rdx(3)

        self%threads_pairs(:)%len = 0
        !$omp parallel do private(i, j, rdx) schedule(dynamic)
        do i = 1, n - 1
            do j = i + 1, n

                call distance2d(self%loc(:, i), self%loc(:, j), rdx(1), rdx(2:3))
#ifdef SERIAL
                if (rdx(1) < radius) call self%threads_pairs(0)%push([i, j], rdx)
#else
                if (rdx(1) < radius) call self%threads_pairs(omp_get_thread_num())%push([i, j], rdx)
#endif

            end do
        end do

#ifndef SERIAL
        if (size(self%threads_pairs) > 1) call self%threads_pairs(0)%merge(self%threads_pairs)
#endif

        associate (len => self%threads_pairs(0)%len)

            pairs => self%threads_pairs(0)%items(1:len*2)
            rdxs => self%threads_pairs(0)%ritems(1:len*3)

        end associate

    end subroutine query

end module nnps_direct2d_module
