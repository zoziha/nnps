!> 3维直接搜索
!> 3D direct neighbor search
module nnps_direct3d_module

    use nnps_kinds, only: wp
    use nnps_vector, only: vector
    use nnps_math, only: distance3d
    use, intrinsic :: iso_fortran_env, only: error_unit
#ifndef SERIAL
    use omp_lib, only: omp_get_thread_num, omp_get_max_threads
#endif
    implicit none

    private
    public :: nnps_direct3d

    !> 三维直接搜索类型
    !> nnps_direct3d
    type nnps_direct3d
        real(wp), pointer :: loc(:, :)                          !! particle 3d coordinate
        type(vector), allocatable, private :: threads_pairs(:)  !! thread local pairs
        integer :: m(2)
    contains
        procedure :: init, query
    end type nnps_direct3d

contains

    !> initialize
    subroutine init(self, loc, m, n)
        class(nnps_direct3d), intent(inout) :: self
        real(wp), intent(in), target :: loc(:, :)
        integer, intent(in) :: m(2)                             !! 粒子的平均粒子对数量区间
        integer, intent(in) :: n                                !! 粒子数量, n = size(loc, 2)
        integer :: thread_num

        self%loc => loc
        self%m = m*n
#ifdef SERIAL
        thread_num = 0
#else
        thread_num = omp_get_max_threads() - 1
#endif
        allocate (self%threads_pairs(0:thread_num))
        call self%threads_pairs(:)%init(3, self%m(1)/(thread_num + 1))

    end subroutine init

    !> 查询粒子对
    !> query
    subroutine query(self, radius, pairs, rdxs, n)
        class(nnps_direct3d), intent(inout), target :: self
        real(wp), intent(in) :: radius              !! query radius
        integer, pointer :: pairs(:)                !! particle pairs
        real(wp), dimension(:), pointer :: rdxs     !! r, dx(3)
        integer, intent(in) :: n                    !! 粒子数量
        integer :: i, j
        real(wp) :: rdx(4)

        self%threads_pairs(:)%len = 0
        !$omp parallel do private(i, j, rdx) schedule(dynamic)
        do i = 1, n - 1
            do j = i + 1, n

                call distance3d(self%loc(:, i), self%loc(:, j), rdx(1), rdx(2:4))
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

            if (len*2 > self%m(2)) then
                write (error_unit, "(a)") "INFO: particle pairs exceed the maximum number"
                stop 99
            end if

            pairs => self%threads_pairs(0)%items(1:len*2)
            rdxs => self%threads_pairs(0)%ritems(1:len*4)

        end associate

    end subroutine query

end module nnps_direct3d_module
