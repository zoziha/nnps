!> 1维直接搜索
!> 1D direct neighbor search
module nnps_direct1d_module

    use nnps_kinds, only: wp
    use nnps_vector, only: vector
    use nnps_math, only: distance1d
    use, intrinsic :: iso_fortran_env, only: error_unit
    implicit none

    private
    public :: nnps_direct1d

    !> 一维直接搜索类型
    !> nnps_direct1d
    type nnps_direct1d
        real(wp), pointer :: loc(:)                 !! 粒子坐标 (指针), particle 1d coordinate
        type(vector) :: pairs                       !! 粒子对, partcile pairs
        integer :: m(2)                             !! 粒子交互粒子对数量 (不含自身) 区间 [m(1), m(2)], m 的存在有利于内存合理分配和及时发出错误提示
                                                    !! @note 数据容器依据 m(1) 分配初始空间, 依据 m(2) 判断粒子对计算是否发生异常
                                                    !! 当粒子的平均粒子对数量大于 m(2) 时, 说明粒子发生坍塌异常, 将抛出异常 "INFO: particle pairs exceed the maximum number"
    contains
        procedure :: init, query
    end type nnps_direct1d

contains

    !> 初始化一维直接搜索
    !> initialize
    subroutine init(self, loc, m, n)
        class(nnps_direct1d), intent(inout) :: self
        real(wp), intent(in), target :: loc(:)      !! 粒子坐标 (地址)
        integer, intent(in) :: m(2)                 !! 粒子的平均粒子对 (不含自身)
        integer, intent(in) :: n                    !! 粒子数量, n = size(loc)

        self%loc => loc
        self%m = m*n
        call self%pairs%init(1, self%m(1))

    end subroutine init

    !> 查询粒子对
    !> query
    subroutine query(self, radius, pairs, rdxs, n)
        class(nnps_direct1d), intent(inout), target :: self
        real(wp), intent(in) :: radius              !! 查询半径, query radius
        integer, dimension(:), pointer :: pairs     !! 粒子对, 粒子在坐标数组中的索引
        real(wp), dimension(:), pointer :: rdxs     !! 粒子对距离 (r, dx), particle pairs distance
        integer, intent(in) :: n                    !! 粒子数量
        integer :: i, j
        real(wp) :: rdx(2)

        associate (len => self%pairs%len)

            len = 0
            do i = 1, n - 1
                do j = i + 1, n

                    call distance1d(self%loc(i), self%loc(j), rdx(1), rdx(2))
                    if (rdx(1) < radius) then
                        call self%pairs%push([i, j], rdx)
                    end if

                end do
            end do

            if (len*2 > self%m(2)) then
                write (error_unit, "(a)") "INFO: particle pairs exceed the maximum number"
                stop 99
            end if

            pairs => self%pairs%items(1:len*2)
            rdxs => self%pairs%ritems(1:len*2)

        end associate

    end subroutine query

end module nnps_direct1d_module
