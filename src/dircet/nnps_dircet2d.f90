!> 直接粒子搜索(2d)
module nnps_direct2d_module

    use nnps_kinds, only: rk
    use nnps_vector, only: vector
    use nnps_math, only: distance2d
    implicit none

    private
    public :: nnps_direct2d

    !> nnps_direct2d
    type nnps_direct2d
        real(rk), pointer :: loc(:, :)  !! particle 2d coordinate
        type(vector) :: pairs  !! partcile pairs
    contains
        procedure :: init, query
    end type nnps_direct2d

contains

    !> initialize
    subroutine init(self, loc)
        class(nnps_direct2d), intent(inout), target :: self
        real(rk), intent(in), target :: loc(:, :)

        self%loc => loc
        call self%pairs%init()

    end subroutine init

    !> query
    subroutine query(self, radius, pairs)
        class(nnps_direct2d), intent(inout), target :: self
        real(rk), intent(in) :: radius  !! query radius
        integer, pointer :: pairs(:)  !! particle pairs
        integer :: i, j
        real(rk) :: r

        self%pairs%len = 0

        do i = 1, size(self%loc, 2) - 1
            do j = i + 1, size(self%loc, 2)

                call distance2d(self%loc(:, i), self%loc(:, j), r)
                if (r < radius) then
                    call self%pairs%push(i)
                    call self%pairs%push(j)
                end if

            end do
        end do

        pairs => self%pairs%items(1:self%pairs%len)

    end subroutine query

end module nnps_direct2d_module
