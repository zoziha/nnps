!> 3D direct neighbor search
module nnps_direct3d_module

    use nnps_kinds, only: rk
    use nnps_vector, only: vector
    use nnps_math, only: distance3d
    implicit none

    private
    public :: nnps_direct3d

    !> nnps_direct3d
    type nnps_direct3d
        real(rk), pointer :: loc(:, :)  !! particle 3d coordinate
        type(vector) :: pairs  !! partcile pairs
    contains
        procedure :: init, query
    end type nnps_direct3d

contains

    !> initialize
    subroutine init(self, loc)
        class(nnps_direct3d), intent(inout) :: self
        real(rk), intent(in), target :: loc(:, :)

        self%loc => loc
        call self%pairs%init()

    end subroutine init

    !> query
    subroutine query(self, radius, pairs)
        class(nnps_direct3d), intent(inout), target :: self
        real(rk), intent(in) :: radius  !! query radius
        integer, pointer :: pairs(:)  !! particle pairs
        integer :: i, j
        real(rk) :: r

        self%pairs%len = 0

        do i = 1, size(self%loc, 2) - 1
            do j = i + 1, size(self%loc, 2)

                call distance3d(self%loc(:, i), self%loc(:, j), r)
                if (r < radius) then
                    call self%pairs%push(i)
                    call self%pairs%push(j)
                end if

            end do
        end do

        pairs => self%pairs%items(1:self%pairs%len)

    end subroutine query

end module nnps_direct3d_module
