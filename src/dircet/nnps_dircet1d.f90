!> 1D direct neighbor search
module nnps_direct1d_module

    use nnps_kinds, only: rk
    use nnps_vector, only: vector
    use nnps_math, only: distance1d
    implicit none

    private
    public :: nnps_direct1d

    !> nnps_direct1d
    type nnps_direct1d
        real(rk), pointer :: loc(:)  !! particle 1d coordinate
        type(vector) :: pairs  !! partcile pairs
    contains
        procedure :: init, query
    end type nnps_direct1d

contains

    !> initialize
    subroutine init(self, loc)
        class(nnps_direct1d), intent(inout) :: self
        real(rk), intent(in), target :: loc(:)

        self%loc => loc
        call self%pairs%init()

    end subroutine init

    !> query
    subroutine query(self, radius, pairs)
        class(nnps_direct1d), intent(inout), target :: self
        real(rk), intent(in) :: radius  !! query radius
        integer, dimension(:), pointer :: pairs
        integer :: i, j
        real(rk) :: r

        self%pairs%len = 0

        do i = 1, size(self%loc) - 1
            do j = i + 1, size(self%loc)

                call distance1d(self%loc(i), self%loc(j), r)
                if (r < radius) then
                    call self%pairs%push(i)
                    call self%pairs%push(j)
                end if

            end do
        end do

        pairs => self%pairs%items(1:self%pairs%len)

    end subroutine query

end module nnps_direct1d_module
