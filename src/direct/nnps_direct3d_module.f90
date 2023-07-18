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
    subroutine init(self, loc, cap)
        class(nnps_direct3d), intent(inout) :: self
        real(rk), intent(in), target :: loc(:, :)
        integer, intent(in), optional :: cap

        self%loc => loc
        call self%pairs%init(3, cap)

    end subroutine init

    !> query
    subroutine query(self, radius, pairs, rdxs)
        class(nnps_direct3d), intent(inout), target :: self
        real(rk), intent(in) :: radius  !! query radius
        integer, pointer :: pairs(:)  !! particle pairs
        real(rk), dimension(:), pointer :: rdxs
        integer :: i, j
        real(rk) :: rdx(4)

        self%pairs%len = 0

        !$omp parallel do private(i, j, rdx)
        do i = 1, size(self%loc, 2) - 1
            do j = i + 1, size(self%loc, 2)

                call distance3d(self%loc(:, i), self%loc(:, j), rdx(1), rdx(2:4))
                if (rdx(1) < radius) then
                    call self%pairs%push([i, j], rdx)
                end if

            end do
        end do

        pairs => self%pairs%items(1:self%pairs%len*2)
        rdxs => self%pairs%ritems(1:self%pairs%len*4)

    end subroutine query

end module nnps_direct3d_module
