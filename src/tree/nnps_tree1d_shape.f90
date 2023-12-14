!> nnps binary tree shape
module nnps_tree1d_shape

    use nnps_kinds, only: wp
    use nnps_math, only: distance1d
    implicit none

    private
    public :: line

    !> line shape
    type line
        real(wp) :: left, right  !! left and right coordinate
    contains
        procedure :: contain, intersect
    end type line

contains

    !> judgy contain
    logical function contain(self, x, rdx)
        class(line), intent(in) :: self
        real(wp), intent(in) :: x
        real(wp), intent(out), optional :: rdx(2)

        if (present(rdx)) call distance1d((self%left + self%right)/2, x, rdx(1), rdx(2))
        contain = x >= self%left .and. x <= self%right

    end function contain

    !> judgy intersect
    pure logical function intersect(self, that)
        class(line), intent(in) :: self
        type(line), intent(in) :: that

        intersect = .not. (self%right <= that%left .or. self%left >= that%right)

    end function intersect

end module nnps_tree1d_shape
