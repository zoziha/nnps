!> nnps binary tree shape
module nnps_tree1d_shape

    use nnps_kinds, only: rk
    implicit none

    private
    public :: line

    type line
        real(rk) :: left, right  !! left and right coordinate
    contains
        procedure :: contain, intersect
    end type line

contains

    !> judgy contain
    pure logical function contain(self, x)
        class(line), intent(in) :: self
        real(rk), intent(in) :: x

        contain = x > self%left .and. x < self%right

    end function contain

    !> judgy intersect
    pure logical function intersect(self, that)
        class(line), intent(in) :: self
        type(line), intent(in) :: that

        intersect = self%contain(that%left) .or. self%contain(that%right)

    end function intersect

end module nnps_tree1d_shape
