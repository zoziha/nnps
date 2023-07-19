!> nnps quadrature tree shape
module nnps_tree2d_shape

    use nnps_kinds, only: rk
    implicit none

    private
    public :: rectangle, circle

    !> rectangle shape
    type rectangle
        real(rk) :: left, right, top, bottom  !! left, right, top and bottom coordinate
    contains
        procedure :: contain => rectangle_contain
    end type rectangle

    !> circle shape
    type circle
        real(rk) :: center(2), radius  !! center and radius
    contains
        procedure :: intersect => circle_intersect
    end type circle

contains

    !> judgy contain
    pure logical function rectangle_contain(self, x)
        class(rectangle), intent(in) :: self
        real(rk), intent(in) :: x(2)

        rectangle_contain = x(1) >= self%left .and. x(1) <= self%right .and. &
                            x(2) >= self%bottom .and. x(2) <= self%top

    end function rectangle_contain

    !> judgy intersect
    pure logical function circle_intersect(self, that)
        class(circle), intent(in) :: self
        type(rectangle), intent(in) :: that

        if (that%contain(self%center)) then
            circle_intersect = .true.
            return
        end if

        associate (d1 => abs(self%center(1) - that%left), &
                   d2 => abs(self%center(1) - that%right), &
                   d3 => abs(self%center(2) - that%top), &
                   d4 => abs(self%center(2) - that%bottom))
            circle_intersect = min(d1, d2, d3, d4) <= self%radius
        end associate

    end function circle_intersect

end module nnps_tree2d_shape
