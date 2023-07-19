!> nnps octave tree shape
module nnps_tree3d_shape

    use nnps_kinds, only: rk
    implicit none

    private
    public :: cuboid, sphere

    !> cuboid shape
    type cuboid
        real(rk) :: left, right, top, bottom, front, back  !! left, right, top and bottom coordinate
    contains
        procedure :: contain => cuboid_contain
    end type cuboid

    !> sphere shape
    type sphere
        real(rk) :: center(3), radius  !! center and radius
    contains
        procedure :: intersect => sphere_intersect
    end type sphere

contains

    !> judgy contain
    pure logical function cuboid_contain(self, x)
        class(cuboid), intent(in) :: self
        real(rk), intent(in) :: x(3)

        cuboid_contain = x(1) >= self%left .and. x(1) <= self%right .and. &
                         x(2) >= self%bottom .and. x(2) <= self%top .and. &
                         x(3) >= self%back .and. x(3) <= self%front

    end function cuboid_contain

    !> judgy intersect
    pure logical function sphere_intersect(self, that)
        class(sphere), intent(in) :: self
        type(cuboid), intent(in) :: that

        if (that%contain(self%center)) then
            sphere_intersect = .true.
            return
        end if

        associate (d1 => abs(self%center(1) - that%left), &
                   d2 => abs(self%center(1) - that%right), &
                   d3 => abs(self%center(2) - that%top), &
                   d4 => abs(self%center(2) - that%bottom), &
                   d5 => abs(self%center(3) - that%front), &
                   d6 => abs(self%center(3) - that%back))
            sphere_intersect = min(d1, d2, d3, d4, d5, d6) <= self%radius
        end associate

    end function sphere_intersect

end module nnps_tree3d_shape
