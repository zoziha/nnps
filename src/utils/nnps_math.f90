!> NNPS math
module nnps_math

    use nnps_kinds, only: rk
    implicit none

    private
    public :: sqrt_eps, distance1d, distance2d, distance3d

    real(rk), parameter :: sqrt_eps = sqrt(epsilon(0.0_rk))  !! eps number

contains

    !> 1d distance
    pure subroutine distance1d(x, y, r, dx)
        real(rk), intent(in) :: x, y
        real(rk), intent(out) :: r, dx

        dx = x - y
        r = abs(dx)

    end subroutine distance1d

    !> 2d distance
    pure subroutine distance2d(x, y, r, dx)
        real(rk), intent(in), dimension(2) :: x, y
        real(rk), intent(out) :: r, dx(2)

        dx = x - y
        r = sqrt(dx(1)*dx(1) + dx(2)*dx(2))

    end subroutine distance2d

    !> 3d distance
    pure subroutine distance3d(x, y, r, dx)
        real(rk), intent(in), dimension(3) :: x, y
        real(rk), intent(out) :: r, dx(3)

        dx = x - y
        r = sqrt(dx(1)*dx(1) + dx(2)*dx(2) + dx(3)*dx(3))

    end subroutine distance3d

end module nnps_math
