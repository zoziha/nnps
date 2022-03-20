! 维度工厂
module ntree_factory_m

    use, intrinsic :: iso_fortran_env, only: rk => real32
    use shape_m, only: point_t, shape_t, cuboid_t, cube_t, square_t, rectangle_t, line_t, circle_t, sphere_t
    use ntree_m, only: ntree_t, make_ntree
    implicit none
    private

    public :: point_t, shape_t, ntree_t
    public :: make_boundary, make_range, make_ntree

contains

    !> 构建边界形状
    subroutine make_boundary(xmin, xmax, boundary, squared)
        real(rk), intent(in) :: xmin(:), xmax(:)
        class(shape_t), intent(out), allocatable :: boundary
        logical, intent(in) :: squared
        select case (size(xmin))
        case (1)
            boundary = line_t((xmin + xmax)*0.5_rk, xmax(1) - xmin(1))
        case (2)
            if (squared) then
                boundary = square_t((xmin + xmax)*0.5_rk, maxval(xmax - xmin))
            else
                boundary = rectangle_t((xmin + xmax)*0.5_rk, xmax - xmin)
            end if
        case (3)
            if (squared) then
                boundary = cube_t((xmin + xmax)*0.5_rk, maxval(xmax - xmin))
            else
                boundary = cuboid_t((xmin + xmax)*0.5_rk, xmax - xmin)
            end if
        end select
    end subroutine make_boundary

    !> 构建查找域边界
    subroutine make_range(point, radius, range)
        type(point_t), intent(in) :: point
        real(rk), intent(in) :: radius
        class(shape_t), intent(out), allocatable :: range
        select case (size(point%x))
        case (1)
            range = line_t(point%x, radius*2.0_rk)
        case (2)
            range = circle_t(point%x, radius)
        case (3)
            range = sphere_t(point%x, radius)
        end select
    end subroutine make_range

end module ntree_factory_m
