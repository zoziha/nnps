! 2 dimension
program demo

    use ntree_factory_m, only: ntree_t, shape_t, make_boundary, make_ntree, make_range, point_t
    implicit none
    class(shape_t), allocatable :: boundary, range
    type (point_t), allocatable :: found(:)
    type (ntree_t) :: quad_tree
    logical :: info
    integer :: i

    call make_boundary([0.0, 0.0], [1.0, 1.0], boundary, squared=.true.)
    call make_ntree(boundary, 2, quad_tree)
    call quad_tree%insert(point_t([0.5, 0.5]), info)
    call quad_tree%insert(point_t([-0.1, 0.1]), info)
    call quad_tree%insert(point_t([0.1, 0.2]), info)
    call quad_tree%insert(point_t([0.9, 0.4]), info)
    call quad_tree%insert(point_t([0.7, 0.6]), info)
    call quad_tree%insert(point_t([0.3, 0.5]), info)
    call quad_tree%insert(point_t([randu(), randu()]), info)
    call make_range(point_t([0.9, 0.9]), 0.5, range)
    call quad_tree%query(range, found)

    if (.not.allocated(found)) stop
    do i = 1, size(found)
        write(*, '("x = [",2(f6.2,1x),"] dist = ",f6.2)') found(i)%x(1), found(i)%x(2), &
            hypot(found(i)%x(1)-0.9, found(i)%x(2)-0.9)
    end do

contains

    function randu() result(r)
        real :: r
        call random_number(r)
    end function randu

end program demo
