! 2 dimension
program demo

    use ntree_factory_m, only: ntree_t, shape_t, make_boundary, make_ntree, make_range, point_t
    use queue_m, only: queue_t
    implicit none
    class(shape_t), allocatable :: boundary, range
    type (queue_t) :: found
    class(*), allocatable :: data
    type(point_t), allocatable :: point
    type (ntree_t) :: quad_tree
    logical :: info
    integer :: i

    call make_boundary([0.0, 0.0], [1.0, 1.0], boundary, squared=.true.)
    call make_ntree(boundary, 2, quad_tree)
    call quad_tree%insert(point_t([0.5, 0.5], 1), info)
    call quad_tree%insert(point_t([-0.1, 0.1], 2), info)
    call quad_tree%insert(point_t([0.1, 0.2], 3), info)
    call quad_tree%insert(point_t([0.9, 0.4], 4), info)
    call quad_tree%insert(point_t([0.7, 0.6], 5), info)
    call quad_tree%insert(point_t([0.3, 0.5], 6), info)
    call quad_tree%insert(point_t([randu(), randu()], 7), info)
    call make_range([0.9, 0.9], 0.5, range)
    call quad_tree%query(range, found)

    if (found%size() == 0) stop
    do i = 1, found%size()
        call found%dequeue(data)
        point = to_point(data)
        write(*, '("x = [",2(f6.2,1x),"] dist = ",f6.2)') point%x(1), point%x(2), &
            hypot(point%x(1)-0.9, point%x(2)-0.9)
    end do

contains

    function randu() result(r)
        real :: r
        call random_number(r)
    end function randu
    
    pure function to_point(x) result(p)
        class(*), intent(in) :: x
        type(point_t) :: p
        select type (x)
        type is (point_t)
            p = x
        end select
    end function to_point

end program demo
