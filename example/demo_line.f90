! 1 dimension
program demo

    use ntree_factory_m, only: ntree_t, shape_t, make_boundary, make_ntree, make_range, point_t
    implicit none
    class(shape_t), allocatable :: boundary, range
    type (point_t), allocatable :: found(:)
    type (ntree_t) :: binarize_tree
    logical :: info
    integer :: i

    call make_boundary([0.0], [1.0], boundary, squared=.true.)
    call make_ntree(boundary, 2, binarize_tree)
    call binarize_tree%insert(point_t([0.5]), info)
    call binarize_tree%insert(point_t([-0.1]), info)
    call binarize_tree%insert(point_t([0.1]), info)
    call binarize_tree%insert(point_t([0.85]), info)
    call binarize_tree%insert(point_t([0.7]), info)
    call binarize_tree%insert(point_t([0.3]), info)
    call binarize_tree%insert(point_t([randu()]), info)
    call make_range(point=point_t([0.9]), radius=0.3, range=range)
    call binarize_tree%query(range, found)

    if (.not.allocated(found)) stop
    do i = 1, size(found)
        write(*, '("x = [",f6.2,"] dist = ",f6.2)') found(i)%x(1), abs(found(i)%x(1) - 0.9)
    end do

contains

    function randu() result(r)
        real :: r
        call random_number(r)
    end function randu

end program demo
