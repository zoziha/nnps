module test_1D

    use, intrinsic :: iso_fortran_env, only: rk => real32
    use testdrive, only: new_unittest, unittest_type, error_type, check
    use queue_m, only: queue_t
    use ntree_factory_m, only: ntree_t, shape_t, make_boundary, make_ntree, make_range, point_t
    use shape_m, only: line_t, point_t
    implicit none
    private

    public :: collect_1D

contains

    !> Collect all exported unit tests
    subroutine collect_1D(testsuite)
        !> Collection of tests
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                    new_unittest("make_boundary_line", test_make_boundary_line), &
                    new_unittest("make_tree", test_make_tree), &
                    new_unittest("tree_insert", test_tree_insert), &
                    new_unittest("make_range", test_make_range), &
                    new_unittest("tree_query", test_tree_query) &
                    ]

    end subroutine collect_1D

    subroutine test_make_boundary_line(error)
        type(error_type), allocatable, intent(out) :: error
        class(shape_t), allocatable :: boundary
        type(line_t) :: line1, line2

        line1 = line_t(0.0_rk, 1.0_rk)
        call check(error, line1%center(1), 0.0_rk, "line center")
        if (allocated(error)) return
        call check(error, line1%length, 1.0_rk, "line length")
        if (allocated(error)) return
        call check(error, line1%contains(point_t([0.0_rk])), .true., "line contains")
        if (allocated(error)) return
        call check(error, line1%contains(point_t([1.0_rk])), .false., "line does not contain")
        if (allocated(error)) return
        call check(error, line1%contains(point_t([-1.5_rk])), .false., "line does not contain")
        if (allocated(error)) return

        line2 = line_t(0.5_rk, 1.0_rk)
        call check(error, line1%intersects(line2), .true., "lines intersect")
        if (allocated(error)) return
        line2 = line_t(0.2_rk, 0.9_rk)
        call check(error, line1%intersects(line2), .true., "lines intersect")
        if (allocated(error)) return

        line2 = line_t(1.0_rk, 1.0_rk)
        call check(error, line1%intersects(line2), .false., "lines do not intersect")
        if (allocated(error)) return
        line2 = line_t(-1.0_rk, 1.0_rk)
        call line2%show()
        call check(error, line1%intersects(line2), .false., "lines do not intersect")
        if (allocated(error)) return

        call make_boundary([-0.5_rk], [0.5_rk], boundary, .true.)
        call boundary%show()
    end subroutine test_make_boundary_line

    subroutine test_make_tree(error)
        type(error_type), allocatable, intent(out) :: error
        type(ntree_t) :: bintree
        type(line_t) line1
        line1 = line_t(0.0_rk, 1.0_rk)
        call make_ntree(line1, 2, bintree)
        call check(error, bintree%capacity, 2, "tree capacity")
        if (allocated(error)) return
        call check(error, bintree%divided, .false., "tree not divided")
        if (allocated(error)) return
        call check(error, allocated(bintree%points), .true., "points allocated")
        if (allocated(error)) return
        call check(error, allocated(bintree%children), .false., "children not allocated")
        if (allocated(error)) return
        call bintree%boundary%show()
    end subroutine test_make_tree

    subroutine test_tree_insert(error)
        type(error_type), allocatable, intent(out) :: error
        type(ntree_t) :: bintree
        logical info
        call make_ntree(line_t(0.0_rk, 1.0_rk), 1, bintree)
        call bintree%insert(point_t([0.5_rk]), info)
        call check(error, size(bintree%points), 0, "tree size")
        if (allocated(error)) return
        call check(error, info, .false., "insert failed")
        if (allocated(error)) return

        call bintree%insert(point_t([-0.5_rk]), info)
        call check(error, size(bintree%points), 0, "tree size")
        if (allocated(error)) return
        call check(error, info, .false., "insert failed")
        if (allocated(error)) return

        call bintree%insert(point_t([0.1_rk]), info)
        call check(error, size(bintree%points), 1, "tree size")
        if (allocated(error)) return
        call check(error, info, .true., "insert succeeded")
        if (allocated(error)) return
        call check(error, bintree%points(1)%x(1), 0.1_rk, "point x")
        if (allocated(error)) return

        call bintree%insert(point_t([0.23_rk]), info)
        call check(error, allocated(bintree%children), .true., "children allocated")
        if (allocated(error)) return
        call check(error, size(bintree%children), 2, "children size")
        if (allocated(error)) return
        call check(error, bintree%children(2)%points(1)%x(1), 0.23_rk, "child point x")
        if (allocated(error)) return

        call bintree%insert(point_t([0.0_rk]), info)
        call check(error, bintree%children(1)%points(1)%x(1), 0.0_rk, "child point x")
        if (allocated(error)) return
        call check(error, bintree%children(2)%divided, .false., "child not divided")

    end subroutine test_tree_insert

    subroutine test_make_range(error)
        type(error_type), allocatable, intent(out) :: error
        class(shape_t), allocatable :: range
        call make_range([-1.0_rk], 1.0_rk, range)
        call range%show()
        select type (range)
        type is (line_t)
            call check(error, range%center(1), -1.0_rk, "range center")
            if (allocated(error)) return
            call check(error, range%length, 2.0_rk, "range length")  
        end select
    end subroutine test_make_range

    subroutine test_tree_query(error)
        type(error_type), allocatable, intent(out) :: error
        type(ntree_t) :: bintree
        type(queue_t) :: queue
        class(shape_t), allocatable :: range
        logical info
        call make_ntree(line_t(0.0_rk, 1.0_rk), 1, bintree)
        call bintree%insert(point_t([0.1_rk]), info)
        call bintree%insert(point_t([0.23_rk]), info)
        call bintree%insert(point_t([0.0_rk]), info)
        call bintree%insert(point_t([0.5_rk]), info)
        call bintree%query(line_t([0.0_rk], 0.5), queue)
        call check(error, queue%size(), 3, "queue size")
        !@todo: check queue contents
    end subroutine test_tree_query      
        

end module test_1D
