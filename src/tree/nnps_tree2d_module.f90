!> 2D quadrature tree search
module nnps_tree2d_module

    use nnps_kinds, only: rk
    use nnps_vector, only: vector
    use nnps_tree2d_quadtree, only: quadtree
    use nnps_tree2d_shape, only: circle
    use omp_lib, only: omp_get_thread_num, omp_get_max_threads
    implicit none

    private
    public :: nnps_quadtree

    !> quadrature tree
    type nnps_quadtree
        real(rk), pointer :: loc(:, :)  !! particle 2d coordinate
        type(vector), allocatable, private :: threads_pairs(:)  !! thread local pairs
        type(vector) :: pairs  !! partcile pairs
        type(quadtree) :: tree  !! data tree
    contains
        procedure :: init, build, query
        procedure, private :: check
    end type nnps_quadtree

contains

    !> initialize
    subroutine init(self, loc, min, max, cap)
        class(nnps_quadtree), intent(inout) :: self
        real(rk), dimension(:, :), intent(in), target :: loc
        real(rk), intent(in), dimension(2) :: min, max
        integer, intent(in), optional :: cap

        self%loc => loc
        allocate (self%threads_pairs(0:omp_get_max_threads() - 1))
        call self%pairs%init(2, cap)
        call self%threads_pairs(:)%init(2, cap)
        call self%tree%init(min(1), max(1), max(2), min(2))

    end subroutine init

    !> build tree
    subroutine build(self)
        class(nnps_quadtree), intent(inout) :: self
        integer :: i
        logical :: done

        call self%check()
        call self%tree%clear()

        do i = 1, size(self%loc, 2)
            done = .false.
            call self%tree%add(self%loc(:, i), i, done)
            if (.not. done) error stop "nnps_tree2d: build failed"
        end do

    end subroutine build

    !> query
    subroutine query(self, radius, pairs, rdxs)
        class(nnps_quadtree), intent(inout), target :: self
        real(rk), intent(in) :: radius
        integer, dimension(:), pointer, intent(out) :: pairs
        real(rk), dimension(:), pointer, intent(out) :: rdxs
        integer :: i

        self%threads_pairs%len = 0

        !$omp parallel do private(i) schedule(dynamic)
        do i = 1, size(self%loc, 2)
            call self%tree%query(self%loc, circle(self%loc(:, i), radius), i, &
                                 self%threads_pairs(omp_get_thread_num()))
        end do

        call self%pairs%merge(self%threads_pairs)

        pairs => self%pairs%items(1:self%pairs%len*2)
        rdxs => self%pairs%ritems(1:self%pairs%len*3)

    end subroutine query

    !> check
    subroutine check(self)
        class(nnps_quadtree), intent(in) :: self

        associate (max => maxval(self%loc, 2), &
                   min => maxval(self%loc, 2))
            if (max(1) > self%tree%boundary%right .or. &
                max(2) > self%tree%boundary%top .or. &
                min(1) < self%tree%boundary%left .or. &
                min(2) < self%tree%boundary%bottom) then
                error stop 'nnps_tree2d: out of range'
            end if
        end associate

    end subroutine check

end module nnps_tree2d_module

