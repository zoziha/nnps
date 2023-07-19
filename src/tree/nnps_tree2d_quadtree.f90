!> quadrature tree
module nnps_tree2d_quadtree

    use nnps_kinds, only: rk
    use nnps_tree2d_shape, only: circle, rectangle
    use nnps_vector, only: vector
    use nnps_int_vector, only: int_vector
    use nnps_math, only: distance2d
    implicit none

    private
    public :: quadtree

    !> quad tree
    type quadtree
        type(rectangle) :: boundary  !! boundary
        type(quadtree), allocatable :: children(:)  !! leef
        type(int_vector) :: points  !! points
    contains
        procedure :: init, add, divide, query, clear
    end type quadtree

contains

    !> initialize
    subroutine init(self, left, right, top, bottom)
        class(quadtree), intent(inout) :: self
        real(rk), intent(in) :: left, right, top, bottom

        call self%points%init(1)
        self%boundary = rectangle(left, right, top, bottom)

    end subroutine init

    !> add point
    recursive subroutine add(self, x, i, done)
        class(quadtree), intent(inout) :: self
        real(rk), intent(in) :: x(2)
        integer, intent(in) :: i
        logical, intent(out) :: done
        integer :: j

        if (.not. self%boundary%contain(x)) return

        if (self%points%len < 1) then
            call self%points%push(i)
            done = .true.
        else
            if (.not. allocated(self%children)) call self%divide()
            do j = 1, 4
                call self%children(j)%add(x, i, done)
                if (done) return
            end do
        end if

    end subroutine add

    !> divide
    subroutine divide(self)
        class(quadtree), intent(inout) :: self

        allocate (self%children(4))
        associate (midx => (self%boundary%left + self%boundary%right)/2, &
                   midy => (self%boundary%top + self%boundary%bottom)/2)
            call self%children(1)%init(self%boundary%left, midx, self%boundary%top, midy)
            call self%children(2)%init(midx, self%boundary%right, self%boundary%top, midy)
            call self%children(3)%init(self%boundary%left, midx, midy, self%boundary%bottom)
            call self%children(4)%init(midx, self%boundary%right, midy, self%boundary%bottom)
        end associate

    end subroutine divide

    !> query
    recursive pure subroutine query(self, loc, range, i, threads_pairs)
        class(quadtree), intent(in) :: self
        real(rk), intent(in) :: loc(:, :)
        type(circle), intent(in) :: range
        integer, intent(in) :: i
        type(vector), intent(inout) :: threads_pairs
        integer :: j
        real(rk) :: rdx(3)

        if (.not. range%intersect(self%boundary)) return

        if (self%points%len > 0) then
            if (self%points%items(1) > i) then
                call distance2d(range%center, loc(:, self%points%items(1)), rdx(1), rdx(2:3))
                if (rdx(1) < range%radius) call threads_pairs%push([i, self%points%items(1)], rdx)
            end if
        end if

        if (allocated(self%children)) then
            do j = 1, 4
                call self%children(j)%query(loc, range, i, threads_pairs)
            end do
        end if

    end subroutine query

    !> clear
    recursive subroutine clear(self)
        class(quadtree), intent(inout) :: self
        integer :: i

        self%points%len = 0

        if (allocated(self%children)) then
            do i = 1, 4
                call self%children(i)%clear()
            end do
        end if

    end subroutine clear

end module nnps_tree2d_quadtree
