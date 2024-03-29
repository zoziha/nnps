!> octave tree
module nnps_tree3d_octree

    use nnps_kinds, only: wp
    use nnps_tree3d_shape, only: sphere, cuboid
    use nnps_vector, only: vector
    use nnps_int_vector, only: int_vector
    use nnps_math, only: distance3d
    implicit none

    private
    public :: octree

    !> octave tree
    type octree
        type(cuboid) :: boundary  !! boundary
        type(octree), allocatable :: children(:)  !! leef
        type(int_vector) :: points  !! points
    contains
        procedure :: init, add, divide, query, clear
    end type octree

contains

    !> initialize
    subroutine init(self, left, right, top, bottom, front, back)
        class(octree), intent(inout) :: self
        real(wp), intent(in) :: left, right, top, bottom, front, back

        self%boundary = cuboid(left, right, top, bottom, front, back)

    end subroutine init

    !> add point
    recursive subroutine add(self, x, i, done)
        class(octree), intent(inout) :: self
        real(wp), intent(in) :: x(3)
        integer, intent(in) :: i
        logical, intent(out) :: done
        integer :: j

        if (.not. self%boundary%contain(x)) return

        if (self%points%len < 1) then
            call self%points%push_back(i)
            done = .true.
        else
            if (.not. allocated(self%children)) call self%divide()
            do j = 1, 8
                call self%children(j)%add(x, i, done)
                if (done) return
            end do
        end if

    end subroutine add

    !> divide
    subroutine divide(self)
        class(octree), intent(inout) :: self

        allocate (self%children(8))
        associate (midy => (self%boundary%left + self%boundary%right)/2, &
                   midz => (self%boundary%top + self%boundary%bottom)/2, &
                   midx => (self%boundary%front + self%boundary%back)/2)
            call self%children(1)%init(self%boundary%left, midy, self%boundary%top, midz, self%boundary%front, midx)
            call self%children(2)%init(midy, self%boundary%right, self%boundary%top, midz, self%boundary%front, midx)
            call self%children(3)%init(self%boundary%left, midy, midz, self%boundary%bottom, self%boundary%front, midx)
            call self%children(4)%init(midy, self%boundary%right, midz, self%boundary%bottom, self%boundary%front, midx)
            call self%children(5)%init(self%boundary%left, midy, self%boundary%top, midz, midx, self%boundary%back)
            call self%children(6)%init(midy, self%boundary%right, self%boundary%top, midz, midx, self%boundary%back)
            call self%children(7)%init(self%boundary%left, midy, midz, self%boundary%bottom, midx, self%boundary%back)
            call self%children(8)%init(midy, self%boundary%right, midz, self%boundary%bottom, midx, self%boundary%back)
        end associate

    end subroutine divide

    !> query
    recursive pure subroutine query(self, loc, range, i, threads_pairs)
        class(octree), intent(in) :: self
        real(wp), intent(in) :: loc(:, :)
        type(sphere), intent(in) :: range
        integer, intent(in) :: i
        type(vector), intent(inout) :: threads_pairs
        integer :: j
        real(wp) :: rdx(4)

        if (.not. range%intersect(self%boundary)) return

        if (self%points%len > 0) then
            if (self%points%items(1) > i) then
                call distance3d(range%center, loc(:, self%points%items(1)), rdx(1), rdx(2:4))
                if (rdx(1) < range%radius) call threads_pairs%push([i, self%points%items(1)], rdx)
            end if
        end if

        if (allocated(self%children)) then
            do j = 1, 8
                call self%children(j)%query(loc, range, i, threads_pairs)
            end do
        end if

    end subroutine query

    !> clear
    recursive subroutine clear(self)
        class(octree), intent(inout) :: self
        integer :: i

        self%points%len = 0

        if (allocated(self%children)) then
            do i = 1, 8
                call self%children(i)%clear()
            end do
        end if

    end subroutine clear

end module nnps_tree3d_octree
