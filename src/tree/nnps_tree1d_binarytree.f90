!> binary tree
module nnps_tree1d_binarytree

    use nnps_kinds, only: rk
    use nnps_tree1d_shape, only: line
    use nnps_vector, only: vector
    implicit none

    private
    public :: binarytree

    !> binary tree
    type binarytree
        type(line) :: boundary  !! boundary
        type(binarytree), allocatable :: children(:)  !! leef
        type(vector) :: points  !! points
    contains
        procedure :: init, add, divide, query, clear
    end type binarytree

contains

    !> initialize
    subroutine init(self, left, right)
        class(binarytree), intent(inout) :: self
        real(rk), intent(in) :: left, right

        call self%points%init(1)
        self%boundary = line(left, right)

    end subroutine init

    !> add point
    recursive subroutine add(self, x, i, done)
        class(binarytree), intent(inout) :: self
        real(rk), intent(in) :: x
        integer, intent(in) :: i
        logical, intent(out) :: done
        integer :: j

        if (.not. self%boundary%contain(x)) return

        if (self%points%len < 1) then
            call self%points%push(i)
            done = .true.
        else
            if (.not. allocated(self%children)) call self%divide()
            do j = 1, 2
                call self%children(j)%add(x, i, done)
                if (done) return
            end do
        end if

    end subroutine add

    !> divide
    subroutine divide(self)
        class(binarytree), intent(inout) :: self

        allocate (self%children(2))
        associate (mid => (self%boundary%left + self%boundary%right)/2)
            call self%children(1)%init(self%boundary%left, mid)
            call self%children(2)%init(mid, self%boundary%right)
        end associate

    end subroutine divide

    !> query
    recursive subroutine query(self, loc, range, i, pairs)
        class(binarytree), intent(in) :: self
        real(rk), intent(in) :: loc(:)
        type(line), intent(in) :: range
        integer, intent(in) :: i
        type(vector), intent(inout) :: pairs
        integer :: j

        if (.not. range%intersect(self%boundary)) return

        if (self%points%len > 0) then
            if (self%points%items(1) > i) then
                if (range%contain(loc(self%points%items(1)))) then
                    call pairs%push(i)
                    call pairs%push(self%points%items(1))
                end if
            end if
        end if

        if (allocated(self%children)) then
            do j = 1, 2
                call self%children(j)%query(loc, range, i, pairs)
            end do
        end if

    end subroutine query

    !> clear
    recursive subroutine clear(self)
        class(binarytree), intent(inout) :: self
        integer :: i

        self%points%len = 0

        if (allocated(self%children)) then
            do i = 1, 2
                call self%children(i)%clear()
            end do
        end if

    end subroutine clear

end module nnps_tree1d_binarytree
