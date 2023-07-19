!> 1D binary tree search
module nnps_tree1d_module

    use nnps_kinds, only: rk
    use nnps_vector, only: vector
    use nnps_math, only: distance1d
    use nnps_tree1d_binarytree, only: binarytree
    use nnps_tree1d_shape, only: line
    implicit none

    private
    public :: nnps_binarytree

    !> binary tree
    type nnps_binarytree
        real(rk), pointer :: loc(:)  !! particle 1d coordinate
        type(vector) :: pairs  !! partcile pairs
        type(binarytree) :: tree  !! data tree
    contains
        procedure :: init, build, query
        procedure, private :: check
    end type nnps_binarytree

contains

    !> initialize
    subroutine init(self, loc, min, max, cap)
        class(nnps_binarytree), intent(inout) :: self
        real(rk), dimension(:), intent(in), target :: loc
        real(rk), intent(in) :: min, max
        integer, intent(in), optional :: cap

        self%loc => loc
        call self%pairs%init(1, cap)
        call self%tree%init(min, max)

    end subroutine init

    !> build tree
    subroutine build(self)
        class(nnps_binarytree), intent(inout) :: self
        integer :: i
        logical :: done

        call self%check()
        call self%tree%clear()

        do i = 1, size(self%loc)
            done = .false.
            call self%tree%add(self%loc(i), i, done)
            if (.not. done) error stop "nnps_tree1d: build failed"
        end do

    end subroutine build

    !> query
    subroutine query(self, radius, pairs, rdxs)
        class(nnps_binarytree), intent(inout), target :: self
        real(rk), intent(in) :: radius
        integer, dimension(:), pointer, intent(out) :: pairs
        real(rk), dimension(:), pointer, intent(out) :: rdxs
        integer :: i

        self%pairs%len = 0

        do i = 1, size(self%loc)
            call self%tree%query(self%loc, line(self%loc(i) - radius, self%loc(i) + radius), i, self%pairs)
        end do

        pairs => self%pairs%items(1:self%pairs%len*2)
        rdxs => self%pairs%ritems(1:self%pairs%len*2)

    end subroutine query

    !> check
    subroutine check(self)
        class(nnps_binarytree), intent(in) :: self

        if (any(self%loc > self%tree%boundary%right) .or. any(self%loc < self%tree%boundary%left)) then
            error stop 'nnps_tree1d: out of range'
        end if

    end subroutine check

end module nnps_tree1d_module

