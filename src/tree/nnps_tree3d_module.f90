!> 3D octave tree search
module nnps_tree3d_module

    use nnps_kinds, only: rk
    use nnps_vector, only: vector
    use nnps_math, only: distance3d
    use nnps_tree3d_octree, only: octree
    use nnps_tree3d_shape, only: sphere
    implicit none

    private
    public :: nnps_octree

    !> octave tree
    type nnps_octree
        real(rk), pointer :: loc(:, :)  !! particle 3d coordinate
        type(vector) :: pairs  !! partcile pairs
        type(octree) :: tree  !! data tree
    contains
        procedure :: init, build, query
        procedure, private :: check
    end type nnps_octree

contains

    !> initialize
    subroutine init(self, loc, min, max, cap)
        class(nnps_octree), intent(inout) :: self
        real(rk), dimension(:, :), intent(in), target :: loc
        real(rk), intent(in), dimension(3) :: min, max
        integer, intent(in), optional :: cap

        self%loc => loc
        call self%pairs%init(3, cap)
        call self%tree%init(min(1), max(1), max(2), min(2), max(3), min(3))

    end subroutine init

    !> build tree
    subroutine build(self)
        class(nnps_octree), intent(inout) :: self
        integer :: i
        logical :: done

        call self%check()
        call self%tree%clear()

        do i = 1, size(self%loc, 2)
            done = .false.
            call self%tree%add(self%loc(:, i), i, done)
            if (.not. done) error stop "nnps_tree3d: build failed"
        end do

    end subroutine build

    !> query
    subroutine query(self, radius, pairs, rdxs)
        class(nnps_octree), intent(inout), target :: self
        real(rk), intent(in) :: radius
        integer, dimension(:), pointer :: pairs
        real(rk), dimension(:), pointer :: rdxs
        integer :: i

        self%pairs%len = 0

        do i = 1, size(self%loc, 2)
            associate (range => sphere(self%loc(:, i), radius))
                call self%tree%query(self%loc, range, i, self%pairs)
            end associate
        end do

        pairs => self%pairs%items(1:self%pairs%len*2)
        rdxs => self%pairs%ritems(1:self%pairs%len*4)

    end subroutine query

    !> check
    subroutine check(self)
        class(nnps_octree), intent(in) :: self

        associate (max => maxval(self%loc, 2), &
                   min => maxval(self%loc, 2))
            if (max(1) > self%tree%boundary%right .or. &
                max(2) > self%tree%boundary%top .or. &
                min(1) < self%tree%boundary%left .or. &
                min(2) < self%tree%boundary%bottom) then
                error stop 'nnps_tree3d: out of range'
            end if
        end associate

    end subroutine check

end module nnps_tree3d_module

