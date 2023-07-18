!> 2D background grid method
module nnps_grid2d_module

    use nnps_kinds, only: rk
    use nnps_vector, only: vector
    use nnps_int_vector, only: int_vector
    use nnps_math, only: distance2d, sqrt_eps
    implicit none

    private
    public :: nnps_grid2d

    !> 2d grid
    type nnps_grid2d
        real(rk), pointer :: loc(:, :)  !! particle 2d coordinate
        type(int_vector), allocatable :: grids(:, :)  !! background grids
        type(vector) :: pairs  !! particle pairs
        type(int_vector), private :: found  !! found particles
        real(rk), dimension(2), private :: min, max
        real(rk), private :: radius
    contains
        procedure :: init, build, query
        procedure, private :: check
    end type nnps_grid2d

contains

    !> initialize: U style
    subroutine init(self, loc, min, max, radius, cap)
        class(nnps_grid2d), intent(inout) :: self
        real(rk), dimension(:, :), intent(in), target :: loc
        real(rk), dimension(2), intent(in) :: min, max
        real(rk), intent(in) :: radius
        integer, intent(in), optional :: cap

        self%loc => loc
        call self%pairs%init(2, cap)
        call self%found%init(8*4)
        self%min = min - sqrt_eps - radius  ! setup empty grids at the boundary
        self%max(1) = max(1) + radius  ! setup empty grids at the boundary
        self%max(2) = max(2)
        self%radius = radius

        associate (ik => ceiling((self%max - self%min)/radius))
            allocate (self%grids(ik(1), ik(2)))
            call self%grids(:, :)%init(8)
        end associate

    end subroutine init

    !> build
    subroutine build(self)
        class(nnps_grid2d), intent(inout) :: self
        integer :: i

        call self%check()
        self%grids%len = 0

        do i = 1, size(self%loc, 2)
            associate (ik => ceiling((self%loc(:, i) - self%min)/self%radius))
                call self%grids(ik(1), ik(2))%push(i)
            end associate
        end do

    end subroutine build

    !> query
    pure subroutine query(self, radius, pairs, rdxs)
        class(nnps_grid2d), intent(inout), target :: self
        real(rk), intent(in) :: radius
        integer, dimension(:), pointer, intent(out) :: pairs
        real(rk), pointer, dimension(:), intent(out) :: rdxs
        integer :: i, j, k, l

        self%pairs%len = 0

        ! U style
        do j = 2, size(self%grids, 2)
            do i = 2, size(self%grids, 1) - 1

                if (self%grids(i, j)%len == 0) cycle
                self%found%len = 0

                ! L style, 4 neighbors                            !          ___________
                call self%found%merge(self%grids(i - 1, j - 1))   !  - - -   | |     | |
                call self%found%merge(self%grids(i, j - 1))       !  x o -   | |     | |
                call self%found%merge(self%grids(i + 1, j - 1))   !  x x x   | |_____| |
                call self%found%merge(self%grids(i - 1, j))       !          |_________|
                call find_nearby_particles(self%grids(i, j), self%found, self%pairs)

            end do
        end do

        pairs => self%pairs%items(1:self%pairs%len*2)
        rdxs => self%pairs%ritems(1:self%pairs%len*3)

    contains

        pure subroutine find_nearby_particles(grid, found, pairs)
            type(int_vector), intent(in) :: grid, found
            type(vector), intent(inout) :: pairs
            integer :: ii, jj
            real(rk) :: rdx(3)

            do ii = 1, grid%len

                do jj = ii + 1, grid%len
                    call distance2d(self%loc(:, grid%items(ii)), &
                                    self%loc(:, grid%items(jj)), rdx(1), rdx(2:3))
                    if (rdx(1) < radius) then
                        call pairs%push([grid%items(ii), grid%items(jj)], rdx)
                    end if
                end do

                do jj = 1, found%len
                    call distance2d(self%loc(:, grid%items(ii)), &
                                    self%loc(:, found%items(jj)), rdx(1), rdx(2:3))
                    if (rdx(1) < radius) then
                        call pairs%push([grid%items(ii), found%items(jj)], rdx)
                    end if
                end do

            end do

        end subroutine find_nearby_particles

    end subroutine query

    !> check
    subroutine check(self)
        class(nnps_grid2d), intent(inout) :: self

        associate (max => maxval(self%loc, 2), &
                   min => maxval(self%loc, 2))
            if (any(max > self%max) .or. any(min < self%min)) then
                error stop 'nnps_grid2d: out of range'
            end if
        end associate

    end subroutine check

end module nnps_grid2d_module
