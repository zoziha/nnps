!> 2D background grid method
module nnps_grid2d_module

    use nnps_kinds, only: rk
    use nnps_vector, only: vector
    use nnps_int_vector, only: int_vector
    use nnps_math, only: distance2d, sqrt_eps
    use omp_lib, only: omp_get_thread_num, omp_get_max_threads
    implicit none

    private
    public :: nnps_grid2d

    !> 2d grid
    type nnps_grid2d
        real(rk), pointer :: loc(:, :)  !! particle 2d coordinate
        type(int_vector), allocatable :: grids(:, :)  !! background grids
        type(vector), allocatable, private :: threads_pairs(:)  !! thread local pairs
        type(vector) :: pairs  !! particle pairs
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
        allocate (self%threads_pairs(0:omp_get_max_threads() - 1))
        call self%pairs%init(2, cap)
        call self%threads_pairs(:)%init(2, cap)
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
    subroutine query(self, radius, pairs, rdxs)
        class(nnps_grid2d), intent(inout), target :: self
        real(rk), intent(in) :: radius
        integer, dimension(:), pointer, intent(out) :: pairs
        real(rk), pointer, dimension(:), intent(out) :: rdxs
        integer :: i, j

        self%threads_pairs%len = 0

        ! U style
        !$omp parallel do private(i, j) schedule(dynamic)
        do j = 2, size(self%grids, 2)
            do i = 2, size(self%grids, 1) - 1

                if (self%grids(i, j)%len == 0) cycle

                ! L style, 4 neighbors
                call find_nearby_particles(self%grids(i, j), &                           !          ___________
                    &[self%grids(i - 1, j - 1)%items(1:self%grids(i - 1, j - 1)%len), &  !  - - -   | |     | |
                    &self%grids(i, j - 1)%items(1:self%grids(i, j - 1)%len), &           !  x o -   | |     | |
                    &self%grids(i + 1, j - 1)%items(1:self%grids(i + 1, j - 1)%len), &   !  x x x   | |_____| |
                    &self%grids(i - 1, j)%items(1:self%grids(i - 1, j)%len)], &          !          |_________|
                    &self%threads_pairs(omp_get_thread_num()))

            end do
        end do

        call self%pairs%merge(self%threads_pairs)

        pairs => self%pairs%items(1:self%pairs%len*2)
        rdxs => self%pairs%ritems(1:self%pairs%len*3)

    contains

        pure subroutine find_nearby_particles(grid, found, threads_pairs)
            type(int_vector), intent(in) :: grid
            integer, intent(in) :: found(:)
            type(vector), intent(inout) :: threads_pairs
            integer :: ii, jj
            real(rk) :: rdx(3)

            do ii = 1, grid%len

                do jj = ii + 1, grid%len
                    call distance2d(self%loc(:, grid%items(ii)), &
                                    self%loc(:, grid%items(jj)), rdx(1), rdx(2:3))
                    if (rdx(1) < radius) call threads_pairs%push([grid%items(ii), grid%items(jj)], rdx)
                end do

                do jj = 1, size(found)
                    call distance2d(self%loc(:, grid%items(ii)), &
                                    self%loc(:, found(jj)), rdx(1), rdx(2:3))
                    if (rdx(1) < radius) call threads_pairs%push([grid%items(ii), found(jj)], rdx)
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
