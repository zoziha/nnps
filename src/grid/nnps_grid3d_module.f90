!> 3D background grid method
module nnps_grid3d_module

    use nnps_kinds, only: rk
    use nnps_vector, only: vector
    use nnps_int_vector, only: int_vector
    use nnps_math, only: distance3d, sqrt_eps
    implicit none

    private
    public :: nnps_grid3d

    !> 3d grid
    type nnps_grid3d
        real(rk), pointer :: loc(:, :)  !! particle 3d coordinate
        type(int_vector), allocatable :: grids(:, :, :)  !! background grids
        type(vector) :: pairs  !! particle pairs
        type(int_vector), private :: found  !! found particles
        real(rk), dimension(3), private :: min, max
        real(rk), private :: radius
    contains
        procedure :: init, build, query
        procedure, private :: check
    end type nnps_grid3d

contains

    !> initialize
    subroutine init(self, loc, min, max, radius, cap)
        class(nnps_grid3d), intent(inout) :: self
        real(rk), dimension(:, :), intent(in), target :: loc
        real(rk), dimension(3), intent(in) :: min, max
        real(rk), intent(in) :: radius
        integer, intent(in), optional :: cap

        self%loc => loc
        call self%pairs%init(3, cap)
        call self%found%init(8*13)
        self%min = min - radius - sqrt_eps    ! setup empty grids at the boundary
        self%max(1:2) = max(1:2) + radius                         ! setup empty grids at the boundary
        self%max(3) = max(3)
        self%radius = radius

        associate (ik => ceiling((self%max - self%min)/radius))
            allocate (self%grids(ik(1), ik(2), ik(3)))
            call self%grids(:, :, :)%init(8)
        end associate

    end subroutine init

    !> build
    subroutine build(self)
        class(nnps_grid3d), intent(inout) :: self
        integer :: i

        call self%check()
        self%grids%len = 0

        do i = 1, size(self%loc, 2)
            associate (ik => ceiling((self%loc(:, i) - self%min)/self%radius))
                call self%grids(ik(1), ik(2), ik(3))%push(i)
            end associate
        end do

    end subroutine build

    !> query
    pure subroutine query(self, radius, pairs, rdxs)
        class(nnps_grid3d), intent(inout), target :: self
        real(rk), intent(in) :: radius
        integer, dimension(:), pointer, intent(out) :: pairs
        real(rk), dimension(:), pointer, intent(out) :: rdxs
        integer :: i, j, k, l, m

        self%pairs%len = 0

        ! 3D U style
        do k = 2, size(self%grids, 3)
            do j = 2, size(self%grids, 2) - 1
                do i = 2, size(self%grids, 1) - 1

                    if (self%grids(i, j, k)%len == 0) cycle
                    self%found%len = 0

                    ! 3D L style, 13 neighbors (9 + 4)
                    do m = j - 1, j + 1
                        do l = i - 1, i + 1
                            call self%found%merge(self%grids(l, m, k - 1))
                        end do
                    end do

                    call self%found%merge(self%grids(i - 1, j - 1, k))
                    call self%found%merge(self%grids(i, j - 1, k))
                    call self%found%merge(self%grids(i + 1, j - 1, k))
                    call self%found%merge(self%grids(i - 1, j, k))

                    call find_nearby_particles(self%grids(i, j, k), self%found, self%pairs)

                end do
            end do
        end do

        pairs => self%pairs%items(1:self%pairs%len*2)
        rdxs => self%pairs%ritems(1:self%pairs%len*4)

    contains

        pure subroutine find_nearby_particles(grid, found, pairs)
            type(int_vector), intent(in) :: grid, found
            type(vector), intent(inout) :: pairs
            integer :: ii, jj
            real(rk) :: rdx(4)

            do ii = 1, grid%len

                do jj = ii + 1, grid%len
                    call distance3d(self%loc(:, grid%items(ii)), &
                                    self%loc(:, grid%items(jj)), rdx(1), rdx(2:4))
                    if (rdx(1) < radius) then
                        call pairs%push([grid%items(ii), grid%items(jj)], rdx)
                    end if
                end do

                do jj = 1, found%len
                    call distance3d(self%loc(:, grid%items(ii)), &
                                    self%loc(:, found%items(jj)), rdx(1), rdx(2:4))
                    if (rdx(1) < radius) then
                        call pairs%push([grid%items(ii), found%items(jj)], rdx)
                    end if
                end do

            end do

        end subroutine find_nearby_particles

    end subroutine query

    !> check
    subroutine check(self)
        class(nnps_grid3d), intent(inout) :: self

        associate (max => maxval(self%loc, 2), &
                   min => maxval(self%loc, 2))
            if (any(max > self%max) .or. any(min < self%min)) then
                error stop 'nnps_grid3d: out of range'
            end if
        end associate

    end subroutine check

end module nnps_grid3d_module
