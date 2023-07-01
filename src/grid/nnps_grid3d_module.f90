!> background grid method
module nnps_grid3d_module

    use nnps_kinds, only: rk
    use nnps_vector, only: vector
    use nnps_math, only: distance3d, sqrt_eps
    implicit none

    private
    public :: nnps_grid3d

    !> 3d grid
    type nnps_grid3d
        real(rk), pointer :: loc(:, :)  !! particle 3d coordinate
        type(vector), allocatable :: grids(:, :, :)  !! background grids
        type(vector) :: pairs  !! particle pairs
        real(rk), dimension(3), private :: min, max
        real(rk), private :: radius
    contains
        procedure :: init, build, query
        procedure, private :: check
    end type nnps_grid3d

contains

    !> initialize
    subroutine init(self, loc, min, max, radius, len)
        class(nnps_grid3d), intent(inout) :: self
        real(rk), dimension(:, :), intent(in), target :: loc
        real(rk), dimension(3), intent(in) :: min, max
        real(rk), intent(in) :: radius
        integer, intent(in), optional :: len

        self%loc => loc
        call self%pairs%init(len)
        self%min(1:2) = min(1:2) - radius - sqrt_eps    ! setup empty grids at the boundary
        self%min(3) = min(3) - sqrt_eps
        self%max = max + radius                         ! setup empty grids at the boundary
        self%radius = radius

        associate (ik => ceiling((self%max - self%min)/radius))
            allocate (self%grids(ik(1), ik(2), ik(3)))
            call self%grids(:, :, :)%init()
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
    subroutine query(self, radius, pairs)
        class(nnps_grid3d), intent(inout), target :: self
        real(rk), intent(in) :: radius
        integer, dimension(:), pointer :: pairs
        integer :: i, j, k, l, m
        real(rk) :: r

        self%pairs%len = 0

        do k = 1, size(self%grids, 3) - 1
            do j = 2, size(self%grids, 2) - 1
                do i = 2, size(self%grids, 1) - 1

                    if (self%grids(i, j, k)%len == 0) cycle
                    do l = 1, self%grids(i, j, k)%len

                        do m = l + 1, self%grids(i, j, k)%len
                            call pairing(i, j, k, i, j, k, l, m, self%pairs)
                        end do

                        do m = 1, self%grids(i - 1, j + 1, k)%len
                            call pairing(i, j, k, i - 1, j + 1, k, l, m, self%pairs)
                        end do
                        do m = 1, self%grids(i, j + 1, k)%len
                            call pairing(i, j, k, i, j + 1, k, l, m, self%pairs)
                        end do
                        do m = 1, self%grids(i + 1, j + 1, k)%len
                            call pairing(i, j, k, i + 1, j + 1, k, l, m, self%pairs)
                        end do
                        do m = 1, self%grids(i + 1, j, k)%len
                            call pairing(i, j, k, i + 1, j, k, l, m, self%pairs)
                        end do

                        do m = 1, self%grids(i - 1, j - 1, k + 1)%len
                            call pairing(i, j, k, i - 1, j - 1, k + 1, l, m, self%pairs)
                        end do
                        do m = 1, self%grids(i, j - 1, k + 1)%len
                            call pairing(i, j, k, i, j - 1, k + 1, l, m, self%pairs)
                        end do
                        do m = 1, self%grids(i + 1, j - 1, k + 1)%len
                            call pairing(i, j, k, i + 1, j - 1, k + 1, l, m, self%pairs)
                        end do
                        do m = 1, self%grids(i - 1, j, k + 1)%len
                            call pairing(i, j, k, i - 1, j, k + 1, l, m, self%pairs)
                        end do
                        do m = 1, self%grids(i, j, k + 1)%len
                            call pairing(i, j, k, i, j, k + 1, l, m, self%pairs)
                        end do
                        do m = 1, self%grids(i + 1, j, k + 1)%len
                            call pairing(i, j, k, i + 1, j, k + 1, l, m, self%pairs)
                        end do
                        do m = 1, self%grids(i - 1, j + 1, k + 1)%len
                            call pairing(i, j, k, i - 1, j + 1, k + 1, l, m, self%pairs)
                        end do
                        do m = 1, self%grids(i, j + 1, k + 1)%len
                            call pairing(i, j, k, i, j + 1, k + 1, l, m, self%pairs)
                        end do
                        do m = 1, self%grids(i + 1, j + 1, k + 1)%len
                            call pairing(i, j, k, i + 1, j + 1, k + 1, l, m, self%pairs)
                        end do

                    end do

                end do
            end do
        end do

        pairs => self%pairs%items(1:self%pairs%len)

    contains

        pure subroutine pairing(i, j, k, ik, jk, kk, l, m, pairs)
            integer, intent(in) :: i, j, k, ik, jk, kk, l, m
            type(vector), intent(inout) :: pairs
            real(rk) :: r

            call distance3d(self%loc(:, self%grids(i, j, k)%items(l)), &
                            self%loc(:, self%grids(ik, jk, kk)%items(m)), r)
            if (r < radius) then
                call pairs%push(self%grids(i, j, k)%items(l))
                call pairs%push(self%grids(ik, jk, kk)%items(m))
            end if

        end subroutine pairing

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
