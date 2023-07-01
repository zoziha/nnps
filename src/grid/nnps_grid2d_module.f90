!> background grid method
module nnps_grid2d_module

    use nnps_kinds, only: rk
    use nnps_vector, only: vector
    use nnps_math, only: distance2d, sqrt_eps
    implicit none

    private
    public :: nnps_grid2d

    !> 2d grid
    type nnps_grid2d
        real(rk), pointer :: loc(:, :)  !! particle 2d coordinate
        type(vector), allocatable :: grids(:, :)  !! background grids
        type(vector) :: pairs  !! particle pairs
        real(rk), dimension(2), private :: min, max
        real(rk), private :: radius
    contains
        procedure :: init, build, query
        procedure, private :: check
    end type nnps_grid2d

contains

    !> initialize
    subroutine init(self, loc, min, max, radius)
        class(nnps_grid2d), intent(inout) :: self
        real(rk), dimension(:, :), intent(in), target :: loc
        real(rk), dimension(2), intent(in) :: min, max
        real(rk), intent(in) :: radius

        self%loc => loc
        call self%pairs%init()
        self%min(1) = min(1) - radius - sqrt_eps    ! setup empty grids at the boundary
        self%min(2) = min(2) - sqrt_eps
        self%max = max + radius                     ! setup empty grids at the boundary
        self%radius = radius

        associate (ik => ceiling((self%max - self%min)/radius))
            allocate (self%grids(ik(1), ik(2)))
            call self%grids(:, :)%init()
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
    subroutine query(self, radius, pairs)
        class(nnps_grid2d), intent(inout), target :: self
        real(rk), intent(in) :: radius
        integer, dimension(:), pointer :: pairs
        integer :: i, j, k, l
        real(rk) :: r

        self%pairs%len = 0

        do l = 1, size(self%grids, 2) - 1
            do i = 2, size(self%grids, 1) - 1

                if (self%grids(i, l)%len == 0) cycle

                do j = 1, self%grids(i, l)%len
                    do k = j + 1, self%grids(i, l)%len
                        call distance2d(self%loc(:, self%grids(i, l)%items(j)), &
                                        self%loc(:, self%grids(i, l)%items(k)), r)
                        if (r < radius) then
                            call self%pairs%push(self%grids(i, l)%items(j))
                            call self%pairs%push(self%grids(i, l)%items(k))
                        end if
                    end do
                    do k = 1, self%grids(i - 1, l)%len
                        call distance2d(self%loc(:, self%grids(i, l)%items(j)), &
                                        self%loc(:, self%grids(i - 1, l)%items(k)), r)
                        if (r < radius) then
                            call self%pairs%push(self%grids(i, l)%items(j))
                            call self%pairs%push(self%grids(i - 1, l)%items(k))
                        end if
                    end do
                end do

            end do
        end do

        pairs => self%pairs%items(1:self%pairs%len)

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
