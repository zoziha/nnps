!> background grid method
module nnps_grid1d_module

    use nnps_kinds, only: rk
    use nnps_vector, only: vector
    use nnps_math, only: distance1d, sqrt_eps
    implicit none

    private
    public :: nnps_grid1d

    !> 1D grid
    type nnps_grid1d
        real(rk), pointer :: loc(:)  !! particle 1d coordinate
        type(vector), allocatable :: grids(:)  !! background grids
        type(vector) :: pairs  !! particle pairs
        real(rk), private :: min, max, radius
    contains
        procedure :: init, build, query
        procedure, private :: check
    end type nnps_grid1d

contains

    !> initialize
    subroutine init(self, loc, min, max, radius)
        class(nnps_grid1d), intent(inout) :: self
        real(rk), dimension(:), intent(in), target :: loc
        real(rk), intent(in) :: min, max, radius

        self%loc => loc
        call self%pairs%init()
        associate (ik => ceiling((max - min)/radius))
            allocate (self%grids(ik))
            call self%grids(:)%init()
        end associate
        self%min = min - sqrt_eps
        self%max = max
        self%radius = radius

    end subroutine init

    !> build
    subroutine build(self)
        class(nnps_grid1d), intent(inout) :: self
        integer :: i

        call self%check()
        self%grids%len = 0

        do i = 1, size(self%loc)
            associate (ik => ceiling((self%loc(i) - self%min)/self%radius))
                call self%grids(ik)%push(i)
            end associate
        end do

    end subroutine build

    !> query
    subroutine query(self, radius, pairs)
        class(nnps_grid1d), intent(inout), target :: self
        real(rk), intent(in) :: radius
        integer, dimension(:), pointer :: pairs
        integer :: i, j, k
        real(rk) :: r

        self%pairs%len = 0

        if (self%grids(1)%len > 1) then
            do j = 1, self%grids(1)%len
                do k = j + 1, self%grids(1)%len
                    call distance1d(self%loc(self%grids(1)%items(j)), self%loc(self%grids(1)%items(k)), r)
                    if (r < radius) then
                        call self%pairs%push(self%grids(1)%items(j))
                        call self%pairs%push(self%grids(1)%items(k))
                    end if
                end do
            end do
        end if

        do i = 2, size(self%grids)
            if (self%grids(i)%len == 0) cycle
            do j = 1, self%grids(i)%len
                do k = j + 1, self%grids(i)%len
                    call distance1d(self%loc(self%grids(i)%items(j)), self%loc(self%grids(i)%items(k)), r)
                    if (r < radius) then
                        call self%pairs%push(self%grids(i)%items(j))
                        call self%pairs%push(self%grids(i)%items(k))
                    end if
                end do
                do k = 1, self%grids(i - 1)%len
                    call distance1d(self%loc(self%grids(i)%items(j)), self%loc(self%grids(i - 1)%items(k)), r)
                    if (r < radius) then
                        call self%pairs%push(self%grids(i)%items(j))
                        call self%pairs%push(self%grids(i - 1)%items(k))
                    end if
                end do
            end do
        end do

        pairs => self%pairs%items(1:self%pairs%len)

    end subroutine query

    !> check
    subroutine check(self)
        class(nnps_grid1d), intent(inout) :: self

        if (any(self%loc > self%max) .or. any(self%loc < self%min)) then
            error stop 'nnps_grid1d: out of range'
        end if

    end subroutine check

end module nnps_grid1d_module
