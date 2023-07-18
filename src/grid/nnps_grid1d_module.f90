!> 1D background grid method
module nnps_grid1d_module

    use nnps_kinds, only: rk
    use nnps_vector, only: vector
    use nnps_int_vector, only: int_vector
    use nnps_math, only: distance1d, sqrt_eps
    implicit none

    private
    public :: nnps_grid1d

    !> 1D grid
    type nnps_grid1d
        real(rk), pointer :: loc(:)  !! particle 1d coordinate
        type(int_vector), allocatable :: grids(:)  !! background grids
        type(vector) :: pairs  !! particle pairs
        real(rk), private :: min, max, radius
    contains
        procedure :: init, build, query
        procedure, private :: check
    end type nnps_grid1d

contains

    !> initialize
    subroutine init(self, loc, min, max, radius, cap)
        class(nnps_grid1d), intent(inout) :: self
        real(rk), dimension(:), intent(in), target :: loc
        real(rk), intent(in) :: min, max, radius
        integer, intent(in), optional :: cap

        self%loc => loc
        call self%pairs%init(1, cap)
        self%min = min - radius - sqrt_eps
        self%max = max
        self%radius = radius

        associate (ik => ceiling((self%max - self%min)/radius))
            allocate (self%grids(ik))
            call self%grids(:)%init(8)
        end associate

    end subroutine init

    !> build
    subroutine build(self)
        class(nnps_grid1d), intent(inout) :: self
        integer :: i

        call self%check()
        self%grids(2:)%len = 0

        do i = 1, size(self%loc)
            associate (ik => ceiling((self%loc(i) - self%min)/self%radius))
                call self%grids(ik)%push(i)
            end associate
        end do

    end subroutine build

    !> query
    pure subroutine query(self, radius, pairs, rdxs)
        class(nnps_grid1d), intent(inout), target :: self
        real(rk), intent(in) :: radius
        integer, dimension(:), pointer, intent(out) :: pairs
        real(rk), pointer, dimension(:), intent(out) :: rdxs
        integer :: i, j, k
        real(rk) :: rdx(2)

        self%pairs%len = 0

        do i = 2, size(self%grids)
            if (self%grids(i)%len == 0) cycle
            do j = 1, self%grids(i)%len
                do k = j + 1, self%grids(i)%len
                    call distance1d(self%loc(self%grids(i)%items(j)), self%loc(self%grids(i)%items(k)), rdx(1), rdx(2))
                    if (rdx(1) < radius) then
                        call self%pairs%push([self%grids(i)%items(j), self%grids(i)%items(k)], rdx)
                    end if
                end do
                do k = 1, self%grids(i - 1)%len
                    call distance1d(self%loc(self%grids(i)%items(j)), self%loc(self%grids(i - 1)%items(k)), rdx(1), rdx(2))
                    if (rdx(1) < radius) then
                        call self%pairs%push([self%grids(i)%items(j), self%grids(i - 1)%items(k)], rdx)
                    end if
                end do
            end do
        end do

        pairs => self%pairs%items(1:self%pairs%len*2)
        rdxs => self%pairs%ritems(1:self%pairs%len*2)

    end subroutine query

    !> check
    subroutine check(self)
        class(nnps_grid1d), intent(inout) :: self

        if (any(self%loc > self%max) .or. any(self%loc < self%min)) then
            error stop 'nnps_grid1d: out of range'
        end if

    end subroutine check

end module nnps_grid1d_module
