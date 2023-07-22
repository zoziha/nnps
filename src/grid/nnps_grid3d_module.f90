!> 3D background grid method
module nnps_grid3d_module

    use nnps_kinds, only: rk
    use nnps_vector, only: vector
    use nnps_int_vector, only: int_vector
    use nnps_math, only: distance3d, sqrt_eps
    use omp_lib, only: omp_get_thread_num, omp_get_max_threads
    implicit none

    private
    public :: nnps_grid3d

    !> 3d grid
    type nnps_grid3d
        real(rk), pointer :: loc(:, :)  !! particle 3d coordinate
        type(int_vector), allocatable :: grids(:, :, :)  !! background grids
        type(vector), allocatable, private :: threads_pairs(:)  !! thread local pairs
        type(vector) :: pairs  !! particle pairs
        real(rk), dimension(3), private :: min, max
        real(rk), private :: radius
    contains
        procedure :: init, build, query, storage
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
        allocate (self%threads_pairs(0:omp_get_max_threads() - 1))
        call self%pairs%init(3, cap)
        call self%threads_pairs(:)%init(3, cap)
        self%min = min - radius - sqrt_eps    ! setup empty grids at the boundary
        self%max(1:2) = max(1:2) + radius                         ! setup empty grids at the boundary
        self%max(3) = max(3)
        self%radius = radius

        associate (ik => ceiling((self%max - self%min)/radius))
            allocate (self%grids(ik(1), ik(2), ik(3)))
        end associate

    end subroutine init

    !> build
    subroutine build(self)
        class(nnps_grid3d), intent(inout) :: self
        integer :: i, ik(3)

        call self%check()
        self%grids%len = 0

        do i = 1, size(self%loc, 2)  ! ifort bug: cannot use `associate` in do loops
            ik = ceiling((self%loc(:, i) - self%min)/self%radius)
            call self%grids(ik(1), ik(2), ik(3))%push_back(i)
        end do

    end subroutine build

    !> query
    subroutine query(self, radius, pairs, rdxs)
        class(nnps_grid3d), intent(inout), target :: self
        real(rk), intent(in) :: radius
        integer, dimension(:), pointer, intent(out) :: pairs
        real(rk), dimension(:), pointer, intent(out) :: rdxs
        integer :: i, j, k, l, m

        self%threads_pairs%len = 0

        ! 3D U style
        !$omp parallel do private(i, j, k, l, m)
        do k = 2, size(self%grids, 3)
            do j = 2, size(self%grids, 2) - 1
                do i = 2, size(self%grids, 1) - 1

                    if (self%grids(i, j, k)%len == 0) cycle

                    ! 3D L style, 13 neighbors (9 + 4)
                    call find_nearby_particles(self%grids(i, j, k), &
                        &[self%grids(i - 1, j - 1, k - 1)%items(1:self%grids(i - 1, j - 1, k - 1)%len), &
                        & self%grids(i, j - 1, k - 1)%items(1:self%grids(i, j - 1, k - 1)%len), &
                        & self%grids(i + 1, j - 1, k - 1)%items(1:self%grids(i + 1, j - 1, k - 1)%len), &
                        & self%grids(i - 1, j, k - 1)%items(1:self%grids(i - 1, j, k - 1)%len), &
                        & self%grids(i, j, k - 1)%items(1:self%grids(i, j, k - 1)%len), &
                        & self%grids(i + 1, j, k - 1)%items(1:self%grids(i + 1, j, k - 1)%len), &
                        & self%grids(i - 1, j + 1, k - 1)%items(1:self%grids(i - 1, j + 1, k - 1)%len), &
                        & self%grids(i, j + 1, k - 1)%items(1:self%grids(i, j + 1, k - 1)%len), &
                        & self%grids(i + 1, j + 1, k - 1)%items(1:self%grids(i + 1, j + 1, k - 1)%len), &
                        & self%grids(i - 1, j - 1, k)%items(1:self%grids(i - 1, j - 1, k)%len), &
                        & self%grids(i, j - 1, k)%items(1:self%grids(i, j - 1, k)%len), &
                        & self%grids(i + 1, j - 1, k)%items(1:self%grids(i + 1, j - 1, k)%len), &
                        & self%grids(i - 1, j, k)%items(1:self%grids(i - 1, j, k)%len) &
                        ], self%threads_pairs(omp_get_thread_num()))

                end do
            end do
        end do

        call self%pairs%merge(self%threads_pairs)

        pairs => self%pairs%items(1:self%pairs%len*2)
        rdxs => self%pairs%ritems(1:self%pairs%len*4)

    contains

        pure subroutine find_nearby_particles(grid, found, threads_pairs)
            type(int_vector), intent(in) :: grid
            integer, intent(in) :: found(:)
            type(vector), intent(inout) :: threads_pairs
            integer :: ii, jj
            real(rk) :: rdx(4)

            do ii = 1, grid%len

                do jj = ii + 1, grid%len
                    call distance3d(self%loc(:, grid%items(ii)), &
                                    self%loc(:, grid%items(jj)), rdx(1), rdx(2:4))
                    if (rdx(1) < radius) call threads_pairs%push([grid%items(ii), grid%items(jj)], rdx)
                end do

                do jj = 1, size(found)
                    call distance3d(self%loc(:, grid%items(ii)), &
                                    self%loc(:, found(jj)), rdx(1), rdx(2:4))
                    if (rdx(1) < radius) call threads_pairs%push([grid%items(ii), found(jj)], rdx)
                end do

            end do

        end subroutine find_nearby_particles

    end subroutine query

    !> check
    subroutine check(self)
        class(nnps_grid3d), intent(inout) :: self
        real(rk) :: max(3), min(3)

        max = maxval(self%loc, 2)
        min = minval(self%loc, 2)
        if (any(max > self%max) .or. any(min < self%min)) then
            error stop 'nnps_grid3d: out of range'
        end if

    end subroutine check

    !> storage
    pure integer function storage(self)
        class(nnps_grid3d), intent(in) :: self
        integer :: i, j, k

        storage = storage_size(self) + storage_size(self%loc) + self%pairs%storage()
        do k = 1, size(self%grids, 3)
            do j = 1, size(self%grids, 2)
                do i = 1, size(self%grids, 1)
                    storage = storage + self%grids(i, j, k)%storage()
                end do
            end do
        end do
        do i = 1, size(self%threads_pairs)
            storage = storage + self%threads_pairs(i)%storage()
        end do

    end function storage

end module nnps_grid3d_module
