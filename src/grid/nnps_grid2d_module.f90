!> 2D background grid method
module nnps_grid2d_module

    use nnps_kinds, only: rk
    use nnps_vector, only: vector
    use nnps_int_vector, only: int_vector
    use nnps_compact_hashing, only: chash_tbl
    use nnps_math, only: distance2d, sqrt_eps
    use omp_lib, only: omp_get_thread_num, omp_get_max_threads
    implicit none

    private
    public :: nnps_grid2d

    !> 2d grid
    type nnps_grid2d
        real(rk), pointer :: loc(:, :)  !! particle 2d coordinate
        type(chash_tbl) :: tbl  !! background grids hash table
        type(vector), allocatable, private :: threads_pairs(:)  !! thread local pairs
        type(vector) :: pairs  !! particle pairs
        real(rk), dimension(2), private :: min, max
        real(rk), private :: radius
    contains
        procedure :: init, query, storage
        procedure, private :: check
    end type nnps_grid2d

contains

    !> initialize: U style
    subroutine init(self, loc, min, max, radius, n)
        class(nnps_grid2d), intent(inout) :: self
        real(rk), dimension(:, :), intent(in), target :: loc
        real(rk), dimension(2), intent(in) :: min, max
        real(rk), intent(in) :: radius
        integer, intent(in) :: n

        self%loc => loc
        allocate (self%threads_pairs(0:omp_get_max_threads() - 1))
        call self%pairs%init(2, 8*n)
        call self%threads_pairs(:)%init(2, 2*n)
        self%min = min - sqrt_eps - radius  ! setup empty grids at the boundary
        self%max(1) = max(1) + radius  ! setup empty grids at the boundary
        self%max(2) = max(2)
        self%radius = radius

        call self%tbl%allocate(m=2*n)

    end subroutine init

    !> query
    subroutine query(self, radius, pairs, rdxs, n)
        class(nnps_grid2d), intent(inout), target :: self
        integer, intent(in) :: n
        real(rk), intent(in) :: radius
        integer, dimension(:), pointer, intent(out) :: pairs
        real(rk), pointer, dimension(:), intent(out) :: rdxs
        integer :: i, j, idx(4), ik(2, n), idxij, max(2)

        call self%check()
        call self%tbl%zeroing()

        !$omp parallel do private(i)
        do i = 1, n
            ik(:, i) = ceiling((self%loc(:, i) - self%min)/self%radius)
        end do

        do i = 1, n  ! cannot use parallel do here
            call self%tbl%push(i=ik(1, i), j=ik(2, i), k=1, index=i)
        end do

        self%threads_pairs%len = 0
        max = maxval(ik, 2)
        associate (grid => self%tbl%buckets)

            ! U style
            !$omp parallel do private(i, j, idx, idxij) schedule(dynamic)
            do j = 2, max(2)
                do i = 2, max(1) - 1

                    idxij = self%tbl%hash(i, j, 1)
                    if (grid(idxij)%len == 0) cycle          !          ___________
                    ! L style, 4 neighbors                   !          |         |
                    idx = [self%tbl%hash(i - 1, j - 1, 1), & !  - - -   | |     | |
                           self%tbl%hash(i, j - 1, 1), &     !  x o -   | |     | |
                           self%tbl%hash(i + 1, j - 1, 1), & !  x x x   | |_____| |
                           self%tbl%hash(i - 1, j, 1)]       !          |_________|

                    call find_nearby_particles(grid(idxij), &
                        &[grid(idx(1))%items(1:grid(idx(1))%len), &
                        &grid(idx(2))%items(1:grid(idx(2))%len), &
                        &grid(idx(3))%items(1:grid(idx(3))%len), &
                        &grid(idx(4))%items(1:grid(idx(4))%len)], &
                        &self%threads_pairs(omp_get_thread_num()))

                end do
            end do

        end associate

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
        real(rk) :: max(2), min(2)

        max = maxval(self%loc, 2)
        min = minval(self%loc, 2)
        if (any(max > self%max) .or. any(min < self%min)) then
            error stop 'nnps_grid2d: out of range'
        end if

    end subroutine check

    !> storage
    integer function storage(self)
        class(nnps_grid2d), intent(in) :: self
        integer :: i, j

        storage = storage_size(self) + storage_size(self%loc) + self%pairs%storage()

        print *, 'storage_size = ', storage
        ! do j = 1, size(self%grids, 2)
        !     do i = 1, size(self%grids, 1)
        !         storage = storage + self%grids(i, j)%storage()
        !     end do
        ! end do
        ! print *, 'storage_size = ', storage
        do i = 0, size(self%threads_pairs) - 1
            storage = storage + self%threads_pairs(i)%storage()
        end do
        print *, 'storage_size = ', storage

    end function storage

end module nnps_grid2d_module
