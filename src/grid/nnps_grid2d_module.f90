!> 2D background grid method
module nnps_grid2d_module

    use nnps_kinds, only: rk
    use nnps_vector, only: vector
    use nnps_int_vector, only: int_vector
    use nnps_spatial_hashing, only: chash_tbl
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
        type(int_vector) :: iks  !! unique keys
        real(rk), private :: radius
    contains
        procedure :: init, query, storage
    end type nnps_grid2d

contains

    !> initialize: U style
    subroutine init(self, loc, radius, n)
        class(nnps_grid2d), intent(inout) :: self  !! nnps_grid2d
        real(rk), dimension(:, :), intent(in), target :: loc  !! particle 2d coordinate
        real(rk), intent(in) :: radius  !! grid length, smoothing length
        integer, intent(in) :: n  !! number of particles

        self%loc => loc
        allocate (self%threads_pairs(0:omp_get_max_threads() - 1))
        call self%pairs%init(2, 8*n)
        call self%threads_pairs(:)%init(2, 2*n)
        self%radius = radius

        call self%tbl%allocate(m=2*n)

    end subroutine init

    !> query
    subroutine query(self, radius, pairs, rdxs, n)
        class(nnps_grid2d), intent(inout), target :: self  !! nnps_grid2d
        integer, intent(in) :: n  !! number of particles
        real(rk), intent(in) :: radius  !! grid length, smoothing length
        integer, dimension(:), pointer, intent(out) :: pairs  !! particle pairs
        real(rk), pointer, dimension(:), intent(out) :: rdxs  !! particle pairs distance
        integer :: i, j, idx(5), ik(3), ijk(3, 5)
        integer, allocatable :: idxs(:)
        integer, pointer :: values(:)
        logical :: lstat

        call self%tbl%zeroing()
        self%min = minval(self%loc, 2) - sqrt_eps

        self%iks%len = 0
        ik(3) = 1
        do i = 1, n  ! cannot use parallel do here
            ik(1:2) = ceiling((self%loc(:, i) - self%min)/radius)
            call self%tbl%set(key=ik, value=i, stat=lstat)
            if (lstat) call self%iks%push_back3(ik)  ! collect unique keys
        end do

        self%threads_pairs%len = 0
        associate (grid => self%tbl%buckets, iks => self%iks%items)

            !$omp parallel do private(i, idx, idxs, ijk, values) schedule(dynamic)
            do i = 1, self%iks%len, 3  ! @todo improve performance, not malloc and realloc frequently

                ijk(:, 1) = [iks(i) - 1, iks(i + 1) - 1, iks(i + 2)]
                ijk(:, 2) = [iks(i), iks(i + 1) - 1, iks(i + 2)]
                ijk(:, 3) = [iks(i) + 1, iks(i + 1) - 1, iks(i + 2)]
                ijk(:, 4) = [iks(i) - 1, iks(i + 1), iks(i + 2)]
                ijk(:, 5) = [iks(i:i + 2)]

                ! U style, L style, 4 neighbors         !          ___________
                idx = [self%tbl%hash(ijk(:, 1)), &      !          |         |
                       self%tbl%hash(ijk(:, 2)), &      !  - - -   | |     | |
                       self%tbl%hash(ijk(:, 3)), &      !  x o -   | |     | |
                       self%tbl%hash(ijk(:, 4)), &      !  x x x   | |_____| |
                       self%tbl%hash(ijk(:, 5))]        !          |_________|

                allocate (idxs(0))

                do j = 1, 4
                    nullify (values)
                    call grid(idx(j))%get_value(ijk(:, j), values)
                    if (associated(values)) idxs = [idxs, values]
                end do

                call grid(idx(5))%get_value(ijk(:, 5), values)
                call find_nearby_particles(values, idxs, &
                                           self%threads_pairs(omp_get_thread_num()))

                deallocate (idxs)
                nullify (values)

            end do

        end associate

        call self%pairs%merge(self%threads_pairs)

        pairs => self%pairs%items(1:self%pairs%len*2)
        rdxs => self%pairs%ritems(1:self%pairs%len*3)

    contains

        pure subroutine find_nearby_particles(main, found, threads_pairs)
            integer, intent(in) :: main(:)
            integer, intent(in) :: found(:)
            type(vector), intent(inout) :: threads_pairs
            integer :: ii, jj
            real(rk) :: rdx(3)

            do ii = 1, size(main)

                do jj = ii + 1, size(main)
                    call distance2d(self%loc(:, main(ii)), &
                                    self%loc(:, main(jj)), rdx(1), rdx(2:3))
                    if (rdx(1) < radius) call threads_pairs%push([main(ii), main(jj)], rdx)
                end do

                do jj = 1, size(found)
                    call distance2d(self%loc(:, main(ii)), &
                                    self%loc(:, found(jj)), rdx(1), rdx(2:3))
                    if (rdx(1) < radius) call threads_pairs%push([main(ii), found(jj)], rdx)
                end do

            end do

        end subroutine find_nearby_particles

    end subroutine query

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
