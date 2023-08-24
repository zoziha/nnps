!> 2D background grid method
module nnps_grid2d_module

    use nnps_kinds, only: rk
    use nnps_vector, only: vector, vector_finalizer
    use nnps_int_vector, only: int_vector, int_vector_finalizer
    use nnps_spatial_hashing, only: shash_tbl, shash_tbl_finalizer
    use nnps_math, only: distance2d, sqrt_eps
#ifdef OMP
    use omp_lib, only: omp_get_thread_num, omp_get_max_threads
#endif
    implicit none

    private
    public :: nnps_grid2d, nnps_grid2d_finalizer

    !> 2d grid
    !> @note shash_tbl/pairs(threads_pairs) are the No.1/No.2 memory consumers
    type nnps_grid2d
        real(rk), pointer :: loc(:, :)  !! particle 2d coordinate
        type(shash_tbl) :: tbl  !! background grids hash table
        type(int_vector), private :: iks  !! unique keys
        type(vector), allocatable, private :: threads_pairs(:)  !! thread local pairs
        type(int_vector), allocatable, private :: threads_idxs(:)  !! thread local indexes
    contains
        procedure :: init, query, storage
    end type nnps_grid2d

contains

    !> finalizer
    pure subroutine nnps_grid2d_finalizer(self)
        type(nnps_grid2d), intent(inout) :: self  !! nnps_grid2d

        if (associated(self%loc)) nullify (self%loc)
        call shash_tbl_finalizer(self%tbl)
        call int_vector_finalizer(self%iks)

        if (allocated(self%threads_pairs)) then
            call vector_finalizer(self%threads_pairs)
            deallocate (self%threads_pairs)
        end if
        if (allocated(self%threads_idxs)) then
            call int_vector_finalizer(self%threads_idxs)
            deallocate (self%threads_idxs)
        end if

    end subroutine nnps_grid2d_finalizer

    !> initialize
    subroutine init(self, loc, n)
        class(nnps_grid2d), intent(inout) :: self  !! nnps_grid2d
        real(rk), dimension(:, :), intent(in), target :: loc  !! particle 2d coordinate
        integer, intent(in) :: n  !! number of particles

        self%loc => loc

#ifdef OMP
        allocate (self%threads_pairs(0:omp_get_max_threads() - 1), &
                  self%threads_idxs(0:omp_get_max_threads() - 1))
#else
        allocate (self%threads_pairs(0:0), self%threads_idxs(0:0))
#endif
        call self%threads_pairs(:)%init(2, n)

        call self%tbl%allocate(m=2*n)

    end subroutine init

    !> query
    subroutine query(self, radius, pairs, rdxs, n)
        class(nnps_grid2d), intent(inout), target :: self  !! nnps_grid2d
        integer, intent(in) :: n  !! number of particles
        real(rk), intent(in) :: radius  !! grid length, smoothing length
        integer, dimension(:), pointer, intent(out) :: pairs  !! particle pairs
        real(rk), pointer, dimension(:), intent(out) :: rdxs  !! particle pairs distance
        integer :: i, j, idx(5), ik(3), ijk(3, 5), thread_id
        integer, pointer :: values(:)
        real(rk) :: min(2)
        logical :: lstat

        call self%tbl%zeroing()
        min = minval(self%loc, 2) - sqrt_eps
        thread_id = 0

        self%iks%len = 0
        ik(3) = 1
        do i = 1, n  !@todo use parallel and sorting here
            ik(1:2) = ceiling((self%loc(:, i) - min)/radius)
            call self%tbl%set(key=ik, value=i, stat=lstat)
            if (lstat) call self%iks%push_back_items(ik, 3)  ! collect unique keys
        end do

        self%threads_pairs%len = 0
        associate (grid => self%tbl%buckets, iks => self%iks%items)

#ifdef OMP
            !$omp parallel do private(i, idx, ijk, values, thread_id) schedule(dynamic)
#endif
            do i = 1, self%iks%len, 3

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

#ifdef OMP
                thread_id = omp_get_thread_num()
#endif
                self%threads_idxs(thread_id)%len = 0

                nullify (values)
                do j = 1, 4
                    call grid(idx(j))%get_value(ijk(:, j), values)
                    if (associated(values)) then
                        call self%threads_idxs(thread_id)%push_back_items(values, size(values))
                        nullify (values)
                    end if
                end do

                call grid(idx(5))%get_value(ijk(:, 5), values)
                if (self%threads_idxs(thread_id)%len == 0) then
                    if (size(values) > 1) call self_grid_neighbors(values, &
                        &self%threads_pairs(thread_id))
                else
                    call adjacent_grid_neighbors(values, &
                    &self%threads_idxs(thread_id)%items(1:self%threads_idxs(thread_id)%len), &
                    &self%threads_pairs(thread_id))
                end if

                nullify (values)

            end do

        end associate

#ifdef OMP
        call self%threads_pairs(0)%merge(self%threads_pairs)
#endif
        pairs => self%threads_pairs(0)%items(1:self%threads_pairs(0)%len*2)
        rdxs => self%threads_pairs(0)%ritems(1:self%threads_pairs(0)%len*3)

    contains

        pure subroutine adjacent_grid_neighbors(main, found, threads_pairs)
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

        end subroutine adjacent_grid_neighbors

        pure subroutine self_grid_neighbors(main, threads_pairs)
            integer, intent(in) :: main(:)
            type(vector), intent(inout) :: threads_pairs
            integer :: ii, jj
            real(rk) :: rdx(3)

            do ii = 1, size(main)

                do jj = ii + 1, size(main)
                    call distance2d(self%loc(:, main(ii)), &
                                    self%loc(:, main(jj)), rdx(1), rdx(2:3))
                    if (rdx(1) < radius) call threads_pairs%push([main(ii), main(jj)], rdx)
                end do

            end do

        end subroutine self_grid_neighbors

    end subroutine query

    !> storage in bits (shash_tbl/all)
    function storage(self)
        class(nnps_grid2d), intent(in) :: self  !! nnps_grid2d
        integer, dimension(2) :: storage
        integer :: i

        storage(1) = self%tbl%storage()
        storage(2) = storage_size(self) + storage(1) + storage_size(self%loc) + &
                     self%iks%storage()

        do i = 0, size(self%threads_pairs) - 1
            storage(2) = storage(2) + self%threads_pairs(i)%storage() + &
                         self%threads_idxs(i)%storage() + &
                         storage_size(self%threads_pairs(i)) + &
                         storage_size(self%threads_idxs(i))
        end do

    end function storage

end module nnps_grid2d_module
