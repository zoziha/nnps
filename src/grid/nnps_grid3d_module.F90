!> 3D background grid method
module nnps_grid3d_module

    use nnps_kinds, only: rk
    use nnps_vector, only: vector, vector_finalizer
    use nnps_int_vector, only: int_vector, int_vector_finalizer
    use nnps_spatial_hashing, only: shash_tbl, shash_tbl_finalizer
    use nnps_math, only: distance3d, sqrt_eps
#ifndef SERIAL
    use omp_lib, only: omp_get_thread_num, omp_get_max_threads
#endif
    implicit none

    private
    public :: nnps_grid3d, nnps_grid3d_finalizer

    !> 3d grid
    !> @note pairs(threads_pairs)/shash_tbl are the No.1/No.2 memory consumers
    type nnps_grid3d
        real(rk), pointer :: loc(:, :)                              !! particle 3d coordinate
        type(shash_tbl) :: tbl                                      !! background grids hash table
        type(int_vector), private :: iks                            !! unique keys

#ifndef SERIAL
        type(int_vector), private :: remains(2)                     !! remains indexes
        type(vector), allocatable, private :: threads_pairs(:)      !! thread local pairs
        type(int_vector), allocatable, private :: threads_idxs(:)   !! thread local indexes
#else
        type(vector), private :: pairs      !! pairs
        type(int_vector), private :: idxs   !! indexes
#endif
    contains
        procedure :: init, query, storage
    end type nnps_grid3d

contains

    !> finalizer
    pure subroutine nnps_grid3d_finalizer(self)
        type(nnps_grid3d), intent(inout) :: self  !! nnps_grid3d

        if (associated(self%loc)) nullify (self%loc)
        call shash_tbl_finalizer(self%tbl)
        call int_vector_finalizer(self%iks)

#ifndef SERIAL
        if (allocated(self%threads_pairs)) then
            call vector_finalizer(self%threads_pairs)
            deallocate (self%threads_pairs)
        end if
        if (allocated(self%threads_idxs)) then
            call int_vector_finalizer(self%threads_idxs)
            deallocate (self%threads_idxs)
        end if
#else
        call vector_finalizer(self%pairs)
        call int_vector_finalizer(self%idxs)
#endif

    end subroutine nnps_grid3d_finalizer

    !> initialize
    subroutine init(self, loc, n)
        class(nnps_grid3d), intent(inout) :: self               !! nnps_grid3d
        real(rk), dimension(:, :), intent(in), target :: loc    !! particle 3d coordinate
        integer, intent(in) :: n                                !! number of particles

        self%loc => loc

#ifndef SERIAL
        allocate (self%threads_pairs(0:omp_get_max_threads() - 1), &  ! allocate by threads number
                  self%threads_idxs(0:omp_get_max_threads() - 1))
        call self%threads_pairs(:)%init(3, n)
#else
        call self%pairs%init(3, n)
#endif

        call self%tbl%allocate(m=n)

    end subroutine init

    !> query
    subroutine query(self, radius, pairs, rdxs, n)
        class(nnps_grid3d), intent(inout), target :: self       !! nnps_grid3d
        integer, intent(in) :: n                                !! number of pairs
        real(rk), intent(in) :: radius                          !! search radius
        integer, dimension(:), pointer, intent(out) :: pairs    !! particle pairs
        real(rk), dimension(:), pointer, intent(out) :: rdxs    !! particle pairs distance
        integer :: i, j, ik(3), ijk(3, 13), istat
        integer, pointer :: values(:)

#ifndef SERIAL
        integer :: thread_id
#endif

        call self%tbl%zeroing()
        self%iks%len = 0

        do i = 1, n  !@todo: parallelize
            ik = ceiling(self%loc(:, i)/radius)
            call self%tbl%set(key=ik, value=i, istat=istat)
            if (istat == 0) call self%iks%push_back_items(ik, 3)  ! collect unique keys
        end do

#ifndef SERIAL

        self%threads_pairs%len = 0
        associate (grid => self%tbl%buckets, iks => self%iks%items)

            !$omp parallel do private(i, ijk, values, thread_id)
            do i = 1, self%iks%len, 3

                ijk(:, 1) = iks(i:i + 2) - 1  ! 3D L style, 13 neighbors (9 + 4)
                ijk(:, 2) = [iks(i), iks(i + 1:i + 2) - 1]
                ijk(:, 3) = [iks(i) + 1, iks(i + 1:i + 2) - 1]
                ijk(:, 4) = [iks(i) - 1, iks(i + 1), iks(i + 2) - 1]
                ijk(:, 5) = [iks(i:i + 1), iks(i + 2) - 1]
                ijk(:, 6) = [iks(i) + 1, iks(i + 1), iks(i + 2) - 1]
                ijk(:, 7) = [iks(i) - 1, iks(i + 1) + 1, iks(i + 2) - 1]
                ijk(:, 8) = [iks(i), iks(i + 1) + 1, iks(i + 2) - 1]
                ijk(:, 9) = [iks(i:i + 1) + 1, iks(i + 2) - 1]
                ijk(:, 10) = [iks(i:i + 1) - 1, iks(i + 2)]
                ijk(:, 11) = [iks(i), iks(i + 1) - 1, iks(i + 2)]
                ijk(:, 12) = [iks(i) + 1, iks(i + 1) - 1, iks(i + 2)]
                ijk(:, 13) = [iks(i) - 1, iks(i + 1:i + 2)]

                thread_id = omp_get_thread_num()
                self%threads_idxs(thread_id)%len = 0

                nullify (values)
                do j = 1, 13
                    call grid(self%tbl%hash(ijk(:, j)))%get_value(ijk(:, j), values)
                    if (associated(values)) then
                        call self%threads_idxs(thread_id)%push_back_items(values, size(values))
                        nullify (values)
                    end if
                end do

                call grid(self%tbl%hash(iks(i:i + 2)))%get_value(iks(i:i + 2), values)
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

        if (size(self%threads_pairs) > 1) call self%threads_pairs(0)%merge(self%threads_pairs)

        pairs => self%threads_pairs(0)%items(1:self%threads_pairs(0)%len*2)
        rdxs => self%threads_pairs(0)%ritems(1:self%threads_pairs(0)%len*4)

#else

        self%pairs%len = 0
        associate (iks => self%iks%items)
            do i = 1, self%iks%len, 3
                ijk(:, 1) = iks(i:i + 2) - 1  ! 3D L style, 13 neighbors (9 + 4)
                ijk(:, 2) = [iks(i), iks(i + 1:i + 2) - 1]
                ijk(:, 3) = [iks(i) + 1, iks(i + 1:i + 2) - 1]
                ijk(:, 4) = [iks(i) - 1, iks(i + 1), iks(i + 2) - 1]
                ijk(:, 5) = [iks(i:i + 1), iks(i + 2) - 1]
                ijk(:, 6) = [iks(i) + 1, iks(i + 1), iks(i + 2) - 1]
                ijk(:, 7) = [iks(i) - 1, iks(i + 1) + 1, iks(i + 2) - 1]
                ijk(:, 8) = [iks(i), iks(i + 1) + 1, iks(i + 2) - 1]
                ijk(:, 9) = [iks(i:i + 1) + 1, iks(i + 2) - 1]
                ijk(:, 10) = [iks(i:i + 1) - 1, iks(i + 2)]
                ijk(:, 11) = [iks(i), iks(i + 1) - 1, iks(i + 2)]
                ijk(:, 12) = [iks(i) + 1, iks(i + 1) - 1, iks(i + 2)]
                ijk(:, 13) = [iks(i) - 1, iks(i + 1:i + 2)]

                nullify (values)
                self%idxs%len = 0
                do j = 1, 13
                    call self%tbl%buckets(self%tbl%hash(ijk(:, j)))%get_value(ijk(:, j), values)
                    if (associated(values)) then
                        call self%idxs%push_back_items(values, size(values))
                        nullify (values)
                    end if
                end do

                call self%tbl%buckets(self%tbl%hash(iks(i:i + 2)))%get_value(iks(i:i + 2), values)
                if (self%idxs%len == 0) then
                    if (size(values) > 1) call self_grid_neighbors(values, self%pairs)
                else
                    call adjacent_grid_neighbors(values, self%idxs%items(1:self%idxs%len), self%pairs)
                end if
                nullify (values)

            end do
        end associate

        pairs => self%pairs%items(1:self%pairs%len*2)
        rdxs => self%pairs%ritems(1:self%pairs%len*4)

#endif

    contains

        pure subroutine adjacent_grid_neighbors(main, found, pairs)
            integer, intent(in) :: main(:)
            integer, intent(in) :: found(:)
            type(vector), intent(inout) :: pairs
            integer :: ii, jj
            real(rk) :: rdx(4)

            do ii = 1, size(main)

                do jj = ii + 1, size(main)
                    call distance3d(self%loc(:, main(ii)), &
                                    self%loc(:, main(jj)), rdx(1), rdx(2:4))
                    if (rdx(1) < radius) call pairs%push([main(ii), main(jj)], rdx)
                end do

                do jj = 1, size(found)
                    call distance3d(self%loc(:, main(ii)), &
                                    self%loc(:, found(jj)), rdx(1), rdx(2:4))
                    if (rdx(1) < radius) call pairs%push([main(ii), found(jj)], rdx)
                end do

            end do

        end subroutine adjacent_grid_neighbors

        pure subroutine self_grid_neighbors(main, pairs)
            integer, intent(in) :: main(:)
            type(vector), intent(inout) :: pairs
            integer :: ii, jj
            real(rk) :: rdx(4)

            do ii = 1, size(main)

                do jj = ii + 1, size(main)
                    call distance3d(self%loc(:, main(ii)), &
                                    self%loc(:, main(jj)), rdx(1), rdx(2:4))
                    if (rdx(1) < radius) call pairs%push([main(ii), main(jj)], rdx)
                end do

            end do

        end subroutine self_grid_neighbors

    end subroutine query

    !> storage in bits (shash_tbl/all)
    function storage(self)
        class(nnps_grid3d), intent(in) :: self  !! nnps_grid3d
        integer, dimension(2) :: storage
        integer :: i

        storage(1) = self%tbl%storage()
#ifndef SERIAL
        storage(2) = storage_size(self) + storage(1) + storage_size(self%loc) + &
                     self%iks%storage()

        do i = 0, size(self%threads_pairs) - 1
            storage(2) = storage(2) + self%threads_pairs(i)%storage() + &
                         self%threads_idxs(i)%storage() + &
                         storage_size(self%threads_pairs(i)) + &
                         storage_size(self%threads_idxs(i))
        end do
#else
        storage(2) = storage_size(self) + storage(1) + storage_size(self%loc) + &
                     self%iks%storage() + storage_size(self%pairs) + &
                     self%idxs%storage()
#endif

    end function storage

end module nnps_grid3d_module
