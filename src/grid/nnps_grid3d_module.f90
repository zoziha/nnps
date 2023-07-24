!> 3D background grid method
module nnps_grid3d_module

    use nnps_kinds, only: rk
    use nnps_vector, only: vector
    use nnps_int_vector, only: int_vector
    use nnps_spatial_hashing, only: shash_tbl
    use nnps_math, only: distance3d, sqrt_eps
    use omp_lib, only: omp_get_thread_num, omp_get_max_threads
    implicit none

    private
    public :: nnps_grid3d

    !> 3d grid
    type nnps_grid3d
        real(rk), pointer :: loc(:, :)  !! particle 3d coordinate
        type(shash_tbl) :: tbl  !! background grid hash table
        type(vector) :: pairs  !! particle pairs
        real(rk), dimension(3), private :: min
        type(int_vector), private :: iks  !! unique keys
        type(vector), allocatable, private :: threads_pairs(:)  !! thread local pairs
        type(int_vector), allocatable, private :: threads_idxs(:)  !! thread local indexes
    contains
        procedure :: init, query, storage
    end type nnps_grid3d

contains

    !> initialize
    subroutine init(self, loc, n)
        class(nnps_grid3d), intent(inout) :: self  !! nnps_grid3d
        real(rk), dimension(:, :), intent(in), target :: loc  !! particle 3d coordinate
        integer, intent(in) :: n  !! number of particles

        self%loc => loc
        allocate (self%threads_pairs(0:omp_get_max_threads() - 1), &
                  self%threads_idxs(0:omp_get_max_threads() - 1))
        call self%pairs%init(3, 12*n)
        call self%threads_pairs(:)%init(3, n)

        call self%tbl%allocate(m=2*n)

    end subroutine init

    !> query
    subroutine query(self, radius, pairs, rdxs, n)
        class(nnps_grid3d), intent(inout), target :: self  !! nnps_grid3d
        integer, intent(in) :: n  !! number of pairs
        real(rk), intent(in) :: radius  !! search radius
        integer, dimension(:), pointer, intent(out) :: pairs  !! particle pairs
        real(rk), dimension(:), pointer, intent(out) :: rdxs  !! particle pairs distance
        integer :: i, j, idx(14), ik(3), ijk(3, 14), thread_id
        integer, pointer :: values(:)
        logical :: lstat

        call self%tbl%zeroing()
        self%min = minval(self%loc, dim=2) - sqrt_eps

        self%iks%len = 0
        do i = 1, n
            ik = ceiling((self%loc(:, i) - self%min)/radius)
            call self%tbl%set(key=ik, value=i, stat=lstat)
            if (lstat) call self%iks%push_back_items(ik, 3)  ! collect unique keys
        end do

        self%threads_pairs%len = 0
        associate (grid => self%tbl%buckets, iks => self%iks%items)

            !$omp parallel do private(i, idx, ijk, values, thread_id) schedule(dynamic)
            do i = 1, self%iks%len, 3

                ijk(:, 1) = [iks(i:i + 2) - 1]  ! 3D L style, 13 neighbors (9 + 4)
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
                ijk(:, 14) = [iks(i:i + 2)]

                do j = 1, 14
                    idx(j) = self%tbl%hash(ijk(:, j))
                end do

                thread_id = omp_get_thread_num()
                self%threads_idxs(thread_id)%len = 0

                nullify (values)
                do j = 1, 13
                    call grid(idx(j))%get_value(ijk(:, j), values)
                    if (associated(values)) then
                        call self%threads_idxs(thread_id)%push_back_items(values, size(values))
                        nullify (values)
                    end if
                end do

                call grid(idx(14))%get_value(ijk(:, 14), values)
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

        call self%pairs%merge(self%threads_pairs)

        pairs => self%pairs%items(1:self%pairs%len*2)
        rdxs => self%pairs%ritems(1:self%pairs%len*4)

    contains

        pure subroutine adjacent_grid_neighbors(main, found, threads_pairs)
            integer, intent(in) :: main(:)
            integer, intent(in) :: found(:)
            type(vector), intent(inout) :: threads_pairs
            integer :: ii, jj
            real(rk) :: rdx(4)

            do ii = 1, size(main)

                do jj = ii + 1, size(main)
                    call distance3d(self%loc(:, main(ii)), &
                                    self%loc(:, main(jj)), rdx(1), rdx(2:4))
                    if (rdx(1) < radius) call threads_pairs%push([main(ii), main(jj)], rdx)
                end do

                do jj = 1, size(found)
                    call distance3d(self%loc(:, main(ii)), &
                                    self%loc(:, found(jj)), rdx(1), rdx(2:4))
                    if (rdx(1) < radius) call threads_pairs%push([main(ii), found(jj)], rdx)
                end do

            end do

        end subroutine adjacent_grid_neighbors

        pure subroutine self_grid_neighbors(main, threads_pairs)
            integer, intent(in) :: main(:)
            type(vector), intent(inout) :: threads_pairs
            integer :: ii, jj
            real(rk) :: rdx(4)

            do ii = 1, size(main)

                do jj = ii + 1, size(main)
                    call distance3d(self%loc(:, main(ii)), &
                                    self%loc(:, main(jj)), rdx(1), rdx(2:4))
                    if (rdx(1) < radius) call threads_pairs%push([main(ii), main(jj)], rdx)
                end do

            end do

        end subroutine self_grid_neighbors

    end subroutine query

    !> storage @todo to complete
    pure integer function storage(self)
        class(nnps_grid3d), intent(in) :: self
        integer :: i, j, k

        storage = storage_size(self) + storage_size(self%loc) + self%pairs%storage()
        ! do k = 1, size(self%grids, 3)
        !     do j = 1, size(self%grids, 2)
        !         do i = 1, size(self%grids, 1)
        !             storage = storage + self%grids(i, j, k)%storage()
        !         end do
        !     end do
        ! end do
        do i = 1, size(self%threads_pairs)
            storage = storage + self%threads_pairs(i)%storage()
        end do

    end function storage

end module nnps_grid3d_module
