!> 二维背景网格搜索
!> 2D background grid method
module nnps_grid2d_module

    use nnps_kinds, only: wp
    use nnps_vector, only: vector, vector_finalizer
    use nnps_int_vector, only: int_vector, int_vector_finalizer
    use nnps_spatial_hashing, only: shash_tbl, shash_tbl_finalizer
    use nnps_math, only: distance2d, sqrt_eps
    use, intrinsic :: iso_fortran_env, only: error_unit
#ifndef SERIAL
    use omp_lib, only: omp_get_thread_num, omp_get_max_threads
#endif
    implicit none

    private
    public :: nnps_grid2d

    !> 二维背景网格搜索
    !> 2d grid
    !> @note pairs(threads_pairs)/shash_tbl are the No.1/No.2 memory consumers
    type nnps_grid2d
        real(wp), pointer :: loc(:, :)                              !! particle 2d coordinate
        type(shash_tbl) :: tbl                                      !! background grids hash table
        type(int_vector), private :: iks                            !! unique keys

#ifdef SERIAL
        type(vector), private :: pairs                              !! pairs
        type(int_vector), private :: idxs                           !! indexes
#else
        type(vector), allocatable, private :: threads_pairs(:)      !! thread local pairs
        type(int_vector), allocatable, private :: threads_idxs(:)   !! thread local indexes
#endif
    contains
        procedure :: init, query, recycle, destroy
    end type nnps_grid2d

contains

    !> 初始化背景网格
    !> initialize
    subroutine init(self, loc, n, load_factor)
        class(nnps_grid2d), intent(inout) :: self               !! nnps_grid2d
        real(wp), dimension(:, :), intent(in), target :: loc    !! particle 2d coordinate
        integer, intent(in) :: n                                !! number of particles
        real(wp), intent(in), optional :: load_factor           !! 负载因子, 默认值为 1.0
        integer :: thread_num

        self%loc => loc
#ifdef SERIAL
        thread_num = 0
        call self%pairs%init(2, n)
#else
        thread_num = omp_get_max_threads() - 1
        allocate (self%threads_pairs(0:thread_num), &           ! allocate by threads number
                  self%threads_idxs(0:thread_num))
        call self%threads_pairs(:)%init(2, n)
#endif

        if (present(load_factor)) then
            call self%tbl%alloc(m=int(n/load_factor))
        else
            call self%tbl%alloc(m=n)
        end if

    end subroutine init

    !> 构建背景网格与查询粒子对
    !> build and query
    subroutine query(self, radius, pairs, rdxs, n)
        class(nnps_grid2d), intent(inout), target :: self       !! nnps_grid2d
        integer, intent(in) :: n                                !! number of particles
        real(wp), intent(in) :: radius                          !! grid length, smoothing length
        integer, dimension(:), pointer, intent(out) :: pairs    !! particle pairs
        real(wp), pointer, dimension(:), intent(out) :: rdxs    !! particle pairs distance
        integer :: i, j, ik(3), ijk(3, 4), istat
        integer, pointer :: values(:)

#ifndef SERIAL
        integer :: thread_id
#endif

        call self%tbl%zeroing()
        ik(3) = 1
        self%iks%len = 0

        do i = 1, n                                             ! 一次查询，将所有粒子的位置信息存入哈希表
            ik(1:2) = ceiling(self%loc(:, i)/radius)
            call self%tbl%set(key=ik, value=i, istat=istat)
            if (istat == 0) call self%iks%push_back_items(ik, 3)! collect unique keys
        end do

#ifndef SERIAL

        self%threads_pairs%len = 0
        associate (grid => self%tbl%buckets, iks => self%iks%items)

            !$omp parallel do private(i, ijk, values, thread_id) schedule(dynamic)
            do i = 1, self%iks%len, 3

                ijk(:, 1) = [iks(i:i + 1) - 1, iks(i + 2)]           !          ___________
                ijk(:, 2) = [iks(i), iks(i + 1) - 1, iks(i + 2)]     !          |         |
                ijk(:, 3) = [iks(i) + 1, iks(i + 1) - 1, iks(i + 2)] !  - - -   | |     | |
                ijk(:, 4) = [iks(i) - 1, iks(i + 1:i + 2)]           !  x o -   | |     | |
                !                                                       x x x   | |_____| |
                ! U style, L style, 4 neighbors                      !          |_________|

                thread_id = omp_get_thread_num()
                self%threads_idxs(thread_id)%len = 0

                nullify (values)
                do j = 1, 4
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

        associate (len => self%threads_pairs(0)%len)
            pairs => self%threads_pairs(0)%items(1:len*2)
            rdxs => self%threads_pairs(0)%ritems(1:len*3)
        end associate
#else

        self%pairs%len = 0
        associate (iks => self%iks%items)
            do i = 1, self%iks%len, 3
                ijk(:, 1) = [iks(i:i + 1) - 1, iks(i + 2)]           !          ___________
                ijk(:, 2) = [iks(i), iks(i + 1) - 1, iks(i + 2)]     !          |         |
                ijk(:, 3) = [iks(i) + 1, iks(i + 1) - 1, iks(i + 2)] !  - - -   | |     | |
                ijk(:, 4) = [iks(i) - 1, iks(i + 1:i + 2)]           !  x o -   | |     | |
                !                                                       x x x   | |_____| |
                ! U style, L style, 4 neighbors                      !          |_________|

                self%idxs%len = 0
                nullify (values)
                do j = 1, 4
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

        associate (len => self%pairs%len)
            pairs => self%pairs%items(1:len*2)
            rdxs => self%pairs%ritems(1:len*3)
        end associate

#endif

    contains

        pure subroutine adjacent_grid_neighbors(main, found, pairs)
            integer, intent(in) :: main(:)
            integer, intent(in) :: found(:)
            type(vector), intent(inout) :: pairs
            integer :: ii, jj
            real(wp) :: rdx(3)

            do ii = 1, size(main)

                do jj = ii + 1, size(main)
                    call distance2d(self%loc(:, main(ii)), &
                                    self%loc(:, main(jj)), rdx(1), rdx(2:3))
                    if (rdx(1) < radius) call pairs%push([main(ii), main(jj)], rdx)
                end do

                do jj = 1, size(found)
                    call distance2d(self%loc(:, main(ii)), &
                                    self%loc(:, found(jj)), rdx(1), rdx(2:3))
                    if (rdx(1) < radius) call pairs%push([main(ii), found(jj)], rdx)
                end do

            end do

        end subroutine adjacent_grid_neighbors

        pure subroutine self_grid_neighbors(main, pairs)
            integer, intent(in) :: main(:)
            type(vector), intent(inout) :: pairs
            integer :: ii, jj
            real(wp) :: rdx(3)

            do ii = 1, size(main)

                do jj = ii + 1, size(main)
                    call distance2d(self%loc(:, main(ii)), &
                                    self%loc(:, main(jj)), rdx(1), rdx(2:3))
                    if (rdx(1) < radius) call pairs%push([main(ii), main(jj)], rdx)
                end do

            end do

        end subroutine self_grid_neighbors

    end subroutine query

    !> 回收惰性空间: 标记
    pure subroutine recycle(self)
        class(nnps_grid2d), intent(inout) :: self

        call self%tbl%recycle()

    end subroutine recycle

    !> finalizer
    pure subroutine destroy(self)
        class(nnps_grid2d), intent(inout) :: self  !! nnps_grid2d

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

    end subroutine destroy

end module nnps_grid2d_module
