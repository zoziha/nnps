!> 2D direct neighbor search
module nnps_direct2d_module

    use nnps_kinds, only: rk
    use nnps_vector, only: vector
    use nnps_math, only: distance2d
    use omp_lib, only: omp_get_thread_num, omp_get_max_threads
    implicit none

    private
    public :: nnps_direct2d

    !> nnps_direct2d
    type nnps_direct2d
        real(rk), pointer :: loc(:, :)  !! particle 2d coordinate
        type(vector), allocatable, private :: threads_pairs(:)  !! thread local pairs
        type(vector) :: pairs  !! partcile pairs
    contains
        procedure :: init, query
    end type nnps_direct2d

contains

    !> initialize
    subroutine init(self, loc, cap)
        class(nnps_direct2d), intent(inout), target :: self
        real(rk), intent(in), target :: loc(:, :)
        integer, intent(in), optional :: cap

        self%loc => loc
        allocate (self%threads_pairs(0:omp_get_max_threads() - 1))
        call self%pairs%init(2, cap)
        call self%threads_pairs(:)%init(2, cap)

    end subroutine init

    !> query
    subroutine query(self, radius, pairs, rdxs)
        class(nnps_direct2d), intent(inout), target :: self
        real(rk), intent(in) :: radius  !! query radius
        integer, pointer :: pairs(:)  !! particle pairs
        real(rk), dimension(:), pointer :: rdxs
        integer :: i, j
        real(rk) :: rdx(3)

        self%threads_pairs%len = 0

        !$omp parallel do private(i, j, rdx) schedule(dynamic)
        do i = 1, size(self%loc, 2) - 1
            do j = i + 1, size(self%loc, 2)

                call distance2d(self%loc(:, i), self%loc(:, j), rdx(1), rdx(2:3))
                if (rdx(1) < radius) call self%threads_pairs(omp_get_thread_num())%push([i, j], rdx)

            end do
        end do

        call self%pairs%merge(self%threads_pairs)

        pairs => self%pairs%items(1:self%pairs%len*2)
        rdxs => self%pairs%ritems(1:self%pairs%len*3)

    end subroutine query

end module nnps_direct2d_module
