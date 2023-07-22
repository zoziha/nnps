!> Compact Hashing
module nnps_compact_hashing

    use nnps_int_vector, only: int_vector
    implicit none

    private
    public :: chash_tbl

    !> Compact Hashing
    type chash_tbl
        type(int_vector), allocatable :: buckets(:)  !! Buckets
    contains
        procedure :: allocate => chash_tbl_allocate, clear, zeroing, push, hash
    end type chash_tbl

contains

    !> Hashing function
    pure integer function hash(self, i, j, k)
        class(chash_tbl), intent(in) :: self
        integer, intent(in) :: i, j, k

        hash = modulo(ieor(ieor(73856093*i, 19349663*j), 83492791*k), size(self%buckets))

    end function hash

    !> Allocate
    pure subroutine chash_tbl_allocate(self, m)
        class(chash_tbl), intent(inout) :: self
        integer, intent(in) :: m

        allocate(self%buckets(0:m-1))

    end subroutine chash_tbl_allocate

    !> Clean up
    pure subroutine clear(self)
        class(chash_tbl), intent(inout) :: self

        if (allocated(self%buckets)) then
            call self%buckets(:)%clear()
            deallocate(self%buckets)
        end if

    end subroutine clear

    !> Zeroing
    pure subroutine zeroing(self)
        class(chash_tbl), intent(inout) :: self

        where (self%buckets(:)%len /= 0) self%buckets(:)%len = 0

    end subroutine zeroing

    !> Push
    pure subroutine push(self, i, j, k, index)
        class(chash_tbl), intent(inout) :: self
        integer, intent(in) :: i, j, k, index

        call self%buckets(self%hash(i, j, k))%push(index)

    end subroutine push

end module nnps_compact_hashing
