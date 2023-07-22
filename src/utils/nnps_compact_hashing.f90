!> Compact Hashing
module nnps_compact_hashing

    use nnps_key_value, only: key_value
    implicit none

    private
    public :: chash_tbl

    !> Compact Hashing
    type chash_tbl
        type(key_value), allocatable :: buckets(:)  !! Buckets
    contains
        procedure :: allocate => chash_tbl_allocate, clear, zeroing, set, hash
    end type chash_tbl

contains

    !> Hashing function
    pure integer function hash(self, key)
        class(chash_tbl), intent(in) :: self
        integer, intent(in) :: key(3)

        hash = modulo(ieor(ieor(73856093*key(1), 19349663*key(2)), 83492791*key(3)), &
                      size(self%buckets))

    end function hash

    !> Allocate
    pure subroutine chash_tbl_allocate(self, m)
        class(chash_tbl), intent(inout) :: self
        integer, intent(in) :: m

        allocate (self%buckets(0:m - 1))

    end subroutine chash_tbl_allocate

    !> Clean up
    pure subroutine clear(self)
        class(chash_tbl), intent(inout) :: self

        if (allocated(self%buckets)) then
            deallocate (self%buckets)
        end if

    end subroutine clear

    !> Zeroing
    pure subroutine zeroing(self)
        class(chash_tbl), intent(inout) :: self

        call self%buckets(:)%zeroing()

    end subroutine zeroing

    !> Push
    pure subroutine set(self, key, value, stat)
        class(chash_tbl), intent(inout) :: self
        integer, intent(in) :: key(3), value
        logical, intent(out) :: stat

        call self%buckets(self%hash(key))%push_back(key, value, stat)

    end subroutine set

end module nnps_compact_hashing
