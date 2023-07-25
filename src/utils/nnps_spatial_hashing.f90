!> spatial Hashing
module nnps_spatial_hashing

    use nnps_key_value, only: key_values, key_values_finalizer
    implicit none

    private
    public :: shash_tbl, shash_tbl_finalizer

    !> spatial Hashing
    type shash_tbl
        type(key_values), allocatable :: buckets(:)  !! Buckets
    contains
        procedure :: allocate => shash_tbl_allocate, zeroing, set, hash
    end type shash_tbl

contains

    !> Hashing function
    pure integer function hash(self, key)
        class(shash_tbl), intent(in) :: self
        integer, intent(in) :: key(3)

        hash = modulo(ieor(ieor(73856093*key(1), 19349663*key(2)), 83492791*key(3)), &
                      size(self%buckets))

    end function hash

    !> Allocate
    pure subroutine shash_tbl_allocate(self, m)
        class(shash_tbl), intent(inout) :: self
        integer, intent(in) :: m

        allocate (self%buckets(0:m - 1))

    end subroutine shash_tbl_allocate

    !> Clean up
    pure subroutine shash_tbl_finalizer(self)
        type(shash_tbl), intent(inout) :: self

        if (allocated(self%buckets)) then
            call key_values_finalizer(self%buckets(:))
            deallocate (self%buckets)
        end if

    end subroutine shash_tbl_finalizer

    !> Zeroing
    pure subroutine zeroing(self)
        class(shash_tbl), intent(inout) :: self

        call self%buckets(:)%zeroing()

    end subroutine zeroing

    !> Push
    pure subroutine set(self, key, value, stat)
        class(shash_tbl), intent(inout) :: self
        integer, intent(in) :: key(3), value
        logical, intent(out) :: stat

        call self%buckets(self%hash(key))%push_back(key, value, stat)

    end subroutine set

end module nnps_spatial_hashing
