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
        procedure :: allocate => shash_tbl_allocate, zeroing, set, hash, &
            storage, activated_buckets
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

    !> Storage
    integer function storage(self)
        class(shash_tbl), intent(in) :: self
        integer :: i

        storage = storage_size(self%buckets)*size(self%buckets)
        do i = 0, size(self%buckets) - 1
            storage = storage + self%buckets(i)%storage()
        end do

    end function storage

    !> Get number of activated buckets (activated buckets, recyclable buckets)
    function activated_buckets(self)
        class(shash_tbl), intent(in) :: self
        integer, dimension(2) :: activated_buckets
        integer :: i

        activated_buckets = 0
        do i = 0, size(self%buckets) - 1
            if (allocated(self%buckets(i)%items)) then
                activated_buckets(1) = activated_buckets(1) + 1
                if (all(self%buckets(i)%items(:)%value%len == 0)) then
                    activated_buckets(2) = activated_buckets(2) + 1
                end if
            end if
        end do

    end function activated_buckets

end module nnps_spatial_hashing
