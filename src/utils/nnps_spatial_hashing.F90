!> 散列哈希表
!> spatial Hashing
module nnps_spatial_hashing

    use nnps_key_value, only: key_values, key_values_finalizer, key_values_recycle => recycle
    implicit none

    private
    public :: shash_tbl, shash_tbl_finalizer

    !> spatial Hashing
    type shash_tbl
        type(key_values), allocatable :: buckets(:)  !! Buckets
    contains
        procedure :: alloc, zeroing, set, hash, recycle
    end type shash_tbl

contains

    !> Hashing function
    pure integer function hash(self, key)
        class(shash_tbl), intent(in) :: self
        integer, intent(in) :: key(3)

        hash = modulo(ieor(ieor(73856093*key(1), 19349663*key(2)), 83492791*key(3)), &
                      size(self%buckets))

    end function hash

    !> 分配哈希桶数量
    pure subroutine alloc(self, m)
        class(shash_tbl), intent(inout) :: self
        integer, intent(in) :: m

        allocate (self%buckets(0:m - 1))

    end subroutine alloc

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
    subroutine set(self, key, value, istat)
        class(shash_tbl), intent(inout) :: self
        integer, intent(in) :: key(3), value
        integer, intent(out) :: istat  !! 0: first, 1: lock, 2: not first and not lock

        call self%buckets(self%hash(key))%push_back(key, value, istat)

    end subroutine set

    !> 回收内存
    pure subroutine recycle(self)
        class(shash_tbl), intent(inout) :: self

        call key_values_recycle(self%buckets)

    end subroutine recycle

end module nnps_spatial_hashing
