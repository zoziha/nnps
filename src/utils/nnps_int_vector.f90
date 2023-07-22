!> int_vector integer vector
module nnps_int_vector

    implicit none

    private
    public :: int_vector

    !> int_vector integer vector
    type int_vector
        integer :: len = 0  !! 有效向量长度
        integer, allocatable :: items(:)  !! 整型数组
    contains
        procedure :: push, storage
        procedure :: clear
    end type int_vector

contains

    !> 向量压入
    pure subroutine push(self, item)
        class(int_vector), intent(inout) :: self
        integer, intent(in) :: item

        if (allocated(self%items)) then
            if (self%len == size(self%items)) then
                self%len = self%len + 1
                self%items = [self%items, item]
            else
                self%len = self%len + 1
                self%items(self%len) = item
            end if
        else
            allocate(self%items(1), source=item)
            self%len = 1
        end if

    end subroutine push

    !> Storage
    pure integer function storage(self)
        class(int_vector), intent(in) :: self

        storage = storage_size(self) + storage_size(self%items)*size(self%items)

    end function storage

    !> 向量清空
    elemental subroutine clear(self)
        class(int_vector), intent(inout) :: self

        deallocate (self%items)

    end subroutine clear

    !> Zeroing
    elemental subroutine zeroing(self)
        class(int_vector), intent(inout) :: self

        self%len = 0

    end subroutine zeroing

end module nnps_int_vector
