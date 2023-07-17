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
        procedure :: init
        procedure :: push, pop
        procedure :: get, set
        procedure :: clear
        procedure, private :: extend
    end type int_vector

contains

    !> 初始化向量
    !> @todo adjust default length
    elemental subroutine init(self, len)
        class(int_vector), intent(inout) :: self
        integer, intent(in), optional :: len

        self%len = 0
        if (.not. allocated(self%items)) then
            if (present(len)) then
                allocate (self%items(len))
            else
                allocate (self%items(64))
            end if
        end if

    end subroutine init

    !> 向量扩容
    pure subroutine extend(self)
        class(int_vector), intent(inout) :: self
        integer, allocatable :: tmp(:)

        allocate (tmp(size(self%items)))
        self%items = [self%items, tmp]

    end subroutine extend

    !> 向量压入
    pure subroutine push(self, item)
        class(int_vector), intent(inout) :: self
        integer, intent(in) :: item

        if (self%len == size(self%items)) call self%extend()
        self%len = self%len + 1
        self%items(self%len) = item

    end subroutine push

    !> 向量弹出
    pure subroutine pop(self, item)
        class(int_vector), intent(inout) :: self
        integer, intent(out), optional :: item

        if (self%len == 0) return
        if (present(item)) item = self%items(self%len)
        self%len = self%len - 1

    end subroutine pop

    !> 向量获取
    pure subroutine get(self, index, item)
        class(int_vector), intent(in) :: self
        integer, intent(in) :: index
        integer, intent(out) :: item

        if (index < 1 .or. index > self%len) return
        item = self%items(index)

    end subroutine get

    !> 向量设置
    pure subroutine set(self, index, item)
        class(int_vector), intent(inout) :: self
        integer, intent(in) :: index
        integer, intent(in) :: item

        if (index < 1 .or. index > self%len) return
        self%items(index) = item

    end subroutine set

    !> 向量清空
    pure subroutine clear(self)
        class(int_vector), intent(inout) :: self

        deallocate (self%items)
        self%len = 0

    end subroutine clear

end module nnps_int_vector
