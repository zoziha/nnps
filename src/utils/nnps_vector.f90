!> Vector data vector
module nnps_vector

    use nnps_kinds, only: rk
    implicit none

    private
    public :: vector

    !> Vector data vector
    type vector
        integer :: len = 0  !! 有效向量长度
        integer, private :: cap = 0  !! 向量容量
        integer, private :: dim  !! 物理场维度
        integer, allocatable :: items(:)  !! 整型数组
        real(rk), allocatable :: ritems(:)  !! 实型数组
    contains
        procedure :: init
        procedure :: push
        procedure :: clear
        procedure, private :: extend
    end type vector

contains

    !> 初始化向量
    !> @todo adjust default length
    elemental subroutine init(self, dim, cap)
        class(vector), intent(inout) :: self
        integer, intent(in) :: dim
        integer, intent(in), optional :: cap
        integer, parameter :: default_len = 32

        self%len = 0
        self%dim = dim
        if (present(cap)) then
            self%cap = cap
        else
            self%cap = default_len
        end if
        if (.not. allocated(self%items)) allocate (self%items(2*self%cap), self%ritems((dim + 1)*self%cap))

    end subroutine init

    !> 向量扩容
    pure subroutine extend(self)
        class(vector), intent(inout) :: self
        integer, allocatable :: tmp(:)
        real(rk), allocatable :: rtmp(:)

        self%cap = 2*self%cap
        allocate (tmp(size(self%items)), rtmp(size(self%ritems)))
        self%items = [self%items, tmp]
        self%ritems = [self%ritems, rtmp]

    end subroutine extend

    !> 向量压入
    pure subroutine push(self, items, ritems)
        class(vector), intent(inout) :: self
        integer, intent(in) :: items(2)
        real(rk), intent(in) :: ritems(:)

        if (self%len == self%cap) call self%extend()
        self%len = self%len + 1
        associate (len2 => 2*self%len, dimlen => (self%dim + 1)*self%len)
            self%items(len2 - 1:len2) = items
            self%ritems(dimlen - self%dim:dimlen) = ritems
        end associate

    end subroutine push

    !> 向量清空
    pure subroutine clear(self)
        class(vector), intent(inout) :: self

        deallocate (self%items, self%ritems)
        self%len = 0

    end subroutine clear

end module nnps_vector
