!> Vector data vector
module nnps_vector

    use nnps_kinds, only: rk
    implicit none

    private
    public :: vector, vector_finalizer

    !> Vector data vector
    type vector
        integer :: len = 0  !! 有效向量长度
        integer, private :: cap = 0  !! 向量容量
        integer, private :: dim  !! 物理场维度
        integer, allocatable :: items(:)  !! 整型数组
        real(rk), allocatable :: ritems(:)  !! 实型数组
    contains
        procedure :: init
        procedure :: push, merge, storage
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

    !> 向量合并
    subroutine merge(self, that)
        class(vector), intent(inout) :: self
        type(vector), intent(in) :: that(1:)
        integer :: idx(size(that)), i

        idx(1) = that(1)%len
        do i = 2, size(that)
            idx(i) = idx(i - 1) + that(i)%len
        end do

        do while (idx(size(that)) > self%cap)
            call self%extend()
        end do

        !$omp parallel do private(i)
        do i = 1, size(that)
            if (that(i)%len == 0) cycle
            self%items((idx(i) - that(i)%len)*2 + 1:idx(i)*2) = &
                that(i)%items(1:that(i)%len*2)
            self%ritems((idx(i) - that(i)%len)*(self%dim + 1) + 1:idx(i)*(self%dim + 1)) = &
                that(i)%ritems(1:(self%dim + 1)*that(i)%len)
        end do
        self%len = idx(size(that))

    end subroutine merge

    !> Storage size
    pure integer function storage(self)
        class(vector), intent(in) :: self

        storage = storage_size(self) + storage_size(self%items)*size(self%items) + &
                  storage_size(self%ritems)*size(self%ritems)

    end function storage

    !> 向量清空
    elemental subroutine vector_finalizer(self)
        type(vector), intent(inout) :: self

        if (allocated(self%items)) deallocate (self%items)
        if (allocated(self%ritems)) deallocate (self%ritems)

    end subroutine vector_finalizer

end module nnps_vector
