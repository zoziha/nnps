!> key-value
module nnps_key_value

    use nnps_int_vector, only: int_vector, int_vector_finalizer
    implicit none

    private
    public :: key_values, key_values_finalizer, recycle

    integer, parameter :: huge_value = huge(0)      !! 最大值, 用于标记 (一般不可能有粒子进入 huge 对应的网格)

    !> 键值对
    !> key-value pair
    type key_value
        integer :: key(3)           !! key
        type(int_vector) :: value   !! value
    end type key_value

    !> key-value pairs
    type key_values
        type(key_value), allocatable :: items(:)    !! key-value pair
        logical :: flag = .false.                   !! 是否已被回收, 存在惰性空间
    contains
        procedure :: push_back, get_value, zeroing
    end type key_values

contains

    !> finalizer
    elemental subroutine key_value_finalizer(self)
        type(key_value), intent(inout) :: self

        call int_vector_finalizer(self%value)

    end subroutine key_value_finalizer

    !> finalizer
    elemental subroutine key_values_finalizer(self)
        type(key_values), intent(inout) :: self

        if (allocated(self%items)) then
            call key_value_finalizer(self%items)
            deallocate (self%items)
        end if

    end subroutine key_values_finalizer

    !> constructor
    pure type(key_value) function key_value_constructor(key, value) result(self)
        integer, intent(in) :: key(3)
        integer, intent(in) :: value

        self%key = key
        call self%value%push_back(value)

    end function key_value_constructor

    !> push back key-value pair
    pure subroutine push_back(self, key, value, istat)
        class(key_values), intent(inout) :: self
        integer, intent(in) :: key(3)
        integer, intent(in) :: value
        integer, intent(out) :: istat  !! 0: first, 1: lock, 2: not first and not lock
        integer :: i

        if (allocated(self%items)) then
            ! 往活跃空间中插入
            do i = 1, size(self%items)
                if (all(self%items(i)%key == key)) then  ! TODO: 耗时
                    call self%items(i)%value%push_back(value)
                    if (self%items(i)%value%len == 1) then
                        istat = 0
                    else
                        istat = 2
                    end if
                    return
                end if
            end do
            ! 往惰性空间中插入
            if (self%flag) then
                do i = 1, size(self%items)
                    if (all(self%items(i)%key == huge_value)) then
                        self%items(i) = key_value_constructor(key, value)  ! 惰性空间被复用
                        istat = 0
                        return
                    end if
                end do
                self%flag = .false.
            end if
            ! 开辟新空间并插入
            self%items = [self%items, key_value_constructor(key, value)]
            istat = 0
        else
            ! 初始化空间并插入
            allocate (self%items(1), source=key_value_constructor(key, value))
            istat = 0
        end if

    end subroutine push_back

    !> zeroing key-value pair
    elemental subroutine zeroing(self)
        class(key_values), intent(inout) :: self

        if (allocated(self%items)) then
            where (self%items(:)%value%len /= 0) self%items(:)%value%len = 0
        end if

    end subroutine zeroing

    !> get value @todo improve performance
    subroutine get_value(self, key, ptr)
        class(key_values), intent(in), target :: self
        integer, intent(in) :: key(3)
        integer, pointer, intent(inout) :: ptr(:)
        integer :: i

        if (.not. allocated(self%items)) return
        do i = 1, size(self%items)
            if (all(self%items(i)%key == key)) then ! TODO: 耗时
                ptr => self%items(i)%value%items(1:self%items(i)%value%len)
                return
            end if
        end do

    end subroutine get_value

    !> 空键值对回收机制: 添加标记 @note 应该低频率实施该回收过程
    elemental subroutine recycle(self)
        type(key_values), intent(inout) :: self
        integer :: i

        if (.not. allocated(self%items)) return
        if (size(self%items) > 10) then                     ! 回收较大桶的内存
            deallocate (self%items)
            return
        end if

        do i = 1, size(self%items)
            if (self%items(i)%value%len == 0) then
                self%items(i)%key = huge_value              ! 使用极大值标记空键值对
                if (.not. self%flag) self%flag = .true.
            end if
        end do

    end subroutine recycle

end module nnps_key_value
