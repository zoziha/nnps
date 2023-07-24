!> key-value
module nnps_key_value

    use nnps_int_vector, only: int_vector
    implicit none

    private
    public :: key_values

    !> key-value pair
    type key_value
        integer :: key(3)  !! key
        type(int_vector) :: value  !! value
    end type key_value

    !> key-value pairs
    type key_values
        type(key_value), allocatable :: items(:)  !! key-value pair
    contains
        procedure :: push_back, get_value, zeroing !, storage, clear
    end type key_values

contains

    !> constructor
    pure type(key_value) function key_value_constructor(key, value) result(self)
        integer, intent(in) :: key(3)
        integer, intent(in) :: value

        self%key = key
        call self%value%push_back(value)

    end function key_value_constructor

    !> push back key-value pair
    pure subroutine push_back(self, key, value, stat)
        class(key_values), intent(inout) :: self
        integer, intent(in) :: key(3)
        integer, intent(in) :: value
        logical, intent(out) :: stat
        integer :: i

        if (allocated(self%items)) then
            do i = 1, size(self%items)
                if (all(self%items(i)%key == key)) then
                    call self%items(i)%value%push_back(value)
                    if (self%items(i)%value%len == 1) then
                        stat = .true.
                    else
                        stat = .false.
                    end if
                    return
                end if
            end do
            self%items = [self%items, key_value_constructor(key, value)]
            stat = .true.
        else
            allocate (self%items(1), source=key_value_constructor(key, value))
            stat = .true.
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
            if (all(self%items(i)%key == key)) then
                ptr => self%items(i)%value%items(1:self%items(i)%value%len)
                return
            end if
        end do

    end subroutine get_value

end module nnps_key_value