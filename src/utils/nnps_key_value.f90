!> key-value
module nnps_key_value

    use nnps_int_vector, only: int_vector
    implicit none

    private
    public :: key_value

    !> key-value pair
    type key_value
        integer, allocatable :: key(:)  !! key
        type(int_vector), allocatable :: value(:)  !! value
    contains
        procedure :: push_back, get_value, zeroing !, storage, clear
    end type key_value

contains

    !> push back key-value pair
    pure subroutine push_back(self, key, value, stat)
        class(key_value), intent(inout) :: self
        integer, intent(in) :: key(3)
        integer, intent(in) :: value
        logical, intent(out) :: stat
        integer :: i

        if (allocated(self%key)) then
            do i = 1, size(self%value)
                if (all(self%key(3*i - 2:3*i) == key)) then
                    call self%value(i)%push_back(value)
                    if (self%value(i)%len == 1) then
                        stat = .true.
                    else
                        stat = .false.
                    end if
                    return
                end if
            end do
            self%key = [self%key, key]
            block
                type(int_vector) :: tmp
                call tmp%push_back(value)
                self%value = [self%value, tmp]
                stat = .true.
            end block
        else
            allocate (self%key, source=key)
            allocate (self%value(1))
            call self%value(1)%push_back(value)
            stat = .true.
        end if

    end subroutine push_back

    !> zeroing key-value pair
    elemental subroutine zeroing(self)
        class(key_value), intent(inout) :: self

        if (allocated(self%key)) then
            where (self%value(:)%len /= 0) self%value(:)%len = 0
        end if

    end subroutine zeroing

    !> get value @todo improve performance
    subroutine get_value(self, key, ptr)
        class(key_value), intent(in), target :: self
        integer, intent(in) :: key(3)
        integer, pointer, intent(inout) :: ptr(:)
        integer :: i

        if (.not. allocated(self%key)) return
        do i = 1, size(self%value)
            if (all(self%key(3*i - 2:3*i) == key)) then
                ptr => self%value(i)%items(1:self%value(i)%len)
                return
            end if
        end do

    end subroutine get_value

end module nnps_key_value
