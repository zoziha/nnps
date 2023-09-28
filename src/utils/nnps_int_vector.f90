!> int_vector integer vector
module nnps_int_vector

    implicit none

    private
    public :: int_vector, int_vector_finalizer

    !> int_vector integer vector
    type int_vector
        integer :: len = 0  !! true vector length
        integer, allocatable :: items(:)  !! integer vector
    contains
        procedure :: push_back, storage, push_back_items, merge
        procedure, private :: extend
    end type int_vector

contains

    !> push back 1 item
    pure subroutine push_back(self, item)
        class(int_vector), intent(inout) :: self
        integer, intent(in) :: item

        if (allocated(self%items)) then
            if (self%len == size(self%items)) call self%extend()
            self%len = self%len + 1
            self%items(self%len) = item
        else
            allocate (self%items(1), source=item)
            self%len = 1
        end if

    end subroutine push_back

    !> push back n items
    pure subroutine push_back_items(self, items, n)
        class(int_vector), intent(inout) :: self
        integer, intent(in) :: items(:)
        integer, intent(in) :: n
        integer :: m

        if (allocated(self%items)) then
            m = self%len + n
            do while (m > size(self%items))
                call self%extend()
            end do
            self%items(self%len + 1:m) = items
            self%len = m
        else
            allocate (self%items(n), source=items)
            self%len = n
        end if

    end subroutine push_back_items

    !> extend
    pure subroutine extend(self)
        class(int_vector), intent(inout) :: self
        integer, allocatable :: tmp(:)

        allocate (tmp(size(self%items)))
        self%items = [self%items, tmp]  ! address of self%items may not be changed, which is good

    end subroutine extend

    !> merge
    subroutine merge(self, that)
        class(int_vector), intent(inout) :: self
        type(int_vector), intent(in) :: that(1:)
        integer :: idx(size(that)), i

        idx(1) = that(1)%len
        do i = 2, size(that)
            idx(i) = idx(i - 1) + that(i)%len
        end do

        do while (idx(size(that)) > size(self%items))
            call self%extend()
        end do

        !$omp parallel do private(i)
        do i = 2, size(that)
            if (that(i)%len == 0) cycle
            self%items(idx(i - 1) + 1:idx(i)) = that(i)%items(1:that(i)%len)
        end do
        self%len = idx(size(that))

    end subroutine merge

    !> Storage
    pure integer function storage(self)
        class(int_vector), intent(in) :: self

        storage = storage_size(self%items)*size(self%items)

    end function storage

    !> 向量清空
    elemental subroutine int_vector_finalizer(self)
        type(int_vector), intent(inout) :: self

        if (allocated(self%items)) deallocate (self%items)

    end subroutine int_vector_finalizer

    !> Zeroing
    elemental subroutine zeroing(self)
        class(int_vector), intent(inout) :: self

        self%len = 0

    end subroutine zeroing

end module nnps_int_vector
