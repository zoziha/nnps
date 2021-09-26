module utils

    implicit none
    
    !> 二维点
    type point_t
        real :: x, y
    end type point_t

    !> 矩形
    type rectangle_t

        real :: x, y        !! 左下角顶点坐标
        real :: w, h        !! 宽度, 高度

    contains

        procedure :: contains
        ! procedure :: intersects

    end type rectangle_t

contains

    !> 查询粒子是否在矩形域内
    logical function contains(self, point)

        class(rectangle_t), intent(inout) :: self
        type(point_t), intent(in) :: point

        associate (x => self%x, &
                   y => self%y, &
                   w => self%w, &
                   h => self%h)

            contains = (point%x >= x - 0.5*w) .and. &
                       (point%x <= x + 0.5*w) .and. &
                       (point%y >= y - 0.5*h) .and. &
                       (point%y <= y + 0.5*h)

        end associate

    end function contains

    ! !> 
    ! logical function intersects(self, range)

    !     class(rectangle_t), intent(inout) :: self
    !     type(rectangle_t), intent(in) :: range

    !     associate (x => self%boundary%x, &
    !                y => slef%boundary%y, &
    !                w => self%boundary%w, &
    !                z => self%boundary%h)

    !         intersects = (range%x - range%w >= x + w) .or. &
    !                      (range%x + w       <= x - w) .or. &
    !                      (range%y - range%h >= y - h) .or. &
    !                      (range%y + h       <= y + h)

    !     end associate

    ! end function intersects

end module utils