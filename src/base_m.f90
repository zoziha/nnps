module base_m

    implicit none
    
    !> 二维点
    type point_t
        real :: x  !! 点的横坐标
        real :: y  !! 点的纵坐标
    end type point_t

    !> 矩形
    type rectangle_t
        
        real :: x  !! 矩形形心横坐标
        real :: y  !! 矩形形心纵坐标
        real :: w  !! 矩形宽度
        real :: h  !! 矩形高度

    contains

        procedure :: contains   => rectangle_t_contains
        procedure :: intersects => rectangle_t_intersects

    end type rectangle_t

    !> 圆形
    type circle_t

        real :: x, y, r

    contains

        procedure :: contains   => circle_t_contains
        procedure :: intersects => circle_t_intersects

    end type circle_t

contains

    !> 查询粒子是否在矩形域内
    logical function rectangle_t_contains(self, point) result(contains)

        class(rectangle_t), intent(in) :: self
        type(point_t), intent(in) :: point

        associate (x => self%x, &
                   y => self%y, &
                   w => self%w, &
                   h => self%h)
            
            !@note: 注意这里和其他地方的等于号！:
            ! - 如果是在碰撞检测中，很可能需要等于号来检测碰撞;
            ! - 如果是在 SPH 粒子法中，很可能不需要等于号，因为在查找域边缘的目标值往往为0.
            contains = (point%x >= x - 0.5*w) .and. &
                       (point%x <= x + 0.5*w) .and. &
                       (point%y >= y - 0.5*h) .and. &
                       (point%y <= y + 0.5*h)

        end associate

    end function rectangle_t_contains

    !> 查询几何形状是否有交集
    logical function rectangle_t_intersects(self, range) result(intersects)

        class(rectangle_t), intent(in) :: self
        type(rectangle_t), intent(in) :: range

        associate (left   => self%x - 0.5*self%w, &
                   right  => self%x + 0.5*self%w, &
                   bottom => self%y - 0.5*self%h, &
                   top    => self%y + 0.5*self%h, &
                   left_   => range%x - 0.5*range%w, &
                   right_  => range%x + 0.5*range%w, &
                   bottom_ => range%y - 0.5*range%h, &
                   top_    => range%y + 0.5*range%h  )

            intersects = (left   <= right_ ) .or. &
                         (right  >= left_  ) .or. &
                         (bottom <= top_   ) .or. &
                         (top    >= bottom_)

        end associate

    end function rectangle_t_intersects

    !> 查询粒子是否在圆形内
    logical function circle_t_contains(self, point) result(contains)

        class(circle_t), intent(in) :: self
        type(point_t), intent(in) :: point

        contains = hypot(self%x - point%x, self%y - point%y) <= self%r

    end function circle_t_contains

    !> 查询几何形状是否有交集
    logical function circle_t_intersects(self, range) result(intersects)

        class(circle_t), intent(in) :: self
        type(rectangle_t), intent(in) :: range

        real :: x_dist, y_dist

        ! 形心距离
        x_dist = abs(range%x - self%x)
        y_dist = abs(range%y - self%y)

        associate(r => self%r , &
                  w => range%w, &
                  h => range%h  )
            associate(edges => hypot(x_dist - 0.5*w, y_dist - 0.5*h))

                !> 圆形边缘在矩形外
                !> no intersection
                if (x_dist > r + 0.5*w .or. y_dist > r + 0.5*h) then
                    intersects = .false.
                    return
                end if

                !> 圆形形心在矩形内
                !> intersection within the circle
                if (x_dist <= 0.5*w .and. y_dist <= 0.5*h) then
                    intersects = .true.
                    return
                end if

                !> 圆形边缘在矩形内
                !> intersection on the edge of the circle
                intersects = edges <= self%r

            end associate
        end associate

    end function circle_t_intersects

end module base_m