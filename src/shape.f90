! 形状 (暂不服从开闭原则)
module shape_m

    use, intrinsic :: iso_fortran_env, only: real32
    implicit none
    private
    
    public :: point_t, shape_t, cuboid_t, cube_t, square_t, rectangle_t, line_t, sphere_t, circle_t
    
    !> 点的坐标
    !> 支持1~3维
    type point_t
        real(real32), allocatable :: x(:) !! 坐标
    end type point_t
    
    ! 形状类型
    type, abstract :: shape_t
    contains
        procedure(shape_t_contains  ), deferred :: contains     !! 查询点是否在形状内
        procedure(shape_t_intersects), deferred :: intersects   !! 查询是否相交
    end type shape_t
    
    abstract interface
        logical function shape_t_contains(self, point)
            import shape_t, point_t
            class(shape_t), intent(in) :: self
            type(point_t) , intent(in) :: point
        end function shape_t_contains 
        logical function shape_t_intersects(self, other)
            import shape_t
            class(shape_t), intent(in) :: self
            class(shape_t), intent(in) :: other
        end function shape_t_intersects
    end interface
    
    type, extends(shape_t) :: line_t
        real(real32) :: center(1)
        real(real32) :: length
    contains
        procedure :: contains   => line_t_contains
        procedure :: intersects => line_t_intersects
    end type line_t
    
    type, extends(shape_t) :: circle_t
        real(real32) :: center(2)
        real(real32) :: radius
    contains
        procedure :: contains   => circle_t_contains
        procedure :: intersects => circle_t_intersects
    end type circle_t
    
    type, extends(shape_t) :: square_t
        real(real32) :: center(2)
        real(real32) :: length
    contains
        procedure :: contains   => square_t_contains
        procedure :: intersects => square_t_intersects
    end type square_t

    type, extends(shape_t) :: rectangle_t
        real(real32) :: center(2)
        real(real32) :: length(2)
    contains
        procedure :: contains   => rectangle_t_contains
        procedure :: intersects => rectangle_t_intersects
    end type rectangle_t
    
    type, extends(shape_t) :: sphere_t
        real(real32) :: center(3)
        real(real32) :: radius
    contains
        procedure :: contains   => sphere_t_contains
        procedure :: intersects => sphere_t_intersects
    end type sphere_t

    type, extends(shape_t) :: cube_t
        real(real32) :: center(3)
        real(real32) :: length
    contains
        procedure :: contains   => cube_t_contains
        procedure :: intersects => cube_t_intersects
    end type cube_t
    
    type, extends(shape_t) :: cuboid_t
        real(real32) :: center(3)
        real(real32) :: length(3)
    contains
        procedure :: contains   => cuboid_t_contains
        procedure :: intersects => cuboid_t_intersects
    end type cuboid_t    
    
contains
            
    pure real(real32) function func1(x, y)
        real(real32), intent(in) :: x, y
        func1 = x - y/2
    end function func1
    
    pure real(real32) function func2(x, y)
        real(real32), intent(in) :: x, y
        func2 = x + y/2
    end function func2

    logical pure function line_t_contains(self, point)
        class(line_t), intent(in) :: self
        type(point_t), intent(in) :: point
        select case (size(point%x))
        case (1)
            associate(x => self%center(1), l => self%length)
            line_t_contains = (point%x(1) >= x - l/2) .and. &
                              (point%x(1) <= x + l/2)
            end associate
        end select
    end function line_t_contains
    
    logical pure function circle_t_contains(self, point)
        class(circle_t), intent(in) :: self
        type (point_t ), intent(in) :: point
        select case (size(point%x))
        case (2)
            associate(x => self%center(1), y => self%center(2), r => self%radius)
            circle_t_contains = hypot(point%x(1) - x, point%x(2) - y) <= r
            end associate
        end select
    end function circle_t_contains
    
    logical pure function square_t_contains(self, point)
        class(square_t), intent(in) :: self
        type (point_t ), intent(in) :: point
        select case (size(point%x))
        case (2)
            associate(x => self%center(1), y => self%center(2), l => self%length)
            square_t_contains = (point%x(1) >= x - l/2) .and. &
                (point%x(1) <= x + l/2) .and. (point%x(2) >= y - l/2) .and. &
                (point%x(2) <= y + l/2)
            end associate
        end select
    end function square_t_contains
    
    logical pure function rectangle_t_contains(self, point)
        class(rectangle_t), intent(in) :: self
        type (point_t    ), intent(in) :: point
        select case (size(point%x))
        case (2)
            associate(x => self%center(1), y => self%center(2), l => self%length(1), w => self%length(2))
            rectangle_t_contains = (point%x(1) >= x - l/2) .and. &
                (point%x(1) <= x + l/2) .and. (point%x(2) >= y - w/2) .and. &
                (point%x(2) <= y + w/2)
            end associate
        end select
    end function rectangle_t_contains
    
    logical pure function sphere_t_contains(self, point)
        class(sphere_t), intent(in) :: self
        type (point_t ), intent(in) :: point
        select case (size(point%x))
        case (3)
            associate(x => self%center(1), y => self%center(2), z => self%center(3), r => self%radius)
            sphere_t_contains = norm2([point%x(1) - x, point%x(2) - y, point%x(3) - z]) <= r
            end associate
        end select
    end function sphere_t_contains
    
    logical pure function cube_t_contains(self, point)
        class(cube_t), intent(in) :: self
        type(point_t), intent(in) :: point
        select case (size(point%x))
        case (3)
            associate(x => self%center(1), y => self%center(2), z => self%center(3), l => self%length)
            cube_t_contains = (point%x(1) >= x - l/2) .and. &
                (point%x(1) <= x + l/2) .and. (point%x(2) >= y - l/2) .and. &
                (point%x(2) <= y + l/2) .and. (point%x(3) >= z - l/2) .and. &
                (point%x(3) <= z + l/2)
            end associate
        end select
    end function cube_t_contains
    
    logical pure function cuboid_t_contains(self, point)
        class(cuboid_t), intent(in) :: self
        type (point_t ), intent(in) :: point
        select case (size(point%x))
        case (3)
            associate(x => self%center(1), y => self%center(2), z => self%center(3), &
                      l => self%length(1), w => self%length(2), h => self%length(3))
            cuboid_t_contains = (point%x(1) >= x - l/2) .and. &
                (point%x(1) <= x + l/2) .and. (point%x(2) >= y - w/2) .and. &
                (point%x(2) <= y + w/2) .and. (point%x(3) >= z - h/2) .and. &
                (point%x(3) <= z + h/2)
            end associate
        end select
    end function cuboid_t_contains
    
    logical pure function line_t_intersects(self, other)
        class(line_t ), intent(in) :: self
        class(shape_t), intent(in) :: other
        select type (other)
        type is (line_t)
            associate(left  => func1(self%center(1) , self%length ), &
                      right => func2(self%center(1) , self%length ), &
                      left_ => func1(other%center(1), other%length), &
                      right_=> func2(other%center(1), other%length))
            line_t_intersects = (left <= right_) .or. (right >= left_)
            end associate
        end select
    end function line_t_intersects
    
    logical pure function circle_t_intersects(self, other)
        class(circle_t), intent(in) :: self
        class(shape_t ), intent(in) :: other
        select type (other)
        type is (square_t)
            associate(x_dist => abs(other%center(1) - self%center(1)), y_dist => abs(other%center(2) - self%center(2)), &
                      r => self%radius, l => other%length)
            if (x_dist > r + l/2 .or. y_dist > r + l/2) then
                circle_t_intersects = .false.
                return
            end if
            if (x_dist <= l/2 .and. y_dist <= l/2) then
                circle_t_intersects = .true.
                return
            end if
            circle_t_intersects = hypot(x_dist - l/2, y_dist - l/2) <= r
            end associate
        type is (rectangle_t)
            associate(x_dist => abs(other%center(1) - self%center(1)), y_dist => abs(other%center(2) - self%center(2)), &
                      r => self%radius, l => other%length(1), w => other%length(2))
            if (x_dist > r + l/2 .or. y_dist > r + w/2) then
                circle_t_intersects = .false.
                return
            end if
            if (x_dist <= l/2 .and. y_dist <= w/2) then
                circle_t_intersects = .true.
                return
            end if
            circle_t_intersects = hypot(x_dist - l/2, y_dist - w/2) <= r
            end associate
        end select
    end function circle_t_intersects
    
    logical pure function square_t_intersects(self, other)
        class(square_t), intent(in) :: self
        class(shape_t ), intent(in) :: other
        select type (other)
        type is (square_t)
            associate(left   => func1(self%center(1) , self%length ), &
                      right  => func2(self%center(1) , self%length ), &
                      bottom => func1(self%center(2) , self%length ), &
                      top    => func2(self%center(2) , self%length ), &
                      left_  => func1(other%center(1), other%length), &
                      right_ => func2(other%center(1), other%length), &
                      bottom_=> func1(other%center(2), other%length), &
                      top_   => func2(other%center(2), other%length))
            square_t_intersects = (left  <= right_ ) &
                             .or. (right >= left_  ) &
                             .or. (bottom<= top_   ) &
                             .or. (top   >= bottom_)
            end associate
        type is (rectangle_t)
            associate(left   => func1(self%center(1) , self%length    ), &
                      right  => func2(self%center(1) , self%length    ), &
                      bottom => func1(self%center(2) , self%length    ), &
                      top    => func2(self%center(2) , self%length    ), &
                      left_  => func1(other%center(1), other%length(1)), &
                      right_ => func2(other%center(1), other%length(1)), &
                      bottom_=> func1(other%center(2), other%length(2)), &
                      top_   => func2(other%center(2), other%length(2)))
            square_t_intersects = (left  <= right_ ) &
                             .or. (right >= left_  ) &
                             .or. (bottom<= top_   ) &
                             .or. (top   >= bottom_)
            end associate
        end select
    end function square_t_intersects
    
    logical pure function rectangle_t_intersects(self, other)
        class(rectangle_t), intent(in) :: self
        class(shape_t    ), intent(in) :: other
        select type (other)
        type is (square_t)
            associate(left   => func1(self%center(1) , self%length(1)), &
                      right  => func2(self%center(1) , self%length(1)), &
                      bottom => func1(self%center(2) , self%length(2)), &
                      top    => func2(self%center(2) , self%length(2)), &
                      left_  => func1(other%center(1), other%length  ), &
                      right_ => func2(other%center(1), other%length  ), &
                      bottom_=> func1(other%center(2), other%length  ), &
                      top_   => func2(other%center(2), other%length  ))
            rectangle_t_intersects = (left  <= right_ ) &
                                .or. (right >= left_  ) &
                                .or. (bottom<= top_   ) &
                                .or. (top   >= bottom_)
            end associate
        type is (rectangle_t)
            associate(left   => func1(self%center(1) , self%length(1) ), &
                      right  => func2(self%center(1) , self%length(1) ), &
                      bottom => func1(self%center(2) , self%length(2) ), &
                      top    => func2(self%center(2) , self%length(2) ), &
                      left_  => func1(other%center(1), other%length(1)), &
                      right_ => func2(other%center(1), other%length(1)), &
                      bottom_=> func1(other%center(2), other%length(2)), &
                      top_   => func2(other%center(2), other%length(2)))
            rectangle_t_intersects = (left  <= right_ ) &
                                .or. (right >= left_  ) &
                                .or. (bottom<= top_   ) &
                                .or. (top   >= bottom_)
            end associate
        end select
    end function rectangle_t_intersects    
    
    logical pure function sphere_t_intersects(self, other)
        class(sphere_t), intent(in) :: self
        class(shape_t ), intent(in) :: other
        select type (other)
        type is (sphere_t)
            associate(r1 => self%radius, r2 => other%radius)
            sphere_t_intersects = norm2(other%center - self%center) <= r1 + r2
            end associate
        type is (cube_t)
            associate(x_dist => abs(other%center(1) - self%center(1)), y_dist => abs(other%center(2) - self%center(2)), &
                      z_dist => abs(other%center(3) - self%center(3)), r => self%radius, l => other%length)
            if (x_dist > r + l/2 .or. y_dist > r + l/2 .or. z_dist > r + l/2) then
                sphere_t_intersects = .false.
                return
            end if
            if (x_dist <= l/2 .and. y_dist <= l/2 .and. z_dist <= l/2) then
                sphere_t_intersects = .true.
                return
            end if
            sphere_t_intersects = norm2([x_dist - l/2, y_dist - l/2, z_dist - l/2]) <= r
            end associate
        type is (cuboid_t)
            associate(x_dist => abs(other%center(1) - self%center(1)), y_dist => abs(other%center(2) - self%center(2)), &
                      z_dist => abs(other%center(3) - self%center(3)), r => self%radius, l => other%length(1), &
                      w => other%length(2), h => other%length(3))
            if (x_dist > r + l/2 .or. y_dist > r + w/2 .or. z_dist > r + h/2) then
                sphere_t_intersects = .false.
                return
            end if
            if (x_dist <= l/2 .and. y_dist <= w/2 .and. z_dist <= h/2) then
                sphere_t_intersects = .true.
                return
            end if
            sphere_t_intersects = norm2([x_dist - l/2, y_dist - w/2, z_dist - h/2]) <= r
            end associate
        end select
    end function sphere_t_intersects
    
    logical pure function cube_t_intersects(self, other)
        class(cube_t ), intent(in) :: self
        class(shape_t), intent(in) :: other
        select type (other)
        type is (cube_t)
            associate(left   => func1(self%center(1) , self%length ), &
                      right  => func2(self%center(1) , self%length ), &
                      bottom => func1(self%center(2) , self%length ), &
                      top    => func2(self%center(2) , self%length ), &
                      back   => func1(self%center(3) , self%length ), &
                      front  => func2(self%center(3) , self%length ), &
                      left_  => func1(other%center(1), other%length), &
                      right_ => func2(other%center(1), other%length), &
                      bottom_=> func1(other%center(2), other%length), &
                      top_   => func2(other%center(2), other%length), &
                      back_  => func1(other%center(3), other%length), &
                      front_ => func2(other%center(3), other%length))
            cube_t_intersects = (left  <= right_ ) &
                           .or. (right >= left_  ) &
                           .or. (bottom<= top_   ) &
                           .or. (top   >= bottom_) &
                           .or. (back  <= front_ ) &
                           .or. (front >= back_  )
            end associate
        type is (cuboid_t)
            associate(left   => func1(self%center(1) , self%length    ), &
                      right  => func2(self%center(1) , self%length    ), &
                      bottom => func1(self%center(2) , self%length    ), &
                      top    => func2(self%center(2) , self%length    ), &
                      back   => func1(self%center(3) , self%length    ), &
                      front  => func2(self%center(3) , self%length    ), &
                      left_  => func1(other%center(1), other%length(1)), &
                      right_ => func2(other%center(1), other%length(1)), &
                      bottom_=> func1(other%center(2), other%length(2)), &
                      top_   => func2(other%center(2), other%length(2)), &
                      back_  => func1(other%center(3), other%length(3)), &
                      front_ => func2(other%center(3), other%length(3)))
            cube_t_intersects = (left  <= right_ ) &
                           .or. (right >= left_  ) &
                           .or. (bottom<= top_   ) &
                           .or. (top   >= bottom_) &
                           .or. (back  <= front_ ) &
                           .or. (front >= back_  )
            end associate
        end select
    end function cube_t_intersects
    
    logical pure function cuboid_t_intersects(self, other)
        class(cuboid_t), intent(in) :: self
        class(shape_t), intent(in) :: other
        select type (other)
        type is (cuboid_t)
            associate(left   => func1(self%center(1) , self%length(1) ), &
                      right  => func2(self%center(1) , self%length(1) ), &
                      bottom => func1(self%center(2) , self%length(2) ), &
                      top    => func2(self%center(2) , self%length(2) ), &
                      back   => func1(self%center(3) , self%length(3) ), &
                      front  => func2(self%center(3) , self%length(3) ), &
                      left_  => func1(other%center(1), other%length(1)), &
                      right_ => func2(other%center(1), other%length(1)), &
                      bottom_=> func1(other%center(2), other%length(2)), &
                      top_   => func2(other%center(2), other%length(2)), &
                      back_  => func1(other%center(3), other%length(3)), &
                      front_ => func2(other%center(3), other%length(3)))
            cuboid_t_intersects = (left  <= right_ ) &
                           .or. (right >= left_  ) &
                           .or. (bottom<= top_   ) &
                           .or. (top   >= bottom_) &
                           .or. (back  <= front_ ) &
                           .or. (front >= back_  )
            end associate
        type is (cube_t)
            associate(left   => func1(self%center(1) , self%length(1)), &
                      right  => func2(self%center(1) , self%length(1)), &
                      bottom => func1(self%center(2) , self%length(2)), &
                      top    => func2(self%center(2) , self%length(2)), &
                      back   => func1(self%center(3) , self%length(3)), &
                      front  => func2(self%center(3) , self%length(3)), &
                      left_  => func1(other%center(1), other%length  ), &
                      right_ => func2(other%center(1), other%length  ), &
                      bottom_=> func1(other%center(2), other%length  ), &
                      top_   => func2(other%center(2), other%length  ), &
                      back_  => func1(other%center(3), other%length  ), &
                      front_ => func2(other%center(3), other%length  ))
            cuboid_t_intersects = (left  <= right_ ) &
                             .or. (right >= left_  ) &
                             .or. (bottom<= top_   ) &
                             .or. (top   >= bottom_) &
                             .or. (back  <= front_ ) &
                             .or. (front >= back_  )
            end associate
        end select
    end function cuboid_t_intersects

end module shape_m