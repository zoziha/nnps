! n叉树
module ntree_m

    use shape_m, only: point_t
    use shape_m, only: shape_t, line_t, rectangle_t, square_t, cube_t, cuboid_t
    implicit none
    private
    
    public :: ntree_t, make_ntree
    
    !> n叉树
    type ntree_t
        class(shape_t), allocatable :: boundary     !! 边界域
        integer :: capacity                         !! 容量 
        logical :: divided                          !! 是否已经分割
        type(point_t), allocatable :: points(:)     !! 点
        type(ntree_t), allocatable :: children(:)   !! 子树
    contains
        procedure :: insert => ntree_t_insert
        procedure :: divide => ntree_t_divide
        procedure :: query  => ntree_t_query
    end type ntree_t
    
contains

    !> 获取树
    pure subroutine make_ntree(boundary, capacity, ntree)
        class(shape_t), intent(in ) :: boundary
        integer       , intent(in ) :: capacity
        type(ntree_t) , intent(out) :: ntree
        ntree%boundary = boundary
        ntree%capacity = capacity
        ntree%divided  = .false.
        allocate(ntree%points(0))
    end subroutine make_ntree
    
    !> 插入点
    recursive subroutine ntree_t_insert(self, point, done)
        class(ntree_t), intent(inout) :: self
        type(point_t) , intent(in )   :: point
        logical       , intent(out), optional :: done
        integer :: i
        
        done = .false.
        if (.not.self%boundary%contains(point)) return
        
        if (size(self%points) < self%capacity) then
            self%points = [self%points, point]
            done = .true.
        else
            if (.not.self%divided) call self%divide()
            do i = 1, size(self%children)
                call self%children(i)%insert(point, done)
                if (done) return
            end do
        end if
        
    end subroutine ntree_t_insert
    
    !> 分割
    pure subroutine ntree_t_divide(self)
        class(ntree_t), intent(inout) :: self
        associate(boundary => self%boundary)
        select type (boundary)
        type is (line_t)
            allocate(self%children(2))
            associate(x => boundary%center(1), &
                      y => boundary%length)
            call make_ntree(line_t([x - y/4], y/2), self%capacity, self%children(1))
            call make_ntree(line_t([x + y/4], y/2), self%capacity, self%children(2))
            end associate
        type is (rectangle_t)
            allocate(self%children(4))
            associate(x => boundary%center(1), &
                      y => boundary%center(2), &
                      l => boundary%length(1), &
                      w => boundary%length(2))
            call make_ntree(rectangle_t([x - l/4, y - w/4], [l/2, w/2]), self%capacity, self%children(1))
            call make_ntree(rectangle_t([x - l/4, y + w/4], [l/2, w/2]), self%capacity, self%children(2))
            call make_ntree(rectangle_t([x + l/4, y - w/4], [l/2, w/2]), self%capacity, self%children(3))
            call make_ntree(rectangle_t([x + l/4, y + w/4], [l/2, w/2]), self%capacity, self%children(4))
            end associate
        type is (square_t)
            allocate(self%children(4))
            associate(x => boundary%center(1), &
                      y => boundary%center(2), &
                      l => boundary%length)
            call make_ntree(square_t([x - l/4, y - l/4], l/2), self%capacity, self%children(1))
            call make_ntree(square_t([x - l/4, y + l/4], l/2), self%capacity, self%children(2))
            call make_ntree(square_t([x + l/4, y - l/4], l/2), self%capacity, self%children(3))
            call make_ntree(square_t([x + l/4, y + l/4], l/2), self%capacity, self%children(4))
            end associate
        type is (cube_t)
            allocate(self%children(8))
            associate(x => boundary%center(1), &
                      y => boundary%center(2), &
                      z => boundary%center(3), &
                      l => boundary%length)
            call make_ntree(cube_t([x - l/4, y - l/4, z - l/4], l/2), self%capacity, self%children(1))
            call make_ntree(cube_t([x - l/4, y - l/4, z + l/4], l/2), self%capacity, self%children(2))
            call make_ntree(cube_t([x - l/4, y + l/4, z - l/4], l/2), self%capacity, self%children(3))
            call make_ntree(cube_t([x - l/4, y + l/4, z + l/4], l/2), self%capacity, self%children(4))
            call make_ntree(cube_t([x + l/4, y - l/4, z - l/4], l/2), self%capacity, self%children(5))
            call make_ntree(cube_t([x + l/4, y - l/4, z + l/4], l/2), self%capacity, self%children(6))
            call make_ntree(cube_t([x + l/4, y + l/4, z - l/4], l/2), self%capacity, self%children(7))
            call make_ntree(cube_t([x + l/4, y + l/4, z + l/4], l/2), self%capacity, self%children(8))
            end associate
        type is (cuboid_t)
            allocate(self%children(8))
            associate(x => boundary%center(1), &
                      y => boundary%center(2), &
                      z => boundary%center(3), &
                      l => boundary%length(1), &
                      w => boundary%length(2), &
                      h => boundary%length(3))
            call make_ntree(cuboid_t([x - l/4, y - l/4, z - w/4], [l/2, w/2, h/2]), self%capacity, self%children(1))
            call make_ntree(cuboid_t([x - l/4, y - l/4, z + w/4], [l/2, w/2, h/2]), self%capacity, self%children(2))
            call make_ntree(cuboid_t([x - l/4, y + l/4, z - w/4], [l/2, w/2, h/2]), self%capacity, self%children(3))
            call make_ntree(cuboid_t([x - l/4, y + l/4, z + w/4], [l/2, w/2, h/2]), self%capacity, self%children(4))
            call make_ntree(cuboid_t([x + l/4, y - l/4, z - w/4], [l/2, w/2, h/2]), self%capacity, self%children(5))
            call make_ntree(cuboid_t([x + l/4, y - l/4, z + w/4], [l/2, w/2, h/2]), self%capacity, self%children(6))
            call make_ntree(cuboid_t([x + l/4, y + l/4, z - w/4], [l/2, w/2, h/2]), self%capacity, self%children(7))
            call make_ntree(cuboid_t([x + l/4, y + l/4, z + w/4], [l/2, w/2, h/2]), self%capacity, self%children(8))
            end associate
        end select
        end associate
        self%divided = .true.
    end subroutine ntree_t_divide

    !> 查询点
    recursive subroutine ntree_t_query(self, range, found)
        class(ntree_t), intent(inout) :: self
        class(shape_t), intent(in   ) :: range
        type(point_t) , intent(inout), allocatable :: found(:)
        integer :: i
        if (.not.range%intersects(self%boundary)) return
        do i = 1, size(self%points)
            if (range%contains(self%points(i))) then
                if (.not.allocated(found)) then
                    allocate(found(1), source=self%points(i))
                else
                    found = [found, self%points(i)]
                end if
            end if
        end do
        if (self%divided) then
            do i = 1, size(self%children)
                call self%children(i)%query(range, found)
            end do
        end if
    end subroutine ntree_t_query

end module ntree_m
