!> 四叉树建模
module quad

    use utils, only: rectangle_t, point_t
    use stdlib_strings, only: to_string
    implicit none
    
    !> 四叉树
    type quad_tree_t

        type(rectangle_t) :: boundary
        integer :: capacity                             !! SPH粒子法中，该变量应该设为1
        type(point_t), allocatable :: points(:)
        logical :: divided
        type(quad_tree_t), allocatable :: children(:)

        character(:), allocatable :: code               !! 调试目标

    contains

        procedure :: constructor
        procedure :: insert
        procedure :: subdivide
        procedure :: show
        procedure :: delete_tree

    end type quad_tree_t

contains

    !> 构建四叉树
    subroutine constructor(self, boundary, n)

        class(quad_tree_t), intent(inout) :: self
        type(rectangle_t), intent(in) :: boundary
        integer, intent(in) :: n
         
        self%boundary = boundary
        self%capacity = n
        allocate(self%points(0))
        ! allocate(character(0) :: self%code)
        self%divided  = .false.

        if (boundary%w /= boundary%h) error stop "*<ERROR>* heigh /= width!"

    end subroutine constructor

    !\TODO: 确认是否需要递归 (是)
    !> 插入点
    recursive logical function insert(self, point, code) result(done)

        class(quad_tree_t), intent(inout) :: self
        type(point_t), intent(in) :: point
        character(*), intent(in)  :: code       !! 调试目标

        integer :: i

        done      = .false.
        self%code = code

        !> 判断粒子是否处于矩形内部
        if (.not.self%boundary%contains(point)) then
            return
        end if


        if (size(self%points) < self%capacity) then

            !> 添加粒子
            self%points = [self%points, point]
            done = .true.
        
        !> 当前节点的矩形最大粒子容量已满，指向下属节点
        else

            !> 查询是否有下属节点，无则创建！
            if (.not.self%divided) then

                allocate( self%children(4) )

                call self%subdivide()

            end if

            !\TODO: 可能无法运行, Fortran语法？(可以使用递归)
            !> 在下属节点依次尝试插入粒子
            do i = 1, size(self%children)

                if (self%children(i)%insert(point, self%code//to_string(i) )) then
                    done = .true.
                    exit
                end if

            end do

        end if

    end function insert

    !> 细分
    subroutine subdivide(self)

        class(quad_tree_t), intent(inout) :: self
        type(rectangle_t) :: nw, ne, sw, se

        associate (x => self%boundary%x, &
                   y => self%boundary%y, &
                   w => self%boundary%w, &
                   h => self%boundary%h  )

            !> children: 1-ne(东北), 2-nw(西北), 3-se(东南), 4-sw(西南)

            ne = rectangle_t(x + 0.25*w, y + 0.25*h, 0.5*w, 0.5*h)
            nw = rectangle_t(x - 0.25*w, y + 0.25*h, 0.5*w, 0.5*h)
            se = rectangle_t(x + 0.25*w, y - 0.25*h, 0.5*w, 0.5*h)
            sw = rectangle_t(x - 0.25*w, y - 0.25*h, 0.5*w, 0.5*h)

            call self%children(1)%constructor(ne, self%capacity)
            call self%children(2)%constructor(nw, self%capacity)
            call self%children(3)%constructor(se, self%capacity)
            call self%children(4)%constructor(sw, self%capacity)

            self%divided = .true.

        end associate

    end subroutine subdivide

    !> 打印四叉树
    recursive subroutine show(self, indent)

        class(quad_tree_t), intent(inout) :: self
        integer, intent(in) :: indent

        integer :: i

        do i = 1, size(self%points)
            if (allocated(self%points)) then
                write(*, "(A,1x,A)", advance="no") repeat("·", indent), self%code
                write(*, "(2(2x,f6.2))") self%points(i)%x, self%points(i)%y
            end if
        end do

        if (allocated(self%children)) then

            do i = 1, size(self%children)
                call self%children(i)%show(indent+1)
            end do

        end if


    end subroutine show

    !> 析构例程
    subroutine delete_tree(self)

        class(quad_tree_t), intent(out) :: self

    end subroutine delete_tree

    !> 粒子查询操作
    ! subroutine query(self, range)

    !     class(quad_tree_t), intent(inout) :: self
    !     type(rectangle_t), intent(in) :: range

    !     if (.not. self%boundary%intersects(range)) then
            
    !     else
    !         associate(p => self%points)
    !             do i = 1, size(p)
    !                 if (range%contains(p(i))) then
    !                     found%push(p(i))
    !                 end if
    !             end do
    !         end associate

    !         if (self%divided) then
    !             if (any(self%children(:)%query(range))) then

    !             end if
    !         end if
        
    !     end if
        

    ! end subroutine query

end module quad