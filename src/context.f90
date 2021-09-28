module context

    use quad, only: quad_tree_t
    use utils, only: rectangle_t, point_t
    use forlab_stats, only: randu
    implicit none
    
contains

    !> 建立场景，构建四叉树。
    subroutine setup()

        type(rectangle_t) :: boundary       !! 目标域
        type(quad_tree_t) :: qt             !! 四叉树
        type(point_t)     :: p              !! 粒子

        type(rectangle_t) :: range          !! 查找域
        type(point_t), allocatable :: found(:)  !! 查找所得粒子

        integer :: i
        logical :: bool

        boundary = rectangle_t(100.0, 100.0, 200.0, 200.0)
        call qt%constructor(boundary, 4)

        ! range = rectangle_t(100.0, 100.0, 200.0, 200.0)
        range = rectangle_t(randu(50.0, 100.0), randu(50.0, 100.0), randu(50.0, 100.0), randu(50.0, 100.0))
        print *, "查找域信息: "//new_line(""), range%x, range%y, range%w, range%h

        !> 布置粒子：随机位置
        do i = 1, 16
            p    = point_t(randu(start=0.0, end=200.0), randu(start=0.0, end=200.0))
            bool = qt%insert(p, code="0")
        end do

        !> 打印四叉树
        print *, "四叉树结构信息: "
        call qt%show(indent=1)

        !> 查询粒子
        call qt%query(range, found)

        !> 打印粒子
        print *, "查询查找域内的粒子(坐标): "
        do i = 1, size(found)
            print "(2(f6.2,2x))", found(i)%x, found(i)%y
        end do
        print *, "查询到粒子数: ", i-1

    end subroutine setup

end module context