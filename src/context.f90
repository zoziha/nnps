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

        integer :: i

        boundary = rectangle_t(100.0, 100.0, 200.0, 200.0)
        call qt%constructor(boundary, 4)

        !> 布置粒子：随机位置
        do i = 1, 6
            p = point_t(randu(start=100.0, end=200.0), randu(start=100.0, end=200.0))
            print *, qt%insert(p, "0")
        end do

        call qt%show(1)

    end subroutine setup

end module context