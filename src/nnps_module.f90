!> Complete interfaces
module nnps_module

    use nnps_kinds, only: wp
    use nnps_direct1d_module, only: nnps_direct1d
    use nnps_direct2d_module, only: nnps_direct2d
    use nnps_direct3d_module, only: nnps_direct3d
    use nnps_tree1d_module, only: nnps_binarytree
    use nnps_tree2d_module, only: nnps_quadtree
    use nnps_tree3d_module, only: nnps_octree
    use nnps_grid1d_module, only: nnps_grid1d
    use nnps_grid2d_module, only: nnps_grid2d, nnps_grid2d_finalizer
    use nnps_grid3d_module, only: nnps_grid3d, nnps_grid3d_finalizer
    implicit none

end module nnps_module
