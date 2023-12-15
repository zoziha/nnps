<div align='center'>

# 近邻粒子搜索算法

![NNPS](https://img.shields.io/badge/NNPS-v1.3.0-blueviolet)
![语言](https://img.shields.io/badge/-Fortran-734f96?logo=fortran&logoColor=white)
[![许可证](https://img.shields.io/badge/License-MIT-pink)](LICENSE)

NNPS 方案在 CFD 粒子法中的实践。

</div>

设计思想：
1. 面向对象;
2. 简单易用;
3. 高效灵活;
4. 仅实施必要的计算, 不引入不必要的CPU计算换时间。

TODO:
1. 背景网格结构回收规则。

## 使用方法

仅支持 [FPM]((https://github.com/fortran-lang/fpm))/Meson/Visual-Studio、
其他编译系统可直接复制源文件、
并测试了 `ifort/ifx` 和 `gfortran` 编译器。

要在`fpm`项目中使用`nnps`，请在`fpm.toml`文件中添加以下几行：

```sh
[dependencies］
nnps = { git = "https://gitee.com/zoziha/nnps" }
```

**注意：Windows-ifort bug，使用 OpenMP 静态编译会导致可分配数组出现问题！请使用 `/libs:dll` 编译选项！**
或者取消 `fpm.toml` 的 OpenMP 依赖，使用串行版本：
```sh
fpm run --flag '/DSERIAL' --profile release --example --all --compiler ifort
```

### 并行

对于 2D/3D NNPS，OpenMP 用于并行加速，并行线程可通过 `OMP_NUM_THREADS` 环境变量设置。

### 示例

```sh
fpm run --example --all # 运行所有示例
```

```fortran
program example_grid2d

    use nnps_module, only: nnps_grid2d, wp
    use display_module, only: display
    implicit none

    type(nnps_grid2d) :: nnps
    real(wp), dimension(2, 4) :: loc = reshape([0.0_wp, 1.0_wp, 2.0_wp, 1.5_wp, &
                                                1.0_wp, 1.0_wp, 0.5_wp, 1.0_wp], [2, 4])
    integer, pointer :: pairs(:)
    real(wp), pointer :: rdxs(:)

    call nnps%init(loc, n=4)
    call nnps%query(0.6_wp, pairs, rdxs, n=4)

    print *, '*** grid find (2D)'
    call display(pairs, 'pairs index:', brief=.false.)
    call display(loc(:, pairs), 'pairs coordinates:', brief=.false.)
    call display(rdxs, 'rdxs:', brief=.false.)

end program example_grid2d
!  *** grid find (2D)
! [vector: 4] pairs index:
! 3, 4, 4, 1
! [matrix: 2*4] pairs coordinates:
!  1.000E+00,  5.000E-01,  5.000E-01,  0.000E+00;
!  1.000E+00,  1.000E+00,  1.000E+00,  1.000E+00
! [vector: 6] rdxs:
!  5.000E-01,  5.000E-01,  0.000E+00,  5.000E-01,  5.000E-01,  0.000E+00
```

这个库最初是研究 SPH 粒子法的 NNPS 算法，现在已经引入到 SPH 代码中，参见 [zoziha/SPH](https://github.com/zoziha/SPH)。

## 链接

+ [wikipedia/Quadtree](https://en.wikipedia.org/wiki/Quadtree)
+ [lewisfish/quad-tree](https://github.com/lewisfish/quad-trees)
+ [bilibili/Dejavu32/[Coding Challenge ]:四叉树的实现 以及 碰撞检测应用（Quadtree Collisions )](https://www.bilibili.com/video/BV1ub411S7N5?spm_id_from=333.999.0.0)
+ [CodingTrain/QuadTree](https://github.com/CodingTrain/QuadTree)