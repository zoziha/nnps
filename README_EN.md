<div align='center'>

# Nearest Neighbor Particle Search

![NNPS](https://img.shields.io/badge/NNPS-v1.3.0-blueviolet)
![Language](https://img.shields.io/badge/-Fortran-734f96?logo=fortran&logoColor=white)
[![license](https://img.shields.io/badge/License-MIT-pink)](LICENSE)

NNPS scheme practice in CFD particle method.

</div>

## Usage

Only [FPM]((https://github.com/fortran-lang/fpm))/Meson/Visual-Studio are supported,
other build systems can copy source files directly,
and `ifort/ifx` and `gfortran` compilers are tested.

To use `nnps` within your `fpm` project, add the following lines to your `fpm.toml` file:

```toml
[dependencies]
nnps = { git="https://github.com/zoziha/nnps" }
```

**Note: Windows-ifort bug, Using static compilation will cause problems with allocatable arrays! Please use the `/libs:dll` compile option!**

### Parallel

For 2D/3D NNPS, OpenMP is used for parallel acceleration, and parallel threads can be set through the `OMP_NUM_THREADS` environment variable.

## Example

```sh
> fpm run --example --all  # run all examples
```

```fortran
program example_grid2d

    use nnps_module, only: nnps_grid2d, rk
    use display_module, only: display
    implicit none

    type(nnps_grid2d) :: nnps
    real(rk), dimension(2, 4) :: loc = reshape([0.0_rk, 1.0_rk, 2.0_rk, 1.5_rk, &
                                                1.0_rk, 1.0_rk, 0.5_rk, 1.0_rk], [2, 4])
    integer, pointer :: pairs(:)
    real(rk), pointer :: rdxs(:)

    call nnps%init(loc, n=4)
    call nnps%query(0.6_rk, pairs, rdxs, n=4)

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

Originally a study of the NNPS algorithm for SPH particle methods, this library has now been introduced into the SPH code, see [zoziha/SPH](https://github.com/zoziha/SPH).

## Links

+ [wikipedia/Quadtree](https://en.wikipedia.org/wiki/Quadtree)
+ [lewisfish/quad-tree](https://github.com/lewisfish/quad-trees)
+ [bilibili/Dejavu32/[Coding Challenge ]:四叉树的实现 以及 碰撞检测应用（Quadtree Collisions )](https://www.bilibili.com/video/BV1ub411S7N5?spm_id_from=333.999.0.0)
+ [CodingTrain/QuadTree](https://github.com/CodingTrain/QuadTree)
