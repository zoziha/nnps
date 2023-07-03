<div align='center'>

# Nearest Neighbor Particle Search

![NNPS](https://img.shields.io/badge/NNPS-v1.0.0-blueviolet)
![Language](https://img.shields.io/badge/-Fortran-734f96?logo=fortran&logoColor=white)
[![license](https://img.shields.io/badge/License-MIT-pink)](LICENSE)

NNPS scheme practice in CFD particle method.

</div>

## Usage

Only [FPM]((https://github.com/fortran-lang/fpm)) and Meson is supported,
other build systems can copy source files directly,
and `ifort/ifx` and `gfortran` compilers are tested.

To use `nnps` within your `fpm` project, add the following lines to your `fpm.toml` file:

```toml
[dependencies]
nnps = { git="https://github.com/zoziha/nnps" }
```

## Example

```sh
> fpm run --example --all  # run all examples
```

```fortran
program example_direct2d

    use nnps_module, only: nnps_direct2d, rk
    use display_module, only: display
    implicit none

    type(nnps_direct2d) :: nnps
    real(rk), dimension(2, 4) :: loc = reshape([0.0_rk, 1.0_rk, 2.0_rk, 1.5_rk, &
                                                1.0_rk, 1.0_rk, 0.5_rk, 1.0_rk], [2, 4])
    integer, pointer :: pairs(:)

    call nnps%init(loc)
    call nnps%query(0.6_rk, pairs)

    call display(real(pairs), 'pairs index:')
    call display(loc(:, pairs), 'pairs coordinates:')

end program example_direct2d
!> [vector: 4] pairs index:
!>  1.000E+00,  4.000E+00,  3.000E+00,  4.000E+00
!> [matrix: 2*4] pairs coordinates:
!>  0.000E+00,  5.000E-01,  1.000E+00,  5.000E-01;
!>  1.000E+00,  1.000E+00,  1.000E+00,  1.000E+00
```

Originally a study of the NNPS algorithm for SPH particle methods, this library has now been introduced into the SPH code, see [zoziha/SPH](https://github.com/zoziha/SPH).

## Links

+ [wikipedia/Quadtree](https://en.wikipedia.org/wiki/Quadtree)
+ [lewisfish/quad-tree](https://github.com/lewisfish/quad-trees)
+ [bilibili/Dejavu32/[Coding Challenge ]:四叉树的实现 以及 碰撞检测应用（Quadtree Collisions )](https://www.bilibili.com/video/BV1ub411S7N5?spm_id_from=333.999.0.0)
+ [CodingTrain/QuadTree](https://github.com/CodingTrain/QuadTree)
