# N-Tree

A collection of commonly used functions for bin-tree/quad-tree/oct-tree for Fortran.

*Suggestions and code contributions are welcome.*

```fortran
use ntree_factory_m, only: ntree_t, shape_t, point_t, &
                           make_ntree, make_boundary, make_range
```

## Build with [fortran-lang/fpm](https://github.com/fortran-lang/fpm)

```sh
fpm run --example --list
```

```toml
[dependencies]
ntree = { git = "https://gitee.com/zoziha/ntree.git" }
```

Originally a study of the NNPS algorithm for SPH particle methods, this library has now been introduced into the SPH code, see [zoziha/SPH](https://github.com/zoziha/SPH).

## Links

+ [wikipedia/Quadtree](https://en.wikipedia.org/wiki/Quadtree)
+ [lewisfish/quad-tree](https://github.com/lewisfish/quad-trees)
+ [bilibili/Dejavu32/[Coding Challenge ]:四叉树的实现 以及 碰撞检测应用（Quadtree Collisions )](https://www.bilibili.com/video/BV1ub411S7N5?spm_id_from=333.999.0.0)
+ [CodingTrain/QuadTree](https://github.com/CodingTrain/QuadTree)
