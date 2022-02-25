# 四叉树模型

四叉树模型在碰撞检测中有很大的应用。

在粒子法流体力学中，时间复杂度为$o(N\lg{N})$将使得粒子搜索非常高效和健壮。

四叉树也是八叉树的基础，有望在后续的工作中添加对八叉树的支持。

| 项目 | 描述 |
| :-: | :-: |
| 版本 | 0.1.0 |
| 许可证 | MIT |
| 版权 | Copyright (c) 2021~2022 quad-tree 贡献者 |

## 开始 [English](README_EN.md)

### 获取代码

```sh
git clone https://github.com/zoziha/quad-tree.git
cd quad-tree
```

### 使用[fortran-lang/fpm](https://github.com/fortran-lang/fpm)构建

FPM是社区驱动的Fortran语言的包管理器和代码构建工具。  
你可以通过提供的`fpm.toml`构建代码：

```sh
fpm build
fpm run --example --list
fpm run --example <example_name>
```

原本本库是作为 SPH 粒子法的 NNPS 算法的学习，现如今已经被引入到了 SPH 代码中，详见：[zoziha/SPH](https://github.com/zoziha/SPH)。

## 链接

+ [wikipedia/Quadtree](https://en.wikipedia.org/wiki/Quadtree)
+ [lewisfish/quad-tree](https://github.com/lewisfish/quad-trees)
+ [bilibili/Dejavu32/[Coding Challenge ]:四叉树的实现 以及 碰撞检测应用（Quadtree Collisions )](https://www.bilibili.com/video/BV1ub411S7N5?spm_id_from=333.999.0.0)
+ [CodingTrain/QuadTree](https://github.com/CodingTrain/QuadTree)