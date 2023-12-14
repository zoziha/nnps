---
title: 监测器
---

# 监测器

本部分将就内存管理设计进行介绍, NNPS 涉及的内存控制主要分为两个部分:

1. 基于粒子总数的现实交互粒子对数量: 如果粒子模拟正常, 
   粒子的交互数量将持续保持为粒子数的确定倍数 (与初始粒子分布有关), 而不是发生发散或者坍塌。
2. 基于粒子总数的现实路径网格存储: 如果粒子持续扩散其背景路径域, 则路径网格存储将持续增加, 对此, 需要评估惰性网格的回收。

三种方法中, 固定光滑长度方案建议采用背景网格法, 变光滑长度方案建议采用树型搜索法。

@todo 如何评价变光滑长度方案的应用场景与在特定场景中的实用价值。