---
title: User Guide
---

This repository contains the source code for nearest neighbor particle search using the methods as below:

* direct search (\\(O(n^2/2)\\));
* kd-tree search (\\(O(nlogn)\\));
* grid search (\\(O(n)\\)).

For optimization purposes: Since distance solving will be involved in NNPS,
and distance solving is needed in other solving processes in CFD particle method,
in order to reduce the amount of computation,the distances obtained from solving will be stored
in the data structure of this warehouse in order to improve the efficiency as much as possible.
The distances stored are the line distance and the axis component of the line distance.

@todo
Adding OMP parallel computing to particle search for large datasets.

@note
The direct search method is relatively simple because it does not require the construction of data structures;
the tree search method may take a larger total time because of the large amount of computation for single-particle searches;
and the background grid method has a low time complexity and a relatively small amount of computation if the background grid is reused.
