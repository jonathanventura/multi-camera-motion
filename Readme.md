### Multi-camera Motion

C++ implementation of solvers from these papers:

Jonathan Ventura, Clemens Arth and Vincent Lepetit.  An efficient minimal solution for multi-camera motion.  International Conference on Computer Vision (ICCV), 2015.  Santiago, Chile.

(Note: these solvers take 6x6 w matrices as input.  Each w matrix is the matrix u*v^T where u and v are corresponding rays in Pluecker coordinates. )

Khaled Alyousefi and Jonathan Ventura.  Multi-camera motion estimation with affine correspondences.  International Conference on Image Analysis and Recognition. Springer, Cham, 2020.

#### Dependencies:

Eigen

[Polynomial](https://github.com/jonathanventura/polynomial) (my polynomial root-finding library)

---

This material is based upon work supported by the National Science Foundation under Grant Number 1464420.
