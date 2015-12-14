Fast multi-camera motion solvers

--
C++ implementation of solvers from the ICCV paper:

Jonathan Ventura, Clemens Arth and Vincent Lepetit.  An efficient minimal solution for multi-camera motion.  International Conference on Computer Vision (ICCV), 2015.  Santiago, Chile.

--
Note: the solvers take 6x6 w matrices as input.  Each w matrix is the matrix u*v^T where u and v are corresponding rays in Pluecker coordinates.
