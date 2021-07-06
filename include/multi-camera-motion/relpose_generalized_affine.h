#pragma once 

#include <Eigen/Core>
#include <vector>

void
relpose_generalized_affine
(
 const Eigen::MatrixXd &x,    // 2xN points in cameras in first rig
 const Eigen::MatrixXd &y,    // 2xN points in cameras in second rig
 const Eigen::MatrixXd &cu,   // 3xN centers of cameras in first rig
 const Eigen::MatrixXd &cv,   // 3xN centers of cameras in second rig
 const std::vector<Eigen::Matrix2d> &Au,  // affine shapes in first rig
 const std::vector<Eigen::Matrix2d> &Av,  // affine shapes in second rig
 std::vector<Eigen::Matrix3d> &Rsolns,  // rotation solutions
 std::vector<Eigen::Vector3d> &tsolns  // translation solutions
);

