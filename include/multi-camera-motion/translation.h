#pragma once

#include <Eigen/Core>

// solve for translation from rotation and observations
Eigen::Vector3d solve_translation(
 const Eigen::MatrixXd &x,    // 2xN points in camera in first rig
 const Eigen::MatrixXd &y,    // 2xN points in camera in second rig
 const Eigen::MatrixXd &cu,   // 3xN camera centers in first rig
 const Eigen::MatrixXd &cv,   // 3xN camera centers in second rig
 const Eigen::Matrix3d &R     // rotation
);
