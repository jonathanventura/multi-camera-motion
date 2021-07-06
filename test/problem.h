#pragma once

#include <Eigen/Core>

void compute_error( const Eigen::Matrix3d &R, const Eigen::Vector3d &t,
                   const Eigen::Matrix3d &Rsoln, const Eigen::Vector3d &tsoln,
                  double &rot_angle_err, double &trans_angle_err, double &trans_scale_err );

