#ifndef APPROX_RELPOSE_AXIAL_H
#define APPROX_RELPOSE_AXIAL_H

#include <Eigen/Core>
#include <vector>

void
approx_relpose_axial
(
    const Eigen::Matrix<double,6,6> &w1,
    const Eigen::Matrix<double,6,6> &w2,
    const Eigen::Matrix<double,6,6> &w3,
    const Eigen::Matrix<double,6,6> &w4,
    const Eigen::Matrix<double,6,6> &w5,
    const Eigen::Matrix<double,6,6> &w6,
    std::vector<Eigen::Vector3d> &rsolns
);

#endif
