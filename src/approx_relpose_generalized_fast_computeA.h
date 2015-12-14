#ifndef APPROX_RELPOSE_GENERALIZED_FAST_COMPUTEA_H
#define APPROX_RELPOSE_GENERALIZED_FAST_COMPUTEA_H

#include <Eigen/Core>

Eigen::Matrix<double,15,35> computeA(
    const Eigen::Matrix<double,6,6> &w1,
    const Eigen::Matrix<double,6,6> &w2,
    const Eigen::Matrix<double,6,6> &w3,
    const Eigen::Matrix<double,6,6> &w4,
    const Eigen::Matrix<double,6,6> &w5,
    const Eigen::Matrix<double,6,6> &w6
);

#endif
