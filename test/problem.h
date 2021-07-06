#pragma once

#include <Eigen/Core>
#include <vector>
#include <algorithm>

void compute_error( const Eigen::Matrix3d &R, const Eigen::Vector3d &t,
                   const Eigen::Matrix3d &Rsoln, const Eigen::Vector3d &tsoln,
                  double &rot_angle_err, double &trans_angle_err, double &trans_scale_err );

struct Ray {
    // Image location in camera coordinates
    Eigen::Vector3d x;
    // Camera center
    Eigen::Vector3d c;
};
typedef std::pair<Ray,Ray> RayPair;
typedef std::vector<RayPair> RayPairList;

struct Problem
{
    Eigen::Matrix3d R;
    Eigen::Vector3d t;
    RayPairList ray_pairs;
};

void generateProblem( double trans_mag, double angle, double noise, Problem &prob, int nrays = 18 );

