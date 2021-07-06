
#pragma once

#include <Eigen/Core>
#include <vector>

struct AffineRay {
    // Image location in camera coordinates
    Eigen::Vector3d x;
    // Camera center
    Eigen::Vector3d c;
    // Affine frame
    Eigen::Matrix2d A;
};
typedef std::pair<AffineRay,AffineRay> AffineRayPair;
typedef std::vector<AffineRayPair> AffineRayPairList;

struct AffineProblem
{
    Eigen::Matrix3d R;
    Eigen::Vector3d t;
    AffineRayPairList ray_pairs;
};

void generateAffineProblem( double trans_mag, double angle, double noise, double affine_noise, AffineProblem &prob, int nrays = 18 );

