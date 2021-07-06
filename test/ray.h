#ifndef RAY_H
#define RAY_H

#include <Eigen/Core>
#include <vector>
#include <algorithm>

typedef Eigen::Matrix<double,6,1> Ray;
typedef std::pair<Ray,Ray> RayPair;
typedef std::vector<RayPair> RayPairList;


#endif

