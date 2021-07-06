
#include <Eigen/Core>
#include "so3.h"

void compute_error( const Eigen::Matrix3d &R, const Eigen::Vector3d &t,
                   const Eigen::Matrix3d &Rsoln, const Eigen::Vector3d &tsoln,
                  double &rot_angle_err, double &trans_angle_err, double &trans_scale_err )
{
    Eigen::Matrix3d Rdiff = R*Rsoln.transpose();
    Eigen::Vector3d rdiff = so3ln(Rdiff);
    rot_angle_err = rdiff.norm();

    trans_angle_err = acos(t.dot(tsoln)/t.norm()/tsoln.norm());

    if ( isnan(trans_angle_err) ) trans_angle_err = 0;
    trans_scale_err = tsoln.norm()/t.norm();
}

