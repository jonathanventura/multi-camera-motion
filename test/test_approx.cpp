
#include <multi-camera-motion/approx_relpose_generalized.h>
#include <multi-camera-motion/translation.h>

#include <Eigen/Geometry>

#include "problem.h"
#include "so3.h"

#include <iostream>

#include <cstdlib>
#include <ctime>

static Eigen::Vector2d project( const Eigen::Vector3d &x )
{
    return x.head(2)/x[2];
}

static Eigen::Matrix<double,6,1> pluecker( const Ray & ray )
{
    // make ray in Pluecker coordinates
    Eigen::Vector3d d = ray.x/ray.x.norm();
    Eigen::Vector3d m = ray.c.cross(d);
    Eigen::Matrix<double,6,1> p;
    p.head(3) = d;
    p.tail(3) = m;
    return p;
}

static Eigen::Matrix<double,6,6> make_w( const RayPair &ray_pair )
{
    Eigen::Matrix<double,6,1> u = pluecker(ray_pair.first);
    Eigen::Matrix<double,6,1> v = pluecker(ray_pair.second);
    return u*v.transpose();
}

int main( int argc, char **argv )
{
    srand(time(NULL));

    const double trans_mag = 1.;
    const double angle = 5.*M_PI/180.;
    const double noise = 0;

    Problem prob;
    generateProblem( trans_mag, angle, noise, prob );

    Eigen::MatrixXd x(2,6);
    Eigen::MatrixXd y(2,6);
    Eigen::MatrixXd cu(3,6);
    Eigen::MatrixXd cv(3,6);
    for ( int i = 0; i < 6; i++ )
    {
        x.col(i) = project(prob.ray_pairs[i].first.x);
        y.col(i) = project(prob.ray_pairs[i].second.x);
        cu.col(i) = prob.ray_pairs[i].first.c;
        cv.col(i) = prob.ray_pairs[i].second.c;
    }
    
    Eigen::Matrix<double,6,6> w1;
    Eigen::Matrix<double,6,6> w2;
    Eigen::Matrix<double,6,6> w3;
    Eigen::Matrix<double,6,6> w4;
    Eigen::Matrix<double,6,6> w5;
    Eigen::Matrix<double,6,6> w6;
    std::vector<Eigen::Vector3d> rsolns;

    w1 = make_w(prob.ray_pairs[0]);
    w2 = make_w(prob.ray_pairs[1]);
    w3 = make_w(prob.ray_pairs[2]);
    w4 = make_w(prob.ray_pairs[3]);
    w5 = make_w(prob.ray_pairs[4]);
    w6 = make_w(prob.ray_pairs[5]);
    approx_relpose_generalized(w1,w2,w3,w4,w5,w6,rsolns);
    
    for ( int i = 0; i < rsolns.size(); i++ )
    {
        Eigen::Matrix3d Rsoln = so3exp(rsolns[i]);
        Eigen::Vector3d tsoln = solve_translation(x,y,cu,cv,Rsoln);
        double rot_angle_err, trans_angle_err, trans_scale_err;
        compute_error( prob.R, prob.t,
                       Rsoln, tsoln,
                       rot_angle_err, trans_angle_err, trans_scale_err );
        std::cout << "solution " << i+1 << " error: ";
        std::cout << "\trotation error: " << rot_angle_err*180./M_PI << " deg\n";
        std::cout << "\ttranslation error: " << trans_angle_err*180./M_PI << " deg\n";
        std::cout << "\ttranslation scale error: " << trans_scale_err << "\n";
    }
        
}
