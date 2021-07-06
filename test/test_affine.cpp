
#include <multi-camera-motion/relpose_generalized_affine.h>

#include "problem.h"
#include "affine_problem.h"

#include <iostream>

#include <cstdlib>
#include <ctime>

static Eigen::Vector2d project( const Eigen::Vector3d &x )
{
    return x.head(2)/x[2];
}

int main( int argc, char **argv )
{
    srand(time(NULL));

    const double trans_mag = 1.;
    const double angle = 5.*M_PI/180.;
    const double noise = 0;
    const double affine_noise = 0;

    AffineProblem prob;
    generateAffineProblem( trans_mag, angle, noise, affine_noise, prob );
    
    Eigen::MatrixXd x(2,6);
    Eigen::MatrixXd y(2,6);
    Eigen::MatrixXd cu(3,6);
    Eigen::MatrixXd cv(3,6);
    std::vector<Eigen::Matrix2d> Au(6);
    std::vector<Eigen::Matrix2d> Av(6);
    std::vector<Eigen::Matrix3d> Rsolns;
    std::vector<Eigen::Vector3d> tsolns;

    for ( int i = 0; i < 6; i++ )
    {
        x.col(i) = project(prob.ray_pairs[i].first.x);
        y.col(i) = project(prob.ray_pairs[i].second.x);
        cu.col(i) = prob.ray_pairs[i].first.c;
        cv.col(i) = prob.ray_pairs[i].second.c;
        Au[i] = prob.ray_pairs[i].first.A;
        Av[i] = prob.ray_pairs[i].second.A;
    }
    
    relpose_generalized_affine(x,y,cu,cv,Au,Av,Rsolns,tsolns);
    
    for ( int i = 0; i < Rsolns.size(); i++ )
    {
        double rot_angle_err, trans_angle_err, trans_scale_err;
        compute_error( prob.R, prob.t,
                       Rsolns[i], tsolns[i],
                       rot_angle_err, trans_angle_err, trans_scale_err );
        std::cout << "solution " << i+1 << " error: ";
        std::cout << "\trotation error: " << rot_angle_err*180./M_PI << " deg\n";
        std::cout << "\ttranslation error: " << trans_angle_err*180./M_PI << " deg\n";
        std::cout << "\ttranslation scale error: " << trans_scale_err << "\n";
    }
        
}
