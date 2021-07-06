
#include <Eigen/Core>
#include "problem.h"
#include "random.h"
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

void generateProblem( double trans_mag, double angle, double noise, Problem &prob, int nrays )
{
    Eigen::Vector3d X;    // 3D points
    Eigen::Vector3d cu;   // camera centers for first multi-camera rig
    Eigen::Vector3d cv;   // camera centers for second multi-camera rig
    Eigen::Vector3d PX;   // 3D points after R,t transformation
    Ray u;                // observations in first multi-camera rig
    Ray v;                // observations in second multi-camera rig

    Eigen::Vector3d w = rand3f();
    w = w/w.norm();
    double theta = angle*M_PI/180;
    prob.R = so3exp(w*theta);                // random rotation
    prob.t = rand3f();                       // random translation
    while ( prob.t.norm() < 1e-10 ) prob.t = rand3f();
    prob.t = prob.t/prob.t.norm()*trans_mag;

    prob.ray_pairs.clear();

    for ( size_t i = 0; i < nrays; i++ ) {
        X = rand3f();    // sample random point
        while ( X.norm() == 0 ) X = rand3f();
        X = X/X.norm();
        double scale_factor = ((double)rand()/RAND_MAX) * 4 + 4;
        X *= scale_factor;

        cu = rand3f();        // sample random camera locations
        cv = rand3f(); // inter-camera
        //cv = cu;    // intra-camera

        PX = prob.R*X+prob.t;    // transform point
        u.x = X-cu;    // make observation rays
        u.c = cu;
        v.x = PX-cv;
        v.c = cv;

        if ( noise > 0. )    // add noise with requested standard deviation
        {
            double unorm = u.x.norm();
            Eigen::Vector3d uproj = u.x/u.x(2);
            uproj.head(2) += grand2f()*noise;
            uproj /= uproj.norm();
            uproj *= unorm;
            if ( (X-cu).dot(uproj) < 0 ) uproj = -uproj;
            u.x = uproj;

            double vnorm = v.x.norm();
            Eigen::Vector3d vproj = v.x/v.x(2);
            vproj.head(2) += grand2f()*noise;
            vproj /= vproj.norm();
            vproj *= vnorm;
            if ( (PX-cv).dot(vproj) < 0 ) vproj = -vproj;
            v.x = vproj;
        }

        prob.ray_pairs.push_back( RayPair( u, v ) );
    }
}

