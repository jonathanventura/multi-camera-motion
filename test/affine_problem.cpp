
#include <Eigen/Core>
#include <Eigen/Geometry>

#include <opencv2/calib3d.hpp>

#include <vector>

#include "random.h"
#include "so3.h"

#include "affine_problem.h"

struct AffineRegion
{
    Eigen::Vector3d X;
    Eigen::Vector3d u1;
    Eigen::Vector3d u2;
    Eigen::Vector3d u3;
    Eigen::Vector3d u4;
};

Eigen::Matrix2d getRelativeAffineMatrix(AffineRegion Au, AffineRegion Av,
                                       Eigen::Vector3d cu, Eigen::Vector3d cv)
{
    // calculate affine matrix from u image to v image
    // solves A*(u-center_u) = (v-center_v)
    // where v is projected point in v image, center_v is projection of PX in v image
    // where u is projected point in u image, center_u is projection of X in u image
    Eigen::Vector3d center_u3 = Au.X-cu;
    Eigen::Vector2d center_u(center_u3[0]/center_u3[2],center_u3[1]/center_u3[2]);

    Eigen::Vector3d u13 = Au.u1 - cu;
    Eigen::Vector3d u23 = Au.u2 - cu;
    Eigen::Vector3d u33 = Au.u3 - cu;
    Eigen::Vector3d u43 = Au.u4 - cu;

    Eigen::Vector2d u1(u13[0]/u13[2],u13[1]/u13[2]);
    Eigen::Vector2d u2(u23[0]/u23[2],u23[1]/u23[2]);
    Eigen::Vector2d u3(u33[0]/u33[2],u33[1]/u33[2]);
    Eigen::Vector2d u4(u43[0]/u43[2],u43[1]/u43[2]);

    Eigen::Vector3d center_v3 = Av.X-cv;
    Eigen::Vector2d center_v(center_v3[0]/center_v3[2],center_v3[1]/center_v3[2]);

    Eigen::Vector3d v13 = Av.u1 - cv;
    Eigen::Vector3d v23 = Av.u2 - cv;
    Eigen::Vector3d v33 = Av.u3 - cv;
    Eigen::Vector3d v43 = Av.u4 - cv;

    Eigen::Vector2d v1(v13[0]/v13[2],v13[1]/v13[2]);
    Eigen::Vector2d v2(v23[0]/v23[2],v23[1]/v23[2]);
    Eigen::Vector2d v3(v33[0]/v33[2],v33[1]/v33[2]);
    Eigen::Vector2d v4(v43[0]/v43[2],v43[1]/v43[2]);

    std::vector<cv::Point2f> srcPoints(4);
    std::vector<cv::Point2f> dstPoints(4);
    srcPoints[0] = cv::Point2f(u1(0),u1(1));
    srcPoints[1] = cv::Point2f(u2(0),u2(1));
    srcPoints[2] = cv::Point2f(u3(0),u3(1));
    srcPoints[3] = cv::Point2f(u4(0),u4(1));
    dstPoints[0] = cv::Point2f(v1(0),v1(1));
    dstPoints[1] = cv::Point2f(v2(0),v2(1));
    dstPoints[2] = cv::Point2f(v3(0),v3(1));
    dstPoints[3] = cv::Point2f(v4(0),v4(1));
    cv::Mat Hmat = cv::findHomography(srcPoints, dstPoints);
    Eigen::MatrixXf H(3,3);
    for ( int i = 0; i < 3; i++ )
        for ( int j = 0; j < 3; j++ )
            H(i,j) = Hmat.at<double>(i,j);

    float s = H(2,0) * center_v(0) + H(2,1) * center_v(1) + H(2,2) ;

    Eigen::Matrix2d affine_matrix;

    affine_matrix << (H(0,0) - H(2,0) * v1(0))/s ,
                     (H(0,1) - H(2,1) * v1(0))/s,
                     (H(1,0) - H(2,0) * v1(1))/s,
                     (H(1,1) - H(2,1) * v1(1))/s;

    return affine_matrix;
}


Eigen::Vector3d add_noise(Eigen::Vector3d p, Eigen::Vector2d noise)
{
    double unorm = p.norm();
    Eigen::Vector3d uproj = p/p(2);
    uproj.head(2) += noise;
    uproj /= uproj.norm();
    uproj *= unorm;
    if ( p.dot(uproj) < 0 ) uproj = -uproj;

    return uproj;
}

Eigen::Vector2d generate_noise(double noise)
{
    return grand2f()*noise;
}

AffineRegion generateAffineRegion(Eigen::Vector3d X)
{
    const double R = -.01;
    AffineRegion a;
    a.X = X;
    a.u1 = offset(X,-R,-R);
    a.u2 = offset(X,R,-R);
    a.u3 = offset(X,-R,R);
    a.u4 = offset(X,R,R);
    return a;
}

void generateAffineProblem( double trans_mag, double angle, double noise, double affine_noise, AffineProblem &prob, int nrays )
{
    Eigen::Vector3d X;    // 3D point
    Eigen::Vector3d cu;   // camera center for first multi-camera rig
    Eigen::Vector3d cv;   // camera center for second multi-camera rig
    Eigen::Vector3d PX;   // 3D point after R,t transformation
    AffineRay u;          // observations in first multi-camera rig (Pluecker line + affine)
    AffineRay v;          // observations in second multi-camera rig

    AffineRegion A1;      // affine region in first camera
    AffineRegion A2;      // affine region in second camera

    Eigen::Vector3d w = rand3f();  // make random axis for rotation
    w = w/w.norm();                // normalize rotation axis
    double theta = angle*M_PI/180; // make random angle for rotation
    prob.R = so3exp(w*theta);                // random rotation
    prob.t = rand3f();                       // random translation
    while ( prob.t.norm() < 1e-10 ) prob.t = rand3f();  // ensure translation is not all zeros

    prob.ray_pairs.clear(); // clear out ray pairs to be safe

    // generate observations
    for ( size_t i = 0; i < nrays; i++ )
    {
        X = rand3f();                          // sample random point
        while ( X.norm() == 0 ) X = rand3f();  // ensure point is not all zeros
        X(2) = 4.;
        X = X/X.norm();                        // normalize point
        double scale_factor = ((double)rand()/RAND_MAX) * 4 + 4;
        X *= scale_factor;                     // stretch point out to distance in range [4,8]

        A1 = generateAffineRegion(X);     // generate affine region

        cu = rand3f();        // sample random camera locations
        cv = rand3f();        // inter-camera
        //cv = cu;            // intra-camera: camera center is same in both rigs

        PX = prob.R*X+prob.t;    // transform point to second rig's frame
        A2.X = PX;               // transform affine region
        A2.u1 = prob.R*A1.u1+prob.t;
        A2.u2 = prob.R*A1.u2+prob.t;
        A2.u3 = prob.R*A1.u3+prob.t;
        A2.u4 = prob.R*A1.u4+prob.t;

        // affine matrix for first rig is identity
        u.x = X-cu;
        u.c = cu;
        u.A << 1, 0,
               0, 1;

        // calculate relative affine matrix for second rig
        v.x = PX-cv;
        v.c = cv;
        v.A = getRelativeAffineMatrix(A1, A2, cu, cv);

        if ( noise > 0. )    // add noise with requested standard deviation
        {
            // add gaussian noise to center points
            u.x = add_noise(u.x, generate_noise(noise));
            v.x = add_noise(v.x, generate_noise(noise));
        }

        if ( affine_noise > 0. )
        {
            // add gaussian noise to affine parameters
            v.A(0,0) += grandf()*affine_noise;
            v.A(0,1) += grandf()*affine_noise;
            v.A(1,0) += grandf()*affine_noise;
            v.A(1,1) += grandf()*affine_noise;
        }

        // store observation
        prob.ray_pairs.push_back( AffineRayPair( u, v ) );
    }
}


