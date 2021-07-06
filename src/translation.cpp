
#include <Eigen/Core>
#include <Eigen/SVD>
#include <Eigen/Jacobi>

static Eigen::Vector4d get_translation_row( const Eigen::Matrix3d &R,
                                           const Eigen::Vector2d &x,
                                           const Eigen::Vector2d &y,
                                           const Eigen::Vector3d &cu,
                                           const Eigen::Vector3d &cv
                                           )
{
    return Eigen::Vector4d( R(1,2) - R(2,2)*y(1) + x(0)*(R(1,0) - R(2,0)*y(1)) + x(1)*(R(1,1) - R(2,1)*y(1)), R(2,2)*y(0) - R(0,2) - x(0)*(R(0,0) - R(2,0)*y(0)) - x(1)*(R(0,1) - R(2,1)*y(0)), R(0,2)*y(1) - R(1,2)*y(0) + x(0)*(R(0,0)*y(1) - R(1,0)*y(0)) + x(1)*(R(0,1)*y(1) - R(1,1)*y(0)), R(2,0)*cu(1) - R(2,1)*cu(0) + R(0,2)*cv(1) - R(1,2)*cv(0) + y(0)*(R(0,0)*cu(1) - R(0,1)*cu(0) + R(1,2)*cv(2) - R(2,2)*cv(1)) + y(1)*(R(1,0)*cu(1) - R(1,1)*cu(0) - R(0,2)*cv(2) + R(2,2)*cv(0)) + x(0)*(R(2,1)*cu(2) - R(2,2)*cu(1) + R(0,0)*cv(1) - R(1,0)*cv(0) + y(0)*(R(0,1)*cu(2) - R(0,2)*cu(1) + R(1,0)*cv(2) - R(2,0)*cv(1)) + y(1)*(R(1,1)*cu(2) - R(1,2)*cu(1) - R(0,0)*cv(2) + R(2,0)*cv(0))) - x(1)*(R(2,0)*cu(2) - R(2,2)*cu(0) - R(0,1)*cv(1) + R(1,1)*cv(0) + y(0)*(R(0,0)*cu(2) - R(0,2)*cu(0) - R(1,1)*cv(2) + R(2,1)*cv(1)) + y(1)*(R(1,0)*cu(2) - R(1,2)*cu(0) + R(0,1)*cv(2) - R(2,1)*cv(0))) );
}

// solve for translation from rotation and observations
Eigen::Vector3d solve_translation(
 const Eigen::MatrixXd &x,    // 2xN points in camera in first rig
 const Eigen::MatrixXd &y,    // 2xN points in camera in second rig
 const Eigen::MatrixXd &cu,   // 3xN camera centers in first rig
 const Eigen::MatrixXd &cv,   // 3xN camera centers in second rig
 const Eigen::Matrix3d &R     // rotation
)
{
    const int N = x.cols();

    // solve for translations
    Eigen::MatrixXd A(N,4);

    for ( int j = 0; j < N; j++ )
    {
        // point
        A.row(j) = get_translation_row( R, x.col(j),  y.col(j), cu.col(j), cv.col(j) );
    }

    Eigen::JacobiSVD<Eigen::MatrixXd> jacobiSvd(A,Eigen::ComputeThinV);
    const Eigen::Matrix<double,4,1> t = jacobiSvd.matrixV().col(3);
    Eigen::Vector3d tsoln = t.head(3)/t[3];
    
    return tsoln;
}

