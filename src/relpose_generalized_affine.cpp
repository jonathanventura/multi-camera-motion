
#include <multi-camera-motion/relpose_generalized_affine.h>

#include <multi-camera-motion/translation.h>

#include <Eigen/LU>
#include <Eigen/Jacobi>
#include <Eigen/SVD>

static Eigen::Matrix<double,1,9> get_AE_row( const Eigen::Vector2d &x,
                                             const Eigen::Vector2d &y,
                                             const Eigen::Vector3d &cu,
                                             const Eigen::Vector3d &cv )
{
    Eigen::Matrix<double,1,9> AE_row;
    AE_row << x(0)*y(0), x(0)*y(1), x(0), x(1)*y(0), x(1)*y(1), x(1), y(0), y(1), 1;
    return AE_row;
}

static Eigen::Matrix<double,1,9> get_AR_row( const Eigen::Vector2d &x,
                                             const Eigen::Vector2d &y,
                                             const Eigen::Vector3d &cu,
                                             const Eigen::Vector3d &cv )
{
    Eigen::Matrix<double,1,9> AR_row;
    const double t6 = cv(2)*y(1);
    const double t2 = cv(1) - t6;
    const double t7 = cv(2)*y(0);
    const double t3 = cv(0) - t7;
    const double t4 = cv(0)*y(1);
    const double t8 = cv(1)*y(0);
    const double t5 = t4 - t8;
    AR_row << cu(1)*y(0) + t2*x(0) - cu(2)*x(1)*y(0), cu(1)*y(1) - t3*x(0) - cu(2)*x(1)*y(1), cu(1) - cu(2)*x(1) + t5*x(0), t2*x(1) - cu(0)*y(0) + cu(2)*x(0)*y(0), cu(2)*x(0)*y(1) - t3*x(1) - cu(0)*y(1), cu(2)*x(0) - cu(0) + t5*x(1), cv(1) - t6 + cu(0)*x(1)*y(0) - cu(1)*x(0)*y(0), t7 - cv(0) + cu(0)*x(1)*y(1) - cu(1)*x(0)*y(1), t4 - t8 + cu(0)*x(1) - cu(1)*x(0);
    return AR_row;
}
static Eigen::Matrix<double,3,9> get_AE_rows( const Eigen::Vector2d &x,
                                            const Eigen::Vector2d &y,
                                            const Eigen::Vector3d &cu,
                                            const Eigen::Vector3d &cv,
                                            double a1,
                                            double a2,
                                            double a3,
                                            double a4 )
{
    Eigen::Matrix<double,3,9> AE_rows;
    AE_rows <<
    x(0)*y(0), x(0)*y(1), x(0), x(1)*y(0), x(1)*y(1), x(1), y(0), y(1), 1,
    a3*x(0), a4*x(0), 0, y(0) + a3*x(1), y(1) + a4*x(1), 1, a3, a4, 0,
    y(0) + a1*x(0), y(1) + a2*x(0), 1, a1*x(1), a2*x(1), 0, a1, a2, 0;
    return AE_rows;
}

static Eigen::Matrix<double,3,9> get_AR_rows( const Eigen::Vector2d &x,
                     const Eigen::Vector2d &y,
                     const Eigen::Vector3d &cu,
                     const Eigen::Vector3d &cv,
                     double a1,
                     double a2,
                     double a3,
                     double a4 )
{
    Eigen::Matrix<double,3,9> AR_rows;
    const double t2 = a4*x(1);
    const double t3 = t2 + y(1);
    const double t4 = a3*x(1);
    const double t5 = t4 + y(0);
    const double t6 = a1*x(0);
    const double t7 = t6 + y(0);
    const double t8 = a2*x(0);
    const double t9 = t8 + y(1);
    AR_rows << cv(1)*x(0) + cu(1)*y(0) - cu(2)*x(1)*y(0) - cv(2)*x(0)*y(1), cu(1)*y(1) - cv(0)*x(0) - cu(2)*x(1)*y(1) + cv(2)*x(0)*y(0), cu(1) - cu(2)*x(1) + cv(0)*x(0)*y(1) - cv(1)*x(0)*y(0), cv(1)*x(1) - cu(0)*y(0) + cu(2)*x(0)*y(0) - cv(2)*x(1)*y(1), cu(2)*x(0)*y(1) - cu(0)*y(1) - cv(0)*x(1) + cv(2)*x(1)*y(0), cu(2)*x(0) - cu(0) + cv(0)*x(1)*y(1) - cv(1)*x(1)*y(0), cv(1) - cv(2)*y(1) + cu(0)*x(1)*y(0) - cu(1)*x(0)*y(0), cv(2)*y(0) - cv(0) + cu(0)*x(1)*y(1) - cu(1)*x(0)*y(1), cu(0)*x(1) - cu(1)*x(0) + cv(0)*y(1) - cv(1)*y(0),
    a3*cu(1) - cu(2)*t5 - a4*cv(2)*x(0), a4*cu(1) - cu(2)*t3 + a3*cv(2)*x(0), a4*cv(0)*x(0) - a3*cv(1)*x(0) - cu(2), cv(1) - a3*cu(0) - cv(2)*t3 + a3*cu(2)*x(0), cv(2)*t5 - a4*cu(0) - cv(0) + a4*cu(2)*x(0), cv(0)*t3 - cv(1)*t5, cu(0)*t5 - a4*cv(2) - a3*cu(1)*x(0), a3*cv(2) + cu(0)*t3 - a4*cu(1)*x(0), cu(0) - a3*cv(1) + a4*cv(0),
    cv(1) + a1*cu(1) - cv(2)*t9 - a1*cu(2)*x(1), a2*cu(1) - cv(0) + cv(2)*t7 - a2*cu(2)*x(1), cv(0)*t9 - cv(1)*t7, cu(2)*t7 - a1*cu(0) - a2*cv(2)*x(1), cu(2)*t9 - a2*cu(0) + a1*cv(2)*x(1), cu(2) - a1*cv(1)*x(1) + a2*cv(0)*x(1), a1*cu(0)*x(1) - cu(1)*t7 - a2*cv(2), a1*cv(2) - cu(1)*t9 + a2*cu(0)*x(1), a2*cv(0) - a1*cv(1) - cu(1);
    return AR_rows;
}

void
relpose_generalized_affine
(
 const Eigen::MatrixXd &x,    // 2xN
 const Eigen::MatrixXd &y,    // 2xN
 const Eigen::MatrixXd &cu,   // 3xN
 const Eigen::MatrixXd &cv,   // 3xN
 const std::vector<Eigen::Matrix2d> &Au,
 const std::vector<Eigen::Matrix2d> &Av,
 std::vector<Eigen::Matrix3d> &Rsolns,
 std::vector<Eigen::Vector3d> &tsolns
)
{
    const int N = x.cols();
    
    Eigen::VectorXd a1(N);
    Eigen::VectorXd a2(N);
    Eigen::VectorXd a3(N);
    Eigen::VectorXd a4(N);
    
    Eigen::MatrixXd AE(3*N,9);
    Eigen::MatrixXd AR(3*N,9);
    
    for ( int i = 0; i < N; i++ )
    {
        Eigen::Matrix2d A = Av[i] * Au[i].inverse();
        a1(i) = A(0,0);
        a2(i) = A(1,0);
        a3(i) = A(0,1);
        a4(i) = A(1,1);
    }
    
    for ( int i = 0; i < N; i++ )
    {
        AE.block<3,9>(i*3,0) = get_AE_rows(x.col(i),y.col(i),cu.col(i),cv.col(i),a1(i),a2(i),a3(i),a4(i));
        AR.block<3,9>(i*3,0) = get_AR_rows(x.col(i),y.col(i),cu.col(i),cv.col(i),a1(i),a2(i),a3(i),a4(i));
    }

    // see Eigen FAQ about pseudo-inverse
    Eigen::JacobiSVD< Eigen::MatrixXd > ARsvd(AR,Eigen::ComputeThinU|Eigen::ComputeThinV);
    Eigen::VectorXd singularValues = ARsvd.singularValues();
    Eigen::VectorXd singularValues_inv = singularValues;

    // see Wikipedia article on pseudo-inverse
    double pinvtoler = Eigen::NumTraits<double>::epsilon() * (3*N) * singularValues(0);
    for ( long i=0; i < 9; ++i ) {
        if ( singularValues(i) > pinvtoler )
            singularValues_inv(i)=1.0/singularValues(i);
    }

    Eigen::MatrixXd ARp = (ARsvd.matrixV() * singularValues_inv.asDiagonal()) * ARsvd.matrixU().transpose();
    
    Eigen::MatrixXd AE2 = (AR*ARp-Eigen::MatrixXd::Identity(3*N,3*N))*AE;
    
    Eigen::Matrix<double,9,1> Evec = AE2.jacobiSvd(Eigen::ComputeThinV).matrixV().col(8);
    Eigen::Matrix3d E;
    E <<
    Evec[0], Evec[3], Evec[6],
    Evec[1], Evec[4], Evec[7],
    Evec[2], Evec[5], Evec[8];
    
    Eigen::JacobiSVD<Eigen::Matrix3d> svdE(E,Eigen::ComputeFullU|Eigen::ComputeFullV);
    
    Eigen::Matrix3d U = svdE.matrixU();
    Eigen::Matrix3d V = svdE.matrixV();
    
    if ( U.determinant() < 0 ) U.col(2) *= -1.0;
    if ( V.determinant() < 0 ) V.col(2) *= -1.0;
    
    Eigen::Matrix3d D;
    D <<
    0,1,0,
    -1,0,0,
    0,0,1;
    
    Eigen::Matrix3d DT;
    DT <<
    0,-1,0,
    1,0,0,
    0,0,1;
    
    Eigen::Matrix3d VT = V.transpose().eval();
    
    Eigen::Matrix3d R1 = U*D*VT;
    Eigen::Matrix3d R2 = U*DT*VT;
    
    Rsolns.resize(2);
    tsolns.resize(2);
    Rsolns[0] = R1;
    Rsolns[1] = R2;
    
    tsolns[0] = solve_translation( x, y, cu, cv, Rsolns[0] );
    tsolns[1] = solve_translation( x, y, cu, cv, Rsolns[1] );
}

