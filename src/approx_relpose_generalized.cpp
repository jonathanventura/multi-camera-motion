
#include "approx_relpose_generalized.h"

#include <Eigen/LU>
#include <Eigen/Jacobi>
#include <Eigen/SVD>
#include <Eigen/Eigenvalues>

#include "approx_relpose_generalized_fast_computeA.h"

#include <iostream>

void
approx_relpose_generalized
(
    const Eigen::Matrix<double,6,6> &w1,
    const Eigen::Matrix<double,6,6> &w2,
    const Eigen::Matrix<double,6,6> &w3,
    const Eigen::Matrix<double,6,6> &w4,
    const Eigen::Matrix<double,6,6> &w5,
    const Eigen::Matrix<double,6,6> &w6,
    std::vector<Eigen::Vector3d> &rsolns
)
{

    const Eigen::Matrix<double,15,35> A = computeA(w1,w2,w3,w4,w5,w6);
    Eigen::Matrix<double,15,35> gbA;
    gbA << A.col(0),A.col(1),A.col(2),A.col(3),A.col(4),A.col(5),A.col(7),A.col(9),A.col(11),A.col(15),A.col(18),A.col(21),A.col(24),A.col(28),A.col(13),A.col(6),A.col(8),A.col(10),A.col(12),A.col(16),A.col(19),A.col(22),A.col(25),A.col(29),A.col(14),A.col(17),A.col(20),A.col(23),A.col(26),A.col(30),A.col(32),A.col(27),A.col(31),A.col(33),A.col(34);

    const Eigen::Matrix<double,15,20> G = gbA.block<15,15>(0,0).lu().solve(gbA.block<15,20>(0,15));

    Eigen::Matrix<double,20,20> M = Eigen::Matrix<double,20,20>::Zero();
    M.block<10,20>(0,0) = -G.block<10,20>(5,0);
    M(10,4) = 1;
    M(11,5) = 1;
    M(12,6) = 1;
    M(13,7) = 1;
    M(14,8) = 1;
    M(15,9) = 1;
    M(16,13) = 1;
    M(17,14) = 1;
    M(18,15) = 1;
    M(19,18) = 1;
    
    const Eigen::EigenSolver< Eigen::Matrix<double,20,20> > eigensolver(M,true);
    const Eigen::EigenSolver< Eigen::Matrix<double,20,20> >::EigenvalueType evalues = eigensolver.eigenvalues();
    const Eigen::EigenSolver< Eigen::Matrix<double,20,20> >::EigenvectorsType evecs = eigensolver.eigenvectors();

    rsolns.clear();
    rsolns.reserve(evalues.size());
    for ( size_t i = 0; i < evalues.size(); i++ )
    {
        if ( evalues[i].imag() != 0 ) continue;
        const double zsoln = evalues(i).real();
        const double xsoln = evecs(16,i).real()/evecs(19,i).real();
        const double ysoln = evecs(17,i).real()/evecs(19,i).real();
        Eigen::Vector3d rsoln;
        rsoln << xsoln, ysoln, zsoln;
        rsolns.push_back(rsoln);
    }

}

