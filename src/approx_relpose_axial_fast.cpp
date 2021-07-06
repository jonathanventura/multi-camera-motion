
#include <multi-camera-motion/approx_relpose_axial_fast.h>

#include <Polynomial/Polynomial.hpp>
using polynomial::Polynomial;

#include <Eigen/LU>
#include <Eigen/Jacobi>
#include <Eigen/SVD>

#include "approx_relpose_axial_fast_computeA.h"

#include <iostream>

void
approx_relpose_axial_fast
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
    const Eigen::Matrix<double,15,35> A = approx_relpose_axial_fast_computeA(w1,w2,w3,w4,w5,w6);
    
    const Eigen::Matrix<double,15,20> G = A.block<15,15>(0,0).lu().solve(A.block<15,20>(0,15));

    const double B11coeffs[] = {-G(6,0), G(5,0) - G(6,1), G(5,1) - G(6,2), G(5,2)};
    const double B12coeffs[] = {-G(6,3), G(5,3) - G(6,4), G(5,4) - G(6,5), G(5,5)};
    const double B13coeffs[] = {-G(6,6), G(5,6) - G(6,7), G(5,7) - G(6,8), G(5,8)};
    const double B14coeffs[] = {-G(6,9), G(5,9) - G(6,10), G(5,10) - G(6,11), G(5,11) - G(6,12), G(5,12)};
    const double B15coeffs[] = {-G(6,13), G(5,13) - G(6,14), G(5,14) - G(6,15), G(5,15) - G(6,16), G(5,16)};
    const double B16coeffs[] = {-G(6,17), G(5,17) - G(6,18), G(5,18) - G(6,19), G(5,19)};
    const double B21coeffs[] = {-G(8,0), G(7,0) - G(8,1), G(7,1) - G(8,2), G(7,2)};
    const double B22coeffs[] = {-G(8,3), G(7,3) - G(8,4), G(7,4) - G(8,5), G(7,5)};
    const double B23coeffs[] = {-G(8,6), G(7,6) - G(8,7), G(7,7) - G(8,8), G(7,8)};
    const double B24coeffs[] = {-G(8,9), G(7,9) - G(8,10), G(7,10) - G(8,11), G(7,11) - G(8,12), G(7,12)};
    const double B25coeffs[] = {-G(8,13), G(7,13) - G(8,14), G(7,14) - G(8,15), G(7,15) - G(8,16), G(7,16)};
    const double B26coeffs[] = {-G(8,17), G(7,17) - G(8,18), G(7,18) - G(8,19), G(7,19)};
    const double B31coeffs[] = {-G(10,0), G(9,0) - G(10,1), G(9,1) - G(10,2), G(9,2)};
    const double B32coeffs[] = {-G(10,3), G(9,3) - G(10,4), G(9,4) - G(10,5), G(9,5)};
    const double B33coeffs[] = {-G(10,6), G(9,6) - G(10,7), G(9,7) - G(10,8), G(9,8)};
    const double B34coeffs[] = {-G(10,9), G(9,9) - G(10,10), G(9,10) - G(10,11), G(9,11) - G(10,12), G(9,12)};
    const double B35coeffs[] = {-G(10,13), G(9,13) - G(10,14), G(9,14) - G(10,15), G(9,15) - G(10,16), G(9,16)};
    const double B36coeffs[] = {-G(10,17), G(9,17) - G(10,18), G(9,18) - G(10,19), G(9,19)};
    const double B41coeffs[] = {-G(12,0), G(11,0) - G(12,1), G(11,1) - G(12,2), G(11,2)};
    const double B42coeffs[] = {-G(12,3), G(11,3) - G(12,4), G(11,4) - G(12,5), G(11,5)};
    const double B43coeffs[] = {-G(12,6), G(11,6) - G(12,7), G(11,7) - G(12,8), G(11,8)};
    const double B44coeffs[] = {-G(12,9), G(11,9) - G(12,10), G(11,10) - G(12,11), G(11,11) - G(12,12), G(11,12)};
    const double B45coeffs[] = {-G(12,13), G(11,13) - G(12,14), G(11,14) - G(12,15), G(11,15) - G(12,16), G(11,16)};
    const double B46coeffs[] = {-G(12,17), G(11,17) - G(12,18), G(11,18) - G(12,19), G(11,19)};
    const double B51coeffs[] = {G(13,0), G(13,1), G(13,2)};
    const double B52coeffs[] = {G(13,3), G(13,4), G(13,5)};
    const double B53coeffs[] = {G(13,6), G(13,7), G(13,8)};
    const double B54coeffs[] = {G(13,9), G(13,10), G(13,11), G(13,12)};
    const double B55coeffs[] = {G(13,13), G(13,14), G(13,15), G(13,16)};
    const double B56coeffs[] = {1, 0, G(13,17), G(13,18), G(13,19)};
    const double B61coeffs[] = {G(14,0), G(14,1), G(14,2)};
    const double B62coeffs[] = {G(14,3), G(14,4), G(14,5)};
    const double B63coeffs[] = {G(14,6), G(14,7), G(14,8)};
    const double B64coeffs[] = {G(14,9), G(14,10), G(14,11), G(14,12)};
    const double B65coeffs[] = {G(14,13), G(14,14), G(14,15), G(14,16)};
    const double B66coeffs[] = {1, G(14,17), G(14,18), G(14,19)};
    const Polynomial<3> B11(B11coeffs);
    const Polynomial<3> B12(B12coeffs);
    const Polynomial<3> B13(B13coeffs);
    const Polynomial<4> B14(B14coeffs);
    const Polynomial<4> B15(B15coeffs);
    const Polynomial<3> B16(B16coeffs);
    const Polynomial<3> B21(B21coeffs);
    const Polynomial<3> B22(B22coeffs);
    const Polynomial<3> B23(B23coeffs);
    const Polynomial<4> B24(B24coeffs);
    const Polynomial<4> B25(B25coeffs);
    const Polynomial<3> B26(B26coeffs);
    const Polynomial<3> B31(B31coeffs);
    const Polynomial<3> B32(B32coeffs);
    const Polynomial<3> B33(B33coeffs);
    const Polynomial<4> B34(B34coeffs);
    const Polynomial<4> B35(B35coeffs);
    const Polynomial<3> B36(B36coeffs);
    const Polynomial<3> B41(B41coeffs);
    const Polynomial<3> B42(B42coeffs);
    const Polynomial<3> B43(B43coeffs);
    const Polynomial<4> B44(B44coeffs);
    const Polynomial<4> B45(B45coeffs);
    const Polynomial<3> B46(B46coeffs);
    const Polynomial<2> B51(B51coeffs);
    const Polynomial<2> B52(B52coeffs);
    const Polynomial<2> B53(B53coeffs);
    const Polynomial<3> B54(B54coeffs);
    const Polynomial<3> B55(B55coeffs);
    const Polynomial<4> B56(B56coeffs);
    const Polynomial<2> B61(B61coeffs);
    const Polynomial<2> B62(B62coeffs);
    const Polynomial<2> B63(B63coeffs);
    const Polynomial<3> B64(B64coeffs);
    const Polynomial<3> B65(B65coeffs);
    const Polynomial<3> B66(B66coeffs);
    const Polynomial<6> t2 = B11*B22;
    const Polynomial<6> t14 = B12*B21;
    const Polynomial<6> t3 = t2 - t14;
    const Polynomial<6> t4 = B11*B32;
    const Polynomial<6> t16 = B12*B31;
    const Polynomial<6> t5 = t4 - t16;
    const Polynomial<6> t6 = B11*B42;
    const Polynomial<6> t27 = B12*B41;
    const Polynomial<6> t7 = t6 - t27;
    const Polynomial<6> t8 = B21*B32;
    const Polynomial<6> t17 = B22*B31;
    const Polynomial<6> t9 = t8 - t17;
    const Polynomial<6> t10 = B21*B42;
    const Polynomial<6> t28 = B22*B41;
    const Polynomial<6> t11 = t10 - t28;
    const Polynomial<6> t12 = B31*B42;
    const Polynomial<6> t39 = B32*B41;
    const Polynomial<6> t13 = t12 - t39;
    const Polynomial<9> t15 = B33*t3;
    const Polynomial<9> t18 = B13*t9;
    const Polynomial<9> t62 = B23*t5;
    const Polynomial<9> t19 = t15 + t18 - t62;
    const Polynomial<5> t20 = B11*B52;
    const Polynomial<5> t32 = B12*B51;
    const Polynomial<5> t21 = t20 - t32;
    const Polynomial<5> t22 = B21*B52;
    const Polynomial<5> t33 = B22*B51;
    const Polynomial<5> t23 = t22 - t33;
    const Polynomial<5> t24 = B31*B52;
    const Polynomial<5> t43 = B32*B51;
    const Polynomial<5> t25 = t24 - t43;
    const Polynomial<9> t26 = B43*t3;
    const Polynomial<9> t29 = B13*t11;
    const Polynomial<9> t129 = B23*t7;
    const Polynomial<9> t30 = t26 + t29 - t129;
    const Polynomial<8> t31 = B53*t3;
    const Polynomial<8> t34 = B13*t23;
    const Polynomial<8> t63 = B23*t21;
    const Polynomial<8> t35 = t31 + t34 - t63;
    const Polynomial<5> t36 = B41*B52;
    const Polynomial<5> t47 = B42*B51;
    const Polynomial<5> t37 = t36 - t47;
    const Polynomial<9> t38 = B43*t5;
    const Polynomial<9> t40 = B13*t13;
    const Polynomial<9> t99 = B33*t7;
    const Polynomial<9> t41 = t38 + t40 - t99;
    const Polynomial<8> t42 = B53*t5;
    const Polynomial<8> t44 = B13*t25;
    const Polynomial<8> t65 = B33*t21;
    const Polynomial<8> t45 = t42 + t44 - t65;
    const Polynomial<8> t46 = B53*t7;
    const Polynomial<8> t48 = B13*t37;
    const Polynomial<8> t101 = B43*t21;
    const Polynomial<8> t49 = t46 + t48 - t101;
    const Polynomial<9> t50 = B43*t9;
    const Polynomial<9> t51 = B23*t13;
    const Polynomial<9> t131 = B33*t11;
    const Polynomial<9> t52 = t50 + t51 - t131;
    const Polynomial<8> t53 = B53*t9;
    const Polynomial<8> t54 = B23*t25;
    const Polynomial<8> t66 = B33*t23;
    const Polynomial<8> t55 = t53 + t54 - t66;
    const Polynomial<8> t56 = B53*t11;
    const Polynomial<8> t57 = B23*t37;
    const Polynomial<8> t146 = B43*t23;
    const Polynomial<8> t58 = t56 + t57 - t146;
    const Polynomial<8> t59 = B53*t13;
    const Polynomial<8> t60 = B33*t37;
    const Polynomial<8> t102 = B43*t25;
    const Polynomial<8> t61 = t59 + t60 - t102;
    const Polynomial<12> t64 = B34*t35;
    const Polynomial<12> t67 = B14*t55;
    const Polynomial<12> t68 = t64 + t67 - B24*t45 - B54*t19;
    const Polynomial<5> t69 = B11*B62;
    const Polynomial<5> t76 = B12*B61;
    const Polynomial<5> t70 = t69 - t76;
    const Polynomial<5> t71 = B21*B62;
    const Polynomial<5> t77 = B22*B61;
    const Polynomial<5> t72 = t71 - t77;
    const Polynomial<5> t73 = B31*B62;
    const Polynomial<5> t83 = B32*B61;
    const Polynomial<5> t74 = t73 - t83;
    const Polynomial<8> t75 = B63*t3;
    const Polynomial<8> t78 = B13*t72;
    const Polynomial<8> t124 = B23*t70;
    const Polynomial<8> t79 = t75 + t78 - t124;
    const Polynomial<4> t80 = B51*B62;
    const Polynomial<4> t87 = B52*B61;
    const Polynomial<4> t81 = t80 - t87;
    const Polynomial<8> t82 = B63*t5;
    const Polynomial<8> t84 = B13*t74;
    const Polynomial<8> t105 = B33*t70;
    const Polynomial<8> t85 = t82 + t84 - t105;
    const Polynomial<7> t86 = B63*t21;
    const Polynomial<7> t88 = B13*t81;
    const Polynomial<7> t109 = B53*t70;
    const Polynomial<7> t89 = t86 + t88 - t109;
    const Polynomial<8> t90 = B63*t9;
    const Polynomial<8> t91 = B23*t74;
    const Polynomial<8> t126 = B33*t72;
    const Polynomial<8> t92 = t90 + t91 - t126;
    const Polynomial<7> t93 = B63*t23;
    const Polynomial<7> t94 = B23*t81;
    const Polynomial<7> t150 = B53*t72;
    const Polynomial<7> t95 = t93 + t94 - t150;
    const Polynomial<7> t96 = B63*t25;
    const Polynomial<7> t97 = B33*t81;
    const Polynomial<7> t111 = B53*t74;
    const Polynomial<7> t98 = t96 + t97 - t111;
    const Polynomial<12> t100 = B44*t45;
    const Polynomial<12> t103 = B14*t61;
    const Polynomial<12> t104 = t100 + t103 - B34*t49 - B54*t41;
    const Polynomial<5> t106 = B41*B62;
    const Polynomial<5> t114 = B42*B61;
    const Polynomial<5> t107 = t106 - t114;
    const Polynomial<11> t108 = B64*t45;
    const Polynomial<11> t110 = B34*t89;
    const Polynomial<11> t112 = t108 + t110 - B14*t98 - B54*t85;
    const Polynomial<8> t113 = B63*t7;
    const Polynomial<8> t115 = B13*t107;
    const Polynomial<8> t133 = B43*t70;
    const Polynomial<8> t116 = t113 + t115 - t133;
    const Polynomial<8> t117 = B63*t13;
    const Polynomial<8> t118 = B33*t107;
    const Polynomial<8> t136 = B43*t74;
    const Polynomial<8> t119 = t117 + t118 - t136;
    const Polynomial<7> t120 = B63*t37;
    const Polynomial<7> t121 = B43*t81;
    const Polynomial<7> t154 = B53*t107;
    const Polynomial<7> t122 = t120 + t121 - t154;
    const Polynomial<12> t123 = B64*t19;
    const Polynomial<12> t125 = B24*t85;
    const Polynomial<12> t127 = t123 + t125 - B14*t92 - B34*t79;
    const Polynomial<13> t128 = B44*t19;
    const Polynomial<13> t130 = B24*t41;
    const Polynomial<13> t132 = t128 + t130 - B34*t30 - B14*t52;
    const Polynomial<12> t134 = B64*t41;
    const Polynomial<12> t135 = B34*t116;
    const Polynomial<12> t137 = t134 + t135 - B44*t85 - B14*t119;
    const Polynomial<8> t138 = B63*t11;
    const Polynomial<8> t139 = B23*t107;
    const Polynomial<8> t143 = B43*t72;
    const Polynomial<8> t140 = t138 + t139 - t143;
    const Polynomial<12> t141 = B64*t30;
    const Polynomial<12> t142 = B24*t116;
    const Polynomial<12> t144 = t141 + t142 - B44*t79 - B14*t140;
    const Polynomial<12> t145 = B44*t35;
    const Polynomial<12> t147 = B14*t58;
    const Polynomial<11> t148 = B64*t35;
    const Polynomial<11> t149 = B24*t89;
    const Polynomial<11> t151 = t148 + t149 - B14*t95 - B54*t79;
    const Polynomial<11> t152 = B64*t49;
    const Polynomial<11> t153 = B44*t89;
    const Polynomial<11> t155 = t152 + t153 - B14*t122 - B54*t116;
    const Polynomial<12> t156 = B64*t52;
    const Polynomial<12> t157 = B34*t140;
    const Polynomial<12> t158 = t156 + t157 - B44*t92 - B24*t119;
    const Polynomial<12> t159 = B44*t55;
    const Polynomial<12> t160 = B24*t61;
    const Polynomial<11> t161 = B64*t55;
    const Polynomial<11> t162 = B34*t95;
    const Polynomial<11> t163 = t161 + t162 - B24*t98 - B54*t92;
    const Polynomial<11> t164 = B64*t58;
    const Polynomial<11> t165 = B44*t95;
    const Polynomial<11> t166 = t164 + t165 - B24*t122 - B54*t140;
    const Polynomial<11> t167 = B64*t61;
    const Polynomial<11> t168 = B44*t98;
    const Polynomial<11> t169 = t167 + t168 - B34*t122 - B54*t119;
    const Polynomial<20> detB = B66*(B45*t68 - B15*(t159 + t160 - B34*t58 - B54*t52) - B35*(t145 + t147 - B24*t49 - B54*t30) + B25*t104 + B55*t132) + B36*(B65*(t145 + t147 - B24*t49 - B54*t30) + B25*t155 - B15*t166 - B45*t151 + B55*t144) + B16*(B65*(t159 + t160 - B34*t58 - B54*t52) - B25*t169 + B35*t166 - B45*t163 + B55*t158) - B46*(B65*t68 + B25*t112 - B15*t163 + B55*t127 - B35*t151) - B26*(B65*t104 - B45*t112 - B15*t169 + B35*t155 + B55*t137) - B56*(B15*t158 - B45*t127 - B25*t137 + B35*t144 + B65*t132);

    std::vector<double> zsolns;
    detB.realRootsSturm(-15.*M_PI/180.,15.*M_PI/180.,zsolns);
    
    rsolns.clear();
    rsolns.reserve(zsolns.size());
    for ( size_t i = 0; i < zsolns.size(); i++ )
    {
        Eigen::Matrix<double,6,6> Bz;
        Bz <<
        B11.eval(zsolns[i]),B12.eval(zsolns[i]),B13.eval(zsolns[i]),B14.eval(zsolns[i]),B15.eval(zsolns[i]),B16.eval(zsolns[i]),
        B21.eval(zsolns[i]),B22.eval(zsolns[i]),B23.eval(zsolns[i]),B24.eval(zsolns[i]),B25.eval(zsolns[i]),B26.eval(zsolns[i]),
        B31.eval(zsolns[i]),B32.eval(zsolns[i]),B33.eval(zsolns[i]),B34.eval(zsolns[i]),B35.eval(zsolns[i]),B36.eval(zsolns[i]),
        B41.eval(zsolns[i]),B42.eval(zsolns[i]),B43.eval(zsolns[i]),B44.eval(zsolns[i]),B45.eval(zsolns[i]),B46.eval(zsolns[i]),
        B51.eval(zsolns[i]),B52.eval(zsolns[i]),B53.eval(zsolns[i]),B54.eval(zsolns[i]),B55.eval(zsolns[i]),B56.eval(zsolns[i]),
        B61.eval(zsolns[i]),B62.eval(zsolns[i]),B63.eval(zsolns[i]),B64.eval(zsolns[i]),B65.eval(zsolns[i]),B66.eval(zsolns[i]);
        Eigen::Matrix<double,6,1> xysoln;
        xysoln << 0, 0, 0, 0, 0, 1;
        xysoln = Bz.householderQr().solve(xysoln);
        const double xsoln = xysoln(3)/xysoln(5);
        const double ysoln = xysoln(4)/xysoln(5);
        Eigen::Matrix<double,3,1> rsoln;
        rsoln << xsoln, ysoln, zsolns[i];
        rsolns.push_back(rsoln);
    }

}

