#pragma once

#include <Eigen/Core>

static inline double randf() { return ((double)rand()/RAND_MAX)*2.-1.; }
static inline Eigen::Vector3d rand3f() { Eigen::Vector3d v; v << randf(),randf(),randf(); return v; }
static inline Eigen::Vector3d offset(Eigen::Vector3d g, double n, double m) { Eigen::Vector3d v; v << g[0]+n,g[1]+m,g[2]; return v; }
static inline
double grandf() {

    double r,v1,v2,fac;

    r=2;
    while (r>=1) {
        v1=(2*((double)rand()/(double)RAND_MAX)-1);
        v2=(2*((double)rand()/(double)RAND_MAX)-1);
        r=v1*v1+v2*v2;
    }
    fac=sqrt(-2*log(r)/r);

    return(v2*fac);

}
static inline Eigen::Vector2d grand2f() { Eigen::Vector2d v; v << grandf(),grandf(); return v; }

