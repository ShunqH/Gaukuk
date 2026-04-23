// C++ headers
#include <cmath>            // sqrt(), abs()
#include <algorithm>        // max()  
#include <stdexcept>        // std::runtime_error()

// Gaukuk dependence
#include "source.hpp"
#include "../utils/utils.hpp" // Config 

namespace Gaukuk{
    
void SourceTerm::UpdateSource(TArray<Real>& cons, const Real dt, const Grid& grid, const Domain& domain){
    if (hasConstAcc){
        ConstGravity(cons, dt, grid); 
    }else if (hasPointG){
        PointGravity(cons, dt, grid, domain); 
    }
}

void SourceTerm::EnrollConstGravityVector(Real gx, Real gy, Real gz){
    gx_ = gx; 
    gy_ = gy; 
    gz_ = gz; 
    if (gx_!=0 || gy_!=0 || gz_!=0){
        sourceEnrolled = true; 
        hasConstAcc = true;  
    } 
} 

void SourceTerm::EnrollPointGravity(Real gm, Real x, Real y, Real z, 
                                    Real vx, Real vy, Real vz, Real rs){
    if (hasBinary){
        throw std::runtime_error("Setup failed: too many point gravity.");
    }
    if (hasPointG) {
        obj1.gm = gm;
        obj1.x = x;
        obj1.y = y;
        obj1.z = z;
        obj1.vx = vx;
        obj1.vy = vy;
        obj1.vz = vz;
        obj1.rs = rs;
        hasBinary = true;
        hasPointG = false;
        sourceEnrolled = true; 
    } else {
        obj0.gm = gm;
        obj0.x = x;
        obj0.y = y;
        obj0.z = z;
        obj0.vx = vx;
        obj0.vy = vy;
        obj0.vz = vz;
        obj0.rs = rs;
        hasPointG = true;
        sourceEnrolled = true; 
    }
}

}