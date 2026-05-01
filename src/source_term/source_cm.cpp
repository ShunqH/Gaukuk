// C++ headers
#include <cmath>            // sqrt(), abs()
#include <algorithm>        // max()  
#include <stdexcept>        // std::runtime_error()

// Gaukuk dependence
#include "source.hpp"
#include "../utils/utils.hpp" // Config 

namespace Gaukuk{
    
void SourceTerm::UpdateSource(TArray<Real>& cons, const Real t, const Real dt, const Grid& grid, const Domain& domain){
    if (hasConstAcc){
        ConstGravity(cons, dt, grid); 
    }
    if (hasPointG){
        PointGravity(cons, dt, grid, domain); 
    }else if (hasBinary){
        BinaryGravity(cons, t, dt, grid, domain); 
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
        ab_ = std::sqrt( (x-obj0.x)*(x-obj0.x) + (y-obj0.y)*(y-obj0.y) + (z-obj0.z)*(z-obj0.z) ); 
        Real gmTot = obj0.gm + gm; 
        omega_ = std::sqrt(gmTot/(ab_*ab_*ab_)); 
        Real r0 = gm * ab_ / gmTot; 
        Real r1 = obj0.gm * ab_ / gmTot; 

        theta1_ = std::atan2(y-obj0.y, x-obj0.x); 
        theta0_ = theta1_ + PI; 

        obj1.gm = gm;
        obj1.x = r1;
        obj1.y = 0;
        obj1.z = 0;
        obj1.vx = 0;
        obj1.vy = 0;
        obj1.vz = 0;
        obj1.rs = rs;

        obj0.x = -r0;
        obj0.y = 0;
        obj0.z = 0;
        obj0.vx = 0;
        obj0.vy = 0;
        obj0.vz = 0;

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