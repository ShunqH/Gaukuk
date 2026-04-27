// C++ headers
#include <cmath>            // sqrt(), abs()
#include <algorithm>        // max()  
#include <iostream>     // std::cout; std::endl; std::cerr

// Gaukuk dependence
#include "source.hpp"

namespace Gaukuk{
    
void SourceTerm::BinaryGravity(TArray<Real>& cons, const Real t, const Real dt, const Grid& grid, const Domain& domain){
    int il = grid.ib;               // first activated cell left side
    int ir = grid.ie;               // last activated cell right side + 1
    int jl = grid.jb;               // first activated cell left side
    int jr = grid.je;               // last activated cell right side + 1
    int kl = grid.kb;               // first activated cell left side
    int kr = grid.ke;               // last activated cell right side + 1

    const Real gm0 = obj0.gm; 
    const Real rs0 = obj0.rs; 
    const Real gm0dt = gm0*dt; 

    const Real gm1 = obj1.gm; 
    const Real rs1 = obj1.rs; 
    const Real gm1dt = gm1*dt; 

    Real theta = omega_ * t + theta0_; 
    obj1.x = ab_ * std::cos(theta);
    obj1.y = ab_ * std::sin(theta);

    Real x0 = obj0.x; 
    Real y0 = obj0.y; 
    Real z0 = obj0.z; 

    Real x1 = obj1.x; 
    Real y1 = obj1.y; 
    Real z1 = obj1.z; 

    // indirect source trem since this is not a barycenter frame
    Real r1Prim = std::sqrt(x1*x1 + y1*y1 + z1*z1); 
    Real r1pInv = 1.0/r1Prim; 
    Real gIndirdtdr1 = - gm1dt * r1pInv * r1pInv * r1pInv; 

    #pragma omp parallel for collapse(2) schedule(static)
    for (int k=kl; k<kr; k++){
        for (int j=jl; j<jr; j++){
    #pragma omp simd 
            for (int i=il; i<ir; i++){
                Real xnow = domain.xc(i); 
                Real ynow = domain.yc(j); 
                Real znow = domain.zc(k); 

                Real delx0 = xnow - x0; 
                Real dely0 = ynow - y0; 
                Real delz0 = znow - z0; 
                Real delx1 = xnow - x1; 
                Real dely1 = ynow - y1; 
                Real delz1 = znow - z1; 

                Real r0Sqr = delx0*delx0 + dely0*dely0 + delz0*delz0; 
                Real r0SmoInv = 1.0 / std::sqrt(r0Sqr + rs0*rs0);
                Real g0Vecdtdr = - gm0dt * r0SmoInv * r0SmoInv * r0SmoInv;

                Real r1Sqr = delx1*delx1 + dely1*dely1 + delz1*delz1; 
                Real r1SmoInv = 1.0 / std::sqrt(r1Sqr + rs1*rs1);
                Real g1Vecdtdr = - gm1dt * r1SmoInv * r1SmoInv * r1SmoInv;

                Real gxdt = g0Vecdtdr * delx0 + g1Vecdtdr * delx1 + gIndirdtdr1 * x1; 
                Real gydt = g0Vecdtdr * dely0 + g1Vecdtdr * dely1 + gIndirdtdr1 * y1; 
                Real gzdt = g0Vecdtdr * delz0 + g1Vecdtdr * delz1 + gIndirdtdr1 * z1; 

                Real& consDen = cons(DEN, k, j, i); 
                Real& consMtx = cons(MTX, k, j, i); 
                Real& consMty = cons(MTY, k, j, i); 
                Real& consMtz = cons(MTZ, k, j, i); 
                Real& consEng = cons(ENG, k, j, i); 
                
                consEng += consMtx*gxdt + consMty*gydt + consMtz*gzdt;
                consMtx += consDen*gxdt;
                consMty += consDen*gydt;
                consMtz += consDen*gzdt; 
            }
        }
    }
}

}