// C++ headers
#include <cmath>            // sqrt(), abs()
#include <algorithm>        // max()  
#include <iostream>     // std::cout; std::endl; std::cerr

// Gaukuk dependence
#include "source.hpp"

namespace Gaukuk{
    
void SourceTerm::PointGravity(TArray<Real>& cons, const Real dt, const Grid& grid, const Domain& domain){
    int il = grid.ib;               // first activated cell left side
    int ir = grid.ie;               // last activated cell right side + 1
    int jl = grid.jb;               // first activated cell left side
    int jr = grid.je;               // last activated cell right side + 1
    int kl = grid.kb;               // first activated cell left side
    int kr = grid.ke;               // last activated cell right side + 1

    Real gm = obj0.gm; 
    Real x0 = obj0.x; 
    Real y0 = obj0.y; 
    Real z0 = obj0.z; 
    Real rs = obj0.rs; 
    Real gmdt = gm*dt; 

    #pragma omp parallel for collapse(2) schedule(static)
    for (int k=kl; k<kr; k++){
        for (int j=jl; j<jr; j++){
    #pragma omp simd 
            for (int i=il; i<ir; i++){
                Real xnow = domain.xc(i); 
                Real ynow = domain.yc(j); 
                Real znow = domain.zc(k); 

                Real delx = xnow - x0; 
                Real dely = ynow - y0; 
                Real delz = znow - z0; 

                Real r2 = delx*delx + dely*dely + delz*delz; 
                Real rSmoInv = 1.0 / std::sqrt(r2 + rs*rs);
                Real gVecdtdr = -gmdt * rSmoInv * rSmoInv * rSmoInv;

                Real gxdt = gVecdtdr * delx; 
                Real gydt = gVecdtdr * dely; 
                Real gzdt = gVecdtdr * delz; 

                Real& consDen = cons(DEN, k, j, i); 
                Real& consMtx = cons(MTX, k, j, i); 
                Real& consMty = cons(MTY, k, j, i); 
                Real& consMtz = cons(MTZ, k, j, i); 
                Real& consEng = cons(ENG, k, j, i); 

                // probably I don't need rk2?
                // ------------------- apply rk2 -------------------
                // Real den0 = consDen;
                // Real mtx0 = consMtx; 
                // Real mty0 = consMty; 
                // Real mtz0 = consMtz; 
                // Real eng0 = consEng; 

                // // half
                // const Real gxdtHalf = gxdt;
                // const Real gydtHalf = gydt;
                // const Real gzdtHalf = gzdt;

                // // second order mid point
                // Real mtxHalf = mtx0 + den0 * gxdtHalf * 0.5;
                // Real mtyHalf = mty0 + den0 * gydtHalf * 0.5;
                // Real mtzHalf = mtz0 + den0 * gzdtHalf * 0.5;

                // // correction
                // consEng = eng0 + (mtxHalf * gxdtHalf + mtyHalf * gydtHalf + mtzHalf * gzdtHalf);
                // consMtx = mtx0 + den0 * gxdtHalf;
                // consMty = mty0 + den0 * gydtHalf;
                // consMtz = mtz0 + den0 * gzdtHalf;
                // -------------------- rk2 end -------------------
  
                consEng += consMtx*gxdt + consMty*gydt + consMtz*gzdt;
                consMtx += consDen*gxdt;
                consMty += consDen*gydt;
                consMtz += consDen*gzdt; 

                // accreation (sink) 
                // Real factor = (r2 < rs*rs) ? 0 : 1; 
                // consEng = (consEng + consMtx*gxdt + consMty*gydt + consMtz*gzdt) * factor; 
                // consMtx = ( consMtx + consDen*gxdt ) * factor;
                // consMty = ( consMty + consDen*gydt ) * factor;
                // consMtz = ( consMtz + consDen*gzdt ) * factor; 
                // consDen *= factor;  
            }
        }
    }
}

}