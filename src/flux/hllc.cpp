// C++ headers 
#include <cmath>            // sqrt
#include <algorithm>        // max 

// Gaukuk dependence
#include "../gaukuk.hpp"    // Real
#include "flux.hpp"         // Sim 
#include "../eos/eos.hpp"   // EquationOfState

namespace Gaukuk{

/* 
---------------------------------------------------------------
    HLLC Riemann Solver 
    Toro E.F.
    Riemann Solvers and Numerical Methods for Fluid Dynamics 
    Chapter 10.4 ~ 10.6
----------------------------------------------------------------
*/

void Flux::RiemannSolver(const TArray<Real>& ul, const TArray<Real>& ur, const int direction, 
                         EquationOfState& eos, TArray<Real>& flux, 
                         int k, int j, int igb, int ige){
    const int IVLX = VLX + (direction - VLX + 0) % 3; 
    const int IVLY = VLX + (direction - VLX + 1) % 3; 
    const int IVLZ = VLX + (direction - VLX + 2) % 3; 
    Real gamma = eos.GetGamma(); 
    Real gmRec = 1/gamma; 

#pragma omp simd
    for (int i=igb; i<ige; i++){
        // load data for SIMD optimization
        Real denl = ul(DEN, i); 
        Real vxl = ul(IVLX, i); 
        Real vyl = ul(IVLY, i); 
        Real vzl = ul(IVLZ, i); 
        Real prel = ul(PRE, i);

        Real denr = ur(DEN, i); 
        Real vxr = ur(IVLX, i); 
        Real vyr = ur(IVLY, i); 
        Real vzr = ur(IVLZ, i); 
        Real prer = ur(PRE, i);

        // Step I pressure estimate. 
        Real al = eos.SoundSpeed(denl, prel); 
        Real ar = eos.SoundSpeed(denr, prer); 
        Real prePVRS = 0.5*(prel + prer) - 0.125*(vxr - vxl)*(denl + denr)*(al + ar); 

        // Step II wave speed estimates. 
        // ql = (prePVRS<=prel) ? 1 : std::sqrt( 1 + 0.5*(1 + gmRec) * (prePVRS/prel - 1) ); 
        // qr = (prePVRS<=prer) ? 1 : std::sqrt( 1 + 0.5*(1 + gmRec) * (prePVRS/prer - 1) ); 
        Real temp = std::max(prePVRS/prel-1, Real(0)); 
        Real ql = std::sqrt( 1 + 0.5*(1 + gmRec) * temp );
        temp = std::max(prePVRS/prer-1, Real(0)); 
        Real qr = std::sqrt( 1 + 0.5*(1 + gmRec) * temp );
        
        Real sl = vxl - al*ql; 
        Real sr = vxr = ar*qr; 
        Real ss = ( prer - prel + denl*vxl*(sl-vxl) - denr*vxr*(sr-vxr) ) / 
                  ( denl*(sl-vxl) - denr*(sr-vxr) ) ;
        
        // Step III HLLC flux
        Real el = eos.EGas(denl, prel) + 0.5*(vxl*vxl+vyl*vyl+vzl*vzl); 
        Real er = eos.EGas(denr, prer) + 0.5*(vxr*vxr+vyr*vyr+vzr*vzr); 

        Real tl = std::min(Real(0), sl)*(ss-vxl)/(sl-ss); 

        Real fl1 = denl*vxl;  
        Real fl2 = denl*vxl*vxl + prel; 
        Real fl3 = denl*vxl*vyl; 
        Real fl4 = denl*vxl*vzl; 
        Real fl5 = (el + prel)*vxl;

        Real fr1 = denr*vxr;  
        Real fr2 = denr*vxr*vxr + prer; 
        Real fr3 = denr*vxr*vyr; 
        Real fr4 = denr*vxr*vzr; 
        Real fr5 = (er + prer)*vxr;


    }

}

} // namespace Gaukuk