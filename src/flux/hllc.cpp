// C++ headers 
#include <cmath>            // std::sqrt()
#include <algorithm>        // std::max() 
#include <iostream>     // std::cout; std::endl; std::cerr

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
    Real gmRec = 1.0/gamma; 

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
        prePVRS = std::max(prePVRS, Real(0));

        // Step II wave speed estimates. 
        // ql = (prePVRS<=prel) ? 1 : std::sqrt( 1 + 0.5*(1 + gmRec) * (prePVRS/prel - 1) ); 
        // qr = (prePVRS<=prer) ? 1 : std::sqrt( 1 + 0.5*(1 + gmRec) * (prePVRS/prer - 1) ); 
        Real temp = std::max(prePVRS/prel-1, Real(0)); 
        Real ql = std::sqrt( 1 + 0.5*(1 + gmRec) * temp );
        temp = std::max(prePVRS/prer-1, Real(0)); 
        Real qr = std::sqrt( 1 + 0.5*(1 + gmRec) * temp );
        
        Real sl = vxl - al*ql; 
        Real sr = vxr + ar*qr; 
        Real ss = ( prer - prel + denl*vxl*(sl-vxl) - denr*vxr*(sr-vxr) ) / 
                  ( denl*(sl-vxl) - denr*(sr-vxr) ) ;
        // Real testTerm = denl*(sl-vxl) - denr*(sr-vxr); 
        // if (std::abs(testTerm)<1e-5){
        //     std::cout<<testTerm<<std::endl;
        // } 
        
        // Step III HLLC flux
        Real el = eos.EGas(denl, prel) + 0.5*denl*(vxl*vxl+vyl*vyl+vzl*vzl); 
        Real er = eos.EGas(denr, prer) + 0.5*denr*(vxr*vxr+vyr*vyr+vzr*vzr); 

        Real cl = std::min(Real(0), sl)*(ss-vxl)/(sl-ss); 
        Real cr = std::max(Real(0), sr)*(ss-vxr)/(sr-ss); 
        Real selectl = (ss >= 0) ? 1.0 : 0.0;
        Real selectr = 1.0 - selectl;

        Real fl1 = denl * (vxl + cl);  
        Real fl2 = denl * (vxl*vxl + cl*sl) + prel; 
        Real fl3 = fl1*vyl;                             // denl*vxl*vyl     + denl*cl*vyl;
        Real fl4 = fl1*vzl;                             // denl*vxl*vzl     + denl*cl*vzl; 
        Real fl5 = (el + prel)*(vxl + cl) + denl*cl*ss*(sl-vxl);

        Real fr1 = denr * (vxr + cr);  
        Real fr2 = denr * (vxr*vxr + cr*sr) + prer; 
        Real fr3 = fr1*vyr;                             // denr*vxr*vyr     + denr*cr*vyr;
        Real fr4 = fr1*vzr;                             // denr*vxr*vzr     + denr*cr*vzr; 
        Real fr5 = (er + prer)*(vxr + cr) + denr*cr*ss*(sr-vxr);

        Real f1 = selectl * fl1 + selectr * fr1; 
        Real f2 = selectl * fl2 + selectr * fr2; 
        Real f3 = selectl * fl3 + selectr * fr3; 
        Real f4 = selectl * fl4 + selectr * fr4; 
        Real f5 = selectl * fl5 + selectr * fr5; 

        flux(DEN, k, j, i) = f1; 
        flux(IVLX, k, j, i) = f2; 
        flux(IVLY, k, j, i) = f3; 
        flux(IVLZ, k, j, i) = f4; 
        flux(ENG, k, j, i) = f5; 
    }
}

} // namespace Gaukuk