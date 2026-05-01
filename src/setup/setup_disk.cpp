// C++ headers
#include <cmath>        //std::sin

// Gaukuk dependence
#include "../gaukuk.hpp" 
#include "../sim.hpp" 

namespace Gaukuk
{
    
void Sim::Setup() {
    Real gamma = Config::getInstance().get("gamma");
    Real GM = Config::getInstance().get("GM", 1);
    Real rs = Config::getInstance().get("rs", 0.2);

    Real r0 = Config::getInstance().get("r0", 1.0); 
    Real rho0 = Config::getInstance().get("rho0", 1.0); 
    Real p = Config::getInstance().get("p", 1.0); 
    Real cs0 = Config::getInstance().get("cs0", 0.1); 
    Real q = Config::getInstance().get("q", 1.0); 

    Real gapWidth = Config::getInstance().get("gapWidth", 0.05); 
    Real gapDepth = Config::getInstance().get("gapDepth", 0.5); 
    Real amp = Config::getInstance().get("amp", 0.05); 
    Real km = Config::getInstance().get("km", 5); 

    Real r0Inv = 1.0/r0; 
    Real deltaInv = 1.0/gapWidth; 
    Real gammaInv = 1.0/gamma; 

    int il = grid.ib;
    int ir = grid.ie;
    int jl = grid.jb;
    int jr = grid.je;
    int kl = grid.kb;
    int kr = grid.ke;

    // enroll the central point source 
    srcTerm.EnrollPointGravity(GM, 0, 0, 0, 0, 0, 0, rs); 

#pragma omp parallel for collapse(2) schedule(static)
    for (int k = kl; k < kr; ++k) {
        for (int j = jl; j < jr; ++j) {
            for (int i = il; i < ir; ++i) {
                Real xNow = domain.xc(i);
                Real yNow = domain.yc(j);
                // Real zNow = domain.zc(k);

                Real r2 = xNow*xNow + yNow*yNow ; 
                Real r = std::sqrt(r2); 
                Real rInv = 1.0/(r + 1e-12); 

                Real rhob = rho0 * std::pow(r*r0Inv, -p);
                Real gaussian = 1 - gapDepth * std::exp(-0.5*(r-r0)*(r-r0)*deltaInv*deltaInv); 
                Real rhoNow = (r<rs) ? 1e-6 : rhob * gaussian; 
                // Real rhoNow = rhob * gaussian;

                Real cs = cs0 * std::pow(r*r0Inv, -0.5*q); 
                Real pressure = gammaInv*cs*cs*rhoNow; 

                Real theta = std::atan2(yNow, xNow); 
                Real vPtb = amp*(std::cos( km*theta )) ;

                Real r2s = r2 + rs*rs;
                Real v_phi = std::sqrt( GM * r2 / (r2s * std::sqrt(r2s))  
                                        + (p + q) * cs * cs ) 
                             + vPtb;
                // Real v_phi = std::sqrt( GM  / std::sqrt(r2s) 
                //                         - (p + q) * cs * cs ) ;
                // vx = -v_phi * sinθ = -v_phi * (y/r)
                // vy =  v_phi * cosθ =  v_phi * (x/r)
                Real vx = -v_phi * (yNow * rInv);
                Real vy =  v_phi * (xNow * rInv);

                // if (r < rs) {
                //     // Real factor = r / rs;
                //     // vx *= factor;
                //     // vy *= factor;
                //     vx = 0;
                //     vy = 0;
                // }

                prim(DEN, k, j, i) = rhoNow;
                prim(VLX, k, j, i) = vx;
                prim(VLY, k, j, i) = vy;
                prim(VLZ, k, j, i) = 0;
                prim(PRE, k, j, i) = pressure;
            }
        }
    }

    // if primitive quantivities is set, make sure you call eos.PrimToCons 
    eos.PrimToCons(prim, cons, grid);
}

} // namespace Gaukuk
