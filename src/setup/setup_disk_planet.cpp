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
    Real rs0 = Config::getInstance().get("rs0", 0.2);
    Real GMp = Config::getInstance().get("GMp", 1);
    Real rp = Config::getInstance().get("rp", 1);
    Real rs1 = Config::getInstance().get("rs1", 0.1);

    Real r0 = Config::getInstance().get("r0", 1.0); 
    Real rho0 = Config::getInstance().get("rho0", 1.0); 
    Real p = Config::getInstance().get("p", 1.0); 
    Real cs0 = Config::getInstance().get("cs0", 0.1); 
    Real q = Config::getInstance().get("q", 1.0); 

    Real xShell = Config::getInstance().get("x0", 0.0);  
    Real yShell = Config::getInstance().get("y0", 0.0);  
    Real zShell = Config::getInstance().get("z0", 0.0);  
    Real rShell = Config::getInstance().get("r_inner"); 

    Real r0Inv = 1.0/r0; 
    Real gammaInv = 1.0/gamma; 

    // int il = grid.ib;
    // int ir = grid.ie;
    // int jl = grid.jb;
    // int jr = grid.je;
    // int kl = grid.kb;
    // int kr = grid.ke;

    int il = grid.igb;
    int ir = grid.ige;
    int jl = grid.jgb;
    int jr = grid.jge;
    int kl = grid.kgb;
    int kr = grid.kge;

    // enroll the central point source 
    srcTerm.EnrollPointGravity(GM, 0, 0, 0, 0, 0, 0, rs0); 
    srcTerm.EnrollPointGravity(GMp, rp, 0, 0, 0, 0, 0, rs1); 

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

                // Real rhoNow = (r<0.4) ? 1e-6 : rho0 * std::pow(r*r0Inv, -p);
                Real rhoNow = rho0 * std::pow(r*r0Inv, -p);

                Real cs = cs0 * std::pow(r*r0Inv, -0.5*q); 
                Real pressure = gammaInv*cs*cs*rhoNow; 

                Real r2s = r2 + rs0*rs0;
                Real v_phi2 = GM * r2 / (r2s * std::sqrt(r2s))  
                              - gammaInv * (p + q) * cs * cs ;
                if (v_phi2 < 0) { v_phi2 = 0; }
                Real v_phi = std::sqrt( v_phi2 ); 
                // Real v_phi = std::sqrt( GM * r2 / (r2s * std::sqrt(r2s)) 
                                        // - (p + q) * cs * cs * gammaInv );
                // Real v_phi = std::sqrt( GM  / std::sqrt(r2s) 
                //                         - gammaInv * (p + q) * cs * cs ;
                // vx = -v_phi * sinθ = -v_phi * (y/r)
                // vy =  v_phi * cosθ =  v_phi * (x/r)
                Real vx = -v_phi * (yNow * rInv);
                Real vy =  v_phi * (xNow * rInv);

                // if (r < rs0) {
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
    boundary.EnrollFixedBoundaryData(cons, grid, true); 
    // boundary.SetupImmersedShell(grid, domain, xShell, yShell, zShell, rShell); 
}

} // namespace Gaukuk
