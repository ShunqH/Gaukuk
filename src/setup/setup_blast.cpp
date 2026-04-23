// C++ headers
#include <cmath>        //std::sin
#include <iostream>     // std::cout; std::endl; std::cerr

// Gaukuk dependence
#include "../gaukuk.hpp" 
#include "../sim.hpp" 

namespace Gaukuk
{
    
void Sim::Setup(){
    // Load parameters from input file
    Real rho0    = Config::getInstance().get("rho0");      // uniform density
    Real engTot = Config::getInstance().get("eng_tot");   // total energy injected at the center region 
    Real pLow    = Config::getInstance().get("pressure_low");   // pressure floor 
    Real radius0 = Config::getInstance().get("r0");        // initial blast radius
    Real gm1Rec = eos.GetGm1Rec(); 

    // Domain geometry
    Real xmin = domain.xmin;
    Real xmax = domain.xmax;
    Real ymin = domain.ymin;
    Real ymax = domain.ymax;
    Real zmin = domain.zmin;
    Real zmax = domain.zmax;

    // Center of the domain
    Real xc = 0.5 * (xmin + xmax);
    Real yc = 0.5 * (ymin + ymax);
    Real zc = 0.5 * (zmin + zmax);

    // activated zone
    int il = grid.ib; 
    int ir = grid.ie; 
    int jl = grid.jb; 
    int jr = grid.je; 
    int kl = grid.kb;
    int kr = grid.ke; 

    int count = 0; 
    #pragma omp parallel for collapse(2) schedule(static) reduction(+:count)
    for (int k=kl; k<kr; k++){
        for (int j=jl; j<jr; j++){
            for (int i=il; i<ir; i++){
                Real xNow = domain.xc(i);
                Real yNow = domain.yc(j);
                Real zNow = domain.zc(k);

                // Distance from center
                Real dx = xNow - xc;
                Real dy = yNow - yc;
                Real dz = zNow - zc;
                Real r = std::sqrt(dx*dx + dy*dy + dz*dz);

                if (r < radius0) {
                    count++; 
                }
            }
        }
    }

    Real V = (4.0/3.0) * PI * radius0 * radius0 * radius0; 
    // Real V = count * domain.dx*domain.dy*domain.dz; 
    Real eHigh = engTot/V; 
    Real eLow = pLow * gm1Rec; 
    // std::cout<< eHigh << "  " << eLow << std::endl; 

    #pragma omp parallel for collapse(2) schedule(static)    
    for (int k=kl; k<kr; k++){
        for (int j=jl; j<jr; j++){
    #pragma omp simd
            for (int i=il; i<ir; i++){
                Real xNow = domain.xc(i);
                Real yNow = domain.yc(j);
                Real zNow = domain.zc(k);

                // Distance from center
                Real dx = xNow - xc;
                Real dy = yNow - yc;
                Real dz = zNow - zc;
                Real r = std::sqrt(dx*dx + dy*dy + dz*dz);

                // energy is high inside sphere, low outside
                Real engNow = (r < radius0) ? eHigh : eLow;

                // Uniform density, zero velocity
                cons(DEN, k, j, i) = rho0;
                cons(MTX, k, j, i) = 0.0;
                cons(MTY, k, j, i) = 0.0;
                cons(MTZ, k, j, i) = 0.0;
                cons(ENG, k, j, i) = engNow;
            }
        }
    }
}

} // namespace Gaukuk
