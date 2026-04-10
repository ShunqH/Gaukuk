// C++ headers
#include <cmath>        //std::sin

// Gaukuk dependence
#include "../gaukuk.hpp" 
#include "../sim.hpp" 

namespace Gaukuk
{
    
void Sim::Setup(){
    // load from setupfile 
    Real rhoIn = Config::getInstance().get("rhoIn"); 
    Real rhoOut = Config::getInstance().get("rhoOut"); 
    Real vxIn = Config::getInstance().get("vxIn"); 
    Real vxOut = Config::getInstance().get("vxOut"); 
    Real pressure = Config::getInstance().get("pressure"); 
    Real amp = Config::getInstance().get("amp"); 

    Real yLayer = domain.ymax/2; 
    Real xmin = domain.xmin; 
    Real Lx = domain.xmax - domain.xmin; 
    Real ymin = domain.ymin; 
    Real Ly = domain.ymax - domain.ymin; 
    const Real PI = 3.1415926535; 

    // activated zone
    int il = grid.ib; 
    int ir = grid.ie; 
    int jl = grid.jb; 
    int jr = grid.je; 
    int kl = grid.kb;
    int kr = grid.ke; 

#pragma omp parallel for collapse(2) schedule(static)    
    for (int k=kl; k<kr; k++){
        for (int j=jl; j<jr; j++){
#pragma omp simd
            for (int i=il; i<ir; i++){
                int iForXc = i - grid.ib; 
                int jForYc = j - grid.jb; 
                Real xNow = domain.xc(iForXc); 
                Real yNow = domain.yc(jForYc); 

                Real rhoNow = (std::abs(yNow)<yLayer) ? rhoIn : rhoOut ; 
                Real vx = (std::abs(yNow)<yLayer) ? vxIn : vxOut ; 

                Real vxRandom = amp * std::sin( 4 * PI * (xNow-xmin) / Lx )
                                    * std::sin( 2 * PI * (yNow-ymin) / Ly ) ; 
                Real vyRandom = amp * std::sin( 4 * PI * (xNow-xmin) / Lx + 0.87*PI)
                                    * std::sin( 2 * PI * (yNow-ymin) / Ly + 1.23*PI) ; 

                // usually you have to setup conservative quantivities (cons) 
                // but you can setup primitive quantivities (cons) 
                // then use eos.PrimToCons to conver 
                prim(DEN, k, j, i) = rhoNow; 
                prim(VLX, k, j, i) = vx + vxRandom; 
                prim(VLY, k, j, i) = vyRandom; 
                prim(VLZ, k, j, i) = 0; 
                prim(PRE, k, j, i) = pressure; 
            }
        }
    }
    // if primitive quantivities is set, make sure you call eos.PrimToCons 
    eos.PrimToCons(prim, cons, grid); 
}

} // namespace Gaukuk
