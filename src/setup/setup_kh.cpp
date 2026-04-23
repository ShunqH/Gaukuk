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
                Real xNow = domain.xc(i); 
                Real yNow = domain.yc(j); 

                Real rhoNow = (std::abs(yNow)<yLayer) ? rhoIn : rhoOut ; 
                Real vx = (std::abs(yNow)<yLayer) ? vxIn : vxOut ; 

                Real vxRandom = amp * std::sin( 4 * PI * (xNow-xmin) / Lx ); 
                Real vyRandom = amp * std::cos( 4 * PI * (xNow-xmin) / Lx );  

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
