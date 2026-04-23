// C++ headers

// Gaukuk dependence
#include "../gaukuk.hpp" 
#include "../sim.hpp" 

namespace Gaukuk
{
    
void Sim::Setup(){
    // load from setupfile 
    Real x0 = Config::getInstance().get("x0"); 
    Real rhoLeft = Config::getInstance().get("rhoLeft"); 
    Real rhoRight = Config::getInstance().get("rhoRight"); 
    Real preLeft = Config::getInstance().get("pressureLeft"); 
    Real preRight = Config::getInstance().get("pressureRight"); 

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
                Real rhoNow = (xNow<x0) ? rhoLeft : rhoRight ; 
                Real preNow = (xNow<x0) ? preLeft : preRight ; 

                // usually you have to setup conservative quantivities (cons) 
                // but you can setup primitive quantivities (cons) 
                // then use eos.PrimToCons to conver 
                prim(DEN, k, j, i) = rhoNow; 
                prim(VLX, k, j, i) = 0; 
                prim(VLY, k, j, i) = 0; 
                prim(VLY, k, j, i) = 0; 
                prim(PRE, k, j, i) = preNow; 
            }
        }
    }
    // if primitive quantivities is set, make sure you call eos.PrimToCons 
    eos.PrimToCons(prim, cons, grid); 
}

} // namespace Gaukuk
