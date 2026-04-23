// C++ Header
#include <algorithm>        //std::min 
// #include <iostream>     // std::cout; std::endl; std::cerr

// Gaukuk dependence
#include "../gaukuk.hpp" 
#include "../sim.hpp"

namespace Gaukuk
{

void Sim::ForwardEuler_(){
    int il = grid.ib; 
    int ir = grid.ie; 
    int jl = grid.jb; 
    int jr = grid.je; 
    int kl = grid.kb;
    int kr = grid.ke; 

    Real dtdx = dt*domain.dxRec; 
    Real dtdy = dt*domain.dyRec; 
    Real dtdz = dt*domain.dzRec; 

    boundary.UpdateBD(cons, grid); 
    eos.ConsToPrim(cons, prim, grid); 
    flux.CalFlux(grid, prim, eos, flx1, flx2, flx3, rcOrder); 
if (grid.nz == 1){ 
    int k = kl; 
    #pragma omp parallel for collapse(2) schedule(static)    
    for (int ivar=DEN; ivar<=ENG; ivar++){
        for (int j=jl; j<jr; j++){
    #pragma omp simd
            for (int i=il; i<ir; i++){
                Real& u = cons(ivar, k, j, i); 

                const Real fl = flx1(ivar, k, j, i); 
                const Real fr = flx1(ivar, k, j, i+1); 
                const Real gl = flx2(ivar, k, j, i); 
                const Real gr = flx2(ivar, k, j+1, i); 

                u = u - dtdx * ( fr - fl ) 
                        - dtdy * ( gr - gl ) ; 
            }
        }
    }  
}else{
    #pragma omp parallel for collapse(3) schedule(static)    
    for (int ivar=DEN; ivar<=ENG; ivar++){
        for (int k=kl; k<kr; k++){
            for (int j=jl; j<jr; j++){
    #pragma omp simd
                for (int i=il; i<ir; i++){
                    Real& u = cons(ivar, k, j, i); 

                    const Real fl = flx1(ivar, k, j, i); 
                    const Real fr = flx1(ivar, k, j, i+1); 
                    const Real gl = flx2(ivar, k, j, i); 
                    const Real gr = flx2(ivar, k, j+1, i); 
                    const Real hl = flx3(ivar, k, j, i); 
                    const Real hr = flx3(ivar, k+1, j, i); 

                    u = u - dtdx * ( fr - fl ) 
                          - dtdy * ( gr - gl ) 
                          - dtdz * ( hr - hl ) ; 
                }
            }
        }
    }
}
}

} // namespace Gaukuk
