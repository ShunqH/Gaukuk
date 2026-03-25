// C++ headers 

// Gaukuk dependence 
#include "../sim.hpp" 
#include "../grid/slice.hpp"

namespace Gaukuk
{
    
void Sim::CalFlux(){
    int il = grid.igb; 
    int ir = grid.ige; 
    int jl = grid.jb; 
    int jr = grid.je; 
    int kl = grid.kb;
    int kr = grid.ke; 
    
    // Flux on x direction 
#pragma omp parallel for
    for (int k=kl; k<kr; k++){
        for (int j=jl; j<jr; j++){
            slice.ExtractXForCalFlux(prim, ul_, ur_, nVar, k, j, il, ir);
            
        }
    }
}

} // namespace Gaukuk
