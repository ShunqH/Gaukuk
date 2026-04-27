// C++ Header
#include <algorithm>        //std::min 

// Gaukuk dependence
#include "../gaukuk.hpp" 
#include "../sim.hpp"

namespace Gaukuk
{

void Sim::RK3_(){
    boundary.UpdateBD(cons, grid); 
    eos.ConsToPrim(cons, prim, grid); 
    flux.CalFlux(grid, prim, eos, flx1, flx2, flx3, rcOrder); 
    UpdateCons(cons, consTemp, 1.0, 1.0); 

    boundary.UpdateBD(consTemp, grid); 
    eos.ConsToPrim(consTemp, prim, grid); 
    flux.CalFlux(grid, prim, eos, flx1, flx2, flx3, rcOrder); 
    UpdateCons(cons, consTemp, 0.75, 0.25, 0.25); 

    boundary.UpdateBD(consTemp, grid); 
    eos.ConsToPrim(consTemp, prim, grid); 
    flux.CalFlux(grid, prim, eos, flx1, flx2, flx3, rcOrder); 
    Real frac13 = 1.0/3.0; 
    UpdateCons(cons, consTemp, frac13, 2*frac13, 2*frac13); 
    cons.Swap(consTemp);  

}

} // namespace Gaukuk
