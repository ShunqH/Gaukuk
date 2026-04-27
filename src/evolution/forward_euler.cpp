// C++ Header
#include <algorithm>        //std::min 

// Gaukuk dependence
#include "../gaukuk.hpp" 
#include "../sim.hpp"

namespace Gaukuk
{

void Sim::ForwardEuler_(){
    boundary.UpdateBD(cons, grid); 
    eos.ConsToPrim(cons, prim, grid); 
    flux.CalFlux(grid, prim, eos, flx1, flx2, flx3, rcOrder); 
    UpdateCons(cons, consTemp, 1.0, 1.0); 
    cons.Swap(consTemp); 
}

} // namespace Gaukuk
