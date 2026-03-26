// C++ headers 

// Gaukuk dependence 
#include "../flux/flux.hpp"     // class Flux   
#include "../grid/grid.hpp"     // class grid
#include "../grid/slice.hpp"    // class slice; void ExtractXForCalFlux

namespace Gaukuk
{
    
void Flux::CalFlux(const Grid& grid, const TArray<Real>& prim, EquationOfState& eos,
                   TArray<Real>& flx1, TArray<Real>& flx2, TArray<Real>& flx3){
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
            slice.ExtractXForCalFlux(prim, ul_, ur_, k, j, il, ir);
            RiemannSolver(ul_, ur_, VLX, eos, flx1, k, j, il, ir); 
        }
    }
}

} // namespace Gaukuk
