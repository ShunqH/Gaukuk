// C++ headers 

// Gaukuk dependence 
#include "../template_array.hpp"// TArray
#include "../flux/flux.hpp"     // class Flux   
#include "../grid.hpp"     // class grid
#include "../reconstruction/reconstruction.hpp"    // class reconstruction; void reconstructXForFlux

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
    #pragma omp parallel
    {
        TArray<Real> ul_(NVar, lenUlUr);
        TArray<Real> ur_(NVar, lenUlUr);
    #pragma omp for collapse(2) schedule(static)
        for (int k = kl; k < kr; k++) {
            for (int j = jl; j < jr; j++) {
                recon.ReconstructXForFlux(prim, ul_, ur_, k, j, il, ir);
                RiemannSolver(ul_, ur_, VLX, eos, flx1, k, j, il, ir);
            }
        }
    }

    // Flux on y direction 
    #pragma omp parallel
    {
        TArray<Real> ul_(NVar, lenUlUr);
        TArray<Real> ur_(NVar, lenUlUr);
    #pragma omp for schedule(static)
        for (int k=kl; k<kr; k++){
            recon.ReconstructYZForFlux(prim, ul_, k, jl-1, il, ir);
            for (int j=jl; j<jr+1; j++){
                recon.ReconstructYZForFlux(prim, ur_, k, j, il, ir);
                RiemannSolver(ul_, ur_, VLY, eos, flx2, k, j, il, ir); 
                ul_.Swap(ur_); 
            }
        }
    }
    
    // Flux on z direction 
    #pragma omp parallel
    {
        TArray<Real> ul_(NVar, lenUlUr);
        TArray<Real> ur_(NVar, lenUlUr);
    #pragma omp for schedule(static)
        for (int j=jl; j<jr; j++){
            recon.ReconstructYZForFlux(prim, ul_, kl-1, j, il, ir);
            for (int k=kl; k<kr+1; k++){
                recon.ReconstructYZForFlux(prim, ur_, k, j, il, ir);
                RiemannSolver(ul_, ur_, VLZ, eos, flx3, k, j, il, ir); 
                ul_.Swap(ur_); 
            }
        }
    }

}


} // namespace Gaukuk
