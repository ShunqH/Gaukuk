// C++ headers 
#include <algorithm>        //std::min() 

// Gaukuk dependence 
#include "../template_array.hpp"// TArray
#include "../flux/flux.hpp"     // class Flux   
#include "../grid.hpp"     // class grid
#include "../reconstruction/reconstruction.hpp"    // class reconstruction; void reconstructXFirstOrder

namespace Gaukuk
{
    
void Flux::CalFlux(const Grid& grid, const TArray<Real>& prim, EquationOfState& eos,
                   TArray<Real>& flx1, TArray<Real>& flx2, TArray<Real>& flx3, const int rcOrder){
    int il = grid.igb; 
    int ir = grid.ige; 
    int jl = grid.jb; 
    int jr = grid.je; 
    int kl = grid.kb;
    int kr = grid.ke; 
    
if (rcOrder == 1){
    // Flux on x direction 
    #pragma omp parallel
    {
        TArray<Real> ul_(NVar, lenUlUr);
        TArray<Real> ur_(NVar, lenUlUr);
    #pragma omp for collapse(2) schedule(static)
        for (int k = kl; k < kr; k++) {
            for (int j = jl; j < jr; j++) {
                recon.ReconstructXFirstOrder(prim, ul_, ur_, k, j, il, ir);
                RiemannSolver(ul_, ur_, VLX, eos, flx1, k, j, il, ir);
            }
        }
    }

  if (grid.nz==1){
    int blockJ = 32; 
    int k = kl; 
    // Flux on y direction 
    #pragma omp parallel
    {
        TArray<Real> ul_(NVar, lenUlUr);
        TArray<Real> ur_(NVar, lenUlUr);
    #pragma omp for schedule(static)
        for (int jj=jl; jj<jr+1; jj+=blockJ){
            recon.ReconstructYZFirstOrder(prim, ul_, k, jj-1, il, ir);
            for (int j=jj; j<std::min(jj+blockJ, jr+1); j++){
                recon.ReconstructYZFirstOrder(prim, ur_, k, j, il, ir);
                RiemannSolver(ul_, ur_, VLY, eos, flx2, k, j, il, ir); 
                ul_.Swap(ur_); 
            }
        }
    }
  }else{
    // Flux on y direction 
    #pragma omp parallel
    {
        TArray<Real> ul_(NVar, lenUlUr);
        TArray<Real> ur_(NVar, lenUlUr);
    #pragma omp for schedule(static)
        for (int k=kl; k<kr; k++){
            recon.ReconstructYZFirstOrder(prim, ul_, k, jl-1, il, ir);
            for (int j=jl; j<jr+1; j++){
                recon.ReconstructYZFirstOrder(prim, ur_, k, j, il, ir);
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
            recon.ReconstructYZFirstOrder(prim, ul_, kl-1, j, il, ir);
            for (int k=kl; k<kr+1; k++){
                recon.ReconstructYZFirstOrder(prim, ur_, k, j, il, ir);
                RiemannSolver(ul_, ur_, VLZ, eos, flx3, k, j, il, ir); 
                ul_.Swap(ur_); 
            }
        }
    }
  }
}else {
    // Flux on x direction 
    #pragma omp parallel
    {
        TArray<Real> ul_(NVar, lenUlUr);
        TArray<Real> ur_(NVar, lenUlUr);
    #pragma omp for collapse(2) schedule(static)
        for (int k = kl; k < kr; k++) {
            for (int j = jl; j < jr; j++) {
                recon.ReconstructXPLM(prim, ul_, ur_, k, j, il, ir);
                RiemannSolver(ul_, ur_, VLX, eos, flx1, k, j, il, ir);
            }
        }
    }

  if (grid.nz==1){
    int blockJ = 32; 
    int k = kl; 
    // Flux on y direction 
    #pragma omp parallel
    {
        TArray<Real> ul_(NVar, lenUlUr);
        TArray<Real> ur_(NVar, lenUlUr);
        TArray<Real> ulNext_(NVar, lenUlUr);
    #pragma omp for schedule(static)
        for (int jj=jl; jj<jr+1; jj+=blockJ){
            recon.ReconstructYZPLM(prim, ur_, ul_, k, jj-1, il, ir);
            for (int j=jj; j<std::min(jj+blockJ, jr+1); j++){
                recon.ReconstructYZPLM(prim, ur_, ulNext_, k, j, il, ir);
                RiemannSolver(ul_, ur_, VLY, eos, flx2, k, j, il, ir); 
                ul_.Swap(ulNext_); 
            }
        }
    }
  }else{
    // Flux on y direction 
    #pragma omp parallel
    {
        TArray<Real> ul_(NVar, lenUlUr);
        TArray<Real> ur_(NVar, lenUlUr);
        TArray<Real> ulNext_(NVar, lenUlUr);
    #pragma omp for schedule(static)
        for (int k=kl; k<kr; k++){
            recon.ReconstructYZPLM(prim, ur_, ul_, k, jl-1, il, ir);
            for (int j=jl; j<jr+1; j++){
                recon.ReconstructYZPLM(prim, ur_, ulNext_, k, j, il, ir);
                RiemannSolver(ul_, ur_, VLY, eos, flx2, k, j, il, ir); 
                ul_.Swap(ulNext_); 
            }
        }
    }

    // Flux on z direction 
    #pragma omp parallel
    {
        TArray<Real> ul_(NVar, lenUlUr);
        TArray<Real> ur_(NVar, lenUlUr);
        TArray<Real> ulNext_(NVar, lenUlUr);
    #pragma omp for schedule(static)
        for (int j=jl; j<jr; j++){
            recon.ReconstructYZPLM(prim, ur_, ul_, kl-1, j, il, ir);
            for (int k=kl; k<kr+1; k++){
                recon.ReconstructYZPLM(prim, ur_, ulNext_, k, j, il, ir);
                RiemannSolver(ul_, ur_, VLZ, eos, flx3, k, j, il, ir); 
                ul_.Swap(ulNext_); 
            }
        }
    }
  }
}

}


} // namespace Gaukuk
