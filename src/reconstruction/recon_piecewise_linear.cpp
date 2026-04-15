// C++ headers 

// Gaukuk dependence
#include "reconstruction.hpp" 

namespace Gaukuk{

void Reconstruction::ReconstructXPLM(const TArray<Real>& prim, 
                         TArray<Real>& ul, TArray<Real>& ur, 
                         int k, int j, int igb, int ige){
    for (int n=0; n<NVar; n++){
    #pragma omp simd 
        for (int i=igb; i<ige-2; i++){
            int iFluxUl = i + 2; 
            int iFluxUr = i + 1; 
            Real u1 = prim(n, k, j, i);
            Real u2 = prim(n, k, j, i+1);
            Real u3 = prim(n, k, j, i+2);
            Real dul = u2-u1; 
            Real dur = u3-u2; 
            // Apply minmod limiter
            // Real sigma = Minmod(dul, dur); 
            // Apply van Leer (VL) limiter
            Real sigma = VanLeer(dul, dur); 
            // Apply monotoniced central (MC) limiter
            // Real sigma = MC(dul, dur); 

            ul(n, iFluxUl) = u2 + 0.5*sigma; 
            ur(n, iFluxUr) = u2 - 0.5*sigma; 
        }
    }
}

void Reconstruction::ReconstructYZPLM(const TArray<Real>& prim, 
                                      TArray<Real>& ur, TArray<Real>& ulNext, 
                                      int k, int j, int igb, int ige){
    for (int n=0; n<NVar; n++){
#pragma omp simd 
        for (int i=igb; i<ige; i++){
            Real u1 = prim(n, k, j-1, i);
            Real u2 = prim(n, k, j, i);
            Real u3 = prim(n, k, j+1, i);
            Real dul = u2-u1; 
            Real dur = u3-u2; 
            // Apply minmod limiter
            // Real sigma = Minmod(dul, dur); 
            // Apply van Leer (VL) limiter
            Real sigma = VanLeer(dul, dur); 
            // Apply monotoniced central (MC) limiter
            // Real sigma = MC(dul, dur); 
            ur(n, i) = u2 - 0.5*sigma; 
            ulNext(n, i) = u2 + 0.5*sigma; 
        }
    }
}


} // namespace Gaukuk