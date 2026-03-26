// C++ headers 

// Gaukuk dependence
#include "slice.hpp" 

namespace Gaukuk{

void Slice::ExtractXForCalFlux(const TArray<Real>& prim, 
                            TArray<Real>& ul, TArray<Real>& ur, 
                            int k, int j, int igb, int ige){
    for (int n=0; n<NVar; n++){
#pragma omp simd 
        for (int i=igb; i<ige; i++){
            Real u = prim(n, k, j, i);
            ul(n, i+1) = u; 
            ur(n, i) = u; 
        }
    }
}

void Slice::ExtractYForCalFlux(const TArray<Real>& prim, TArray<Real>& ur, 
                            int k, int j, int igb, int ige){
    for (int n=0; n<NVar; n++){
#pragma omp simd 
        for (int i=igb; i<ige; i++){
            ur(n, i) = prim(n, k, j, i);
        }
    }
}

void Slice::ExtractZForCalFlux(const TArray<Real>& prim, TArray<Real>& ur, 
                            int k, int j, int igb, int ige){
    for (int n=0; n<NVar; n++){
#pragma omp simd 
        for (int i=igb; i<ige; i++){
            ur(n, i) = prim(n, k, j, i);
        }
    }
}

} // namespace Gaukuk