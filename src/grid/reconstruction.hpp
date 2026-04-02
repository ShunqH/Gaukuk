#pragma once 

// C++ headers 

// Gaukuk dependence
#include "../gaukuk.hpp" 
#include "../template_array.hpp"

namespace Gaukuk{

class Reconstruction{
// friend class Sim; 
public:
    void ReconstructXForFlux(const TArray<Real>& prim, 
                             TArray<Real>& ul, TArray<Real>& ur, 
                             int k, int j, int igb, int ige); 
    void ReconstructYZForFlux(const TArray<Real>& prim, TArray<Real>& ur, 
                              int k, int j, int igb, int ige); 
};

} // namespace Gaukuk