#pragma once 

// C++ headers 

// Gaukuk dependence
#include "../gaukuk.hpp" 
#include "../template_array.hpp"

namespace Gaukuk{

class Slice{
// friend class Sim; 
public:
    void ExtractXForCalFlux(const TArray<Real>& prim, 
                            TArray<Real>& ul, TArray<Real>& ur, 
                            int nVar, int k, int j, int igb, int ige); 
    void ExtractYForCalFlux(const TArray<Real>& prim, 
                            TArray<Real>& ul, TArray<Real>& ur, 
                            int nVar, int k, int j, int igb, int ige); 
    void ExtractZForCalFlux(const TArray<Real>& prim, 
                            TArray<Real>& ul, TArray<Real>& ur, 
                            int nVar, int k, int j, int igb, int ige);                         

};

} // namespace Gaukuk