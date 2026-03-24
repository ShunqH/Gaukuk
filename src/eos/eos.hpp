#pragma once 

#include "../gaukuk.hpp"
#include "../template_array.hpp"

namespace Gaukuk{

class EquationOfState{
public:
    EquationOfState(); 
    void ConsToPrim(TArray<Real>& cons, TArray<Real>& prim, 
                    int ib, int ie, int jb, int je, int kb, int ke); 
    void PrimToCons(const TArray<Real>& prim, TArray<Real>& cons, 
                    int ib, int ie, int jb, int je, int kb, int ke); 
    Real SoundSpeed(const Real den, const Real pre); 
    
private:
    Real gamma_, densityMin_, pressureMin_; 
};


} // namespace Gaukuk