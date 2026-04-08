#pragma once 

#include "../gaukuk.hpp"
#include "../template_array.hpp"
#include "../grid/grid.hpp"

namespace Gaukuk{

class EquationOfState{
public:
    EquationOfState(); 
    void ConsToPrim(TArray<Real>& cons, TArray<Real>& prim, 
                    const Grid& grid); 
    void ConsToPrim(TArray<Real>& cons, TArray<Real>& prim, 
                    const Grid& grid, Real& cmax); 
    void PrimToCons(const TArray<Real>& prim, TArray<Real>& cons, 
                    const Grid& grid); 
    Real SoundSpeed(const Real den, const Real pre); 

    const Real GetGamma(){ return gamma_; }
    const Real GetGm1Rec(){ return gm1Rec_; }
    Real EGas(const Real den, const Real pre); 
    
private:
    Real gamma_, densityMin_, pressureMin_, gm1Rec_; 
};


} // namespace Gaukuk