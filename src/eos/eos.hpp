#pragma once 

// C++ headers
#include <cmath>            // sqrt()

// Gaukuk dependence
#include "../gaukuk.hpp"
#include "../template_array.hpp"
#include "../grid.hpp"

namespace Gaukuk{

class EquationOfState{
public:
    EquationOfState(); 
    void ConsToPrim(const TArray<Real>& cons, TArray<Real>& prim, 
                    const Grid& grid); 
    void ConsToPrim(const TArray<Real>& cons, TArray<Real>& prim, 
                    const Grid& grid, Real& cmax); 
    void CalCmax(const TArray<Real>& cons, const Grid& grid, Real& cmax);  
    void PrimToCons(const TArray<Real>& prim, TArray<Real>& cons, 
                    const Grid& grid); 
    
    Real GetGamma() const { return gamma_; }
    Real GetGm1Rec() const { return gm1Rec_; }

    inline Real SoundSpeed(const Real den, const Real pre) const {
        return std::sqrt(gamma_*pre/den); 
    }
    inline Real EGas(const Real den, const Real pre) const {
        return pre * gm1Rec_; 
    }
    
private:
    Real gamma_, gm1Rec_; 
};


} // namespace Gaukuk