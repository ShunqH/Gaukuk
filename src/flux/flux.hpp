#pragma once 

// C++ headers

// Gaukuk dependence
#include "../gaukuk.hpp"            //Real 
#include "../template_array.hpp"    // TArray
#include "../grid/grid.hpp"         // class grid
#include "../grid/reconstruction.hpp"        // class slice; void ExtractXForCalFlux 
#include "../eos/eos.hpp"   // EquationOfState

namespace Gaukuk{

class Flux{
friend class Sim; 
public:
    Flux(const int lenx): lenUlUr(lenx+1){
        // ul_.NewArray(NVar, lenx+1);
        // ur_.NewArray(NVar, lenx+1);  
    }
    void CalFlux(const Grid& grid, const TArray<Real>& prim, EquationOfState& eos,
                 TArray<Real>& flx1, TArray<Real>& flx2, TArray<Real>& flx3); 
    void RiemannSolver(const TArray<Real>& ul, const TArray<Real>& ur, const int direction, 
                       EquationOfState& eos, TArray<Real>& flux, 
                       int k, int j, int igb, int ige); 
private:
    Reconstruction recon; 
    const int lenUlUr; 
    // TArray<Real> ul_, ur_;
};
    
} // namespace Gaukuk