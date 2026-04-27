#pragma once 

// C++ headers
#include <cstddef>  // size_t

// Gaukuk dependence
#include "gaukuk.hpp" 
#include "template_array.hpp"
#include "grid.hpp"
#include "reconstruction/reconstruction.hpp"
#include "eos/eos.hpp"
#include "flux/flux.hpp"
#include "boundary/boundary.hpp" 
#include "source_term/source.hpp"

namespace Gaukuk{

enum class DataType {
    Prim,
    Cons
};

class Sim{
public:
    Sim(); 
    bool isContinue; 
    const Grid grid; 
    const Domain domain; 
    Flux flux; 

    TArray<Real> cons, prim;
    TArray<Real> flx1, flx2, flx3; 

    EquationOfState eos; 
    Boundary boundary; 
    SourceTerm srcTerm; 
    
    void Setup(); 
    void Advance(Real dtoutput);
    
    void WriteData(const int outputID, DataType dType); 

    Real GetTime(){ return t; }
    Real Getdt(){ return dt; }
    int GetStep() { return step; }
private:
using VoidFunc = void (Sim::*)();
    int step, stepMax; 
    Real t, dt, dtUntilOutput, cmax, CFL; 
    int rcOrder, integratorType; 
    VoidFunc hydroIntegrator_; 
    void UpdateCons(const TArray<Real>& cons_, TArray<Real>& consTemp_, const Real coef1, const Real coef2); 
    void UpdateCons(const TArray<Real>& cons_, TArray<Real>& consTemp_, const Real coef1, const Real coef2, const Real coef3); 
    void ForwardEuler_(); 
    void RK2_(); 
    void RK3_(); 
    TArray<Real> consTemp, consTemp2; 
}; 

} // namespace Gaukuk