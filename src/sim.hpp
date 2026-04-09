#pragma once 

// C++ headers
#include <cstddef>  // size_t

// Gaukuk dependence
#include "gaukuk.hpp" 
#include "template_array.hpp"
#include "grid/grid.hpp"
#include "grid/reconstruction.hpp" 
#include "eos/eos.hpp"
#include "flux/flux.hpp"
#include "boundary/boundary.hpp" 

namespace Gaukuk{

enum class DataType {
    Prim,
    Cons
};

class Domain{
public:
friend class Sim; 
    Domain(const Grid& grid); 

private:
    Real xmin, xmax, ymin, ymax, zmin, zmax; 
    Real dx, dy, dz, drmin; 
    TArray<Real> xc, yc, zc;       // cell center 
};

class Sim{
public:

    Sim(); 
    const Grid grid; 
    const Domain domain; 
    Flux flux; 

    TArray<Real> cons, prim;
    TArray<Real> flx1, flx2, flx3; 

    EquationOfState eos; 
    Boundary boundary; 
    
    void Setup(); 
    void Advance(Real dtoutput);
    
    void WriteData(const int outputID, DataType dType); 

    Real GetTime(){ return t; };
    Real Getdt(){ return dt; };
private:
using IntegratorFunc = void (Sim::*)();
    Real t, dt, dtUntilOutput, cmax, CFL; 
    IntegratorFunc integrator_; 
    void ForwardEuler_(); 
    // void RK2(); 
    TArray<Real> consTemp1_, consTemp2_; 
}; 

} // namespace Gaukuk