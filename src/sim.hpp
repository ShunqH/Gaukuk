#pragma once 

// C++ headers
#include <cstddef>  // size_t

// Gaukuk dependence
#include "template_array.hpp"
#include "grid/grid.hpp"
#include "grid/slice.hpp" 
#include "eos/eos.hpp"
#include "gaukuk.hpp" 

namespace Gaukuk{

class Domain{
public:
friend class Sim; 
    Domain(const Grid& grid); 

private:
    Real xmin, xmax, ymin, ymax, zmin, zmax; 
    Real dx, dy, dz; 
    TArray<Real> xGrid, yGrid, zGrid; 
};

class Sim{
public:
    const Grid grid; 
    const Domain domain; 
    const int nVar; 
    Real t, dt, cmax, CFL;  
    
    TArray<Real> cons, prim;
    TArray<Real> flx1, flx2, flx3; 

    Sim(); 

    EquationOfState eos; 
    Slice slice; 

    void CalFlux(); 
    void RiemannSolver(); 

    void UpdateBD(); 
    void UpdateCons(); 
private:
    TArray<Real> ul_, ur_;
    
}; 

} // namespace Gaukuk