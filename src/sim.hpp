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
    Flux flux; 
    Real t, dt, cmax, CFL;  
    
    TArray<Real> cons, prim;
    TArray<Real> flx1, flx2, flx3; 

    Sim(); 

    EquationOfState eos; 
    Boundary boundary; 
    
    void UpdateCons(); 
private:
    
}; 

} // namespace Gaukuk