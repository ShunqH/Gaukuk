#pragma once 

// C++ headers
#include <cstddef>  // size_t

// Gaukuk dependence
#include "template_array.hpp"
#include "eos/eos.hpp" 
#include "gaukuk.hpp" 

namespace Gaukuk{

class Grid{
public:
friend class Sim; 
    Grid(); 
    int nx, ny, nz, nGhost, lenx, leny, lenz, lenArr; 
};

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

    void UpdateBD(); 
    void UpdateFlux(); 
    void RiemannSolver(); 
    void UpdateCons(); 
// private:
    
}; 

} // namespace Gaukuk