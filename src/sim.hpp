#pragma once 

// C++ headers
#include <cstddef>  // size_t

// Gaukuk dependence
#include "template_array.hpp"
#include "gaukuk.hpp" 

namespace Gaukuk{

class Domain{
public:
friend class Gaukuk; 
    void SetDomain(const int nx, const int ny, const int nz); 

private:
    Real xmin, xmax, ymin, ymax, zmin, zmax; 
    Real dx, dy, dz; 
    TArray<Real> xGrid, yGrid, zGrid; 
};


class Gaukuk{
public:
    Domain domain; 
    Real t, dt, cmax, CFL;  
    int nx, ny, nz, nGhost, lenx, leny, lenz, lenArr; 

    TArray<Real> cons, prim;
    TArray<Real> flx1, flx2, flx3; 

    Gaukuk(); 

    void UpdateBD(); 
    void UpdateFlux(); 
    void RiemannSolver(); 
    void UpdateCons(); 
}; 

} // namespace Gaukuk