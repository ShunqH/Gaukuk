#pragma once 

// C++ headers
#include <cstddef>  // size_t

// Gaukuk dependence
#include "template_array.hpp"
#include "gaukuk.hpp" 

namespace Gaukuk{

class Sim{
public:
    int nx, ny, nz, nGhost, lenx, leny, lenz, lenArr; 

    TArray<Real> cons;
    TArray<Real> flux1, flux2, flux3; 
}; 

} // namespace Gaukuk