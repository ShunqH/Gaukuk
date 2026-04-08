#pragma once 

// C++ headers

// Gaukuk dependence
#include "../gaukuk.hpp" 
#include "../template_array.hpp" 
#include "../grid/grid.hpp"

namespace Gaukuk
{
    
class Boundary{
using BoundaryFunc = void (*)(TArray<Real>&, const Grid&);
public:
friend class Sim; 
    Boundary(); 
    BoundaryFunc Bdxl; 
    BoundaryFunc Bdxr; 
    BoundaryFunc Bdyl; 
    BoundaryFunc Bdyr; 
    BoundaryFunc Bdzl; 
    BoundaryFunc Bdzr; 
    void UpdateBD(TArray<Real>& cons, const Grid& grid); 

    // simple copy boundary condition
    static void OutflowCopyXL(TArray<Real>& cons, const Grid& grid); 
    static void OutflowCopyXR(TArray<Real>& cons, const Grid& grid); 
    static void OutflowCopyYL(TArray<Real>& cons, const Grid& grid); 
    static void OutflowCopyYR(TArray<Real>& cons, const Grid& grid); 
    static void OutflowCopyZL(TArray<Real>& cons, const Grid& grid); 
    static void OutflowCopyZR(TArray<Real>& cons, const Grid& grid); 

    static void SelfDefineCopyXL(TArray<Real>& cons, const Grid& grid); 
    static void SelfDefineCopyXR(TArray<Real>& cons, const Grid& grid); 
};


} // namespace Gaukuk
