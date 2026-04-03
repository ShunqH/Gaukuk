#pragma once 

// C++ headers

// Gaukuk dependence
#include "../gaukuk.hpp" 
#include "../template_array.hpp" 
#include "../grid/grid.hpp"

namespace Gaukuk
{
    
class Boundary{
using BoundaryFunc = void (*)(TArray<Real>&, Grid&);
public:
friend class Sim; 
    Boundary(); 
    BoundaryFunc bdxl; 
    BoundaryFunc bdxr; 
    BoundaryFunc bdyl; 
    BoundaryFunc bdyr; 
    BoundaryFunc bdzl; 
    BoundaryFunc bdzr; 
    void UpdateBD(TArray<Real>& cons, Grid& grid); 

    // simple copy boundary condition
    static void OutflowCopyXL(TArray<Real>& cons, Grid& grid); 
    static void OutflowCopyXR(TArray<Real>& cons, Grid& grid); 
    static void OutflowCopyYL(TArray<Real>& cons, Grid& grid); 
    static void OutflowCopyYR(TArray<Real>& cons, Grid& grid); 
    static void OutflowCopyZL(TArray<Real>& cons, Grid& grid); 
    static void OutflowCopyZR(TArray<Real>& cons, Grid& grid); 

    static void SelfDefineCopyXL(TArray<Real>& cons, Grid& grid); 
    static void SelfDefineCopyXR(TArray<Real>& cons, Grid& grid); 
};


} // namespace Gaukuk
