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
    BoundaryFunc bdxr; 
    BoundaryFunc bdxl; 
    BoundaryFunc bdyr; 
    BoundaryFunc bdyl; 
    BoundaryFunc bdzr; 
    BoundaryFunc bdzl; 

    // simple copy boundary condition
    void OutflowCopyXL(TArray<Real>& cons, Grid& grid); 
    void OutflowCopyXR(TArray<Real>& cons, Grid& grid); 
    void OutflowCopyYL(TArray<Real>& cons, Grid& grid); 
    void OutflowCopyYR(TArray<Real>& cons, Grid& grid); 
    void OutflowCopyZL(TArray<Real>& cons, Grid& grid); 
    void OutflowCopyZR(TArray<Real>& cons, Grid& grid); 

    void SelfDefineCopyXL(TArray<Real>& cons, Grid& grid); 
    void SelfDefineCopyXR(TArray<Real>& cons, Grid& grid); 
};


} // namespace Gaukuk
