#pragma once 

// C++ headers

// Gaukuk dependence
#include "../gaukuk.hpp" 
#include "../template_array.hpp" 
#include "../grid.hpp"

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
    void EnrollSelfDefineBDXL(){
        Bdxl = &SelfDefineBDXL;
    } 
    void EnrollSelfDefineBDXR(){
        Bdxr = &SelfDefineBDXR;
    } 
    void EnrollSelfDefineBDYL(){
        Bdyl = &SelfDefineBDYL;
    } 
    void EnrollSelfDefineBDYR(){
        Bdyr = &SelfDefineBDYR;
    } 
    void EnrollSelfDefineBDZL(){
        Bdzl = &SelfDefineBDZL;
    } 
    void EnrollSelfDefineBDZR(){
        Bdzr = &SelfDefineBDZR;
    } 


    // simple copy boundary condition
    static void OutflowCopyXL(TArray<Real>& cons, const Grid& grid); 
    static void OutflowCopyXR(TArray<Real>& cons, const Grid& grid); 
    static void OutflowCopyYL(TArray<Real>& cons, const Grid& grid); 
    static void OutflowCopyYR(TArray<Real>& cons, const Grid& grid); 
    static void OutflowCopyZL(TArray<Real>& cons, const Grid& grid); 
    static void OutflowCopyZR(TArray<Real>& cons, const Grid& grid); 

    // periodic boundary condition
    static void PeriodicXL(TArray<Real>& cons, const Grid& grid); 
    static void PeriodicXR(TArray<Real>& cons, const Grid& grid); 
    static void PeriodicYL(TArray<Real>& cons, const Grid& grid); 
    static void PeriodicYR(TArray<Real>& cons, const Grid& grid); 
    static void PeriodicZL(TArray<Real>& cons, const Grid& grid); 
    static void PeriodicZR(TArray<Real>& cons, const Grid& grid); 

    // reflective boundary condition
    static void ReflectiveXL(TArray<Real>& cons, const Grid& grid); 
    static void ReflectiveXR(TArray<Real>& cons, const Grid& grid); 
    static void ReflectiveYL(TArray<Real>& cons, const Grid& grid); 
    static void ReflectiveYR(TArray<Real>& cons, const Grid& grid); 
    static void ReflectiveZL(TArray<Real>& cons, const Grid& grid); 
    static void ReflectiveZR(TArray<Real>& cons, const Grid& grid); 

    static void SelfDefineBDXL(TArray<Real>& cons, const Grid& grid); 
    static void SelfDefineBDXR(TArray<Real>& cons, const Grid& grid); 
    static void SelfDefineBDYL(TArray<Real>& cons, const Grid& grid); 
    static void SelfDefineBDYR(TArray<Real>& cons, const Grid& grid); 
    static void SelfDefineBDZL(TArray<Real>& cons, const Grid& grid); 
    static void SelfDefineBDZR(TArray<Real>& cons, const Grid& grid); 
};


} // namespace Gaukuk
