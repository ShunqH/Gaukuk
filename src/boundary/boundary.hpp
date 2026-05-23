#pragma once 

// C++ headers
#include <functional>

// Gaukuk dependence
#include "../gaukuk.hpp" 
#include "../template_array.hpp" 
#include "../grid.hpp"
#include "../eos/eos.hpp"

namespace Gaukuk
{
    
class Boundary{
// using BoundaryFunc = void (*)(TArray<Real>&, const Grid&, const EquationOfState& eos);
using BoundaryFunc = std::function<void(TArray<Real>&, const Grid&, const EquationOfState&)>;
public:
// friend class Sim; 
    Boundary(); 
    BoundaryFunc Bdxl; 
    BoundaryFunc Bdxr; 
    BoundaryFunc Bdyl; 
    BoundaryFunc Bdyr; 
    BoundaryFunc Bdzl; 
    BoundaryFunc Bdzr; 
    void UpdateBD(TArray<Real>& cons, const Grid& grid, const Domain& domain, const EquationOfState& eos); 
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

    void SetupImmersedShell(const Grid& grid, const Domain& domain, 
                            Real xShell, Real yShell, Real zShell, Real rShell); 
    void EnrollFixedBoundaryData(const TArray<Real>& cons, const Grid& grid, bool fixInner = false); 

private:
    // simple copy boundary condition
    static void OutflowCopyXL(TArray<Real>& cons, const Grid& grid, const EquationOfState& eos); 
    static void OutflowCopyXR(TArray<Real>& cons, const Grid& grid, const EquationOfState& eos); 
    static void OutflowCopyYL(TArray<Real>& cons, const Grid& grid, const EquationOfState& eos); 
    static void OutflowCopyYR(TArray<Real>& cons, const Grid& grid, const EquationOfState& eos); 
    static void OutflowCopyZL(TArray<Real>& cons, const Grid& grid, const EquationOfState& eos); 
    static void OutflowCopyZR(TArray<Real>& cons, const Grid& grid, const EquationOfState& eos); 

    // periodic boundary condition
    static void PeriodicXL(TArray<Real>& cons, const Grid& grid, const EquationOfState& eos); 
    static void PeriodicXR(TArray<Real>& cons, const Grid& grid, const EquationOfState& eos); 
    static void PeriodicYL(TArray<Real>& cons, const Grid& grid, const EquationOfState& eos); 
    static void PeriodicYR(TArray<Real>& cons, const Grid& grid, const EquationOfState& eos); 
    static void PeriodicZL(TArray<Real>& cons, const Grid& grid, const EquationOfState& eos); 
    static void PeriodicZR(TArray<Real>& cons, const Grid& grid, const EquationOfState& eos); 

    // reflective boundary condition
    static void ReflectiveXL(TArray<Real>& cons, const Grid& grid, const EquationOfState& eos); 
    static void ReflectiveXR(TArray<Real>& cons, const Grid& grid, const EquationOfState& eos); 
    static void ReflectiveYL(TArray<Real>& cons, const Grid& grid, const EquationOfState& eos); 
    static void ReflectiveYR(TArray<Real>& cons, const Grid& grid, const EquationOfState& eos); 
    static void ReflectiveZL(TArray<Real>& cons, const Grid& grid, const EquationOfState& eos); 
    static void ReflectiveZR(TArray<Real>& cons, const Grid& grid, const EquationOfState& eos); 

    // outflow boundary condition
    static void OutflowXL(TArray<Real>& cons, const Grid& grid, const EquationOfState& eos); 
    static void OutflowXR(TArray<Real>& cons, const Grid& grid, const EquationOfState& eos); 
    static void OutflowYL(TArray<Real>& cons, const Grid& grid, const EquationOfState& eos); 
    static void OutflowYR(TArray<Real>& cons, const Grid& grid, const EquationOfState& eos); 
    static void OutflowZL(TArray<Real>& cons, const Grid& grid, const EquationOfState& eos); 
    static void OutflowZR(TArray<Real>& cons, const Grid& grid, const EquationOfState& eos); 

    // user define boundary condition 
    static void SelfDefineBDXL(TArray<Real>& cons, const Grid& grid, const EquationOfState& eos); 
    static void SelfDefineBDXR(TArray<Real>& cons, const Grid& grid, const EquationOfState& eos); 
    static void SelfDefineBDYL(TArray<Real>& cons, const Grid& grid, const EquationOfState& eos); 
    static void SelfDefineBDYR(TArray<Real>& cons, const Grid& grid, const EquationOfState& eos); 
    static void SelfDefineBDZL(TArray<Real>& cons, const Grid& grid, const EquationOfState& eos); 
    static void SelfDefineBDZR(TArray<Real>& cons, const Grid& grid, const EquationOfState& eos);
    
    // user define boundary condition 
    void FixedBDXL(TArray<Real>& cons, const Grid& grid, const EquationOfState& eos); 
    void FixedBDXR(TArray<Real>& cons, const Grid& grid, const EquationOfState& eos); 
    void FixedBDYL(TArray<Real>& cons, const Grid& grid, const EquationOfState& eos); 
    void FixedBDYR(TArray<Real>& cons, const Grid& grid, const EquationOfState& eos); 
    void FixedBDZL(TArray<Real>& cons, const Grid& grid, const EquationOfState& eos); 
    void FixedBDZR(TArray<Real>& cons, const Grid& grid, const EquationOfState& eos);
    void FixedBDInner(TArray<Real>& cons, const Grid& grid, const Domain& domain, const EquationOfState& eos);
    // static/fixed boundary condition 
    bool isFixedBdEnrolled, isFixedInner; 
    TArray<Real> consFixedInner; 
    TArray<Real> consFixedBdxl, consFixedBdxr, consFixedBdyl, consFixedBdyr, consFixedBdzl, consFixedBdzr; 

    // immersed reflective boundary; solid body
    void ImmersedReflective(TArray<Real>& cons, const EquationOfState& eos); 
    // for immersed boundary condition, like a solid body
    bool isImmersed; 
    const int nTarMax = 8; 
    int nShell; 
    // if there are more than one body, I can expand the array here, 
    // because the cons update function should be the same to all immersed bd cells
    TArray<int> shellIDs, nTars; 
    TArray<Real> shellDirection, disGamma; 

    
};


} // namespace Gaukuk
