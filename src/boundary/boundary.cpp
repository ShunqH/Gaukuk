// C++ Headers
#include <iostream>     // std::cout; std::endl; std::cerr

// Gaukuk dependence
#include "../gaukuk.hpp" 
#include "../sim.hpp"
#include "boundary.hpp"
#include "../utils/utils.hpp"

namespace Gaukuk
{

Boundary::Boundary() : isFixedBdEnrolled(false), isFixedInner(false), isImmersed(false), nShell(0) {
    int typeBDXL = static_cast<int>(Config::getInstance().get("xleft")); 
    int typeBDXR = static_cast<int>(Config::getInstance().get("xright")); 
    int typeBDYL = static_cast<int>(Config::getInstance().get("yleft")); 
    int typeBDYR = static_cast<int>(Config::getInstance().get("yright")); 
    int typeBDZL = static_cast<int>(Config::getInstance().get("zleft")); 
    int typeBDZR = static_cast<int>(Config::getInstance().get("zright"));
    int enableImmersed = static_cast<int>(Config::getInstance().get("enableImmersed", 0));

    // X left boundary registration
    if (typeBDXL == 0){
        Bdxl = &OutflowCopyXL;
    }else if (typeBDXL == 1){
        Bdxl = &PeriodicXL;
    }else if (typeBDXL == 2){
        Bdxl = &ReflectiveXL;
    }else if (typeBDXL == 3){
        Bdxl = &OutflowXL;
    }else if (typeBDXL == 4){
        Bdxl = [this](TArray<Real>& cons, const Grid& grid, const EquationOfState& eos) {
            this->FixedBDXL(cons, grid, eos);
        };
    }else{
        Bdxl = &OutflowCopyXL;
    }

    // X right boundary registration
    if (typeBDXR == 0){
        Bdxr = &OutflowCopyXR;
    }else if (typeBDXR == 1){
        Bdxr = &PeriodicXR;
    }else if (typeBDXR == 2){
        Bdxr = &ReflectiveXR;
    }else if (typeBDXR == 3){
        Bdxr = &OutflowXR;
    }else if (typeBDXR == 4){
        Bdxr = [this](TArray<Real>& cons, const Grid& grid, const EquationOfState& eos) {
            this->FixedBDXR(cons, grid, eos);
        };
    }else{
        Bdxr = &OutflowCopyXR;
    }

    // y left boundary registration
    if (typeBDYL == 0){
        Bdyl = &OutflowCopyYL;
    }else if (typeBDYL == 1){
        Bdyl = &PeriodicYL;
    }else if (typeBDYL == 2){
        Bdyl = &ReflectiveYL;
    }else if (typeBDYL == 3){
        Bdyl = &OutflowYL;
    }else if (typeBDYL == 4){
        Bdyl = [this](TArray<Real>& cons, const Grid& grid, const EquationOfState& eos) {
            this->FixedBDYL(cons, grid, eos);
        };
    }else{
        Bdyl = &OutflowCopyYL;
    }

    // y right boundary registration
    if (typeBDYR == 0){
        Bdyr = &OutflowCopyYR;
    }else if (typeBDYR == 1){
        Bdyr = &PeriodicYR;
    }else if (typeBDYR == 2){
        Bdyr = &ReflectiveYR;
    }else if (typeBDYR == 3){
        Bdyr = &OutflowYR;
    }else if (typeBDYR == 4){
        Bdyr = [this](TArray<Real>& cons, const Grid& grid, const EquationOfState& eos) {
            this->FixedBDYR(cons, grid, eos);
        };
    }else{
        Bdyr = &OutflowCopyYR;
    }

    // z left boundary registration
    if (typeBDZL == 0){
        Bdzl = &OutflowCopyYL;
    }else if (typeBDZL == 1){
        Bdzl = &PeriodicZL;
    }else if (typeBDZL == 2){
        Bdzl = &ReflectiveZL;
    }else if (typeBDZL == 3){
        Bdzl = &OutflowZL;
    }else if (typeBDZL == 4){
        Bdzl = [this](TArray<Real>& cons, const Grid& grid, const EquationOfState& eos) {
            this->FixedBDZL(cons, grid, eos);
        };
    }else{
        Bdzl = &OutflowCopyYL;
    }

    // z right boundary registration
    if (typeBDZR == 0){
        Bdzr = &OutflowCopyYR;
    }else if (typeBDZR == 1){
        Bdzr = &PeriodicZR;
    }else if (typeBDZR == 2){
        Bdzr = &ReflectiveZR;
    }else if (typeBDZR == 3){
        Bdzr = &OutflowZR;
    }else if (typeBDZR == 4){
        Bdzr = [this](TArray<Real>& cons, const Grid& grid, const EquationOfState& eos) {
            this->FixedBDZR(cons, grid, eos);
        };
    }else{
        Bdzr = &OutflowCopyYR;
    }

    if (enableImmersed != 0){ isImmersed = true; }
}

void Boundary::UpdateBD(TArray<Real>& cons, const Grid& grid, const Domain& domain, const EquationOfState& eos){
    Bdxl(cons, grid, eos); 
    Bdxr(cons, grid, eos); 
    Bdyl(cons, grid, eos); 
    Bdyr(cons, grid, eos); 
if (grid.nz>1){
    Bdzl(cons, grid, eos); 
    Bdzr(cons, grid, eos); 
}
if (isImmersed){
    ImmersedReflective(cons, eos); 
}
if (isFixedInner){
    FixedBDInner(cons, grid, domain, eos); 
}
}

} // namespace Gaukuk