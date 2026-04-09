// C++ Headers
#include <iostream>     // std::cout; std::endl; std::cerr

// Gaukuk dependence
#include "../gaukuk.hpp" 
#include "../sim.hpp"
#include "boundary.hpp"
#include "../utils/utils.hpp"

namespace Gaukuk
{

Boundary::Boundary(){
    int typeBDXL = static_cast<int>(Config::getInstance().get("xleft")); 
    int typeBDXR = static_cast<int>(Config::getInstance().get("xright")); 
    int typeBDYL = static_cast<int>(Config::getInstance().get("yleft")); 
    int typeBDYR = static_cast<int>(Config::getInstance().get("yright")); 
    int typeBDZL = static_cast<int>(Config::getInstance().get("zleft")); 
    int typeBDZR = static_cast<int>(Config::getInstance().get("zright")); 
    // X left boundary registration
    if (typeBDXL == 0){
        Bdxl = &OutflowCopyXL;
    }else if (typeBDXL == 1){
        Bdxl = &PeriodicXL;
    }else{
        Bdxl = &OutflowCopyXL;
    }

    // X right boundary registration
    if (typeBDXR == 0){
        Bdxr = &OutflowCopyXR;
    }else if (typeBDXR == 1){
        Bdxr = &PeriodicXR;
    }else{
        Bdxr = &OutflowCopyXR;
    }

    // y left boundary registration
    if (typeBDYL == 0){
        Bdyl = &OutflowCopyYL;
    }else if (typeBDYL == 1){
        Bdyl = &PeriodicYL;
    }else{
        Bdyl = &OutflowCopyYL;
    }

    // y right boundary registration
    if (typeBDYR == 0){
        Bdyr = &OutflowCopyYR;
    }else if (typeBDYR == 1){
        Bdyr = &PeriodicYR;
    }else{
        Bdyr = &OutflowCopyYR;
    }

    // z left boundary registration
    if (typeBDZL == 0){
        Bdzl = &OutflowCopyYL;
    }else if (typeBDZL == 1){
        Bdzl = &PeriodicZL;
    }else{
        Bdzl = &OutflowCopyYL;
    }

    // z right boundary registration
    if (typeBDZR == 0){
        Bdzr = &OutflowCopyYR;
    }else if (typeBDZR == 1){
        Bdzr = &PeriodicZR;
    }else{
        Bdzr = &OutflowCopyYR;
    }
}

void Boundary::UpdateBD(TArray<Real>& cons, const Grid& grid){
    Bdxl(cons, grid); 
    Bdxr(cons, grid); 
    Bdyl(cons, grid); 
    Bdyr(cons, grid); 
    Bdzl(cons, grid); 
    Bdzr(cons, grid); 
}

} // namespace Gaukuk