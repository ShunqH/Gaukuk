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
    }else if (typeBDXL == 2){
        Bdxl = &ReflectiveXL;
    }else if (typeBDXL == 3){
        Bdxl = &OutflowXL;
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
    }else{
        Bdzr = &OutflowCopyYR;
    }
}

void Boundary::UpdateBD(TArray<Real>& cons, const Grid& grid){
    Bdxl(cons, grid); 
    Bdxr(cons, grid); 
    Bdyl(cons, grid); 
    Bdyr(cons, grid); 
if (grid.nz>1){
    Bdzl(cons, grid); 
    Bdzr(cons, grid); 
}
}

} // namespace Gaukuk