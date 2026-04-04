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
        bdxr = &OutflowCopyXL;
    }else{
        bdxr = &OutflowCopyXL;
    }

    // X right boundary registration
    if (typeBDXR == 0){
        bdxr = &OutflowCopyXR;
    }else{
        bdxr = &OutflowCopyXR;
    }

    // y left boundary registration
    if (typeBDYL == 0){
        bdyr = &OutflowCopyYL;
    }else{
        bdyr = &OutflowCopyYL;
    }

    // y right boundary registration
    if (typeBDYR == 0){
        bdyr = &OutflowCopyYR;
    }else{
        bdyr = &OutflowCopyYR;
    }

    // z left boundary registration
    if (typeBDYL == 0){
        bdyr = &OutflowCopyYL;
    }else{
        bdyr = &OutflowCopyYL;
    }

    // z right boundary registration
    if (typeBDYR == 0){
        bdyr = &OutflowCopyYR;
    }else{
        bdyr = &OutflowCopyYR;
    }
}

void Boundary::UpdateBD(TArray<Real>& cons, const Grid& grid){
    bdxl(cons, grid); 
    bdxr(cons, grid); 
    bdyl(cons, grid); 
    bdyr(cons, grid); 
    bdzl(cons, grid); 
    bdzr(cons, grid); 
}

} // namespace Gaukuk