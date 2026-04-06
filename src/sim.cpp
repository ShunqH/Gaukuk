// CPP header
#include <stdexcept>    // runtime_error
#include <cmath>        // std::min()

#include "sim.hpp"
#include "utils/utils.hpp"

namespace Gaukuk{

Domain::Domain(const Grid& grid) : 
    xmin(Config::getInstance().get("xmin")), 
    xmax(Config::getInstance().get("xmax")), 
    ymin(Config::getInstance().get("ymin")), 
    ymax(Config::getInstance().get("ymax")), 
    zmin(Config::getInstance().get("zmin")), 
    zmax(Config::getInstance().get("zmax")) {
    dx = ( xmax - xmin ) / grid.nx; 
    dy = ( ymax - ymin ) / grid.ny; 
    dz = ( zmax - zmin ) / grid.nz; 
    drmin = std::min(dx, dy, dz); 
    xGrid.NewArray(grid.nx); 
    yGrid.NewArray(grid.ny); 
    zGrid.NewArray(grid.nz); 
    for (int i=0; i<grid.nx; i++) { xGrid(i) = xmin + 0.5*dx + i*dx; }
    for (int j=0; j<grid.ny; j++) { yGrid(j) = ymin + 0.5*dy + j*dy; }
    for (int k=0; k<grid.nz; k++) { zGrid(k) = zmin + 0.5*dz + k*dz; }
}


Sim::Sim(): domain(grid), flux(grid.lenx){
    t = 0;
    CFL = Config::getInstance().get("CFL"); 
    cmax = 1; 

    cons.NewArray(NVar, grid.lenz, grid.leny, grid.lenx);
    prim.NewArray(NVar, grid.lenz, grid.leny, grid.lenx);
    flx1.NewArray(NVar, grid.lenz, grid.leny, grid.lenx+1); 
    flx2.NewArray(NVar, grid.lenz, grid.leny+1, grid.lenx); 
    flx3.NewArray(NVar, grid.lenz+1, grid.leny, grid.lenx); 

    if (static_cast<int>(Config::getInstance().get("rk_order")) >= 2) {
        consTemp1_.NewArray(NVar, grid.lenz, grid.leny, grid.lenx); 
    }
    if (static_cast<int>(Config::getInstance().get("rk_order")) >= 3) {
        consTemp2_.NewArray(NVar, grid.lenz, grid.leny, grid.lenx); 
    }
}

} // namespace Gaukuk