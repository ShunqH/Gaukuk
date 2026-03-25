// CPP header
#include <stdexcept>    // runtime_error

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
    xGrid.NewArray(grid.nx); 
    yGrid.NewArray(grid.ny); 
    zGrid.NewArray(grid.nz); 
    for (int i=0; i<grid.nx; i++) { xGrid(i) = xmin + 0.5*dx + i*dx; }
    for (int j=0; j<grid.ny; j++) { yGrid(j) = ymin + 0.5*dy + j*dy; }
    for (int k=0; k<grid.nz; k++) { zGrid(k) = zmin + 0.5*dz + k*dz; }
}


Sim::Sim(): domain(grid), nVar(5){
    t = 0;
    CFL = Config::getInstance().get("CFL"); 
    cmax = 1; 

    cons.NewArray(nVar, grid.lenz, grid.leny, grid.lenx);
    prim.NewArray(nVar, grid.lenz, grid.leny, grid.lenx);
    flx1.NewArray(nVar, grid.lenz, grid.leny, grid.lenx+1); 
    flx2.NewArray(nVar, grid.lenz, grid.leny+1, grid.lenx); 
    flx3.NewArray(nVar, grid.lenz+1, grid.leny, grid.lenx); 

    ul_.NewArray(nVar, grid.lenx+1);
    ur_.NewArray(nVar, grid.lenx+1);  
}

} // namespace Gaukuk