// CPP header
#include <stdexcept>    // runtime_error
#include <algorithm>    // std::min()
#include <iostream>     // std::cout; std::endl; std::cerr

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
    dxRec = 1./dx; 
    dyRec = 1./dy; 
    dzRec = 1./dz; 
    if (grid.ny<=1){
        dyRec = 0; 
    }
    if (grid.nz<=1){
        dzRec = 0; 
    }
    drmin = std::min(std::min(dx, dy), dz); 
    xc.NewArray(grid.nx); 
    yc.NewArray(grid.ny); 
    zc.NewArray(grid.nz); 
    for (int i=0; i<grid.nx; i++) { xc(i) = xmin + 0.5*dx + i*dx; }
    for (int j=0; j<grid.ny; j++) { yc(j) = ymin + 0.5*dy + j*dy; }
    for (int k=0; k<grid.nz; k++) { zc(k) = zmin + 0.5*dz + k*dz; }
}


Sim::Sim(): domain(grid), flux(grid.lenx), 
            step(0), t(0), dt(1e10), dtUntilOutput(1e10), cmax(1e-16){
    CFL = Config::getInstance().get("CFL"); 

    cons.NewArray(NVar, grid.lenz, grid.leny, grid.lenx);
    prim.NewArray(NVar, grid.lenz, grid.leny, grid.lenx);
    flx1.NewArray(NVar, grid.lenz, grid.leny, grid.lenx+1); 
    flx2.NewArray(NVar, grid.lenz, grid.leny+1, grid.lenx); 
    flx3.NewArray(NVar, grid.lenz+1, grid.leny, grid.lenx); 
    
    integratorType = static_cast<int>(Config::getInstance().get("integrator", 2)); 
    rcOrder = static_cast<int>(Config::getInstance().get("RcOrder", 1)); 
    
    if (integratorType == 1) {
        integrator_ = &Sim::ForwardEuler_; 
    }else if (integratorType == 3) {
        consTemp_.NewArray(NVar, grid.lenz, grid.leny, grid.lenx); 
        integrator_ = &Sim::RK3_; 
    }else {
        consTemp_.NewArray(NVar, grid.lenz, grid.leny, grid.lenx); 
        integrator_ = &Sim::RK2_; 
    }
}

void Sim::Advance(Real dtoutput){
    Real tNext = t + dtoutput; 
    dtUntilOutput = dtoutput; 
    while (t<tNext){
        (this->*integrator_)(); 
        t += dt; 
        step ++; 
        dtUntilOutput = tNext - t; 
        std::cout << "step = " << step 
                  << ", t = " << t 
                  << ", dt = " << dt 
                  << std::endl; 
    }
}

} // namespace Gaukuk