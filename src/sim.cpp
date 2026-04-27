// CPP header
#include <stdexcept>    // runtime_error
#include <algorithm>    // std::min()
#include <iostream>     // std::cout; std::endl; std::cerr

// #include <chrono>       // std::chrono::high_resolution_clock
// #include <ctime>        // clock_t
// #include <iomanip>      // std::setw()

// Gaukuk dependence
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
    xc.NewArray(grid.lenx); 
    yc.NewArray(grid.leny); 
    zc.NewArray(grid.lenz); 
    for (int i=0; i<grid.lenx; i++) { xc(i) = xmin - grid.nGhost*dx + 0.5*dx + i*dx; }
    for (int j=0; j<grid.leny; j++) { yc(j) = ymin - grid.nGhost*dy + 0.5*dy + j*dy; }
    for (int k=0; k<grid.lenz; k++) { zc(k) = zmin - grid.nGhost*dz + 0.5*dz + k*dz; }
}


Sim::Sim(): isContinue(true), domain(grid), flux(grid.lenx), 
            step(0), t(0), dt(1e10), dtUntilOutput(1e10), cmax(1e-16){
    CFL = Config::getInstance().get("CFL"); 

    cons.NewArray(NVar, grid.lenz, grid.leny, grid.lenx);
    prim.NewArray(NVar, grid.lenz, grid.leny, grid.lenx);
    flx1.NewArray(NVar, grid.lenz, grid.leny, grid.lenx+1); 
    flx2.NewArray(NVar, grid.lenz, grid.leny+1, grid.lenx); 
    flx3.NewArray(NVar, grid.lenz+1, grid.leny, grid.lenx); 
    consTemp.NewArray(NVar, grid.lenz, grid.leny, grid.lenx); 
    
    integratorType = static_cast<int>(Config::getInstance().get("integrator", 2)); 
    rcOrder = static_cast<int>(Config::getInstance().get("RcOrder", 1));
    stepMax = Config::getInstance().get("stepmax", -1); 
    
    if (integratorType == 1) {
        hydroIntegrator_ = &Sim::ForwardEuler_; 
    }else if (integratorType == 3) {
        hydroIntegrator_ = &Sim::RK3_; 
    }else {
        hydroIntegrator_ = &Sim::RK2_; 
    }
}

void Sim::Advance(Real dtoutput){
    Real tNext = t + dtoutput; 
    dtUntilOutput = dtoutput; 
    while (std::abs(tNext - t) > 1e-12 && isContinue){

        // auto time0 = std::chrono::high_resolution_clock::now();
        // clock_t cputime0 = clock(); 

        cmax = 1e-16;
        eos.CalCmax(cons, grid, cmax); 
        dt = CFL * domain.drmin / cmax; 
        // dt = std::max(dt, 1e-7);
        dt = std::min(dt, dtUntilOutput); 

        // Strang Splitting 
        // source(0.5*dt) -> hydro(dt) -> source(0.5*dt)
        if (srcTerm.sourceEnrolled){
            srcTerm.UpdateSource(cons, t+0.5*dt, 0.5*dt, grid, domain); 
        }
        (this->*hydroIntegrator_)(); 
        if (srcTerm.sourceEnrolled){
            srcTerm.UpdateSource(cons, t + dt, 0.5*dt, grid, domain); 
        }

        t += dt; 
        step ++; 
        dtUntilOutput = tNext - t; 
        if (stepMax>0 && step>=stepMax){
            isContinue = false; 
        }

        std::cout << "step = " << step 
                  << ", t = " << t 
                  << ", dt = " << dt 
                  << std::endl; 

        // auto time1 = std::chrono::high_resolution_clock::now();
        // clock_t cputime1 = clock();
        // Real walltime = std::chrono::duration<Real>(time1 - time0).count(); 
        // Real cputime = double(cputime1 - cputime0) / CLOCKS_PER_SEC; 
        // std::cout << std::right
        //         << "step time: " 
        //         << "\n    walltime = " << std::setw(20) << std::fixed << std::setprecision(15) << walltime
        //         << " s, \n    cputime = " << std::setw(20) << std::fixed << std::setprecision(15) << cputime << " s"
        //         << "\n-----------------------------------------------------------"
        //         << std::endl;
        // std::cout<<cmax<<std::endl;
    }
}

} // namespace Gaukuk