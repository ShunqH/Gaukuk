#include "sim.hpp"

namespace Gaukuk{

void Domain::SetDomain(const int nx, const int ny, const int nz){
    dx = ( xmax - xmin ) / nx; 
    dy = ( ymax - ymin ) / ny; 
    dz = ( zmax - zmin ) / nz; 
    xGrid.NewArray(nx); 
    yGrid.NewArray(ny); 
    zGrid.NewArray(nz); 
    for (int i=0; i<nx; i++) { xGrid(i) = xmin + 0.5*dx + i*dx; }
    for (int j=0; j<ny; j++) { yGrid(j) = ymin + 0.5*dy + j*dy; }
    for (int k=0; k<nz; k++) { zGrid(k) = zmin + 0.5*dz + k*dz; }
}

Sim::Sim(){
    domain.xmax = 0.5; 
    domain.xmin = -0.5; 
    domain.ymax = 0.5; 
    domain.ymin = -0.5; 
    domain.zmax = 0.5; 
    domain.zmin = -0.5; 
    t = 0;
    CFL = 0.5; 
    
    nVar = 5;
    nx = 1001; 
    ny = 101; 
    nz = 101; 
    nGhost = 1; 

    lenx = nx + 2*nGhost; 
    leny = ny + 2*nGhost; 
    lenz = nz + 2*nGhost; 
    lenArr = lenx*leny*lenz; 

    cmax = 1; 

    domain.SetDomain(nx, ny, nz); 
    cons.NewArray(nVar, lenz, leny, lenx);
    prim.NewArray(nVar, lenz, leny, lenx);
    flx1.NewArray(nVar, nz, ny, nx+1); 
    flx2.NewArray(nVar, nz, ny+1, nx); 
    flx3.NewArray(nVar, nz+1, ny, nx); 
}

} // namespace Gaukuk