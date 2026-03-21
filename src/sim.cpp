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
    for (int j=0; j<nx; j++) { yGrid(j) = ymin + 0.5*dy + j*dy; }
    for (int k=0; k<nx; k++) { zGrid(k) = zmin + 0.5*dz + k*dz; }
}

Gaukuk::Gaukuk(){
    domain.xmax = 0.5; 
    domain.xmin = -0.5; 
    domain.ymax = 0.5; 
    domain.ymin = -0.5; 
    domain.zmax = 0.5; 
    domain.zmin = -0.5; 
    t = 0;
    CFL = 0.5; 
    
    nx = 101; 
    ny = 1; 
    nz = 1; 
    nGhost = 1; 

    lenx = nx + 2*nGhost; 
    leny = ny + 2*nGhost; 
    lenz = nz + 2*nGhost; 
    lenArr = lenx*leny*lenz; 

    cmax = 1; 

    domain.SetDomain(nx, ny, nz); 
    cons.NewArray(5, lenz, leny, lenx);
    prim.NewArray(5, lenz, leny, lenx);
    flx1.NewArray(5, nz, ny, nx+1); 
    flx2.NewArray(5, nz, ny+1, nx); 
    flx3.NewArray(5, nz+1, ny, nx); 
}

} // namespace Gaukuk