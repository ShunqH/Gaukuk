// C++ Headers
// #include <iostream>     // std::cout; std::endl; std::cerr
#include <cmath>            // std::sqrt, std::abs
#include <algorithm>        // std::max(), std::min()  
#include <stdexcept>
#include <vector>       

// Gaukuk dependence
#include "../gaukuk.hpp" 
#include "../sim.hpp"
#include "boundary.hpp"
// #include "../utils/utils.hpp"

namespace Gaukuk
{

namespace{

int FindLeftId(const TArray<Real>& xc, int ib, int ie, Real xShell, Real rShell, Real xmin){
    // int il = ib; 
    Real xb = xShell - 1.2*rShell; 
    if (xb<xmin){
        throw std::invalid_argument("immersed boundary out of domain");
    }
    for (int i=ib; i<ie; i++){ 
        Real xnow = xc(i); 
        if (xnow>=xb){
            return i; 
        }
    }
    return ie; 
}

int FindRightId(const TArray<Real>& xc, int ib, int ie, Real xShell, Real rShell, Real xmax){
    // int il = ib; 
    Real xe = xShell + 1.2*rShell; 
    if (xe>xmax){
        throw std::invalid_argument("immersed boundary out of domain");
    }
    for (int i=ie-1; i>=ib; i--){ 
        Real xnow = xc(i); 
        if (xnow<=xe){
            return i; 
        }
    }
    return ib; 
}

bool IsShellCell(Real delx, Real dely, Real delz, Real rShell, Real dxShell, Real dyShell, Real dzShell){
    Real rnow = std::sqrt(delx*delx + dely*dely + delz*delz); 
    if (rnow <rShell){
        rnow=std::sqrt((delx+dxShell)*(delx+dxShell) + dely*dely + delz*delz); 
        if (rnow > rShell) { return true;  }
        rnow=std::sqrt((delx-dxShell)*(delx-dxShell) + dely*dely + delz*delz); 
        if (rnow > rShell) { return true;  }

        rnow=std::sqrt(delx*delx + (dely+dyShell)*(dely+dyShell) + delz*delz); 
        if (rnow > rShell) { return true;  }
        rnow=std::sqrt(delx*delx + (dely-dyShell)*(dely-dyShell) + delz*delz); 
        if (rnow > rShell) { return true;  }

        rnow=std::sqrt(delx*delx + dely*dely + (delz+dzShell)*(delz+dzShell)); 
        if (rnow > rShell) { return true;  }
        rnow=std::sqrt(delx*delx + dely*dely + (delz-dzShell)*(delz-dzShell)); 
        if (rnow > rShell) { return true;  }
    }
    return false; 
}

// transform Cartesian velocity to v_perpendicular and v_shear
// vp: v_perpendicular; vs: v_shear
void CarToLocal(const Real v1, const Real v2, const Real v3, const Real cosp, const Real sinp, 
        const Real cost, const Real sint, Real& vr, Real& vp, Real& vt){
  vr = v1*cosp*sint + v2*sinp*sint + v3*cost;
  vp = - v1*sinp + v2*cosp; 
  vt = v1*cosp*cost + v2*sinp*cost - v3*sint; 
}

// invert of CarToLocal, transform back to Cartesian
void LocalToCar(Real& v1, Real& v2, Real& v3, const Real cosp, const Real sinp, 
        const Real cost, const Real sint, const Real vr, const Real vp, const Real vt){
  v1 = vr*sint*cosp + vt*cost*cosp - vp*sinp; 
  v2 = vr*sint*sinp + vt*cost*sinp + vp*cosp;
  v3 = vr*cost - vt*sint; 
}

}

void Boundary::SetupImmersedShell(const Grid& grid, const Domain& domain, Real xShell, Real yShell, Real zShell, Real rShell){
    if (grid.nx <=1 || grid.ny <=1 ) {
        throw std::invalid_argument("immersed boundary has to be at least 2D. ");
    }
    isImmersed = true; 
    int il = FindLeftId(domain.xc, grid.ib, grid.ie, xShell, rShell, domain.xmin); 
    int ir = FindRightId(domain.xc, grid.ib, grid.ie, xShell, rShell, domain.xmax); 
    int jl = FindLeftId(domain.yc, grid.jb, grid.je, yShell, rShell, domain.ymin); 
    int jr = FindRightId(domain.yc, grid.jb, grid.je, yShell, rShell, domain.ymax); 
    int kl=grid.kb; 
    int kr = grid.ke;
    if (grid.nz>1){
        kl = FindLeftId(domain.zc, grid.kb, grid.ke, zShell, rShell, domain.zmin); 
        kr = FindRightId(domain.zc, grid.kb, grid.ke, zShell, rShell, domain.zmax); 
    }

    if (il > ir || jl > jr || kl > kr) {
        nShell = 0;
        throw std::invalid_argument("immersed boundary out of domain");
        return;
    }
    nShell = 0; 
    std::vector<int> shellIDx;
    std::vector<int> shellIDy;
    std::vector<int> shellIDz;
    #pragma omp parallel
    {
        std::vector<int> shellIDxLocal; // Thread-local storage
        std::vector<int> shellIDyLocal;
        std::vector<int> shellIDzLocal;
        #pragma omp for collapse(2) schedule(static)
        for (int k=kl; k<kr; k++){
            for (int j=jl; j<jr; j++){
                for (int i=il; i<ir; i++){
                    Real delx = domain.xc(i) - xShell; 
                    Real dely = domain.yc(j) - yShell; 
                    Real delz = domain.zc(k) - zShell; 
                    Real dxShell = grid.nGhost*domain.dx; 
                    Real dyShell = grid.nGhost*domain.dy; 
                    Real dzShell = grid.nGhost*domain.dz; 
                    if ( IsShellCell(delx, dely, delz, rShell, dxShell, dyShell, dzShell) ){ 
                        shellIDxLocal.push_back(i); 
                        shellIDyLocal.push_back(j); 
                        shellIDzLocal.push_back(k); 
                    }
                }
            }
        }
        #pragma omp critical
        {
            shellIDx.insert(shellIDx.end(), shellIDxLocal.begin(), shellIDxLocal.end());
            shellIDy.insert(shellIDy.end(), shellIDyLocal.begin(), shellIDyLocal.end());
            shellIDz.insert(shellIDz.end(), shellIDzLocal.begin(), shellIDzLocal.end());
        }
    }
    nShell = shellIDx.size(); 

    shellIDs.NewArray(3, nShell, nTarMax+1); 
    nTars.NewArray(nShell); 
    // targetIDs.NewArray(3, nShell, nTarMax+1); 
    shellDirection.NewArray(4, nShell); 
    disGamma.NewArray(nShell, nTarMax+1); 

    int nxSearch = 2*grid.nGhost, nySearch = 2*grid.nGhost, nzSearch = 0;
    if (grid.nz>1) { nzSearch = 2*grid.nGhost; } 

    #pragma omp for schedule(static)
    for (int iShell=0; iShell<nShell; iShell++){
        int i = shellIDx[iShell]; 
        int j = shellIDy[iShell]; 
        int k = shellIDz[iShell]; 
        shellIDs(0, iShell, 0) = i; 
        shellIDs(1, iShell, 0) = j; 
        shellIDs(2, iShell, 0) = k; 
        Real delx = domain.xc(i) - xShell; 
        Real dely = domain.yc(j) - yShell; 
        Real delz = domain.zc(k) - zShell; 
        Real rnow = std::sqrt(delx*delx + dely*dely + delz*delz); 

        Real cosp, sinp, cost, sint; 
        Real rxy = std::sqrt(delx*delx + dely*dely);
        if (rxy < 1e-15) {
            cosp = 1.0; sinp = 0.0;
        } else {
            cosp = delx / rxy;
            sinp = dely / rxy;
        }
        if (rnow < 1e-15) {
            cost = 1.0; sint = 0.0;
        } else {
            cost = delz / rnow;
            sint = rxy / rnow; 
        }
        shellDirection(0, iShell) = cosp; 
        shellDirection(1, iShell) = sinp; 
        shellDirection(2, iShell) = cost; 
        shellDirection(3, iShell) = sint; 

        Real xImg = sint*cosp*(2*rShell - rnow) + xShell;
        Real yImg = sint*sinp*(2*rShell - rnow) + yShell;
        Real zImg = cost*(2*rShell - rnow) + zShell;

        int nTarNow = 0; 
        Real disDel[nTarMax+1]; 
        disDel[0] = 2*(rShell - rnow); 
        Real xTar, yTar, zTar; 
        for (int ii= -nxSearch; ii < nxSearch+1; ++ii){
            for (int jj= -nySearch; jj < nySearch+1; ++jj){
                for (int kk= -nzSearch; kk < nzSearch+1; ++kk){
                    int iTar = i+ii; 
                    int jTar = j+jj; 
                    int kTar = k+kk; 
                    if (iTar<grid.ib || iTar>=grid.ie || 
                        jTar<grid.jb || jTar>=grid.je || 
                        kTar<grid.kb || kTar>=grid.ke   ){
                            throw std::invalid_argument("immersed body is too close to the domain boundary, neighboor search out of range.");
                        }
                    xTar = domain.xc(iTar);
                    yTar = domain.yc(jTar); 
                    zTar = domain.zc(kTar); 
                    if ( std::abs(xTar - xImg) <= domain.dx && 
                         std::abs(yTar - yImg) <= domain.dy && 
                         std::abs(zTar - zImg) <= domain.dz    ){
                        Real rTar = std::sqrt( (xTar - xShell)*(xTar - xShell) 
                                             + (yTar - yShell)*(yTar - yShell) 
                                             + (zTar - zShell)*(zTar - zShell)); 
                        if ( rTar > rShell){
                            nTarNow++; 
                            if (nTarNow > nTarMax) {
                                // usually won't happen, just a protection
                                throw std::logic_error("n of imagine points is more than nTarMax");
                            }
                            shellIDs(0, iShell, nTarNow) = iTar; 
                            shellIDs(1, iShell, nTarNow) = jTar; 
                            shellIDs(2, iShell, nTarNow) = kTar; 
                            disDel[nTarNow] = std::sqrt( (xTar - xImg)*(xTar - xImg) 
                                                       + (yTar - yImg)*(yTar - yImg) 
                                                       + (zTar - zImg)*(zTar - zImg)); 
                        }
                    }
                    
                }
            }
        }
        nTars(iShell) = nTarNow; 
        // if (nTarNow == 0) {
        //     throw std::logic_error("A shell cell has no image target (nTarNow==0).");
        // }
        Real SumDel = 0.;
        for (int n=0; n<nTarNow+1; n++) SumDel = SumDel + 1./(disDel[n]*disDel[n]) ;
        for (int n=0; n<nTarNow+1; n++) disGamma(iShell, n) = 1./(disDel[n]*disDel[n]*SumDel); 
    }
}

void Boundary::ImmersedReflective(TArray<Real>& cons, const EquationOfState& eos){
    Real dmin = DENSITY_FLOOR;
    Real pmin = PRESSURE_FLOOR;
    Real gm1 = eos.GetGamma() - 1; 
    #pragma omp for schedule(static)
    for (int iShell=0; iShell<nShell; iShell++){
        int iNow = shellIDs(0, iShell, 0); 
        int jNow = shellIDs(1, iShell, 0); 
        int kNow = shellIDs(2, iShell, 0); 
        Real cosp = shellDirection(0, iShell); 
        Real sinp = shellDirection(1, iShell); 
        Real cost = shellDirection(2, iShell); 
        Real sint = shellDirection(3, iShell); 

        Real vrloc, vploc, vtloc;

        Real dens = 0; 
        Real vr = 0; 
        Real vp = 0;
        Real vt = 0;
        Real pres = 0; 
        int nTarNow = nTars(iShell); 
        for (int n=1; n<=nTarNow; n++){
            int iTar = shellIDs(0, iShell, n); 
            int jTar = shellIDs(1, iShell, n); 
            int kTar = shellIDs(2, iShell, n); 
            Real denTar = cons(DEN, kTar, jTar, iTar); 
            denTar = std::max(denTar , dmin); 
            Real denInv = 1.0/denTar;
            dens += 1./(1-disGamma(iShell, 0)) * disGamma(iShell, n)*denTar; 
            Real vx = cons(MTX, kTar, jTar, iTar)*denInv; 
            Real vy = cons(MTY, kTar, jTar, iTar)*denInv; 
            Real vz = cons(MTZ, kTar, jTar, iTar)*denInv; 

            CarToLocal(vx, vy, vz, cosp, sinp, cost, sint, vrloc, vploc, vtloc);

            // vr += 1./(-1.-disGamma(iShell, 0)) * disGamma(iShell, n)*vrloc;
            vr += 1./(1.-disGamma(iShell, 0)) * disGamma(iShell, n)*vrloc;
            vp += 1./(1.-disGamma(iShell, 0)) * disGamma(iShell, n)*vploc;
            vt += 1./(1.-disGamma(iShell, 0)) * disGamma(iShell, n)*vtloc;

            Real engTar = cons(ENG, kTar, jTar, iTar); 
            Real engKin = 0.5 * denTar*( vx*vx + vy*vy + vz*vz ); 
            Real preTar = gm1 * (engTar - engKin); 
            preTar = (preTar > pmin) ? preTar : pmin;
            pres += 1./(1.-disGamma(iShell, 0)) * disGamma(iShell, n)*preTar; 
        }

        vr = std::min(0.0, vr); 

        Real vx, vy, vz; 
        LocalToCar(vx, vy, vz, cosp, sinp, cost, sint, vr, vp, vt); 
        cons(DEN, kNow, jNow, iNow) = dens; 
        cons(MTX, kNow, jNow, iNow) = dens*vx; 
        cons(MTY, kNow, jNow, iNow) = dens*vy; 
        cons(MTZ, kNow, jNow, iNow) = dens*vz; 
        Real engKin = 0.5 * dens * ( vx*vx + vy*vy + vz*vz ); 
        Real engGas = eos.EGas(dens, pres); 
        cons(ENG, kNow, jNow, iNow) = engKin + engGas; 
    }
}

}