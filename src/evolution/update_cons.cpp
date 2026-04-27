// C++ Header
#include <algorithm>        //std::min 

// Gaukuk dependence
#include "../gaukuk.hpp" 
#include "../sim.hpp"

namespace Gaukuk
{

void Sim::UpdateCons(const TArray<Real>& cons_, TArray<Real>& consTemp_, const Real coef1, const Real coef2){
    int il = grid.ib; 
    int ir = grid.ie; 
    int jl = grid.jb; 
    int jr = grid.je; 
    int kl = grid.kb;
    int kr = grid.ke; 

    Real dtdx = dt*domain.dxRec; 
    Real dtdy = dt*domain.dyRec; 
    Real dtdz = dt*domain.dzRec; 

    const Real dmin = DENSITY_FLOOR;
    const Real pmin = PRESSURE_FLOOR;
    const Real gm1 = eos.GetGamma() - 1; 
if (grid.nz == 1){ 
    int k = kl; 
    #pragma omp parallel for schedule(static)      
    for (int j = jl; j < jr; ++j) {
        #pragma omp simd
        for (int i = il; i < ir; ++i) {
            const Real u1 = cons_(DEN, k, j, i);
            const Real u2 = cons_(MTX, k, j, i);
            const Real u3 = cons_(MTY, k, j, i);
            const Real u4 = cons_(MTZ, k, j, i);
            const Real u5 = cons_(ENG, k, j, i);

            const Real dF1 = (flx1(DEN, k, j, i+1) - flx1(DEN, k, j, i)) * dtdx
                           + (flx2(DEN, k, j+1, i) - flx2(DEN, k, j, i)) * dtdy;
            const Real dF2 = (flx1(MTX, k, j, i+1) - flx1(MTX, k, j, i)) * dtdx
                           + (flx2(MTX, k, j+1, i) - flx2(MTX, k, j, i)) * dtdy;
            const Real dF3 = (flx1(MTY, k, j, i+1) - flx1(MTY, k, j, i)) * dtdx
                           + (flx2(MTY, k, j+1, i) - flx2(MTY, k, j, i)) * dtdy;
            const Real dF4 = (flx1(MTZ, k, j, i+1) - flx1(MTZ, k, j, i)) * dtdx
                           + (flx2(MTZ, k, j+1, i) - flx2(MTZ, k, j, i)) * dtdy;
            const Real dF5 = (flx1(ENG, k, j, i+1) - flx1(ENG, k, j, i)) * dtdx
                           + (flx2(ENG, k, j+1, i) - flx2(ENG, k, j, i)) * dtdy;

            Real den = coef1*u1 - coef2*dF1;
            Real mx  = coef1*u2 - coef2*dF2;
            Real my  = coef1*u3 - coef2*dF3;
            Real mz  = coef1*u4 - coef2*dF4;
            Real eng = coef1*u5 - coef2*dF5;

            if (den < dmin) {
                den = dmin;
                mx  = 0.0;
                my  = 0.0;
                mz  = 0.0;
                eng = pmin / gm1 + 0.0;   // ke = 0
            } else {
                Real inv_den = 1.0 / den;
                Real ke = 0.5 * inv_den * (mx*mx + my*my + mz*mz);
                Real press = gm1 * (eng - ke);
                if (press < pmin) {
                    eng = pmin / gm1 + ke;   
                }
            }
            consTemp_(DEN, k, j, i) = den;
            consTemp_(MTX, k, j, i) = mx;
            consTemp_(MTY, k, j, i) = my;
            consTemp_(MTZ, k, j, i) = mz;
            consTemp_(ENG, k, j, i) = eng;
        }
    }
}else{
    #pragma omp parallel for collapse(2) schedule(static)    
    for (int k=kl; k<kr; k++){
        for (int j=jl; j<jr; j++){
    #pragma omp simd
            for (int i=il; i<ir; i++){
                const Real u1 = cons_(DEN, k, j, i);
                const Real u2 = cons_(MTX, k, j, i);
                const Real u3 = cons_(MTY, k, j, i);
                const Real u4 = cons_(MTZ, k, j, i);
                const Real u5 = cons_(ENG, k, j, i);

                const Real dF1 = (flx1(DEN, k, j, i+1) - flx1(DEN, k, j, i)) * dtdx
                               + (flx2(DEN, k, j+1, i) - flx2(DEN, k, j, i)) * dtdy
                               + (flx3(DEN, k+1, j, i) - flx3(DEN, k, j, i)) * dtdz; 
                const Real dF2 = (flx1(MTX, k, j, i+1) - flx1(MTX, k, j, i)) * dtdx
                               + (flx2(MTX, k, j+1, i) - flx2(MTX, k, j, i)) * dtdy
                               + (flx3(MTX, k+1, j, i) - flx3(MTX, k, j, i)) * dtdz;
                const Real dF3 = (flx1(MTY, k, j, i+1) - flx1(MTY, k, j, i)) * dtdx
                               + (flx2(MTY, k, j+1, i) - flx2(MTY, k, j, i)) * dtdy
                               + (flx3(MTY, k+1, j, i) - flx3(MTY, k, j, i)) * dtdz;
                const Real dF4 = (flx1(MTZ, k, j, i+1) - flx1(MTZ, k, j, i)) * dtdx
                               + (flx2(MTZ, k, j+1, i) - flx2(MTZ, k, j, i)) * dtdy
                               + (flx3(MTZ, k+1, j, i) - flx3(MTZ, k, j, i)) * dtdz;
                const Real dF5 = (flx1(ENG, k, j, i+1) - flx1(ENG, k, j, i)) * dtdx
                               + (flx2(ENG, k, j+1, i) - flx2(ENG, k, j, i)) * dtdy
                               + (flx3(ENG, k+1, j, i) - flx3(ENG, k, j, i)) * dtdz;

                Real den = coef1*u1 - coef2*dF1;
                Real mx  = coef1*u2 - coef2*dF2;
                Real my  = coef1*u3 - coef2*dF3;
                Real mz  = coef1*u4 - coef2*dF4;
                Real eng = coef1*u5 - coef2*dF5;

                if (den < dmin) {
                    den = dmin;
                    mx  = 0.0;
                    my  = 0.0;
                    mz  = 0.0;
                    eng = pmin / gm1 + 0.0;   // ke = 0
                } else {
                    Real inv_den = 1.0 / den;
                    Real ke = 0.5 * inv_den * (mx*mx + my*my + mz*mz);
                    Real press = gm1 * (eng - ke);
                    if (press < pmin) {
                        eng = pmin / gm1 + ke;   
                    }
                }
                consTemp_(DEN, k, j, i) = den;
                consTemp_(MTX, k, j, i) = mx;
                consTemp_(MTY, k, j, i) = my;
                consTemp_(MTZ, k, j, i) = mz;
                consTemp_(ENG, k, j, i) = eng;
            }
        }
    }
}
}

void Sim::UpdateCons(const TArray<Real>& cons_, TArray<Real>& consTemp_, const Real coef1, const Real coef2, const Real coef3){
    int il = grid.ib; 
    int ir = grid.ie; 
    int jl = grid.jb; 
    int jr = grid.je; 
    int kl = grid.kb;
    int kr = grid.ke; 

    Real dtdx = dt*domain.dxRec; 
    Real dtdy = dt*domain.dyRec; 
    Real dtdz = dt*domain.dzRec; 

    const Real dmin = DENSITY_FLOOR;
    const Real pmin = PRESSURE_FLOOR;
    const Real gm1 = eos.GetGamma() - 1; 
if (grid.nz == 1){ 
    int k = kl; 
    #pragma omp parallel for schedule(static)    
    for (int j=jl; j<jr; j++){
    #pragma omp simd
        for (int i=il; i<ir; i++){
            const Real u1 = cons_(DEN, k, j, i);
            const Real u2 = cons_(MTX, k, j, i);
            const Real u3 = cons_(MTY, k, j, i);
            const Real u4 = cons_(MTZ, k, j, i);
            const Real u5 = cons_(ENG, k, j, i);
            Real& ut1 = consTemp_(DEN, k, j, i);
            Real& ut2 = consTemp_(MTX, k, j, i);
            Real& ut3 = consTemp_(MTY, k, j, i);
            Real& ut4 = consTemp_(MTZ, k, j, i);
            Real& ut5 = consTemp_(ENG, k, j, i);

            const Real dF1 = (flx1(DEN, k, j, i+1) - flx1(DEN, k, j, i)) * dtdx
                           + (flx2(DEN, k, j+1, i) - flx2(DEN, k, j, i)) * dtdy;
            const Real dF2 = (flx1(MTX, k, j, i+1) - flx1(MTX, k, j, i)) * dtdx
                           + (flx2(MTX, k, j+1, i) - flx2(MTX, k, j, i)) * dtdy;
            const Real dF3 = (flx1(MTY, k, j, i+1) - flx1(MTY, k, j, i)) * dtdx
                           + (flx2(MTY, k, j+1, i) - flx2(MTY, k, j, i)) * dtdy;
            const Real dF4 = (flx1(MTZ, k, j, i+1) - flx1(MTZ, k, j, i)) * dtdx
                           + (flx2(MTZ, k, j+1, i) - flx2(MTZ, k, j, i)) * dtdy;
            const Real dF5 = (flx1(ENG, k, j, i+1) - flx1(ENG, k, j, i)) * dtdx
                           + (flx2(ENG, k, j+1, i) - flx2(ENG, k, j, i)) * dtdy;

            Real den = coef1*u1 + coef2*ut1 - coef3*dF1;
            Real mx  = coef1*u2 + coef2*ut2 - coef3*dF2;
            Real my  = coef1*u3 + coef2*ut3 - coef3*dF3;
            Real mz  = coef1*u4 + coef2*ut4 - coef3*dF4;
            Real eng = coef1*u5 + coef2*ut5 - coef3*dF5;

            if (den < dmin) {
                den = dmin;
                mx  = 0.0;
                my  = 0.0;
                mz  = 0.0;
                eng = pmin / gm1 + 0.0;   // ke = 0
            } else {
                Real inv_den = 1.0 / den;
                Real ke = 0.5 * inv_den * (mx*mx + my*my + mz*mz);
                Real press = gm1 * (eng - ke);
                if (press < pmin) {
                    eng = pmin / gm1 + ke;   
                }
            }

            ut1 = den;
            ut2 = mx;
            ut3 = my;
            ut4 = mz;
            ut5 = eng;
        }
    }
}else{
    #pragma omp parallel for collapse(2) schedule(static)    
    for (int k=kl; k<kr; k++){
        for (int j=jl; j<jr; j++){
    #pragma omp simd
            for (int i=il; i<ir; i++){
                const Real u1 = cons_(DEN, k, j, i);
                const Real u2 = cons_(MTX, k, j, i);
                const Real u3 = cons_(MTY, k, j, i);
                const Real u4 = cons_(MTZ, k, j, i);
                const Real u5 = cons_(ENG, k, j, i);
                Real& ut1 = consTemp_(DEN, k, j, i);
                Real& ut2 = consTemp_(MTX, k, j, i);
                Real& ut3 = consTemp_(MTY, k, j, i);
                Real& ut4 = consTemp_(MTZ, k, j, i);
                Real& ut5 = consTemp_(ENG, k, j, i);

                const Real dF1 = (flx1(DEN, k, j, i+1) - flx1(DEN, k, j, i)) * dtdx
                                + (flx2(DEN, k, j+1, i) - flx2(DEN, k, j, i)) * dtdy
                                + (flx3(DEN, k+1, j, i) - flx3(DEN, k, j, i)) * dtdz; 
                const Real dF2 = (flx1(MTX, k, j, i+1) - flx1(MTX, k, j, i)) * dtdx
                                + (flx2(MTX, k, j+1, i) - flx2(MTX, k, j, i)) * dtdy
                                + (flx3(MTX, k+1, j, i) - flx3(MTX, k, j, i)) * dtdz;
                const Real dF3 = (flx1(MTY, k, j, i+1) - flx1(MTY, k, j, i)) * dtdx
                                + (flx2(MTY, k, j+1, i) - flx2(MTY, k, j, i)) * dtdy
                                + (flx3(MTY, k+1, j, i) - flx3(MTY, k, j, i)) * dtdz;
                const Real dF4 = (flx1(MTZ, k, j, i+1) - flx1(MTZ, k, j, i)) * dtdx
                                + (flx2(MTZ, k, j+1, i) - flx2(MTZ, k, j, i)) * dtdy
                                + (flx3(MTZ, k+1, j, i) - flx3(MTZ, k, j, i)) * dtdz;
                const Real dF5 = (flx1(ENG, k, j, i+1) - flx1(ENG, k, j, i)) * dtdx
                                + (flx2(ENG, k, j+1, i) - flx2(ENG, k, j, i)) * dtdy
                                + (flx3(ENG, k+1, j, i) - flx3(ENG, k, j, i)) * dtdz;
                
                Real den = coef1*u1 + coef2*ut1 - coef3*dF1;
                Real mx  = coef1*u2 + coef2*ut2 - coef3*dF2;
                Real my  = coef1*u3 + coef2*ut3 - coef3*dF3;
                Real mz  = coef1*u4 + coef2*ut4 - coef3*dF4;
                Real eng = coef1*u5 + coef2*ut5 - coef3*dF5;

                if (den < dmin) {
                    den = dmin;
                    mx  = 0.0;
                    my  = 0.0;
                    mz  = 0.0;
                    eng = pmin / gm1 + 0.0;   // ke = 0
                } else {
                    Real inv_den = 1.0 / den;
                    Real ke = 0.5 * inv_den * (mx*mx + my*my + mz*mz);
                    Real press = gm1 * (eng - ke);
                    if (press < pmin) {
                        eng = pmin / gm1 + ke;   
                    }
                }

                ut1 = den;
                ut2 = mx;
                ut3 = my;
                ut4 = mz;
                ut5 = eng;
            }
        }
    }
}
}

} // namespace Gaukuk
