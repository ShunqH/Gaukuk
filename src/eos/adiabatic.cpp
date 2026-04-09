#include <cmath>            // sqrt(), abs()
#include <algorithm>        // max()  

#include "eos.hpp"
#include "../utils/utils.hpp" // Config 

namespace Gaukuk{
    
// adiabatic equation of state
EquationOfState::EquationOfState(){
    // read adiabatic index gamma from input file
    gamma_ = Config::getInstance().get("gamma") ; 
    // density minimum (floor)
    densityMin_ = Config::getInstance().get("rho_floor", 1e-16) ; 
    // pressure minimum (floor)
    pressureMin_ = Config::getInstance().get("pressure_floor", 1e-16) ; 
    // 1/(gamma-1)
    gm1Rec_ = 1.0 / (gamma_ - 1.0);
}

void EquationOfState::ConsToPrim(TArray<Real>& cons, TArray<Real>& prim, const Grid& grid){
    int il = grid.igb;                      // first ghost cell left side
    int ir = grid.ige;                      // last ghost cell right side + 1
    int jl = grid.jgb;                      // first ghost cell left side
    int jr = grid.jge;                      // last ghost cell right side + 1
    int kl = grid.kgb;                      // first ghost cell left side
    int kr = grid.kge;                      // last ghost cell right side + 1
    Real gm1 = gamma_ - 1.0;
#pragma omp parallel for
    for (int k=kl; k<kr; k++){
        for (int j=jl; j<jr; j++){
#pragma omp simd 
            for (int i=il; i<ir; i++){
                Real& consDen = cons(DEN, k, j, i); 
                Real& consMtx = cons(MTX, k, j, i); 
                Real& consMty = cons(MTY, k, j, i); 
                Real& consMtz = cons(MTZ, k, j, i); 
                Real& consEng = cons(ENG, k, j, i); 

                Real& primDen = prim(DEN, k, j, i); 
                Real& primVlx = prim(VLX, k, j, i); 
                Real& primVly = prim(VLY, k, j, i); 
                Real& primVlz = prim(VLZ, k, j, i); 
                Real& primPre = prim(PRE, k, j, i); 

                consDen = (consDen > densityMin_) ? consDen : densityMin_; 
                primDen = consDen; 
                Real denInv = 1.0/consDen; 
                primVlx = consMtx * denInv; 
                primVly = consMty * denInv; 
                primVlz = consMtz * denInv; 
                Real engKin = 0.5 * denInv * (consMtx*consMtx + consMty*consMty + consMtz*consMtz); 
                primPre = gm1 * (consEng - engKin); 
                primPre = (primPre > pressureMin_) ? primPre : pressureMin_; 
                consEng = (primPre > pressureMin_) ? consEng : ( pressureMin_/gm1 + engKin ); 
            }
        }
    }
}

void EquationOfState::ConsToPrim(TArray<Real>& cons, TArray<Real>& prim, const Grid& grid, Real& cmax){
    int il = grid.igb;                      // first ghost cell left side
    int ir = grid.ige;                      // last ghost cell right side + 1
    int jl = grid.jgb;                      // first ghost cell left side
    int jr = grid.jge;                      // last ghost cell right side + 1
    int kl = grid.kgb;                      // first ghost cell left side
    int kr = grid.kge;                      // last ghost cell right side + 1
    Real gm1 = gamma_ - 1.0;
#pragma omp parallel for collapse(2) reduction(max: cmax)
    for (int k=kl; k<kr; k++){
        for (int j=jl; j<jr; j++){
            Real cmaxRow = 1e-16;
#pragma omp simd reduction(max: cmaxRow)
            for (int i=il; i<ir; i++){
                Real& consDen = cons(DEN, k, j, i); 
                Real& consMtx = cons(MTX, k, j, i); 
                Real& consMty = cons(MTY, k, j, i); 
                Real& consMtz = cons(MTZ, k, j, i); 
                Real& consEng = cons(ENG, k, j, i); 

                Real& primDen = prim(DEN, k, j, i); 
                Real& primVlx = prim(VLX, k, j, i); 
                Real& primVly = prim(VLY, k, j, i); 
                Real& primVlz = prim(VLZ, k, j, i); 
                Real& primPre = prim(PRE, k, j, i); 

                consDen = (consDen > densityMin_) ? consDen : densityMin_; 
                primDen = consDen; 
                Real denInv = 1.0/consDen; 
                primVlx = consMtx * denInv; 
                primVly = consMty * denInv; 
                primVlz = consMtz * denInv; 
                Real engKin = 0.5 * denInv * (consMtx*consMtx + consMty*consMty + consMtz*consMtz); 
                primPre = gm1 * (consEng - engKin); 
                primPre = (primPre > pressureMin_) ? primPre : pressureMin_; 
                consEng = (primPre > pressureMin_) ? consEng : ( pressureMin_/gm1 + engKin ); 

                // also calculate the |v|+cs maximum
                Real cs = std::sqrt(gamma_ * primPre * denInv);
                Real local = std::max({
                    std::abs(primVlx) + cs,
                    std::abs(primVly) + cs,
                    std::abs(primVlz) + cs
                });
                cmaxRow = std::max(cmaxRow, local);
            }
            cmax = cmaxRow;
        }
    }
}

void EquationOfState::PrimToCons(const TArray<Real>& prim, TArray<Real>& cons, const Grid& grid){
    int il = grid.igb;                      // first ghost cell left side
    int ir = grid.ige;                      // last ghost cell right side + 1
    int jl = grid.jgb;                      // first ghost cell left side
    int jr = grid.jge;                      // last ghost cell right side + 1
    int kl = grid.kgb;                      // first ghost cell left side
    int kr = grid.kge;                      // last ghost cell right side + 1
#pragma omp parallel for
    for (int k=kl; k<kr; k++){
        for (int j=jl; j<jr; j++){
#pragma omp simd 
            for (int i=il; i<ir; i++){
                Real& consDen = cons(DEN, k, j, i); 
                Real& consMtx = cons(MTX, k, j, i); 
                Real& consMty = cons(MTY, k, j, i); 
                Real& consMtz = cons(MTZ, k, j, i); 
                Real& consEng = cons(ENG, k, j, i); 

                const Real& primDen = prim(DEN, k, j, i); 
                const Real& primVlx = prim(VLX, k, j, i); 
                const Real& primVly = prim(VLY, k, j, i); 
                const Real& primVlz = prim(VLZ, k, j, i); 
                const Real& primPre = prim(PRE, k, j, i); 

                consDen = primDen; 
                consMtx = primDen * primVlx; 
                consMty = primDen * primVly; 
                consMtz = primDen * primVlz; 
                consEng = primPre * gm1Rec_ + 0.5 * primDen * ( primVlx*primVlx + primVly*primVly + primVlz*primVlz );
            }
        }
    }
}

} // namespace Gaukuk
