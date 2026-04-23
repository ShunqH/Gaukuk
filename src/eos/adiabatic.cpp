// C++ headers
#include <cmath>            // sqrt(), abs()
#include <algorithm>        // max()  

// Gaukuk dependence
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
    if (grid.nz==1){
        kl = grid.kb;
        kr = grid.ke;
    }
    Real gm1 = gamma_ - 1.0;
    Real dmin = densityMin_;
    Real pmin = pressureMin_;
#pragma omp parallel for collapse(2) schedule(static)
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

                consDen = (consDen > dmin) ? consDen : dmin; 
                primDen = consDen; 
                Real denInv = 1.0/consDen; 
                primVlx = consMtx * denInv; 
                primVly = consMty * denInv; 
                primVlz = consMtz * denInv; 
                Real engKin = 0.5 * denInv * (consMtx*consMtx + consMty*consMty + consMtz*consMtz); 
                primPre = gm1 * (consEng - engKin); 
                primPre = (primPre > pmin) ? primPre : pmin; 
                consEng = (primPre > pmin) ? consEng : ( pmin/gm1 + engKin ); 
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
    if (grid.nz==1){
        kl = grid.kb;
        kr = grid.ke;
    }
    Real gm1 = gamma_ - 1.0;
    Real dmin = densityMin_;
    Real pmin = pressureMin_;
#pragma omp parallel for collapse(2) reduction(max: cmax) schedule(static)
    for (int k=kl; k<kr; k++){
        for (int j=jl; j<jr; j++){
#pragma omp simd reduction(max: cmax)
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

                consDen = (consDen > dmin) ? consDen : dmin; 
                primDen = consDen; 
                Real denInv = 1.0/consDen; 
                primVlx = consMtx * denInv; 
                primVly = consMty * denInv; 
                primVlz = consMtz * denInv; 
                Real engKin = 0.5 * denInv * (consMtx*consMtx + consMty*consMty + consMtz*consMtz); 
                primPre = gm1 * (consEng - engKin); 
                primPre = (primPre > pmin) ? primPre : pmin; 
                consEng = (primPre > pmin) ? consEng : ( pmin/gm1 + engKin ); 

                // also calculate the |v|+cs maximum
                Real cs = std::sqrt(gamma_ * primPre * denInv);
                Real local = std::max({
                    std::abs(primVlx) + cs,
                    std::abs(primVly) + cs,
                    std::abs(primVlz) + cs
                });
                cmax = std::max(cmax, local);
            }
        }
    }
}

void EquationOfState::CalCmax(const TArray<Real>& cons, const Grid& grid, Real& cmax){
    int il = grid.ib;               // first activated cell left side
    int ir = grid.ie;               // last activated cell right side + 1
    int jl = grid.jb;               // first activated cell left side
    int jr = grid.je;               // last activated cell right side + 1
    int kl = grid.kb;               // first activated cell left side
    int kr = grid.ke;               // last activated cell right side + 1
    if (grid.nz==1){
        kl = grid.kb;
        kr = grid.ke;
    }
    Real gm1 = gamma_ - 1.0;
    Real dmin = densityMin_;
    Real pmin = pressureMin_;

#pragma omp parallel for collapse(2) reduction(max: cmax) schedule(static)
    for (int k=kl; k<kr; k++){
        for (int j=jl; j<jr; j++){
#pragma omp simd reduction(max: cmax)
            for (int i=il; i<ir; i++){
                Real consDen = cons(DEN, k, j, i); 
                Real consMtx = cons(MTX, k, j, i); 
                Real consMty = cons(MTY, k, j, i); 
                Real consMtz = cons(MTZ, k, j, i); 
                Real consEng = cons(ENG, k, j, i); 

                // Clamp density to local variable
                consDen = (consDen > dmin) ? consDen : dmin;
                Real denInv = 1.0 / consDen;
                Real vx = consMtx * denInv;
                Real vy = consMty * denInv;
                Real vz = consMtz * denInv;
                // Kinetic energy
                Real engKin = 0.5 * (vx * consMtx + vy * consMty + vz * consMtz);
                Real pres = gm1 * (consEng - engKin);
                pres = (pres > pmin) ? pres : pmin;
                Real cs = std::sqrt(gamma_ * pres * denInv);
                Real local = std::abs(vx) + cs;
                local = std::max(local, std::abs(vy) + cs);
                local = std::max(local, std::abs(vz) + cs);
                
                cmax = std::max(cmax, local);
            }
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
    if (grid.nz==1){
        kl = grid.kb;
        kr = grid.ke;
    }
#pragma omp parallel for collapse(2) schedule(static)
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
