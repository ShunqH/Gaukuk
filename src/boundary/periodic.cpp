// C++ Headers
// #include <iostream>     // std::cout; std::endl; std::cerr

// Gaukuk dependence
#include "../gaukuk.hpp" 
#include "../sim.hpp"
#include "boundary.hpp"

namespace Gaukuk
{

//   *** periodic boundary condition ***
//
// copy the edge cell of the activated zone 
// to the other side of the ghost cells
//
//------------------------------------------------------------
// X direction, left side 
void Boundary::PeriodicXL(TArray<Real>& cons, const Grid& grid){
    int il = grid.igb;                      // first ghost cell left side
    int ir = grid.ib;                       // first activated cell 
    int jl = grid.jb;                       // first activated cell 
    int jr = grid.je;                       // first ghost cell right side
    int kl = grid.kb;                       // first activated cell 
    int kr = grid.ke;                       // first ghost cell right side
    int iDist = grid.nx;                    // distance between the activated and ghost cell 
#pragma omp parallel for collapse(3) schedule(static)
    for (int ivar=DEN; ivar<=ENG; ivar++){
        // loop z activated zone  
        for (int k=kl; k<kr; k++){
            // loop y activated zone 
            for (int j=jl; j<jr; j++){
#pragma omp simd
                // loop x ghost zone
                for (int i=il; i<ir; i++){
                    int iTarget = i + iDist; 
                    const Real& uTarget = cons(ivar, k, j, iTarget); 
                    Real& u = cons(ivar, k, j, i); 
                    u = uTarget;
                }
            }
        }
    }
}

//------------------------------------------------------------
// X direction, right side 
void Boundary::PeriodicXR(TArray<Real>& cons, const Grid& grid){
    int il = grid.ie;                       // first ghost cell right side
    int ir = grid.ige;                      // last ghost cell right side + 1
    int jl = grid.jb;                       // first activated cell 
    int jr = grid.je;                       // first ghost cell right side
    int kl = grid.kb;                       // first activated cell 
    int kr = grid.ke;                       // first ghost cell right side
    int iDist = grid.nx;                    // distance between the activated and ghost cell 
#pragma omp parallel for collapse(3) schedule(static)
    for (int ivar=DEN; ivar<=ENG; ivar++){
        // loop z activated zone  
        for (int k=kl; k<kr; k++){
            // loop y activated zone 
            for (int j=jl; j<jr; j++){
#pragma omp simd
                // loop x ghost zone
                for (int i=il; i<ir; i++){
                    int iTarget = i - iDist; 
                    const Real& uTarget = cons(ivar, k, j, iTarget); 
                    Real& u = cons(ivar, k, j, i); 
                    u = uTarget;
                }
            }
        }
    }
}

//------------------------------------------------------------
// Y direction, left side 
void Boundary::PeriodicYL(TArray<Real>& cons, const Grid& grid){
    int il = grid.ib;                       // first activated cell 
    int ir = grid.ie;                       // first ghost cell right side
    int jl = grid.jgb;                      // first ghost cell left side 
    int jr = grid.jb;                       // first activated cell 
    int kl = grid.kb;                       // first activated cell 
    int kr = grid.ke;                       // first ghost cell right side
    int jDist = grid.ny;                    // distance between the activated and ghost cell 
#pragma omp parallel for collapse(3) schedule(static)
    for (int ivar=DEN; ivar<=ENG; ivar++){
        // loop z activated zone  
        for (int k=kl; k<kr; k++){
            // loop y ghost zone 
            for (int j=jl; j<jr; j++){
#pragma omp simd
                // loop x activated zone
                for (int i=il; i<ir; i++){
                    int jTarget = j + jDist; 
                    const Real& uTarget = cons(ivar, k, jTarget, i); 
                    Real& u = cons(ivar, k, j, i); 
                    u = uTarget;
                }
            }
        }
    }
}

//------------------------------------------------------------
// Y direction, right side 
void Boundary::PeriodicYR(TArray<Real>& cons, const Grid& grid){
    int il = grid.ib;                       // first activated cell 
    int ir = grid.ie;                       // first ghost cell right side
    int jl = grid.je;                       // first ghost cell right side 
    int jr = grid.jge;                      // last ghost cell right side + 1
    int kl = grid.kb;                       // first activated cell 
    int kr = grid.ke;                       // first ghost cell right side
    int jDist = grid.ny;                    // distance between the activated and ghost cell 
#pragma omp parallel for collapse(3) schedule(static)
    for (int ivar=DEN; ivar<=ENG; ivar++){
        // loop z activated zone  
        for (int k=kl; k<kr; k++){
            // loop y ghost zone 
            for (int j=jl; j<jr; j++){
#pragma omp simd
                // loop x activated zone
                for (int i=il; i<ir; i++){
                    int jTarget = j - jDist; 
                    const Real& uTarget = cons(ivar, k, jTarget, i); 
                    Real& u = cons(ivar, k, j, i); 
                    u = uTarget;
                }
            }
        }
    }
}

//------------------------------------------------------------
// Z direction, left side 
void Boundary::PeriodicZL(TArray<Real>& cons, const Grid& grid){
    int il = grid.ib;                       // first activated cell 
    int ir = grid.ie;                       // first ghost cell right side
    int jl = grid.jb;                       // first activated cell 
    int jr = grid.je;                       // first ghost cell right side
    int kl = grid.kgb;                      // first ghost cell left side 
    int kr = grid.kb;                       // first activated cell 
    int kDist = grid.nz;                    // distance between the activated and ghost cell 
#pragma omp parallel for collapse(3) schedule(static)
    for (int ivar=DEN; ivar<=ENG; ivar++){
        // loop z activated zone  
        for (int k=kl; k<kr; k++){
            // loop y ghost zone 
            for (int j=jl; j<jr; j++){
#pragma omp simd
                // loop x activated zone
                for (int i=il; i<ir; i++){
                    int kTarget = k + kDist; 
                    const Real& uTarget = cons(ivar, kTarget, j, i); 
                    Real& u = cons(ivar, k, j, i); 
                    u = uTarget;
                }
            }
        }
    }
}

//------------------------------------------------------------
// Z direction, right side 
void Boundary::PeriodicZR(TArray<Real>& cons, const Grid& grid){
    int il = grid.ib;                       // first activated cell 
    int ir = grid.ie;                       // first ghost cell right side
    int jl = grid.jb;                       // first activated cell 
    int jr = grid.je;                       // first ghost cell right side
    int kl = grid.ke;                       // first ghost cell right side  
    int kr = grid.kge;                      // last ghost cell right side + 1
    int kDist = grid.nz;                    // distance between the activated and ghost cell 
#pragma omp parallel for collapse(3) schedule(static)
    for (int ivar=DEN; ivar<=ENG; ivar++){
        // loop z activated zone  
        for (int k=kl; k<kr; k++){
            // loop y ghost zone 
            for (int j=jl; j<jr; j++){
#pragma omp simd
                // loop x activated zone
                for (int i=il; i<ir; i++){
                    int kTarget = k - kDist; 
                    const Real& uTarget = cons(ivar, kTarget, j, i); 
                    Real& u = cons(ivar, k, j, i); 
                    u = uTarget;
                }
            }
        }
    }
}

//------------------------------------------------------------
} // namespace Gaukuk