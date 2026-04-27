// Gaukuk dependence
#include "../gaukuk.hpp" 
#include "../sim.hpp"
#include "boundary.hpp"

namespace Gaukuk
{

// *** Simple copy boundary condition ***
//
// copy the edge cell of the activated zone 
// to all the ghost cell next to it
//
//------------------------------------------------------------
// X direction, left side 
void Boundary::OutflowCopyXL(TArray<Real>& cons, const Grid& grid){
    int il = grid.igb;                      // first ghost cell left side
    int ir = grid.ib;                       // first activated cell 
    int jl = grid.jb;                       // first activated cell 
    int jr = grid.je;                       // first ghost cell right side
    int kl = grid.kb;                       // first activated cell 
    int kr = grid.ke;                       // first ghost cell right side
    int iTarget = ir;                       // copy cell's id
#pragma omp parallel for collapse(3) schedule(static)
    for (int ivar=DEN; ivar<=ENG; ivar++){
        // loop z activated zone  
        for (int k=kl; k<kr; k++){
            // loop y activated zone 
            for (int j=jl; j<jr; j++){
#pragma omp simd
                // loop x ghost zone
                for (int i=il; i<ir; i++){
                    Real uTarget = cons(ivar, k, j, iTarget); 
                    Real& u = cons(ivar, k, j, i); 
                    u = uTarget;
                }
            }
        }
    }
}

//------------------------------------------------------------
// X direction, right side 
void Boundary::OutflowCopyXR(TArray<Real>& cons, const Grid& grid){
    int il = grid.ie;                       // first ghost cell right side
    int ir = grid.ige;                      // last ghost cell right side + 1
    int jl = grid.jb;                       // first activated cell 
    int jr = grid.je;                       // first ghost cell right side
    int kl = grid.kb;                       // first activated cell 
    int kr = grid.ke;                       // first ghost cell right side
    int iTarget = il - 1;                   // copy cell's id
#pragma omp parallel for collapse(3) schedule(static)
    for (int ivar=DEN; ivar<=ENG; ivar++){
        // loop z activated zone  
        for (int k=kl; k<kr; k++){
            // loop y activated zone 
            for (int j=jl; j<jr; j++){
#pragma omp simd
                // loop x ghost zone
                for (int i=il; i<ir; i++){
                    Real uTarget = cons(ivar, k, j, iTarget); 
                    Real& u = cons(ivar, k, j, i); 
                    u = uTarget;
                }
            }
        }
    }
}

//------------------------------------------------------------
// Y direction, left side 
void Boundary::OutflowCopyYL(TArray<Real>& cons, const Grid& grid){
    int il = grid.ib;                       // first activated cell 
    int ir = grid.ie;                       // first ghost cell right side
    int jl = grid.jgb;                      // first ghost cell left side 
    int jr = grid.jb;                       // first activated cell 
    int kl = grid.kb;                       // first activated cell 
    int kr = grid.ke;                       // first ghost cell right side
    int jTarget = jr;                       // copy cell's id
#pragma omp parallel for collapse(3) schedule(static)
    for (int ivar=DEN; ivar<=ENG; ivar++){
        // loop z activated zone  
        for (int k=kl; k<kr; k++){
            // loop y ghost zone 
            for (int j=jl; j<jr; j++){
#pragma omp simd
                // loop x activated zone
                for (int i=il; i<ir; i++){
                    Real uTarget = cons(ivar, k, jTarget, i); 
                    Real& u = cons(ivar, k, j, i); 
                    u = uTarget;
                }
            }
        }
    }
}

//------------------------------------------------------------
// Y direction, right side 
void Boundary::OutflowCopyYR(TArray<Real>& cons, const Grid& grid){
    int il = grid.ib;                       // first activated cell 
    int ir = grid.ie;                       // first ghost cell right side
    int jl = grid.je;                       // first ghost cell right side 
    int jr = grid.jge;                      // last ghost cell right side + 1
    int kl = grid.kb;                       // first activated cell 
    int kr = grid.ke;                       // first ghost cell right side
    int jTarget = jl - 1;                   // copy cell's id
#pragma omp parallel for collapse(3) schedule(static)
    for (int ivar=DEN; ivar<=ENG; ivar++){
        // loop z activated zone  
        for (int k=kl; k<kr; k++){
            // loop y ghost zone 
            for (int j=jl; j<jr; j++){
#pragma omp simd
                // loop x activated zone
                for (int i=il; i<ir; i++){
                    Real uTarget = cons(ivar, k, jTarget, i); 
                    Real& u = cons(ivar, k, j, i); 
                    u = uTarget;
                }
            }
        }
    }
}

//------------------------------------------------------------
// Z direction, left side 
void Boundary::OutflowCopyZL(TArray<Real>& cons, const Grid& grid){
    int il = grid.ib;                       // first activated cell 
    int ir = grid.ie;                       // first ghost cell right side
    int jl = grid.jb;                       // first activated cell 
    int jr = grid.je;                       // first ghost cell right side
    int kl = grid.kgb;                      // first ghost cell left side 
    int kr = grid.kb;                       // first activated cell 
    int kTarget = kr;                       // copy cell's id
#pragma omp parallel for collapse(3) schedule(static)
    for (int ivar=DEN; ivar<=ENG; ivar++){
        // loop z ghost zone  
        for (int k=kl; k<kr; k++){
            // loop y activated zone 
            for (int j=jl; j<jr; j++){
#pragma omp simd
                // loop x activated zone
                for (int i=il; i<ir; i++){
                    Real uTarget = cons(ivar, kTarget, j, i); 
                    Real& u = cons(ivar, k, j, i); 
                    u = uTarget;
                }
            }
        }
    }
}

//------------------------------------------------------------
// Z direction, right side 
void Boundary::OutflowCopyZR(TArray<Real>& cons, const Grid& grid){
    int il = grid.ib;                       // first activated cell 
    int ir = grid.ie;                       // first ghost cell right side
    int jl = grid.jb;                       // first activated cell 
    int jr = grid.je;                       // first ghost cell right side
    int kl = grid.ke;                       // first ghost cell right side  
    int kr = grid.kge;                      // last ghost cell right side + 1
    int kTarget = kl - 1;                   // copy cell's id
#pragma omp parallel for collapse(3) schedule(static)
    for (int ivar=DEN; ivar<=ENG; ivar++){
        // loop z ghost zone  
        for (int k=kl; k<kr; k++){
            // loop y activated zone 
            for (int j=jl; j<jr; j++){
#pragma omp simd
                // loop x activated zone
                for (int i=il; i<ir; i++){
                    Real uTarget = cons(ivar, kTarget, j, i); 
                    Real& u = cons(ivar, k, j, i); 
                    u = uTarget;
                }
            }
        }
    }
}

//------------------------------------------------------------
} // namespace Gaukuk