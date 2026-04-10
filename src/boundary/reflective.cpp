// C++ Headers
// #include <iostream>     // std::cout; std::endl; std::cerr

// Gaukuk dependence
#include "../gaukuk.hpp" 
#include "../sim.hpp"
#include "boundary.hpp"

namespace Gaukuk
{

// *** reflective boundary condition ***
//
// mirror the activated cell to the ghost cells 
// along the boundary
//
//------------------------------------------------------------
// X direction, left side 
void Boundary::ReflectiveXL(TArray<Real>& cons, const Grid& grid){
    int il = grid.igb;                      // first ghost cell left side
    int ir = grid.ib;                       // first activated cell 
    int jl = grid.jb;                       // first activated cell 
    int jr = grid.je;                       // first ghost cell right side
    int kl = grid.kb;                       // first activated cell 
    int kr = grid.ke;                       // first ghost cell right side
#pragma omp parallel for collapse(2) schedule(static)
    for (int k=kl; k<kr; k++){
        // loop y activated zone 
        for (int j=jl; j<jr; j++){
#pragma omp simd
            // loop x ghost zone
            for (int i=il; i<ir; i++){
                int iTarget = ir + ir - i - 1; 
                const Real& rhoTarget = cons(DEN, k, j, iTarget); 
                const Real& mtxTarget = cons(MTX, k, j, iTarget); 
                const Real& mtyTarget = cons(MTY, k, j, iTarget); 
                const Real& mtzTarget = cons(MTZ, k, j, iTarget); 
                const Real& engTarget = cons(ENG, k, j, iTarget); 
                Real& rho = cons(DEN, k, j, i); 
                Real& mtx = cons(MTX, k, j, i); 
                Real& mty = cons(MTY, k, j, i); 
                Real& mtz = cons(MTZ, k, j, i); 
                Real& eng = cons(ENG, k, j, i); 
                rho = rhoTarget;
                mtx = - mtxTarget;          // reflect momentum in x 
                mty = mtyTarget;
                mtz = mtzTarget;
                eng = engTarget;
            }
        }
    }
}

//------------------------------------------------------------
// X direction, right side 
void Boundary::ReflectiveXR(TArray<Real>& cons, const Grid& grid){
    int il = grid.ie;                       // first ghost cell right side
    int ir = grid.ige;                      // last ghost cell right side + 1
    int jl = grid.jb;                       // first activated cell 
    int jr = grid.je;                       // first ghost cell right side
    int kl = grid.kb;                       // first activated cell 
    int kr = grid.ke;                       // first ghost cell right side
#pragma omp parallel for collapse(2) schedule(static)
    // loop z activated zone  
    for (int k=kl; k<kr; k++){
        // loop y activated zone 
        for (int j=jl; j<jr; j++){
#pragma omp simd
            // loop x ghost zone
            for (int i=il; i<ir; i++){
                int iTarget = il + il - i - 1; 
                const Real& rhoTarget = cons(DEN, k, j, iTarget); 
                const Real& mtxTarget = cons(MTX, k, j, iTarget); 
                const Real& mtyTarget = cons(MTY, k, j, iTarget); 
                const Real& mtzTarget = cons(MTZ, k, j, iTarget); 
                const Real& engTarget = cons(ENG, k, j, iTarget); 
                Real& rho = cons(DEN, k, j, i); 
                Real& mtx = cons(MTX, k, j, i); 
                Real& mty = cons(MTY, k, j, i); 
                Real& mtz = cons(MTZ, k, j, i); 
                Real& eng = cons(ENG, k, j, i); 
                rho = rhoTarget;
                mtx = - mtxTarget;          // reflect momentum in x 
                mty = mtyTarget;
                mtz = mtzTarget;
                eng = engTarget;
            }
        }
    }
}

//------------------------------------------------------------
// Y direction, left side 
void Boundary::ReflectiveYL(TArray<Real>& cons, const Grid& grid){
    int il = grid.ib;                       // first activated cell 
    int ir = grid.ie;                       // first ghost cell right side
    int jl = grid.jgb;                      // first ghost cell left side 
    int jr = grid.jb;                       // first activated cell 
    int kl = grid.kb;                       // first activated cell 
    int kr = grid.ke;                       // first ghost cell right side
#pragma omp parallel for collapse(2) schedule(static)
    // loop z activated zone  
    for (int k=kl; k<kr; k++){
        // loop y ghost zone 
        for (int j=jl; j<jr; j++){
#pragma omp simd
            // loop x activated zone
            for (int i=il; i<ir; i++){
                int jTarget = jr + jr - j - 1; 
                const Real& rhoTarget = cons(DEN, k, jTarget, i); 
                const Real& mtxTarget = cons(MTX, k, jTarget, i); 
                const Real& mtyTarget = cons(MTY, k, jTarget, i); 
                const Real& mtzTarget = cons(MTZ, k, jTarget, i); 
                const Real& engTarget = cons(ENG, k, jTarget, i); 
                Real& rho = cons(DEN, k, j, i); 
                Real& mtx = cons(MTX, k, j, i); 
                Real& mty = cons(MTY, k, j, i); 
                Real& mtz = cons(MTZ, k, j, i); 
                Real& eng = cons(ENG, k, j, i); 
                rho = rhoTarget;
                mtx = mtxTarget;           
                mty = - mtyTarget;          // reflect momentum in y
                mtz = mtzTarget;
                eng = engTarget;
            }
        }
    }
}

//------------------------------------------------------------
// Y direction, right side 
void Boundary::ReflectiveYR(TArray<Real>& cons, const Grid& grid){
    int il = grid.ib;                       // first activated cell 
    int ir = grid.ie;                       // first ghost cell right side
    int jl = grid.je;                       // first ghost cell right side 
    int jr = grid.jge;                      // last ghost cell right side + 1
    int kl = grid.kb;                       // first activated cell 
    int kr = grid.ke;                       // first ghost cell right side
#pragma omp parallel for collapse(2) schedule(static)
    // loop z activated zone  
    for (int k=kl; k<kr; k++){
        // loop y ghost zone 
        for (int j=jl; j<jr; j++){
#pragma omp simd
            // loop x activated zone
            for (int i=il; i<ir; i++){
                int jTarget = jl + jl - j - 1; 
                const Real& rhoTarget = cons(DEN, k, jTarget, i); 
                const Real& mtxTarget = cons(MTX, k, jTarget, i); 
                const Real& mtyTarget = cons(MTY, k, jTarget, i); 
                const Real& mtzTarget = cons(MTZ, k, jTarget, i); 
                const Real& engTarget = cons(ENG, k, jTarget, i); 
                Real& rho = cons(DEN, k, j, i); 
                Real& mtx = cons(MTX, k, j, i); 
                Real& mty = cons(MTY, k, j, i); 
                Real& mtz = cons(MTZ, k, j, i); 
                Real& eng = cons(ENG, k, j, i); 
                rho = rhoTarget;
                mtx = mtxTarget;           
                mty = - mtyTarget;          // reflect momentum in y
                mtz = mtzTarget;
                eng = engTarget;
            }
        }
    }
}

//------------------------------------------------------------
// Z direction, left side 
void Boundary::ReflectiveZL(TArray<Real>& cons, const Grid& grid){
    int il = grid.ib;                       // first activated cell 
    int ir = grid.ie;                       // first ghost cell right side
    int jl = grid.jb;                       // first activated cell 
    int jr = grid.je;                       // first ghost cell right side
    int kl = grid.kgb;                      // first ghost cell left side 
    int kr = grid.kb;                       // first activated cell 
#pragma omp parallel for collapse(2) schedule(static)
    // loop z ghost zone  
    for (int k=kl; k<kr; k++){
        // loop y activated zone 
        for (int j=jl; j<jr; j++){
#pragma omp simd
            // loop x activated zone
            for (int i=il; i<ir; i++){
                int kTarget = kr + kr - k - 1; 
                const Real& rhoTarget = cons(DEN, kTarget, j, i); 
                const Real& mtxTarget = cons(MTX, kTarget, j, i); 
                const Real& mtyTarget = cons(MTY, kTarget, j, i); 
                const Real& mtzTarget = cons(MTZ, kTarget, j, i); 
                const Real& engTarget = cons(ENG, kTarget, j, i); 
                Real& rho = cons(DEN, k, j, i); 
                Real& mtx = cons(MTX, k, j, i); 
                Real& mty = cons(MTY, k, j, i); 
                Real& mtz = cons(MTZ, k, j, i); 
                Real& eng = cons(ENG, k, j, i); 
                rho = rhoTarget;
                mtx = mtxTarget;           
                mty = mtyTarget;          
                mtz = - mtzTarget;          // reflect momentum in z
                eng = engTarget;
            }
        }
    }
}

//------------------------------------------------------------
// Z direction, right side 
void Boundary::ReflectiveZR(TArray<Real>& cons, const Grid& grid){
    int il = grid.ib;                       // first activated cell 
    int ir = grid.ie;                       // first ghost cell right side
    int jl = grid.jb;                       // first activated cell 
    int jr = grid.je;                       // first ghost cell right side
    int kl = grid.ke;                       // first ghost cell right side  
    int kr = grid.kge;                      // last ghost cell right side + 1
#pragma omp parallel for collapse(2) schedule(static)
    // loop z ghost zone  
    for (int k=kl; k<kr; k++){
        // loop y activated zone 
        for (int j=jl; j<jr; j++){
#pragma omp simd
            // loop x activated zone
            for (int i=il; i<ir; i++){
                int kTarget = kl + kl - k - 1; 
                const Real& rhoTarget = cons(DEN, kTarget, j, i); 
                const Real& mtxTarget = cons(MTX, kTarget, j, i); 
                const Real& mtyTarget = cons(MTY, kTarget, j, i); 
                const Real& mtzTarget = cons(MTZ, kTarget, j, i); 
                const Real& engTarget = cons(ENG, kTarget, j, i); 
                Real& rho = cons(DEN, k, j, i); 
                Real& mtx = cons(MTX, k, j, i); 
                Real& mty = cons(MTY, k, j, i); 
                Real& mtz = cons(MTZ, k, j, i); 
                Real& eng = cons(ENG, k, j, i); 
                rho = rhoTarget;
                mtx = mtxTarget;           
                mty = mtyTarget;          
                mtz = - mtzTarget;          // reflect momentum in z
                eng = engTarget;
            }
        }
    }
}

//------------------------------------------------------------
} // namespace Gaukuk