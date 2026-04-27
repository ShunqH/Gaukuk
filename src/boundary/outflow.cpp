// C++ headers
#include <algorithm>                // std::min(), std::max()

// Gaukuk dependence
#include "../gaukuk.hpp" 
#include "../sim.hpp"
#include "boundary.hpp"

namespace Gaukuk
{

// ***** Outflow boundary condition *****
//
// copy the edge cell of the activated zone 
// to all the ghost cell next to it. 
// if v_n is inflow, set to 0.  
//
//------------------------------------------------------------
// X direction, left side 
void Boundary::OutflowXL(TArray<Real>& cons, const Grid& grid){
    int il = grid.igb;                      // first ghost cell left side
    int ir = grid.ib;                       // first activated cell 
    int jl = grid.jb;                       // first activated cell 
    int jr = grid.je;                       // first ghost cell right side
    int kl = grid.kb;                       // first activated cell 
    int kr = grid.ke;                       // first ghost cell right side
    const Real denMin = DENSITY_FLOOR; 
#pragma omp parallel for collapse(2) schedule(static)
    // loop z activated zone  
    for (int k=kl; k<kr; k++){
        // loop y activated zone 
        for (int j=jl; j<jr; j++){
#pragma omp simd
            // loop x ghost zone
            for (int i=il; i<ir; i++){
                int iTarget = ir + ir - i - 1; 
                Real rhoTarget = cons(DEN, k, j, iTarget); 
                Real mtxTarget = cons(MTX, k, j, iTarget); 
                Real mtyTarget = cons(MTY, k, j, iTarget); 
                Real mtzTarget = cons(MTZ, k, j, iTarget); 
                Real engTarget = cons(ENG, k, j, iTarget); 
                Real& rho = cons(DEN, k, j, i); 
                Real& mtx = cons(MTX, k, j, i); 
                Real& mty = cons(MTY, k, j, i); 
                Real& mtz = cons(MTZ, k, j, i); 
                Real& eng = cons(ENG, k, j, i); 

                Real mtxClamp = std::min(mtxTarget, 0.0); 
                rhoTarget = std::max(rhoTarget, denMin); 
                Real rhoInv = 1.0/rhoTarget; 
                Real keOld = 0.5*rhoInv*mtxTarget*mtxTarget; 
                Real keNew = 0.5*rhoInv*mtxClamp*mtxClamp; 
                Real eint = engTarget - keOld;
                eint = std::max(eint, 1e-16);
                Real engNew = eint + keNew;

                rho = rhoTarget; 
                mtx = mtxClamp; 
                mty = mtyTarget; 
                mtz = mtzTarget; 
                eng = engNew; 
            }
        }
    }
}

//------------------------------------------------------------
// X direction, right side 
void Boundary::OutflowXR(TArray<Real>& cons, const Grid& grid){
    int il = grid.ie;                       // first ghost cell right side
    int ir = grid.ige;                      // last ghost cell right side + 1
    int jl = grid.jb;                       // first activated cell 
    int jr = grid.je;                       // first ghost cell right side
    int kl = grid.kb;                       // first activated cell 
    int kr = grid.ke;                       // first ghost cell right side
    const Real denMin = DENSITY_FLOOR; 
#pragma omp parallel for collapse(2) schedule(static)
    // loop z activated zone  
    for (int k=kl; k<kr; k++){
        // loop y activated zone 
        for (int j=jl; j<jr; j++){
#pragma omp simd
            // loop x ghost zone
            for (int i=il; i<ir; i++){
                int iTarget = il + il - i - 1; 
                Real rhoTarget = cons(DEN, k, j, iTarget); 
                Real mtxTarget = cons(MTX, k, j, iTarget); 
                Real mtyTarget = cons(MTY, k, j, iTarget); 
                Real mtzTarget = cons(MTZ, k, j, iTarget); 
                Real engTarget = cons(ENG, k, j, iTarget); 
                Real& rho = cons(DEN, k, j, i); 
                Real& mtx = cons(MTX, k, j, i); 
                Real& mty = cons(MTY, k, j, i); 
                Real& mtz = cons(MTZ, k, j, i); 
                Real& eng = cons(ENG, k, j, i); 

                Real mtxClamp = std::max(mtxTarget, 0.0); 
                rhoTarget = std::max(rhoTarget, denMin); 
                Real rhoInv = 1.0/rhoTarget; 
                Real keOld = 0.5*rhoInv*mtxTarget*mtxTarget; 
                Real keNew = 0.5*rhoInv*mtxClamp*mtxClamp; 
                Real eint = engTarget - keOld;
                eint = std::max(eint, 1e-16);
                Real engNew = eint + keNew;

                rho = rhoTarget;
                mtx = mtxClamp;     
                mty = mtyTarget;
                mtz = mtzTarget;
                eng = engNew;
            }
        }
    }
}

//------------------------------------------------------------
// Y direction, left side 
void Boundary::OutflowYL(TArray<Real>& cons, const Grid& grid){
    int il = grid.ib;                       // first activated cell 
    int ir = grid.ie;                       // first ghost cell right side
    int jl = grid.jgb;                      // first ghost cell left side 
    int jr = grid.jb;                       // first activated cell 
    int kl = grid.kb;                       // first activated cell 
    int kr = grid.ke;                       // first ghost cell right side
    const Real denMin = DENSITY_FLOOR; 
#pragma omp parallel for collapse(2) schedule(static)
    // loop z activated zone  
    for (int k=kl; k<kr; k++){
        // loop y ghost zone 
        for (int j=jl; j<jr; j++){
#pragma omp simd
            // loop x activated zone
            for (int i=il; i<ir; i++){
                int jTarget = jr + jr - j - 1; 
                Real rhoTarget = cons(DEN, k, jTarget, i); 
                Real mtxTarget = cons(MTX, k, jTarget, i); 
                Real mtyTarget = cons(MTY, k, jTarget, i); 
                Real mtzTarget = cons(MTZ, k, jTarget, i); 
                Real engTarget = cons(ENG, k, jTarget, i); 
                Real& rho = cons(DEN, k, j, i); 
                Real& mtx = cons(MTX, k, j, i); 
                Real& mty = cons(MTY, k, j, i); 
                Real& mtz = cons(MTZ, k, j, i); 
                Real& eng = cons(ENG, k, j, i); 

                Real mtyClamp = std::min(mtyTarget, 0.0); 
                rhoTarget = std::max(rhoTarget, denMin); 
                Real rhoInv = 1.0/rhoTarget; 
                Real keOld = 0.5*rhoInv*mtyTarget*mtyTarget; 
                Real keNew = 0.5*rhoInv*mtyClamp*mtyClamp; 
                Real eint = engTarget - keOld;
                eint = std::max(eint, 1e-16);
                Real engNew = eint + keNew;

                rho = rhoTarget;
                mtx = mtxTarget;           
                mty = mtyClamp;     
                mtz = mtzTarget;
                eng = engNew;
            }
        }
    }
}

//------------------------------------------------------------
// Y direction, right side 
void Boundary::OutflowYR(TArray<Real>& cons, const Grid& grid){
    int il = grid.ib;                       // first activated cell 
    int ir = grid.ie;                       // first ghost cell right side
    int jl = grid.je;                       // first ghost cell right side 
    int jr = grid.jge;                      // last ghost cell right side + 1
    int kl = grid.kb;                       // first activated cell 
    int kr = grid.ke;                       // first ghost cell right side
    const Real denMin = DENSITY_FLOOR; 
#pragma omp parallel for collapse(2) schedule(static)
    // loop z activated zone  
    for (int k=kl; k<kr; k++){
        // loop y ghost zone 
        for (int j=jl; j<jr; j++){
#pragma omp simd
            // loop x activated zone
            for (int i=il; i<ir; i++){
                int jTarget = jl + jl - j - 1; 
                Real rhoTarget = cons(DEN, k, jTarget, i); 
                Real mtxTarget = cons(MTX, k, jTarget, i); 
                Real mtyTarget = cons(MTY, k, jTarget, i); 
                Real mtzTarget = cons(MTZ, k, jTarget, i); 
                Real engTarget = cons(ENG, k, jTarget, i); 
                Real& rho = cons(DEN, k, j, i); 
                Real& mtx = cons(MTX, k, j, i); 
                Real& mty = cons(MTY, k, j, i); 
                Real& mtz = cons(MTZ, k, j, i); 
                Real& eng = cons(ENG, k, j, i); 

                Real mtyClamp = std::max(mtyTarget, 0.0); 
                rhoTarget = std::max(rhoTarget, denMin); 
                Real rhoInv = 1.0/rhoTarget; 
                Real keOld = 0.5*rhoInv*mtyTarget*mtyTarget; 
                Real keNew = 0.5*rhoInv*mtyClamp*mtyClamp; 
                Real eint = engTarget - keOld;
                eint = std::max(eint, 1e-16);
                Real engNew = eint + keNew;

                rho = rhoTarget;
                mtx = mtxTarget;           
                mty = mtyClamp;        
                mtz = mtzTarget;
                eng = engNew;
            }
        }
    }
}

//------------------------------------------------------------
// Z direction, left side 
void Boundary::OutflowZL(TArray<Real>& cons, const Grid& grid){
    int il = grid.ib;                       // first activated cell 
    int ir = grid.ie;                       // first ghost cell right side
    int jl = grid.jb;                       // first activated cell 
    int jr = grid.je;                       // first ghost cell right side
    int kl = grid.kgb;                      // first ghost cell left side 
    int kr = grid.kb;                       // first activated cell 
    const Real denMin = DENSITY_FLOOR; 
#pragma omp parallel for collapse(2) schedule(static)
    // loop z ghost zone  
    for (int k=kl; k<kr; k++){
        // loop y activated zone 
        for (int j=jl; j<jr; j++){
#pragma omp simd
            // loop x activated zone
            for (int i=il; i<ir; i++){
                int kTarget = kr + kr - k - 1; 
                Real rhoTarget = cons(DEN, kTarget, j, i); 
                Real mtxTarget = cons(MTX, kTarget, j, i); 
                Real mtyTarget = cons(MTY, kTarget, j, i); 
                Real mtzTarget = cons(MTZ, kTarget, j, i); 
                Real engTarget = cons(ENG, kTarget, j, i); 
                Real& rho = cons(DEN, k, j, i); 
                Real& mtx = cons(MTX, k, j, i); 
                Real& mty = cons(MTY, k, j, i); 
                Real& mtz = cons(MTZ, k, j, i); 
                Real& eng = cons(ENG, k, j, i); 

                Real mtzClamp = std::min(mtzTarget, 0.0);
                rhoTarget = std::max(rhoTarget, denMin);  
                Real rhoInv = 1.0/rhoTarget; 
                Real keOld = 0.5*rhoInv*mtzTarget*mtzTarget; 
                Real keNew = 0.5*rhoInv*mtzClamp*mtzClamp; 
                Real eint = engTarget - keOld;
                eint = std::max(eint, 1e-16);
                Real engNew = eint + keNew;

                rho = rhoTarget;
                mtx = mtxTarget;           
                mty = mtyTarget;          
                mtz = mtzClamp;      
                eng = engNew;
            }
        }
    }
}

//------------------------------------------------------------
// Z direction, right side 
void Boundary::OutflowZR(TArray<Real>& cons, const Grid& grid){
    int il = grid.ib;                       // first activated cell 
    int ir = grid.ie;                       // first ghost cell right side
    int jl = grid.jb;                       // first activated cell 
    int jr = grid.je;                       // first ghost cell right side
    int kl = grid.ke;                       // first ghost cell right side  
    int kr = grid.kge;                      // last ghost cell right side + 1
    const Real denMin = DENSITY_FLOOR; 
#pragma omp parallel for collapse(2) schedule(static)
    // loop z ghost zone  
    for (int k=kl; k<kr; k++){
        // loop y activated zone 
        for (int j=jl; j<jr; j++){
#pragma omp simd
            // loop x activated zone
            for (int i=il; i<ir; i++){
                int kTarget = kl + kl - k - 1; 
                Real rhoTarget = cons(DEN, kTarget, j, i); 
                Real mtxTarget = cons(MTX, kTarget, j, i); 
                Real mtyTarget = cons(MTY, kTarget, j, i); 
                Real mtzTarget = cons(MTZ, kTarget, j, i); 
                Real engTarget = cons(ENG, kTarget, j, i); 
                Real& rho = cons(DEN, k, j, i); 
                Real& mtx = cons(MTX, k, j, i); 
                Real& mty = cons(MTY, k, j, i); 
                Real& mtz = cons(MTZ, k, j, i); 
                Real& eng = cons(ENG, k, j, i); 

                Real mtzClamp = std::max(mtzTarget, 0.0); 
                rhoTarget = std::max(rhoTarget, denMin); 
                Real rhoInv = 1.0/rhoTarget; 
                Real keOld = 0.5*rhoInv*mtzTarget*mtzTarget; 
                Real keNew = 0.5*rhoInv*mtzClamp*mtzClamp; 
                Real eint = engTarget - keOld;
                eint = std::max(eint, 1e-16);
                Real engNew = eint + keNew;

                rho = rhoTarget;
                mtx = mtxTarget;           
                mty = mtyTarget;          
                mtz = mtzClamp;         
                eng = engNew;
            }
        }
    }
}

//------------------------------------------------------------
} // namespace Gaukuk