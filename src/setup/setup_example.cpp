// C++ headers
#include <cmath>        //std::sin

// Gaukuk dependence
#include "../gaukuk.hpp" 
#include "../sim.hpp" 

namespace Gaukuk
{
    
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
// ***** problem setup example *****
//
// I used the kh instability as an example 
// to show how to setup the problem 
//--------------------------------------------------------------------------
void Sim::Setup(){
    // load parameter from setupfile 
    Real rhoIn = Config::getInstance().get("rhoIn"); 
    Real rhoOut = Config::getInstance().get("rhoOut"); 
    Real vxIn = Config::getInstance().get("vxIn"); 
    Real vxOut = Config::getInstance().get("vxOut"); 
    Real pressure = Config::getInstance().get("pressure"); 
    Real amp = Config::getInstance().get("amp"); 

    /* //------------------------------------------------------
    sim.domain includes the information below
    class Domain{
        ----- maximum and minimum of the domain ---------
        Real xmin, xmax, ymin, ymax, zmin, zmax 
        ----- cell size ---------------------------------
        Real dx, dy, dz, drmin
        ----- x, y, z axis, located at cell center ------
        TArray<Real> xc, yc, zc 
    }
    */ //------------------------------------------------------
    Real yLayer = domain.ymax/2; 
    Real xmin = domain.xmin; 
    Real Lx = domain.xmax - domain.xmin; 
    Real ymin = domain.ymin; 
    Real Ly = domain.ymax - domain.ymin; 
    
    /* //------------------------------------------------------
    sim.grid includes the grid information below
    class Grid{
        ----- array lengths / resolutions ---------------
        int nx, ny, nz -> activated cell lengths
        int nGhost     -> number of ghost cell
        int lenx, leny, lenz    
                       -> length including ghost cells
                          lenx = nx + 2 * nGhost
        int lenArr     -> cons and prim length
                          lenArr = lenx*leny*lenz
        ----- index -------------------------------------
        int igb, ige, jgb, jge, kgb, kge
        int ib,  ie,  jb,  je,  kb,  ke
            igb     -> i ghost begin
            ige     -> i ghost end 
            ib      -> i begin 
            ie      -> i end 
        alignment: 
         igb + nGhost -> ib + nx  -> ie + nGhost   -> ige
        |     nGhost     |    nx    |     nGhost     |
    }
    */ //------------------------------------------------------
    // loop for activated zone 
    int il = grid.ib; 
    int ir = grid.ie; 
    int jl = grid.jb; 
    int jr = grid.je; 
    int kl = grid.kb;
    int kr = grid.ke; 

    // setup simulation
#pragma omp parallel for collapse(2) schedule(static)    
    for (int k=kl; k<kr; k++){
        for (int j=jl; j<jr; j++){
#pragma omp simd
            for (int i=il; i<ir; i++){
                Real xNow = domain.xc(i); 
                Real yNow = domain.yc(j); 

                // setup shear stream with perturbations
                Real rhoNow = (std::abs(yNow)<yLayer) ? rhoIn : rhoOut ; 
                Real vx = (std::abs(yNow)<yLayer) ? vxIn : vxOut ; 

                Real vxRandom = amp * std::sin( 4 * PI * (xNow-xmin) / Lx )
                                    * std::sin( 2 * PI * (yNow-ymin) / Ly ) ; 
                Real vyRandom = amp * std::sin( 4 * PI * (xNow-xmin) / Lx + 0.87*PI)
                                    * std::sin( 2 * PI * (yNow-ymin) / Ly + 1.23*PI) ; 

                // usually you have to setup conservative quantivities (cons) 
                // but you can setup primitive quantivities (cons) 
                // then use eos.PrimToCons to conver 
                prim(DEN, k, j, i) = rhoNow; 
                prim(VLX, k, j, i) = vx + vxRandom; 
                prim(VLY, k, j, i) = vyRandom; 
                prim(VLZ, k, j, i) = 0; 
                prim(PRE, k, j, i) = pressure; 
            }
        }
    }
    // if primitive quantivities is set, make sure you call eos.PrimToCons 
    eos.PrimToCons(prim, cons, grid); 

    // Make sure you enroll the boundary condition if you use selfdefine
    // boundary.EnrollSelfDefineBDXL();
    // boundary.EnrollSelfDefineBDXR();
    // boundary.EnrollSelfDefineBDYL();
    // boundary.EnrollSelfDefineBDYR();
    // boundary.EnrollSelfDefineBDZL();
    // boundary.EnrollSelfDefineBDZR();
}

//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
// ***** user self define boundary condition *****
//
// You can define the boundary condition yourself 
// Implement the function below and enroll in the setup function
// Once the bd function is enroll, the option in input file will exprie
// Below is the copy outflow boundary as a template 
// Copy the edge cell of the activated zone 
// to all the ghost cell next to it
//--------------------------------------------------------------------------
// X direction, left side 
void Boundary::SelfDefineBDXL(TArray<Real>& cons, const Grid& grid){
    int il = grid.igb;                      // first ghost cell left side
    int ir = grid.ib;                       // first activated cell 
    int jl = grid.jb;                       // first activated cell 
    int jr = grid.je;                       // first ghost cell right side
    int kl = grid.kb;                       // first activated cell 
    int kr = grid.ke;                       // first ghost cell right side
    int iAct = ir;                          // copy cell's id
#pragma omp parallel for collapse(2) schedule(static)
    for (int k=kl; k<kr; k++){
        // loop y activated zone 
        for (int j=jl; j<jr; j++){
#pragma omp simd
            // loop x ghost zone
            for (int i=il; i<ir; i++){
                const Real& rhoTarget = cons(DEN, k, j, iAct); 
                const Real& mtxTarget = cons(MTX, k, j, iAct); 
                const Real& mtyTarget = cons(MTY, k, j, iAct); 
                const Real& mtzTarget = cons(MTZ, k, j, iAct); 
                const Real& engTarget = cons(ENG, k, j, iAct); 
                Real& rho = cons(DEN, k, j, i); 
                Real& mtx = cons(MTX, k, j, i); 
                Real& mty = cons(MTY, k, j, i); 
                Real& mtz = cons(MTZ, k, j, i); 
                Real& eng = cons(ENG, k, j, i); 
                rho = rhoTarget;
                mtx = mtxTarget;          
                mty = mtyTarget;
                mtz = mtzTarget;
                eng = engTarget;
            }
        }
    }
}

//-------------------------------------------------------------------------
// X direction, right side 
void Boundary::SelfDefineBDXR(TArray<Real>& cons, const Grid& grid){
    int il = grid.ie;                       // first ghost cell right side
    int ir = grid.ige;                      // last ghost cell right side + 1
    int jl = grid.jb;                       // first activated cell 
    int jr = grid.je;                       // first ghost cell right side
    int kl = grid.kb;                       // first activated cell 
    int kr = grid.ke;                       // first ghost cell right side
    int iAct = il - 1;                      // copy cell's id
#pragma omp parallel for collapse(2) schedule(static)
    // loop z activated zone  
    for (int k=kl; k<kr; k++){
        // loop y activated zone 
        for (int j=jl; j<jr; j++){
#pragma omp simd
            // loop x ghost zone
            for (int i=il; i<ir; i++){
                const Real& rhoTarget = cons(DEN, k, j, iAct); 
                const Real& mtxTarget = cons(MTX, k, j, iAct); 
                const Real& mtyTarget = cons(MTY, k, j, iAct); 
                const Real& mtzTarget = cons(MTZ, k, j, iAct); 
                const Real& engTarget = cons(ENG, k, j, iAct); 
                Real& rho = cons(DEN, k, j, i); 
                Real& mtx = cons(MTX, k, j, i); 
                Real& mty = cons(MTY, k, j, i); 
                Real& mtz = cons(MTZ, k, j, i); 
                Real& eng = cons(ENG, k, j, i); 
                rho = rhoTarget;
                mtx = mtxTarget;          
                mty = mtyTarget;
                mtz = mtzTarget;
                eng = engTarget;
            }
        }
    }
}

//-------------------------------------------------------------------------
// Y direction, left side 
void Boundary::SelfDefineBDYL(TArray<Real>& cons, const Grid& grid){
    int il = grid.ib;                       // first activated cell 
    int ir = grid.ie;                       // first ghost cell right side
    int jl = grid.jgb;                      // first ghost cell left side 
    int jr = grid.jb;                       // first activated cell 
    int kl = grid.kb;                       // first activated cell 
    int kr = grid.ke;                       // first ghost cell right side
    int jAct = jr;                          // copy cell's id
#pragma omp parallel for collapse(2) schedule(static)
    // loop z activated zone  
    for (int k=kl; k<kr; k++){
        // loop y ghost zone 
        for (int j=jl; j<jr; j++){
#pragma omp simd
            // loop x activated zone
            for (int i=il; i<ir; i++){
                const Real& rhoTarget = cons(DEN, k, jAct, i); 
                const Real& mtxTarget = cons(MTX, k, jAct, i); 
                const Real& mtyTarget = cons(MTY, k, jAct, i); 
                const Real& mtzTarget = cons(MTZ, k, jAct, i); 
                const Real& engTarget = cons(ENG, k, jAct, i); 
                Real& rho = cons(DEN, k, j, i); 
                Real& mtx = cons(MTX, k, j, i); 
                Real& mty = cons(MTY, k, j, i); 
                Real& mtz = cons(MTZ, k, j, i); 
                Real& eng = cons(ENG, k, j, i); 
                rho = rhoTarget;
                mtx = mtxTarget;           
                mty = mtyTarget;         
                mtz = mtzTarget;
                eng = engTarget;
            }
        }
    }
}

//-------------------------------------------------------------------------
// Y direction, right side 
void Boundary::SelfDefineBDYR(TArray<Real>& cons, const Grid& grid){
    int il = grid.ib;                       // first activated cell 
    int ir = grid.ie;                       // first ghost cell right side
    int jl = grid.je;                       // first ghost cell right side 
    int jr = grid.jge;                      // last ghost cell right side + 1
    int kl = grid.kb;                       // first activated cell 
    int kr = grid.ke;                       // first ghost cell right side
    int jAct = jl - 1;                      // copy cell's id
#pragma omp parallel for collapse(2) schedule(static)
    // loop z activated zone  
    for (int k=kl; k<kr; k++){
        // loop y ghost zone 
        for (int j=jl; j<jr; j++){
#pragma omp simd
            // loop x activated zone
            for (int i=il; i<ir; i++){
                const Real& rhoTarget = cons(DEN, k, jAct, i); 
                const Real& mtxTarget = cons(MTX, k, jAct, i); 
                const Real& mtyTarget = cons(MTY, k, jAct, i); 
                const Real& mtzTarget = cons(MTZ, k, jAct, i); 
                const Real& engTarget = cons(ENG, k, jAct, i); 
                Real& rho = cons(DEN, k, j, i); 
                Real& mtx = cons(MTX, k, j, i); 
                Real& mty = cons(MTY, k, j, i); 
                Real& mtz = cons(MTZ, k, j, i); 
                Real& eng = cons(ENG, k, j, i); 
                rho = rhoTarget;
                mtx = mtxTarget;           
                mty = mtyTarget;     
                mtz = mtzTarget;
                eng = engTarget;
            }
        }
    }
}

//-------------------------------------------------------------------------
// Z direction, left side 
void Boundary::SelfDefineBDZL(TArray<Real>& cons, const Grid& grid){
    int il = grid.ib;                       // first activated cell 
    int ir = grid.ie;                       // first ghost cell right side
    int jl = grid.jb;                       // first activated cell 
    int jr = grid.je;                       // first ghost cell right side
    int kl = grid.kgb;                      // first ghost cell left side 
    int kr = grid.kb;                       // first activated cell 
    int kAct = kr;                          // copy cell's id
#pragma omp parallel for collapse(2) schedule(static)
    // loop z ghost zone  
    for (int k=kl; k<kr; k++){
        // loop y activated zone 
        for (int j=jl; j<jr; j++){
#pragma omp simd
            // loop x activated zone
            for (int i=il; i<ir; i++){
                const Real& rhoTarget = cons(DEN, kAct, j, i); 
                const Real& mtxTarget = cons(MTX, kAct, j, i); 
                const Real& mtyTarget = cons(MTY, kAct, j, i); 
                const Real& mtzTarget = cons(MTZ, kAct, j, i); 
                const Real& engTarget = cons(ENG, kAct, j, i); 
                Real& rho = cons(DEN, k, j, i); 
                Real& mtx = cons(MTX, k, j, i); 
                Real& mty = cons(MTY, k, j, i); 
                Real& mtz = cons(MTZ, k, j, i); 
                Real& eng = cons(ENG, k, j, i); 
                rho = rhoTarget;
                mtx = mtxTarget;           
                mty = mtyTarget;          
                mtz = mtzTarget;       
                eng = engTarget;
            }
        }
    }
}

//-------------------------------------------------------------------------
// Z direction, right side 
void Boundary::SelfDefineBDZR(TArray<Real>& cons, const Grid& grid){
    int il = grid.ib;                       // first activated cell 
    int ir = grid.ie;                       // first ghost cell right side
    int jl = grid.jb;                       // first activated cell 
    int jr = grid.je;                       // first ghost cell right side
    int kl = grid.ke;                       // first ghost cell right side  
    int kr = grid.kge;                      // last ghost cell right side + 1
    int kAct = kl - 1;                      // copy cell's id
#pragma omp parallel for collapse(2) schedule(static)
    // loop z ghost zone  
    for (int k=kl; k<kr; k++){
        // loop y activated zone 
        for (int j=jl; j<jr; j++){
#pragma omp simd
            // loop x activated zone
            for (int i=il; i<ir; i++){
                const Real& rhoTarget = cons(DEN, kAct, j, i); 
                const Real& mtxTarget = cons(MTX, kAct, j, i); 
                const Real& mtyTarget = cons(MTY, kAct, j, i); 
                const Real& mtzTarget = cons(MTZ, kAct, j, i); 
                const Real& engTarget = cons(ENG, kAct, j, i); 
                Real& rho = cons(DEN, k, j, i); 
                Real& mtx = cons(MTX, k, j, i); 
                Real& mty = cons(MTY, k, j, i); 
                Real& mtz = cons(MTZ, k, j, i); 
                Real& eng = cons(ENG, k, j, i); 
                rho = rhoTarget;
                mtx = mtxTarget;           
                mty = mtyTarget;          
                mtz = mtzTarget;         
                eng = engTarget;
            }
        }
    }
}

} // namespace Gaukuk
