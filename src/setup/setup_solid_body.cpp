// C++ headers
#include <cmath>        // std::exp, M_PI

// Gaukuk dependence
#include "../gaukuk.hpp"
#include "../sim.hpp"

namespace Gaukuk
{

void Sim::Setup()
{
    // background 
    Real rho    = Config::getInstance().get("rho");   
    Real pre    = Config::getInstance().get("p");    
    Real vxIn   = Config::getInstance().get("vxIn"); 
    Real xShell = Config::getInstance().get("x0", 0.0);  
    Real yShell = Config::getInstance().get("y0", 0.0);  
    Real zShell = Config::getInstance().get("z0", 0.0);  
    Real rShell = Config::getInstance().get("r_solid"); 

    // sound speed

    int il = grid.ib;
    int ir = grid.ie;
    int jl = grid.jb;
    int jr = grid.je;
    int kl = grid.kb;
    int kr = grid.ke;

    #pragma omp parallel for collapse(2) schedule(static)
    for (int k = kl; k < kr; ++k) {
        for (int j = jl; j < jr; ++j) {
            #pragma omp simd
            for (int i = il; i < ir; ++i) {
                prim(DEN, k, j, i) = rho;
                prim(VLX, k, j, i) = vxIn;
                prim(VLY, k, j, i) = 0;
                prim(VLZ, k, j, i) = 0;
                prim(PRE, k, j, i) = pre;
            }
        }
    }

    eos.PrimToCons(prim, cons, grid);
    boundary.EnrollSelfDefineBDXL();
    boundary.SetupImmersedShell(grid, domain, xShell, yShell, zShell, rShell); 
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
void Boundary::SelfDefineBDXL(TArray<Real>& cons, const Grid& grid, const EquationOfState& eos){
    int il = grid.igb;                      // first ghost cell left side
    int ir = grid.ib;                       // first activated cell 
    int jl = grid.jb;                       // first activated cell 
    int jr = grid.je;                       // first ghost cell right side
    int kl = grid.kb;                       // first activated cell 
    int kr = grid.ke;                       // first ghost cell right side
    Real rho    = Config::getInstance().get("rho");   
    Real pre    = Config::getInstance().get("p");    
    Real vxIn   = Config::getInstance().get("vxIn"); 
    Real engKin = 0.5*rho* ( vxIn*vxIn ); 
    Real engGas = eos.EGas(rho, pre); 
    Real eng = engKin + engGas; 

#pragma omp parallel for collapse(2) schedule(static)
    for (int k=kl; k<kr; k++){
        // loop y activated zone 
        for (int j=jl; j<jr; j++){
#pragma omp simd
            // loop x ghost zone
            for (int i=il; i<ir; i++){
                Real& consDen = cons(DEN, k, j, i); 
                Real& consMtx = cons(MTX, k, j, i); 
                Real& consMty = cons(MTY, k, j, i); 
                Real& consMtz = cons(MTZ, k, j, i); 
                Real& consEng = cons(ENG, k, j, i); 
                consDen = rho;
                consMtx = rho*vxIn;          
                consMty = 0;
                consMtz = 0;
                consEng = eng;
            }
        }
    }
}

} // namespace Gaukuk