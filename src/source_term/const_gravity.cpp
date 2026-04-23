// C++ headers
#include <cmath>            // sqrt(), abs()
#include <algorithm>        // max()  

// Gaukuk dependence
#include "source.hpp"

namespace Gaukuk{
    
void SourceTerm::ConstGravity(TArray<Real>& cons, const Real dt, const Grid& grid){
    int il = grid.ib;               // first activated cell left side
    int ir = grid.ie;               // last activated cell right side + 1
    int jl = grid.jb;               // first activated cell left side
    int jr = grid.je;               // last activated cell right side + 1
    int kl = grid.kb;               // first activated cell left side
    int kr = grid.ke;               // last activated cell right side + 1

    if (gx_!=0){
    #pragma omp parallel for collapse(2) schedule(static)
        for (int k=kl; k<kr; k++){
            for (int j=jl; j<jr; j++){
    #pragma omp simd 
                for (int i=il; i<ir; i++){
                    Real& consDen = cons(DEN, k, j, i); 
                    Real& consMtx = cons(MTX, k, j, i); 
                    Real& consEng = cons(ENG, k, j, i); 

                    consEng += consMtx*gx_*dt;
                    consMtx += consDen*gx_*dt; 
                }
            }
        }
    }

    if (gy_!=0){
    #pragma omp parallel for collapse(2) schedule(static)
        for (int k=kl; k<kr; k++){
            for (int j=jl; j<jr; j++){
    #pragma omp simd 
                for (int i=il; i<ir; i++){
                    Real& consDen = cons(DEN, k, j, i); 
                    Real& consMty = cons(MTY, k, j, i); 
                    Real& consEng = cons(ENG, k, j, i); 

                    consEng += consMty*gy_*dt;
                    consMty += consDen*gy_*dt; 
                }
            }
        }
    }

    if (gz_!=0){
    #pragma omp parallel for collapse(2) schedule(static)
        for (int k=kl; k<kr; k++){
            for (int j=jl; j<jr; j++){
    #pragma omp simd 
                for (int i=il; i<ir; i++){
                    Real& consDen = cons(DEN, k, j, i); 
                    Real& consMtz = cons(MTZ, k, j, i); 
                    Real& consEng = cons(ENG, k, j, i); 

                    consEng += consMtz*gz_*dt;
                    consMtz += consDen*gz_*dt; 
                }
            }
        }
    }
}

}