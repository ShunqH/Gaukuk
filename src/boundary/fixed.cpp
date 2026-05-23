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
void Boundary::FixedBDXL(TArray<Real>& cons, const Grid& grid, const EquationOfState& eos)
{
    if (!isFixedBdEnrolled)
        throw std::invalid_argument("fixed boundary is not initialized. call 'Boundary::EnrollFixedBoundaryData(const TArray<Real>& cons, const Grid& grid)' in your setup");

    const int il = grid.igb;   // 0
    const int ir = grid.ib;    // nGhost
    const int jl = grid.jb;    // nGhost
    const int jr = grid.je;    // ny + nGhost
    const int kl = grid.kb;    // nGhost
    const int kr = grid.ke;    // nz + nGhost

    #pragma omp parallel for collapse(3) schedule(static)
    for (int ivar = DEN; ivar <= ENG; ++ivar) {
        for (int k = kl; k < kr; ++k) {
            for (int j = jl; j < jr; ++j) {
                #pragma omp simd
                for (int i = il; i < ir; ++i) {
                    cons(ivar, k, j, i) = consFixedBdxl(ivar, k, j, i);
                }
            }
        }
    }
}

//------------------------------------------------------------
// X direction, right side
void Boundary::FixedBDXR(TArray<Real>& cons, const Grid& grid, const EquationOfState& eos)
{
    if (!isFixedBdEnrolled)
        throw std::invalid_argument("fixed boundary is not initialized. call 'Boundary::EnrollFixedBoundaryData(const TArray<Real>& cons, const Grid& grid)' in your setup");

    const int il = grid.ie;    // nx + nGhost
    const int ir = grid.ige;   // nx + 2*nGhost
    const int jl = grid.jb;
    const int jr = grid.je;
    const int kl = grid.kb;
    const int kr = grid.ke;

    #pragma omp parallel for collapse(3) schedule(static)
    for (int ivar = DEN; ivar <= ENG; ++ivar) {
        for (int k = kl; k < kr; ++k) {
            for (int j = jl; j < jr; ++j) {
                #pragma omp simd
                for (int i = il; i < ir; ++i) {
                    const int iLocal = i - grid.ie;   // 0 .. nGhost-1
                    cons(ivar, k, j, i) = consFixedBdxr(ivar, k, j, iLocal);
                }
            }
        }
    }
}

//------------------------------------------------------------
// Y direction, left side
void Boundary::FixedBDYL(TArray<Real>& cons, const Grid& grid, const EquationOfState& eos)
{
    if (!isFixedBdEnrolled)
        throw std::invalid_argument("fixed boundary is not initialized. call 'Boundary::EnrollFixedBoundaryData(const TArray<Real>& cons, const Grid& grid)' in your setup");

    const int il = grid.ib;    // nGhost
    const int ir = grid.ie;    // nx + nGhost
    const int jl = grid.jgb;   // 0
    const int jr = grid.jb;    // nGhost
    const int kl = grid.kb;
    const int kr = grid.ke;

    #pragma omp parallel for collapse(3) schedule(static)
    for (int ivar = DEN; ivar <= ENG; ++ivar) {
        for (int k = kl; k < kr; ++k) {
            for (int j = jl; j < jr; ++j) {
                #pragma omp simd
                for (int i = il; i < ir; ++i) {
                    cons(ivar, k, j, i) = consFixedBdyl(ivar, k, j, i);
                }
            }
        }
    }
}

//------------------------------------------------------------
// Y direction, right side
void Boundary::FixedBDYR(TArray<Real>& cons, const Grid& grid, const EquationOfState& eos)
{
    if (!isFixedBdEnrolled)
        throw std::invalid_argument("fixed boundary is not initialized. call 'Boundary::EnrollFixedBoundaryData(const TArray<Real>& cons, const Grid& grid)' in your setup");

    const int il = grid.ib;
    const int ir = grid.ie;
    const int jl = grid.je;    // ny + nGhost
    const int jr = grid.jge;   // ny + 2*nGhost
    const int kl = grid.kb;
    const int kr = grid.ke;

    #pragma omp parallel for collapse(3) schedule(static)
    for (int ivar = DEN; ivar <= ENG; ++ivar) {
        for (int k = kl; k < kr; ++k) {
            for (int j = jl; j < jr; ++j) {
                #pragma omp simd
                for (int i = il; i < ir; ++i) {
                    const int jLocal = j - grid.je;   // 0 .. nGhost-1
                    cons(ivar, k, j, i) = consFixedBdyr(ivar, k, jLocal, i);
                }
            }
        }
    }
}

//------------------------------------------------------------
// Z direction, left side
void Boundary::FixedBDZL(TArray<Real>& cons, const Grid& grid, const EquationOfState& eos)
{
    if (!isFixedBdEnrolled)
        throw std::invalid_argument("fixed boundary is not initialized. call 'Boundary::EnrollFixedBoundaryData(const TArray<Real>& cons, const Grid& grid)' in your setup");

    const int il = grid.ib;
    const int ir = grid.ie;
    const int jl = grid.jb;
    const int jr = grid.je;
    const int kl = grid.kgb;   // 0
    const int kr = grid.kb;    // nGhost

    #pragma omp parallel for collapse(3) schedule(static)
    for (int ivar = DEN; ivar <= ENG; ++ivar) {
        for (int k = kl; k < kr; ++k) {
            for (int j = jl; j < jr; ++j) {
                #pragma omp simd
                for (int i = il; i < ir; ++i) {
                    cons(ivar, k, j, i) = consFixedBdzl(ivar, k, j, i);
                }
            }
        }
    }
}

//------------------------------------------------------------
// Z direction, right side
void Boundary::FixedBDZR(TArray<Real>& cons, const Grid& grid, const EquationOfState& eos)
{
    if (!isFixedBdEnrolled)
        throw std::invalid_argument("fixed boundary is not initialized. call 'Boundary::EnrollFixedBoundaryData(const TArray<Real>& cons, const Grid& grid)' in your setup");

    const int il = grid.ib;
    const int ir = grid.ie;
    const int jl = grid.jb;
    const int jr = grid.je;
    const int kl = grid.ke;    // nz + nGhost
    const int kr = grid.kge;   // nz + 2*nGhost

    #pragma omp parallel for collapse(3) schedule(static)
    for (int ivar = DEN; ivar <= ENG; ++ivar) {
        for (int k = kl; k < kr; ++k) {
            for (int j = jl; j < jr; ++j) {
                #pragma omp simd
                for (int i = il; i < ir; ++i) {
                    const int kLocal = k - grid.ke;   // 0 .. nGhost-1
                    cons(ivar, k, j, i) = consFixedBdzr(ivar, kLocal, j, i);
                }
            }
        }
    }
}

//------------------------------------------------------------
// fix inner 
void Boundary::FixedBDInner(TArray<Real>& cons, const Grid& grid, const Domain& domain, const EquationOfState& eos)
{
    const int il = static_cast<int>(grid.ib + 0.25*(grid.ie - grid.ib));
    const int ir = static_cast<int>(grid.ib + 0.75*(grid.ie - grid.ib));
    const int jl = static_cast<int>(grid.jb + 0.25*(grid.je - grid.jb));
    const int jr = static_cast<int>(grid.jb + 0.75*(grid.je - grid.jb));
    const int kl = grid.kb;    // nz + nGhost
    const int kr = grid.ke;   // nz + 2*nGhost

    #pragma omp parallel for collapse(3) schedule(static)
    for (int ivar = DEN; ivar <= ENG; ++ivar) {
        for (int k = kl; k < kr; ++k) {
            for (int j = jl; j < jr; ++j) {
                #pragma omp simd
                for (int i = il; i < ir; ++i) {
                    Real xnow = domain.xc(i); 
                    Real ynow = domain.yc(j); 
                    Real znow = domain.zc(k); 

                    Real r2 = xnow*xnow + ynow*ynow + znow*znow; 
                    if (r2<0.5*0.5){
                        cons(ivar, k, j, i) = consFixedInner(ivar, k, j, i);
                    }
                }
            }
        }
    }
}
//------------------------------------------------------------


namespace{

void CopyBD(const TArray<Real>& cons, TArray<Real>& consFixedBd, int il, int ir, int jl, int jr, int kl, int kr, int n, int axis, int side){
    #pragma omp parallel for collapse(3) schedule(static)
    for (int ivar=DEN; ivar<=ENG; ivar++){
        // loop z ghost zone  
        for (int k=kl; k<kr; k++){
            // loop y activated zone 
            for (int j=jl; j<jr; j++){
                // loop x activated zone
                for (int i=il; i<ir; i++){
                    int iTar = i; 
                    int jTar = j; 
                    int kTar = k; 
                    if (side == 1){
                        iTar = (axis==0) ? i+n : i; 
                        jTar = (axis==1) ? j+n : j; 
                        kTar = (axis==2) ? k+n : k; 
                    }
                    Real uTarget = cons(ivar, kTar, jTar, iTar); 
                    Real& u = consFixedBd(ivar, k, j, i); 
                    u = uTarget;
                }
            }
        }
    }
}

}

void Boundary::EnrollFixedBoundaryData(const TArray<Real>& cons, const Grid& grid, bool fixInner){
    isFixedBdEnrolled = true; 
    consFixedBdxl.NewArray(NVar, grid.lenz, grid.leny, grid.nGhost);
    consFixedBdxr.NewArray(NVar, grid.lenz, grid.leny, grid.nGhost);
    consFixedBdyl.NewArray(NVar, grid.lenz, grid.nGhost, grid.lenx);
    consFixedBdyr.NewArray(NVar, grid.lenz, grid.nGhost, grid.lenx);
    if (grid.nz>1){
        consFixedBdzl.NewArray(NVar, grid.nGhost, grid.leny, grid.lenx);
        consFixedBdzr.NewArray(NVar, grid.nGhost, grid.leny, grid.lenx);
    }

    CopyBD(cons, consFixedBdxl, grid.igb, grid.ib, grid.jgb, grid.jge, grid.kgb, grid.kge, grid.nx+grid.nGhost, 0, 0); 
    CopyBD(cons, consFixedBdxr, grid.igb, grid.ib, grid.jgb, grid.jge, grid.kgb, grid.kge, grid.nx+grid.nGhost, 0, 1); 

    CopyBD(cons, consFixedBdyl, grid.igb, grid.ige, grid.jgb, grid.jb, grid.kgb, grid.kge, grid.ny+grid.nGhost, 1, 0); 
    CopyBD(cons, consFixedBdyr, grid.igb, grid.ige, grid.jgb, grid.jb, grid.kgb, grid.kge, grid.ny+grid.nGhost, 1, 1); 

    if (grid.nz>1){
        CopyBD(cons, consFixedBdzl, grid.igb, grid.ige, grid.jgb, grid.jge, grid.kgb, grid.kb, grid.nz+grid.nGhost, 2, 0); 
        CopyBD(cons, consFixedBdzr, grid.igb, grid.ige, grid.jgb, grid.jge, grid.kgb, grid.kb, grid.nz+grid.nGhost, 2, 1); 
    }

    if (fixInner){
        isFixedInner = true; 
        consFixedInner = cons;  
    }

} 


} // namespace Gaukuk