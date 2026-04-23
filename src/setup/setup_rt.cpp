// C++ headers
#include <cmath>        //std::sin

// Gaukuk dependence
#include "../gaukuk.hpp" 
#include "../sim.hpp" 

namespace Gaukuk
{
    
void Sim::Setup() {
    Real gamma = Config::getInstance().get("gamma");
    Real rhoUp = Config::getInstance().get("rhoUp");
    Real rhoDown = Config::getInstance().get("rhoDown");
    Real amp = Config::getInstance().get("amp");
    Real gravAccY = Config::getInstance().get("gravAccY"); 
    Real nxPtb = Config::getInstance().get("nx_ptb", 2.0); 
    Real nyPtb = Config::getInstance().get("ny_ptb", 2.0); 
    Real gy = Config::getInstance().get("gravAccY", -0.1); 

    Real xmin = domain.xmin;
    Real xmax = domain.xmax;
    Real ymin = domain.ymin;
    Real ymax = domain.ymax;
    Real Lx = xmax - xmin;
    Real Ly = ymax - ymin;

    Real kx = nxPtb * PI / Lx;
    Real ky = nyPtb * PI / Ly;

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
                Real xNow = domain.xc(i);
                Real yNow = domain.yc(j);

                Real rhoNow = (yNow > 0.0) ? rhoUp : rhoDown;

                Real shape = (1.0 + std::cos(kx * xNow)) *
                             (1.0 + std::cos(ky * yNow)) / 4.0;
                Real vx = 0.0;
                Real vy = amp * shape;
                Real vz = 0.0;

                Real pressure = 1.0 / gamma + gravAccY * rhoNow * yNow;

                prim(DEN, k, j, i) = rhoNow;
                prim(VLX, k, j, i) = vx;
                prim(VLY, k, j, i) = vy;
                prim(VLZ, k, j, i) = vz;
                prim(PRE, k, j, i) = pressure;
            }
        }
    }
    // Enroll Const Gravity as a Vector
    // void SourceTerm::EnrollConstGravityVector(Real gx, Real gy, Real gz)
    // This will enable the constant gravity source term 
    srcTerm.EnrollConstGravityVector(0, gy, 0); 
    // if primitive quantivities is set, make sure you call eos.PrimToCons 
    eos.PrimToCons(prim, cons, grid);
}

} // namespace Gaukuk
