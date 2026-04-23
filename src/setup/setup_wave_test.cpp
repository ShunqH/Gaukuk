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
    Real rho0    = Config::getInstance().get("rho0");   
    Real p0      = Config::getInstance().get("p0");    
    // gaussian pulse 
    Real amp     = Config::getInstance().get("amp");   
    Real sigma   = Config::getInstance().get("sigma");   
    Real x0      = Config::getInstance().get("x0");    
    // option
    int  addVel   = Config::getInstance().get("addVel");   

    // sound speed
    Real gamma = Config::getInstance().get("gamma");
    Real c0 = std::sqrt(gamma * p0 / rho0);

    // domain
    Real xmin = domain.xmin;
    Real xmax = domain.xmax;
    Real ymin = domain.ymin;
    Real ymax = domain.ymax;
    Real zmin = domain.zmin;
    Real zmax = domain.zmax;

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
                Real x = domain.xc(i);
                Real y = domain.yc(j);
                Real z = domain.zc(k);

                // gaussian perturbation
                Real dx = x - x0;
                Real perturbation = amp * p0 * std::exp(-dx*dx / (2.0*sigma*sigma));

                // (rho/rho0) = (p/p0)^{1/gamma}
                // rho = rho0 * (1 + perturbation/(gamma*p0))
                Real p = p0 + perturbation;
                Real rho = rho0 * std::pow(p / p0, 1.0/gamma);

                // velocity field, if addVel is on 
                // u = p' / (rho0 * c0)
                Real vx = 0.0;
                if (addVel == 1) {
                    vx = perturbation / (rho0 * c0); 
                }
                Real vy = 0.0;
                Real vz = 0.0;

                prim(DEN, k, j, i) = rho;
                prim(VLX, k, j, i) = vx;
                prim(VLY, k, j, i) = vy;
                prim(VLZ, k, j, i) = vz;
                prim(PRE, k, j, i) = p;
            }
        }
    }

    eos.PrimToCons(prim, cons, grid);
}

} // namespace Gaukuk