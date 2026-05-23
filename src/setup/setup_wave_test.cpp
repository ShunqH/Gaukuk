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

    // period lengths
    Real Lx = xmax - xmin;
    Real Ly = ymax - ymin;

    // center in y (domain mid-point, can be changed)
    Real y0 = 0.5 * (ymin + ymax);

    // direction cosines for 45° upper-right (1,1,0)
    const Real inv_sqrt2 = 1.0 / std::sqrt(2.0);

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

                // accumulate contributions from periodic mirror images
                // (here we assume sigma is small enough that i,j in {-1,0,1} is sufficient;
                //  increase the range if sigma is comparable to Lx or Ly)
                Real perturbation = 0.0;
                Real vx = 0.0;
                Real vy = 0.0;
                Real vz = 0.0;

                for (int iper = -1; iper <= 1; ++iper) {
                    for (int jper = -1; jper <= 1; ++jper) {
                        Real xc = x0 + iper * Lx;
                        Real yc = y0 + jper * Ly;

                        // distance along the 45° direction
                        Real dr = (x - xc) * inv_sqrt2 + (y - yc) * inv_sqrt2;

                        // gaussian perturbation
                        Real pert_ij = amp * p0 * std::exp(-dr*dr / (2.0*sigma*sigma));
                        perturbation += pert_ij;

                        if (addVel == 1) {
                            Real v_amp = pert_ij / (rho0 * c0);
                            vx += v_amp * inv_sqrt2;
                            vy += v_amp * inv_sqrt2;
                        }
                    }
                }

                // pressure and density (isentropic relation)
                Real p = p0 + perturbation;
                Real rho = rho0 * std::pow(p / p0, 1.0/gamma);

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