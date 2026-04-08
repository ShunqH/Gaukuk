#pragma once 

// CPP header
#include <stdexcept>    // runtime_error

// Gaukuk dependence
#include "../utils/utils.hpp"   // Config

namespace Gaukuk{

class Grid{
public:
friend class Sim; 
    int nx, ny, nz, nGhost, lenx, leny, lenz, lenArr; 
    int igb, ige, jgb, jge, kgb, kge;   // begin and end ids including ghost cells 
    int ib, ie, jb, je, kb, ke;         // activate cell begin and end ids 
    int iib, iie, jjb, jje, kkb, kke;   // flux (half) cell begin and end ids

    Grid() {
        nx = Config::getInstance().get("NX"); 
        ny = Config::getInstance().get("NY"); 
        nz = Config::getInstance().get("NZ"); 
        nGhost = Config::getInstance().get("NGhost"); 
        if (nx <= 0 || ny <= 0 || nz <= 0 || nGhost < 0){
            throw std::runtime_error("Invalid grid parameters from config");
        }
        lenx = nx + 2 * nGhost;
        leny = ny + 2 * nGhost;
        lenz = nz + 2 * nGhost;
        lenArr = lenx * leny * lenz;

        igb = jgb = kgb = 0; 
        ige = lenx; 
        jge = leny; 
        kge = lenz;
        ib = igb + nGhost; 
        ie = ib + nx; 
        jb = jgb + nGhost; 
        je = jb + ny; 
        kb = kgb + nGhost; 
        ke = kb + nz; 

        iib = jjb = kkb = 0; 
        iie = lenx-1; 
        jje = leny-1; 
        kke = lenz-1;  
    }
};

} // namespace Gaukuk