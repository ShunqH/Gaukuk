#pragma once 

// C++ headers 
#include <iostream>     // std::cout; std::endl; std::cerr

// Gaukuk dependence
#include "../gaukuk.hpp" 
#include "../template_array.hpp"

namespace Gaukuk{

class Reconstruction{
// friend class Sim; 
public:
    void ReconstructXFirstOrder(const TArray<Real>& prim, 
                             TArray<Real>& ul, TArray<Real>& ur, 
                             int k, int j, int igb, int ige); 
    void ReconstructYZFirstOrder(const TArray<Real>& prim, TArray<Real>& ur, 
                              int k, int j, int igb, int ige); 

    // has to be uniform grid 
    void ReconstructXPLM(const TArray<Real>& prim, 
                         TArray<Real>& ul, TArray<Real>& ur, 
                         int k, int j, int igb, int ige);
    void ReconstructYZPLM(const TArray<Real>& prim,
                          TArray<Real>& ur, TArray<Real>& ulNext, 
                          int k, int j, int igb, int ige); 

private:
    inline static Real Minmod(Real a, Real b){
        Real multi = a*b; 
        Real sigma = (std::abs(a) < std::abs(b)) ? a : b;
        return (multi <= 0.0) ? 0.0 : sigma;
    }

    inline Real VanLeer(Real a, Real b) {
        Real multi = a*b; 
        Real sigma = 2.0 * multi / (a + b); 
        return (multi <= 0.0) ? 0.0 : sigma;
        // if (a*b <= 0.0) return 0.0;
        // return (2.0 * a * b) / (a + b);
    }

    inline Real MC(Real a, Real b) {
        return Minmod(Minmod(2*a, 2*b), 0.5*(a+b));
    }
};

} // namespace Gaukuk