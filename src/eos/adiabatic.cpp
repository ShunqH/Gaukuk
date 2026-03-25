#include <cmath> // sqrt()

#include "eos.hpp"
#include "../utils/utils.hpp" // Config 

namespace Gaukuk{
    
// adiabatic equation of state
EquationOfState::EquationOfState(){
    // read adiabatic index gamma from input file
    gamma_ = Config::getInstance().get("gamma") ; 
    // density minimum (floor)
    densityMin_ = Config::getInstance().get("rho_floor", 1e-10) ; 
    // pressure minimum (floor)
    pressureMin_ = Config::getInstance().get("pressure_floor", 1e-10) ; 
}

void EquationOfState::ConsToPrim(TArray<Real>& cons, TArray<Real>& prim, 
                int il, int ir, int jl, int jr, int kl, int kr){
    Real gm1 = gamma_ - 1.0;
#pragma omp parallel for
    for (int k=kl; k<kr; k++){
        for (int j=jl; j<jr; j++){
#pragma omp simd 
            for (int i=il; i<ir; i++){
                Real& consDen = cons(DEN, k, j, i); 
                Real& consMt1 = cons(MT1, k, j, i); 
                Real& consMt2 = cons(MT2, k, j, i); 
                Real& consMt3 = cons(MT3, k, j, i); 
                Real& consEng = cons(ENG, k, j, i); 

                Real& primDen = prim(DEN, k, j, i); 
                Real& primVl1 = prim(VL1, k, j, i); 
                Real& primVl2 = prim(VL2, k, j, i); 
                Real& primVl3 = prim(VL3, k, j, i); 
                Real& primPre = prim(PRE, k, j, i); 

                consDen = (consDen > densityMin_) ? consDen : densityMin_; 
                primDen = consDen; 
                Real denInv = 1.0/consDen; 
                primVl1 = consMt1 * denInv; 
                primVl2 = consMt2 * denInv; 
                primVl3 = consMt3 * denInv; 
                Real engKin = 0.5 * denInv * (consMt1*consMt1 + consMt2*consMt2 + consMt3*consMt3); 
                primPre = gm1 * (consEng - engKin); 
                primPre = (primPre > pressureMin_) ? primPre : pressureMin_; 
                consEng = (primPre > pressureMin_) ? consEng : ( pressureMin_/gm1 + engKin ); 
            }
        }
    }
}

void EquationOfState::PrimToCons(const TArray<Real>& prim, TArray<Real>& cons, 
                int il, int ir, int jl, int jr, int kl, int kr){
    Real gm1Inv = 1.0 / (gamma_ - 1.0);
#pragma omp parallel for
    for (int k=kl; k<kr; k++){
        for (int j=jl; j<jr; j++){
#pragma omp simd 
            for (int i=il; i<ir; i++){
                Real& consDen = cons(DEN, k, j, i); 
                Real& consMt1 = cons(MT1, k, j, i); 
                Real& consMt2 = cons(MT2, k, j, i); 
                Real& consMt3 = cons(MT3, k, j, i); 
                Real& consEng = cons(ENG, k, j, i); 

                const Real& primDen = prim(DEN, k, j, i); 
                const Real& primVl1 = prim(VL1, k, j, i); 
                const Real& primVl2 = prim(VL2, k, j, i); 
                const Real& primVl3 = prim(VL3, k, j, i); 
                const Real& primPre = prim(PRE, k, j, i); 

                consDen = primDen; 
                consMt1 = primDen * primVl1; 
                consMt2 = primDen * primVl2; 
                consMt3 = primDen * primVl3; 
                consEng = primPre * gm1Inv + 0.5 * primDen * ( primVl1*primVl1 + primVl2*primVl2 + primVl3*primVl3 );
            }
        }
    }
}

Real EquationOfState::SoundSpeed(const Real den, const Real pre){
    return std::sqrt(gamma_*pre/den); 
}

} // namespace Gaukuk
