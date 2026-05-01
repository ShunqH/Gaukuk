#pragma once 

// C++ headers

// Gaukuk dependence
#include "../gaukuk.hpp"
#include "../template_array.hpp"
#include "../grid.hpp"

namespace Gaukuk{

struct GravitySource
{
    Real gm; 
    Real x, y, z; 
    Real vx, vy, vz; 
    Real rs; 
    GravitySource(): gm(0), x(0), y(0), z(0), 
                     vx(0), vy(0), vz(0), rs(0) {}
};


class SourceTerm{
public:
    SourceTerm():sourceEnrolled(false), 
                 hasConstAcc(false), 
                 hasPointG(false), 
                 hasBinary(false), 
                 gx_(0.0), gy_(0.0), gz_(0.0){} 
    bool sourceEnrolled; 

    void UpdateSource(TArray<Real>& cons, const Real t,  const Real dt, const Grid& grid, const Domain& domain); 
    void ConstGravity(TArray<Real>& cons, const Real dt, const Grid& grid);
    void PointGravity(TArray<Real>& cons, const Real dt, const Grid& grid, const Domain& domain); 
    void BinaryGravity(TArray<Real>& cons, const Real t, const Real dt, const Grid& grid, const Domain& domain); 
    void InnerWaveKilling(TArray<Real>& cons, const Real dt, const Grid& grid, const Domain& domain);
    
    void EnrollConstGravityVector(Real gx, Real gy, Real gz); 
    void EnrollPointGravity(Real gm, Real x, Real y, Real z, 
                            Real vx, Real vy, Real vz, Real rs); 
private:
    // SourceType srcType;
    bool hasConstAcc, hasPointG, hasBinary; 
    // for constant source
    Real gx_, gy_, gz_; 
    // for binary / star-planet
    GravitySource obj0, obj1; 
    Real omega_, ab_, theta0_, theta1_; 

};


} // namespace Gaukuk