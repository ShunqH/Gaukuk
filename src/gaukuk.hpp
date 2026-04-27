#pragma once 

namespace Gaukuk{

#ifdef GAUKUK_USE_FLOAT
    using Real = float;
#else
    using Real = double;
#endif

#ifndef GAUKUK_DENSITY_FLOOR
#define GAUKUK_DENSITY_FLOOR 1e-16
#endif
#ifndef GAUKUK_PRESSURE_FLOOR
#define GAUKUK_PRESSURE_FLOOR 1e-16
#endif

constexpr int NVar = 5;
constexpr Real PI = 3.1415926535; 
constexpr Real DENSITY_FLOOR      = GAUKUK_DENSITY_FLOOR;
constexpr Real PRESSURE_FLOOR = GAUKUK_PRESSURE_FLOOR;

enum ConsIDs {DEN = 0, MTX = 1, MTY = 2, MTZ = 3, ENG = 4, }; 
enum PrimIDs {VLX = 1, VLY = 2, VLZ = 3, PRE = 4, }; 

} // namespace Gaukuk