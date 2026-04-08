#pragma once 

namespace Gaukuk{

using Real = double; 

constexpr int NVar = 5;

enum ConsIDs {DEN = 0, MTX = 1, MTY = 2, MTZ = 3, ENG = 4, }; 
enum PrimIDs {VLX = 1, VLY = 2, VLZ = 3, PRE = 4, }; 

} // namespace Gaukuk