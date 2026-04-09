#pragma once 

// Gaukuk dependence
#include "gaukuk.hpp"
#include "../template_array.hpp"

namespace Gaukuk
{
    
void CheckNaN(const TArray<Real>& arr, const char* name); 
void WriteTarray(const TArray<Real>&arr, const char* filename, const int ID); 

} // namespace Gaukuk
