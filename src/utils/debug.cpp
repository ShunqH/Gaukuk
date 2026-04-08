// C++ Header 
#include <cmath>
#include <stdexcept>
#include <iostream>

// Gaukuk dependence
#include "debug.hpp"

namespace Gaukuk
{
    
void CheckNaN(const TArray<Real>& arr, const char* name) {
    int nArr = arr.GetSize(); 
    for (int i=0; i<nArr; i++){
        Real v = arr(nArr);
            if (std::isnan(v) || std::isinf(v)) {
                std::cerr << "NaN detected in " << name
                            << " value = " << v << std::endl;
                throw std::runtime_error("NaN detected");
            }
    }
}

} // namespace Gaukuk
