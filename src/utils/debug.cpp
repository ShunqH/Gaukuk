// C++ Header 
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <string>           // std::string 
#include <fstream>          // std::ofstream
#include <sstream>          // std::ostringstream; .str(); 
#include <iomanip>          // std::setw(); std::setfill()

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

void WriteTarray(const TArray<Real>& arr, const char* name, const int ID){
    std::ostringstream filename_stream;
    filename_stream << name << "_"
                    << std::setw(5) << std::setfill('0') << ID;

    std::string filename = filename_stream.str();

    std::ofstream outFile(filename, std::ios::binary);
    if (!outFile) {
        throw std::runtime_error("Failed to open " + filename);
    }

    int elem_size = sizeof(Real);
    int size = arr.GetSize();

    outFile.write(reinterpret_cast<const char*>(&elem_size), sizeof(int));
    outFile.write(reinterpret_cast<const char*>(&size), sizeof(int));
    outFile.write(reinterpret_cast<const char*>(arr.data()),
                  arr.GetSizeInBytes());
}

} // namespace Gaukuk
