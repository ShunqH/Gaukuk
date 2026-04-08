// C++ headers
#include <string>           // std::string 
#include <fstream>          // std::ofstream
#include <sstream>          // std::ostringstream; .str(); 
#include <iomanip>          // std::setw(); std::setfill()

// Gaukuk dependence
#include "../gaukuk.hpp" 
#include "../sim.hpp" 

namespace Gaukuk
{
    
void Sim::WriteCons(const int outputID){
    std::ostringstream filename_stream;
    filename_stream << "output_" << std::setw(5) << std::setfill('0') << outputID;  
    std::string filename = filename_stream.str(); 

    std::ofstream outFile(filename, std::ios::binary); 
    if (!outFile) {
        throw std::runtime_error("Failed to open " + filename);
    }

    // write frame information
    outFile.write(reinterpret_cast<const char*>(&t), sizeof(Real));

    // write number of quantity and size of x, y, z, and ghost
    outFile.write(reinterpret_cast<const char*>(&NVar), sizeof(int));
    outFile.write(reinterpret_cast<const char*>(&grid.nx), sizeof(int));
    outFile.write(reinterpret_cast<const char*>(&grid.ny), sizeof(int));
    outFile.write(reinterpret_cast<const char*>(&grid.nz), sizeof(int));
    outFile.write(reinterpret_cast<const char*>(&grid.nGhost), sizeof(int));

    // write mesh: x, y, z 
    outFile.write(reinterpret_cast<const char*>(domain.xc.data()), domain.xc.GetSizeInBytes());
    outFile.write(reinterpret_cast<const char*>(domain.yc.data()), domain.yc.GetSizeInBytes());
    outFile.write(reinterpret_cast<const char*>(domain.zc.data()), domain.zc.GetSizeInBytes());

    // write data
    outFile.write(reinterpret_cast<const char*>(cons.data()), cons.GetSizeInBytes()); 
}

} // namespace Gaukuk
