// CPP header
#include <algorithm> 
#include <iostream>     // std::cout; std::endl; std::cerr

// Gaukuk dependence
#include "gaukuk.hpp" 
#include "sim.hpp" 

int main(int argc, char* argv[]){
    using namespace Gaukuk; 

    if (argc < 3 || std::string(argv[1]) != "-i") {
        std::cerr << "Usage: " << argv[0] << " -i input.in" << std::endl;
        return 1;
    }
    Config::getInstance().loadFromFile(argv[2]);

    Real tmax = Config::getInstance().get("tmax"); 
    Real dtoutput = Config::getInstance().get("dtoutput"); 
    int stepmax = Config::getInstance().get("stepmax"); 

    Sim sim; 
    sim.Setup(); 
    sim.boundary.UpdateBD(sim.cons, sim.grid);

    int step = 0; 
    Real tnow = sim.GetTime(); 

    DataType outputType = DataType::Cons; 
    int oTypeFromConfig = static_cast<int>(Config::getInstance().get("dataType")); 
    if (oTypeFromConfig == 0){
        outputType = DataType::Cons; 
    }else if (oTypeFromConfig == 1){
        outputType = DataType::Prim; 
    }

    sim.WriteData(step, outputType); 
    while (tnow<tmax && (stepmax<0 || step<stepmax)){
        dtoutput = std::min(dtoutput, tmax-tnow); 
        sim.Advance(dtoutput); 
        tnow = sim.GetTime();
        step ++;  
        sim.WriteData(step, outputType); 
        std::cout<<"output "<< step << std::endl; 
        // WriteTarray(sim.flx1, "flx1", step); 
    }
    // sim.WriteData(step, outputType); 
    std::cout<<"finished! "<< std::endl; 
}