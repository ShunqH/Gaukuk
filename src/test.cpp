// g++ -fopenmp -O3 -march=native -fopt-info-vec test.cpp adiabatic.cpp sim.cpp

#include <iostream> 
#include <iomanip> 
#include <chrono> 
#include <ctime>

#include "gaukuk.hpp"
#include "sim.hpp"
#include "template_array.hpp"

#include "eos/eos.hpp"
#include "utils/utils.hpp"

int main(int argc, char* argv[]){
    using namespace Gaukuk;

    if (argc < 3 || std::string(argv[1]) != "-i") {
        std::cerr << "Usage: " << argv[0] << " -i input.in" << std::endl;
        return 1;
    }
    Config::getInstance().loadFromFile(argv[2]);
    
    Sim sim; 
    Real den = Config::getInstance().get("rho"); 
    Real vlx = Config::getInstance().get("vx"); 
    Real vly = Config::getInstance().get("vy"); 
    Real vlz = Config::getInstance().get("vz"); 
    Real pre = Config::getInstance().get("pressure"); 
    Real gamma = Config::getInstance().get("gamma"); 

    Real engInt = pre/(gamma - 1); 
    Real engKin = 0.5 * den * ( vlx*vlx + vly*vly + vlz*vlz ); 
    Real eng = engInt + engKin; 
#pragma omp parallel for
    for (int k = 0; k<sim.grid.lenz; k++){
        for (int j = 0; j<sim.grid.leny; j++){
            #pragma omp simd 
            for (int i = 0; i<sim.grid.lenx; i++){
                Real& consDen = sim.cons(DEN, k, j, i); 
                Real& consMTX = sim.cons(MTX, k, j, i); 
                Real& consMTY = sim.cons(MTY, k, j, i); 
                Real& consMTZ = sim.cons(MTZ, k, j, i); 
                Real& consEng = sim.cons(ENG, k, j, i); 

                consDen = den; 
                consMTX = den * vlx; 
                consMTY = den * vly; 
                consMTZ = den * vlz; 
                consEng = eng; 
            }
        }
    }

    
    int steps = 100; 
    auto time1 = std::chrono::high_resolution_clock::now();
    auto time2 = std::chrono::high_resolution_clock::now();
    clock_t cputime1 = clock();
    clock_t cputime2 = clock();
    for (int n=0; n<steps; n++){
        sim.eos.ConsToPrim(sim.cons, sim.prim, sim.grid);
        sim.eos.PrimToCons(sim.prim, sim.cons, sim.grid); 
    }

    // std::cout<<sim.cons.GetN1()<<std::endl; 
    // std::cout<<sim.cons.GetN2()<<std::endl; 
    // std::cout<<sim.cons.GetN3()<<std::endl; 
    // std::cout<<sim.cons.GetN4()<<std::endl; 
    // std::cout<<sim.cons.GetSize()/5<<std::endl; 
    sim.eos.ConsToPrim(sim.cons, sim.prim, sim.grid); 
    // std::cout<<sim.cons(DEN, 5, 5, 5)<<std::endl; 
    // std::cout<<sim.cons(MTX, 5, 5, 5)<<std::endl; 
    // std::cout<<sim.cons(ENG, 5, 5, 5)<<std::endl; 
    // std::cout<<eng<<std::endl; 
    std::cout<<pre<<std::endl; 
    std::cout<<(gamma-1)*(eng - engKin)<<std::endl; 
    std::cout<<sim.prim(PRE, 5, 5, 5)<<std::endl; 

    std::cout<<eng<<std::endl; 
    std::cout<<sim.cons(ENG, 5, 5, 5)<<std::endl; 
    sim.eos.PrimToCons(sim.prim, sim.cons, sim.grid); 
    std::cout<<sim.cons(ENG, 5, 5, 5)<<std::endl; 

    std::cout<<gamma*pre/den<<std::endl; 
    Real cs = sim.eos.SoundSpeed(sim.prim(DEN, 5, 5, 5), sim.prim(PRE, 5, 5, 5)); 
    std::cout<<cs*cs<<std::endl; 

    time2 = std::chrono::high_resolution_clock::now();
    cputime2 = clock(); 
    Real walltime = std::chrono::duration<double>(time2 - time1).count();
    Real cputime = double(cputime2 - cputime1) / CLOCKS_PER_SEC;

    std::cout << std::right
        << "Finish! " 
        << "\n    walltime = " << std::setw(8) << std::fixed << std::setprecision(3) << walltime
        << " s, \n    cputime = " << std::setw(8) << std::fixed << std::setprecision(3) << cputime << " s"
        << "\n-----------------------------------------------------------"
        << std::endl;
    return 0; 
}



/*
    int NVar = 5; 
    int nx = 4; 
    int ny = 4; 
    int nz = 4; 
    TArray<double> arr1(NVar, nz, ny, nx);
    // TArray<double> arr2(NVar, ny, nx);
    TArray<double> arr2; 
    arr2.NewArray(NVar, ny, nx); 
    for (size_t i = 0; i < arr1.GetSize(); i++){
        arr1(i) = i; 
    } 
    for (size_t i = 0; i < arr2.GetSize(); i++){
        arr2(i) = i + arr1.GetSize(); 
    } 
    // std::cout<<arr1(2, 3, 2, 1)<<std::endl; 
    std::cout<<arr1.GetSize()<<std::endl; 
    std::cout<<arr1.GetSizeInBytes()<<std::endl; 
    std::cout<<arr2.GetSize()<<std::endl; 
    std::cout<<arr2.GetSizeInBytes()<<std::endl; 
    std::cout<<arr1(20)<<std::endl; 
    std::cout<<arr2(20)<<std::endl; 
    std::cout<<"swap array"<<std::endl; 
    arr1.Swap(arr2); 
    std::cout<<arr1.GetSize()<<std::endl; 
    std::cout<<arr1.GetSizeInBytes()<<std::endl; 
    std::cout<<arr2.GetSize()<<std::endl; 
    std::cout<<arr2.GetSizeInBytes()<<std::endl; 
    std::cout<<arr1(20)<<std::endl; 
    std::cout<<arr2(20)<<std::endl; 
*/
    