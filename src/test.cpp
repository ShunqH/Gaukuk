#include <iostream>
#include "TArray.hpp"

int main(){
    using namespace Gaukuk;
    int nVar = 5; 
    int nx = 2; 
    int ny = 3; 
    int nz = 4; 
    TArray<double> arr1(nVar, nz, ny, nx);
    for (size_t i = 0; i < arr1.Size(); i++){
        arr1(i) = i; 
    } 
    std::cout<<arr1(2)<<std::endl; 
    std::cout<<arr1(2, 2)<<std::endl; 
    std::cout<<arr1(2, 3, 2)<<std::endl; 
    std::cout<<arr1(2, 3, 2, 1)<<std::endl; 
    return 0; 
}