#include <iostream>
#include "template_array.hpp"

int main(){
    using namespace Gaukuk;
    int nVar = 5; 
    int nx = 2; 
    int ny = 3; 
    int nz = 4; 
    TArray<double> arr1(nVar, nz, ny, nx);
    // TArray<double> arr2(nVar, ny, nx);
    TArray<double> arr2; 
    arr2.NewArray(nVar, ny, nx); 
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

    return 0; 
}