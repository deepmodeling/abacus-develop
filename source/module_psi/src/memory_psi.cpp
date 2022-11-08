#include <iostream>
#include <complex>
#include <string.h>
#include "module_psi/include/types.h"
#include "module_psi/include/memory.h"

namespace psi{
namespace memory{

template <typename FPTYPE> 
struct resize_memory_op<FPTYPE, psi::DEVICE_CPU> {
  void operator()(const psi::DEVICE_CPU* dev, FPTYPE*& arr, const size_t size) {
    if (arr != nullptr) {
      free(arr);
    }
    arr = (FPTYPE*) malloc(sizeof(FPTYPE) * size);
  }
};

template <typename FPTYPE> 
struct set_memory_op<FPTYPE, psi::DEVICE_CPU> {
  void operator()(const psi::DEVICE_CPU* dev, FPTYPE* arr, const int var, const size_t size) {
    memset(arr, var, sizeof(FPTYPE) * size);
  }
};

template <typename FPTYPE> 
struct synchronize_memory_op<FPTYPE, psi::DEVICE_CPU, psi::DEVICE_CPU> {
  void operator()(const psi::DEVICE_CPU* dev_out, 
                  const psi::DEVICE_CPU* dev_in, 
                  FPTYPE* arr_out, 
                  const FPTYPE* arr_in, 
                  const size_t size) {
    memcpy(arr_out, arr_in, sizeof(FPTYPE) * size);
  }
};

template <typename FPTYPE>
struct delete_memory_op<FPTYPE, psi::DEVICE_CPU> {
  void operator()(const psi::DEVICE_CPU* dev, FPTYPE* arr) {
    free(arr);
  }
};

template struct resize_memory_op<int, psi::DEVICE_CPU>;
template struct resize_memory_op<double, psi::DEVICE_CPU>;
template struct resize_memory_op<std::complex<double>, psi::DEVICE_CPU>;

template struct set_memory_op<int, psi::DEVICE_CPU>;
template struct set_memory_op<double, psi::DEVICE_CPU>;
template struct set_memory_op<std::complex<double>, psi::DEVICE_CPU>;

template struct synchronize_memory_op<int, psi::DEVICE_CPU, psi::DEVICE_CPU>;
template struct synchronize_memory_op<double, psi::DEVICE_CPU, psi::DEVICE_CPU>;
template struct synchronize_memory_op<std::complex<double>, psi::DEVICE_CPU, psi::DEVICE_CPU>;

template struct delete_memory_op<int, psi::DEVICE_CPU>;
template struct delete_memory_op<double, psi::DEVICE_CPU>;
template struct delete_memory_op<std::complex<double>, psi::DEVICE_CPU>;

}
}