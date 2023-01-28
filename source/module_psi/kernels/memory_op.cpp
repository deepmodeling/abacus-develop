#include <iostream>
#include <complex>
#include <string.h>
#include "module_psi/kernels/types.h"
#include "module_psi/kernels/memory_op.h"
#include "module_base/memory.h"

namespace psi{
namespace memory{

template <typename FPTYPE> 
struct resize_memory_op<FPTYPE, psi::DEVICE_CPU> {
  void operator()(const psi::DEVICE_CPU* dev, FPTYPE*& arr, const size_t size, const char* record_in) {
    if (arr != nullptr) {
      free(arr);
    }
    arr = (FPTYPE*) malloc(sizeof(FPTYPE) * size);
    std::string record_string;
    if(record_in != nullptr)
    {
      record_string = record_in;
    }
    else
    {
      record_string = "no_record";
    }
    
    if(record_string != "no_record" )
    {
      ModuleBase::Memory::record(record_string , sizeof(FPTYPE) * size);
    }
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

template <typename FPTYPE_out, typename FPTYPE_in>
struct cast_memory_op<FPTYPE_out, FPTYPE_in, psi::DEVICE_CPU, psi::DEVICE_CPU> {
    void operator()(const psi::DEVICE_CPU* dev_out,
                    const psi::DEVICE_CPU* dev_in,
                    FPTYPE_out* arr_out,
                    const FPTYPE_in* arr_in,
                    const size_t size) {
        for (int ii = 0; ii < size; ii++) {
            arr_out[ii] = static_cast<FPTYPE_out>(arr_in[ii]);
        }
    }
};

template <typename FPTYPE>
struct delete_memory_op<FPTYPE, psi::DEVICE_CPU> {
  void operator()(const psi::DEVICE_CPU* dev, FPTYPE* arr) {
    free(arr);
  }
};

template struct resize_memory_op<int, psi::DEVICE_CPU>;
template struct resize_memory_op<float, psi::DEVICE_CPU>;
template struct resize_memory_op<double, psi::DEVICE_CPU>;
template struct resize_memory_op<std::complex<float>, psi::DEVICE_CPU>;
template struct resize_memory_op<std::complex<double>, psi::DEVICE_CPU>;

template struct set_memory_op<int, psi::DEVICE_CPU>;
template struct set_memory_op<float, psi::DEVICE_CPU>;
template struct set_memory_op<double, psi::DEVICE_CPU>;
template struct set_memory_op<std::complex<float>, psi::DEVICE_CPU>;
template struct set_memory_op<std::complex<double>, psi::DEVICE_CPU>;

template struct synchronize_memory_op<int, psi::DEVICE_CPU, psi::DEVICE_CPU>;
template struct synchronize_memory_op<float, psi::DEVICE_CPU, psi::DEVICE_CPU>;
template struct synchronize_memory_op<double, psi::DEVICE_CPU, psi::DEVICE_CPU>;
template struct synchronize_memory_op<std::complex<float>, psi::DEVICE_CPU, psi::DEVICE_CPU>;
template struct synchronize_memory_op<std::complex<double>, psi::DEVICE_CPU, psi::DEVICE_CPU>;

template struct cast_memory_op<float, float, psi::DEVICE_CPU, psi::DEVICE_CPU>;
template struct cast_memory_op<double, double, psi::DEVICE_CPU, psi::DEVICE_CPU>;
template struct cast_memory_op<float, double, psi::DEVICE_CPU, psi::DEVICE_CPU>;
template struct cast_memory_op<double, float, psi::DEVICE_CPU, psi::DEVICE_CPU>;
template struct cast_memory_op<std::complex<float>, std::complex<float>, psi::DEVICE_CPU, psi::DEVICE_CPU>;
template struct cast_memory_op<std::complex<double>, std::complex<double>, psi::DEVICE_CPU, psi::DEVICE_CPU>;
template struct cast_memory_op<std::complex<float>, std::complex<double>, psi::DEVICE_CPU, psi::DEVICE_CPU>;
template struct cast_memory_op<std::complex<double>, std::complex<float>, psi::DEVICE_CPU, psi::DEVICE_CPU>;

template struct delete_memory_op<int, psi::DEVICE_CPU>;
template struct delete_memory_op<float, psi::DEVICE_CPU>;
template struct delete_memory_op<double, psi::DEVICE_CPU>;
template struct delete_memory_op<std::complex<float>, psi::DEVICE_CPU>;
template struct delete_memory_op<std::complex<double>, psi::DEVICE_CPU>;

#if !(defined(__CUDA) || defined(__ROCM))
template <typename FPTYPE>
struct resize_memory_op<FPTYPE, psi::DEVICE_GPU> {
    void operator()(const psi::DEVICE_GPU* dev, FPTYPE*& arr, const size_t size, const char* record_in = nullptr) {}
};

template <typename FPTYPE>
struct set_memory_op<FPTYPE, psi::DEVICE_GPU> {
    void operator()(const psi::DEVICE_GPU* dev, FPTYPE* arr, const int var, const size_t size) {}
};

template <typename FPTYPE>
struct synchronize_memory_op<FPTYPE, psi::DEVICE_GPU, psi::DEVICE_GPU> {
    void operator()(const psi::DEVICE_GPU* dev_out,
                    const psi::DEVICE_GPU* dev_in,
                    FPTYPE* arr_out,
                    const FPTYPE* arr_in,
                    const size_t size) {}
};

template <typename FPTYPE>
struct synchronize_memory_op<FPTYPE, psi::DEVICE_GPU, psi::DEVICE_CPU> {
    void operator()(const psi::DEVICE_GPU* dev_out,
                    const psi::DEVICE_CPU* dev_in,
                    FPTYPE* arr_out,
                    const FPTYPE* arr_in,
                    const size_t size) {}
};

template <typename FPTYPE>
struct synchronize_memory_op<FPTYPE, psi::DEVICE_CPU, psi::DEVICE_GPU> {
    void operator()(const psi::DEVICE_CPU* dev_out,
                    const psi::DEVICE_GPU* dev_in,
                    FPTYPE* arr_out,
                    const FPTYPE* arr_in,
                    const size_t size) {}
};

template <typename FPTYPE_out, typename FPTYPE_in>
struct cast_memory_op<FPTYPE_out, FPTYPE_in, psi::DEVICE_GPU, psi::DEVICE_GPU> {
    void operator()(const psi::DEVICE_GPU* dev_out,
                    const psi::DEVICE_GPU* dev_in,
                    FPTYPE_out* arr_out,
                    const FPTYPE_in* arr_in,
                    const size_t size) {}
};

template <typename FPTYPE_out, typename FPTYPE_in>
struct cast_memory_op<FPTYPE_out, FPTYPE_in, psi::DEVICE_GPU, psi::DEVICE_CPU> {
    void operator()(const psi::DEVICE_GPU* dev_out,
                    const psi::DEVICE_CPU* dev_in,
                    FPTYPE_out* arr_out,
                    const FPTYPE_in* arr_in,
                    const size_t size) {}
};

template <typename FPTYPE_out, typename FPTYPE_in>
struct cast_memory_op<FPTYPE_out, FPTYPE_in, psi::DEVICE_CPU, psi::DEVICE_GPU> {
    void operator()(const psi::DEVICE_CPU* dev_out,
                    const psi::DEVICE_GPU* dev_in,
                    FPTYPE_out* arr_out,
                    const FPTYPE_in* arr_in,
                    const size_t size) {}
};

template <typename FPTYPE>
struct delete_memory_op<FPTYPE, psi::DEVICE_GPU> {
    void operator()(const psi::DEVICE_GPU* dev, FPTYPE* arr) {}
};

template struct resize_memory_op<int, psi::DEVICE_GPU>;
template struct resize_memory_op<float, psi::DEVICE_GPU>;
template struct resize_memory_op<double, psi::DEVICE_GPU>;
template struct resize_memory_op<std::complex<float>, psi::DEVICE_GPU>;
template struct resize_memory_op<std::complex<double>, psi::DEVICE_GPU>;

template struct set_memory_op<int, psi::DEVICE_GPU>;
template struct set_memory_op<float, psi::DEVICE_GPU>;
template struct set_memory_op<double, psi::DEVICE_GPU>;
template struct set_memory_op<std::complex<float>, psi::DEVICE_GPU>;
template struct set_memory_op<std::complex<double>, psi::DEVICE_GPU>;

template struct synchronize_memory_op<int, psi::DEVICE_CPU, psi::DEVICE_GPU>;
template struct synchronize_memory_op<int, psi::DEVICE_GPU, psi::DEVICE_CPU>;
template struct synchronize_memory_op<int, psi::DEVICE_GPU, psi::DEVICE_GPU>;
template struct synchronize_memory_op<float, psi::DEVICE_CPU, psi::DEVICE_GPU>;
template struct synchronize_memory_op<float, psi::DEVICE_GPU, psi::DEVICE_CPU>;
template struct synchronize_memory_op<float, psi::DEVICE_GPU, psi::DEVICE_GPU>;
template struct synchronize_memory_op<double, psi::DEVICE_CPU, psi::DEVICE_GPU>;
template struct synchronize_memory_op<double, psi::DEVICE_GPU, psi::DEVICE_CPU>;
template struct synchronize_memory_op<double, psi::DEVICE_GPU, psi::DEVICE_GPU>;
template struct synchronize_memory_op<std::complex<float>, psi::DEVICE_CPU, psi::DEVICE_GPU>;
template struct synchronize_memory_op<std::complex<float>, psi::DEVICE_GPU, psi::DEVICE_CPU>;
template struct synchronize_memory_op<std::complex<float>, psi::DEVICE_GPU, psi::DEVICE_GPU>;
template struct synchronize_memory_op<std::complex<double>, psi::DEVICE_CPU, psi::DEVICE_GPU>;
template struct synchronize_memory_op<std::complex<double>, psi::DEVICE_GPU, psi::DEVICE_CPU>;
template struct synchronize_memory_op<std::complex<double>, psi::DEVICE_GPU, psi::DEVICE_GPU>;

template struct cast_memory_op<float, float, psi::DEVICE_GPU, psi::DEVICE_GPU>;
template struct cast_memory_op<double, double, psi::DEVICE_GPU, psi::DEVICE_GPU>;
template struct cast_memory_op<float, double, psi::DEVICE_GPU, psi::DEVICE_GPU>;
template struct cast_memory_op<double, float, psi::DEVICE_GPU, psi::DEVICE_GPU>;
template struct cast_memory_op<std::complex<float>, std::complex<float>, psi::DEVICE_GPU, psi::DEVICE_GPU>;
template struct cast_memory_op<std::complex<double>, std::complex<double>, psi::DEVICE_GPU, psi::DEVICE_GPU>;
template struct cast_memory_op<std::complex<float>, std::complex<double>, psi::DEVICE_GPU, psi::DEVICE_GPU>;
template struct cast_memory_op<std::complex<double>, std::complex<float>, psi::DEVICE_GPU, psi::DEVICE_GPU>;
template struct cast_memory_op<float, float, psi::DEVICE_GPU, psi::DEVICE_CPU>;
template struct cast_memory_op<double, double, psi::DEVICE_GPU, psi::DEVICE_CPU>;
template struct cast_memory_op<float, double, psi::DEVICE_GPU, psi::DEVICE_CPU>;
template struct cast_memory_op<double, float, psi::DEVICE_GPU, psi::DEVICE_CPU>;
template struct cast_memory_op<std::complex<float>, std::complex<float>, psi::DEVICE_GPU, psi::DEVICE_CPU>;
template struct cast_memory_op<std::complex<double>, std::complex<double>, psi::DEVICE_GPU, psi::DEVICE_CPU>;
template struct cast_memory_op<std::complex<float>, std::complex<double>, psi::DEVICE_GPU, psi::DEVICE_CPU>;
template struct cast_memory_op<std::complex<double>, std::complex<float>, psi::DEVICE_GPU, psi::DEVICE_CPU>;
template struct cast_memory_op<float, float, psi::DEVICE_CPU, psi::DEVICE_GPU>;
template struct cast_memory_op<double, double, psi::DEVICE_CPU, psi::DEVICE_GPU>;
template struct cast_memory_op<float, double, psi::DEVICE_CPU, psi::DEVICE_GPU>;
template struct cast_memory_op<double, float, psi::DEVICE_CPU, psi::DEVICE_GPU>;
template struct cast_memory_op<std::complex<float>, std::complex<float>, psi::DEVICE_CPU, psi::DEVICE_GPU>;
template struct cast_memory_op<std::complex<double>, std::complex<double>, psi::DEVICE_CPU, psi::DEVICE_GPU>;
template struct cast_memory_op<std::complex<float>, std::complex<double>, psi::DEVICE_CPU, psi::DEVICE_GPU>;
template struct cast_memory_op<std::complex<double>, std::complex<float>, psi::DEVICE_CPU, psi::DEVICE_GPU>;

template struct delete_memory_op<int, psi::DEVICE_GPU>;
template struct delete_memory_op<float, psi::DEVICE_GPU>;
template struct delete_memory_op<double, psi::DEVICE_GPU>;
template struct delete_memory_op<std::complex<float>, psi::DEVICE_GPU>;
template struct delete_memory_op<std::complex<double>, psi::DEVICE_GPU>;
#endif

}
}