#include <iostream>
#include <complex>
#include <string.h>
#include <ATen/kernels/memory_op.h>

namespace container {
namespace op {

template <typename T>
struct resize_memory_op<T, DEVICE_CPU> {
  void operator()(T*& arr, const size_t size, const char* /*record_in*/) {
    if (arr != nullptr) {
      free(arr);
    }
    arr = (T*) malloc(sizeof(T) * size);
  }
};

template <typename T>
struct set_memory_op<T, DEVICE_CPU> {
  void operator()(T* arr, const T& var, const size_t& size) {
    for (size_t ii = 0; ii < size; ii++) {
      arr[ii] = var;
    }
  }
};

template <typename T>
struct synchronize_memory_op<T, DEVICE_CPU, DEVICE_CPU> {
  void operator()(T* arr_out,
                  const T* arr_in,
                  const size_t size) {
    memcpy(arr_out, arr_in, sizeof(T) * size);
  }
};

template <typename FPTYPE_out, typename FPTYPE_in>
struct cast_memory_op<FPTYPE_out, FPTYPE_in, DEVICE_CPU, DEVICE_CPU> {
    void operator()(FPTYPE_out* arr_out,
                    const FPTYPE_in* arr_in,
                    const size_t size) {
        for (int ii = 0; ii < size; ii++) {
            arr_out[ii] = static_cast<FPTYPE_out>(arr_in[ii]);
        }
    }
};

template <typename T>
struct delete_memory_op<T, DEVICE_CPU> {
  void operator()(T* arr) {
    free(arr);
  }
};

template struct resize_memory_op<int, DEVICE_CPU>;
template struct resize_memory_op<int64_t, DEVICE_CPU>;
template struct resize_memory_op<float, DEVICE_CPU>;
template struct resize_memory_op<double, DEVICE_CPU>;
template struct resize_memory_op<std::complex<float>, DEVICE_CPU>;
template struct resize_memory_op<std::complex<double>, DEVICE_CPU>;

template struct set_memory_op<int, DEVICE_CPU>;
template struct set_memory_op<int64_t, DEVICE_CPU>;
template struct set_memory_op<float, DEVICE_CPU>;
template struct set_memory_op<double, DEVICE_CPU>;
template struct set_memory_op<std::complex<float>, DEVICE_CPU>;
template struct set_memory_op<std::complex<double>, DEVICE_CPU>;

template struct synchronize_memory_op<int, DEVICE_CPU, DEVICE_CPU>;
template struct synchronize_memory_op<int64_t, DEVICE_CPU, DEVICE_CPU>;
template struct synchronize_memory_op<float, DEVICE_CPU, DEVICE_CPU>;
template struct synchronize_memory_op<double, DEVICE_CPU, DEVICE_CPU>;
template struct synchronize_memory_op<std::complex<float>, DEVICE_CPU, DEVICE_CPU>;
template struct synchronize_memory_op<std::complex<double>, DEVICE_CPU, DEVICE_CPU>;

template struct cast_memory_op<float, float, DEVICE_CPU, DEVICE_CPU>;
template struct cast_memory_op<double, double, DEVICE_CPU, DEVICE_CPU>;
template struct cast_memory_op<float, double, DEVICE_CPU, DEVICE_CPU>;
template struct cast_memory_op<double, float, DEVICE_CPU, DEVICE_CPU>;
template struct cast_memory_op<std::complex<float>, std::complex<float>, DEVICE_CPU, DEVICE_CPU>;
template struct cast_memory_op<std::complex<double>, std::complex<double>, DEVICE_CPU, DEVICE_CPU>;
template struct cast_memory_op<std::complex<float>, std::complex<double>, DEVICE_CPU, DEVICE_CPU>;
template struct cast_memory_op<std::complex<double>, std::complex<float>, DEVICE_CPU, DEVICE_CPU>;

template struct delete_memory_op<int, DEVICE_CPU>;
template struct delete_memory_op<int64_t, DEVICE_CPU>;
template struct delete_memory_op<float, DEVICE_CPU>;
template struct delete_memory_op<double, DEVICE_CPU>;
template struct delete_memory_op<std::complex<float>, DEVICE_CPU>;
template struct delete_memory_op<std::complex<double>, DEVICE_CPU>;

#if !(defined(__CUDA) || defined(__ROCM))
template <typename T>
struct resize_memory_op<T, DEVICE_GPU> {
    void operator()(T*& arr, const size_t size, const char* record_in = nullptr) {}
};

template <typename T>
struct set_memory_op<T, DEVICE_GPU> {
    void operator()(T* arr, const int var, const size_t size) {}
};

template <typename T>
struct synchronize_memory_op<T, DEVICE_GPU, DEVICE_GPU> {
    void operator()(T* arr_out, const T* arr_in, const size_t size) {}
};

template <typename T>
struct synchronize_memory_op<T, DEVICE_GPU, DEVICE_CPU> {
    void operator()(T* arr_out, const T* arr_in, const size_t size) {}
};

template <typename T>
struct synchronize_memory_op<T, DEVICE_CPU, DEVICE_GPU> {
    void operator()(T* arr_out, const T* arr_in, const size_t size) {}
};

template <typename FPTYPE_out, typename FPTYPE_in>
struct cast_memory_op<FPTYPE_out, FPTYPE_in, DEVICE_GPU, DEVICE_GPU> {
    void operator()(FPTYPE_out* arr_out, const FPTYPE_in* arr_in, const size_t size) {}
};

template <typename FPTYPE_out, typename FPTYPE_in>
struct cast_memory_op<FPTYPE_out, FPTYPE_in, DEVICE_GPU, DEVICE_CPU> {
    void operator()(FPTYPE_out* arr_out, const FPTYPE_in* arr_in, const size_t size) {}
};

template <typename FPTYPE_out, typename FPTYPE_in>
struct cast_memory_op<FPTYPE_out, FPTYPE_in, DEVICE_CPU, DEVICE_GPU> {
    void operator()(FPTYPE_out* arr_out, const FPTYPE_in* arr_in, const size_t size) {}
};

template <typename T>
struct delete_memory_op<T, DEVICE_GPU> {
    void operator()(T* arr) {}
};

template struct resize_memory_op<int, DEVICE_GPU>;
template struct resize_memory_op<int64_t, DEVICE_GPU>;
template struct resize_memory_op<float, DEVICE_GPU>;
template struct resize_memory_op<double, DEVICE_GPU>;
template struct resize_memory_op<std::complex<float>, DEVICE_GPU>;
template struct resize_memory_op<std::complex<double>, DEVICE_GPU>;

template struct set_memory_op<int, DEVICE_GPU>;
template struct set_memory_op<int64_t, DEVICE_GPU>;
template struct set_memory_op<float, DEVICE_GPU>;
template struct set_memory_op<double, DEVICE_GPU>;
template struct set_memory_op<std::complex<float>, DEVICE_GPU>;
template struct set_memory_op<std::complex<double>, DEVICE_GPU>;

template struct synchronize_memory_op<int, DEVICE_CPU, DEVICE_GPU>;
template struct synchronize_memory_op<int, DEVICE_GPU, DEVICE_CPU>;
template struct synchronize_memory_op<int, DEVICE_GPU, DEVICE_GPU>;
template struct synchronize_memory_op<int64_t, DEVICE_CPU, DEVICE_GPU>;
template struct synchronize_memory_op<int64_t, DEVICE_GPU, DEVICE_CPU>;
template struct synchronize_memory_op<int64_t, DEVICE_GPU, DEVICE_GPU>;
template struct synchronize_memory_op<float, DEVICE_CPU, DEVICE_GPU>;
template struct synchronize_memory_op<float, DEVICE_GPU, DEVICE_CPU>;
template struct synchronize_memory_op<float, DEVICE_GPU, DEVICE_GPU>;
template struct synchronize_memory_op<double, DEVICE_CPU, DEVICE_GPU>;
template struct synchronize_memory_op<double, DEVICE_GPU, DEVICE_CPU>;
template struct synchronize_memory_op<double, DEVICE_GPU, DEVICE_GPU>;
template struct synchronize_memory_op<std::complex<float>, DEVICE_CPU, DEVICE_GPU>;
template struct synchronize_memory_op<std::complex<float>, DEVICE_GPU, DEVICE_CPU>;
template struct synchronize_memory_op<std::complex<float>, DEVICE_GPU, DEVICE_GPU>;
template struct synchronize_memory_op<std::complex<double>, DEVICE_CPU, DEVICE_GPU>;
template struct synchronize_memory_op<std::complex<double>, DEVICE_GPU, DEVICE_CPU>;
template struct synchronize_memory_op<std::complex<double>, DEVICE_GPU, DEVICE_GPU>;

template struct cast_memory_op<float, float, DEVICE_GPU, DEVICE_GPU>;
template struct cast_memory_op<double, double, DEVICE_GPU, DEVICE_GPU>;
template struct cast_memory_op<float, double, DEVICE_GPU, DEVICE_GPU>;
template struct cast_memory_op<double, float, DEVICE_GPU, DEVICE_GPU>;
template struct cast_memory_op<std::complex<float>, std::complex<float>, DEVICE_GPU, DEVICE_GPU>;
template struct cast_memory_op<std::complex<double>, std::complex<double>, DEVICE_GPU, DEVICE_GPU>;
template struct cast_memory_op<std::complex<float>, std::complex<double>, DEVICE_GPU, DEVICE_GPU>;
template struct cast_memory_op<std::complex<double>, std::complex<float>, DEVICE_GPU, DEVICE_GPU>;
template struct cast_memory_op<float, float, DEVICE_GPU, DEVICE_CPU>;
template struct cast_memory_op<double, double, DEVICE_GPU, DEVICE_CPU>;
template struct cast_memory_op<float, double, DEVICE_GPU, DEVICE_CPU>;
template struct cast_memory_op<double, float, DEVICE_GPU, DEVICE_CPU>;
template struct cast_memory_op<std::complex<float>, std::complex<float>, DEVICE_GPU, DEVICE_CPU>;
template struct cast_memory_op<std::complex<double>, std::complex<double>, DEVICE_GPU, DEVICE_CPU>;
template struct cast_memory_op<std::complex<float>, std::complex<double>, DEVICE_GPU, DEVICE_CPU>;
template struct cast_memory_op<std::complex<double>, std::complex<float>, DEVICE_GPU, DEVICE_CPU>;
template struct cast_memory_op<float, float, DEVICE_CPU, DEVICE_GPU>;
template struct cast_memory_op<double, double, DEVICE_CPU, DEVICE_GPU>;
template struct cast_memory_op<float, double, DEVICE_CPU, DEVICE_GPU>;
template struct cast_memory_op<double, float, DEVICE_CPU, DEVICE_GPU>;
template struct cast_memory_op<std::complex<float>, std::complex<float>, DEVICE_CPU, DEVICE_GPU>;
template struct cast_memory_op<std::complex<double>, std::complex<double>, DEVICE_CPU, DEVICE_GPU>;
template struct cast_memory_op<std::complex<float>, std::complex<double>, DEVICE_CPU, DEVICE_GPU>;
template struct cast_memory_op<std::complex<double>, std::complex<float>, DEVICE_CPU, DEVICE_GPU>;

template struct delete_memory_op<int, DEVICE_GPU>;
template struct delete_memory_op<int64_t, DEVICE_GPU>;
template struct delete_memory_op<float, DEVICE_GPU>;
template struct delete_memory_op<double, DEVICE_GPU>;
template struct delete_memory_op<std::complex<float>, DEVICE_GPU>;
template struct delete_memory_op<std::complex<double>, DEVICE_GPU>;
#endif

}
}