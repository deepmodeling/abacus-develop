#include <iostream>
#include <complex>
#include <string.h>

#include <ATen/kernels/memory.h>

namespace container {
namespace kernels {

template <typename T>
struct resize_memory<T, DEVICE_CPU> {
    void operator()(T*& arr, const size_t& size, const char* /*record_in*/) {
        if (arr != nullptr) {
            free(arr);
        }
        arr = (T*) malloc(sizeof(T) * size);
    }
};

template <typename T>
struct set_memory<T, DEVICE_CPU> {
    void operator()(T* arr, const T& var, const size_t& size) {
        for (size_t ii = 0; ii < size; ii++) {
            arr[ii] = var;
        }
    }
};

template <typename T>
struct synchronize_memory<T, DEVICE_CPU, DEVICE_CPU> {
    void operator()(
        T* arr_out,
        const T* arr_in,
        const size_t& size)
    {
      memcpy(arr_out, arr_in, sizeof(T) * size);
    }
};

template <typename T_out, typename T_in>
struct cast_memory<T_out, T_in, DEVICE_CPU, DEVICE_CPU> {
    void operator()(
        T_out* arr_out,
        const T_in* arr_in,
        const size_t& size)
        {
            for (int ii = 0; ii < size; ii++) {
                arr_out[ii] = static_cast<T_out>(arr_in[ii]);
            }
        }
};

template <typename T>
struct delete_memory<T, DEVICE_CPU> {
    void operator()(T* arr) {
        free(arr);
    }
};

template struct resize_memory<int, DEVICE_CPU>;
template struct resize_memory<int64_t, DEVICE_CPU>;
template struct resize_memory<float, DEVICE_CPU>;
template struct resize_memory<double, DEVICE_CPU>;
template struct resize_memory<std::complex<float>, DEVICE_CPU>;
template struct resize_memory<std::complex<double>, DEVICE_CPU>;

template struct set_memory<int, DEVICE_CPU>;
template struct set_memory<int64_t, DEVICE_CPU>;
template struct set_memory<float, DEVICE_CPU>;
template struct set_memory<double, DEVICE_CPU>;
template struct set_memory<std::complex<float>, DEVICE_CPU>;
template struct set_memory<std::complex<double>, DEVICE_CPU>;

template struct synchronize_memory<int, DEVICE_CPU, DEVICE_CPU>;
template struct synchronize_memory<int64_t, DEVICE_CPU, DEVICE_CPU>;
template struct synchronize_memory<float, DEVICE_CPU, DEVICE_CPU>;
template struct synchronize_memory<double, DEVICE_CPU, DEVICE_CPU>;
template struct synchronize_memory<std::complex<float>, DEVICE_CPU, DEVICE_CPU>;
template struct synchronize_memory<std::complex<double>, DEVICE_CPU, DEVICE_CPU>;

template struct cast_memory<float, float, DEVICE_CPU, DEVICE_CPU>;
template struct cast_memory<double, double, DEVICE_CPU, DEVICE_CPU>;
template struct cast_memory<float, double, DEVICE_CPU, DEVICE_CPU>;
template struct cast_memory<double, float, DEVICE_CPU, DEVICE_CPU>;
template struct cast_memory<std::complex<float>, std::complex<float>, DEVICE_CPU, DEVICE_CPU>;
template struct cast_memory<std::complex<double>, std::complex<double>, DEVICE_CPU, DEVICE_CPU>;
template struct cast_memory<std::complex<float>, std::complex<double>, DEVICE_CPU, DEVICE_CPU>;
template struct cast_memory<std::complex<double>, std::complex<float>, DEVICE_CPU, DEVICE_CPU>;

template struct delete_memory<int, DEVICE_CPU>;
template struct delete_memory<int64_t, DEVICE_CPU>;
template struct delete_memory<float, DEVICE_CPU>;
template struct delete_memory<double, DEVICE_CPU>;
template struct delete_memory<std::complex<float>, DEVICE_CPU>;
template struct delete_memory<std::complex<double>, DEVICE_CPU>;

#if !(defined(__CUDA) || defined(__ROCM))
template <typename T>
struct resize_memory<T, DEVICE_GPU> {
    void operator()(T*& arr, const size_t& size, const char* record_in = nullptr) {}
};

template <typename T>
struct set_memory<T, DEVICE_GPU> {
    void operator()(T* arr, const int var, const size_t& size) {}
};

template <typename T>
struct synchronize_memory<T, DEVICE_GPU, DEVICE_GPU> {
    void operator()(T* arr_out, const T* arr_in, const size_t& size) {}
};

template <typename T>
struct synchronize_memory<T, DEVICE_GPU, DEVICE_CPU> {
    void operator()(T* arr_out, const T* arr_in, const size_t& size) {}
};

template <typename T>
struct synchronize_memory<T, DEVICE_CPU, DEVICE_GPU> {
    void operator()(T* arr_out, const T* arr_in, const size_t& size) {}
};

template <typename T_out, typename T_in>
struct cast_memory<T_out, T_in, DEVICE_GPU, DEVICE_GPU> {
    void operator()(T_out* arr_out, const T_in* arr_in, const size_t& size) {}
};

template <typename T_out, typename T_in>
struct cast_memory<T_out, T_in, DEVICE_GPU, DEVICE_CPU> {
    void operator()(T_out* arr_out, const T_in* arr_in, const size_t& size) {}
};

template <typename T_out, typename T_in>
struct cast_memory<T_out, T_in, DEVICE_CPU, DEVICE_GPU> {
    void operator()(T_out* arr_out, const T_in* arr_in, const size_t& size) {}
};

template <typename T>
struct delete_memory<T, DEVICE_GPU> {
    void operator()(T* arr) {}
};

template struct resize_memory<int, DEVICE_GPU>;
template struct resize_memory<int64_t, DEVICE_GPU>;
template struct resize_memory<float, DEVICE_GPU>;
template struct resize_memory<double, DEVICE_GPU>;
template struct resize_memory<std::complex<float>, DEVICE_GPU>;
template struct resize_memory<std::complex<double>, DEVICE_GPU>;

template struct set_memory<int, DEVICE_GPU>;
template struct set_memory<int64_t, DEVICE_GPU>;
template struct set_memory<float, DEVICE_GPU>;
template struct set_memory<double, DEVICE_GPU>;
template struct set_memory<std::complex<float>, DEVICE_GPU>;
template struct set_memory<std::complex<double>, DEVICE_GPU>;

template struct synchronize_memory<int, DEVICE_CPU, DEVICE_GPU>;
template struct synchronize_memory<int, DEVICE_GPU, DEVICE_CPU>;
template struct synchronize_memory<int, DEVICE_GPU, DEVICE_GPU>;
template struct synchronize_memory<int64_t, DEVICE_CPU, DEVICE_GPU>;
template struct synchronize_memory<int64_t, DEVICE_GPU, DEVICE_CPU>;
template struct synchronize_memory<int64_t, DEVICE_GPU, DEVICE_GPU>;
template struct synchronize_memory<float, DEVICE_CPU, DEVICE_GPU>;
template struct synchronize_memory<float, DEVICE_GPU, DEVICE_CPU>;
template struct synchronize_memory<float, DEVICE_GPU, DEVICE_GPU>;
template struct synchronize_memory<double, DEVICE_CPU, DEVICE_GPU>;
template struct synchronize_memory<double, DEVICE_GPU, DEVICE_CPU>;
template struct synchronize_memory<double, DEVICE_GPU, DEVICE_GPU>;
template struct synchronize_memory<std::complex<float>, DEVICE_CPU, DEVICE_GPU>;
template struct synchronize_memory<std::complex<float>, DEVICE_GPU, DEVICE_CPU>;
template struct synchronize_memory<std::complex<float>, DEVICE_GPU, DEVICE_GPU>;
template struct synchronize_memory<std::complex<double>, DEVICE_CPU, DEVICE_GPU>;
template struct synchronize_memory<std::complex<double>, DEVICE_GPU, DEVICE_CPU>;
template struct synchronize_memory<std::complex<double>, DEVICE_GPU, DEVICE_GPU>;

template struct cast_memory<float, float, container::DEVICE_GPU, container::DEVICE_GPU>;
template struct cast_memory<double, double, container::DEVICE_GPU, container::DEVICE_GPU>;
template struct cast_memory<float, double, container::DEVICE_GPU, container::DEVICE_GPU>;
template struct cast_memory<double, float, container::DEVICE_GPU, container::DEVICE_GPU>;
template struct cast_memory<std::complex<float>, std::complex<float>, container::DEVICE_GPU, container::DEVICE_GPU>;
template struct cast_memory<std::complex<double>, std::complex<double>, container::DEVICE_GPU, container::DEVICE_GPU>;
template struct cast_memory<std::complex<float>, std::complex<double>, container::DEVICE_GPU, container::DEVICE_GPU>;
template struct cast_memory<std::complex<double>, std::complex<float>, container::DEVICE_GPU, container::DEVICE_GPU>;
template struct cast_memory<float, float, container::DEVICE_GPU, container::DEVICE_CPU>;
template struct cast_memory<double, double, container::DEVICE_GPU, container::DEVICE_CPU>;
template struct cast_memory<float, double, container::DEVICE_GPU, container::DEVICE_CPU>;
template struct cast_memory<double, float, container::DEVICE_GPU, container::DEVICE_CPU>;
template struct cast_memory<std::complex<float>, std::complex<float>, container::DEVICE_GPU, container::DEVICE_CPU>;
template struct cast_memory<std::complex<double>, std::complex<double>, container::DEVICE_GPU, container::DEVICE_CPU>;
template struct cast_memory<std::complex<float>, std::complex<double>, container::DEVICE_GPU, container::DEVICE_CPU>;
template struct cast_memory<std::complex<double>, std::complex<float>, container::DEVICE_GPU, container::DEVICE_CPU>;
template struct cast_memory<float, float, container::DEVICE_CPU, container::DEVICE_GPU>;
template struct cast_memory<double, double, container::DEVICE_CPU, container::DEVICE_GPU>;
template struct cast_memory<float, double, container::DEVICE_CPU, container::DEVICE_GPU>;
template struct cast_memory<double, float, container::DEVICE_CPU, container::DEVICE_GPU>;
template struct cast_memory<std::complex<float>, std::complex<float>, container::DEVICE_CPU, container::DEVICE_GPU>;
template struct cast_memory<std::complex<double>, std::complex<double>, container::DEVICE_CPU, container::DEVICE_GPU>;
template struct cast_memory<std::complex<float>, std::complex<double>, container::DEVICE_CPU, container::DEVICE_GPU>;
template struct cast_memory<std::complex<double>, std::complex<float>, container::DEVICE_CPU, container::DEVICE_GPU>;

template struct delete_memory<int, DEVICE_GPU>;
template struct delete_memory<int64_t, DEVICE_GPU>;
template struct delete_memory<float, DEVICE_GPU>;
template struct delete_memory<double, DEVICE_GPU>;
template struct delete_memory<std::complex<float>, DEVICE_GPU>;
template struct delete_memory<std::complex<double>, DEVICE_GPU>;
#endif

} // namespace kernels
} // namespace container