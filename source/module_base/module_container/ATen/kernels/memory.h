#ifndef ATEN_KERNELS_MEMORY_H_
#define ATEN_KERNELS_MEMORY_H_

#include <vector>
#include <complex>

#include <base/macros/macros.h>
#include <ATen/core/tensor_types.h>

namespace container {
namespace kernels {

/**
 * @brief A functor to resize memory allocation.
 * @tparam T Floating-point type of the allocated memory.
 * @tparam Device Device type where the memory will be allocated.
 */
template <typename T, typename Device>
struct resize_memory {
    /**
     * @brief Resize memory allocation.
     *
     * @param dev Device where the memory will be allocated.
     * @param arr Pointer to the allocated memory.
     * @param size New size of the allocated memory.
     * @param record_in Optional message to record the resize operation.
     */
    void operator()(T*& arr, const size_t& size, const char* record_in = nullptr);
};

/**
 * @brief A functor to set memory to a constant value.
 * @tparam T Floating-point type of the memory.
 * @tparam Device Device type where the memory is allocated.
 */
template <typename T, typename Device>
struct set_memory {
    /**
     * @brief Set memory to a constant value.
     *
     * @param arr Pointer to the memory.
     * @param var Constant value to set.
     * @param size Size of the memory to set.
     */
    void operator()(T* arr, const T& var, const size_t& size);
};

/**
 * @brief Synchronizes memory between devices.
 *
 * This class synchronizes memory between two different devices.
 *
 * @tparam T The type of data in the arrays.
 * @tparam Device_out The output device.
 * @tparam Device_in The input device.
 */
template <typename T, typename Device_out, typename Device_in>
struct synchronize_memory {
    /**
     * @brief Synchronizes memory between devices.
     *
     * This method synchronizes memory between two different devices.
     *
     * @param dev_out The output device.
     * @param dev_in The input device.
     * @param arr_out The output array.
     * @param arr_in The input array.
     * @param size The size of the array.
     */
    void operator()(
        T* arr_out,
        const T* arr_in,
        const size_t& size);
};

template <typename T, typename Device_out, typename Device_in>
struct synchronize_memory_stride {
    void operator()(
        T* arr_out,
        const T* arr_in,
        const std::vector<int64_t>& out_size,
        const std::vector<int64_t>& in_size)
    {
        REQUIRES_OK(in_size.size() == out_size.size() && in_size.size() <= 2);
        if (in_size.size() == 1) {
            synchronize_memory<T, Device_out, Device_in>()(arr_out, arr_in, in_size[0]);
        }
        else {
            for (int64_t ii = 0; ii < out_size[0]; ii++) {
                synchronize_memory<T, Device_out, Device_in>()(arr_out + ii * out_size[1], arr_in + ii * in_size[1], in_size[1]);
            }
        }
    }
};

/**
 * @brief Casts memory between devices.
 *
 * This class casts memory between two different devices.
 *
 * @tparam T_out The output data type.
 * @tparam T_in The input data type.
 * @tparam Device_out The output device.
 * @tparam Device_in The input device.
 */
template <typename T_out, typename T_in, typename Device_out, typename Device_in>
struct cast_memory {
    /**
     * @brief Casts memory between devices.
     *
     * This method casts memory between two different devices.
     *
     * @param dev_out The output device.
     * @param dev_in The input device.
     * @param arr_out The output array.
     * @param arr_in The input array.
     * @param size The size of the array.
     */
    void operator()(
        T_out* arr_out,
        const T_in* arr_in,
        const size_t& size);
};


/**
 * @brief Deletes memory on a device.
 *
 * This class deletes memory on a device.
 *
 * @tparam T The type of data in the array.
 * @tparam Device The device.
 */
template <typename T, typename Device>
struct delete_memory {
    /**
     * @brief Deletes memory on a device.
     *
     * This method deletes memory on a device.
     *
     * @param dev The device.
     * @param arr The array to be deleted.
     */
    void operator()(T* arr);
};

} // namespace kernels
} // namespace container

#endif // ATEN_KERNELS_MEMORY_H_