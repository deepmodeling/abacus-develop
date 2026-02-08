/**
 * @file pybind_utils.h
 * @brief Utility functions for pybind11 array operations in PyABACUS
 *
 * This header provides common utilities for:
 * - Array validation (dimension, size checks)
 * - Array pointer access
 * - Array conversion (copy and zero-copy views)
 * - Safe Python function calls
 * - Buffer helper class
 */

#ifndef PYABACUS_UTILS_PYBIND_UTILS_H
#define PYABACUS_UTILS_PYBIND_UTILS_H

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <complex>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#include <cstring>

namespace py = pybind11;

namespace pyabacus {
namespace utils {

// ============================================================================
// Array Validation Functions
// ============================================================================

/**
 * @brief Check that a numpy array is 1-dimensional
 * @param arr The numpy array to check
 * @param name Name of the array for error messages
 * @throws std::runtime_error if array is not 1D
 */
template <typename T>
void check_1d_array(const py::array_t<T>& arr, const std::string& name = "array")
{
    py::buffer_info info = arr.request();
    if (info.ndim != 1)
    {
        throw std::runtime_error(name + " must be 1-dimensional, got " +
                                 std::to_string(info.ndim) + " dimensions");
    }
}

/**
 * @brief Check that a numpy array is 2-dimensional
 * @param arr The numpy array to check
 * @param name Name of the array for error messages
 * @throws std::runtime_error if array is not 2D
 */
template <typename T>
void check_2d_array(const py::array_t<T>& arr, const std::string& name = "array")
{
    py::buffer_info info = arr.request();
    if (info.ndim != 2)
    {
        throw std::runtime_error(name + " must be 2-dimensional, got " +
                                 std::to_string(info.ndim) + " dimensions");
    }
}

/**
 * @brief Check that a numpy array has the expected size
 * @param arr The numpy array to check
 * @param expected_size Expected total number of elements
 * @param name Name of the array for error messages
 * @throws std::runtime_error if size doesn't match
 */
template <typename T>
void check_array_size(const py::array_t<T>& arr, size_t expected_size,
                      const std::string& name = "array")
{
    if (static_cast<size_t>(arr.size()) != expected_size)
    {
        throw std::runtime_error(name + " must have size " +
                                 std::to_string(expected_size) + ", got " +
                                 std::to_string(arr.size()));
    }
}

/**
 * @brief Check that a numpy array has the expected shape
 * @param arr The numpy array to check
 * @param expected_shape Expected shape as vector
 * @param name Name of the array for error messages
 * @throws std::runtime_error if shape doesn't match
 */
template <typename T>
void check_array_shape(const py::array_t<T>& arr,
                       const std::vector<ssize_t>& expected_shape,
                       const std::string& name = "array")
{
    py::buffer_info info = arr.request();
    if (info.ndim != static_cast<ssize_t>(expected_shape.size()))
    {
        throw std::runtime_error(name + " has wrong number of dimensions");
    }
    for (size_t i = 0; i < expected_shape.size(); ++i)
    {
        if (info.shape[i] != expected_shape[i])
        {
            throw std::runtime_error(name + " has wrong shape at dimension " +
                                     std::to_string(i));
        }
    }
}

/**
 * @brief Check that two arrays have the same size
 * @param arr1 First array
 * @param arr2 Second array
 * @param name1 Name of first array
 * @param name2 Name of second array
 * @throws std::runtime_error if sizes don't match
 */
template <typename T1, typename T2>
void check_same_size(const py::array_t<T1>& arr1, const py::array_t<T2>& arr2,
                     const std::string& name1 = "array1",
                     const std::string& name2 = "array2")
{
    if (arr1.size() != arr2.size())
    {
        throw std::runtime_error(name1 + " and " + name2 +
                                 " must have the same size");
    }
}

// ============================================================================
// Array Pointer Access
// ============================================================================

/**
 * @brief Get mutable pointer to array data
 * @param arr The numpy array
 * @return Pointer to the array data
 */
template <typename T>
T* get_array_ptr(py::array_t<T>& arr)
{
    py::buffer_info info = arr.request();
    return static_cast<T*>(info.ptr);
}

/**
 * @brief Get const pointer to array data
 * @param arr The numpy array
 * @return Const pointer to the array data
 */
template <typename T>
const T* get_array_ptr(const py::array_t<T>& arr)
{
    py::buffer_info info = arr.request();
    return static_cast<const T*>(info.ptr);
}

/**
 * @brief Get pointer to array data with 1D validation
 * @param arr The numpy array
 * @param name Name for error messages
 * @return Pointer to the array data
 */
template <typename T>
T* get_1d_array_ptr(py::array_t<T>& arr, const std::string& name = "array")
{
    check_1d_array(arr, name);
    return get_array_ptr(arr);
}

/**
 * @brief Get const pointer to array data with 1D validation
 * @param arr The numpy array
 * @param name Name for error messages
 * @return Const pointer to the array data
 */
template <typename T>
const T* get_1d_array_ptr(const py::array_t<T>& arr, const std::string& name = "array")
{
    check_1d_array(arr, name);
    return get_array_ptr(arr);
}

// ============================================================================
// Array Conversion (Copy)
// ============================================================================

/**
 * @brief Create numpy array from std::vector (copy)
 * @param vec Source vector
 * @return New numpy array with copied data
 */
template <typename T>
py::array_t<T> numpy_from_vector_copy(const std::vector<T>& vec)
{
    py::array_t<T> result(vec.size());
    std::memcpy(result.mutable_data(), vec.data(), vec.size() * sizeof(T));
    return result;
}

/**
 * @brief Create 1D numpy array from raw pointer (copy)
 * @param ptr Source pointer
 * @param size Number of elements
 * @return New numpy array with copied data
 */
template <typename T>
py::array_t<T> numpy_from_ptr_copy(const T* ptr, size_t size)
{
    py::array_t<T> result(size);
    std::memcpy(result.mutable_data(), ptr, size * sizeof(T));
    return result;
}

/**
 * @brief Create 2D numpy array from raw pointer (copy, row-major)
 * @param ptr Source pointer
 * @param nrow Number of rows
 * @param ncol Number of columns
 * @return New numpy array with copied data
 */
template <typename T>
py::array_t<T> numpy_from_ptr_copy_2d(const T* ptr, size_t nrow, size_t ncol)
{
    std::vector<ssize_t> shape = {static_cast<ssize_t>(nrow),
                                   static_cast<ssize_t>(ncol)};
    py::array_t<T> result(shape);
    std::memcpy(result.mutable_data(), ptr, nrow * ncol * sizeof(T));
    return result;
}

/**
 * @brief Create numpy array with specified shape from raw pointer (copy)
 * @param ptr Source pointer
 * @param shape Shape of the output array
 * @return New numpy array with copied data
 */
template <typename T>
py::array_t<T> numpy_from_ptr_copy_nd(const T* ptr, const std::vector<ssize_t>& shape)
{
    py::array_t<T> result(shape);
    size_t total_size = 1;
    for (auto dim : shape)
    {
        total_size *= static_cast<size_t>(dim);
    }
    std::memcpy(result.mutable_data(), ptr, total_size * sizeof(T));
    return result;
}

// ============================================================================
// Zero-Copy Views with Lifetime Management
// ============================================================================

/**
 * @brief Create numpy array view with keepalive (zero-copy)
 *
 * Creates a numpy array that references existing memory without copying.
 * The owner_ptr shared_ptr ensures the underlying data stays alive.
 *
 * @param ptr Pointer to the data
 * @param shape Shape of the array
 * @param owner_ptr Shared pointer to the owner object (keeps data alive)
 * @return Numpy array view
 */
template <typename T, typename Owner>
py::array_t<T> numpy_view_with_keepalive(T* ptr,
                                          const std::vector<ssize_t>& shape,
                                          std::shared_ptr<Owner> owner_ptr)
{
    // Calculate strides for C-contiguous array
    std::vector<ssize_t> strides(shape.size());
    ssize_t stride = sizeof(T);
    for (int i = static_cast<int>(shape.size()) - 1; i >= 0; --i)
    {
        strides[i] = stride;
        stride *= shape[i];
    }

    // Create capsule to prevent deallocation
    py::capsule free_when_done(owner_ptr.get(), [](void*) {
        // The shared_ptr prevents deallocation, capsule just holds reference
    });

    return py::array_t<T>(shape, strides, ptr, free_when_done);
}

/**
 * @brief Create 1D numpy array view (zero-copy, no ownership management)
 *
 * WARNING: Caller must ensure the data outlives the returned array!
 *
 * @param ptr Pointer to the data
 * @param size Number of elements
 * @return Numpy array view
 */
template <typename T>
py::array_t<T> numpy_view_1d(T* ptr, size_t size)
{
    return py::array_t<T>({static_cast<ssize_t>(size)},
                          {sizeof(T)},
                          ptr,
                          py::none());
}

/**
 * @brief Create 2D numpy array view (zero-copy, no ownership management)
 *
 * WARNING: Caller must ensure the data outlives the returned array!
 *
 * @param ptr Pointer to the data
 * @param nrow Number of rows
 * @param ncol Number of columns
 * @return Numpy array view (row-major)
 */
template <typename T>
py::array_t<T> numpy_view_2d(T* ptr, size_t nrow, size_t ncol)
{
    return py::array_t<T>(
        {static_cast<ssize_t>(nrow), static_cast<ssize_t>(ncol)},
        {static_cast<ssize_t>(ncol * sizeof(T)), sizeof(T)},
        ptr,
        py::none()
    );
}

// ============================================================================
// Safe Python Function Calls
// ============================================================================

/**
 * @brief Call a Python function safely with exception handling
 *
 * Wraps Python function calls to provide better error messages
 * and handle GIL properly.
 *
 * @param func Python function to call
 * @param args Arguments to pass to the function
 * @return Return value from the function
 */
template <typename Ret, typename... Args>
Ret call_python_safe(const py::function& func, Args&&... args)
{
    try
    {
        return func(std::forward<Args>(args)...).template cast<Ret>();
    }
    catch (const py::error_already_set& e)
    {
        throw std::runtime_error(std::string("Python callback error: ") + e.what());
    }
    catch (const py::cast_error& e)
    {
        throw std::runtime_error(std::string("Python return type error: ") + e.what());
    }
}

/**
 * @brief Call a Python function that returns void
 * @param func Python function to call
 * @param args Arguments to pass to the function
 */
template <typename... Args>
void call_python_safe_void(const py::function& func, Args&&... args)
{
    try
    {
        func(std::forward<Args>(args)...);
    }
    catch (const py::error_already_set& e)
    {
        throw std::runtime_error(std::string("Python callback error: ") + e.what());
    }
}

// ============================================================================
// Buffer Helper Class
// ============================================================================

/**
 * @brief Helper class for working with numpy buffer info
 *
 * Provides convenient access to buffer properties and validation methods.
 */
template <typename T>
struct BufferHelper
{
    py::buffer_info info;
    T* ptr;
    size_t size;

    /**
     * @brief Construct from numpy array
     * @param arr The numpy array
     */
    explicit BufferHelper(py::array_t<T>& arr)
        : info(arr.request())
        , ptr(static_cast<T*>(info.ptr))
        , size(static_cast<size_t>(arr.size()))
    {
    }

    /**
     * @brief Construct from const numpy array
     * @param arr The numpy array
     */
    explicit BufferHelper(const py::array_t<T>& arr)
        : info(arr.request())
        , ptr(static_cast<T*>(info.ptr))
        , size(static_cast<size_t>(arr.size()))
    {
    }

    /**
     * @brief Require array to be 1-dimensional
     * @param name Name for error messages
     * @throws std::runtime_error if not 1D
     */
    void require_1d(const std::string& name = "array") const
    {
        if (info.ndim != 1)
        {
            throw std::runtime_error(name + " must be 1-dimensional");
        }
    }

    /**
     * @brief Require array to be 2-dimensional
     * @param name Name for error messages
     * @throws std::runtime_error if not 2D
     */
    void require_2d(const std::string& name = "array") const
    {
        if (info.ndim != 2)
        {
            throw std::runtime_error(name + " must be 2-dimensional");
        }
    }

    /**
     * @brief Require array to have specific size
     * @param expected Expected size
     * @param name Name for error messages
     * @throws std::runtime_error if size doesn't match
     */
    void require_size(size_t expected, const std::string& name = "array") const
    {
        if (size != expected)
        {
            throw std::runtime_error(name + " must have size " +
                                     std::to_string(expected) + ", got " +
                                     std::to_string(size));
        }
    }

    /**
     * @brief Require array to have specific shape
     * @param expected_shape Expected shape
     * @param name Name for error messages
     * @throws std::runtime_error if shape doesn't match
     */
    void require_shape(const std::vector<ssize_t>& expected_shape,
                       const std::string& name = "array") const
    {
        if (info.ndim != static_cast<ssize_t>(expected_shape.size()))
        {
            throw std::runtime_error(name + " has wrong number of dimensions");
        }
        for (size_t i = 0; i < expected_shape.size(); ++i)
        {
            if (info.shape[i] != expected_shape[i])
            {
                throw std::runtime_error(name + " has wrong shape at dimension " +
                                         std::to_string(i));
            }
        }
    }

    /**
     * @brief Get number of dimensions
     * @return Number of dimensions
     */
    ssize_t ndim() const { return info.ndim; }

    /**
     * @brief Get shape at dimension i
     * @param i Dimension index
     * @return Size at dimension i
     */
    ssize_t shape(size_t i) const { return info.shape[i]; }

    /**
     * @brief Get number of rows (for 2D arrays)
     * @return Number of rows
     */
    ssize_t nrow() const { return info.shape[0]; }

    /**
     * @brief Get number of columns (for 2D arrays)
     * @return Number of columns
     */
    ssize_t ncol() const { return info.shape[1]; }
};

// ============================================================================
// Convenience Type Aliases
// ============================================================================

using DoubleBuffer = BufferHelper<double>;
using ComplexBuffer = BufferHelper<std::complex<double>>;
using IntBuffer = BufferHelper<int>;

} // namespace utils
} // namespace pyabacus

#endif // PYABACUS_UTILS_PYBIND_UTILS_H
