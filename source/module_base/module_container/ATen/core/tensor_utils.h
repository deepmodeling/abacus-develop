#ifndef ATEN_CORE_TENSOR_UTILS_H_
#define ATEN_CORE_TENSOR_UTILS_H_

#include <ATen/core/tensor.h>
#include <ATen/core/tensor_shape.h>

namespace container {

/**
 *
 * @brief Removes trailing zeros from a string.
 *
 * This function removes any trailing zeros from a given string.
 *
 * @param str The string to remove trailing zeros from.
 *
 * @return A string without trailing zeros.
 */
__inline__
std::string removeTrailingZeros(std::string str) {
    int i = static_cast<int>(str.length()) - 1;
    while (i >= 0 && str[i] == '0') {
        i--;
    }
    if (i == -1) {
        return "0";
    }
    return str.substr(0, i + 1);
}

/**
 *
 * @brief Calculates the length of the longest integer and fractional part of an array of numbers.
 *
 * This function takes an array of numbers and determines the length of the longest integer and fractional part.
 *
 * @tparam T The type of the array.
 *
 * @param arr The array of numbers.
 * @param size The size of the array.
 * @param integer_count The length of the longest integer part.
 * @param fraction_count The length of the longest fractional part.
 *
 * @return The total length of the longest integer and fractional part.
 */

template<typename T>
__inline__
int _get_digit_places(
        const T* arr,
        int size,
        int& integer_count,
        int& fraction_count)
{
    integer_count = 0;
    fraction_count = 0;

    for (int i = 0; i < size; i++) {
        int digits = 0;
        if (arr[i] < 0) {
            digits = log10(-arr[i]) + 1;
            if (digits + 1 > integer_count) {
                integer_count = digits + 1;
            }
        }
        else {
            digits = log10(arr[i]) + 1;
            if (digits > integer_count) {
                integer_count = digits;
            }
        }

        T fraction = arr[i] - std::floor(arr[i]);
        if (fraction == 0) {
            continue;
        }
        std::string str = removeTrailingZeros(std::to_string(fraction));
        digits = str.length() - str.find('.');
        if (digits > fraction_count) {
            fraction_count = digits;
        }
    }

    return integer_count + fraction_count;
}

/**
 *
 * @brief Overloaded function to calculate the length of the longest integer and fractional part of an array of complex numbers.
 *
 * This function is an overloaded version of _get_digit_places for an array of complex numbers.
 *
 * @tparam T The type of the array.
 *
 * @param arr The array of numbers.
 * @param size The size of the array.
 * @param integer_count The length of the longest integer part.
 * @param fraction_count The length of the longest fractional part.
 *
 * @return The total length of the longest integer and fractional part.
 */
template<typename T>
__inline__
int _get_digit_places(
        const std::complex<T>* arr,
        int size,
        int& integer_count,
        int& fraction_count)
{
    return _get_digit_places<T>(reinterpret_cast<const T*>(arr), size * 2, integer_count, fraction_count);
}

/**
 *
 * @brief Output wrapper for a data value with given formatting parameters.
 *
 * @tparam T The type of data to output.
 *
 * @param os The output stream to write the data to.
 * @param data The data value to output.
 * @param digit_width The minimum width for the output.
 * @param fraction_count The number of digits to display after the decimal point.
 */
template <typename T>
__inline__
void _output_wrapper(
        std::ostream& os,
        const T data,
        const int& digit_width,
        const int& fraction_count)
{
    os << std::setw(digit_width) \
       << std::setprecision(fraction_count) << std::fixed << data;
}

/**
 *
 * @brief Output wrapper for a complex data value with given formatting parameters.
 *
 * @tparam T The type of data to output.
 *
 * @param os The output stream to write the data to.
 * @param data The data value to output.
 * @param digit_width The minimum width for the output.
 * @param fraction_count The number of digits to display after the decimal point.
 */
template <typename T>
__inline__
void _output_wrapper(
        std::ostream& os,
        const std::complex<T> data,
        const int& digit_width,
        const int& fraction_count)
{
    // Write the real and imaginary parts of the complex value to the output stream
    // with the specified formatting.
    os << "{";
    os << std::setw(digit_width) \
       << std::setprecision(fraction_count) << std::fixed
       << data.real();
    os << ", ";
    os << std::setw(digit_width) \
       << std::setprecision(fraction_count) << std::fixed
       << data.imag();
    os << "}";
}

/**
 *
 * @brief Output wrapper for an integer data value with given formatting parameters.
 *
 * @param os The output stream to write the data to.
 * @param data The data value to output.
 * @param digit_width The minimum width for the output.
 * @param fraction_count The number of digits to display after the decimal point.
 */
template <>
__inline__
void _output_wrapper(
        std::ostream& os,
        const int data,
        const int& digit_width,
        const int& fraction_count)
{
    os << std::setw(digit_width - 1) \
       << std::setprecision(fraction_count) << std::fixed << data;
}

/**
 *
 * @brief Outputs tensor data to a given output stream.
 *
 * This function outputs tensor data to the specified output stream. It determines the format of the tensor
 * data and prints it accordingly. The format is determined based on the number of dimensions of the tensor.
 *
 * @tparam T The data type of the tensor.
 *
 * @param os The output stream to which the tensor data is to be printed.
 * @param data A pointer to the tensor data.
 * @param shape The shape of the tensor.
 * @param num_elements The total number of elements in the tensor.
*/
template <typename T>
__inline__
void _internal_output(
        std::ostream& os,
        const T * data,
        const TensorShape& shape,
        const int64_t& num_elements)
{
    int integer_count = 0, fraction_count = 0;
    int digit_width = _get_digit_places(data, num_elements, integer_count, fraction_count) + 1;
    if (shape.ndim() == 1) {
        os << "[";
        for (int i = 0; i < num_elements; ++i) {
            _output_wrapper(os, data[i], digit_width, fraction_count);
            if (i != num_elements - 1) {
                os << ",";
            }
        }
        os << "]";
    }
    else if (shape.ndim() == 2) {
        os << "[";
        for (int i = 0; i < shape.dim_size(0); ++i) {
            if (i != 0) os << "       ";
            os << "[";
            for (int j = 0; j < shape.dim_size(1); ++j) {
                _output_wrapper(os, data[i * shape.dim_size(1) + j], digit_width, fraction_count);
                if (j != shape.dim_size(1) - 1) {
                    os << ", ";
                }
            }
            os << "]";
            if (i != shape.dim_size(0) - 1) os << ",\n";
        }
        os << "]";
    }
    else if (shape.ndim() == 3) {
        os << "[";
        for (int i = 0; i < shape.dim_size(0); ++i) {
            if (i != 0) os << "       ";
            os << "[";
            for (int j = 0; j < shape.dim_size(1); ++j) {
                if (j != 0) os << "        ";
                os << "[";
                for (int k = 0; k < shape.dim_size(2); ++k) {
                    _output_wrapper(os, data[i * shape.dim_size(1) * shape.dim_size(2) + j * shape.dim_size(2) + k], digit_width, fraction_count);
                    if (k != shape.dim_size(2) - 1) {
                        os << ", ";
                    }
                }
                os << "]";
                if (j != shape.dim_size(1) - 1) os << ",\n";
            }
            os << "]";
            if (i != shape.dim_size(0) - 1) os << ",\n\n";
        }
        os << "]";
    }
    else {
        for (int64_t i = 0; i < num_elements; ++i) {
            _output_wrapper(os, data[i], 0, 0);
            if (i < num_elements - 1) {
                os << ", ";
            }
        }
    }
}

template <typename T>
T extract(const container::Tensor& tensor) {
    return reinterpret_cast<T*>(tensor.data())[0];
}

} // namespace container

#endif // ATEN_CORE_TENSOR_UTILS_H_