#include <ATen/core/tensor_types.h>

namespace container {

// Overloaded operator<< for the Tensor class.
// Prints the data type of the enum type DataType.
std::ostream& operator<<(std::ostream& os, const DataType& data_type) {
    switch (data_type) {
        case DataType::DT_FLOAT:
            os << "float";
            break;
        case DataType::DT_DOUBLE:
            os << "float64";
            break;
        case DataType::DT_INT:
            os << "int32";
            break;
        case DataType::DT_INT64:
            os << "int64";
            break;
        case DataType::DT_COMPLEX:
            os << "complex<float>";
            break;
        case DataType::DT_COMPLEX_DOUBLE:
            os << "complex<double>";
            break;
    }
    return os;
}

// Overloaded operator<< for the Tensor class.
// Prints the memory type of the enum type DeviceType.
std::ostream& operator<<(std::ostream& os, const DeviceType& device_type) {
    switch (device_type) {
        case DeviceType::CpuDevice:
            os << "cpu";
            break;
    #if __CUDA || __ROCM
        case DeviceType::GpuDevice:
            os << "gpu";
            break;
    #endif
        default:
            os << "unknown";
            break;
    }
    return os;
}

} // namespace container