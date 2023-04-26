#include <cmath>
#include <iostream>
#include <complex>

#include "tensor.h"
#include "tensor_utils.h"
#include "cpu_allocator.h"
#if defined(__CUDA)
#include "gpu_allocator.h"
#endif // __CUDA
namespace container {

// Constructor that creates a tensor with the given data type and shape using the default allocator.
Tensor::Tensor(DataType data_type, const TensorShape& shape)
        : data_type_(data_type),
          shape_(shape),
          device_(DeviceType::CpuDevice),
          allocator_(GetAllocator(device_)),
          buffer_(allocator_, allocator_->allocate(shape.NumElements() * SizeOfType(data_type))) {}

// Constructor that creates a tensor with the given data pointer, data type, device type and shape.
Tensor::Tensor(void *data, DataType data_type, DeviceType device, const TensorShape &shape)
        : data_type_(data_type),
          shape_(shape),
          device_(device),
          allocator_(GetAllocator(device_)),
          buffer_(data) {}

// Construct a new Tensor object with the given data type and shape.
Tensor::Tensor(DataType data_type, DeviceType device, const TensorShape& shape)
        : data_type_(data_type),
          shape_(shape),
          device_(device),
          allocator_(GetAllocator(device_)),
          buffer_(allocator_, allocator_->allocate(shape.NumElements() * SizeOfType(data_type))) {}

// Construct a new Tensor object by copying another Tensor.
Tensor::Tensor(const Tensor& other)
        : data_type_(other.data_type_),
          shape_(other.shape_),
          device_(other.device_),
          allocator_(other.allocator_),
          buffer_(allocator_, allocator_->allocate(shape_.NumElements() * SizeOfType(data_type_)))
{
    TEMPLATE_ALL_2(data_type_, device_,
            op::synchronize_memory_op<T_, DEVICE_, DEVICE_>()(
                    this->data<T_>(), other.data<T_>(), this->NumElements()))
}

// Get the data type of the tensor.
DataType Tensor::data_type() const { return data_type_; }

// Get the device type of the tensor.
DeviceType Tensor::device_type() const { return device_; }

// Get the shape of the tensor.
const TensorShape& Tensor::shape() const { return shape_; }

// Get the total number of elements in the tensor.
int64_t Tensor::NumElements() const { return shape_.NumElements(); }

// Get a pointer to the data buffer of the tensor.
void* Tensor::data() const { return buffer_.data(); }

// Get the TensorBuffer object that holds the data of the tensor.
const TensorBuffer& Tensor::buffer() const { return buffer_; }

// Get the Allocator object according to the given device type.
Allocator* Tensor::GetAllocator(DeviceType device) {
    Allocator * allocator;
    if (device == DeviceType::CpuDevice) {
        allocator = new CPUAllocator();
    }
#if defined(__CUDA)
    else if (device == DeviceType::GpuDevice) {
        allocator = new GPUAllocator();
    }
#endif // __CUDA
    else {
        std::cerr << "Tensor device type " << device << " does not match requested type." << std::endl;
        exit(EXIT_FAILURE);
    }
    return allocator;
}

// Set the tensor to zero
void Tensor::zero() {
    TEMPLATE_ALL_2(this->data_type_, this->device_,
            op::set_memory_op<T_, DEVICE_>()(this->data<T_>(), 0, this->NumElements()))
}

// Reshape the current tensor
void Tensor::reshape(TensorShape shape) {
    // check the -1 dimension
    int num = 1, auto_shape = 0, dim_count = -1, dim_idx = -1;

    for (int dim : shape.dims()) {
        dim_count++;
        if (dim < 1 && dim != -1) {
            throw std::invalid_argument("Invalid shape, dim of tensor must >= 1 or equal to -1(auto shape).");
        }
        if (dim == -1) {
            auto_shape++;
            dim_idx = dim_count;
        }
        num *= dim;
    }
    // more than one -1 dimension.
    if (auto_shape > 1) {
        throw std::invalid_argument("Invalid shape, there can be only one -1 dim in TensorShape object.");
    }
    // auto reshape
    if (auto_shape == 1) {
        int dim_ = static_cast<int>(this->NumElements()) / (-num);
        if (dim_ < 1 || -dim_ * num != this->NumElements()) {
            throw std::invalid_argument("Invalid shape, total number of elements does not match!");
        }
        shape.set_dim_size(dim_idx, dim_);
    }
    else {
        if (num != this->NumElements()) {
            throw std::invalid_argument("Invalid shape, total number of elements does not match!");
        }
    }
    this->shape_ = shape;
}

// Slice the current tensor object.
Tensor Tensor::slice(const std::vector<int> &start, const std::vector<int> &size) const {
    // check the ndim of input shape
    if (start.size() > 3 || size.size() > 3) {
        throw std::invalid_argument("TensorSlice: The slice method only supports tensor ranks that are less than or equal to 2.");
    }
    // check the dimension size
    if (start.size() != shape_.ndim() || size.size() != shape_.ndim()) {
        throw std::invalid_argument("TensorSlice: start and size vectors must have same length as number of dimensions");
    }

    // check the boundary
    for (int i = 0; i < start.size(); i++) {
        if (start[i] < 0 || start[i] >= shape_.dim_size(i)) {
            throw std::invalid_argument("TensorSlice: start index is out of bounds");
        }
        if (size[i] < 0 || start[i] + size[i] > shape_.dim_size(i)) {
            throw std::invalid_argument("TensorSlice: size is out of bounds");
        }
    }

    // set the output shape of the current tensor.
    TensorShape output_shape = shape_;
    for (int i = 0; i < start.size(); i++) {
        output_shape.set_dim_size(i, size[i]);
    }
    Tensor output(this->data_type_, this->device_, output_shape);

    // TODO: implement the data copy.
    // copy the data from the input tensor to the output tensor
    unsigned int ndim = shape_.ndim();
    if (ndim == 1) {
        TEMPLATE_ALL_2(this->data_type_, this->device_,
                       op::synchronize_memory_op<T_, DEVICE_, DEVICE_>()(
                               output.data<T_>(), this->data<T_>() + start[0], size[0]))
    }
    else if (ndim == 2) {
        for (int i = 0; i < size[0]; i++) {
            int offset = (start[0] + i) * shape_.dim_size(1) + start[1];
            int offset_out = i * size[1];
            TEMPLATE_ALL_2(this->data_type_, this->device_,
                           op::synchronize_memory_op<T_, DEVICE_, DEVICE_>()(
                                   output.data<T_>() + offset_out, this->data<T_>() + offset, size[1]))
        }
    }
    else if (ndim == 3) {
        for (int i = 0; i < size[0]; i++) {
            for (int j = 0; j < size[1]; j++) {
                int offset = (i + start[0]) * shape_.dim_size(1) * shape_.dim_size(2) +
                        (j + start[1]) * shape_.dim_size(2) + start[2];
                int offset_out = i * size[1] * size[2] + j * size[2];
                TEMPLATE_ALL_2(this->data_type_, this->device_,
                               op::synchronize_memory_op<T_, DEVICE_, DEVICE_>()(
                                       output.data<T_>() + offset_out, this->data<T_>() + offset, size[1]))
            }
        }
    }
    return output;
}

// Resize tensor object with the given tensor_shape
void Tensor::resize(const TensorShape& new_shape) {
    if (shape_ == new_shape) {
        return;
    }
    buffer_.resize(new_shape.NumElements() * Tensor::SizeOfType(data_type_));
    shape_ = new_shape;
    this->zero();
}

// Overloaded operator<< for the Tensor class.
std::ostream& operator<<(std::ostream& os, const Tensor& tensor) {
    std::ios::fmtflags flag(os.flags());
    const int64_t num_elements = tensor.NumElements();
    const DataType data_type = tensor.data_type();
    const DeviceType device_type = tensor.device_type();
    const TensorShape& shape = tensor.shape();

    // Copy the data from device to host for output
    auto * data_ = tensor.data();
#if __CUDA || __ROCM
    if (device_type != DeviceType::CpuDevice) {
        data_ = malloc(num_elements * Tensor::SizeOfType(data_type));
        // Copy data to a specified device
        TEMPLATE_ALL_2(data_type, device_type,
                       container::op::synchronize_memory_op<T_, DEVICE_CPU, DEVICE_>()(
                               reinterpret_cast<T_ *>(data_), tensor.data<T_>(), num_elements))
    }
#endif

    os << "Tensor(";
    os << "shape=[";
    for (int i = 0; i < shape.ndim(); ++i) {
        os << shape.dim_size(i);
        if (i < shape.ndim() - 1) {
            os << ",";
        }
    }
    os << "], data_type=" << data_type;
    os << ", device_type=" << device_type;
    os << ", owns_memory=" << tensor.buffer().OwnsMemory();
    os << ", buffer=\narray(";
    switch (data_type) {
        case DataType::DT_FLOAT: {
            const auto* data = static_cast<const float*>(data_);
            _internal_output(os, data, shape, num_elements);
            break;
        }
        case DataType::DT_DOUBLE: {
            const auto* data = static_cast<const double*>(data_);
            _internal_output(os, data, shape, num_elements);
            break;
        }
        case DataType::DT_INT: {
            const auto* data = static_cast<const int*>(data_);
            _internal_output(os, data, shape, num_elements);
            break;
        }
        case DataType::DT_INT64: {
            const auto* data = static_cast<const int64_t*>(data_);
            _internal_output(os, data, shape, num_elements);
            break;
        }
        case DataType::DT_COMPLEX: {
            const auto* data = static_cast<const std::complex<float>*>(data_);
            _internal_output(os, data, shape, num_elements);
            break;
        }
        case DataType::DT_COMPLEX_DOUBLE: {
            const auto* data = static_cast<const std::complex<double>*>(data_);
            _internal_output(os, data, shape, num_elements);
            break;
        }
    }
    os << "))\n";

#if __CUDA || __ROCM
    // delete the temporary data
    if (device_type != DeviceType::CpuDevice) {
        free(data_);
    }
#endif
    // restore the os settings
    os.flags(flag);
    return os;
}

} // namespace container
