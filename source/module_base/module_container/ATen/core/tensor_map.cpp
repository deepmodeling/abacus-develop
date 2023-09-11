#include <ATen/core/tensor_map.h>

namespace container {
// Constructor that creates a tensor with the given data pointer, data type, device type and shape.
TensorMap::TensorMap(void *data, DataType data_type, DeviceType device, const TensorShape &shape)
        : Tensor(data_type, device, {})
{
    this->shape_ = shape;
    this->buffer_ = new TensorBuffer(data);
}

void TensorMap::resize(const TensorShape &new_shape) {
    throw std::logic_error("TensorMap object does not support the resize method.");
}

} // namespace container