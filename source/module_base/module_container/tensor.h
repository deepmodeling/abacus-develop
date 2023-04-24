#ifndef CONTAINER_TENSOR_H
#define CONTAINER_TENSOR_H

#include <complex>

#include "allocator.h"
#include "tensor_types.h"
#include "tensor_shape.h"
#include "tensor_buffer.h"
#include "kernels/memory_op.h"

namespace container {

/**
 * @brief A multi-dimensional array of elements of a single data type.
 *
 * This class represents a Tensor, which is a fundamental concept in container.
 * A Tensor has a data type, shape, and memory buffer that stores the actual data.
 *
 * This class is not thread-safe and should not be accessed by multiple threads
 * concurrently.
 */
class Tensor {
  public:
    /**
     * @brief Constructor that creates a tensor with the given data type and shape using the default allocator.
     *
     * @param data_type The data type of the tensor.
     * @param shape The shape of the tensor.
     */
    Tensor(DataType data_type, const TensorShape& shape);

    /**
     * @brief Construct a new Tensor object with the given data type, shape and device type.
     *
     * The memory for the tensor is allocated according to the given device type.
     *
     * @param data_type The data type of the tensor.
     * @param shape The shape of the tensor.
     * @param device The data type of the tensor.
     */
    Tensor(DataType data_type, DeviceType device, const TensorShape& shape);

    /**
     * @brief Construct a new Tensor object by copying another Tensor.
     *
     * This constructor performs a deep copy of the data buffer of the other tensor.
     *
     * @param other The tensor to copy from.
     */
    Tensor(const Tensor& other);

    /**
     * @brief Get the data type of the tensor.
     *
     * @return The data type of the tensor.
     */
    DataType data_type() const;

    /**
     * @brief Get the data type of the tensor.
     *
     * @return The data type of the tensor.
     */
    DeviceType device_type() const;

    /**
     * @brief Get the shape of the tensor.
     *
     * @return The shape of the tensor.
     */
    const TensorShape& shape() const;

    /**
     * @brief Get the total number of elements in the tensor.
     *
     * @return The total number of elements in the tensor.
     */
    int64_t NumElements() const;

    /**
     * @brief Get a pointer to the data buffer of the tensor.
     *
     * @return A void pointer to the data buffer of the tensor.
     */
    void* data() const;

    /**
     * @brief Get a typed pointer to the data buffer of the tensor.
     *
     * @tparam T The data type of the pointer to return.
     *
     * @return A typed pointer to the data buffer of the tensor.
     *
     * @note The template parameter `T` must match the data type of the tensor. If `T`
     * does not match the data type of the tensor, the behavior is undefined.
     *
     * @note This function returns a pointer to the first element in the data buffer
     * of the tensor. If the tensor is empty, the behavior is undefined.
     */
    template <typename T>
    T* data() const {
        if ((std::is_same<T, float>::value && data_type_ != DataType::DT_FLOAT) ||
            (std::is_same<T, int>::value && data_type_ != DataType::DT_INT) ||
            (std::is_same<T, int64_t>::value && data_type_ != DataType::DT_INT64) ||
            (std::is_same<T, double>::value && data_type_ != DataType::DT_DOUBLE) ||
            (std::is_same<T, std::complex<float>>::value && data_type_ != DataType::DT_COMPLEX) ||
            (std::is_same<T, std::complex<double>>::value && data_type_ != DataType::DT_COMPLEX_DOUBLE))
        {
            std::cerr << "Tensor data type does not match requested type." << std::endl;
            exit(EXIT_FAILURE);
        }
        return buffer_.base<T>();
    }

    /**
     * @brief Returns the size of a single element for a given data type.
     *
     * @param data_type The data type.
     *
     * @return The size of a single element.
     *
     * @note This function is not thread-safe.
     * @note If an unsupported data type is passed in, an error message will be printed and the program will exit.
     * @note This function returns the size in bytes of a single element of the given data type.
     * If an unsupported data type is provided, the function outputs an error message and exits the program.
     * Supported data types:
     * DT_FLOAT: 4 bytes
     * DT_INT32: 4 bytes
     * DT_DOUBLE: 8 bytes
     * DT_COMPLEX: 8 bytes (2 floats)
     * DT_COMPLEX_DOUBLE: 16 bytes (2 doubles)
     */
    static size_t SizeOfType(DataType data_type) {
        switch (data_type) {
            case DataType::DT_FLOAT:
                return sizeof(float);
            case DataType::DT_INT:
                return sizeof(int32_t);
            case DataType::DT_INT64:
                return sizeof(int64_t);
            case DataType::DT_DOUBLE:
                return sizeof(double);
            case DataType::DT_COMPLEX:
                return sizeof(std::complex<float>);
            case DataType::DT_COMPLEX_DOUBLE:
                return sizeof(std::complex<double>);
            default:
                std::cerr << "Unsupported data type!" << std::endl;
                exit(EXIT_FAILURE);
        }
    }

    /**
     * @brief Get the TensorBuffer object that holds the data of the tensor.
     *
     * @return The TensorBuffer object that holds the data of the tensor.
     */
    const TensorBuffer& buffer() const;

    /**
     * @brief Method to transform data from a given tensor object to the output tensor with a given device type
     *
     * @tparam DEVICE The device type of the returned tensor.
     *
     * @return Tensor A tensor object with data transformed to the output tensor
     */
    template <typename DEVICE>
    Tensor to_device() const {
        // Create output tensor on device
        Tensor output(this->data_type_, DeviceTypeToEnum<DEVICE>::value, this->shape_);

        // Copy data to a specified device
        // TODO: move the memory operator into the tensor_buff class.
        TEMPLATE_ALL_2(this->data_type_, this->device_,
                   op::synchronize_memory_op<T_, DEVICE, DEVICE_>()(
                           output.data<T_>(), this->data<T_>(), this->NumElements()))

        return output;
    }

    /**
     * @brief Method to transform data from a given tensor object to the output tensor with a given data type
     *
     * @tparam T The data type of the returned tensor.
     *
     * @return Tensor A tensor object with data transformed to the output tensor
     */
    template <typename T>
    Tensor cast() const {
        // Create output tensor on device
        Tensor output(DataTypeToEnum<T>::value, this->device_, this->shape_);

        // TODO: error handle of cast memory
        // TODO: move the memory operator into the tensor_buff class.
        // Copy data to a specified device
        TEMPLATE_CZ_2(this->data_type_, this->device_,
                   op::cast_memory_op<T, T_, DEVICE_, DEVICE_>()(
                           output.data<T>(), this->data<T_>(), this->NumElements()))

        return output;
    }

    /**
     * @brief Set all elements in current tensor object to zero.
     */
    void zero();

    /**
     * @brief Set all elements in current tensor object to zero.
     *
     * @param shape The new shape of the tensor.
     *
     * @note There can be one -1 dimension in the input shape, indicates the auto reshape.
     */
    void reshape(TensorShape shape);

    /**
     * @brief Return a new Tensor slice starting at the specified indices with the given size.
     *
     * @param start A vector of integers representing the starting indices of the slice.
     * @param size A vector of integers representing the size of the slice along each dimension.
     *
     * @return A new Tensor slice.
     *
     * @note Currently, this method only supports tensors with a ndim of less than or equal to 3.
     */
    Tensor slice(const std::vector<int>& start, const std::vector<int>& size) const;

    /**
     * @brief Resize the tensor to the new shape.
     *
     * The tensor will be resized to the new shape, and its data buffer will be reallocated
     * if necessary. If the new shape has a different number of elements than the current
     * shape, the data buffer will be reallocated. If the new shape has the same number of
     * elements as the current shape, the data buffer will not be reallocated, but the
     * shape of the tensor will be updated.
     *
     * @param new_shape The new shape of the tensor.
     *
     * @note This method will automatically zero the resized tensor object.
     */
    virtual void resize(const TensorShape& new_shape);

    /**
     * @brief Constructor that creates a tensor with the given data pointer,
     * data type, device type and shape.
     *
     * This tensor does not own memory.
     *
     * @param data The data pointer.
     * @param data_type The data type of the tensor.
     * @param device The data type of the tensor.
     * @param shape The shape of the tensor.
     */
    Tensor(void * data, DataType data_type, DeviceType device, const TensorShape& shape);

    /**
     * @brief Get the Allocator object according to the given device type.
     *
     * @param device The device type.
     *
     * @return The related Allocator class pointer.
     */
    static Allocator* GetAllocator(DeviceType device);

private:

    /**
     * @brief The data type of the tensor.
     */
    DataType data_type_;

    /**
     * @brief The device type of the tensor.
     */
    DeviceType device_;

    /**
     * @brief The shape of the tensor.
     */
    TensorShape shape_;

    /**
     * @brief The allocator used to allocate the memory for the tensor.
     */
    Allocator* allocator_;

    /**
     * @brief The TensorBuffer object that holds the data of the tensor.
     */
    TensorBuffer buffer_;


};

/**
 * @brief Overloaded operator<< for the Tensor class.
 *
 * Prints the data of the Tensor object to the given output stream.
 *
 * @param os The output stream to write to.
 * @param tensor The Tensor object to print.
 *
 * @return The output stream.
 */
std::ostream& operator<<(std::ostream& os, const Tensor& tensor);

} // namespace container

#endif // CONTAINER_TENSOR_H
