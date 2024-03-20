#ifndef ATEN_CORE_TENSOR_MAP_H_
#define ATEN_CORE_TENSOR_MAP_H_

#include <complex>

#include <ATen/core/tensor.h>

namespace container {

/**
 * @brief A multi-dimensional reference array of elements of a single data type.
 *
 * This class represents a Tensor, which is a fundamental concept in container.
 * A Tensor has a data type, shape, and memory buffer that stores the actual data.
 *
 * This class is not thread-safe and should not be accessed by multiple threads
 * concurrently.
 */
 class TensorMap : public Tensor {
 public:

    /**
     * @brief Constructor that map the given data pointer to a tensor object with the given
     * data type, device type and shape.
     *
     * This tensor does not own memory.
     *
     * @param data The data pointer.
     * @param data_type The data type of the tensor.
     * @param device The data type of the tensor.
     * @param shape The shape of the tensor.
     */
    TensorMap(void *data, DataType data_type, DeviceType device, const TensorShape &shape);

    /**
     * @brief Constructs a TensorMap from the provided data, with the specified shape.
     *
     * This constructor creates a TensorMap that references the data pointed to by the
     * 'data' parameter, and associates it with the given Tensor 'other' and TensorShape 'shape'.
     * The TensorMap allows access to the data with the specified shape without owning the data.
     *
     * @param data A pointer to the data to be referenced by the TensorMap.
     * @param other The Tensor object to be associated with the TensorMap.
     * @param shape The shape of the data referenced by the TensorMap.
     */
    TensorMap(void *data, const Tensor& other, const TensorShape& shape);

    /**
     * @brief Constructs a TensorMap from the provided data, using the shape of the provided Tensor.
     *
     * This constructor creates a TensorMap that references the data pointed to by the 'data' parameter,
     * and associates it with the given Tensor 'other'. The shape of the data is determined by the shape
     * of the 'other' Tensor. The TensorMap allows access to the data with the shape of 'other' without
     * owning the data.
     *
     * @param data A pointer to the data to be referenced by the TensorMap.
     * @param other The Tensor object to be associated with the TensorMap, which defines the shape.
     */
    TensorMap(void *data, const Tensor& other);

};

} // namespace container

#endif // ATEN_CORE_TENSOR_MAP_H_
