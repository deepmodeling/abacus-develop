#ifndef ATEN_CORE_TENSOR_SHAPE_H_
#define ATEN_CORE_TENSOR_SHAPE_H_

#include <vector>
#include <iostream>
#include <initializer_list>

namespace container {

/**
 * @brief A class for representing the shape of a tensor.
 */
class TensorShape {
public:
    /**
     * @brief Default constructor.
     */
    TensorShape();

    /**
     * @brief Constructor with an initializer list of integers.
     * @param dims An initializer list of integers representing the dimensions of the tensor.
     */
    TensorShape(std::initializer_list<int64_t> dims);

    /**
     * @brief Constructor with a vector of integers.
     * @param dims A vector of integers representing the dimensions of the tensor.
     */
    TensorShape(const std::vector<int64_t>& dims);

    /**
     * @brief Copy constructor.
     * @param other The TensorShape object to be copied.
     */
    TensorShape(const TensorShape& other);

    /**
     * @brief Get the size of a dimension in the tensor.
     * @param dim The index of the dimension.
     * @return The size of the specified dimension.
     */
    int64_t dim_size(int dim) const;

    /**
     * @brief Get all dimension sizes in the tensor.
     * @return A const reference to the vector of dimension sizes.
     */
    const std::vector<int64_t>& dims() const;

    const std::vector<int64_t>& strides() const;

    /**
     * @brief Get the ndim of the tensor.
     * @return The number of dimensions in the tensor.
     */
    unsigned int ndim() const;

    /**
     * @brief Modify the size of a dimension in the tensor.
     * @param dim The index of the dimension to be modified.
     * @param size The new size of the dimension.
     */
    void set_dim_size(int dim, int64_t size);

    /**
     * @brief Add a new dimension to the tensor.
     * @param size The size of the new dimension.
     */
    void add_dim(int64_t size);

    /**
     * @brief Remove a dimension from the tensor.
     * @param dim The index of the dimension to be removed.
     */
    void remove_dim(int dim);

    /**
    * @brief Returns the total number of elements in the shape.
    *
    * @return int64_t The number of elements.
    */
    int64_t NumElements() const;

    /**
     * @brief Overload the == operator to compare two TensorShape objects.
     * @param other The other TensorShape object to be compared.
     * @return True if the two objects have the same dimensions, false otherwise.
     */
    bool operator==(const TensorShape& other) const;

    /**
     * @brief Overload the != operator to compare two TensorShape objects.
     * @param other The other TensorShape object to be compared.
     * @return True if the two objects have different dimensions, false otherwise.
     */
    bool operator!=(const TensorShape& other) const;

private:
    std::vector<int64_t> dims_ = {};     // Save dimension sizes of the tensor
    // Note: strides are not always equals to the dimension sizes.
    // The strides specifies the number of elements to step in each dimension when traversing a tensor.
    // There could be some sparse region in the tensor, and the strides will be larger than the dimension sizes.
    // For example, given a 2D tensor with shape [3, 4], 
    // and the strides could be [6, 1] if the actual data is stored in a 1D array with size 18 [3, 6].
    // The strides could also be [12, 3] if the actual data is stored in a 1D array with size 36 [3, 12].
    // strides can only be modified by the TensorMap object.
    std::vector<int64_t> strides_ = {};  // Save dimension strides of the tensor

    /**
     * @brief Compute the strides of the tensor.
     */
    std::vector<int64_t> get_strides_(const std::vector<int64_t>& dim);
};

/**
 * @brief Overload the << operator to print TensorShape objects.
 * @param os The output stream to which the tensor shape is printed.
 * @param shape The TensorShape object to be printed.
 * @return A reference to the output stream.
 */
std::ostream& operator<<(std::ostream& os, const TensorShape& shape);

} // container

#endif  // ATEN_CORE_TENSOR_SHAPE_H_
