#include "tensor_shape.h"

namespace container {
/**
 * @brief Namespace containing constants for default constructor
 */
namespace {
    /**
     * @brief Default size of a dimension
     */
    constexpr int kDefaultDimSize = 0;
} // namespace

// Default constructor for TensorShape class
// Initializes TensorShape with default dimensions
TensorShape::TensorShape() : dims_(kDefaultDimSize) {}

// Constructor for TensorShape class
TensorShape::TensorShape(std::initializer_list<int> dims) : dims_(dims) {}

// Constructor for TensorShape class
TensorShape::TensorShape(const std::vector<int>& dims) : dims_(dims) {}

// Copy constructor for TensorShape class
TensorShape::TensorShape(const TensorShape& other) = default;

// Get size of a specific dimension in the tensor
int TensorShape::dim_size(int dim) const {
    return dims_[dim];
}

// Get all dimension sizes in the tensor
const std::vector<int>& TensorShape::dims() const {
    return dims_;
}

// Get ndim of the tensor, i.e., number of dimensions
unsigned int TensorShape::ndim() const {
    return dims_.size();
}

// Returns the total number of elements in the shape.
int64_t TensorShape::NumElements() const {
    int64_t num_elements = 1;
    for (int i = 0; i < this->ndim(); ++i) {
        num_elements *= dims_[i];
    }
    return num_elements;
}

// Modify size of a specific dimension in the tensor
void TensorShape::set_dim_size(int dim, int size) {
    dims_[dim] = size;
}

// Add a new dimension to the tensor
void TensorShape::add_dim(int size) {
    dims_.push_back(size);
}

// Remove a dimension from the tensor
void TensorShape::remove_dim(int dim) {
    dims_.erase(dims_.begin() + dim);
}

// Overload the == operator to compare two TensorShape objects
bool TensorShape::operator==(const TensorShape& other) const {
    return dims_ == other.dims_;
}

// Overload the != operator to compare two TensorShape objects
bool TensorShape::operator!=(const TensorShape& other) const {
    return dims_ != other.dims_;
}

// Overload the << operator to print the tensor shape
std::ostream& operator<<(std::ostream& os, const TensorShape& shape) {
    os << "[";
    for (int i = 0; i < shape.ndim(); ++i) {
        os << shape.dims()[i];
        if (i < shape.ndim() - 1) {
            os << ",";
        }
    }
    os << "]";
    return os;
}

} // namespace container