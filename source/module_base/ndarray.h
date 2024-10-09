/**
 * @file NDArray.h
 * @author your name (you@domain.com)
 * @brief under the restriction of C++11, a simple alternative to std::vector<T> + std::mdspan. In module_base/module_container/ATen/tensor.h, tensor class provides a cross-device container, but std::string is not supported. Therefore, this class is to provide a general (but CPU-only) container for multi-dimensional data. It can easily convert to ontainer::Tensor.
 * @version 0.1
 * @date 2024-04-24
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef NDARRAY_H
#define NDARRAY_H

#include <vector>
#include <cassert>
#include <iostream>
#include <algorithm>
#include <type_traits>
#include <numeric>
// for heterogeneous computing, we can use ATen::Tensor
//#include "./module_container/ATen/tensor.h"

/**
 * @brief under the restriction of C++11, a simple alternative to std::vector<T> + std::mdspan. In module_base/module_container/ATen/tensor.h, tensor class provides a cross-device container, but std::string is not supported. Therefore, this class is to provide a general (but CPU-only) container for multi-dimensional data. It can easily convert to container::Tensor.
 * 
 * @tparam T 
 */
template<typename T>
class NDArray
{
    // align with STL container implementation, there are several functions compulsory to be implemented
    // constructor: default, copy, move, initializer_list
    // operator: =, ==, !=, <, <=, >, >=
    // iterator: begin, cbegin, end, cend
    // capacity: size, empty, max_size, reserve, shrink_to_fit
    // element access: [], at, front, back, data
    // modifiers: clear, insert, erase, push_back, pop_back, resize, swap
    // allocator: get_allocator
public:
    // constructors
    /**
     * @brief Construct a new NDArray object
     * 
     */
    NDArray()= delete;
    // initializer_list constructor
    NDArray(std::initializer_list<size_t> il) : shape_(il), data_(std::accumulate(shape_.begin(), shape_.end(), 1, std::multiplies<size_t>())) {}
    NDArray(std::initializer_list<int> il) : shape_(il.begin(), il.end()), data_(std::accumulate(shape_.begin(), shape_.end(), 1, std::multiplies<size_t>())) {}
    // variadic template constructor, (delegate constructor)
    template<typename... Args> NDArray(const size_t idx, Args... args) : NDArray({idx, static_cast<size_t>(args)...}) {}
    template<typename... Args> NDArray(const int& idx, Args... args) : NDArray({idx, static_cast<int>(args)...}) {} // not happy with this because size_t can have larger range
    // copy constructor
    NDArray(const NDArray& other) : data_(other.data_), shape_(other.shape_) {}
    // move constructor
    NDArray(NDArray&& other) : data_(std::move(other.data_)), shape_(std::move(other.shape_)) {}

    // destructor
    ~NDArray() {}

    // operators
    /**
     * @brief = operator, copy assignment
     * 
     * @param other 
     * @return NDArray& 
     */
    NDArray& operator=(const NDArray& other)
    {
        if (this != &other)
        {
            data_ = other.data_;
            shape_ = other.shape_;
        }
        return *this;
    }
    /**
     * @brief = operator, move assignment
    */
    NDArray& operator=(NDArray&& other)
    {
        if (this != &other)
        {
            data_ = std::move(other.data_);
            shape_ = std::move(other.shape_);
        }
        return *this;
    }

    /**
     * @brief == operator
     * 
     * @param other 
     * @return true if the data and shape are the same
     * @return false otherwise
     */
    bool operator==(const NDArray& other) const { return data_ == other.data_ && shape_ == other.shape_; }
    /**
     * @brief != operator
     * 
     * @param other 
     * @return true if the data and shape are different
     * @return false otherwise
     */
    bool operator!=(const NDArray& other) const { return !(*this == other); }
    // other operators are not generally supported

    // element access
    /**
     * @brief at function
     * 
     * @tparam Args 
     * @param args indices of the element
     * @return T& or const T&
     */
    template<typename... Args> T& at(const size_t idx, Args... args) { return data_[index(idx, args...)]; }
    template<typename... Args> const T& at(const size_t idx, Args... args) const { return data_[index(idx, args...)]; }
    /**
     * @brief [] operator
     * 
     * @tparam Args 
     * @param args indices of the element
     * @return T& or const T&
     */
    template<typename... Args> T& operator()(const size_t idx, Args... args) { return data_[index(idx, args...)]; }
    template<typename... Args> const T& operator()(const size_t idx, Args... args) const { return data_[index(idx, args...)]; }
    // front
    T& front() { return data_.front(); }
    const T& front() const { return data_.front(); }
    // back
    T& back() { return data_.back(); }
    const T& back() const { return data_.back(); }
    // data
    T* data() { return data_.data(); }
    const T* data() const { return data_.data(); }

    // iterators
    // iterators on the whole data
    T* begin() { return data_.data(); }
    T* end() { return data_.data() + data_.size(); }
    const T* cbegin() const { return data_.data(); }
    const T* cend() const { return data_.data() + data_.size(); }
    // iterators on different dimensions
    
    // capacity
    // size
    size_t size() const { return data_.size(); }
    size_t size(const size_t& dim) const { return shape_.at(dim); }
    // empty
    bool empty() const { return data_.empty(); }
    // multi-dimensional specific
    // shape
    const std::vector<size_t>& shape() const { return shape_; }
    // reshape
    template<typename... Args>
    void reshape(Args... args)
    {
        // DEVELP WARNING: what if arg = -2? :)
        // save args into a vector
        //std::vector<int64_t> dims = {static_cast<int64_t>(args)...};
        std::vector<int64_t> dims = {args...};
        // assert number of -1 in dims is at most 1
        // -1 is not type-safe!!!
        size_t count = std::count_if(dims.begin(), dims.end(), [](size_t i) { return i == -1; });
        assert(count <= 1);
        // if there is -1, calculate the size
        if (count == 1)
        {
            size_t size = 1;
            for (size_t i = 0; i < dims.size(); ++i)
            {
                if (dims[i] != -1)
                {
                    size *= dims[i];
                }
            }
            size_t idx = std::find(dims.begin(), dims.end(), -1) - dims.begin();
            dims[idx] = data_.size() / size;
        }
        // calculate the size
        size_t size = std::accumulate(dims.begin(), dims.end(), 1, std::multiplies<size_t>());
        // assert size is the same
        assert(size == data_.size());
        // assign dims to shape_
        std::copy(dims.begin(), dims.end(), shape_.begin());
    }

    // interface to ATen::Tensor, but constraint to int, double, float, std::complex<float>, std::complex<double>
    /**
     * @brief SFINAE (Substitution Failure Is Not An Error) to_tensor function, only if T is int, double, float, std::complex<float>, std::complex<double>, otherwise there is no such function
     * 
     * @return container::Tensor, only if T is int, double, float, std::complex<float>, std::complex<double>
     */
    // std::enable_if<
    //     std::is_same<T, int>::value 
    //  || std::is_same<T, double>::value 
    //  || std::is_same<T, float>::value 
    //  || std::is_same<T, std::complex<float>>::value 
    //  || std::is_same<T, std::complex<double>>::value, container::Tensor
    // >::type to_tensor() const
    // {
    //     container::TensorShape shape(shape_);
    //     container::Tensor result = container::Tensor(container::DataTypeToEnum<T>::value, shape);
    //     std::memcpy(result.data<T>(), data_.data(), data_.size() * sizeof(T));
    //     return result;
    // }
    template<typename... Args>
    size_t index(const size_t idx, Args... args) const
    {
        assert(sizeof...(args) == shape_.size() - 1); // assert the indices are the same as the shape
        size_t indices[] = {idx, static_cast<size_t>(args)...};
        size_t index = 0;
        for (size_t i = 0; i < shape_.size(); ++i)
        {
            index += indices[i] * std::accumulate(shape_.begin() + i + 1, shape_.end(), 1, std::multiplies<size_t>());
        }
        assert(index < data_.size()); // assert the index is within the data
        return index;
    }
private:
    std::vector<size_t> shape_;
    // for GPU-compatible data container, will be replaced by raw pointer
    std::vector<T> data_;
};

#endif // NDARRAY_H