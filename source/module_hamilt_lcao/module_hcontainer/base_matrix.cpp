#include "base_matrix.h"
#include <complex>
#include <vector>

#include "module_base/global_function.h"

namespace hamilt
{

template <typename T>
BaseMatrix<T>::BaseMatrix(const int& nrow_, const int& ncol_, T* data_existed)
{
    nrow_local = nrow_;
    ncol_local = ncol_;
    this->allocated = false;
    value_begin = data_existed;
}

// move constructor
template <typename T>
BaseMatrix<T>::BaseMatrix(BaseMatrix<T>&& matrix)
{
    this->nrow_local = matrix.nrow_local;
    this->ncol_local = matrix.ncol_local;
    this->value_begin = matrix.value_begin;
    this->allocated = matrix.allocated;
    if (matrix.allocated)
    {
        matrix.allocated = false;
    }
    matrix.value_begin = nullptr;
}

// copy constructor
template <typename T>
BaseMatrix<T>::BaseMatrix(const BaseMatrix<T>& matrix)
{
    this->nrow_local = matrix.nrow_local;
    this->ncol_local = matrix.ncol_local;
    if (matrix.allocated)
    {
        this->value_begin = new T[nrow_local * ncol_local];
        ModuleBase::GlobalFunc::ZEROS(this->value_begin, nrow_local * ncol_local);
        this->allocated = true;
        for (int i = 0; i < nrow_local * ncol_local; i++)
        {
            this->value_begin[i] = matrix.value_begin[i];
        }
    }
    else
    {
        this->value_begin = matrix.value_begin;
        this->allocated = false;
    }
}

template <typename T>
BaseMatrix<T>::~BaseMatrix()
{
    if (this->allocated)
    {
        delete[] value_begin;
    }
}

// allocate
template <typename T>
void BaseMatrix<T>::allocate(T* data_array, bool if_zero)
{
#ifdef __DEBUG
assert(nrow_local*ncol_local>0);
#endif
    if(data_array != nullptr && !this->allocated)
    {
        this->value_begin = data_array;
    }
    else if(data_array != nullptr && this->allocated)
    {
        delete[] this->value_begin;
        this->value_begin = data_array;
        this->allocated = false;
    }
    else if(data_array == nullptr && !this->allocated)
    {
        this->value_begin = new T[nrow_local * ncol_local];
        this->allocated = true;
    }
    else
    {
        // do nothing
    }
    if(if_zero) 
    {
        this->set_zero();
    }
}

// zeros
template <typename T>
void BaseMatrix<T>::set_zero()
{
#ifdef __DEBUG
assert(this->value_begin != nullptr);
#endif
    ModuleBase::GlobalFunc::ZEROS(this->value_begin, nrow_local * ncol_local);
}

// add_array
template <typename T>
void BaseMatrix<T>::add_array(T* array)
{
    // if allocated, save data from array into matrix
    // if whole matrix and 2d-block format, save data from array into matrix either
    for (int i = 0; i < nrow_local * ncol_local; ++i)
    {
        value_begin[i] += array[i];
    }
}

template <typename T>
void BaseMatrix<T>::add_element(int mu, int nu, const T& value)
{
    int index = mu * this->ncol_local + nu;
    value_begin[index] += value;
}

template <typename T>
T& BaseMatrix<T>::get_value(const size_t& i_row, const size_t& j_col) const
{
    int index = i_row * this->ncol_local + j_col;
    return value_begin[index];
}

template <typename T>
T* BaseMatrix<T>::get_pointer() const
{
    return value_begin;
}

// operator= for copy assignment
template <typename T>
BaseMatrix<T>& BaseMatrix<T>::operator=(const BaseMatrix<T>& other)
{
    if (this != &other)
    {
        this->nrow_local = other.nrow_local;
        this->ncol_local = other.ncol_local;
        if (other.allocated)
        {
            this->value_begin = new T[nrow_local * ncol_local];
            ModuleBase::GlobalFunc::ZEROS(this->value_begin, nrow_local * ncol_local);
            this->allocated = true;
            for (int i = 0; i < nrow_local * ncol_local; ++i)
            {
                this->value_begin[i] = other.value_begin[i];
            }
        }
        else
        {
            this->value_begin = other.value_begin;
            this->allocated = false;
        }
    }
    return *this;
}

// operator= for move assignment
template <typename T>
BaseMatrix<T>& BaseMatrix<T>::operator=(BaseMatrix<T>&& other) noexcept
{
    if (this != &other)
    {
        this->nrow_local = other.nrow_local;
        this->ncol_local = other.ncol_local;
        this->value_begin = other.value_begin;
        this->allocated = other.allocated;
        if (other.allocated)
        {
            other.allocated = false;
        }
        other.value_begin = nullptr;
    }
    return *this;
}

// get_memory_size
template <typename T>
size_t BaseMatrix<T>::get_memory_size() const
{
    size_t memory_size = sizeof(*this);
    if(this->allocated)
    {
        memory_size += nrow_local * ncol_local * sizeof(T);
    }
    return memory_size;
}

// T of BaseMatrix can be double or complex<double>
template class BaseMatrix<double>;
template class BaseMatrix<std::complex<double>>;

} // namespace hamilt