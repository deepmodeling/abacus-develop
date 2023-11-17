#ifndef BASE_MATRIX_H
#define BASE_MATRIX_H

#include <iostream>

namespace hamilt
{
/**
 * class: BaseMatrix
 * used to store a matrix for atom-pair local Hamiltonian with specific R-index
 * T can be double or complex<double>
 * It has two ways to arrange data:
 * 1. allocate data itself
 * 2. only access data but be arranged by other class
 */
template <typename T>
class BaseMatrix
{
  public:
    // Constructor of class BaseMatrix
    BaseMatrix(const int& nrow_, const int& ncol_, T* data_existed = nullptr);
    // copy constructor
    BaseMatrix(const BaseMatrix<T>& matrix);
    // move constructor
    BaseMatrix(BaseMatrix<T>&& matrix);
    // Destructor of class BaseMatrix
    ~BaseMatrix();

    /**
     * @brief allocate memory for the matrix
     * if this->value_begin is not nullptr, it will be neglected
     * if this->value_begin is nullptr, it will allocate memory with size nrow_local * ncol_local
    */
    void allocate(T* data_array = nullptr, bool if_zero = false);

    /**
     * @brief set value in the matrix to zero
    */
    void set_zero();

    /**
     * @brief save an array to the matrix
     *
     * @param array array to be saved
     */
    void add_array(T* array);
    /**
     * @brief add a single element to the matrix
     *
     * @param mu row index
     * @param nu column index
     * @param value value to be added
     */
    void add_element(int mu, int nu, const T& value);
    // for inside matrix
    /**
     * @brief get value from a whole matrix
     *
     * @param i_row row index
     * @param j_col column index
     * @return T&
     */
    T& get_value(const size_t& i_row, const size_t& j_col) const;
    /**
     * @brief get pointer of value from a submatrix
     */
    T* get_pointer() const;

    // operator= for copy assignment
    BaseMatrix& operator=(const BaseMatrix& other);

    // operator= for move assignment
    BaseMatrix& operator=(BaseMatrix&& other) noexcept;

    /**
     * @brief get total memory size of BaseMatrix
    */
    size_t get_memory_size() const;

  private:
    bool allocated = false;

    // pointer for accessing data
    // two ways to arrange data:
    // 1. allocate data itself
    // 2. only access data but be arranged by RealSparseHamiltonian
    T* value_begin = nullptr;

    // int current_multiple = 0;

    // number of rows and columns
    int nrow_local = 0;
    int ncol_local = 0;
};

} // namespace hamilt

#endif