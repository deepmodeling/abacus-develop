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
 * It has four ways to access data:
 * 1. whole matrix
 * 2. 2d-block
 * 3. submatrix in whole matrix
 * 4. sparse matrix in whole matrix , not implemented yet
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
    void allocate(bool if_zero = false);

    /**
     * @brief set value in the matrix to zero
     * if memory_type == 1 , it will set all value in the matrix to zero
     * if memory_type == 2 , it will set all value in the submatrix to zero
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
     * for memory_type = 0 or 1, ncol_local will be used to calculate the index
     * for memory_type = 2, ldc will be used to calculate the index
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

    void set_ldc(const int& ldc_in);

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

    // memory type, choose how to access value via pointer of array
    // 1 is dense matrix
    // 2 is submatrix in whole matrix
    // 3 is sparse matrix in whole matrix , not implemented yet
    int memory_type = 1;

    // leading dimension of matrix, used with memory_type = 2
    int ldc = 0;
};

} // namespace hamilt

#endif