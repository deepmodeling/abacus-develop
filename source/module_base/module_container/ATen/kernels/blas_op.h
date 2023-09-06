#ifndef ATEN_KERNELS_BLAS_OP_H_
#define ATEN_KERNELS_BLAS_OP_H_

#include <ATen/core/tensor.h>
#include <ATen/core/tensor_types.h>

#include <base/third_party/blas.h>

namespace container {
namespace op {

/**
 * @brief A functor representing the dot product operation using BLAS(CPU)/cuBLAS(GPU) library.
 *
 * This functor defines the dot product operation between two vectors x and y.
 * The dot product of two vectors x and y is the sum of the element-wise products
 * of their corresponding elements.
 *
 * @tparam T The type of the elements in the vectors (e.g., float, double, etc.).
 * @tparam Device The device where the vectors reside (e.g., CPU, GPU, etc.).
 */
template <typename T, typename Device>
struct blas_dot {
    /**
     * @brief Compute the dot product of two vectors.
     *
     * This function calculates the dot product of the input vectors x and y.
     * The dot product is the sum of the element-wise products of their corresponding elements.
     *
     * @param n The number of elements in the vectors x and y.
     * @param x A pointer to the first input vector x.
     * @param incx The increment between consecutive elements of vector x.
     * @param y A pointer to the second input vector y.
     * @param incy The increment between consecutive elements of vector y.
     * @param result A pointer to store the result of the dot product.
     *
     * @note The lengths of the x and y arrays should be at least n * incx and n * incy, respectively.
     * @note The dot product is calculated as: result = sum(x[i] * y[i]), for i = 0, 1, ..., n-1.
     */
    void operator()(
        const int& n,
        const T* x,
        const int& incx,
        const T* y,
        const int& incy,
        T* result);
};


/**
 * @brief A functor representing a scaling operation on a vector using BLAS(CPU)/cuBLAS(GPU) library.
 *
 * This struct defines the operation of scaling a vector with a scalar.
 * The scaling operation modifies the input vector by multiplying each element
 * of the vector by the given scalar alpha.
 *
 * @tparam T The type of the elements in the vector (e.g., float, double, etc.).
 * @tparam Device The device where the vector resides (e.g., CPU, GPU, etc.).
 */
template <typename T, typename Device>
struct blas_scal {
    /**
     * @brief Perform the scaling operation on the input vector.
     *
     * This function applies the scaling operation to the input vector X.
     * It multiplies each element of the vector by the given scalar alpha.
     *
     * @param n The number of elements in the vector.
     * @param alpha A pointer to the scalar used for scaling.
     * @param x A pointer to the input vector to be scaled.
     * @param incx The increment between consecutive elements of the vector X.
     *
     * @note The length of the vector X should be at least (1 + (N - 1) * abs(incx)).
     */
    void operator()(
        const int& n,
        const T* alpha,
        T* x,
        const int& incx);
};


/**
 * @brief A functor representing the axpy operation on vectors using BLAS(CPU)/cuBLAS(GPU) library.
 *
 * This strfunctor defines the axpy (Y = alpha * X + Y) operation on two vectors X and Y,
 * where alpha is a scalar. The axpy operation modifies the elements of the output vector Y
 * by adding the scaled elements of the input vector X to it.
 *
 * @tparam T The type of the elements in the vectors (e.g., float, double, etc.).
 * @tparam Device The device where the vectors reside (e.g., CPU, GPU, etc.).
 */
template <typename T, typename Device>
struct blas_axpy {
    /**
     * @brief Perform the axpy operation on the input vectors.
     *
     * This function applies the axpy operation to the input vectors X and Y.
     * The axpy operation modifies the elements of vector Y by adding the scaled
     * elements of vector X to it. The scaling factor alpha is a scalar.
     *
     * @param n The number of elements in the vectors.
     * @param alpha A pointer to the scalar used for scaling vector X.
     * @param x A pointer to the input vector X to be scaled and added to vector Y.
     * @param incx The increment between consecutive elements of vector X.
     * @param y A pointer to the output vector Y to be updated with the axpy result.
     * @param incy The increment between consecutive elements of vector Y.
     *
     * @note The length of vectors X and Y should be at least (1 + (N - 1) * abs(incX)) and
     *       (1 + (N - 1) * abs(incY)), respectively.
     */
    void operator()(
        const int& n,
        const T* alpha,
        const T* x,
        const int& incx,
        T* y,
        const int& incy);
};


/**
 * @brief A functor representing the matrix-vector multiplication operation using BLAS(CPU)/cuBLAS(GPU) library.
 *
 * This functor defines the general matrix-vector multiplication operation:
 * Y = alpha * op(A) * X + beta * Y, where A is a matrix, X and Y are vectors,
 * alpha and beta are scalars, and op(A) depends on the value of 'trans'.
 *
 * @tparam T The type of the elements in the matrix and vectors (e.g., float, double, etc.).
 * @tparam Device The device where the matrix and vectors reside (e.g., CPU, GPU, etc.).
 */
template <typename T, typename Device>
struct blas_gemv {
    /**
     * @brief Perform the general matrix-vector multiplication operation.
     *
     * This function applies the general matrix-vector multiplication operation
     * to the input matrix A and vectors X and Y. The result is stored in vector Y.
     *
     * @param trans Specifies whether to use the transpose or conjugate transpose of matrix A.
     *              - 'N' or 'n': No transpose (op(A) = A).
     *              - 'T' or 't': Transpose (op(A) = transpose of A).
     *              - 'C' or 'c': Conjugate transpose (op(A) = conjugate transpose of A).
     * @param m The number of rows in matrix A.
     * @param n The number of columns in matrix A.
     * @param alpha A pointer to the scalar used to scale the matrix-vector product.
     * @param A A pointer to the input matrix A.
     * @param lda The leading dimension (stride) of matrix A (the number of elements between
     *            consecutive columns for row-major storage or consecutive rows for column-major storage).
     * @param x A pointer to the input vector X.
     * @param incx The increment between consecutive elements of vector X.
     * @param beta A pointer to the scalar used to scale vector Y before adding the matrix-vector product.
     * @param y A pointer to the output vector Y to store the result.
     * @param incy The increment between consecutive elements of vector Y.
     *
     * @note The length of vector X should be at least (1 + (n - 1) * abs(incx)).
     * @note The length of vector Y should be at least (1 + (m - 1) * abs(incy)).
     * @note The lda parameter should be at least max(1, m) for row-major storage and max(1, n) for column-major storage.
     */
    void operator()(
        const char& trans,
        const int& m,
        const int& n,
        const T* alpha,
        const T* A,
        const int& lda,
        const T* x,
        const int& incx,
        const T* beta,
        T* y,
        const int& incy);
};


/**
 * @brief A functor representing the batched matrix-vector multiplication operation using BLAS(CPU)/cuBLAS(GPU) library.
 *
 * This functor defines the batched matrix-vector multiplication operation:
 * y[i] = alpha * op(A[i]) * x[i] + beta * y[i], for i = 0, 1, ..., batch_size-1,
 * where A[i] is a matrix, x[i] and y[i] are vectors, alpha and beta are scalars,
 * and op(A[i]) depends on the value of 'trans'.
 *
 * @tparam T The type of the elements in the matrix and vectors (e.g., float, double, etc.).
 * @tparam Device The device where the matrix and vectors reside (e.g., CPU, GPU, etc.).
 */
template <typename T, typename Device>
struct blas_gemv_batched {
    /**
     * @brief Perform the batched matrix-vector multiplication operation.
     *
     * This function applies the batched matrix-vector multiplication operation
     * to the input batches of matrices A and vectors x. The result is stored in the batches of vectors y.
     *
     * @param trans Specifies whether to use the transpose or conjugate transpose of matrices A[i].
     *              - 'N' or 'n': No transpose (op(A[i]) = A[i]).
     *              - 'T' or 't': Transpose (op(A[i]) = transpose of A[i]).
     *              - 'C' or 'c': Conjugate transpose (op(A[i]) = conjugate transpose of A[i]).
     * @param m The number of rows in each matrix A[i].
     * @param n The number of columns in each matrix A[i].
     * @param alpha A pointer to the scalar used to scale the matrix-vector products.
     * @param A Array of pointers to the input matrices A[i].
     * @param lda The leading dimension (stride) of each matrix A[i].
     * @param x Array of pointers to the input vectors x[i].
     * @param incx The increment between consecutive elements of vectors x[i].
     * @param beta A pointer to the scalar used to scale each vector y[i] before adding the matrix-vector product.
     * @param y Array of pointers to the output vectors y[i] to store the results.
     * @param incy The increment between consecutive elements of vectors y[i].
     * @param batch_size The number of batches.
     *
     * @note The leading dimensions lda should be at least max(1, m) for row-major storage
     *       and max(1, n) for column-major storage.
     * @note The length of the A, x, and y arrays should be at least batch_size.
     */
    void operator()(
        const char& trans,
        const int& m,
        const int& n,
        const T* alpha,
        T** A,
        const int& lda,
        T** x,
        const int& incx,
        const T* beta,
        T** y,
        const int& incy,
        const int& batch_size);
};


/**
 * @brief A functor representing the strided batched matrix-vector multiplication operation using BLAS(CPU)/cuBLAS(GPU) library.
 *
 * This functor defines the strided batched matrix-vector multiplication operation:
 * y[i] = alpha * op(A[i]) * x[i] + beta * y[i], for i = 0, 1, ..., batch_size-1,
 * where A[i] is a matrix with a strided memory layout, x[i] and y[i] are vectors with
 * strided memory layouts, alpha and beta are scalars, and op(A[i]) depends on the value of 'trans'.
 *
 * @tparam T The type of the elements in the matrix and vectors (e.g., float, double, etc.).
 * @tparam Device The device where the matrix and vectors reside (e.g., CPU, GPU, etc.).
 */
template <typename T, typename Device>
struct blas_gemv_batched_strided {
    /**
     * @brief Perform the strided batched matrix-vector multiplication operation.
     *
     * This function applies the strided batched matrix-vector multiplication operation
     * to the input batches of matrices A and vectors x. The result is stored in the batches of vectors y.
     * The strided memory layout allows each matrix and vector to have a different stride, allowing for efficient
     * handling of non-contiguous data.
     *
     * @param trans Specifies whether to use the transpose or conjugate transpose of matrices A[i].
     *              - 'N' or 'n': No transpose (op(A[i]) = A[i]).
     *              - 'T' or 't': Transpose (op(A[i]) = transpose of A[i]).
     *              - 'C' or 'c': Conjugate transpose (op(A[i]) = conjugate transpose of A[i]).
     * @param m The number of rows in each matrix A[i].
     * @param n The number of columns in each matrix A[i].
     * @param alpha A pointer to the scalar used to scale the matrix-vector products.
     * @param A A pointer to the input batched matrices A with a strided memory layout.
     * @param lda The leading dimension (stride) of each matrix A[i].
     * @param stride_a The stride between the starting elements of consecutive matrices in A.
     * @param x A pointer to the input batched vectors x with a strided memory layout.
     * @param incx The increment between consecutive elements of vectors x[i].
     * @param stride_x The stride between the starting elements of consecutive vectors in x.
     * @param beta A pointer to the scalar used to scale each vector y[i] before adding the matrix-vector product.
     * @param y A pointer to the output batched vectors y to store the results with a strided memory layout.
     * @param incy The increment between consecutive elements of vectors y[i].
     * @param stride_y The stride between the starting elements of consecutive vectors in y.
     * @param batch_size The number of batches.
     *
     * @note The leading dimension lda should be at least max(1, m) for row-major storage
     *       and max(1, n) for column-major storage.
     * @note The strides stride_a, stride_x, and stride_y are in units of elements (not bytes).
     * @note The lengths of the A, x, and y arrays should be at least batch_size * lda, batch_size * incx,
     *       and batch_size * incy, respectively.
     */
    void operator()(
        const char& trans,
        const int& m,
        const int& n,
        const T* alpha,
        const T* A,
        const int& lda,
        const int64_t& stride_a,
        const T* x,
        const int& incx,
        const int64_t& stride_x,
        const T* beta,
        T* y,
        const int& incy,
        const int64_t& stride_y,
        const int& batch_size);
};


/**
 * @brief A functor representing the matrix-matrix multiplication operation using BLAS(CPU)/cuBLAS(GPU) library.
 *
 * This functor defines the general matrix-matrix multiplication operation:
 * C = alpha * op(A) * op(B) + beta * C, where A and B are matrices, C is a matrix,
 * alpha and beta are scalars, and op(A) and op(B) depend on the values of 'transa' and 'transb'.
 *
 * @tparam T The type of the elements in the matrices (e.g., float, double, etc.).
 * @tparam Device The device where the matrices reside (e.g., CPU, GPU, etc.).
 */
template <typename T, typename Device>
struct blas_gemm {
    /**
     * @brief Perform the general matrix-matrix multiplication operation.
     *
     * This function applies the general matrix-matrix multiplication operation
     * to the input matrices A and B. The result is stored in matrix C.
     *
     * @param transa Specifies whether to use the transpose or conjugate transpose of matrix A.
     *               - 'N' or 'n': No transpose (op(A) = A).
     *               - 'T' or 't': Transpose (op(A) = transpose of A).
     *               - 'C' or 'c': Conjugate transpose (op(A) = conjugate transpose of A).
     * @param transb Specifies whether to use the transpose or conjugate transpose of matrix B.
     *               - 'N' or 'n': No transpose (op(B) = B).
     *               - 'T' or 't': Transpose (op(B) = transpose of B).
     *               - 'C' or 'c': Conjugate transpose (op(B) = conjugate transpose of B).
     * @param m The number of rows in matrix C.
     * @param n The number of columns in matrix C.
     * @param k The number of columns in matrix A (number of rows in matrix B).
     * @param alpha A pointer to the scalar used to scale the matrix-matrix product.
     * @param A A pointer to the input matrix A.
     * @param lda The leading dimension (stride) of matrix A (the number of elements between
     *            consecutive columns for row-major storage or consecutive rows for column-major storage).
     * @param B A pointer to the input matrix B.
     * @param ldb The leading dimension (stride) of matrix B (the number of elements between
     *            consecutive columns for row-major storage or consecutive rows for column-major storage).
     * @param beta A pointer to the scalar used to scale matrix C before adding the matrix-matrix product.
     * @param C A pointer to the output matrix C to store the result.
     * @param ldc The leading dimension (stride) of matrix C (the number of elements between
     *            consecutive columns for row-major storage or consecutive rows for column-major storage).
     *
     * @note The leading dimensions lda, ldb, and ldc should be at least max(1, m) for row-major storage
     *       and max(1, k) for column-major storage.
     */
    void operator()(
        const char& transa,
        const char& transb,
        const int& m,
        const int& n,
        const int& k,
        const T* alpha,
        const T* A,
        const int& lda,
        const T* B,
        const int& ldb,
        const T* beta,
        T* C,
        const int& ldc);
};

/**
 * @brief A functor representing the batched matrix-matrix multiplication operation using BLAS(CPU)/cuBLAS(GPU) library.
 *
 * This functor defines the batched matrix-matrix multiplication operation:
 * C[i] = alpha * op(A[i]) * op(B[i]) + beta * C[i], for i = 0, 1, ..., batch_size-1,
 * where A[i] and B[i] are matrices, C[i] is a matrix, alpha and beta are scalars,
 * and op(A[i]) and op(B[i]) depend on the values of 'transa' and 'transb'.
 *
 * @tparam T The type of the elements in the matrices (e.g., float, double, etc.).
 * @tparam Device The device where the matrices reside (e.g., CPU, GPU, etc.).
 */
template <typename T, typename Device>
struct blas_gemm_batched {
    /**
     * @brief Perform the batched matrix-matrix multiplication operation.
     *
     * This function applies the batched matrix-matrix multiplication operation
     * to the input batches of matrices A and B. The result is stored in the batches of matrices C.
     *
     * @param transa Specifies whether to use the transpose or conjugate transpose of matrices A[i].
     *               - 'N' or 'n': No transpose (op(A[i]) = A[i]).
     *               - 'T' or 't': Transpose (op(A[i]) = transpose of A[i]).
     *               - 'C' or 'c': Conjugate transpose (op(A[i]) = conjugate transpose of A[i]).
     * @param transb Specifies whether to use the transpose or conjugate transpose of matrices B[i].
     *               - 'N' or 'n': No transpose (op(B[i]) = B[i]).
     *               - 'T' or 't': Transpose (op(B[i]) = transpose of B[i]).
     *               - 'C' or 'c': Conjugate transpose (op(B[i]) = conjugate transpose of B[i]).
     * @param m The number of rows in each matrix C[i].
     * @param n The number of columns in each matrix C[i].
     * @param k The number of columns in each matrix A[i] (number of rows in each matrix B[i]).
     * @param alpha A pointer to the scalar used to scale the matrix-matrix products.
     * @param A Array of pointers to the input matrices A[i].
     * @param lda The leading dimension (stride) of each matrix A[i].
     * @param B Array of pointers to the input matrices B[i].
     * @param ldb The leading dimension (stride) of each matrix B[i].
     * @param beta A pointer to the scalar used to scale each matrix C[i] before adding the matrix-matrix product.
     * @param C Array of pointers to the output matrices C[i] to store the results.
     * @param ldc The leading dimension (stride) of each matrix C[i].
     * @param batch_size The number of batches.
     *
     * @note The leading dimensions lda, ldb, and ldc should be at least max(1, m) for row-major storage
     *       and max(1, k) for column-major storage.
     * @note The length of the A, B, and C arrays should be at least batch_size.
     */
    void operator()(
        const char& transa,
        const char& transb,
        const int& m,
        const int& n,
        const int& k,
        const T* alpha,
        T** A,
        const int& lda,
        T** B,
        const int& ldb,
        const T* beta,
        T** C,
        const int& ldc,
        const int& batch_size);
};


/**
 * @brief A functor representing the strided batched matrix-matrix multiplication operation using BLAS(CPU)/cuBLAS(GPU) library.
 *
 * This functor defines the strided batched matrix-matrix multiplication operation:
 * C[i] = alpha * op(A[i]) * op(B[i]) + beta * C[i], for i = 0, 1, ..., batch_size-1,
 * where A[i] and B[i] are matrices with strided memory layout, C[i] is a matrix with
 * a strided memory layout, alpha and beta are scalars, and op(A[i]) and op(B[i]) depend on
 * the values of 'transa' and 'transb'.
 *
 * @tparam T The type of the elements in the matrices (e.g., float, double, etc.).
 * @tparam Device The device where the matrices reside (e.g., CPU, GPU, etc.).
 */
template <typename T, typename Device>
struct blas_gemm_batched_strided {
    /**
     * @brief Perform the strided batched matrix-matrix multiplication operation.
     *
     * This function applies the strided batched matrix-matrix multiplication operation
     * to the input batches of matrices A and B. The result is stored in the batches of matrices C.
     * The strided memory layout allows each matrix to have a different stride, allowing for efficient
     * handling of non-contiguous data.
     *
     * @param transa Specifies whether to use the transpose or conjugate transpose of matrices A[i].
     *               - 'N' or 'n': No transpose (op(A[i]) = A[i]).
     *               - 'T' or 't': Transpose (op(A[i]) = transpose of A[i]).
     *               - 'C' or 'c': Conjugate transpose (op(A[i]) = conjugate transpose of A[i]).
     * @param transb Specifies whether to use the transpose or conjugate transpose of matrices B[i].
     *               - 'N' or 'n': No transpose (op(B[i]) = B[i]).
     *               - 'T' or 't': Transpose (op(B[i]) = transpose of B[i]).
     *               - 'C' or 'c': Conjugate transpose (op(B[i]) = conjugate transpose of B[i]).
     * @param m The number of rows in each matrix C[i].
     * @param n The number of columns in each matrix C[i].
     * @param k The number of columns in each matrix A[i] (number of rows in each matrix B[i]).
     * @param alpha A pointer to the scalar used to scale the matrix-matrix products.
     * @param A A pointer to the input batched matrices A.
     * @param lda The leading dimension (stride) of each matrix A[i].
     * @param stride_a The stride between the starting elements of consecutive matrices in A.
     * @param B A pointer to the input batched matrices B.
     * @param ldb The leading dimension (stride) of each matrix B[i].
     * @param stride_b The stride between the starting elements of consecutive matrices in B.
     * @param beta A pointer to the scalar used to scale each matrix C[i] before adding the matrix-matrix product.
     * @param C A pointer to the output batched matrices C to store the results.
     * @param ldc The leading dimension (stride) of each matrix C[i].
     * @param stride_c The stride between the starting elements of consecutive matrices in C.
     * @param batch_size The number of batches.
     *
     * @note The leading dimensions lda, ldb, and ldc should be at least max(1, m) for row-major storage
     *       and max(1, k) for column-major storage.
     * @note The length of the A, B, and C arrays should be at least batch_size * stride_a, batch_size * stride_b,
     *       and batch_size * stride_c, respectively.
     * @note The strides stride_a, stride_b, and stride_c are in units of elements (not bytes).
     */
    void operator()(
        const char& transa,
        const char& transb,
        const int& m,
        const int& n,
        const int& k,
        const T* alpha,
        const T* A,
        const int& lda,
        const int& stride_a,
        const T* B,
        const int& ldb,
        const int& stride_b,
        const T* beta,
        T* C,
        const int& ldc,
        const int& stride_c,
        const int& batch_size);
};

#if __CUDA || __UT_USE_CUDA || __ROCM || __UT_USE_ROCM
void createBlasHandle();  // create blas handle
void destroyBlasHandle(); // destory blas handle
#endif // __CUDA || __UT_USE_CUDA

} // namespace op
} // namespace container

#endif // ATEN_KERNELS_BLAS_OP_H_