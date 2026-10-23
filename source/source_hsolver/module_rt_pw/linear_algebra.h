#ifndef HSOLVER_LINEAR_ALGEBRA_H
#define HSOLVER_LINEAR_ALGEBRA_H

#include "source_base/macros.h"
#include "source_base/module_container/ATen/core/tensor.h"
#include "source_base/module_device/memory_op.h"
#include "source_hsolver/diag_comm_info.h"

#include <algorithm>
#include <complex>
#include <cstdint>
#include <vector>

namespace hsolver
{
/** @brief Device allocation helper, including type and empty-partition checks. */
template <typename T, typename Device>
void linear_buffer(ct::Tensor* buffer, const int64_t size)
{
    using CtDevice = typename ct::PsiToContainer<Device>::type;
    const ct::DeviceType device = ct::DeviceTypeToEnum<CtDevice>::value;
    const int64_t elements = std::max<int64_t>(1, size);
    if (buffer->NumElements() < elements || buffer->data_type() != ct::DataTypeToEnum<T>::value || buffer->device_type() != device)
    {
        *buffer = ct::Tensor(ct::DataTypeToEnum<T>::value, device, {elements});
    }
}

/** @brief Partial-pivoted factorization of a small column-major complex matrix. */
class LinearSmallLU
{
  private:
    int size_ = 0;
    std::vector<std::complex<double>> lu_;
    std::vector<int> pivots_;

  public:
    bool factor(const std::vector<std::complex<double>>& matrix, int n);
    bool solve(std::vector<std::complex<double>>* rhs, int columns) const;
};

/** @brief Rank-revealing coefficients for a small Gram matrix, with a norm cutoff. */
std::vector<std::complex<double>> linear_gram_basis(const std::vector<std::complex<double>>& gram, int n, double cutoff, int* rank);

/** @brief Collective block algebra; only small coefficients cross the device boundary. */
template <typename T, typename Device>
class LinearAlgebra
{
  private:
    using Wide = std::complex<double>;
    const diag_comm_info comm_;
    ct::Tensor products_;
    ct::Tensor coefficients_;
    ct::Tensor left_;
    ct::Tensor right_;
    ct::Tensor native_products_;

    void reduce(std::vector<Wide>* values) const;

  public:
    explicit LinearAlgebra(const diag_comm_info& comm) : comm_(comm)
    {
    }

    /** @brief Projection products and reductions in wavefunction precision; return double coefficients. */
    std::vector<Wide> projection_cross(int ld, int dim, int nx, int ny, const T* x, const T* y);

    /** @brief Arnoldi products in wavefunction precision; use dots() for reliable norms. */
    std::vector<Wide> arnoldi_dots(int ld, int dim, int bands, int count, int stride, const T* basis, const T* x);

    /** @brief Products between corresponding bands and several Arnoldi blocks. */
    std::vector<Wide> dots(int ld, int dim, int bands, int count, int stride, const T* basis, const T* x);

    /** @brief Global X-adjoint times Y, with double accumulation for both precisions. */
    std::vector<Wide> cross(int ld, int dim, int nx, int ny, const T* x, const T* y);

    /** @brief Global X-adjoint times X, reusing one FP64 conversion of the valid rows. */
    std::vector<Wide> gram(int ld, int dim, int bands, const T* input);

    /** @brief Y = X*C + beta*Y, using a small host coefficient matrix. */
    void expand(int ld, int dim, int nx, int ny, const T* x, const std::vector<Wide>& c, T* y, T beta);
};
} // namespace hsolver
#endif
