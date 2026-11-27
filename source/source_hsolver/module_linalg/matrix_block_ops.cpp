#include "matrix_block_ops.h"

#include "source_base/kernels/math_kernel_op.h"
#include "source_base/matrix_block.h"
#include "source_base/module_external/blas_connector.h"
#ifdef __MPI
#include "source_base/module_external/scalapack_connector.h"
#endif

#include <complex>

namespace hsolver
{
template <typename T>
void multiply_matrix_blocks(char transa, char transb, int m, int n, int k,
                            const ModuleBase::MatrixBlock<T>& a,
                            const ModuleBase::MatrixBlock<T>& b,
                            const ModuleBase::MatrixBlock<T>& c)
{
#ifdef __MPI
    if (c.desc != nullptr)
    {
        const int one = 1;
        ScalapackConnector::gemm(transa, transb, m, n, k, T(1), a.p, one, one, a.desc,
                                 b.p, one, one, b.desc, T(0), c.p, one, one, c.desc);
        return;
    }
#endif
    BlasConnector::gemm_cm(transa, transb, m, n, k, T(1), a.p, a.row, b.p, b.row, T(0), c.p, c.row);
}

template <typename T>
void multiply_hermitian_upper(int m, int n,
                              const ModuleBase::MatrixBlock<T>& a,
                              const ModuleBase::MatrixBlock<T>& b,
                              const ModuleBase::MatrixBlock<T>& c)
{
#ifdef __MPI
    if (c.desc != nullptr)
    {
        ScalapackConnector::hemm('L', 'U', m, n, T(1), a.p, a.desc, b.p, b.desc, T(0), c.p, c.desc);
        return;
    }
#endif
    BlasConnector::hemm_cm('L', 'U', m, n, T(1), a.p, a.row, b.p, b.row, T(0), c.p, c.row);
}

template <typename T>
void copy_hermitian_upper(int order,
                          const ModuleBase::MatrixBlock<T>& source,
                          const ModuleBase::MatrixBlock<T>& target)
{
#ifdef __MPI
    if (target.desc != nullptr)
    {
        const int one = 1;
        ScalapackConnector::tranc(order, order, T(1), source.p, one, one, source.desc,
                                  T(0), target.p, one, one, target.desc);
        ScalapackConnector::lacpy('U', order, order, source.p, source.desc, target.p, target.desc);
        return;
    }
#endif
    for (int col = 0; col < order; ++col)
    {
        target.p[col + col * target.row] = source.p[col + col * source.row];
        for (int row = 0; row < col; ++row)
        {
            const T value = source.p[row + col * source.row];
            target.p[row + col * target.row] = value;
            target.p[col + row * target.row] = ModuleBase::get_conj(value);
        }
    }
}

template void multiply_matrix_blocks(char, char, int, int, int,
                                     const ModuleBase::MatrixBlock<double>&,
                                     const ModuleBase::MatrixBlock<double>&,
                                     const ModuleBase::MatrixBlock<double>&);
template void multiply_matrix_blocks(char, char, int, int, int,
                                     const ModuleBase::MatrixBlock<std::complex<double>>&,
                                     const ModuleBase::MatrixBlock<std::complex<double>>&,
                                     const ModuleBase::MatrixBlock<std::complex<double>>&);

template void multiply_hermitian_upper(int, int,
                                       const ModuleBase::MatrixBlock<double>&,
                                       const ModuleBase::MatrixBlock<double>&,
                                       const ModuleBase::MatrixBlock<double>&);
template void multiply_hermitian_upper(int, int,
                                       const ModuleBase::MatrixBlock<std::complex<double>>&,
                                       const ModuleBase::MatrixBlock<std::complex<double>>&,
                                       const ModuleBase::MatrixBlock<std::complex<double>>&);

template void copy_hermitian_upper(int, const ModuleBase::MatrixBlock<double>&, const ModuleBase::MatrixBlock<double>&);
template void copy_hermitian_upper(int,
                                  const ModuleBase::MatrixBlock<std::complex<double>>&,
                                  const ModuleBase::MatrixBlock<std::complex<double>>&);
} // namespace hsolver
