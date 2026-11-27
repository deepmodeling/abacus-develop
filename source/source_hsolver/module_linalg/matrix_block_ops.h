#ifndef HSOLVER_MATRIX_BLOCK_OPS_H
#define HSOLVER_MATRIX_BLOCK_OPS_H

namespace ModuleBase
{
template <typename T> struct MatrixBlock;
}

namespace hsolver
{
// Operands use column-major host storage. Null descriptors select local BLAS.
// Distributed operands must satisfy the corresponding PBLAS alignment rules.

// C = op(A) op(B) for the leading subblocks of the given dimensions.
template <typename T>
void multiply_matrix_blocks(char transa, char transb, int m, int n, int k,
                            const ModuleBase::MatrixBlock<T>& a,
                            const ModuleBase::MatrixBlock<T>& b,
                            const ModuleBase::MatrixBlock<T>& c);

// C = A B, where A is m by m Hermitian and only its upper triangle is read.
template <typename T>
void multiply_hermitian_upper(int m, int n,
                              const ModuleBase::MatrixBlock<T>& a,
                              const ModuleBase::MatrixBlock<T>& b,
                              const ModuleBase::MatrixBlock<T>& c);

// Reconstruct a Hermitian matrix from the source's upper triangle.
// Source and destination must be distinct and have the same block-cyclic layout.
template <typename T>
void copy_hermitian_upper(int order,
                          const ModuleBase::MatrixBlock<T>& source,
                          const ModuleBase::MatrixBlock<T>& target);
} // namespace hsolver

#endif
