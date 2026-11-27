#ifndef EDM_TDDFT_H
#define EDM_TDDFT_H

#include "source_base/complexmatrix.h"
#include "source_base/matrix_block.h"
#include "source_base/parallel_blkcyc.h"

namespace module_dm
{
// Compute one k point from host matrices in the supplied block-cyclic layout.
// H, S and DMK are borrowed, read-only inputs; EDM keeps the DMK storage convention.
void edm_tddft(const Parallel_Blkcyclic& layout,
              const ModuleBase::MatrixBlock<std::complex<double>>& h_mat,
              const ModuleBase::MatrixBlock<std::complex<double>>& s_mat,
              std::complex<double>* dmk,
              ModuleBase::ComplexMatrix& edmk);

template <typename Device>
void edm_tddft_lapack(const Parallel_Blkcyclic& layout,
                     const ModuleBase::MatrixBlock<std::complex<double>>& h_mat,
                     const ModuleBase::MatrixBlock<std::complex<double>>& s_mat,
                     std::complex<double>* dmk,
                     ModuleBase::ComplexMatrix& edmk);
} // namespace module_dm
#endif // EDM_TDDFT_H
