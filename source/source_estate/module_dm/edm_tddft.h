#ifndef EDM_TDDFT_H
#define EDM_TDDFT_H

#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_cell/klist.h"
#include "source_estate/module_dm/dm_holder.h"
#include "source_hamilt/hamilt.h"

namespace module_dm
{
void edm_tddft(Parallel_Orbitals& pv,
                   Setup_DM<std::complex<double>>& dmat,
                   K_Vectors& kv,
                   hamilt::Hamilt<std::complex<double>>* p_hamilt);

template <typename Device>
void edm_tddft_lapack(Parallel_Orbitals& pv,
                                 Setup_DM<std::complex<double>>& dmat,
                                 K_Vectors& kv,
                                 hamilt::Hamilt<std::complex<double>>* p_hamilt);
} // namespace module_dm
#endif // EDM_TDDFT_H
