#ifndef CAL_EDM_TDDFT_H
#define CAL_EDM_TDDFT_H

#include "module_basis/module_ao/parallel_orbitals.h"
#include "module_cell/klist.h"
#include "module_elecstate/elecstate_lcao.h"
#include "module_hamilt_general/hamilt.h"

namespace elecstate
{
void cal_edm_tddft(Parallel_Orbitals& pv,
                   elecstate::ElecState* pelec,
                   K_Vectors& kv,
                   hamilt::Hamilt<std::complex<double>>* p_hamilt);
} // namespace elecstate
#endif // CAL_EDM_TDDFT_H
