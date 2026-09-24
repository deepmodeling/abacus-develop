#ifndef FORCE_RT_OVERLAP_H
#define FORCE_RT_OVERLAP_H

#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_cell/klist.h"
#include "source_lcao/allocate_dm.h"
#include "source_hamilt/hamilt.h"

template <typename T>
void cal_foverlap_rt(ModuleBase::matrix& foverlap,
                     const module_dm::Setup_DM<T>& dmat,
                     hamilt::Hamilt<T>* p_hamilt,
                     const K_Vectors& kv,
                     Parallel_Orbitals& pv,
                     UnitCell& ucell);

#endif