#ifndef FORCE_RT_OVERLAP_H
#define FORCE_RT_OVERLAP_H

#include "source_cell/klist.h"
#include "source_hamilt/hamilt.h"
#include "source_lcao/allocate_dm.h"
#include "source_lcao/module_rt/td_info.h"

class BasisFlattenedIndexing;
class BlockCyclicDistribution;

template <typename T>
void cal_foverlap_rt(ModuleBase::matrix& foverlap,
                     const module_dm::Setup_DM<T>& dmat,
                     hamilt::Hamilt<T>* p_hamilt,
                     const K_Vectors& kv,
                     const BasisFlattenedIndexing& basis,
                     const BlockCyclicDistribution& layout,
                     UnitCell& ucell,
                     TD_info& td_state);

#endif
