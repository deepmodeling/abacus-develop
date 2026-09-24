#ifndef LCAO_ALLOCATE_DM_H
#define LCAO_ALLOCATE_DM_H

#include "source_estate/module_dm/dm_holder.h"

class K_Vectors;
class Parallel_Orbitals;

namespace LCAO_domain
{
// Allocate the DensityMatrix held by dmat, based on the k-point list and
// orbital parallel layout. Free function so that module_dm itself does not
// depend on LCAO-layer modules (K_Vectors, Parallel_Orbitals).
template <typename TK>
void allocate_dm(module_dm::Setup_DM<TK>& dmat,
                 const K_Vectors* kv,
                 const Parallel_Orbitals* pv,
                 const int nspin);
} // namespace LCAO_domain

#endif
