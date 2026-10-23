#ifndef LIBRPA_BZ_SAMPLING_H
#define LIBRPA_BZ_SAMPLING_H

#include "source_base/vector3.h"

namespace RpaLriDetail
{
struct LibRpaStoredQIndex
{
    int coulomb_irreducible_index;
    int representative_scf_index;
};

ModuleBase::Vector3<double> librpa_cartesian_kvec(const ModuleBase::Vector3<double>& kvec_c,
                                                  double lat0_bohr);

int librpa_stored_coulomb_q_count(int stored_kpoint_count);

LibRpaStoredQIndex librpa_stored_q_index(int stored_kpoint_index, int stored_kpoint_count);
} // namespace RpaLriDetail

#endif
