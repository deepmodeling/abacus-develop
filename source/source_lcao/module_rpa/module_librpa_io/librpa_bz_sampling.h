#ifndef LIBRPA_BZ_SAMPLING_H
#define LIBRPA_BZ_SAMPLING_H

#include "source_base/vector3.h"

namespace RpaLriDetail
{
ModuleBase::Vector3<double> librpa_cartesian_kvec(const ModuleBase::Vector3<double>& kvec_c,
                                                  double lat0_bohr);
} // namespace RpaLriDetail

#endif
