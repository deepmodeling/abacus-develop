#ifndef SOURCE_MD_MODULE_DP_DP_APARAM_H
#define SOURCE_MD_MODULE_DP_DP_APARAM_H

#include <cstdint>
#include <vector>

struct LocalAtom;

namespace ModuleMD
{
namespace DPDetail
{
std::vector<double> arrange_aparam(
    const std::vector<double> &input,
    const std::vector<LocalAtom> &owned,
    const std::vector<LocalAtom> &ghost,
    const std::vector<std::int64_t> &type_counts,
    int dimension,
    bool include_ghost);
}
} // namespace ModuleMD

#endif
