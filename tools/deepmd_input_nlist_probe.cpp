#if defined(__DPMDC)
#include <deepmd/deepmd.hpp>
using DPProbeInputNlist = deepmd::hpp::InputNlist;
#else
#include <deepmd/DeepPot.h>
using DPProbeInputNlist = deepmd::InputNlist;
#endif

int main()
{
#if defined(DPMD_PROBE_COMM_13)
    DPProbeInputNlist nlist(0,
                            nullptr,
                            nullptr,
                            nullptr,
                            0,
                            nullptr,
                            nullptr,
                            nullptr,
                            nullptr,
                            nullptr,
                            nullptr,
                            nullptr,
                            1);
#elif defined(DPMD_PROBE_COMM_12)
    DPProbeInputNlist nlist(0,
                            nullptr,
                            nullptr,
                            nullptr,
                            0,
                            nullptr,
                            nullptr,
                            nullptr,
                            nullptr,
                            nullptr,
                            nullptr,
                            nullptr);
#elif defined(DPMD_PROBE_MAPPING)
    DPProbeInputNlist nlist(0, nullptr, nullptr, nullptr);
    nlist.set_mapping(nullptr);
#else
#error "Select one DeePMD InputNlist probe"
#endif
    return 0;
}
