#include "density_matrix.h"

#include "source_base/libm/libm.h"
#include "source_base/tool_title.h"
#include "source_base/tool_quit.h"
#include "source_base/constants.h"
#include "source_base/timer.h"
#include "source_cell/klist.h"

namespace module_dm
{

template <typename TK, typename TR_in, typename TR_out>
void cal_dmr_td(
    DensityMatrix<TK, TR_in>& dm,
    std::vector<hamilt::HContainer<TR_out>*>& dmR_out,
    const std::map<ModuleBase::Vector3<int>, std::complex<double>>& phase_hybrid,
    const ModuleBase::Vector3<double> At,
    const int ik_in)
{
    ModuleBase::TITLE("DensityMatrix", "cal_dmr_td");
    assert(dmR_out.size() == dm.spin_mult && "DMR has not been initialized!");

    // validate ik_in: either -1 (all k-points) or a valid index
    if (ik_in < -1 || ik_in >= dm._nk)
    {
        ModuleBase::WARNING_QUIT("module_dm::cal_dmr_td",
                                 "ik_in out of range: must be -1 (all k) or 0 <= ik_in < nk");
    }

    ModuleBase::timer::start("DensityMatrix", "cal_dmr_td");
    accumulate_dmr(dm, dmR_out, phase_hybrid, ik_in, "module_dm::cal_dmr_td");
    dm._dmr_ready = true;
    ModuleBase::timer::end("DensityMatrix", "cal_dmr_td");
}

template <>
void DensityMatrix<double, double>::cal_dmr_td(
    const std::map<ModuleBase::Vector3<int>, std::complex<double>>& phase_hybrid,
    const ModuleBase::Vector3<double> At,
    const int ik_in)
{
    return;
}
template <>
void DensityMatrix<std::complex<double>, double>::cal_dmr_td(
    const std::map<ModuleBase::Vector3<int>, std::complex<double>>& phase_hybrid,
    const ModuleBase::Vector3<double> At,
    const int ik_in)
{
    module_dm::cal_dmr_td(*this, this->dmr, phase_hybrid, At, ik_in);
}

template <>
void DensityMatrix<std::complex<double>, std::complex<double>>::cal_dmr_td(
    const std::map<ModuleBase::Vector3<int>, std::complex<double>>& phase_hybrid,
    const ModuleBase::Vector3<double> At,
    const int ik_in)
{
    module_dm::cal_dmr_td(*this, this->dmr, phase_hybrid, At, ik_in);
}

} // namespace module_dm
