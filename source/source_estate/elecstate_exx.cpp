#include "source_estate/elecstate.h"
#include "source_base/tool_quit.h"
#include <complex> // use std::complex

namespace elecstate
{

/// @brief calculation if converged
/// @date Peize Lin add 2016-12-03
void ElecState::set_exx(const double& Eexx, const bool cal_exx, const double hybrid_alpha)
{
    ModuleBase::TITLE("energy", "set_exx");

    if (cal_exx)
    {
        this->f_en.exx = hybrid_alpha * Eexx;
    }
    return;
}

void ElecState::set_exx(const std::complex<double>& Eexx, const bool cal_exx, const double hybrid_alpha)
{
    ModuleBase::WARNING_QUIT("ElecState::set_exx",
                             "std::complex<double> version is not implemented yet");
}

}
