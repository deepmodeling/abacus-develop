#include "module_hamilt_pw/hamilt_pwdft/global.h"

namespace elecstate
{

#ifdef __EXX
#ifdef __LCAO
/// @brief calculation if converged
/// @date Peize Lin add 2016-12-03
void ElecState::set_exx(const double& Eexx)
{
    ModuleBase::TITLE("energy", "set_exx");

    if (GlobalC::exx_info.info_global.cal_exx)
    {
        this->f_en.exx = GlobalC::exx_info.info_global.hybrid_alpha * Eexx;
    }
    return;
}
void ElecState::set_exx(const std::complex<double>& Eexx)
{
    ModuleBase::TITLE("energy", "set_exx");

    if (GlobalC::exx_info.info_global.cal_exx)
    {
        this->f_en.exx = GlobalC::exx_info.info_global.hybrid_alpha * std::real(Eexx);
    }
    return;
}
#endif //__LCAO
#endif //__EXX

}