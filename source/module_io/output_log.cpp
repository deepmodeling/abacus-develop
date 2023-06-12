#include "output_log.h"
#include "module_base/global_variable.h"
#include "module_base/constants.h"

namespace ModuleIO
{
void output_convergence_after_scf(bool& convergence, double& energy, std::ofstream& ofs_running)
{
    if (convergence)
    {
        ofs_running << "\n charge density convergence is achieved" << std::endl;
        ofs_running << " final etot is " << std::setprecision(11) << energy * ModuleBase::Ry_to_eV << " eV" << std::endl;
    }
    else
    {
        ofs_running << " !! convergence has not been achieved @_@" << std::endl;
        std::cout << " !! CONVERGENCE HAS NOT BEEN ACHIEVED !!" << std::endl;
    }
}

void output_efermi(bool& convergence, double& efermi, std::ofstream& ofs_running)
{
    if (convergence && GlobalV::OUT_LEVEL != "m")
    {
        ofs_running << std::setprecision(16);
        ofs_running << " EFERMI = " << std::setprecision(11) << efermi * ModuleBase::Ry_to_eV << " eV" << std::endl;
    }
}

}// namespace ModuleIO