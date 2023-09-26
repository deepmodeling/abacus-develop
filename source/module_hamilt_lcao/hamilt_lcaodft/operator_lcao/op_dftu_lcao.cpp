#include "op_dftu_lcao.h"
#include "module_base/timer.h"
#include "module_base/tool_title.h"
#include "module_hamilt_lcao/module_dftu/dftu.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

namespace hamilt
{

template class OperatorDFTU<OperatorLCAO<double, double>>;

template class OperatorDFTU<OperatorLCAO<std::complex<double>, double>>;

template class OperatorDFTU<OperatorLCAO<std::complex<double>, std::complex<double>>>;

template<typename TK, typename TR>
void OperatorDFTU<OperatorLCAO<TK, TR>>::contributeHR()
{
    //no calculation of HR yet for DFTU operator
    return;
}

template<>
void OperatorDFTU<OperatorLCAO<double, double>>::contributeHk(int ik)
{
    ModuleBase::TITLE("OperatorDFTU", "contributeHk");
    ModuleBase::timer::tick("OperatorDFTU", "contributeHk");
    // Effective potential of DFT+U is added to total Hamiltonian here; Quxin adds on 20201029
    std::vector<double> eff_pot(this->LM->ParaV->nloc);
    GlobalC::dftu.cal_eff_pot_mat_real(ik, &eff_pot[0], isk);

    for (int irc = 0; irc < this->LM->ParaV->nloc; irc++)
    {
        this->LM->Hloc[irc] += eff_pot[irc];
    }

    ModuleBase::timer::tick("OperatorDFTU", "contributeHk");
}

template<>
void OperatorDFTU<OperatorLCAO<std::complex<double>, double>>::contributeHk(int ik)
{
    ModuleBase::TITLE("OperatorDFTU", "contributeHk");
    ModuleBase::timer::tick("OperatorDFTU", "contributeHk");
    // Effective potential of DFT+U is added to total Hamiltonian here; Quxin adds on 20201029
    std::vector<std::complex<double>> eff_pot(this->LM->ParaV->nloc);
    GlobalC::dftu.cal_eff_pot_mat_complex(ik, &eff_pot[0], isk);

    for (int irc = 0; irc < this->LM->ParaV->nloc; irc++)
    {
        this->LM->Hloc2[irc] += eff_pot[irc];
    }

    ModuleBase::timer::tick("OperatorDFTU", "contributeHk");
}

template<>
void OperatorDFTU<OperatorLCAO<std::complex<double>, std::complex<double>>>::contributeHk(int ik)
{
    ModuleBase::TITLE("OperatorDFTU", "contributeHk");
    ModuleBase::timer::tick("OperatorDFTU", "contributeHk");
    // Effective potential of DFT+U is added to total Hamiltonian here; Quxin adds on 20201029
    std::vector<std::complex<double>> eff_pot(this->LM->ParaV->nloc);
    GlobalC::dftu.cal_eff_pot_mat_complex(ik, &eff_pot[0], isk);

    for (int irc = 0; irc < this->LM->ParaV->nloc; irc++)
    {
        this->LM->Hloc2[irc] += eff_pot[irc];
    }

    ModuleBase::timer::tick("OperatorDFTU", "contributeHk");
}

}