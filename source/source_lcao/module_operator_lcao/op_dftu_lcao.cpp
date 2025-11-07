#include "op_dftu_lcao.h"
#include "source_base/timer.h"
#include "source_base/tool_title.h"
#include "source_lcao/module_dftu/dftu.h"
#include "source_pw/module_pwdft/global.h"

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
    std::vector<double> eff_pot(this->hsk->get_pv()->nloc);

    this->dftu->cal_eff_pot_mat_real(ik, &eff_pot[0], isk, this->hsk->get_sk());

    double* hk = this->hsk->get_hk();

    for (int irc = 0; irc < this->hsk->get_pv()->nloc; irc++)
    {
        hk[irc] += eff_pot[irc];
    }

    ModuleBase::timer::tick("OperatorDFTU", "contributeHk");
}

template<>
void OperatorDFTU<OperatorLCAO<std::complex<double>, double>>::contributeHk(int ik)
{
    ModuleBase::TITLE("OperatorDFTU", "contributeHk");
    ModuleBase::timer::tick("OperatorDFTU", "contributeHk");

    // Effective potential of DFT+U is added to total Hamiltonian here; Quxin adds on 20201029
    std::vector<std::complex<double>> eff_pot(this->hsk->get_pv()->nloc);

    this->dftu->cal_eff_pot_mat_complex(ik, &eff_pot[0], isk, this->hsk->get_sk());

    std::complex<double>* hk = this->hsk->get_hk();

    for (int irc = 0; irc < this->hsk->get_pv()->nloc; irc++)
    {
        hk[irc] += eff_pot[irc];
    }

    ModuleBase::timer::tick("OperatorDFTU", "contributeHk");
}

template<>
void OperatorDFTU<OperatorLCAO<std::complex<double>, std::complex<double>>>::contributeHk(int ik)
{
    ModuleBase::TITLE("OperatorDFTU", "contributeHk");
    ModuleBase::timer::tick("OperatorDFTU", "contributeHk");
    // Effective potential of DFT+U is added to total Hamiltonian here; Quxin adds on 20201029
    std::vector<std::complex<double>> eff_pot(this->hsk->get_pv()->nloc);

    this->dftu->cal_eff_pot_mat_complex(ik, &eff_pot[0], isk, this->hsk->get_sk());

    std::complex<double>* hk = this->hsk->get_hk();
    for (int irc = 0; irc < this->hsk->get_pv()->nloc; irc++)
    {
        hk[irc] += eff_pot[irc];
    }

    ModuleBase::timer::tick("OperatorDFTU", "contributeHk");
}

}
