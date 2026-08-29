#include "dftu_lcao_op_legacy.h"
#include "dftu_hamilt.h"
#include "source_base/timer.h"
#include "source_base/tool_title.h"

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
    ModuleBase::timer::start("OperatorDFTU", "contributeHk");
    // Effective potential of DFT+U is added to total Hamiltonian here; Quxin adds on 20201029
    std::vector<double> pot_uterm(this->hsk->get_pv()->nloc);

    DFTU_LCAO::pot_uterm_real(*this->dftu, ik, &pot_uterm[0], isk, this->hsk->get_sk(), this->npol);

    double* hk = this->hsk->get_hk();

    for (int irc = 0; irc < this->hsk->get_pv()->nloc; irc++)
    {
        hk[irc] += pot_uterm[irc];
    }

    ModuleBase::timer::end("OperatorDFTU", "contributeHk");
}

template<>
void OperatorDFTU<OperatorLCAO<std::complex<double>, double>>::contributeHk(int ik)
{
    ModuleBase::TITLE("OperatorDFTU", "contributeHk");
    ModuleBase::timer::start("OperatorDFTU", "contributeHk");

    // Effective potential of DFT+U is added to total Hamiltonian here; Quxin adds on 20201029
    std::vector<std::complex<double>> pot_uterm(this->hsk->get_pv()->nloc);

    DFTU_LCAO::pot_uterm_complex(*this->dftu, ik, &pot_uterm[0], isk, this->hsk->get_sk(), this->npol);

    std::complex<double>* hk = this->hsk->get_hk();

    for (int irc = 0; irc < this->hsk->get_pv()->nloc; irc++)
    {
        hk[irc] += pot_uterm[irc];
    }

    ModuleBase::timer::end("OperatorDFTU", "contributeHk");
}

template<>
void OperatorDFTU<OperatorLCAO<std::complex<double>, std::complex<double>>>::contributeHk(int ik)
{
    ModuleBase::TITLE("OperatorDFTU", "contributeHk");
    ModuleBase::timer::start("OperatorDFTU", "contributeHk");
    // Effective potential of DFT+U is added to total Hamiltonian here; Quxin adds on 20201029
    std::vector<std::complex<double>> pot_uterm(this->hsk->get_pv()->nloc);

    DFTU_LCAO::pot_uterm_complex(*this->dftu, ik, &pot_uterm[0], isk, this->hsk->get_sk(), this->npol);

    std::complex<double>* hk = this->hsk->get_hk();
    for (int irc = 0; irc < this->hsk->get_pv()->nloc; irc++)
    {
        hk[irc] += pot_uterm[irc];
    }

    ModuleBase::timer::end("OperatorDFTU", "contributeHk");
}

}
