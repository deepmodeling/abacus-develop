#include "dftu_nao_op_legacy.h"
#include "dftu_nao_pots.h"
#include "source_base/timer.h"
#include "source_base/tool_title.h"

namespace hamilt
{

template class DFTU_firstzeta<OperatorLCAO<double, double>>;

template class DFTU_firstzeta<OperatorLCAO<std::complex<double>, double>>;

template class DFTU_firstzeta<OperatorLCAO<std::complex<double>, std::complex<double>>>;

template<typename TK, typename TR>
void DFTU_firstzeta<OperatorLCAO<TK, TR>>::contributeHR()
{
    //no calculation of HR yet for DFTU operator
    return;
}

template<>
void DFTU_firstzeta<OperatorLCAO<double, double>>::contributeHk(int ik)
{
    ModuleBase::TITLE("DFTU_firstzeta", "contributeHk");
    if (!this->dftu->is_occmat_ready())
    {
        return;
    }
    ModuleBase::timer::start("DFTU_firstzeta", "contributeHk");
    // Effective potential of DFT+U is added to total Hamiltonian here; Quxin adds on 20201029
    std::vector<double> pot_uterm(this->hsk->get_pv()->nloc);

    DFTU_LCAO::cal_pot_uterm(*this->dftu, *this->ucell, this->hsk->get_pv(), isk[ik], &pot_uterm[0], this->hsk->get_sk());

    double* hk = this->hsk->get_hk();

    for (int irc = 0; irc < this->hsk->get_pv()->nloc; irc++)
    {
        hk[irc] += pot_uterm[irc];
    }

    ModuleBase::timer::end("DFTU_firstzeta", "contributeHk");
}

template<>
void DFTU_firstzeta<OperatorLCAO<std::complex<double>, double>>::contributeHk(int ik)
{
    ModuleBase::TITLE("DFTU_firstzeta", "contributeHk");
    if (!this->dftu->is_occmat_ready())
    {
        return;
    }
    ModuleBase::timer::start("DFTU_firstzeta", "contributeHk");

    // Effective potential of DFT+U is added to total Hamiltonian here; Quxin adds on 20201029
    std::vector<std::complex<double>> pot_uterm(this->hsk->get_pv()->nloc);

    DFTU_LCAO::cal_pot_uterm(*this->dftu, *this->ucell, this->hsk->get_pv(), isk[ik], &pot_uterm[0], this->hsk->get_sk());

    std::complex<double>* hk = this->hsk->get_hk();

    for (int irc = 0; irc < this->hsk->get_pv()->nloc; irc++)
    {
        hk[irc] += pot_uterm[irc];
    }

    ModuleBase::timer::end("DFTU_firstzeta", "contributeHk");
}

template<>
void DFTU_firstzeta<OperatorLCAO<std::complex<double>, std::complex<double>>>::contributeHk(int ik)
{
    ModuleBase::TITLE("DFTU_firstzeta", "contributeHk");
    if (!this->dftu->is_occmat_ready())
    {
        return;
    }
    ModuleBase::timer::start("DFTU_firstzeta", "contributeHk");
    // Effective potential of DFT+U is added to total Hamiltonian here; Quxin adds on 20201029
    std::vector<std::complex<double>> pot_uterm(this->hsk->get_pv()->nloc);

    DFTU_LCAO::cal_pot_uterm(*this->dftu, *this->ucell, this->hsk->get_pv(), isk[ik], &pot_uterm[0], this->hsk->get_sk());

    std::complex<double>* hk = this->hsk->get_hk();
    for (int irc = 0; irc < this->hsk->get_pv()->nloc; irc++)
    {
        hk[irc] += pot_uterm[irc];
    }

    ModuleBase::timer::end("DFTU_firstzeta", "contributeHk");
}

}
