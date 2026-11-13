#include "dftu_nao_op_legacy.h"
#include "dftu_nao_pots.h"
#include "source_base/parallel_blkcyc.h"
#include "source_base/timer.h"
#include "source_base/tool_quit.h"
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

template<typename TK, typename TR>
void DFTU_firstzeta<OperatorLCAO<TK, TR>>::contributeHk(int ik)
{
    ModuleBase::TITLE("DFTU_firstzeta", "contributeHk");
    if (!this->dftu->is_occmat_ready())
    {
        return;
    }
    ModuleBase::timer::start("DFTU_firstzeta", "contributeHk");
    // Effective potential of DFT+U is added to total Hamiltonian here; Quxin adds on 20201029
    const Parallel_Blkcyclic* dist = this->hR->matrix_distribution();
    if (dist == nullptr)
    {
        ModuleBase::WARNING_QUIT("DFTU_firstzeta::contributeHk", "matrix distribution is not initialized");
    }
    std::vector<TK> pot_uterm(dist->size());

    DFTU_LCAO::cal_pot_uterm(
        *this->dftu, *this->ucell, *dist, isk[ik], pot_uterm.data(), this->hsk->get_sk());

    TK* hk = this->hsk->get_hk();

    for (std::size_t irc = 0; irc < dist->size(); ++irc)
    {
        hk[irc] += pot_uterm[irc];
    }

    ModuleBase::timer::end("DFTU_firstzeta", "contributeHk");
}

}
