#ifdef __EXX
#include "op_exx_lcao.h"
#include "module_base/timer.h"
#include "module_base/tool_title.h"
#include "module_hamilt_general/module_xc/xc_functional.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_ri/RI_2D_Comm.h"
#include "module_base/blacs_connector.h"

namespace hamilt
{

template class OperatorEXX<OperatorLCAO<double, double>>;

template class OperatorEXX<OperatorLCAO<std::complex<double>, double>>;

template class OperatorEXX<OperatorLCAO<std::complex<double>, std::complex<double>>>;

template<typename TK, typename TR>
void OperatorEXX<OperatorLCAO<TK, TR>>::contributeHR()
{

}

template<>
void OperatorEXX<OperatorLCAO<double, double>>::add_loaded_Hexx(const int ik)
{
    BlasConnector::axpy(this->LM->ParaV->get_local_size(), 1.0, this->LM->Hexxd_k_load[ik].data(), 1, this->LM->Hloc.data(), 1);
}
template<>
void OperatorEXX<OperatorLCAO<std::complex<double>, double>>::add_loaded_Hexx(const int ik)
{

    BlasConnector::axpy(this->LM->ParaV->get_local_size(), 1.0, this->LM->Hexxc_k_load[ik].data(), 1, this->LM->Hloc2.data(), 1);
}
template<>
void OperatorEXX<OperatorLCAO<std::complex<double>, std::complex<double>>>::add_loaded_Hexx(const int ik)
{
    BlasConnector::axpy(this->LM->ParaV->get_local_size(), 1.0, this->LM->Hexxc_k_load[ik].data(), 1, this->LM->Hloc2.data(), 1);
}


template<typename TK, typename TR>
void OperatorEXX<OperatorLCAO<TK, TR>>::contributeHk(int ik)
{
    // Peize Lin add 2016-12-03
    if (this->two_level_step != nullptr && *this->two_level_step == 0 && !this->restart) return;  //in the non-exx loop, do nothing 
    if (XC_Functional::get_func_type() == 4 || XC_Functional::get_func_type() == 5)
    {
        if (this->restart && this->two_level_step != nullptr)
        {
            if (*this->two_level_step == 0)
            {
                add_loaded_Hexx(ik);
                return;
            }
            else // clear loaded Hexx and release memory
            {
                if (this->LM->Hexxd_k_load.size() > 0)
                {
                    this->LM->Hexxd_k_load.clear();
                    this->LM->Hexxd_k_load.shrink_to_fit();
                }
                else if (this->LM->Hexxc_k_load.size() > 0)
                {
                    this->LM->Hexxc_k_load.clear();
                    this->LM->Hexxc_k_load.shrink_to_fit();
                }
            }
        }
        // cal H(k) from H(R) normally

        if (GlobalC::exx_info.info_ri.real_number)
            RI_2D_Comm::add_Hexx(
                kv,
                ik,
                GlobalC::exx_info.info_global.hybrid_alpha,
                *this->Hexxd == nullptr ? *this->LM->Hexxd : *this->Hexxd,
                *this->LM->ParaV,
                *this->hK);
        else
            RI_2D_Comm::add_Hexx(
                kv,
                ik,
                GlobalC::exx_info.info_global.hybrid_alpha,
                this->Hexxc == nullptr ? *this->LM->Hexxc : *this->Hexxc,
                *this->LM->ParaV,
                *this->hK);
    }
}

} // namespace hamilt
#endif