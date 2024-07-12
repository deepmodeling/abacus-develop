#ifdef __EXX
#include "op_exx_lcao.h"
#include "module_base/blacs_connector.h"

namespace hamilt
{

template<>
void OperatorEXX<OperatorLCAO<double, double>>::add_loaded_Hexx(const int ik)
{
    auto* paraV = this->hsk->get_pv();// get parallel orbitals from HK
    BlasConnector::axpy(paraV->get_local_size(), 1.0, this->LM->Hexxd_k_load[ik].data(), 1, this->hsk->get_hk(), 1);
}
template<>
void OperatorEXX<OperatorLCAO<std::complex<double>, double>>::add_loaded_Hexx(const int ik)
{
    auto* paraV = this->hsk->get_pv();// get parallel orbitals from HK
    BlasConnector::axpy(paraV->get_local_size(), 1.0, this->LM->Hexxc_k_load[ik].data(), 1, this->hsk->get_hk(), 1);
}
template<>
void OperatorEXX<OperatorLCAO<std::complex<double>, std::complex<double>>>::add_loaded_Hexx(const int ik)
{
    auto* paraV = this->hsk->get_pv();// get parallel orbitals from HK
    BlasConnector::axpy(paraV->get_local_size(), 1.0, this->LM->Hexxc_k_load[ik].data(), 1, this->hsk->get_hk(), 1);
}

} // namespace hamilt
#endif