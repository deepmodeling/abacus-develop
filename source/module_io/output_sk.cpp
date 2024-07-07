#include "output_sk.h"

#include "module_base/tool_quit.h"

namespace ModuleIO
{

template <typename TK>
Output_Sk<TK>::Output_Sk(LCAO_Matrix* LM, hamilt::Hamilt<TK>* p_hamilt, Parallel_Orbitals* ParaV, int nspin, int nks)
    : LM_(LM), p_hamilt_(p_hamilt), ParaV_(ParaV), nspin_(nspin), nks_(nks)
{
}

template <>
double* Output_Sk<double>::get_Sk(int ik)
{
    if (ik < 0 || ik >= this->nks_)
    {
        ModuleBase::WARNING_QUIT("Output_Sk::get_sk", "ik out of range");
    }
    return dynamic_cast<hamilt::HamiltLCAO<double, double>*>(this->p_hamilt_)->getSk();
}

template <>
std::complex<double>* Output_Sk<std::complex<double>>::get_Sk(int ik)
{
    if (ik < 0 || ik >= this->nks_)
    {
        ModuleBase::WARNING_QUIT("Output_Sk::get_sk", "ik out of range");
    }
    if (this->nspin_ == 4)
    {
        dynamic_cast<hamilt::HamiltLCAO<std::complex<double>, std::complex<double>>*>(this->p_hamilt_)
            ->updateSk(ik, 1);
        return dynamic_cast<hamilt::HamiltLCAO<std::complex<double>, std::complex<double>>*>(this->p_hamilt_)->getSk();
    }
    else
    {
        dynamic_cast<hamilt::HamiltLCAO<std::complex<double>, double>*>(this->p_hamilt_)->updateSk(ik, 1);
        return dynamic_cast<hamilt::HamiltLCAO<std::complex<double>, double>*>(this->p_hamilt_)->getSk();
    }
    return nullptr;
}

template class Output_Sk<double>;
template class Output_Sk<std::complex<double>>;

} // namespace ModuleIO