#include "esolver_ks_lcao.h"
#include "module_hamilt_lcao/module_dftu/dftu.h"
namespace ModuleESolver
{
    template <>
    void ESolver_KS_LCAO<double, double>::dftu_cal_occup_m(int& iter, const std::vector<std::vector<double>>& dm)const
    {
        GlobalC::dftu.cal_occup_m_gamma(iter, dm, this->p_chgmix->get_mixing_beta());
    }

    template <>
    void ESolver_KS_LCAO<std::complex<double>, double>::dftu_cal_occup_m(int& iter, const std::vector<std::vector<std::complex<double>>>& dm)const
    {
        GlobalC::dftu.cal_occup_m_k(iter, dm, this->kv, this->p_chgmix->get_mixing_beta(), this->p_hamilt);
    }
    template <>
    void ESolver_KS_LCAO<std::complex<double>, std::complex<double>>::dftu_cal_occup_m(int& iter, const std::vector<std::vector<std::complex<double>>>& dm)const
    {
        GlobalC::dftu.cal_occup_m_k(iter, dm, this->kv, this->p_chgmix->get_mixing_beta(), this->p_hamilt);
    }
}