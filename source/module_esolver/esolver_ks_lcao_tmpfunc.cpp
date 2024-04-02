#include "esolver_ks_lcao.h"
#include "module_hamilt_lcao/module_dftu/dftu.h"
#ifdef __DEEPKS
#include "module_hamilt_lcao/module_deepks/LCAO_deepks.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#endif
namespace ModuleESolver
{
    template <>
    void ESolver_KS_LCAO<double, double>::dftu_cal_occup_m(const int& iter, const std::vector<std::vector<double>>& dm)const
    {
        GlobalC::dftu.cal_occup_m_gamma(iter, dm, this->p_chgmix->get_mixing_beta());
    }

    template <>
    void ESolver_KS_LCAO<std::complex<double>, double>::dftu_cal_occup_m(const int& iter, const std::vector<std::vector<std::complex<double>>>& dm)const
    {
        GlobalC::dftu.cal_occup_m_k(iter, dm, this->kv, this->p_chgmix->get_mixing_beta(), this->p_hamilt);
    }
    template <>
    void ESolver_KS_LCAO<std::complex<double>, std::complex<double>>::dftu_cal_occup_m(const int& iter, const std::vector<std::vector<std::complex<double>>>& dm)const
    {
        GlobalC::dftu.cal_occup_m_k(iter, dm, this->kv, this->p_chgmix->get_mixing_beta(), this->p_hamilt);
    }
#ifdef __DEEPKS
    template<>
    void ESolver_KS_LCAO<double, double>::dpks_cal_e_delta_band(const std::vector<std::vector<double>>& dm)const
    {
        GlobalC::ld.cal_e_delta_band(dm);
    }
    template<>
    void ESolver_KS_LCAO<std::complex<double>, double>::dpks_cal_e_delta_band(const std::vector<std::vector<std::complex<double>>>& dm)const
    {
        GlobalC::ld.cal_e_delta_band_k(dm, this->kv.nks);
    }
    template<>
    void ESolver_KS_LCAO<std::complex<double>, std::complex<double>>::dpks_cal_e_delta_band(const std::vector<std::vector<std::complex<double>>>& dm)const
    {
        GlobalC::ld.cal_e_delta_band_k(dm, this->kv.nks);
    }
    template<>
    void ESolver_KS_LCAO<double, double>::dpks_cal_projected_DM(const elecstate::DensityMatrix<double, double>* dm)const
    {
        GlobalC::ld.cal_projected_DM(dm, //this->LOC.dm_gamma,
            GlobalC::ucell,
            GlobalC::ORB,
            GlobalC::GridD);
    }
    template<>
    void ESolver_KS_LCAO<std::complex<double>, double>::dpks_cal_projected_DM(const elecstate::DensityMatrix<std::complex<double>, double>* dm)const
    {
        GlobalC::ld.cal_projected_DM_k(dm, //this->LOC.dm_k,
            GlobalC::ucell,
            GlobalC::ORB,
            GlobalC::GridD);
    }
    template<>
    void ESolver_KS_LCAO<std::complex<double>, std::complex<double>>::dpks_cal_projected_DM(const elecstate::DensityMatrix<std::complex<double>, double>* dm)const
    {
        GlobalC::ld.cal_projected_DM_k(dm, //this->LOC.dm_k,
            GlobalC::ucell,
            GlobalC::ORB,
            GlobalC::GridD);
    }
#endif
}