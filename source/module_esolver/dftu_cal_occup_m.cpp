#include "esolver_ks_lcao.h"
#include "module_hamilt_lcao/module_dftu/dftu.h"

namespace ModuleESolver
{

    using namespace std;

    //! dftu occupation matrix for gamma only using dm(double)
    template <>
    void ESolver_KS_LCAO<double, double>::dftu_cal_occup_m(
         const int& iter, 
         const vector<vector<double>>& dm)const
    {
		GlobalC::dftu.cal_occup_m_gamma(
				iter, 
				dm, 
				this->p_chgmix->get_mixing_beta(),
                this->p_hamilt);
    }

    //! dftu occupation matrix for multiple k-points using dm(complex)
    template <>
    void ESolver_KS_LCAO<complex<double>, double>::dftu_cal_occup_m(
       const int& iter, 
       const vector<vector<complex<double>>>& dm)const
    {
        GlobalC::dftu.cal_occup_m_k(
                 iter, 
                 dm, 
                 this->kv, 
                 this->p_chgmix->get_mixing_beta(), 
                 this->p_hamilt);
    }

    //! dftu occupation matrix
    template <>
    void ESolver_KS_LCAO<complex<double>, complex<double>>::dftu_cal_occup_m(
       const int& iter, 
       const vector<vector<complex<double>>>& dm)const
    {
        GlobalC::dftu.cal_occup_m_k(
                 iter, 
                 dm, 
                 this->kv, 
                 this->p_chgmix->get_mixing_beta(), 
                 this->p_hamilt);
    }

}
