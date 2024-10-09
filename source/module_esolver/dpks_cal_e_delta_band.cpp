#include "esolver_ks_lcao.h"
#ifdef __DEEPKS
#include "module_hamilt_lcao/module_deepks/LCAO_deepks.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#endif

namespace ModuleESolver
{

    using namespace std;

#ifdef __DEEPKS
    template<>
    void ESolver_KS_LCAO<double, double>::dpks_cal_e_delta_band(
       const vector<vector<double>>& dm)const
    {
        GlobalC::ld.cal_e_delta_band(dm);
    }


    template<>
    void ESolver_KS_LCAO<complex<double>, double>::dpks_cal_e_delta_band(
        const vector<vector<complex<double>>>& dm)const
    {
        GlobalC::ld.cal_e_delta_band_k(dm, this->kv.get_nks());
    }


    template<>
    void ESolver_KS_LCAO<complex<double>, complex<double>>::dpks_cal_e_delta_band(
        const vector<vector<complex<double>>>& dm)const
    {
        GlobalC::ld.cal_e_delta_band_k(dm, this->kv.get_nks());
    }


    template<>
    void ESolver_KS_LCAO<double, double>::dpks_cal_projected_DM(
        const elecstate::DensityMatrix<double, double>* dm)const
    {
        GlobalC::ld.cal_projected_DM(dm,
            GlobalC::ucell,
            orb_,
            GlobalC::GridD);
    }


    template<>
    void ESolver_KS_LCAO<complex<double>, double>::dpks_cal_projected_DM(
          const elecstate::DensityMatrix<complex<double>, double>* dm)const
    {
        GlobalC::ld.cal_projected_DM_k(dm,
            GlobalC::ucell,
            orb_,
            GlobalC::GridD);
    }


    template<>
    void ESolver_KS_LCAO<complex<double>, complex<double>>::dpks_cal_projected_DM(
          const elecstate::DensityMatrix<complex<double>, double>* dm)const
    {
        GlobalC::ld.cal_projected_DM_k(dm,
            GlobalC::ucell,
            orb_,
            GlobalC::GridD);
    }
#endif
}
