#ifndef ELECOND_H
#define ELECOND_H

#include "module_base/matrix.h"
#include "module_basis/module_pw/pw_basis_k.h"
#include "module_cell/klist.h"
#include "module_cell/unitcell.h"
#include "module_elecstate/elecstate.h"
#include "module_hamilt_pw/hamilt_pwdft/VNL_in_pw.h"
#include "module_hamilt_pw/hamilt_pwdft/operator_pw/velocity_pw.h"

class EleCond
{
  public:
    EleCond(UnitCell* p_ucell_in, K_Vectors* p_kv_in, elecstate::ElecState* p_elec_in, ModulePW::PW_Basis_K* p_wfcpw_in,
            psi::Psi<std::complex<double>>* p_psi_in, pseudopot_cell_vnl* p_ppcell_in);
    ~EleCond(){};

    /**
     * @brief calculate Onsager coefficients Lmn(\omega) and conductivities with Kubo-Greenwood formula
     *
     * @param fwhmin FWHM for delta function
     * @param smear_type 1: Gaussian, 2: Lorentzian
     * @param wcut cutoff \omega for Lmn(\omega)
     * @param dw_in \omega step
     * @param dt_in time step
     * @param nonlocal whether to include the nonlocal potential corrections for velocity operator
     * @param wg wg(ik,ib) occupation for the ib-th band in the ik-th kpoint
     */
    void KG(const int& smear_type, const double& fwhmin, const double& wcut, const double& dw_in, const double& dt_in,
            const bool& nonlocal, ModuleBase::matrix& wg);

  protected:
    pseudopot_cell_vnl* p_ppcell = nullptr;          ///< pointer to the pseudopotential
    UnitCell* p_ucell = nullptr;                     ///< pointer to the unit cell
    ModulePW::PW_Basis_K* p_wfcpw = nullptr;         ///< pointer to the plane wave basis
    K_Vectors* p_kv = nullptr;                       ///< pointer to the k vectors
    elecstate::ElecState* p_elec = nullptr;          ///< pointer to the electronic state
    psi::Psi<std::complex<double>>* p_psi = nullptr; ///< pointer to the wavefunction

  protected:
    /**
     * @brief calculate the response function Cmn(t) for currents
     *
     * @param ik k point
     * @param nt number of steps of time
     * @param dt time step
     * @param decut ignore dE which is larger than decut
     * @param wg wg(ik,ib) occupation for the ib-th band in the ik-th kpoint
     * @param velop velocity operator
     * @param ct11 C11(t)
     * @param ct12 C12(t)
     * @param ct22 C22(t)
     */
    void jjresponse_ks(const int ik, const int nt, const double dt, const double decut, ModuleBase::matrix& wg,
                       hamilt::Velocity& velop, double* ct11, double* ct12, double* ct22);
    /**
     * @brief Calculate the conductivity using the response function
     *
     * @param nt number of time steps
     * @param dt time step
     * @param smear_type smearing type 1: gaussian, 2: lorentzian
     * @param fwhmin full width at half maximum of the smearing function
     * @param wcut cutoff frequency
     * @param dw_in frequency step
     * @param ct11 C11 component of the response function
     * @param ct12 C12 component of the response function
     * @param ct22 C22 component of the response function
     */
    void calcondw(const int nt, const double dt, const int& smear_type, const double fwhmin, const double wcut,
                  const double dw_in, double* ct11, double* ct12, double* ct22);
};

#endif // ELECOND_H