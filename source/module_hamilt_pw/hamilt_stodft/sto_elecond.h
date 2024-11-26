#ifndef STOELECOND_H
#define STOELECOND_H

#include "module_hamilt_general/hamilt.h"
#include "module_hamilt_pw/hamilt_pwdft/elecond.h"
#include "module_hamilt_pw/hamilt_stodft/sto_wf.h"
#include "module_hsolver/hsolver_pw_sdft.h"

class Sto_EleCond : protected EleCond
{
  public:
    Sto_EleCond(UnitCell* p_ucell_in,
                K_Vectors* p_kv_in,
                elecstate::ElecState* p_elec_in,
                ModulePW::PW_Basis_K* p_wfcpw_in,
                psi::Psi<std::complex<double>>* p_psi_in,
                pseudopot_cell_vnl* p_ppcell_in,
                hamilt::Hamilt<std::complex<double>>* p_hamilt_in,
                StoChe<double>& stoche,
                Stochastic_WF<std::complex<double>, base_device::DEVICE_CPU>* p_stowf_in);
    ~Sto_EleCond(){};
    /**
     * @brief Set the N order of Chebyshev expansion for conductivities
     *        It will change class member : fd_nche, cond_nche
     *
     * @param dt t step
     * @param nbatch number of t batch
     * @param cond_thr threshold of errors for conductivities
     * @param fd_nche  N order of Chebyshev for Fermi-Dirac function
     * @param try_emin trial Emin
     * @param try_emax trial Emax
     *
     */
    void decide_nche(const double dt, const double cond_thr, const int& fd_nche, double try_emin, double try_emax);
    /**
     * @brief calculate Stochastic Kubo-Greenwood
     *
     * @param fwhmin FWHM
     * @param smear_type 1: Gaussian, 2: Lorentzian
     * @param wcut cutoff omega
     * @param dw_in omega step
     * @param dt_in t step
     * @param nonlocal whether to include the nonlocal potential corrections for velocity operator
     * @param nbatch t step batch
     * @param npart_sto number stochastic wavefunctions parts to evalution simultaneously
     */
    void sKG(const int& smear_type,
             const double& fwhmin,
             const double& wcut,
             const double& dw_in,
             const double& dt_in,
             const bool& nonlocal,
             const int& npart_sto);

  protected:
    int nbands_ks = 0;                              ///< number of KS bands
    int nbands_sto = 0;                             ///< number of stochastic bands
    int cond_nche = 0;                              ///< number of Chebyshev orders for conductivities
    int fd_nche = 0;                                ///< number of Chebyshev orders for Fermi-Dirac function
    int cond_dtbatch = 0;                           ///< number of time steps in a batch
    hamilt::Hamilt<std::complex<double>>* p_hamilt; ///< pointer to the Hamiltonian
    Stochastic_WF<std::complex<double>, base_device::DEVICE_CPU>* p_stowf
        = nullptr;            ///< pointer to the stochastic wavefunctions
    Sto_Func<double> stofunc; ///< functions

    hamilt::HamiltSdftPW<std::complex<double>>* p_hamilt_sto = nullptr; ///< pointer to the Hamiltonian for sDFT

  protected:
    /**
     * @brief calculate Jmatrix  <leftv|J|rightv>
     *
     */
    void cal_jmatrix(const psi::Psi<std::complex<float>>& kspsi_all,
                     const psi::Psi<std::complex<float>>& vkspsi,
                     const double* en,
                     const double* en_all,
                     std::complex<double>* leftfact,
                     std::complex<double>* rightfact,
                     const psi::Psi<std::complex<double>>& leftchi,
                     psi::Psi<std::complex<double>>& rightchi,
                     psi::Psi<std::complex<double>>& left_hchi,
                     psi::Psi<std::complex<double>>& batch_vchi,
                     psi::Psi<std::complex<double>>& batch_vhchi,
#ifdef __MPI
                     psi::Psi<std::complex<float>>& chi_all,
                     psi::Psi<std::complex<float>>& hchi_all,
                     void* gatherinfo_ks,
                     void* gatherinfo_sto,
#endif
                     const int& bsize_psi,
                     std::vector<std::complex<float>>& j1,
                     std::vector<std::complex<float>>& j2,
                     hamilt::Velocity& velop,
                     const int& ik,
                     const std::complex<double>& factor,
                     const int bandinfo[6]);
};
#endif // ELECOND_H