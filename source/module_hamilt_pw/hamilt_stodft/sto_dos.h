#ifndef STO_DOS
#define STO_DOS
#include "module_hamilt_pw/hamilt_stodft/sto_wf.h"
#include "module_hsolver/hsolver_pw_sdft.h"

class Sto_DOS
{
  public:
    Sto_DOS(ModulePW::PW_Basis_K* p_wfcpw_in, K_Vectors* p_kv_in, elecstate::ElecState* p_elec_in,
            psi::Psi<std::complex<double>>* p_psi_in, hamilt::Hamilt<std::complex<double>>* p_hamilt_in,
            hsolver::HSolverPW_SDFT* p_hsol_in, Stochastic_WF* p_stowf_in);
    /**
     * @brief decide the parameters for the DOS calculation
     *
     * @param dos_nche Number of Chebyshev orders
     * @param emin_sto Emin input for sdft
     * @param emax_sto Emax input for sdft
     * @param dos_setemin whether to set the minimum energy
     * @param dos_setemax whether to set the maximum energy
     * @param dos_emin_ev Emin input for DOS
     * @param dos_emax_ev Emax input for DOS
     * @param dos_scale dos_scale input for DOS
     */
    void decide_param(const int& dos_nche, const double& emin_sto, const double& emax_sto, const bool& dos_setemin,
                      const bool& dos_setemax, const double& dos_emin_ev, const double& dos_emax_ev,
                      const double& dos_scale);
    /**
     * @brief Calculate DOS using stochastic wavefunctions
     *
     * @param sigmain  sigma for the gaussian broadening
     * @param de       energy step
     * @param npart    number of parts to reduce the memory usage
     */
    void caldos(const double sigmain, const double de, const int npart);

  protected:
    int nbands_ks = 0;                               ///< number of KS bands
    int nbands_sto = 0;                              ///< number of stochastic bands
    int dos_nche = 0;                                ///< number of Chebyshev orders for DOS
    double emax = 0.0;                               ///< maximum energy
    double emin = 0.0;                               ///< minimum energy
    ModulePW::PW_Basis_K* p_wfcpw = nullptr;         ///< pointer to the plane wave basis
    K_Vectors* p_kv = nullptr;                       ///< pointer to the k vectors
    elecstate::ElecState* p_elec = nullptr;          ///< pointer to the electronic state
    psi::Psi<std::complex<double>>* p_psi = nullptr; ///< pointer to the wavefunction
    hamilt::Hamilt<std::complex<double>>* p_hamilt;  ///< pointer to the Hamiltonian
    hsolver::HSolverPW_SDFT* p_hsol = nullptr;       ///< pointer to the Hamiltonian solver
    Stochastic_WF* p_stowf = nullptr;                ///< pointer to the stochastic wavefunctions
};

#endif // STO_DOS