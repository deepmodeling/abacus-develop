#ifndef HAMILTSDFTPW_H
#define HAMILTSDFTPW_H

#include "module_hamilt_pw/hamilt_pwdft/hamilt_pw.h"

namespace hamilt
{

template <typename T, typename Device = base_device::DEVICE_CPU>
class HamiltSdftPW : public HamiltPW<T, Device>
{
  public:
    /**
     * @brief Construct a new Hamilt Sdft P W object
     *
     * @param pot_in potential
     * @param wfc_basis pw basis for wave functions
     * @param p_kv k vectors
     * @param npol the length of wave function is npol * npwk_max
     * @param emin_in Emin of the Hamiltonian
     * @param emax_in Emax of the Hamiltonian
     */
    HamiltSdftPW(elecstate::Potential* pot_in,
                 ModulePW::PW_Basis_K* wfc_basis,
                 K_Vectors* p_kv,
                 pseudopot_cell_vnl* nlpp,
                 const UnitCell* ucell,
                 const int& npol,
                 double* emin_in,
                 double* emax_in);
    /**
     * @brief Destroy the Hamilt Sdft P W object
     *
     */
    ~HamiltSdftPW(){};

    // void update_emin_emax(const double& emin, const double& emax)
    // {
    //     this->emin = &emin;
    //     this->emax = &emax;
    // }

    /**
     * @brief Calculate \hat{H}|psi>
     *
     * @param psi_in input wave function
     * @param hpsi output wave function
     * @param nbands number of bands
     */
    void hPsi(const T* psi_in, T* hpsi, const int& nbands = 1);

    /**
     * @brief Calculate \hat{H}|psi> and normalize it
     *
     * @param psi_in input wave function
     * @param hpsi output wave function
     * @param nbands number of bands
     */
    void hPsi_norm(const T* psi_in, T* hpsi, const int& nbands = 1);

    double* emin = nullptr; ///< Emin of the Hamiltonian
    double* emax = nullptr; ///< Emax of the Hamiltonian

  private:
    int npwk_max = 0;      ///< maximum number of plane waves
    int npol = 0;          ///< number of polarizations
    std::vector<int>& ngk; ///< number of G vectors
};

} // namespace hamilt

#endif