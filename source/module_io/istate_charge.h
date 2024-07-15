#ifndef ISTATE_CHARGE_H
#define ISTATE_CHARGE_H
#include "module_basis/module_pw/pw_basis.h"
#include "module_elecstate/module_dm/density_matrix.h"
#include "module_hamilt_lcao/module_gint/gint.h"
#include "module_hamilt_lcao/module_gint/gint_gamma.h"
#include "module_psi/psi.h"

#include <module_base/complexmatrix.h>
#include <module_base/matrix.h>
#include <stdexcept>
#include <vector>

/**
 * @brief Manages the computation of the charge density for different bands.
 *
 * This class is responsible for initializing and managing the
 * charge state computation process, offering functionality to
 * calculate and plot the decomposed charge density below and above
 * the Fermi surface based on specified bands.
 */
class IState_Charge
{
  public:
    IState_Charge(psi::Psi<double>* psi_gamma_in, const Parallel_Orbitals* ParaV_in);
    IState_Charge(psi::Psi<std::complex<double>>* psi_k_in, const Parallel_Orbitals* ParaV_in)
    {
        throw std::logic_error("IState_Charge for multi-k is not implemented.");
    };

    ~IState_Charge();

    void begin(Gint_Gamma& gg,
               double** rho,
               const ModuleBase::matrix& wg,
               const std::vector<double>& ef_all_spin,
               const int rhopw_nrxx,
               const int rhopw_nplane,
               const int rhopw_startz_current,
               const int rhopw_nx,
               const int rhopw_ny,
               const int rhopw_nz,
               const int bigpw_bz,
               const int bigpw_nbz,
               const bool gamma_only_local,
               const int nbands_istate,
               const std::vector<int>& out_band_kb,
               const int nbands,
               const double nelec,
               const int nspin,
               const int nlocal,
               const std::string& global_out_dir,
               const int my_rank,
               std::ofstream& ofs_warning,
               const UnitCell* ucell_in,
               Grid_Driver* GridD_in);

  private:
    std::vector<int> bands_picked_;

#ifdef __MPI
    /**
     * @brief Calculates the density matrix for a given band and spin.
     *
     * This method calculates the density matrix for a given band and spin using the wave function coefficients.
     * It adjusts the coefficients based on the Fermi energy and performs a matrix multiplication to produce the density
     * matrix.
     *
     * @param ib Band index.
     * @param nspin Number of spin channels.
     * @param nelec Total number of electrons.
     * @param nlocal Number of local orbitals.
     * @param wg Weight matrix for bands and spins.
     * @param DM Density matrix to be calculated.
     */
    void idmatrix(const int& ib,
                  const int nspin,
                  const double& nelec,
                  const int nlocal,
                  const ModuleBase::matrix& wg,
                  elecstate::DensityMatrix<double, double>& DM);
#endif
    psi::Psi<double>* psi_gamma;
    const Parallel_Orbitals* ParaV;
};
#endif
