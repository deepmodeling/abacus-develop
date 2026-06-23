#ifndef ELECSTATE_H
#define ELECSTATE_H

#include "fp_energy.h"
#include "source_base/matrix.h"
#include "source_psi/psi.h" // Psi<T> appears in virtual method signatures; psi.h is light (own_fan=2)
#include "source_io/module_parameter/parameter.h" // PARAM/Input_para; kept (own_fan=4, cheap) to avoid churning ~120 downstream PARAM users

// Base-layer leaf utilities (own_fan~0) that the old heavy includes used to
// pull in transitively. Kept here so the ~100 downstream files using TITLE /
// WARNING_QUIT / timer / GlobalFunc / Memory do not each need a new include;
// global_function.h also transitively provides tool_title.h and tool_quit.h.
#include "source_base/global_function.h"
#include "source_base/timer.h"
#include "source_base/memory_recorder.h"

#include <complex>
#include <string>
#include <vector>

// The heavy data headers (potential_new own_fan=58, klist=42, charge=37) are
// forward-declared instead of included; their full definitions are pulled in by
// elecstate.cpp. parameter.h is kept above because dropping it saves only 4 fan
// but would force ~120 downstream files to add their own include.
class Charge;
class K_Vectors;
class UnitCell;
class Parallel_Grid;
namespace ModuleBase
{
class ComplexMatrix;
}
namespace ModulePW
{
class PW_Basis;
class PW_Basis_Big;
} // namespace ModulePW

namespace elecstate
{

class Potential;

class ElecState
{
  public:
    ElecState()
    {
    }
    ElecState(Charge* chr_in, ModulePW::PW_Basis* rhopw_in, ModulePW::PW_Basis_Big* bigpw_in);
    virtual ~ElecState();
    void init_ks(Charge* chr_in, // pointer for class Charge
                 const K_Vectors* klist_in,
                 int nk_in, // number of k points
                 const ModulePW::PW_Basis_Big* bigpw_in);

    // return current electronic density rho, as a input for constructing Hamiltonian
    virtual const double* getRho(int spin) const;

    // calculate electronic charge density on grid points or density matrix in real space
    // the consequence charge density rho saved into rho_out, preparing for charge mixing.
    virtual void psiToRho(const psi::Psi<std::complex<double>>& psi)
    {
        return;
    }
    virtual void psiToRho(const psi::Psi<double>& psi)
    {
        return;
    }
    virtual void cal_tau(const psi::Psi<std::complex<double>>& psi)
    {
        return;
    }
    virtual void cal_tau(const psi::Psi<double>& psi)
    {
        return;
    }
    virtual void cal_tau(const psi::Psi<std::complex<float>>& psi)
    {
        return;
    }

    // update charge density for next scf step
    // in this function, 1. input rho for construct Hamilt and 2. calculated rho from Psi will mix to 3. new charge
    // density rho among these rho,
    // 1. input rho would be store to file for restart
    // 2. calculated rho should be near with input rho when convergence has achieved
    // 3. new rho should be input rho for next scf step.
    virtual void getNewRho()
    {
        return;
    }

    // use occupied weights from INPUT and skip calculate_weights
    // mohan updated on 2024-06-08
    
    // if nupdown is not 0(TWO_EFERMI case),
    // nelec_spin will be fixed and weights will be constrained
    void init_nelec_spin();
    // used to record number of electrons per spin index
    // for NSPIN=2, it will record number of spin up and number of spin down
    // for NSPIN=4, it will record total number, magnetization for x, y, z direction
    std::vector<double> nelec_spin;

    virtual void print_psi(const psi::Psi<double>& psi_in, const int istep = -1)
    {
        return;
    }
    virtual void print_psi(const psi::Psi<std::complex<double>>& psi_in, const int istep = -1)
    {
        return;
    }



    std::string classname = "elecstate";

    int iter = 0;                                  ///< scf iteration
    Potential* pot = nullptr;                      ///< pointer to potential
    Charge* charge = nullptr;                      ///< pointer to charge density
    const K_Vectors* klist = nullptr;              ///< pointer to k points lists
    const ModulePW::PW_Basis_Big* bigpw = nullptr; ///< bigpw will be removed later

  public: // something aboud energies. See elecstate_energy.cpp
    void cal_bandgap();
    void cal_bandgap_updw();

    double cal_delta_eband(const UnitCell& ucell) const;
    double cal_delta_escf() const;

    ModuleBase::matrix vnew;
    bool vnew_exist = false;
    void cal_converged();
    void cal_energies(const int type);
    void set_exx(const double& Eexx);
    void set_exx(const std::complex<double>& Eexx);

    double get_hartree_energy();
    double get_etot_efield();
    double get_etot_gatefield();

    double get_solvent_model_Ael();
    double get_solvent_model_Acav();

    virtual double get_spin_constrain_energy()
    {
        return 0.0;
    }

    double get_dftu_energy();
    double get_local_pp_energy();

    fenergy f_en; ///< energies contribute to the total free energy
    Efermi eferm; ///< fermi energies

    // below defines the bandgap:

    double bandgap = 0.0;    ///< bandgap = E_{lumo} - E_{homo}
    double bandgap_up = 0.0; ///< spin up bandgap
    double bandgap_dw = 0.0; ///< spin down bandgap

    ModuleBase::matrix ekb; ///< band energy at each k point, each band.
    ModuleBase::matrix wg;  ///< occupation weight for each k-point and band

  public:

    bool skip_weights = false;
};

/**
 * @brief Init rho_core, init rho, renormalize rho, init pot
 *
 * @param ucell unit cell
 * @param pgrid parallel grid
 * @param strucfac structure factor
 * @param numeric numeric flag
 * @param istep ionic step index
 * @param out_dir output directory
 * @param inp input parameters
 * @param pelec pointer to ElecState
 */
void init_scf(const UnitCell& ucell,
              const Parallel_Grid& pgrid,
              const ModuleBase::ComplexMatrix& strucfac,
              const bool* numeric,
              const int istep,
              const std::string& out_dir,
              const Input_para& inp,
              ElecState* pelec);

} // namespace elecstate
#endif
