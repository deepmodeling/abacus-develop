#ifndef ESOLVER_FP_H
#define ESOLVER_FP_H

#include "esolver.h"

#ifndef __MPI
#include <chrono>
#endif

//! plane wave basis
#include "source_basis/module_pw/pw_basis.h"

//! symmetry analysis
#include "source_cell/module_symmetry/symmetry.h"

//! electronic states
#include "source_estate/elecstate.h"

//! charge extrapolation
#include "source_estate/module_charge/charge_extra.h"

//! solvation model
#include "source_hamilt/module_surchem/surchem.h"

//! local pseudopotential
#include "source_pw/hamilt_pwdft/VL_in_pw.h"

//! structure factor related to plane wave basis
#include "source_pw/hamilt_pwdft/structure_factor.h"

#include <fstream>


//! The First-Principles (FP) Energy Solver Class
/**
 * This class represents components that needed in
 * first-principles energy solver, such as the plane
 * wave basis, the structure factors, and the k points.
 *
 */

namespace ModuleESolver
{
class ESolver_FP: public ESolver
{
  public:
    //! Constructor
    ESolver_FP();

    //! Deconstructor
    virtual ~ESolver_FP();

    //! Initialize of the first-principels energy solver
    virtual void before_all_runners(UnitCell& ucell, const Input_para& inp) override;

    virtual void after_all_runners(UnitCell& ucell) override;

  protected:
    //! Something to do before SCF iterations.
    virtual void before_scf(UnitCell& ucell, const int istep);

    //! Something to do after SCF iterations when SCF is converged or comes to the max iter step.
    virtual void after_scf(UnitCell& ucell, const int istep, const bool conv_esolver);

    //! Something to do after hamilt2rho function in each iter loop.
    virtual void iter_finish(UnitCell& ucell, const int istep, int& iter, bool &conv_esolver);

    //! ------------------------------------------------------------------------------
    //! These pointers will be deleted in the free_pointers() function every ion step.
    //! ------------------------------------------------------------------------------
    elecstate::ElecState* pelec = nullptr; ///< Electronic states

    //! K points in Brillouin zone
    K_Vectors kv;

    //! Electorn charge density
    Charge chr;

    //! pw_rho: Plane-wave basis set for charge density
    //! pw_rhod: same as pw_rho for NCPP. Here 'd' stands for 'dense',
    //!          dense grid for for uspp, used for ultrasoft augmented charge density.
    //!          charge density and potential are defined on dense grids,
    //!          but effective potential needs to be interpolated on smooth grids in order to compute Veff|psi>
    ModulePW::PW_Basis* pw_rho;
    ModulePW::PW_Basis* pw_rhod;    //! dense grid for USPP
    ModulePW::PW_Basis_Big* pw_big; ///< [temp] pw_basis_big class

    //! parallel for rho grid
    Parallel_Grid Pgrid;

    //! Structure factors that used with plane-wave basis set
    Structure_Factor sf;

    //! local pseudopotentials
    pseudopot_cell_vl locpp;

    //! charge extrapolation method
    Charge_Extra CE;

    //! solvent model
    surchem solvent;

    int pw_rho_flag  = false; ///< flag for pw_rho, 0: not initialized, 1: initialized

    //! the start time of scf iteration
    #ifdef __MPI
        double iter_time;
    #else
        std::chrono::system_clock::time_point iter_time;
    #endif
};
} // namespace ModuleESolver

#endif
