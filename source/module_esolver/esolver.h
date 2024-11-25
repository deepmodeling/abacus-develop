#ifndef ESOLVER_H
#define ESOLVER_H

#include "module_base/matrix.h"
#include "module_cell/unitcell.h"
#include "module_parameter/parameter.h"

namespace ModuleESolver
{
class ESolver
{
  public:
    ESolver()
    {
        classname = "ESolver";
    }

    virtual ~ESolver()
    {
    }

    //! initialize the energy solver by using input parameters and cell modules
    virtual void before_all_runners(UnitCell& ucell, const Input_para& inp) = 0;

    //! run energy solver
    virtual void runner(UnitCell& cell, const int istep) = 0;

    //! perform post processing calculations
    virtual void after_all_runners(UnitCell& ucell){};

    //! deal with exx and other calculation than scf/md/relax/cell-relax:
    //! such as nscf, get_wf and get_pchg
    virtual void others(UnitCell& ucell, const int istep) {};

    //! calculate total energy of a given system
    virtual double cal_energy() = 0;

    //! calcualte forces for the atoms in the given cell
    virtual void cal_force(UnitCell& ucell, ModuleBase::matrix& force) = 0;

    //! calcualte stress of given cell
    virtual void cal_stress(UnitCell& ucell, ModuleBase::matrix& stress) = 0;

    bool conv_esolver = true; // whether esolver is converged

    std::string classname;
};

/**
 * @brief A subrutine called in init_esolver()
 *        This function returns type of ESolver
 *        Based on PARAM.inp.basis_type and PARAM.inp.esolver_type
 * 
 * @return [out] std::string The type of ESolver
 */
std::string determine_type();

/**
 * @brief Determine and initialize an ESolver based on input information.
 *
 * This function determines the type of ESolver to create based on input information and initializes
 * the corresponding ESolver child class. It supports various ESolver types including ksdft_pw,
 * ksdft_lcao, ksdft_lcao_tddft, sdft_pw, ofdft, lj_pot, and dp_pot.
 *
 * @return [out] A pointer to an ESolver object that will be initialized.
 */
ESolver* init_esolver(const Input_para& inp, UnitCell& ucell);

void clean_esolver(ESolver*& pesolver, const bool lcao_cblacs_exit = false);

} // namespace ModuleESolver

#endif
