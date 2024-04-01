#ifndef ESOLVER_H
#define ESOLVER_H

#include "module_base/matrix.h"
#include "module_cell/unitcell.h"
#include "module_io/input.h"

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
    virtual void init(Input& inp, UnitCell& cell) = 0;

    //! run energy solver
    virtual void run(int istep, UnitCell& cell) = 0;

    //! deal with exx and other calculation than scf/md/relax:
    //! such as nscf, get_wf and get_pchg
    virtual void others(const int istep){};

    //! calculate total energy of a given system
    virtual double cal_energy() = 0;

    //! calcualte forces for the atoms in the given cell
    virtual void cal_force(ModuleBase::matrix& force) = 0;

    //! calcualte stress of given cell
    virtual void cal_stress(ModuleBase::matrix& stress) = 0;

    //! perform post processing calculations
    virtual void post_process(){};

    // Print current classname.
    void printname();

    // temporarily
    // get iterstep used in current scf
    virtual int getniter()
    {
        return 0;
    }
    std::string classname;
};

/**
 * @brief A subrutine called in init_esolver()
 *        This function returns type of ESolver
 *        Based on GlobalV::BASIS_TYPE and GlobalV::ESOLVER_TYPE
 * 
 * @return [out] std::string The type of ESolver
 */
std::string determine_type(void);

/**
 * @brief Determine and initialize an ESolver based on input information.
 *
 * This function determines the type of ESolver to create based on input information and initializes
 * the corresponding ESolver child class. It supports various ESolver types including ksdft_pw,
 * ksdft_lcao, ksdft_lcao_tddft, sdft_pw, ofdft, lj_pot, and dp_pot.
 *
 * @param [in, out] p_esolver A pointer to an ESolver object that will be initialized.
 */
void init_esolver(ESolver*& p_esolver);

void clean_esolver(ESolver*& pesolver);

} // namespace ModuleESolver

#endif
