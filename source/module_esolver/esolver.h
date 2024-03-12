#ifndef ESOLVER_H
#define ESOLVER_H

#include "module_base/matrix.h"
#include "module_cell/unitcell.h"
#include "module_io/input.h"

namespace ModuleESolver
{
class ESolver
{
    // protected:
    //     ModuleBase::matrix lattice_v;
  public:
    ESolver()
    {
        classname = "ESolver";
    }

    virtual ~ESolver()
    {
    }

    // virtual void Init(Input_EnSolver &inp, matrix &lattice_v)=0
    virtual void Init(Input& inp, UnitCell& cell) = 0;

    // They shoud be add after atom class is refactored
    // virtual void UpdateLatAtom(ModuleBase::matrix &lat_in, Atom &atom_in);
    // virtual void UpdateLat(ModuleBase::matrix &lat_in);
    // virtual void UpdateAtom(Atom &atom_in);

    virtual void Run(int istep, UnitCell& cell) = 0;

    // Deal with exx and other calculation than scf/md/relax:
    //  such as nscf, get_wf and get_pchg
    virtual void othercalculation(const int istep){};

    virtual double cal_Energy() = 0;
    virtual void cal_Force(ModuleBase::matrix& force) = 0;
    virtual void cal_Stress(ModuleBase::matrix& stress) = 0;
    virtual void postprocess(){};

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
std::string determine_type();

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
