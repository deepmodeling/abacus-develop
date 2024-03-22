#ifndef OUTPUT_INFO_H
#define OUTPUT_INFO_H

#include "module_cell/module_symmetry/symmetry.h"
#include "module_cell/atom_spec.h"
#include "module_cell/unitcell.h"
#include "module_base/matrix.h"


/**
* @brief In this part of the code to complete the output part of the json tree.
* @param ucell: ucell for reading json parameters
*/
namespace Json
{
#ifdef __RAPIDJSON

    void init_output_array_obj();

    void add_output_cell_coo_stress_force(
        UnitCell *ucell,
        ModuleBase::matrix force, double fac,
        ModuleBase::matrix stress, double unit_transform
    );

    void add_output_efermi_energy_converge(double efermi, double energy ,bool scf_converge );

    void add_output_scf_mag(
        double total_mag, double absolute_mag,
        double energy, double ediff, double drho,double time
    );

#endif
}
#endif