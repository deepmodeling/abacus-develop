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
        const UnitCell *ucell,
        const ModuleBase::matrix force, const double fac,
        const ModuleBase::matrix stress, const double unit_transform
    );

    void add_output_efermi_converge(const double efermi, const bool scf_converge );
    void add_output_energy(const double energy );

    void add_output_scf_mag(
        const double total_mag, const double absolute_mag,
        const double energy, const double ediff, const double drho,const double time
    );

#endif
}
#endif