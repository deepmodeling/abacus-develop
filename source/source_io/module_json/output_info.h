#ifndef OUTPUT_INFO_H
#define OUTPUT_INFO_H

class UnitCell;

namespace ModuleBase
{
class matrix;
}

/**
 * @brief Generate the output section of the JSON tree.
 */
namespace Json
{
#ifdef __JSON

void init_output_array_obj();

void add_output_cell_coo_stress_force(const UnitCell& ucell,
                                      const ModuleBase::matrix& force,
                                      double fac,
                                      const ModuleBase::matrix& stress,
                                      double unit_transform,
                                      bool cal_force,
                                      bool cal_stress);

void add_output_efermi_converge(double efermi, bool scf_converge);
void add_output_energy(double energy);

void add_output_scf_mag(double total_mag,
                        double absolute_mag,
                        double energy,
                        double ediff,
                        double drho,
                        double time);

#endif // __JSON
} // namespace Json

#endif
