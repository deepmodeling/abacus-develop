#ifndef OUTPUT_LOG
#define OUTPUT_LOG

#include <fstream>

#include "module_base/global_variable.h"
#include "module_base/matrix.h"
#include "module_cell/unitcell.h"

namespace ModuleIO
{

/// @brief output if is convergence and energy after scf
/// @param convergence if is convergence
/// @param energy the total energy in Ry
/// @param ofs_running the output stream
void output_convergence_after_scf(bool& convergence, double& energy, std::ofstream& ofs_running = GlobalV::ofs_running);

/// @brief output the fermi energy
/// @param convergence if is convergence
/// @param efermi
/// @param ofs_running the output stream
void output_efermi(bool& convergence, double& efermi, std::ofstream& ofs_running = GlobalV::ofs_running);

/// @brief output atomic forces
/// @param ofs the output stream
/// @param cell the unitcell
/// @param name force term name
/// @param force atomic forces
/// @param ry true if the unit of force is a.u.
void print_force(std::ofstream& ofs_running,
                 const UnitCell& cell,
                 const std::string& name,
                 const ModuleBase::matrix& force,
                 bool ry = true);

/// @brief output stress components
/// @param name stress term name
/// @param f stress components
/// @param ry true if the unit of force is a.u.
void print_stress(const std::string& name, const ModuleBase::matrix& f, const bool screen, bool ry);

/// @brief output total stress
/// @param scs total stress
/// @param ry true if the unit of force is a.u.
void printstress_total(const ModuleBase::matrix& scs, bool ry);

} // namespace ModuleIO

#endif