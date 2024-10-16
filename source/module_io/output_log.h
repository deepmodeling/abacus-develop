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

/// @brief output after relaxation
/// @param conv_ion if is convergence for ions
/// @param conv_esolver if is convergence for electrons
/// @param ofs_running the output stream
void output_after_relax(bool conv_ion, bool conv_esolver, std::ofstream& ofs_running = GlobalV::ofs_running);

/// @brief output the fermi energy
/// @param convergence if is convergence
/// @param efermi
/// @param ofs_running the output stream
void output_efermi(bool& convergence, double& efermi, std::ofstream& ofs_running = GlobalV::ofs_running);

/// @brief calculate and output the vacuum level
/// We first determine the vacuum direction, then get the vacuum position based on the minimum of charge density,
/// finally output the vacuum level, i.e., the electrostatic potential at the vacuum position.
///
/// @param ucell the unitcell
/// @param rho charge density
/// @param v_elecstat electrostatic potential
/// @param ofs_running the output stream
void output_vacuum_level(const UnitCell* ucell,
                         const double* const* rho,
                         const double* v_elecstat,
                         const int& nx,
                         const int& ny,
                         const int& nz,
                         const int& nxyz,
                         const int& nrxx,
                         const int& nplane,
                         const int& startz_current,
                         std::ofstream& ofs_running = GlobalV::ofs_running);

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
void print_stress(const std::string& name, const ModuleBase::matrix& scs, const bool screen, const bool ry);

/// @brief write head for scf iteration
/// @param ofs_running output stream
/// @param istep the ion step
/// @param iter the scf iteration step
/// @param basisname basis set name
void write_head(std::ofstream& ofs_running, const int& istep, const int& iter, const std::string& basisname);

} // namespace ModuleIO

#endif