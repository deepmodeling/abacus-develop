#ifndef OUTPUT_LOG
#define OUTPUT_LOG
#include <fstream>
#include "module_base/global_variable.h"

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

} // namespace ModuleIO

#endif