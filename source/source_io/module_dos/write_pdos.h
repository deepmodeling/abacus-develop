#ifndef WRITE_PDOS_H
#define WRITE_PDOS_H

#include <map>
#include <string>
#include <utility>
#include <vector>

class UnitCell;
namespace ModuleBase { class matrix; }

namespace ModuleIO
{
// Each (l, real-harmonic column) maps to the zeta sum, ordered by spin and energy.
using PdosBlock = std::map<std::pair<int, int>, std::vector<double>>;
PdosBlock atom_pdos(const UnitCell& cell, int iat, const ModuleBase::matrix* pdos, int spins, bool spinor);

std::string write_pdos_npy(const UnitCell& cell, const std::vector<ModuleBase::matrix>& pdos,
                           int nspin, double emin, double de, const std::string& out_dir,
                           int istep, int out_app_flag);
}
#endif
