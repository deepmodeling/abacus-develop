#ifndef EXPERIMENTAL_GAUSSIAN_H
#define EXPERIMENTAL_GAUSSIAN_H
#include <vector>
class UnitCell;
namespace ModulePW { class PW_Basis; }
namespace ModuleBase { class matrix; }
namespace ModuleSccs {
std::vector<double> gaussian_ionic_density(const UnitCell&, const ModulePW::PW_Basis&, double);
ModuleBase::matrix gaussian_ionic_force(const UnitCell&, const ModulePW::PW_Basis&, double,
                                       const std::vector<double>&);
}
#endif
