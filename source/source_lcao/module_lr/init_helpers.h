#ifndef ABACUS_LR_INIT_HELPERS_H
#define ABACUS_LR_INIT_HELPERS_H
#include <iosfwd>
#include <string>
namespace ModuleBase
{
class matrix;
}
struct Input_para;
class TwoCenterBundle;
class LCAO_Orbitals;
class UnitCell;
namespace LR
{
/// Prepare nonlocal beta projectors before tabulating two-center integrals.
void setup_2center_table(TwoCenterBundle&, LCAO_Orbitals&, UnitCell&, const Input_para&, int, std::ofstream&);
/// Validate LR solver, functional, spin count and optical gauge without global input access.
void check_lr_parameters(const Input_para&, const std::string& xc_kernel, int nspin);
}
#endif
