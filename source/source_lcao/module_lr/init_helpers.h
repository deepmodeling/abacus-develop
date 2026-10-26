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
/// Initialize nonlocal projector metadata used by neighbor searches, independently of optical output.
void setup_nonlocal_info(LCAO_Orbitals&, UnitCell&, const Input_para&, int, std::ofstream&);
/// Build beta integrals and tabulate two-center integrals after setup_nonlocal_info.
void setup_2center_table(TwoCenterBundle&, UnitCell&, const Input_para&);
/// Validate LR solver, functional, spin count and optical gauge without global input access.
void check_lr_parameters(const Input_para&, const std::string& xc_kernel, int nspin);
}
#endif
