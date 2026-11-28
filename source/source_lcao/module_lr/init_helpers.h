#ifndef ABACUS_LR_INIT_HELPERS_H
#define ABACUS_LR_INIT_HELPERS_H
#include <iosfwd>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
namespace ModuleBase
{
class matrix;
}
struct Input_para;
class TwoCenterBundle;
class LCAO_Orbitals;
class UnitCell;
#ifdef __EXX
template <typename T>
class Exx_LRI;
class K_Vectors;
#endif
namespace LR
{
/// Initialize nonlocal projector metadata used by neighbor searches, independently of optical output.
void setup_nonlocal_info(LCAO_Orbitals&, UnitCell&, const Input_para&, int, std::ofstream&);
/// Build beta integrals and tabulate two-center integrals after setup_nonlocal_info.
void setup_2center_table(TwoCenterBundle&, UnitCell&, const Input_para&);
/// Validate LR solver, functional, spin count and optical gauge without global input access.
void check_lr_parameters(const Input_para&, const std::string& xc_kernel, int nspin);
#ifdef __EXX
/// Share the ground-state Exx_LRI as an LR workspace (the KS solver keeps using the original).
template <typename T>
inline std::shared_ptr<Exx_LRI<T>> make_lr_exx(const std::shared_ptr<Exx_LRI<T>>& exx_ks, const UnitCell& ucell,
                                               const K_Vectors& kv)
{
    return exx_ks->make_lr_workspace(ucell, kv);
}
/// Mismatched Exx_LRI types cannot back this LR solver.
template <typename T, typename U, typename std::enable_if<!std::is_same<T, U>::value, int>::type = 0>
inline std::shared_ptr<Exx_LRI<T>> make_lr_exx(const std::shared_ptr<Exx_LRI<U>>&, const UnitCell&, const K_Vectors&)
{
    throw std::runtime_error("LR::make_lr_exx: cannot share different Exx_LRI types");
}
#endif
}
#endif
