#ifndef ABACUS_LR_KERNEL_POTS_H
#define ABACUS_LR_KERNEL_POTS_H
#include <memory>
#include <string>
#include <vector>
struct Input_para;
class UnitCell;
class Charge;
class Parallel_Grid;
namespace ModulePW
{
class PW_Basis;
}
namespace LR
{
class PotHxcLR;
class KernelXC;
/// Own the LR spin potentials and the shared coefficient set prepared for this geometry.
struct KernelPotentials
{
    std::vector<std::shared_ptr<PotHxcLR>> pot;
    std::shared_ptr<const KernelXC> kernel;
    int gxc = 0;
};
/// Prepare LR potential views. At nspin=1, lr_target_spin selects a single channel;
/// triplets use a polarized kernel evaluated at rho_up=rho_down=rho/2.
/// nspin=2 retains both closed-shell views or the unrestricted spin blocks.
KernelPotentials prepare_lr_potentials(const Input_para&, const std::string& xc_kernel,
    const ModulePW::PW_Basis&, const UnitCell&, const Charge&, const Parallel_Grid&, int nspin, bool openshell);
/// Reuse LR coefficients for an identical GS functional and sufficient derivative order.
std::shared_ptr<PotHxcLR> prepare_gs_potential(const Input_para&, const std::string& xc_kernel,
    const ModulePW::PW_Basis&, const UnitCell&, const Charge&, const Parallel_Grid&, int nspin, bool openshell,
    const KernelPotentials&);
}
#endif
