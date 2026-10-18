#include "kernel_pots.h"
#include "potentials/pot_hxc_lrtd.h"
#include "utils/lr_util.h"
#include "source_io/module_parameter/input_parameter.h"
#include "source_cell/unitcell.h"
namespace LR
{
KernelPotentials prepare_lr_potentials(const Input_para& input, const std::string& xc_kernel,
    const ModulePW::PW_Basis& rho_basis, const UnitCell& ucell, const Charge& chg_gs,
    const Parallel_Grid& grid, const int nspin, const bool openshell)
{
    using ST = PotHxcLR::SpinType;
    using GX = LR::KernelXC::GxcSpin;
    KernelPotentials result;
    result.pot.resize(nspin);

    // The singlet and triplet potentials evaluate the *same* kernel arrays and differ only in which
    // spin combination of them they read, so they share one `KernelXC` to save memory.
    const bool oshell = (nspin == 2) && openshell;
    // $g^{xc}$ (third-order) is only ever needed by the LR gradient, and only for the spin
    // combinations that are actually going to be requested.
    const int gxc_lr = (!input.cal_force || !LR_Util::has_local_xc(xc_kernel)) ? GX::NoGxc
        : ((nspin == 1) ? GX::Singlet : GX::BothSpins);

    result.kernel = LR::build_xc_kernel(
        xc_kernel, rho_basis, ucell, chg_gs, grid, nspin, oshell, gxc_lr, input.lr_init_xc_kernel);
    switch (nspin)
    {
    case 1:
        result.pot[0] = std::make_shared<PotHxcLR>(result.kernel, rho_basis, ucell, chg_gs.nrxx, nspin, ST::S1);
        break;
    case 2:
    {
        const ST first_spin = oshell ? ST::S2_updown : ST::S2_singlet;
        const ST second_spin = oshell ? ST::S2_updown : ST::S2_triplet;
        result.pot[0] = std::make_shared<PotHxcLR>(result.kernel, rho_basis, ucell, chg_gs.nrxx, nspin, first_spin);
        result.pot[1] = std::make_shared<PotHxcLR>(result.kernel, rho_basis, ucell, chg_gs.nrxx, nspin, second_spin);
        break;
    }
    default:
        throw std::invalid_argument("ESolver_LR: nspin must be 1 or 2");
    }
    result.gxc = gxc_lr;
    return result;
}
std::shared_ptr<PotHxcLR> prepare_gs_potential(const Input_para& input, const std::string& xc_kernel,
    const ModulePW::PW_Basis& rho_basis, const UnitCell& ucell, const Charge& chg_gs,
    const Parallel_Grid& grid, const int nspin, const bool openshell,
    const KernelPotentials& prepared)
{
    using ST = PotHxcLR::SpinType;
    using GX = KernelXC::GxcSpin;
    const bool oshell = nspin == 2 && openshell;
    // `dft_functional == "default"` leaves the raw INPUT string unresolved (it never gets
    // overwritten to the actual functional in use); the functional actually read from the
    // pseudopotential lives in `ucell.atoms[i].ncpp.xc_func` instead. Comparing against the
    // literal "default" string here would always disagree with `xc_kernel`, even when the
    // ground state and LR use the same functional. Resolve it before comparing the names
    // or constructing a separate GS kernel.
    const std::string xc_kernel_gs = (input.dft_functional == "default")
        ? LR_Util::tolower(ucell.atoms[0].ncpp.xc_func)
        : LR_Util::tolower(input.dft_functional);
    // `ST::S1` is only correct when nspin=1. `KernelXC` is built with
    // the input `nspin`, so at nspin=2 the kernel arrays carry 3 spin components per grid point
    // while the S1 integrand indexes them as if there were 1 -- it does not even read a
    // consistent spin combination. Use `ST::S2_gs` there, which is exactly half of S2_singlet,
    // matching the `K_Hxc(singlet) = 2 * pot_hxc_gs` convention of the gradient operators.
    const ST st_gs = (nspin == 1) ? ST::S1_gs : (oshell ? ST::S2_updown : ST::S2_gs);
    // GS response terms retain this kernel. The g^xc terms originating from the LR K
    // contribution to W and the force use the corresponding LR potential instead.
    const int gxc_gs = (!LR_Util::has_local_xc(xc_kernel) || !LR_Util::has_local_xc(xc_kernel_gs)) ? GX::NoGxc
        : ((nspin == 1) ? GX::Singlet : GX::BothSpins);
    // When the LR kernel *is* the ground-state functional -- the usual TDDFT case -- the two
    // `KernelXC` are bit-for-bit identical, so reuse the one just built.
    const bool share_lr = (xc_kernel_gs == xc_kernel) && ((prepared.gxc & gxc_gs) == gxc_gs);
    std::shared_ptr<const LR::KernelXC> kernel_gs = share_lr ? prepared.kernel
        : LR::build_xc_kernel(xc_kernel_gs, rho_basis, ucell, chg_gs, grid, nspin, oshell, gxc_gs, input.lr_init_xc_kernel);
    return std::make_shared<LR::PotHxcLR>(kernel_gs, rho_basis, ucell, chg_gs.nrxx, nspin, st_gs);
}
}
