#include "electronic_control.h"
#include "source_esolver/esolver.h"
#include "source_estate/electronic_control.h"
#include "source_io/module_parameter/parameter.h"
#include "source_hamilt/module_xc/xc_functional.h"

#include <cmath>
#include <stdexcept>

namespace relax
{
void run_electronic_step(ModuleESolver::ESolver& solver, UnitCell& cell,
                         Parameter& parameters, const elecstate::ElectronicControl& control,
                         int istep)
{
    const auto& inp = parameters.inp;
    if (inp.assume_isolated == "makov-payne")
    {
        throw std::runtime_error(
            "constant-potential electron control is incompatible with the energy-only Makov-Payne correction");
    }
    if (inp.esolver_type != "ksdft" || (inp.basis_type != "pw" && inp.basis_type != "lcao")
        || inp.device != "cpu" || inp.precision != "double" || parameters.globalv.use_uspp
        || inp.sc_mag_switch || inp.deepks_scf || inp.rpa || inp.rdmft || inp.dft_plus_dmft)
    {
        throw std::runtime_error(
            "constant-potential control requires CPU double-precision PW/LCAO KSDFT with "
            "norm-conserving pseudopotentials; specialized solver corrections are not supported");
    }
    const int functional_type = XC_Functional::get_func_type();
    if (functional_type < 1 || functional_type > 5)
    {
        throw std::runtime_error(
            "constant-potential control requires a supported KS exchange-correlation functional");
    }
    if (inp.dft_plus_u != 0 && inp.dft_plus_u != 1)
    {
        throw std::runtime_error(
            "constant-potential forces require dft_plus_u=1; "
            "the legacy dft_plus_u=2 force implementation is unavailable");
    }
    if (inp.dft_plus_u && (inp.init_occ_mat < 0 || inp.init_occ_mat > 1
                          || inp.yukawa_potential || parameters.globalv.uramping > 0.0))
    {
        throw std::runtime_error(
            "constant-potential DFT+U requires fixed Hubbard U and self-consistent occupations "
            "(init_occ_mat=0/1, yukawa_potential=false, uramping disabled)");
    }
    if (inp.init_chg == "dm_no_renormalize")
    { throw std::runtime_error("electron-number control requires charge-density renormalization"); }
    if (inp.nspin != 1 && inp.nspin != 2)
    { throw std::runtime_error("constant-potential control supports only nspin=1 or 2"); }
    if (inp.ocp || inp.nbands <= 0 || !std::isfinite(control.nelect) || control.nelect <= 0.0)
    { throw std::runtime_error("invalid total electron number, band capacity, or fixed occupation array"); }
    const bool two_fermi = parameters.globalv.two_fermi;
    if (two_fermi && inp.nspin != 2)
    { throw std::runtime_error("two_fermi electron control requires nspin=2"); }
    const double up = 0.5 * (control.nelect + inp.nupdown);
    const double down = 0.5 * (control.nelect - inp.nupdown);
    if ((two_fermi && (!std::isfinite(inp.nupdown) || up <= 0.0 || down <= 0.0
                      || up > inp.nbands || down > inp.nbands))
        || control.nelect > 2.0 * inp.nbands)
    {
        throw std::runtime_error(
            "electron number exceeds band capacity or gives invalid spin populations; increase nbands");
    }
    if (inp.smearing_method == "fixed"
        && (std::abs(control.nelect - std::round(control.nelect)) > 1.0e-8
            || (two_fermi && (std::abs(up-std::round(up)) > 1.0e-8 || std::abs(down-std::round(down)) > 1.0e-8))))
    { throw std::runtime_error("fractional electron populations require a smearing occupation method"); }
    if (!solver.apply_electronic_control(parameters, control))
    {
        throw std::runtime_error("selected ESolver does not support constant-potential electron control");
    }
    solver.runner(cell, istep);
}
}
