#include "esolver_fp.h"
#include "electronic_control.h"
#include "electronic_occupations.h"
#include "source_estate/param_update.h"
#include "source_estate/module_pot/electrostatic_potential.h"
#include "source_io/module_parameter/parameter.h"
#include "source_hamilt/module_xc/xc_functional.h"

#include <cmath>
#include <stdexcept>

namespace ModuleESolver
{
bool ESolver_FP::apply_electronic_control(Parameter& parameters, const ElectronicControl& control)
{
    const auto& inp = parameters.inp;
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
    if (inp.dft_plus_u && (inp.occ_mat_ctrl < 0 || inp.occ_mat_ctrl > 1
                          || inp.yukawa_potential || parameters.globalv.uramping > 0.0))
    {
        throw std::runtime_error(
            "constant-potential DFT+U requires fixed Hubbard U and self-consistent occupations "
            "(omc=0/1, yukawa_potential=false, uramping disabled)");
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
    if (!pelec || !pelec->charge || !pelec->pot)
    { throw std::runtime_error("electronic state is not initialized for electron control"); }

    elecstate::ParamUpdater::update_electron_count(control.nelect, parameters);
    pelec->eferm.two_efermi = two_fermi;
    pelec->init_nelec_spin();
    return true;
}

bool ESolver_FP::electronic_observables(const UnitCell& ucell, const Parameter& parameters,
                                       const PotentialWindow& window, ElectronicObservables& result) const
{
    if (!this->conv_esolver || !pelec || !pelec->charge || !pelec->pot)
    { throw std::runtime_error("electronic observables require a converged SCF state"); }
    const int npool = parameters.inp.kpar * parameters.inp.bndpar;
    result.nelect = checked_occupation_electron_count(pelec->wg, parameters.inp.nelec,
                                                      npool, pw_rhod->poolnproc);
    result.fermi_ry = pelec->eferm.ef;
    if (parameters.globalv.two_fermi)
    {
        // The total-N derivative at fixed N_up-N_down is the mean spin chemical potential.
        result.fermi_ry = 0.5 * (pelec->eferm.ef_up + pelec->eferm.ef_dw);
    }
    if (!std::isfinite(result.fermi_ry)) { throw std::runtime_error("non-finite SCF chemical potential"); }
    if (window.axis >= 0)
    {
        const double* fixed = pw_rhod->nrxx > 0 ? pelec->pot->get_fixed_v() : nullptr;
        const auto* pcc_potential = pelec->pot->pcc_electron_potential();
        const auto potential = elecstate::electrostatic_potential(ucell, *pw_rhod, *pelec->charge,
                                                                 fixed, solvent, parameters.inp, pcc_potential);
        result.reference_ry = elecstate::planar_potential_average(ucell, *pw_rhod, potential, window.axis,
                                                                 window.minimum_A, window.maximum_A);
    }
    return true;
}
}
