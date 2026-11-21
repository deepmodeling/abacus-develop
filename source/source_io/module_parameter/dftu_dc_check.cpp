#include "dftu_dc_check.h"
#include "input_parameter.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <stdexcept>

namespace ModuleIO
{
double parse_dftu_alpha(const std::string& token)
{
    std::size_t consumed = 0;
    const double alpha = std::stod(token, &consumed);
    if (consumed != token.size() || !std::isfinite(alpha) || alpha < 0.0 || alpha > 1.0)
    {
        throw std::invalid_argument("dftu_alpha must be a finite number in [0,1]");
    }
    return alpha;
}

namespace
{
bool output_enabled(const std::vector<int>& setting)
{
    return !setting.empty() && setting.front() != 0;
}
} // namespace

std::string dftu_dc_error(const Input_para& input, const bool alpha_was_read, const int kpar_lcao)
{
    if (input.dftu_dc != "fll" && input.dftu_dc != "amf" && input.dftu_dc != "mixed")
    {
        return "dftu_dc must be fll, amf or mixed";
    }
    if (input.dftu_dc == "mixed")
    {
        if (!alpha_was_read || !std::isfinite(input.dftu_alpha)
            || input.dftu_alpha < 0.0 || input.dftu_alpha > 1.0)
        {
            return "dftu_dc mixed requires an explicit finite dftu_alpha in [0,1]";
        }
    }
    else if (alpha_was_read)
    {
        return "dftu_alpha is only used with dftu_dc mixed";
    }
    if (input.dftu_dc == "fll")
    {
        return {}; // No new restrictions on legacy inputs.
    }
    if (input.dft_plus_u != 1 || input.basis_type != "lcao")
    {
        return "AMF/mixed currently requires Full dft_plus_u=1 and basis_type=lcao";
    }
    if (input.device != "cpu" || input.precision != "double" || input.gint_precision != "double"
        || input.esolver_type != "ksdft" || input.calculation != "scf" || input.socket_driver
        || input.cal_force || input.cal_stress)
    {
        return "AMF/mixed requires double-precision CPU fixed-structure KS SCF without force/stress";
    }
    if ((input.nspin != 1 && input.nspin != 2) || input.lspinorb || input.noncolin
        || input.sc_mag_switch || input.nupdown != 0.0)
    {
        return "AMF/mixed requires collinear nspin=1/2 without spin constraints";
    }
    if (input.kpar != 1 || kpar_lcao != 1 || input.bndpar != 1
        || (input.symmetry != "0" && input.symmetry != "-1"))
    {
        return "AMF/mixed requires one k pool and symmetry=0 or -1";
    }
    if (input.ks_solver != "genelpa" && input.ks_solver != "elpa"
        && input.ks_solver != "lapack" && input.ks_solver != "scalapack_gvx")
    {
        return "AMF/mixed requires genelpa, elpa, lapack or scalapack_gvx";
    }
    if (input.init_chg != "atomic" || input.init_wfc != "atomic" || input.init_occ_mat != 0
        || input.omc != 0 || input.restart_load || input.restart_save || input.dm_to_rho)
    {
        return "AMF/mixed requires atomic cold start without occupation/restart files";
    }
    if (input.uramping_eV > 0.0 || input.yukawa_potential || input.mixing_dftu || input.mixing_dmr)
    {
        return "AMF/mixed does not yet support U ramping, Yukawa U or occupation/DMR mixing";
    }
    std::string functional = input.dft_functional;
    std::transform(functional.begin(), functional.end(), functional.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    if (functional != "pbe" || !input.exx_fock_alpha.empty() || !input.exx_erfc_alpha.empty()
        || input.rpa || input.ml_exx || input.deepks_scf || input.dft_plus_dmft)
    {
        return "AMF/mixed currently requires explicit PBE without EXX, RPA, DeePKS or DMFT";
    }
    const std::vector<const std::vector<int>*> matrix_outputs = {
        &input.out_hsk, &input.out_hsr, &input.out_mat_hs, &input.out_mat_hs2,
        &input.out_mat_dh, &input.out_dmr, &input.out_dmk, &input.out_mat_h_t,
        &input.out_mat_h_vnl, &input.out_mat_h_vl, &input.out_mat_h_vh, &input.out_mat_h_vxc,
        &input.out_mat_h_exx, &input.out_mat_dh_t, &input.out_mat_dh_vl, &input.out_mat_dh_vnl,
        &input.out_mat_dh_vh, &input.out_mat_dh_vxc, &input.out_mat_ds};
    for (const auto* setting : matrix_outputs)
    {
        if (output_enabled(*setting))
        {
            return "AMF/mixed Hamiltonian, derivative and density matrix exports are not supported yet";
        }
    }
    if (input.out_hr_npz || input.out_hsr_npz || input.out_hsr_npz_compat || input.out_dm_npz
        || input.out_eband_terms || input.towannier90 || input.qo_switch)
    {
        return "AMF/mixed matrix export, Wannier and QO combinations are not supported yet";
    }
    return {};
}
} // namespace ModuleIO
