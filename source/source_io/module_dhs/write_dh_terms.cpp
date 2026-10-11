#include "source_base/timer.h"
#include "source_io/module_hs/hsr_writer.h"
#include "source_cell/ucell_io.h"
#include "source_io/module_parameter/parameter.h"
#include "source_hamilt/module_hcontainer/hcontainer_funcs.h"
#include "source_hamilt/module_hcontainer/output_hcontainer.h"
#include "source_lcao/module_operator_lcao/ekinetic.h"
#include "source_lcao/module_operator_lcao/nonlocal.h"
#include "source_lcao/module_operator_lcao/operator_fs_utils.h"
#include "source_lcao/module_operator_lcao/veff_lcao.h"
#include "source_hamilt/module_gint/gint_interface.h"
#include "source_lcao/module_lr/utils/lr_util_xc.hpp"
#include "source_base/global_variable.h"
#include "source_base/parallel_reduce.h"
#include "write_dh.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

namespace ModuleIO
{

namespace
{

    // RAII holder for the per-atom-I dH containers: one HContainer per (atom I, direction d).
    // g[d] are raw-pointer views into the owned containers, ready to hand to an
// operator's cal_dH(...) and then to write_dh_perI(...).
struct PerIContainers
{
    std::array<std::vector<std::unique_ptr<hamilt::HContainer<double>>>, 3> owned;
    std::array<std::vector<hamilt::HContainer<double>*>, 3> g;

    PerIContainers(const Parallel_Orbitals& pv, int nat)
    {
        for (int d = 0; d < 3; ++d)
        {
            owned[d].reserve(nat);
            g[d].reserve(nat);
            for (int iat = 0; iat < nat; ++iat)
            {
                owned[d].push_back(make_unique<hamilt::HContainer<double>>(&pv));
                g[d].push_back(owned[d].back().get());
            }
        }
    }
};

#ifdef __DEBUG
// Self-validation for cal_gint_drho (otherwise untested). For a symmetric DM the product rule
// gives  grad(rho) = sum_{K,L} D[K,L]( (grad phi_K) phi_L + phi_K (grad phi_L) )
//                  = 2 * sum_{K,L} D[K,L] (grad phi_K) phi_L = 2 * cal_gint_drho(D).
// We compare 2*cal_gint_drho(D) against the FFT gradient of cal_gint_rho(D) (same density,
// independent operator) and log the max relative deviation under a grep-able key so the
// integrate harness can assert it stays ~0. Cheap: one rho + one drho + one FFT grad.
void validate_gint_drho(const UnitCell& ucell,
                        elecstate::Potential* pot,
                        const hamilt::HContainer<double>* dmR)
{
    if (dmR == nullptr)
        return;
    const ModulePW::PW_Basis* rho_basis = pot->get_rho_basis();
    const int nrxx = rho_basis->nrxx;
    std::vector<hamilt::HContainer<double>*> dm_vec = {const_cast<hamilt::HContainer<double>*>(dmR)};

    // rho via Gint
    std::vector<double> rho(nrxx, 0.0);
    double* rho_p[1] = {rho.data()};
    ModuleGint::cal_gint_rho(dm_vec, 1, rho_p, false);

    // grad rho via FFT
    std::vector<ModuleBase::Vector3<double>> gradrho(nrxx);
    LR_Util::grad(rho.data(), gradrho.data(), *rho_basis, ucell.tpiba);

    // drho via Gint (gradient on the first/row orbital)
    std::vector<double> dx(nrxx, 0.0), dy(nrxx, 0.0), dz(nrxx, 0.0);
    double* dxp[1] = {dx.data()};
    double* dyp[1] = {dy.data()};
    double* dzp[1] = {dz.data()};
    ModuleGint::cal_gint_drho(dm_vec, 1, dxp, dyp, dzp);

    double maxdev = 0.0, maxref = 0.0;
    for (int ir = 0; ir < nrxx; ++ir)
    {
        const double r[3] = {2.0 * dx[ir], 2.0 * dy[ir], 2.0 * dz[ir]};
        const double g[3] = {gradrho[ir].x, gradrho[ir].y, gradrho[ir].z};
        for (int d = 0; d < 3; ++d)
        {
            maxdev = std::max(maxdev, std::abs(r[d] - g[d]));
            maxref = std::max(maxref, std::abs(g[d]));
        }
    }
#ifdef __MPI
    Parallel_Reduce::reduce_all(maxdev);
    Parallel_Reduce::reduce_all(maxref);
#endif
    const double reldev = (maxref > 1e-30) ? (maxdev / maxref) : maxdev;
    GlobalV::ofs_running << " GINT_DRHO_MAXDEV_REL " << reldev << std::endl;
}
#endif

// Per-(spin) fillers: build one term's per-atom-I dH containers, no file output. Shared by the
// individual term writers below and by write_dH_sum (which accumulates them). Keeping the
// operator construction in one place avoids duplicating it in the summation path.
void fill_dH_t(WriteDHParams& params, PerIContainers& c)
{
    const UnitCell& ucell = *params.ucell;
    const Grid_Driver& gd = *params.gd;
    const TwoCenterBundle& two_center_bundle = *params.two_center_bundle;
    const std::vector<double>& orb_cutoff = params.orb->cutoffs();

    hamilt::EKinetic<hamilt::OperatorLCAO<double, double>> tmp_ekinetic(nullptr,
                                                                        params.kv->kvec_d,
                                                                        nullptr,
                                                                        &ucell,
                                                                        orb_cutoff,
                                                                        &gd,
                                                                        two_center_bundle.kinetic_orb.get());

    tmp_ekinetic.cal_dH(c.g);
}

void fill_dH_vnl(WriteDHParams& params, PerIContainers& c)
{
    const UnitCell& ucell = *params.ucell;
    const Grid_Driver& gd = *params.gd;
    const TwoCenterBundle& two_center_bundle = *params.two_center_bundle;
    const std::vector<double>& orb_cutoff = params.orb->cutoffs();

    hamilt::Nonlocal<hamilt::OperatorLCAO<double, double>> tmp_nonlocal(nullptr,
                                                                        params.kv->kvec_d,
                                                                        nullptr,
                                                                        &ucell,
                                                                        orb_cutoff,
                                                                        &gd,
                                                                        two_center_bundle.overlap_orb_beta.get());

    tmp_nonlocal.cal_dH(c.g);
}

void fill_dH_veff(WriteDHParams& params,
                  elecstate::Potential* pot,
                  const std::string& hf_type,
                  int ispin,
                  PerIContainers& c)
{
    const UnitCell& ucell = *params.ucell;
    const Grid_Driver& gd = *params.gd;
    const Parallel_Orbitals& pv = *params.pv;
    const std::vector<double>& orb_cutoff = params.orb->cutoffs();
    const int nspin = params.nspin;

    hamilt::HContainer<double> hR_dummy(const_cast<Parallel_Orbitals*>(&pv));

    hamilt::Veff<hamilt::OperatorLCAO<double, double>> veff(nullptr,
                                                            params.kv->kvec_d,
                                                            pot,
                                                            &hR_dummy,
                                                            &ucell,
                                                            orb_cutoff,
                                                            &gd,
                                                            nspin);

    veff.cal_dH(c.g, hf_type, params.dmR, params.chg, ispin);
}

// Shared driver for the Veff-based terms (V^L, V^H, V^XC), which differ only in the
// Hellmann-Feynman type passed to Veff::cal_dH and in the output prefixes/label.
bool write_dH_veff_term(WriteDHParams& params,
                        elecstate::Potential* pot,
                        const std::string& hf_type,
                        const std::string& rprefix,
                        const std::string& kprefix,
                        const std::string& label,
                        const std::string& term_name,
                        const std::vector<int>& atom_filter = {})
{
    const UnitCell& ucell = *params.ucell;
    const Parallel_Orbitals& pv = *params.pv;
    const int nat = ucell.nat;
    const int nspin = params.nspin;

#ifdef __DEBUG
    // Validate cal_gint_drho once (it underpins the V^H Hellmann-Feynman term).
    if (hf_type == "hartree")
        validate_gint_drho(ucell, pot, params.dmR.empty() ? nullptr : params.dmR[0]);
#endif

    for (int ispin = 0; ispin < (nspin == 2 ? 2 : 1); ispin++)
    {
        PerIContainers c(pv, nat);

        fill_dH_veff(params, pot, hf_type, ispin, c);

        ModuleIO::write_dh_perI(params, ispin, rprefix, kprefix, label, term_name, c.g, atom_filter);
    }
    return true;
}

} // namespace

bool write_dH_t(WriteDHParams& params)
{
    ModuleBase::TITLE("ModuleIO", "write_dH_t");
    ModuleBase::timer::start("ModuleIO", "write_dH_t");

    const Parallel_Orbitals& pv = *params.pv;
    const int nat = params.ucell->nat;
    const int nspin = params.nspin;

    const std::vector<int> af_t = dh_atom_filter(PARAM.inp.out_mat_dh_t);
    for (int ispin = 0; ispin < (nspin == 2 ? 2 : 1); ispin++)
    {
        // per-atom-I containers: dT_*[iat] = d<phi|T|phi>/dtau_iat
        PerIContainers c(pv, nat);

        fill_dH_t(params, c);

        ModuleIO::write_dh_perI(params, ispin, "dtr", "dtk", "dT", "kinetic", c.g, af_t);
    }

    ModuleBase::timer::end("ModuleIO", "write_dH_t");
    return true;
}

bool write_dH_vnl(WriteDHParams& params)
{
    ModuleBase::TITLE("ModuleIO", "write_dH_vnl");
    ModuleBase::timer::start("ModuleIO", "write_dH_vnl");

    const Parallel_Orbitals& pv = *params.pv;
    const int nat = params.ucell->nat;
    const int nspin = params.nspin;

    const std::vector<int> af_vnl = dh_atom_filter(PARAM.inp.out_mat_dh_vnl);
    for (int ispin = 0; ispin < (nspin == 2 ? 2 : 1); ispin++)
    {
        PerIContainers c(pv, nat);

        fill_dH_vnl(params, c);

        ModuleIO::write_dh_perI(params, ispin, "dvnlr", "dvnlk", "dV^NL", "nonlocal", c.g, af_vnl);
    }

    ModuleBase::timer::end("ModuleIO", "write_dH_vnl");
    return true;
}

bool write_dH_vl(WriteDHParams& params)
{
    ModuleBase::TITLE("ModuleIO", "write_dH_vl");
    ModuleBase::timer::start("ModuleIO", "write_dH_vl");

    const std::vector<int> af_vl = dh_atom_filter(PARAM.inp.out_mat_dh_vl);
    const bool ok = write_dH_veff_term(params, params.pot_vl, "vl", "dvlr", "dvlk", "dV^L", "local", af_vl);

    ModuleBase::timer::end("ModuleIO", "write_dH_vl");
    return ok;
}

bool write_dH_vh(WriteDHParams& params)
{
    ModuleBase::TITLE("ModuleIO", "write_dH_vh");
    ModuleBase::timer::start("ModuleIO", "write_dH_vh");

    const std::vector<int> af_vh = dh_atom_filter(PARAM.inp.out_mat_dh_vh);
    const bool ok = write_dH_veff_term(params, params.pot_vh, "hartree", "dvhr", "dvhk", "dV^H", "hartree", af_vh);

    ModuleBase::timer::end("ModuleIO", "write_dH_vh");
    return ok;
}

bool write_dH_vh_pulay(WriteDHParams& params)
{
    ModuleBase::TITLE("ModuleIO", "write_dH_vh_pulay");
    ModuleBase::timer::start("ModuleIO", "write_dH_vh_pulay");

    const std::vector<int> af_vh_pulay = dh_atom_filter(PARAM.inp.out_mat_dh_vh);
    const bool ok = write_dH_veff_term(params, params.pot_vh, "none", "dvhr_pulay_", "dvhk_pulay_", "dV^H (Pulay)", "hartree Pulay", af_vh_pulay);

    ModuleBase::timer::end("ModuleIO", "write_dH_vh_pulay");
    return ok;
}

bool write_dH_vxc(WriteDHParams& params)
{
    ModuleBase::TITLE("ModuleIO", "write_dH_vxc");
    ModuleBase::timer::start("ModuleIO", "write_dH_vxc");

    const std::vector<int> af_vxc = dh_atom_filter(PARAM.inp.out_mat_dh_vxc);
    const bool ok = write_dH_veff_term(params, params.pot_vxc,
                                       params.chg ? "xc" : "none",
                                       "dvxcr", "dvxck", "dV^XC", "xc", af_vxc);

    ModuleBase::timer::end("ModuleIO", "write_dH_vxc");
    return ok;
}

bool write_dH_vxc_pulay(WriteDHParams& params)
{
    ModuleBase::TITLE("ModuleIO", "write_dH_vxc_pulay");
    ModuleBase::timer::start("ModuleIO", "write_dH_vxc_pulay");

    const std::vector<int> af_vxc_pulay = dh_atom_filter(PARAM.inp.out_mat_dh_vxc);
    const bool ok = write_dH_veff_term(params, params.pot_vxc, "none", "dvxcr_pulay_", "dvxck_pulay_", "dV^XC (Pulay)", "xc Pulay", af_vxc_pulay);

    ModuleBase::timer::end("ModuleIO", "write_dH_vxc_pulay");
    return ok;
}

// Total dH = sum of ALL dH terms (dT + dV^NL + dV^L + dV^H + dV^XC), independent of the
// per-component out_mat_dh_* flags: out_mat_dh on its own yields the full sum. Each term is
// built into its own per-atom-I containers (via the same fillers the per-term writers use) and
// accumulated with HContainer::add_value_union, which unions the (generally different)
// sparsities and sums values. Each term already carries its own sign.
//
// Note: the EXX (dV^EXX/dR) term is intentionally NOT included: it requires the LibRI
// cal_dHs/dHs_HF APIs that are only available in the unmerged LibRI PR#10 (the former
// EXX_DEV path), so it cannot be computed in the current develop branch. When EXX is active
// and out_mat_dh is requested at gamma-only, the resulting dH sum therefore omits the EXX
// contribution.
bool write_dH_sum(WriteDHParams& params)
{
    ModuleBase::TITLE("ModuleIO", "write_dH_sum");
    ModuleBase::timer::start("ModuleIO", "write_dH_sum");

    const Parallel_Orbitals& pv = *params.pv;
    const int nat = params.ucell->nat;
    const int nspin = params.nspin;

    for (int ispin = 0; ispin < (nspin == 2 ? 2 : 1); ++ispin)
    {
        PerIContainers sum(pv, nat);
        auto accumulate = [&](PerIContainers& term) {
            for (int d = 0; d < 3; ++d)
                for (int iat = 0; iat < nat; ++iat)
                    sum.g[d][iat]->add_value_union(*term.g[d][iat]);
        };

        // dT (kinetic) and dV^NL (nonlocal) need no potential.
        {
            PerIContainers c(pv, nat);
            fill_dH_t(params, c);
            accumulate(c);
        }
        {
            PerIContainers c(pv, nat);
            fill_dH_vnl(params, c);
            accumulate(c);
        }
        // veff terms (potentials are allocated whenever out_mat_dh is on; guard defensively).
        if (params.pot_vl != nullptr)
        {
            PerIContainers c(pv, nat);
            fill_dH_veff(params, params.pot_vl, "vl", ispin, c);
            accumulate(c);
        }
        if (params.pot_vh != nullptr)
        {
            // total dV^H = Hellmann-Feynman part + Pulay part
            PerIContainers c(pv, nat);
            fill_dH_veff(params, params.pot_vh, "hartree", ispin, c);
            accumulate(c);
        }
        if (params.pot_vxc != nullptr)
        {
            // total dV^XC = Hellmann-Feynman part + Pulay part
            PerIContainers c(pv, nat);
            fill_dH_veff(params, params.pot_vxc, params.chg ? "xc" : "none", ispin, c);
            accumulate(c);
        }

        ModuleIO::write_dh_perI(params, ispin, "dhr", "dhk", "dH", "total", sum.g, dh_atom_filter(PARAM.inp.out_mat_dh));
    }

    ModuleBase::timer::end("ModuleIO", "write_dH_sum");
    return true;
}

} // namespace ModuleIO
