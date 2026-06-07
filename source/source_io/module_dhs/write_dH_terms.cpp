#include "source_base/timer.h"
#include "source_io/module_hs/write_HS_R.h"
#include "source_io/module_output/ucell_io.h"
#include "source_io/module_parameter/parameter.h"
#include "source_lcao/module_hcontainer/hcontainer_funcs.h"
#include "source_lcao/module_hcontainer/output_hcontainer.h"
#include "source_lcao/module_operator_lcao/ekinetic.h"
#include "source_lcao/module_operator_lcao/nonlocal.h"
#include "source_lcao/module_operator_lcao/operator_force_stress_utils.h"
#include "source_lcao/module_operator_lcao/veff_lcao.h"
#include "write_dH.h"
#ifdef __EXX
#include "source_lcao/module_operator_lcao/op_exx_lcao.h"
#include "source_lcao/module_ri/Exx_LRI_interface.hpp"
#endif

#include <memory>
#include <string>

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
                owned[d].push_back(std::make_unique<hamilt::HContainer<double>>(&pv));
                g[d].push_back(owned[d].back().get());
            }
        }
    }
};

// Shared driver for the Veff-based terms (V^L, V^H, V^XC), which differ only in the
// Hellmann-Feynman type passed to Veff::cal_dH and in the output prefixes/label.
bool write_dH_veff_term(WriteDHParams& params,
                        elecstate::Potential* pot,
                        const std::string& hf_type,
                        const std::string& rprefix,
                        const std::string& kprefix,
                        const std::string& label)
{
    const UnitCell& ucell = *params.ucell;
    const Grid_Driver& gd = *params.gd;
    const Parallel_Orbitals& pv = *params.pv;
    const std::vector<double>& orb_cutoff = params.orb->cutoffs();
    const int nat = ucell.nat;
    const int nspin = params.nspin;

    for (int ispin = 0; ispin < (nspin == 2 ? 2 : 1); ispin++)
    {
        hamilt::HContainer<double> hR_dummy(const_cast<Parallel_Orbitals*>(&pv));

        hamilt::Veff<hamilt::OperatorLCAO<double, double>> veff(nullptr,
                                                                params.kv->kvec_d,
                                                                pot,
                                                                &hR_dummy,
                                                                &ucell,
                                                                orb_cutoff,
                                                                &gd,
                                                                nspin);

        PerIContainers c(pv, nat);

        veff.cal_dH(c.g, hf_type, params.dmR);

        ModuleIO::write_dh_perI(params, ispin, rprefix, kprefix, label, c.g);
    }
    return true;
}

#ifdef __EXX
// Shared driver for the EXX dH term, templated on the Hexx tensor data type (double for the
// real interface exd, std::complex<double> for the complex interface exc). The per-atom-I dH
// is always written into real HContainer<double> (add_HexxR converts Tdata -> double).
template <typename Tdata>
void write_dH_exx_impl(WriteDHParams& params, Exx_LRI_Interface<double, Tdata>* ex)
{
    const UnitCell& ucell = *params.ucell;
    const Parallel_Orbitals& pv = *params.pv;
    const int nat = ucell.nat;
    const int nspin = params.nspin;

    // 1+2. build the exx-form per-direction/atom/spin dH (dHexxs) from the current mixed DM
    ex->cal_exx_dHs(ucell, pv, nspin);

    // 3+4. convert dHexxs to per-atom-I HContainers and write, one spin channel at a time
    for (int ispin = 0; ispin < (nspin == 2 ? 2 : 1); ++ispin)
    {
        PerIContainers c(pv, nat);

        // OperatorEXX dereferences hR_in in its constructor and reallocates it, so pass a
        // throwaway container (its cell_nearest is built from kv and reused for dhR below).
        hamilt::HContainer<double> hR_dummy(const_cast<Parallel_Orbitals*>(&pv));
        hamilt::OperatorEXX<hamilt::OperatorLCAO<double, double>> op_exx(nullptr, &hR_dummy, ucell, *params.kv);

        op_exx.cal_dH(ispin, c.g, ex->get_dHexxs());

        ModuleIO::write_dh_perI(params, ispin, "dvexxr", "dvexxk", "dV^EXX", c.g);
    }
}
#endif

} // namespace

bool write_dH_t(WriteDHParams& params)
{
    ModuleBase::TITLE("ModuleIO", "write_dH_t");
    ModuleBase::timer::start("ModuleIO", "write_dH_t");

    const UnitCell& ucell = *params.ucell;
    const Grid_Driver& gd = *params.gd;
    const Parallel_Orbitals& pv = *params.pv;
    const TwoCenterBundle& two_center_bundle = *params.two_center_bundle;
    const LCAO_Orbitals& orb = *params.orb;

    const int nat = ucell.nat;
    const int nspin = params.nspin;
    const std::vector<double>& orb_cutoff = orb.cutoffs();

    for (int ispin = 0; ispin < (nspin == 2 ? 2 : 1); ispin++)
    {
        // per-atom-I containers: dT_*[iat] = d<phi|T|phi>/dtau_iat
        PerIContainers c(pv, nat);

        hamilt::EKinetic<hamilt::OperatorLCAO<double, double>> tmp_ekinetic(nullptr,
                                                                            params.kv->kvec_d,
                                                                            nullptr,
                                                                            &ucell,
                                                                            orb_cutoff,
                                                                            &gd,
                                                                            two_center_bundle.kinetic_orb.get());

        tmp_ekinetic.cal_dH(c.g);

        ModuleIO::write_dh_perI(params, ispin, "dtr", "dtk", "dT", c.g);
    }

    ModuleBase::timer::end("ModuleIO", "write_dH_t");
    return true;
}

bool write_dH_vnl(WriteDHParams& params)
{
    ModuleBase::TITLE("ModuleIO", "write_dH_vnl");
    ModuleBase::timer::start("ModuleIO", "write_dH_vnl");

    const UnitCell& ucell = *params.ucell;
    const Grid_Driver& gd = *params.gd;
    const Parallel_Orbitals& pv = *params.pv;
    const TwoCenterBundle& two_center_bundle = *params.two_center_bundle;
    const int nat = ucell.nat;
    const int nspin = params.nspin;
    const std::vector<double>& orb_cutoff = params.orb->cutoffs();

    for (int ispin = 0; ispin < (nspin == 2 ? 2 : 1); ispin++)
    {
        PerIContainers c(pv, nat);

        hamilt::Nonlocal<hamilt::OperatorLCAO<double, double>> tmp_nonlocal(nullptr,
                                                                            params.kv->kvec_d,
                                                                            nullptr,
                                                                            &ucell,
                                                                            orb_cutoff,
                                                                            &gd,
                                                                            two_center_bundle.overlap_orb_beta.get());

        tmp_nonlocal.cal_dH(c.g);

        ModuleIO::write_dh_perI(params, ispin, "dvnlr", "dvnlk", "dV^NL", c.g);
    }

    ModuleBase::timer::end("ModuleIO", "write_dH_vnl");
    return true;
}

bool write_dH_vl(WriteDHParams& params)
{
    ModuleBase::TITLE("ModuleIO", "write_dH_vl");
    ModuleBase::timer::start("ModuleIO", "write_dH_vl");

    const bool ok = write_dH_veff_term(params, params.pot_vl, "vl", "dvlr", "dvlk", "dV^L");

    ModuleBase::timer::end("ModuleIO", "write_dH_vl");
    return ok;
}

bool write_dH_vh(WriteDHParams& params)
{
    ModuleBase::TITLE("ModuleIO", "write_dH_vh");
    ModuleBase::timer::start("ModuleIO", "write_dH_vh");

    const bool ok = write_dH_veff_term(params, params.pot_vh, "hartree", "dvhr", "dvhk", "dV^H");

    ModuleBase::timer::end("ModuleIO", "write_dH_vh");
    return ok;
}

bool write_dH_vh_pulay(WriteDHParams& params)
{
    ModuleBase::TITLE("ModuleIO", "write_dH_vh_pulay");
    ModuleBase::timer::start("ModuleIO", "write_dH_vh_pulay");

    const bool ok = write_dH_veff_term(params, params.pot_vh, "none", "dvhr_pulay_", "dvhk_pulay_", "dV^H (Pulay)");

    ModuleBase::timer::end("ModuleIO", "write_dH_vh_pulay");
    return ok;
}

bool write_dH_vxc(WriteDHParams& params)
{
    ModuleBase::TITLE("ModuleIO", "write_dH_vxc");
    ModuleBase::timer::start("ModuleIO", "write_dH_vxc");

    const bool ok = write_dH_veff_term(params, params.pot_vxc, "none", "dvxcr", "dvxck", "dV^XC");

    ModuleBase::timer::end("ModuleIO", "write_dH_vxc");
    return ok;
}

#ifdef __EXX
bool write_dH_exx(WriteDHParams& params)
{
    ModuleBase::TITLE("ModuleIO", "write_dH_exx");
    ModuleBase::timer::start("ModuleIO", "write_dH_exx");

    bool ok = false;
    // exd (real Hexx) and exc (complex Hexx) are mutually exclusive; pick by real_number.
    if (GlobalC::exx_info.info_ri.real_number)
    {
        if (params.exd != nullptr)
        {
            write_dH_exx_impl(params, params.exd);
            ok = true;
        }
    }
    else
    {
        if (params.exc != nullptr)
        {
            write_dH_exx_impl(params, params.exc);
            ok = true;
        }
    }

    ModuleBase::timer::end("ModuleIO", "write_dH_exx");
    return ok;
}
#endif

} // namespace ModuleIO
