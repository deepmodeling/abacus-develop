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

#include <memory>
#include <string>

namespace ModuleIO
{

namespace
{

// RAII holder for the per-atom-I dH containers: one (dx, dy, dz) HContainer per atom I.
// gx/gy/gz are raw-pointer views into the owned containers, ready to hand to an
// operator's cal_dH(...) and then to write_dh_perI(...).
struct PerIContainers
{
    std::vector<std::unique_ptr<hamilt::HContainer<double>>> ox, oy, oz;
    std::vector<hamilt::HContainer<double>*> gx, gy, gz;

    PerIContainers(const Parallel_Orbitals& pv, int nat)
    {
        ox.reserve(nat);
        oy.reserve(nat);
        oz.reserve(nat);
        gx.reserve(nat);
        gy.reserve(nat);
        gz.reserve(nat);
        for (int iat = 0; iat < nat; ++iat)
        {
            ox.push_back(std::make_unique<hamilt::HContainer<double>>(&pv));
            oy.push_back(std::make_unique<hamilt::HContainer<double>>(&pv));
            oz.push_back(std::make_unique<hamilt::HContainer<double>>(&pv));
            gx.push_back(ox.back().get());
            gy.push_back(oy.back().get());
            gz.push_back(oz.back().get());
        }
    }
};

// Shared driver for the Veff-based terms (V^L, V^H, V^XC), which differ only in the
// Hellmann-Feynman type passed to Veff::cal_dH and in the output prefixes/label.
bool write_dH_veff_term(WriteDHParams& params,
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
                                                                params.pot,
                                                                &hR_dummy,
                                                                &ucell,
                                                                orb_cutoff,
                                                                &gd,
                                                                nspin);

        PerIContainers c(pv, nat);

        veff.cal_dH(c.gx, c.gy, c.gz, hf_type);

        ModuleIO::write_dh_perI(params, ispin, rprefix, kprefix, label, c.gx, c.gy, c.gz);
    }
    return true;
}

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

        tmp_ekinetic.cal_dH(c.gx, c.gy, c.gz);

        ModuleIO::write_dh_perI(params, ispin, "dtr", "dtk", "dT", c.gx, c.gy, c.gz);
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

        tmp_nonlocal.cal_dH(c.gx, c.gy, c.gz);

        ModuleIO::write_dh_perI(params, ispin, "dvnlr", "dvnlk", "dV^NL", c.gx, c.gy, c.gz);
    }

    ModuleBase::timer::end("ModuleIO", "write_dH_vnl");
    return true;
}

bool write_dH_vl(WriteDHParams& params)
{
    ModuleBase::TITLE("ModuleIO", "write_dH_vl");
    ModuleBase::timer::start("ModuleIO", "write_dH_vl");

    const bool ok = write_dH_veff_term(params, "vl", "dvlr", "dvlk", "dV^L");

    ModuleBase::timer::end("ModuleIO", "write_dH_vl");
    return ok;
}

bool write_dH_vh(WriteDHParams& params)
{
    ModuleBase::TITLE("ModuleIO", "write_dH_vh");
    ModuleBase::timer::start("ModuleIO", "write_dH_vh");

    const bool ok = write_dH_veff_term(params, "hartree", "dvhr", "dvhk", "dV^H");

    ModuleBase::timer::end("ModuleIO", "write_dH_vh");
    return ok;
}

bool write_dH_vxc(WriteDHParams& params)
{
    ModuleBase::TITLE("ModuleIO", "write_dH_vxc");
    ModuleBase::timer::start("ModuleIO", "write_dH_vxc");

    const bool ok = write_dH_veff_term(params, "none", "dvxcr", "dvxck", "dV^XC");

    ModuleBase::timer::end("ModuleIO", "write_dH_vxc");
    return ok;
}

} // namespace ModuleIO
