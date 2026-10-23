#include "chg_tau.h"
#include "chg_uspp.h"

#include "source_base/timer.h"
#include "source_base/tool_quit.h"

namespace module_charge {
namespace detail {

void mix_tau_recip(Charge* chr,
                   const int nspin,
                   const bool double_grid,
                   ModulePW::PW_Basis* rhopw,
                   ModulePW::PW_Basis* rhodpw,
                   Base_Mixing::Mixing* mixing,
                   Base_Mixing::Mixing_Data& tau_mdata,
                   Base_Mixing::Plain_Mixing* mixing_highf)
{
    ModuleBase::TITLE("module_charge", "mix_tau_recip");
    ModuleBase::timer::start("module_charge", "mix_tau_recip");

    if (chr == nullptr)
    {
        ModuleBase::WARNING_QUIT("mix_tau_recip", "chr is null");
    }
    if (rhopw == nullptr || rhodpw == nullptr)
    {
        ModuleBase::WARNING_QUIT("mix_tau_recip", "grid pointer is null");
    }
    if (mixing == nullptr)
    {
        ModuleBase::WARNING_QUIT("mix_tau_recip", "mixing is null");
    }
    if (nspin < 1)
    {
        ModuleBase::WARNING_QUIT("mix_tau_recip", "nspin must be >= 1");
    }
    if (double_grid && mixing_highf == nullptr)
    {
        ModuleBase::WARNING_QUIT("mix_tau_recip", "mixing_highf is null when double_grid is on");
    }

    std::vector<std::complex<double>> kin_g(nspin * rhodpw->npw);
    std::vector<std::complex<double>> kin_g_save(nspin * rhodpw->npw);
    // FFT to get kin_g and kin_g_save
    for (int is = 0; is < nspin; ++is)
    {
        rhodpw->real2recip(chr->kin_r[is], &kin_g[is * rhodpw->npw]);
        rhodpw->real2recip(chr->kin_r_save[is], &kin_g_save[is * rhodpw->npw]);
    }

    // RAII owners for the smooth / high-frequency parts on the double grid;
    // raw pointers below alias these vectors when double_grid is on, or
    // alias kin_g[_save] directly when double_grid is off so the mixing
    // mutates the dense buffer in place.
    std::vector<std::complex<double>> tau_sg_in;
    std::vector<std::complex<double>> tau_sg_out;
    std::vector<std::complex<double>> tau_hf_in;
    std::vector<std::complex<double>> tau_hf_out;
    std::complex<double>* taugs_in = nullptr;
    std::complex<double>* taugs_out = nullptr;
    std::complex<double>* taughf_in = nullptr;
    std::complex<double>* taughf_out = nullptr;

    if (double_grid)
    {
        const int npw_smooth = rhopw->npw;
        const int npw_dense = rhodpw->npw;
        tau_sg_in.resize(nspin * npw_smooth);
        tau_hf_in.resize(nspin * (npw_dense - npw_smooth));
        tau_sg_out.resize(nspin * npw_smooth);
        tau_hf_out.resize(nspin * (npw_dense - npw_smooth));
        module_charge::split_dgrid(kin_g_save.data(), tau_sg_in, tau_hf_in,
                                   nspin, npw_smooth, npw_dense);
        module_charge::split_dgrid(kin_g.data(), tau_sg_out, tau_hf_out,
                                   nspin, npw_smooth, npw_dense);
        taugs_in = tau_sg_in.data();
        taughf_in = tau_hf_in.data();
        taugs_out = tau_sg_out.data();
        taughf_out = tau_hf_out.data();
    }
    else
    {
        taugs_in = kin_g_save.data();
        taugs_out = kin_g.data();
    }

    // Note: there is no kerker modification for tau because I'm not sure
    // if we should have it. If necessary we can try it in the future.
    mixing->push_data(tau_mdata, taugs_in, taugs_out, nullptr, false);
    mixing->mix_data(tau_mdata, taugs_out);

    if (double_grid)
    {
        // simple mixing for high_frequencies
        const int ndimhf = (rhodpw->npw - rhopw->npw) * nspin;
        mixing_highf->plain_mix(taughf_out, taughf_in, taughf_out, ndimhf, nullptr);

        // combine smooth part and high_frequency part
        module_charge::merge_dgrid(kin_g.data(), tau_sg_out, tau_hf_out,
                                   nspin, rhopw->npw, rhodpw->npw);
    }

    // kin_g to kin_r
    for (int is = 0; is < nspin; is++)
    {
        rhodpw->recip2real(&kin_g[is * rhodpw->npw], chr->kin_r[is]);
    }

    ModuleBase::timer::end("module_charge", "mix_tau_recip");
}

} // namespace detail
} // namespace module_charge
