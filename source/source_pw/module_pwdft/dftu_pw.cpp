#include "source_pw/module_pwdft/dftu_pw.h"

#include "source_pw/module_pwdft/dftu_pw_tools.h"
#include "source_pw/module_pwdft/dftu_base.h"
#include "source_pw/module_pwdft/dftu_base_io.h"
#include "source_pw/module_proj/onsite_proj.h"
#include "source_cell/unitcell.h"
#include "source_estate/module_charge/chg_mix.h"
#include "source_base/timer.h"

namespace DFTU_BASE {

template <typename Device>
void accumulate_occ_one_k(const void* psi_in,
                          const ModuleBase::matrix& wg_in,
                          const UnitCell& cell,
                          const int* isk,
                          const int nspin,
                          const std::vector<int>& l_channel,
                          OccupationMatrix& occmat)
{
    projectors::OnsiteProjector<double, Device>* onsite_p = projectors::OnsiteProjector<double, Device>::get_instance();
    const psi::Psi<std::complex<double>, Device>* psi_p =
        (const psi::Psi<std::complex<double>, Device>*)psi_in;
    const int nbands = psi_p->get_nbands();
    const int npol = psi_p->get_npol();
    for(int ik = 0; ik < psi_p->get_nk(); ik++)
    {
        int is = (nspin == 2) ? isk[ik] : 0;
        psi_p->fix_k(ik);
        onsite_p->tabulate_atomic(ik);

        onsite_p->overlap_proj_psi(nbands*npol, psi_p->get_pointer());
        const std::complex<double>* becp = onsite_p->get_h_becp();
        int nkb = onsite_p->get_size_becp() / nbands / npol;

        int begin_ih = 0;
        for(int iat = 0; iat < cell.nat; iat++)
        {
            const int it = cell.iat2it[iat];
            const int nh = onsite_p->get_nh(iat);
            const int target_l = l_channel[it];
            if(target_l == -1)
            {
                begin_ih += nh;
                continue;
            }
            const int m_begin = target_l * target_l;
            const int tlp1 = 2 * target_l + 1;
            if(nspin == 4)
            {
                pw::accumulate_occ_spinor(
                    occmat.mat(iat, target_l, 0).c,
                    becp, nbands, npol, nkb, begin_ih, m_begin, tlp1,
                    wg_in, ik);
            }
            else // nspin=1 or nspin=2
            {
                pw::accumulate_occ_scalar(
                    occmat.mat(iat, target_l, is).c,
                    becp, nbands, nkb, begin_ih, m_begin, tlp1,
                    wg_in, ik);
            }
            begin_ih += nh;
        }
    }
}

void cal_occ_pw(const void* psi_in,
                const ModuleBase::matrix& wg_in,
                const UnitCell& cell,
                Charge_Mixing* p_chgmix,
                const int* isk,
                const int kpar,
                const int nspin,
                const std::string& device,
                const std::vector<int>& l_channel,
                const std::vector<double>& u_current,
                const std::vector<int>& uterm_mat_index,
                OccupationMatrix& occmat,
                OccMatMixer* occ_mixer,
                std::vector<std::complex<double>>& uterm_mat,
                double& energy_u)
{
    ModuleBase::timer::start("Plus_U_Base", "cal_occ_pw");
    occmat.copy_to_save(cell, l_channel);
    if (occ_mixer != nullptr)
    {
        occ_mixer->begin_iter(occmat);
    }
    occmat.zero(cell, l_channel);

    if (device == "cpu")
    {
        DFTU_BASE::accumulate_occ_one_k<base_device::DEVICE_CPU>(
            psi_in, wg_in, cell, isk, nspin, l_channel, occmat);
    }
#if defined(__CUDA) || defined(__ROCM)
    else
    {
        DFTU_BASE::accumulate_occ_one_k<base_device::DEVICE_GPU>(
            psi_in, wg_in, cell, isk, nspin, l_channel, occmat);
    }
#endif

    // reduce occ_mat across k-pools
    pw::reduce_occ_mat(cell, nspin, kpar, l_channel, occmat);

    // mixing: flatten the fresh occ, mix against the saved one, write back
    if (occ_mixer != nullptr && p_chgmix != nullptr)
    {
        occ_mixer->collect(occmat);
        p_chgmix->mix_uom(occ_mixer->uom(), occ_mixer->uom_save());
        occ_mixer->write_back(occmat);
    }

    pw::compute_pot_uterm_and_energy(cell, nspin,
        u_current, l_channel, uterm_mat_index,
        occmat, uterm_mat, energy_u);

    ModuleBase::timer::end("Plus_U_Base", "cal_occ_pw");
}

template <typename Device>
std::vector<ScfDiagnostics::ScfOccupationBlock> observe_occ_pw(
    const void* psi_in, const ModuleBase::matrix& wg, const UnitCell& cell,
    const int* isk, const int kpar, const Plus_U_Base& dftu, const bool input_ready)
{
    OccupationMatrix raw;
    const std::vector<int>& channels = dftu.get_l_channel_vec();
    raw.init(cell, channels, 2, 1);
    raw.zero(cell, channels);
    const auto* wavefunctions = static_cast<const psi::Psi<std::complex<double>, Device>*>(psi_in);
    const int saved_k = wavefunctions->get_current_k();
    accumulate_occ_one_k<Device>(psi_in, wg, cell, isk, 2, channels, raw);
    pw::reduce_occ_mat(cell, 2, kpar, channels, raw);
    wavefunctions->fix_k(saved_k);
    auto* projector = projectors::OnsiteProjector<double, Device>::get_instance();
    projector->tabulate_atomic(saved_k);
    const int nbands = wavefunctions->get_nbands();
    projector->overlap_proj_psi(nbands, wavefunctions->get_pointer());
    std::vector<ScfDiagnostics::ScfOccupationBlock> blocks;
    for (int atom = 0; atom < cell.nat; ++atom)
    {
        const int type = cell.iat2it[atom];
        const int l = channels[type];
        if (l < 0) { continue; }
        const int size = (2*l+1)*(2*l+1);
        for (int spin = 0; spin < 2; ++spin)
        {
            ScfDiagnostics::ScfOccupationBlock block;
            block.atom = atom+1;
            block.spin = spin;
            block.l = l;
            block.u_eV = dftu.get_u_current(type)*ModuleBase::Ry_to_eV;
            const double* output = raw.mat(atom,l,spin).c;
            block.raw_output.assign(output,output+size);
            if (input_ready)
            {
                const double* input = dftu.occmat().mat(atom,l,spin).c;
                block.input.assign(input,input+size);
            }
            blocks.push_back(block);
        }
    }
    return blocks;
}

} // namespace DFTU_BASE

// explicit instantiations
template void DFTU_BASE::accumulate_occ_one_k<base_device::DEVICE_CPU>(
    const void*, const ModuleBase::matrix&, const UnitCell&, const int*,
    const int, const std::vector<int>&, OccupationMatrix&);
#if defined(__CUDA) || defined(__ROCM)
template void DFTU_BASE::accumulate_occ_one_k<base_device::DEVICE_GPU>(
    const void*, const ModuleBase::matrix&, const UnitCell&, const int*,
    const int, const std::vector<int>&, OccupationMatrix&);
#endif

template std::vector<ScfDiagnostics::ScfOccupationBlock> DFTU_BASE::observe_occ_pw<base_device::DEVICE_CPU>(
    const void*, const ModuleBase::matrix&, const UnitCell&, const int*, const int, const Plus_U_Base&, const bool);
#if defined(__CUDA) || defined(__ROCM)
template std::vector<ScfDiagnostics::ScfOccupationBlock> DFTU_BASE::observe_occ_pw<base_device::DEVICE_GPU>(
    const void*, const ModuleBase::matrix&, const UnitCell&, const int*, const int, const Plus_U_Base&, const bool);
#endif
