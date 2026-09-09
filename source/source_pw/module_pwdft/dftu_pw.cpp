#include "source_pw/module_pwdft/dftu_pw.h"

#include "source_pw/module_pwdft/dftu_pw_tools.h"
#include "source_pw/module_pwdft/dftu_base.h"
#include "source_pw/module_pwdft/dftu_base_io.h"
#include "source_pw/module_pwdft/onsite_proj.h"
#include "source_cell/unitcell.h"
#include "source_estate/module_charge/charge_mixing.h"
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
    auto* onsite_p = projectors::OnsiteProjector<double, Device>::get_instance();
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
                    occmat.mat(iat, target_l, 0, 0).c,
                    becp, nbands, npol, nkb, begin_ih, m_begin, tlp1,
                    wg_in, ik);
            }
            else // nspin=1 or nspin=2
            {
                pw::accumulate_occ_scalar(
                    occmat.mat(iat, target_l, 0, is).c,
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
