#include "source_pw/module_pwdft/dftu_base.h"
#include "source_pw/module_pwdft/dftu_base_io.h"
#include "source_pw/module_pwdft/dftu_base_tools.h"
#include "source_pw/module_pwdft/onsite_proj.h"
#include "source_cell/unitcell.h"
#include "source_estate/module_charge/charge_mixing.h"
#include "source_base/parallel_reduce.h"
#include "source_base/global_variable.h"
#include "source_base/timer.h"
#include "source_base/parallel_global.h"




/// calculate occupation matrix for DFT+U (PW basis)
///
/// nspin=1 (npol=1): single spin channel; occ_mat[iat][l][n][0] only;
///   pot_uterm_pw has one block of tlp1^2 per atom.
///
/// nspin=2 (npol=1): two spin channels stored separately:
///   occ_mat[iat][l][n][0] = spin-up, occ_mat[iat][l][n][1] = spin-down;
///   becp indices: ib*nkb + begin_ih + m (same formula for both spins);
///   spin channel selected by `isk[ik]` (not ik >= nk/2, which fails for kpar>1);
///
/// nspin=4 (npol=2): spinor calculation;
///   occ_mat has a single matrix of size (2*tlp1) x (2*tlp1) per atom
///   storing all 4 Pauli blocks contiguously.
void Plus_U_Base::cal_occ_pw(const void* psi_in,
            const ModuleBase::matrix& wg_in,
            const UnitCell& cell,
            Charge_Mixing* p_chgmix,
            const int* isk)
{
    ModuleBase::timer::start("Plus_U_Base", "cal_occ_pw");
    this->occmat_.copy_to_save(cell, this->orbital_corr);
    this->occmat_.write_save_to_flat(cell, this->orbital_corr,
                                     this->pot_uterm_pw_index, this->uom_save);
    this->occmat_.zero(cell, this->orbital_corr);

    if(this->device == "cpu")
    {
        DFTU_BASE::accumulate_occ_one_k<base_device::DEVICE_CPU>(
            psi_in, wg_in, cell, isk, this->nspin, this->orbital_corr, this->occmat_);
    }
#if defined(__CUDA) || defined(__ROCM)
    else
    {
        DFTU_BASE::accumulate_occ_one_k<base_device::DEVICE_GPU>(
            psi_in, wg_in, cell, isk, this->nspin, this->orbital_corr, this->occmat_);
    }
#endif

    // reduce occ_mat across k-pools, then copy to uom_array for mixing
    DFTU_BASE::reduce_occ_mat(cell, this->nspin, this->kpar,
                              this->orbital_corr, this->occmat_);
    this->occmat_.write_to_flat(cell, this->orbital_corr,
                                this->pot_uterm_pw_index, this->uom_array);

    // mixing
    if(is_mixing_enabled() && p_chgmix != nullptr)
    {
        p_chgmix->mix_uom(this->uom_array, this->uom_save);
        this->occmat_.read_from_flat(cell, this->orbital_corr,
                                     this->pot_uterm_pw_index, this->uom_array);
    }

    DFTU_BASE::compute_pot_uterm_and_energy(cell, this->nspin,
        this->u_current, this->orbital_corr, this->pot_uterm_pw_index,
        this->occmat_, this->pot_uterm_pw, this->energy_u);

    ModuleBase::timer::end("Plus_U_Base", "cal_occ_pw");
}

namespace DFTU_BASE {

void reduce_occ_mat(const UnitCell& cell,
                    const int nspin,
                    const int kpar,
                    const std::vector<int>& orbital_corr,
                    OccupationMatrix& occmat)
{
    for(int iat = 0; iat < cell.nat; iat++)
    {
        const int it = cell.iat2it[iat];
        const int target_l = orbital_corr[it];
        if(target_l == -1)
        {
            continue;
        }
        const int size = (2 * target_l + 1) * (2 * target_l + 1);

        if(nspin != 4)
        {
            Parallel_Reduce::reduce_double_allpool(kpar,
                    GlobalV::NPROC_IN_POOL,
                    occmat.mat(iat, target_l, 0, 0).c,
                    size);
            if(nspin == 2)
            {
                Parallel_Reduce::reduce_double_allpool(kpar,
                        GlobalV::NPROC_IN_POOL,
                        occmat.mat(iat, target_l, 0, 1).c,
                        size);
            }
        }
        else
        {
            Parallel_Reduce::reduce_double_allpool(kpar,
                    GlobalV::NPROC_IN_POOL,
                    occmat.mat(iat, target_l, 0, 0).c,
                    size * 4);
        }
    }
}

void compute_pot_uterm_and_energy(const UnitCell& cell,
                                  const int nspin,
                                  const std::vector<double>& u_current,
                                  const std::vector<int>& orbital_corr,
                                  const std::vector<int>& pot_uterm_pw_index,
                                  const OccupationMatrix& occmat,
                                  std::vector<std::complex<double>>& pot_uterm_pw,
                                  double& energy_u)
{
    energy_u = 0.0;
    const double weight_eu = (nspin == 1) ? 1.0 : (nspin == 2) ? 0.5 : 0.25;
    const double diag_coeff = (nspin == 4) ? 1.0 : 0.5;
    // calculate pot_onsite and energy (occ_mat already reduced above)
    for(int iat = 0; iat < cell.nat; iat++)
    {
        const int it = cell.iat2it[iat];
        const int target_l = orbital_corr[it];
        if(target_l == -1)
        {
            continue;
        }
        const int size = (2 * target_l + 1) * (2 * target_l + 1);

        //update effective potential
        const double u_value = u_current[it];
        std::complex<double>* pot_onsite_iat = &(pot_uterm_pw[pot_uterm_pw_index[iat]]);
        const int m_size = 2 * target_l + 1;

        if(nspin == 4)
        {
            // pot_onsite is stored as 4 contiguous Pauli blocks per atom:
            //   is=0: charge channel (identity), Hubbard U contributes the
            //         diagonal term diag_coeff*delta(m1,m2)
            //   is=1,2,3: spin channels (sigma_x/y/z), no U diagonal term
            // The occupation matrix occ_mat[...][0][0].c packs all 4 blocks
            // contiguously, each of size m_size*m_size.
            energy_u += compute_pot_onsite_spinor(
                pot_onsite_iat,
                occmat.mat(iat, target_l, 0, 0).c,
                u_value, diag_coeff, weight_eu, m_size);
        }
        else // nspin=1 or nspin=2
        {
            // spin-up channel
            energy_u += compute_pot_onsite_scalar(
                pot_onsite_iat,
                occmat.mat(iat, target_l, 0, 0).c,
                u_value, diag_coeff, weight_eu, m_size);
            // spin-down channel for nspin=2
            if(nspin == 2)
            {
                std::complex<double>* pot_onsite_iat1 = &(pot_uterm_pw[pot_uterm_pw.size()/2 + pot_uterm_pw_index[iat]]);
                energy_u += compute_pot_onsite_scalar(
                    pot_onsite_iat1,
                    occmat.mat(iat, target_l, 0, 1).c,
                    u_value, diag_coeff, weight_eu, m_size);
            }
        }
    }
}

} // namespace DFTU_BASE

namespace DFTU_BASE {

template <typename Device>
void accumulate_occ_one_k(const void* psi_in,
                          const ModuleBase::matrix& wg_in,
                          const UnitCell& cell,
                          const int* isk,
                          const int nspin,
                          const std::vector<int>& orbital_corr,
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
            const int target_l = orbital_corr[it];
            if(target_l == -1)
            {
                begin_ih += nh;
                continue;
            }
            const int m_begin = target_l * target_l;
            const int tlp1 = 2 * target_l + 1;
            if(nspin == 4)
            {
                accumulate_occ_spinor(
                    occmat.mat(iat, target_l, 0, 0).c,
                    becp, nbands, npol, nkb, begin_ih, m_begin, tlp1,
                    wg_in, ik);
            }
            else // nspin=1 or nspin=2
            {
                accumulate_occ_scalar(
                    occmat.mat(iat, target_l, 0, is).c,
                    becp, nbands, nkb, begin_ih, m_begin, tlp1,
                    wg_in, ik);
            }
            begin_ih += nh;
        }
    }
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
