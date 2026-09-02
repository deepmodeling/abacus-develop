#include "source_pw/module_pwdft/dftu_base.h"
#include "source_pw/module_pwdft/dftu_tools_pw.h"
#include "source_pw/module_pwdft/onsite_proj.h"
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
    this->copy_occ_mat(cell);
    this->zero_occ_mat(cell);

    if(this->device == "cpu")
    {
        this->accumulate_occ_one_k<base_device::DEVICE_CPU>(psi_in, wg_in, cell, isk);
    }
#if defined(__CUDA) || defined(__ROCM)
    else
    {
        this->accumulate_occ_one_k<base_device::DEVICE_GPU>(psi_in, wg_in, cell, isk);
    }
#endif

    // reduce occ_mat across k-pools, then copy to uom_array for mixing
    this->reduce_occ_mat(cell);
    this->sync_occ_to_uom(cell);

    // mixing
    if(is_mixing_enabled() && p_chgmix != nullptr)
    {
        p_chgmix->mix_uom(this->uom_array, this->uom_save);
        this->set_occ_mat(cell);
    }

    this->compute_eff_pot_and_energy(cell);

    ModuleBase::timer::end("Plus_U_Base", "cal_occ_pw");
}

void Plus_U_Base::update_eff_pot_pw(const UnitCell& cell)
{
    ModuleBase::TITLE("Plus_U_Base", "update_eff_pot_pw");
    ModuleBase::timer::start("Plus_U_Base", "update_eff_pot_pw");

    if (this->is_occ_mat_initialized())
    {
        this->compute_eff_pot_and_energy(cell);
    }

    ModuleBase::timer::end("Plus_U_Base", "update_eff_pot_pw");
}

/// reduce occ_mat across all k-pools.
///
/// Each k-pool only accumulates occ_mat contributions from the k-points it
/// owns; this sums them across pools so occ_mat holds the full result.
/// nspin=1: single channel, size elements
/// nspin=2: two channels (spin-up/down) reduced separately
/// nspin=4: 4 Pauli blocks packed contiguously, reduced in one shot
void Plus_U_Base::reduce_occ_mat(const UnitCell& cell)
{
    for(int iat = 0; iat < cell.nat; iat++)
    {
        const int it = cell.iat2it[iat];
        const int target_l = get_orbital_corr(it);
        if(!has_correlated_orbital(it))
        {
            continue;
        }
        const int size = (2 * target_l + 1) * (2 * target_l + 1);

        if(this->nspin != 4)
        {
            Parallel_Reduce::reduce_double_allpool(this->kpar,
                    GlobalV::NPROC_IN_POOL,
                    this->occ_mat[iat][target_l][0][0].c,
                    size);
            if(this->nspin == 2)
            {
                Parallel_Reduce::reduce_double_allpool(this->kpar,
                        GlobalV::NPROC_IN_POOL,
                        this->occ_mat[iat][target_l][0][1].c,
                        size);
            }
        }
        else
        {
            Parallel_Reduce::reduce_double_allpool(this->kpar,
                    GlobalV::NPROC_IN_POOL,
                    this->occ_mat[iat][target_l][0][0].c,
                    size * 4);
        }
    }
}

/// copy occ_mat to uom_array for mixing.
///
/// Layout:
///   nspin=1: uom_array[pot_uterm_pw_index[iat] + mm] = occ_mat[...][0][0]
///   nspin=2: split layout [all_up | all_dn], each atom's spin-up in the
///            first half and spin-down in the second half, both indexed by
///            pot_uterm_pw_index[iat]
///   nspin=4: all 4 Pauli blocks are packed contiguously per atom
void Plus_U_Base::sync_occ_to_uom(const UnitCell& cell)
{
    if(this->uom_array.size() == 0)
    {
        return;
    }
    for(int iat = 0; iat < cell.nat; iat++)
    {
        const int it = cell.iat2it[iat];
        const int target_l = get_orbital_corr(it);
        if(!has_correlated_orbital(it))
        {
            continue;
        }
        const int size = (2 * target_l + 1) * (2 * target_l + 1);
        const int spin_fold = (this->nspin == 4) ? 4 : 1;

        for (int mm = 0; mm < size * spin_fold; mm++)
        {
            this->uom_array[pot_uterm_pw_index[iat] + mm] =
                this->occ_mat[iat][target_l][0][0].c[mm];
        }
        if(this->nspin == 2)
        {
            const int half_size = this->uom_array.size() / 2;
            for(int mm = 0; mm < size; mm++)
            {
                this->uom_array[half_size + pot_uterm_pw_index[iat] + mm] =
                    this->occ_mat[iat][target_l][0][1].c[mm];
            }
        }
    }
}

/// compute effective potential pot_onsite and DFT+U energy from occ_mat.
///
/// Preconditions:
///   - occ_mat has been accumulated from psi and reduced across k-pools
///     (cal_occ_pw calls this after the reduce + mixing steps).
///
/// Outputs:
///   - pot_uterm_pw: pot_onsite = U * (diag*delta - occ) written per atom
///     nspin=4: 4 Pauli blocks per atom, then transformed to spin basis
///     nspin=1: single channel
///     nspin=2: two channels in split layout [all_up | all_dn]
///   - energy_u: E_U = sum U * weight_eu * occ(m2,m1) * occ(m1,m2)
void Plus_U_Base::compute_eff_pot_and_energy(const UnitCell& cell)
{
    this->energy_u = 0.0;
    const double weight_eu = (this->nspin == 1) ? 1.0 : (this->nspin == 2) ? 0.5 : 0.25;
    const double diag_coeff = (this->nspin == 4) ? 1.0 : 0.5;
    // calculate pot_onsite and energy (occ_mat already reduced above)
    for(int iat = 0; iat < cell.nat; iat++)
    {
        const int it = cell.iat2it[iat];
        const int target_l = get_orbital_corr(it);
        if(!has_correlated_orbital(it))
        {
            continue;
        }
        const int size = (2 * target_l + 1) * (2 * target_l + 1);

        //update effective potential
        const double u_value = this->u_current[it];
        std::complex<double>* pot_onsite_iat = &(this->pot_uterm_pw[this->pot_uterm_pw_index[iat]]);
        const int m_size = 2 * target_l + 1;

        if(this->nspin == 4)
        {
            // pot_onsite is stored as 4 contiguous Pauli blocks per atom:
            //   is=0: charge channel (identity), Hubbard U contributes the
            //         diagonal term diag_coeff*delta(m1,m2)
            //   is=1,2,3: spin channels (sigma_x/y/z), no U diagonal term
            // The occupation matrix occ_mat[...][0][0].c packs all 4 blocks
            // contiguously, each of size m_size*m_size.
            this->energy_u += dftu_pw::compute_pot_onsite_spinor(
                pot_onsite_iat,
                this->occ_mat[iat][target_l][0][0].c,
                u_value, diag_coeff, weight_eu, m_size);
        }
        else // nspin=1 or nspin=2
        {
            // spin-up channel
            this->energy_u += dftu_pw::compute_pot_onsite_scalar(
                pot_onsite_iat,
                this->occ_mat[iat][target_l][0][0].c,
                u_value, diag_coeff, weight_eu, m_size);
            // spin-down channel for nspin=2
            if(this->nspin == 2)
            {
                std::complex<double>* pot_onsite_iat1 = &(this->pot_uterm_pw[this->pot_uterm_pw.size()/2 + this->pot_uterm_pw_index[iat]]);
                this->energy_u += dftu_pw::compute_pot_onsite_scalar(
                    pot_onsite_iat1,
                    this->occ_mat[iat][target_l][0][1].c,
                    u_value, diag_coeff, weight_eu, m_size);
            }
        }
    }
}

template <typename Device>
void Plus_U_Base::accumulate_occ_one_k(const void* psi_in,
                                       const ModuleBase::matrix& wg_in,
                                       const UnitCell& cell,
                                       const int* isk)
{
    auto* onsite_p = projectors::OnsiteProjector<double, Device>::get_instance();
    const psi::Psi<std::complex<double>, Device>* psi_p =
        (const psi::Psi<std::complex<double>, Device>*)psi_in;
    const int nbands = psi_p->get_nbands();
    const int npol = psi_p->get_npol();
    for(int ik = 0; ik < psi_p->get_nk(); ik++)
    {
        int is = (this->nspin == 2) ? isk[ik] : 0;
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
            const int target_l = get_orbital_corr(it);
            if(!has_correlated_orbital(it))
            {
                begin_ih += nh;
                continue;
            }
            const int m_begin = target_l * target_l;
            const int tlp1 = 2 * target_l + 1;
            if(this->nspin == 4)
            {
                dftu_pw::accumulate_occ_spinor(
                    this->occ_mat[iat][target_l][0][0].c,
                    becp, nbands, npol, nkb, begin_ih, m_begin, tlp1,
                    wg_in, ik);
            }
            else // nspin=1 or nspin=2
            {
                dftu_pw::accumulate_occ_scalar(
                    this->occ_mat[iat][target_l][0][is].c,
                    becp, nbands, nkb, begin_ih, m_begin, tlp1,
                    wg_in, ik);
            }
            begin_ih += nh;
        }
    }
}

// explicit instantiations
template void Plus_U_Base::accumulate_occ_one_k<base_device::DEVICE_CPU>(
    const void*, const ModuleBase::matrix&, const UnitCell&, const int*);
#if defined(__CUDA) || defined(__ROCM)
template void Plus_U_Base::accumulate_occ_one_k<base_device::DEVICE_GPU>(
    const void*, const ModuleBase::matrix&, const UnitCell&, const int*);
#endif
