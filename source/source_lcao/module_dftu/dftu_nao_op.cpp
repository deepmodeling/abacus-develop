#include "dftu_nao_op.h"

#include "source_base/timer.h"
#include "source_base/tool_title.h"
#include "source_basis/module_nao/two_center_integrator.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "source_cell/unitcell.h"
#include "source_estate/module_dm/density_matrix.h"
#include "source_lcao/module_operator_lcao/operator_lcao.h"
#include "source_pw/module_pwdft/dftu_base.h"
#include "source_base/parallel_reduce.h"
#include "source_cell/klist.h"
#include "source_cell/module_symmetry/symmetry.h"
#include "source_lcao/module_ri/ri_util.h"

#include <memory>

#include "dftu_nao_adj.h"
#include "dftu_nao_fs_r.h"
#include "dftu_nao_ijr.h"
#include "dftu_nao_pots.h"

template <typename TK, typename TR>
hamilt::DFTU<hamilt::OperatorLCAO<TK, TR>>::DFTU(HS_Matrix_K<TK>* hsk_in,
                                                 const std::vector<ModuleBase::Vector3<double>>& kvec_d_in,
                                                 hamilt::HContainer<TR>* hR_in,
                                                 const UnitCell& ucell_in,
                                                 const Grid_Driver* GridD_in,
                                                 const TwoCenterIntegrator* intor,
                                                 const std::vector<double>& orb_cutoff,
                                                 Plus_U_Base* p_dftu,
                                                 const int nspin_in,
                                                 const double onsite_radius,
                                                 const elecstate::DensityMatrix<TK, double>* dm_in)
    : hamilt::OperatorLCAO<TK, TR>(hsk_in, kvec_d_in, hR_in),
      ucell(&ucell_in),
      dftu(p_dftu),
      dm_(dm_in),
      intor_(intor),
      orb_cutoff_(orb_cutoff),
      nspin(nspin_in)
{
    ModuleBase::timer::start("DFTU", "DFTU");
    this->cal_type = calculation_type::lcao_dftu;

    assert(this->ucell != nullptr);
    assert(this->dm_ != nullptr);

    // structure snapshot: both members depend only on the atomic structure
    // and are computed once here; the operator is rebuilt every ionic step.
    // Kept in the constructor body (not the initializer list) so failures
    // inside these calls are easy to debug.
    this->adjs_all = DFTU_LCAO::build_adjacent_atoms(this->ucell, this->dftu, GridD_in, this->orb_cutoff_, onsite_radius);
    const Parallel_Orbitals* pv = this->hR->get_atom_pair(0).get_paraV();
    this->nlm_tot = DFTU_LCAO::cal_nlm_all(*this->ucell, *this->dftu, *this->intor_, this->adjs_all, *pv);

    ModuleBase::timer::end("DFTU", "DFTU");
}

// contributeHR()
/**
 * @brief Contribute DFT+U Hamiltonian to real-space HR matrix
 * 
 * @details This function handles different scenarios based on:
 * 1. Whether occ_mat (occupation matrix) is read from file (is_occmat_ready)
 * 2. Spin configuration (nspin=1, 2, or 4)
 * 3. SCF iteration stage (first vs subsequent iterations)
 * 
 * Case 1: Occ_mat NOT ready (!is_occmat_ready)
 *   - First electronic iteration: calculates occupation matrix from density matrix (DMR)
 *     * Fetches the real-space DMR via dm_->get_DMR_pointer()
 *     * Accumulates contributions from all atom pairs via DFTU_LCAO::cal_occ_ijr()
 *     * Performs MPI reduction to sum occ across processes
 *     * Stores result via set_occ_mat_flat() for use in pot_onsite calculation
 *     * For nspin=1: occ is scaled by 0.5 (since only one spin channel computed)
 *   - Subsequent iterations: occ_mat is computed fresh each iteration from updated DMR
 * 
 * Case 2: Occ_mat IS ready (is_occmat_ready, i.e., read from dm_onsite.txt file)
 *   - First electronic iteration: uses pre-read occ_mat directly without DMR calculation
 *     * Skips DMR-based occ calculation entirely
 *     * Reads occ_mat from stored data via get_occ_mat()
 *     * Different indexing for nspin=4 vs nspin=1/2 (see below)
 *   - After first iteration: set_occmat_stale() is called to force recomputation
 * 
 * Spin configurations:
 *   nspin=1 (non-spin-polarized):
 *     - Single spin channel, occ computed once
 *     - Energy correction doubled at end (set_double_energy)
 *     - current_spin always 0
 *   
 *   nspin=2 (collinear spin-polarized):
 *     - Two separate spin channels (spin-up: 0, spin-down: 1)
 *     - current_spin toggles between 0 and 1 across iterations
 *     - set_occmat_stale() called when current_spin == 1 (last spin)
 *     - HR accumulated separately for each spin
 *   
 *   nspin=4 (non-collinear/SOC):
 *     - Single 4x4 Pauli matrix representation per atom
 *     - occ has 4*(2l+1)^2 elements (spin_fold=4)
 *     - get_occ_mat uses spin=0, ipol indices for Pauli blocks
 *     - set_occmat_stale() always called (current_spin check always true)
 *     - No current_spin toggling (all spins handled simultaneously)
 * 
 * @warning THREAD SAFETY: DFTU_LCAO::cal_hr_ijr() updates shared HR matrix entries.
 *          Different iat0 may contribute to same HR(iat1, iat2, R), requiring
 *          critical section protection for multithreaded correctness.
 *          TODO: Consider refactoring to atom_row_list pattern (see nonlocal.cpp)
 *          for better parallel performance instead of critical section.
 */
template <typename TK, typename TR>
void hamilt::DFTU<hamilt::OperatorLCAO<TK, TR>>::contributeHR()
{
    ModuleBase::TITLE("DFTU", "contributeHR");
    // Early exit: DMR not available (first SCF iteration before the first
    // diagonalization) AND occ_mat not yet initialized
    const bool dmr_null = (this->dm_ == nullptr || !this->dm_->is_dmr_ready());
    const bool occ_mat_not_init = !this->dftu->is_occmat_ready();

    if (dmr_null && occ_mat_not_init)
    {
        return;
    }
    if (this->current_spin == 0)
    {
        this->dftu->set_energy(0.0);
    }
    ModuleBase::timer::start("DFTU", "contributeHR");

    const Parallel_Orbitals* pv = this->hR->get_atom_pair(0).get_paraV();
    // nlm_tot is precomputed in the constructor (structure snapshot)

    // (symmetry) when crystal symmetry reduces the k-mesh, dm_->get_DMR_pointer()
    // was Fourier-transformed from the irreducible k-points only and is not
    // actually symmetric; reconstruct the full-BZ DMR once here (reused by every
    // atom below) via the same D(k) restoration EXX already uses for its own
    // real-space density matrix (ModuleSymmetry::Symmetry_rotation::restore_dm).
    std::unique_ptr<elecstate::DensityMatrix<TK, double>> dmr_sym;
    if (!this->dftu->is_occmat_ready() && this->kv_ != nullptr && ModuleSymmetry::Symmetry::symm_flag == 1
        && !this->kv_->kstars.empty())
    {
        if (!this->symrot_built_)
        {
            const std::array<int, 3>& period = RI_Util::get_Born_vonKarmen_period(*this->kv_);
            // for return_lattice to calculate Ms
            this->symrot_.find_irreducible_sector(this->ucell->symm, this->ucell->atoms, this->ucell->st,
                RI_Util::get_Born_von_Karmen_cells(period), period, this->ucell->lat);
            this->symrot_.cal_Ms(*this->kv_, *this->ucell, *pv);
            this->symrot_built_ = true;
        }
        const int nspin0 = (this->nspin == 2) ? 2 : 1;
        std::vector<ModuleBase::Vector3<double>> kvec_d_full;
        for (int ik_ibz = 0; ik_ibz < static_cast<int>(this->kv_->kstars.size()); ++ik_ibz)
        {
            for (auto& isym_kvd : this->kv_->kstars[ik_ibz]) { kvec_d_full.push_back(isym_kvd.second); }
        }
        const std::vector<std::vector<TK>> dmk_full = this->symrot_.restore_dm(*this->kv_, this->dm_->get_DMK_vector(), *pv);
        dmr_sym.reset(new elecstate::DensityMatrix<TK, double>(pv, nspin0, kvec_d_full, static_cast<int>(kvec_d_full.size())));
        dmr_sym->init_DMR(*this->dm_->get_DMR_pointer(1));
        dmr_sym->get_DMK_vector() = dmk_full;
        dmr_sym->cal_DMR();
    }

    // loop over all Hubbard-projector center atoms (iat0)
    int atom_index = 0;
    for (int iat0 = 0; iat0 < this->ucell->nat; iat0++)
    {
        int T0 = 0;
        int I0 = 0;
        ucell->iat2iait(iat0, &I0, &T0);
        if (!this->dftu->has_l_channel(T0))
        {
            continue;
        }
        const int target_L = this->dftu->get_l_channel(T0);
        const int tlp1 = 2 * target_L + 1;
        AdjacentAtomInfo& adjs = this->adjs_all[atom_index++];

        const int spin_fold = (this->nspin == 4) ? 4 : 1;
        std::vector<double> occ(tlp1 * tlp1 * spin_fold, 0.0);

        // compute or load occupation matrix
        if (!this->dftu->is_occmat_ready())
        {
            // DMR is guaranteed ready here: otherwise the early exit above
            // would have returned. DMR index is 1-based, hence +1.
            const hamilt::HContainer<double>* dmr = this->dm_->get_DMR_pointer(this->current_spin + 1);
            if (dmr_sym) { dmr = dmr_sym->get_DMR_pointer(this->current_spin + 1); }
            DFTU_LCAO::compute_occ_from_dmr(*this->ucell,
                                            *this->dftu,
                                            iat0,
                                            target_L,
                                            this->current_spin,
                                            this->nspin,
                                            adjs,
                                            *pv,
                                            this->nlm_tot,
                                            *dmr,
                                            occ);
        }
        else
        {
            DFTU_LCAO::load_occ_from_file(*this->dftu, iat0, target_L, this->nspin, this->current_spin, occ);
        }

        // compute Hubbard potential and energy
        const double u_value = this->dftu->get_u_current(T0);
        std::vector<double> pot_onsite_tmp(occ.size());
        double u_energy = this->dftu->get_energy();
        DFTU_LCAO::cal_pot_onsite(occ, tlp1, u_value, pot_onsite_tmp.data(), u_energy);
        this->dftu->set_energy(u_energy);

        std::vector<TR> pot_onsite(occ.size());
        DFTU_LCAO::transfer_pot_onsite(pot_onsite_tmp, pot_onsite);

        // accumulate HR contributions from neighbor pairs
        DFTU_LCAO::accumulate_hr_for_iat0<TR>(*this->ucell, this->hR, this->nlm_tot, iat0, adjs, *pv, pot_onsite);
    }

    // post-processing: energy doubling for nspin=1
    if (this->nspin == 1)
    {
        this->dftu->set_double_energy();
    }
    // mark occ_mat stale for next iteration
    if (this->current_spin == this->nspin - 1 || this->nspin == 4)
    {
        this->dftu->set_occmat_stale();
    }
    // toggle spin channel for nspin=2
    if (this->nspin == 2)
    {
        this->current_spin = 1 - this->current_spin;
    }

    ModuleBase::timer::end("DFTU", "contributeHR");
}

template class hamilt::DFTU<hamilt::OperatorLCAO<double, double>>;
template class hamilt::DFTU<hamilt::OperatorLCAO<std::complex<double>, double>>;
template class hamilt::DFTU<hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>>;
