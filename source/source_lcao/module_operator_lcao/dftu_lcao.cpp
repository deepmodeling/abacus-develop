#include "dftu_lcao.h"

#include "source_base/timer.h"
#include "source_base/tool_title.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "source_lcao/module_operator_lcao/operator_lcao.h"
#include "source_lcao/module_hcontainer/hcontainer_funcs.h"
#include "source_io/module_parameter/parameter.h"
#ifdef _OPENMP
#include <unordered_set>
#endif
#include "source_base/parallel_reduce.h"

template <typename TK, typename TR>
hamilt::DFTU<hamilt::OperatorLCAO<TK, TR>>::DFTU(HS_Matrix_K<TK>* hsk_in,
                                                 const std::vector<ModuleBase::Vector3<double>>& kvec_d_in,
                                                 hamilt::HContainer<TR>* hR_in,
                                                 const UnitCell& ucell_in,
                                                 const Grid_Driver* GridD_in,
                                                 const TwoCenterIntegrator* intor,
                                                 const std::vector<double>& orb_cutoff,
                                                 Plus_U* p_dftu)
    : hamilt::OperatorLCAO<TK, TR>(hsk_in, kvec_d_in, hR_in), intor_(intor), orb_cutoff_(orb_cutoff)
{
    this->cal_type = calculation_type::lcao_dftu;
    this->ucell = &ucell_in;
    this->dftu = p_dftu;
#ifdef __DEBUG
    assert(this->ucell != nullptr);
#endif
    // initialize HR to allocate sparse Nonlocal matrix memory
    this->initialize_HR(GridD_in);
    // set nspin
    this->nspin = PARAM.inp.nspin;
}

// destructor
template <typename TK, typename TR>
hamilt::DFTU<hamilt::OperatorLCAO<TK, TR>>::~DFTU()
{
}

// initialize_HR()
template <typename TK, typename TR>
void hamilt::DFTU<hamilt::OperatorLCAO<TK, TR>>::initialize_HR(const Grid_Driver* GridD)
{
    ModuleBase::TITLE("DFTU", "initialize_HR");
    ModuleBase::timer::start("DFTU", "initialize_HR");

    this->adjs_all.clear();
    this->adjs_all.reserve(this->ucell->nat);
    for (int iat0 = 0; iat0 < ucell->nat; iat0++)
    {
        auto tau0 = ucell->get_tau(iat0);
        int T0=0;
        int I0=0;
        ucell->iat2iait(iat0, &I0, &T0);
        if (!this->dftu->has_correlated_orbital(T0))
        {
            continue;
        }
        const int target_L = this->dftu->get_orbital_corr(T0);

        AdjacentAtomInfo adjs;
        GridD->Find_atom(*ucell, tau0, T0, I0, &adjs);
        std::vector<bool> is_adj(adjs.adj_num + 1, false);
        for (int ad1 = 0; ad1 < adjs.adj_num + 1; ++ad1)
        {
            const int T1 = adjs.ntype[ad1];
            const int I1 = adjs.natom[ad1];
            const int iat1 = ucell->itia2iat(T1, I1);
            const ModuleBase::Vector3<double>& tau1 = adjs.adjacent_tau[ad1];
            const ModuleBase::Vector3<int>& R_index1 = adjs.box[ad1];
            // choose the real adjacent atoms
            // Note: the distance of atoms should less than the cutoff radius,
            // When equal, the theoretical value of matrix element is zero,
            // but the calculated value is not zero due to the numerical error, which would lead to result changes.
            if (this->ucell->cal_dtau(iat0, iat1, R_index1).norm() * this->ucell->lat0
                < orb_cutoff_[T1] + PARAM.inp.onsite_radius)
            {
                is_adj[ad1] = true;
            }
        }
        filter_adjs(is_adj, adjs);
        this->adjs_all.push_back(adjs);
    }

    ModuleBase::timer::end("DFTU", "initialize_HR");
}

template <typename TK, typename TR>
void hamilt::DFTU<hamilt::OperatorLCAO<TK, TR>>::cal_nlm_all(const Parallel_Orbitals* paraV)
{
    ModuleBase::TITLE("DFTU", "cal_nlm_all");
	if (this->precal_nlm_done) 
	{
		return;
	}

    ModuleBase::timer::start("DFTU", "cal_nlm_all");
    nlm_tot.resize(this->ucell->nat);
    const int npol = this->ucell->get_npol();
    int atom_index = 0;
    for (int iat0 = 0; iat0 < ucell->nat; iat0++)
    {
        auto tau0 = ucell->get_tau(iat0);
        int T0=0;
        int I0=0;
        ucell->iat2iait(iat0, &I0, &T0);
        if (!this->dftu->has_correlated_orbital(T0))
        {
            continue;
        }
        const int target_L = this->dftu->get_orbital_corr(T0);
        const int tlp1 = 2 * target_L + 1;
        AdjacentAtomInfo& adjs = this->adjs_all[atom_index++];

        // calculate and save the table of two-center integrals
        nlm_tot[iat0].resize(adjs.adj_num + 1);

        for (int ad = 0; ad < adjs.adj_num + 1; ++ad)
        {
            const int T1 = adjs.ntype[ad];
            const int I1 = adjs.natom[ad];
            const int iat1 = ucell->itia2iat(T1, I1);
            const ModuleBase::Vector3<double>& tau1 = adjs.adjacent_tau[ad];
            const Atom* atom1 = &ucell->atoms[T1];

            auto all_indexes = paraV->get_indexes_row(iat1);
            auto col_indexes = paraV->get_indexes_col(iat1);
            // insert col_indexes into all_indexes to get universal set with no repeat elements
            all_indexes.insert(all_indexes.end(), col_indexes.begin(), col_indexes.end());
            std::sort(all_indexes.begin(), all_indexes.end());
            all_indexes.erase(std::unique(all_indexes.begin(), all_indexes.end()), all_indexes.end());
            for (int iw1l = 0; iw1l < all_indexes.size(); iw1l += npol)
            {
                const int iw1 = all_indexes[iw1l] / npol;
                // only first zeta orbitals in target L of atom iat0 are needed
                std::vector<double> nlm_target(tlp1);
                const int L1 = atom1->iw2l[iw1];
                const int N1 = atom1->iw2n[iw1];
                const int m1 = atom1->iw2m[iw1];
                std::vector<std::vector<double>> nlm;
                // nlm is a vector of vectors, but size of outer vector is only 1 here
                // If we are calculating force, we need also to store the gradient
                // and size of outer vector is then 4
                // inner loop : all projectors (L0,M0)

                // convert m (0,1,...2l) to M (-l, -l+1, ..., l-1, l)
                const int M1 = (m1 % 2 == 0) ? -m1 / 2 : (m1 + 1) / 2;

                ModuleBase::Vector3<double> dtau = tau0 - tau1;
                intor_->snap(T1, L1, N1, M1, T0, dtau * this->ucell->lat0, false /*cal_deri*/, nlm);
                // select the elements of nlm with target_L
                for (int iw = 0; iw < this->ucell->atoms[T0].nw; iw++)
                {
                    const int L0 = this->ucell->atoms[T0].iw2l[iw];
                    if (L0 == target_L)
                    {
                        for (int m = 0; m < 2 * L0 + 1; m++)
                        {
                            nlm_target[m] = nlm[0][iw + m];
                        }
                        break;
                    }
                }
                nlm_tot[iat0][ad].insert({all_indexes[iw1l], nlm_target});
            }
        }
    }
    this->precal_nlm_done = true;
    ModuleBase::timer::end("DFTU", "cal_nlm_all");
}

// contributeHR()
/**
 * @brief Contribute DFT+U Hamiltonian to real-space HR matrix
 * 
 * @details This function handles different scenarios based on:
 * 1. Whether locale (occupation matrix) is read from file (is_locale_initialized)
 * 2. Spin configuration (nspin=1, 2, or 4)
 * 3. SCF iteration stage (first vs subsequent iterations)
 * 
 * Case 1: Locale NOT initialized (!is_locale_initialized)
 *   - First electronic iteration: calculates occupation matrix from density matrix (DMR)
 *     * Uses get_dmr(current_spin) to get real-space density matrix
 *     * Accumulates contributions from all atom pairs via cal_occ()
 *     * Performs MPI reduction to sum occ across processes
 *     * Stores result via set_locale_flat() for use in VU calculation
 *     * For nspin=1: occ is scaled by 0.5 (since only one spin channel computed)
 *   - Subsequent iterations: locale is computed fresh each iteration from updated DMR
 * 
 * Case 2: Locale IS initialized (is_locale_initialized, i.e., read from onsite.dm file)
 *   - First electronic iteration: uses pre-read locale directly without DMR calculation
 *     * Skips DMR-based occ calculation entirely
 *     * Reads locale from stored data via get_locale()
 *     * Different indexing for nspin=4 vs nspin=1/2 (see below)
 *   - After first iteration: mark_locale_dirty() is called to force recomputation
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
 *     - mark_locale_dirty() called when current_spin == 1 (last spin)
 *     - HR accumulated separately for each spin
 *   
 *   nspin=4 (non-collinear/SOC):
 *     - Single 4x4 Pauli matrix representation per atom
 *     - occ has 4*(2l+1)^2 elements (spin_fold=4)
 *     - get_locale uses spin=0, ipol indices for Pauli blocks
 *     - mark_locale_dirty() always called (current_spin check always true)
 *     - No current_spin toggling (all spins handled simultaneously)
 * 
 * @warning THREAD SAFETY: cal_HR_IJR() updates shared HR matrix entries.
 *          Different iat0 may contribute to same HR(iat1, iat2, R), requiring
 *          critical section protection for multithreaded correctness.
 *          TODO: Consider refactoring to atom_row_list pattern (see nonlocal.cpp)
 *          for better parallel performance instead of critical section.
 */
template <typename TK, typename TR>
void hamilt::DFTU<hamilt::OperatorLCAO<TK, TR>>::contributeHR()
{
    ModuleBase::TITLE("DFTU", "contributeHR");
    // Early exit conditions:
    // - get_dmr(0) == nullptr: DMR not available (typical in first iteration without file input)
    // - !is_locale_initialized(): locale not read from file AND not yet computed from DMR
    // When both true, skip DFT+U contribution entirely (first iteration, no file input)
    const bool dmr_null = (this->dftu->get_dmr(0) == nullptr);
    const bool locale_not_init = !this->dftu->is_locale_initialized();

    if (dmr_null && locale_not_init)
    {
        return;
    }
    else
    {
        // Reset DFT+U energy at start of each spin cycle
        // For nspin=2: reset when current_spin==0 (start of spin-up calculation)
        // For nspin=1/4: reset once (current_spin always 0)
		if (this->current_spin == 0) 
		{
            this->dftu->set_energy(0.0);
		}
	}
    ModuleBase::timer::start("DFTU", "contributeHR");

    const Parallel_Orbitals* paraV = this->hR->get_atom_pair(0).get_paraV();
    const int npol = this->ucell->get_npol();
    // 1. Calculate <psi|alpha> two-center integrals for all atom pairs
    //    This is reused in both occ and HR calculations
    this->cal_nlm_all(paraV);

    // 2. Loop over all Hubbard-projector center atoms (iat0)
    int atom_index = 0;
    for (int iat0 = 0; iat0 < this->ucell->nat; iat0++)
    {
        auto tau0 = ucell->get_tau(iat0);
        int T0, I0;
        ucell->iat2iait(iat0, &I0, &T0);
        if (!this->dftu->has_correlated_orbital(T0))
        {
            continue;
        }
        const int target_L = this->dftu->get_orbital_corr(T0);
        const int tlp1 = 2 * target_L + 1;
        AdjacentAtomInfo& adjs = this->adjs_all[atom_index++];

        ModuleBase::timer::start("DFTU", "cal_occ");
        // spin_fold: number of spin components in occ array
        // nspin=4: 4 (Pauli matrix blocks), nspin=1/2: 1 (single spin channel)
        const int spin_fold = (this->nspin == 4) ? 4 : 1;
        std::vector<double> occ(tlp1 * tlp1 * spin_fold, 0.0);
        
        // ============================================================
        // BRANCH 1: Locale NOT initialized (compute from DMR)
        // ============================================================
        // This branch is taken when:
        // - is_locale_initialized() == false (no file read or omc != 0)
        // - DMR is available (get_dmr() != nullptr)
        // Typical scenario: normal SCF iterations after first step
        if (!this->dftu->is_locale_initialized())
        {
            // TODO: UNSAFE - get_dmr(current_spin) assumes DMR has correct spin indexing.
            // For nspin=2, current_spin must be correctly toggled (0 then 1).
            // If current_spin is wrong, wrong spin channel's DMR is used.
            const hamilt::HContainer<double>* dmR_current = this->dftu->get_dmr(this->current_spin);
            for (int ad1 = 0; ad1 < adjs.adj_num + 1; ++ad1)
            {
                const int T1 = adjs.ntype[ad1];
                const int I1 = adjs.natom[ad1];
                const int iat1 = ucell->itia2iat(T1, I1);
                ModuleBase::Vector3<int>& R_index1 = adjs.box[ad1];
                const std::unordered_map<int, std::vector<double>>& nlm1 = nlm_tot[iat0][ad1];
                for (int ad2 = 0; ad2 < adjs.adj_num + 1; ++ad2)
                {
                    const int T2 = adjs.ntype[ad2];
                    const int I2 = adjs.natom[ad2];
                    const int iat2 = ucell->itia2iat(T2, I2);
                    const std::unordered_map<int, std::vector<double>>& nlm2 = nlm_tot[iat0][ad2];
                    ModuleBase::Vector3<int>& R_index2 = adjs.box[ad2];
                    ModuleBase::Vector3<int> R_vector(R_index2[0] - R_index1[0],
                                                      R_index2[1] - R_index1[1],
                                                      R_index2[2] - R_index1[2]);
                    const hamilt::BaseMatrix<double>* tmp
                        = dmR_current->find_matrix(iat1, iat2, R_vector[0], R_vector[1], R_vector[2]);
                    if (tmp != nullptr)
                    {
                        this->cal_occ(iat1, iat2, paraV, nlm1, nlm2, tmp->get_pointer(), occ);
                    }
                }
            }
#ifdef __MPI
            // CRITICAL: MPI reduction required for distributed DMR calculations.
            // Each process computes partial occ from its local DMR blocks.
            // Without this, occ would be incomplete and DFT+U potential wrong.
            // TODO: Verify that occ size is consistent across processes before reduction.
            // TODO: Consider using MPI_IN_PLACE to avoid extra buffer allocation.
            Parallel_Reduce::reduce_all(occ.data(), occ.size());
#endif
            // For nspin=1: occ computed from single spin channel, but should represent
            // total occupation (both spins). Scale by 0.5 to account for this.
            if (this->nspin == 1)
            {
                for (auto& v : occ) { v *= 0.5; }
            }
            this->dftu->set_locale_flat(iat0, target_L, this->current_spin, occ);
        }
        // ============================================================
        // BRANCH 2: Locale IS initialized (use pre-read data)
        // ============================================================
        // This branch is taken when:
        // - is_locale_initialized() == true (locale read from onsite.dm file)
        // - OR omc != 0 (occupation matrix control with initial_onsite.dm)
        // Typical scenario: first SCF iteration with file input, or restart calculation
        else
        {
            // nspin=4: Non-collinear case with Pauli matrix representation
            // Locale stored as single 4x4 block per atom, with spin indices embedded
            // in the matrix indices (ipol0, ipol1 for Pauli block indices)
            if (this->nspin == 4)
            {
                // For nspin=4, locale is stored as 4 stacked tlp1^2 blocks
                // at offsets 0, tlp1^2, 2*tlp1^2, 3*tlp1^2 for the 4 Pauli channels.
                // Use get_locale_flat to read the stacked blocks directly
                this->dftu->get_locale_flat(iat0, target_L, occ);
            }
            // nspin=1 or nspin=2: Collinear spin case
            // Locale stored separately for each spin channel
            else
            {
                for (int i = 0; i < static_cast<int>(occ.size()); i++)
                {
                    // TODO: UNSAFE - current_spin must be correct for nspin=2.
                    // If current_spin is not toggled properly, wrong spin channel's locale is read.
                    // This can happen if contributeHR() is called out of expected order.
                    occ[i] = this->dftu->get_locale(iat0, target_L, 0, this->current_spin,
                                                      i / (2 * target_L + 1), i % (2 * target_L + 1));
                }
            }
        }
        ModuleBase::timer::end("DFTU", "cal_occ");

        // 3. Calculate Hubbard potential VU from occupation matrix
        // VU = U * (1/2 * delta(m,m') - occ(m,m')) for each spin channel
        // Energy: EU = U * 1/2 * occ(m,m') * occ(m',m)
        ModuleBase::timer::start("DFTU", "cal_vu");
        const double u_value = this->dftu->U[T0];
        std::vector<double> VU_tmp(occ.size());

        // TODO: GLOBAL STATE - Plus_U::get_energy()/set_energy() uses static member variable.
        // This is NOT thread-safe for parallel SCF calculations.
        // TODO: Refactor to use instance member or pass energy by reference.
        double u_energy = Plus_U::get_energy();
        this->cal_v_of_u(occ, tlp1, u_value, VU_tmp.data(), u_energy);
        Plus_U::set_energy(u_energy);

        // 4. Convert VU to appropriate data type (real or complex)
        // For nspin=4 with complex Hamiltonian, VU needs Pauli matrix transformation
        std::vector<TR> VU(occ.size());
        this->transfer_vu(VU_tmp, VU);

        // 5. Second iteration: Calculate Hamiltonian matrix contribution
        // HR += <psi_I|beta_m> * VU(m,m') * <beta_m'|psi_{J,R}>
        // for all atom pairs <I,J,R> within cutoff
        // Note: different iat0 may contribute to the same HR(iat1, iat2, R), so we need to protect the update
        // to avoid race conditions in multithreading. Reference: nonlocal.cpp for the atom_row_list pattern.
        // TODO: CRITICAL SECTION PERFORMANCE - This critical section serializes HR updates.
        // For systems with many Hubbard atoms, this becomes a bottleneck.
        // Consider refactoring to atom_row_list pattern (see nonlocal.cpp lines 127-220):
        //   1. Use #pragma omp for to distribute iat0 across threads
        //   2. Each thread records its assigned iat0 in thread-local atom_row_list
        //   3. When updating HR(iat1, iat2, R), skip if iat1 not in thread's atom_row_list
        //   4. This eliminates race conditions without critical section
        for (int ad1 = 0; ad1 < adjs.adj_num + 1; ++ad1)
        {
            const int T1 = adjs.ntype[ad1];
            const int I1 = adjs.natom[ad1];
            const int iat1 = ucell->itia2iat(T1, I1);
            ModuleBase::Vector3<int>& R_index1 = adjs.box[ad1];
            const std::unordered_map<int, std::vector<double>>& nlm1 = nlm_tot[iat0][ad1];
            for (int ad2 = 0; ad2 < adjs.adj_num + 1; ++ad2)
            {
                const int T2 = adjs.ntype[ad2];
                const int I2 = adjs.natom[ad2];
                const int iat2 = ucell->itia2iat(T2, I2);
                const std::unordered_map<int, std::vector<double>>& nlm2 = nlm_tot[iat0][ad2];
                ModuleBase::Vector3<int>& R_index2 = adjs.box[ad2];
                ModuleBase::Vector3<int> R_vector(R_index2[0] - R_index1[0],
                                                  R_index2[1] - R_index1[1],
                                                  R_index2[2] - R_index1[2]);
                hamilt::BaseMatrix<TR>* tmp = this->hR->find_matrix(iat1, iat2, R_vector[0], R_vector[1], R_vector[2]);
                if (tmp != nullptr)
                {
#ifdef _OPENMP
#pragma omp critical(dftu_hr_update)
#endif
                    {
                        this->cal_HR_IJR(iat1, iat2, paraV, nlm1, nlm2, VU, tmp->get_pointer());
                    }
                }
            }
        }
        ModuleBase::timer::end("DFTU", "cal_vu");
    }

    // 6. Post-processing: Energy correction and locale state management
    // For nspin=1: DFT+U energy computed for single spin channel, but should count both spins
    // set_double_energy() doubles the energy to account for degenerate spin-up/down
	if (this->nspin == 1) 
	{
        this->dftu->set_double_energy();
	}
	
    // 7. Mark locale as dirty to force recomputation in next iteration
    // This is called when:
    // - nspin=4: Always (all spins handled simultaneously, current_spin==0==nspin-1)
    // - nspin=2: When current_spin==1 (after spin-down calculation, last spin channel)
    // - nspin=1: When current_spin==0==nspin-1 (always called)
    // 
    // Purpose: Ensure locale is recomputed from updated DMR in next SCF iteration,
    // rather than using stale pre-read data from file.
    // TODO: This logic is confusing. Consider explicit variable like `is_last_spin_channel`.
	if (this->current_spin == this->nspin - 1 || this->nspin == 4) 
	{
		this->dftu->mark_locale_dirty();
	}

    // 8. Spin channel toggling for nspin=2
    // nspin=2 requires separate HR updates for spin-up (current_spin=0) and spin-down (current_spin=1)
    // The HR matrix is updated twice per SCF iteration, once for each spin channel
    // current_spin toggles: 0 -> 1 -> 0 -> 1 ...
    // For nspin=1: current_spin always 0 (no toggling needed)
    // For nspin=4: current_spin always 0 (all spins handled simultaneously via Pauli matrices)
    // TODO: UNSAFE - This assumes contributeHR() is called in strict alternating order.
    // If called out of order (e.g., due to parallel k-point distribution), current_spin may be wrong.
    // TODO: Consider deriving current_spin from ik or explicit parameter instead of toggling.
	if (this->nspin == 2) 
	{
		this->current_spin = 1 - this->current_spin;
	}

    ModuleBase::timer::end("DFTU", "contributeHR");
}

// cal_HR_IJR()
template <typename TK, typename TR>
void hamilt::DFTU<hamilt::OperatorLCAO<TK, TR>>::cal_HR_IJR(
    const int& iat1,
    const int& iat2,
    const Parallel_Orbitals* paraV,
    const std::unordered_map<int, std::vector<double>>& nlm1_all,
    const std::unordered_map<int, std::vector<double>>& nlm2_all,
    const std::vector<TR>& VU,
    TR* data_pointer)
{

    // npol is the number of polarizations,
    // 1 for non-magnetic (one Hamiltonian matrix only has spin-up or spin-down),
    // 2 for magnetic (one Hamiltonian matrix has both spin-up and spin-down)
    const int npol = this->ucell->get_npol();
    // ---------------------------------------------
    // calculate the Nonlocal matrix for each pair of orbitals
    // ---------------------------------------------
    auto row_indexes = paraV->get_indexes_row(iat1);
    auto col_indexes = paraV->get_indexes_col(iat2);
    const int m_size = int(sqrt(VU.size()) / npol);
    // step_trace = 0 for NSPIN=1,2; ={0, 1, local_col, local_col+1} for NSPIN=4
    std::vector<int> step_trace(npol * npol, 0);
    for (int is = 0; is < npol; is++)
    {
        for (int is2 = 0; is2 < npol; is2++)
        {
            step_trace[is * npol + is2] = paraV->get_ncol_atom(iat2) * is + is2;
        }
    }
    // calculate the local matrix
    const TR* tmp_d = nullptr;
    for (int iw1l = 0; iw1l < row_indexes.size(); iw1l += npol)
    {
        const std::vector<double>& nlm1 = nlm1_all.find(row_indexes[iw1l])->second;
        for (int iw2l = 0; iw2l < col_indexes.size(); iw2l += npol)
        {
            const std::vector<double>& nlm2 = nlm2_all.find(col_indexes[iw2l])->second;
#ifdef __DEBUG
            assert(nlm1.size() == nlm2.size());
#endif
            for (int is = 0; is < npol * npol; ++is)
            {
                int start = is * m_size * m_size;
                TR nlm_tmp = TR(0);
                for (int m1 = 0; m1 < m_size; m1++)
                {
                    for (int m2 = 0; m2 < m_size; m2++)
                    {
                        nlm_tmp += nlm1[m1] * nlm2[m2] * VU[m1 * m_size + m2 + start];
                    }
                }
                data_pointer[step_trace[is]] += nlm_tmp;
            }
            data_pointer += npol;
        }
        data_pointer += (npol - 1) * col_indexes.size();
    }
}

template <typename TK, typename TR>
void hamilt::DFTU<hamilt::OperatorLCAO<TK, TR>>::cal_occ(const int& iat1,
                                                         const int& iat2,
                                                         const Parallel_Orbitals* paraV,
                                                         const std::unordered_map<int, std::vector<double>>& nlm1_all,
                                                         const std::unordered_map<int, std::vector<double>>& nlm2_all,
                                                         const double* dm_pointer,
                                                         std::vector<double>& occ)
{

    // npol is the number of polarizations,
    // 1 for non-magnetic (one Hamiltonian matrix only has spin-up or spin-down),
    // 2 for magnetic (one Hamiltonian matrix has both spin-up and spin-down)
    const int npol = this->ucell->get_npol();
    // ---------------------------------------------
    // calculate the Nonlocal matrix for each pair of orbitals
    // ---------------------------------------------
    auto row_indexes = paraV->get_indexes_row(iat1);
    auto col_indexes = paraV->get_indexes_col(iat2);
    const int m_size = int(sqrt(occ.size()) / npol);
    const int m_size2 = m_size * m_size;
#ifdef __DEBUG
    assert(m_size * m_size == occ.size());
#endif
    // step_trace = 0 for NSPIN=1,2; ={0, 1, local_col, local_col+1} for NSPIN=4
    std::vector<int> step_trace(npol * npol, 0);
    for (int is = 0; is < npol; is++)
    {
        for (int is2 = 0; is2 < npol; is2++)
        {
            step_trace[is * npol + is2] = paraV->get_ncol_atom(iat2) * is + is2;
        }
    }
    // calculate the local matrix
    for (int iw1l = 0; iw1l < row_indexes.size(); iw1l += npol)
    {
        const std::vector<double>& nlm1 = nlm1_all.find(row_indexes[iw1l])->second;
        for (int iw2l = 0; iw2l < col_indexes.size(); iw2l += npol)
        {
            const std::vector<double>& nlm2 = nlm2_all.find(col_indexes[iw2l])->second;
#ifdef __DEBUG
            assert(nlm1.size() == nlm2.size());
#endif
            for (int is1 = 0; is1 < npol; ++is1)
            {
                for (int is2 = 0; is2 < npol; ++is2)
                {
                    for (int m1 = 0; m1 < m_size; ++m1)
                    {
                        for (int m2 = 0; m2 < m_size; ++m2)
                        {
                            occ[m1 * m_size + m2 + (is1 * npol + is2) * m_size2]
                                += nlm1[m1] * nlm2[m2] * dm_pointer[step_trace[is1 * npol + is2]];
                        }
                    }
                }
            }
            dm_pointer += npol;
        }
        dm_pointer += (npol - 1) * col_indexes.size();
    }
}

template <typename TK, typename TR>
void hamilt::DFTU<hamilt::OperatorLCAO<TK, TR>>::transfer_vu(std::vector<double>& vu_tmp, std::vector<TR>& vu)
{
#ifdef __DEBUG
    assert(vu.size() == vu_tmp.size());
#endif
    for (int i = 0; i < vu_tmp.size(); i++)
    {
        vu[i] = vu_tmp[i];
    }
}

template <>
void hamilt::DFTU<hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>>::transfer_vu(
    std::vector<double>& vu_tmp,
    std::vector<std::complex<double>>& vu)
{
#ifdef __DEBUG
    assert(vu.size() == vu_tmp.size());
#endif

    // Pauli-to-spinor conversion for DFT+U potential:
    // V = V_0*I + V_x*sigma_x + V_y*sigma_y + V_z*sigma_z
    // sigma_y = [[0,-i],[i,0]], so:
    //   V_{up,up}   = 0.5*(V_0 + V_z)
    //   V_{down,down} = 0.5*(V_0 - V_z)
    //   V_{up,down}  = 0.5*(V_x - i*V_y)  <- note: minus sign from sigma_y
    //   V_{down,up}  = 0.5*(V_x + i*V_y)  <- note: plus sign from sigma_y
    // This is consistent with the convention in gint_common.cpp merge_hr_part_to_hR().
    const int m_size = int(sqrt(vu.size()) / 2);
    const int m_size2 = m_size * m_size;
    vu.resize(vu_tmp.size());
    for (int m1 = 0; m1 < m_size; m1++)
    {
        for (int m2 = 0; m2 < m_size; m2++)
        {
            int index[4];
            index[0] = m1 * m_size + m2;
            index[1] = m1 * m_size + m2 + m_size2;
            index[2] = m2 * m_size + m1 + m_size2 * 2;
            index[3] = m2 * m_size + m1 + m_size2 * 3;
            vu[index[0]] = 0.5 * (vu_tmp[index[0]] + vu_tmp[index[3]]);
            vu[index[3]] = 0.5 * (vu_tmp[index[0]] - vu_tmp[index[3]]);
            vu[index[1]] = 0.5 * (vu_tmp[index[1]] - std::complex<double>(0.0, 1.0) * vu_tmp[index[2]]);
            vu[index[2]] = 0.5 * (vu_tmp[index[1]] + std::complex<double>(0.0, 1.0) * vu_tmp[index[2]]);
        }
    }
}

template <typename TK, typename TR>
void hamilt::DFTU<hamilt::OperatorLCAO<TK, TR>>::cal_v_of_u(const std::vector<double>& occ,
                                                            const int m_size,
                                                            const double u_value,
                                                            double* vu,
                                                            double& eu)
{
    // calculate the local matrix
    int spin_fold = occ.size() / m_size / m_size;
    if (spin_fold < 4) {
        for (int is = 0; is < spin_fold; ++is)
        {
            int start = is * m_size * m_size;
            for (int m1 = 0; m1 < m_size; m1++)
            {
                for (int m2 = 0; m2 < m_size; m2++)
                {
                    vu[start + m1 * m_size + m2] = u_value * (0.5 * (m1 == m2) - occ[start + m2 * m_size + m1]);
                    eu += u_value * 0.5 * occ[start + m2 * m_size + m1] * occ[start + m1 * m_size + m2];
                }
            }
        }
    } else
    {
        for (int m1 = 0; m1 < m_size; m1++)
        {
            for (int m2 = 0; m2 < m_size; m2++)
            {
                vu[m1 * m_size + m2] = u_value * (1.0 * (m1 == m2) - occ[m2 * m_size + m1]);
                eu += u_value * 0.25 * occ[m2 * m_size + m1] * occ[m1 * m_size + m2];
            }
        }
        for (int is = 1; is < spin_fold; ++is)
        {
            int start = is * m_size * m_size;
            for (int m1 = 0; m1 < m_size; m1++)
            {
                for (int m2 = 0; m2 < m_size; m2++)
                {
                    vu[start + m1 * m_size + m2] = u_value * (0 - occ[start + m2 * m_size + m1]);
                    eu += u_value * 0.25 * occ[start + m2 * m_size + m1] * occ[start + m1 * m_size + m2];
                }
            }
        }
    }
}

#include "dftu_force_stress.hpp"

template class hamilt::DFTU<hamilt::OperatorLCAO<double, double>>;
template class hamilt::DFTU<hamilt::OperatorLCAO<std::complex<double>, double>>;
template class hamilt::DFTU<hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>>;
