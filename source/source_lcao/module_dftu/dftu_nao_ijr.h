#ifndef DFTU_LCAO_IJR_H
#define DFTU_LCAO_IJR_H

#include "source_base/parallel_reduce.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h" // AdjacentAtomInfo (complete type needed)
#include "source_cell/unitcell.h"
#include "source_hamilt/module_hcontainer/hcontainer.h"
#include "source_pw/module_pwdft/dftu_base.h"

#include <cassert>
#include <cmath>
#include <complex>
#include <unordered_map>
#include <vector>

namespace DFTU_LCAO
{

/// The <phi|alpha^I> overlap values for all
/// [atoms][neighbors][orb_index(iw) in NAOs][m of target_l in Projectors]
using NlmTot = std::vector<std::vector<std::unordered_map<int, std::vector<double>>>>;

/**
 * @brief load the flattened occupation matrix of one Hubbard atom from
 *        pre-read occ_mat data.
 *
 * @param dftu         Plus_U state providing the stored occupation matrix
 * @param iat0         global atom index of the Hubbard atom
 * @param target_L     angular momentum channel of the correlated shell
 * @param nspin        number of spin components
 * @param current_spin active spin channel for nspin=2
 * @param occ          output flattened occupation matrix
 */
inline void load_occ_from_file(const Plus_U_Base& dftu,
                               const int iat0,
                               const int target_L,
                               const int nspin,
                               const int current_spin,
                               std::vector<double>& occ)
{
    if (nspin == 4)
    {
        dftu.occmat().get_flat(iat0, target_L, occ);
    }
    else
    {
        for (int i = 0; i < static_cast<int>(occ.size()); i++)
        {
            occ[i] = dftu.occmat().get(iat0, target_L, current_spin,
                                       i / (2 * target_L + 1), i % (2 * target_L + 1));
        }
    }
}

/**
 * @brief accumulate one real-space HR atom-pair block for DFT+U:
 *        HR += <psi_I|beta_m> * pot_onsite(m,m') * <beta_m'|psi_{J,R}>
 *
 * @tparam TR           real-space Hamiltonian scalar type: double for
 *                      NSPIN=1,2 and std::complex<double> for NSPIN=4
 * @param iat1          global atom index of the row atom I
 * @param iat2          global atom index of the column atom J
 * @param npol          number of polarizations: 1 for NSPIN=1,2 and 2 for
 *                      the non-collinear case (NSPIN=4)
 * @param pv           parallel-orbitals descriptor providing local index maps
 * @param nlm1_all     <psi_I|beta_m> overlap values keyed by local row index
 * @param nlm2_all     <beta_m'|psi_J> overlap values keyed by local column index
 * @param pot_onsite   onsite potential matrix packed in npol*npol spin blocks
 * @param data_pointer pointer to the local HR matrix block; updated in place
 */
template <typename TR>
void cal_hr_ijr(const int iat1,
                const int iat2,
                const int npol,
                const Parallel_Orbitals& pv,
                const std::unordered_map<int, std::vector<double>>& nlm1_all,
                const std::unordered_map<int, std::vector<double>>& nlm2_all,
                const std::vector<TR>& pot_onsite,
                TR* data_pointer)
{
    assert(iat1 >= 0);
    assert(iat2 >= 0);
    assert(npol > 0);
    assert(data_pointer != nullptr);
    std::vector<int> row_indexes = pv.get_indexes_row(iat1);
    std::vector<int> col_indexes = pv.get_indexes_col(iat2);
    const int m_size = int(sqrt(pot_onsite.size()) / npol);
    // step_trace = 0 for NSPIN=1,2; ={0, 1, local_col, local_col+1} for NSPIN=4
    std::vector<int> step_trace(npol * npol, 0);
    for (int is = 0; is < npol; is++)
    {
        for (int is2 = 0; is2 < npol; is2++)
        {
            step_trace[is * npol + is2] = pv.get_ncol_atom(iat2) * is + is2;
        }
    }
    // calculate the local matrix
    for (int iw1l = 0; iw1l < int(row_indexes.size()); iw1l += npol)
    {
        const std::vector<double>& nlm1 = nlm1_all.find(row_indexes[iw1l])->second;
        for (int iw2l = 0; iw2l < int(col_indexes.size()); iw2l += npol)
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
                        nlm_tmp += nlm1[m1] * nlm2[m2] * pot_onsite[m1 * m_size + m2 + start];
                    }
                }
                data_pointer[step_trace[is]] += nlm_tmp;
            }
            data_pointer += npol;
        }
        data_pointer += (npol - 1) * col_indexes.size();
    }
}

/**
 * @brief accumulate one atom-pair contribution to the DFT+U occupation matrix:
 *        occ_mm' += sum_R DMR(I,J,R) * <phi_0|alpha^I_m> * <alpha^J_m'|phi_R>
 *
 * @param iat1        global atom index of the row atom I
 * @param iat2        global atom index of the column atom J
 * @param npol        number of polarizations: 1 for NSPIN=1,2 and 2 for
 *                    the non-collinear case (NSPIN=4)
 * @param pv          parallel-orbitals descriptor providing local index maps
 * @param nlm1_all    <phi_0|alpha^I_m> overlap values keyed by local row index
 * @param nlm2_all    <alpha^J_m'|phi_R> overlap values keyed by local column index
 * @param dm_pointer  pointer to the local real-space DMR block of (I,J,R)
 * @param occ         occupation matrix packed in npol*npol spin blocks;
 *                    updated in place
 */
inline void cal_occ_ijr(const int iat1,
                        const int iat2,
                        const int npol,
                        const Parallel_Orbitals& pv,
                        const std::unordered_map<int, std::vector<double>>& nlm1_all,
                        const std::unordered_map<int, std::vector<double>>& nlm2_all,
                        const double* dm_pointer,
                        std::vector<double>& occ)
{
    assert(iat1 >= 0);
    assert(iat2 >= 0);
    assert(npol > 0);
    assert(dm_pointer != nullptr);
    std::vector<int> row_indexes = pv.get_indexes_row(iat1);
    std::vector<int> col_indexes = pv.get_indexes_col(iat2);
    const int m_size = int(sqrt(occ.size()) / npol);
    const int m_size2 = m_size * m_size;
#ifdef __DEBUG
    assert(m_size2 * npol * npol == int(occ.size()));
#endif
    // step_trace = 0 for NSPIN=1,2; ={0, 1, local_col, local_col+1} for NSPIN=4
    std::vector<int> step_trace(npol * npol, 0);
    for (int is = 0; is < npol; is++)
    {
        for (int is2 = 0; is2 < npol; is2++)
        {
            step_trace[is * npol + is2] = pv.get_ncol_atom(iat2) * is + is2;
        }
    }
    for (int iw1l = 0; iw1l < int(row_indexes.size()); iw1l += npol)
    {
        const std::vector<double>& nlm1 = nlm1_all.find(row_indexes[iw1l])->second;
        for (int iw2l = 0; iw2l < int(col_indexes.size()); iw2l += npol)
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

/**
 * @brief walk all (ad1, ad2) adjacent-atom pairs of one Hubbard atom iat0
 *        and invoke @p body for each pair that has a matching matrix block.
 *
 * This is the common iteration skeleton shared by the HR accumulation
 * (accumulate_hr_for_iat0), occupation-matrix computation
 * (compute_occ_from_dmr), and the force/stress real-space implementation
 * (cal_fs_nao_r_impl). The body receives (iat1, iat2, R_vector, nlm1, nlm2)
 * and is responsible for the actual physics.
 *
 * @param ucell   [in] unit cell (atom index maps)
 * @param iat0    [in] global atom index of the Hubbard atom
 * @param adjs    [in] adjacent atom info of the Hubbard atom
 * @param nlm_tot [in] <phi|alpha^I> overlap table for all atoms
 * @param body    [in] callable with signature
 *        void(int iat1, int iat2,
 *             const ModuleBase::Vector3<int>& R_vector,
 *             const std::unordered_map<int, std::vector<double>>& nlm1,
 *             const std::unordered_map<int, std::vector<double>>& nlm2)
 */
template <typename Body>
inline void for_adj_pair(const UnitCell& ucell,
                         const int iat0,
                         const AdjacentAtomInfo& adjs,
                         const NlmTot& nlm_tot,
                         Body&& body)
{
    for (int ad1 = 0; ad1 < adjs.adj_num + 1; ++ad1)
    {
        const int T1 = adjs.ntype[ad1];
        const int I1 = adjs.natom[ad1];
        const int iat1 = ucell.itia2iat(T1, I1);
        const ModuleBase::Vector3<int>& R_index1 = adjs.box[ad1];
        const std::unordered_map<int, std::vector<double>>& nlm1 = nlm_tot[iat0][ad1];
        for (int ad2 = 0; ad2 < adjs.adj_num + 1; ++ad2)
        {
            const int T2 = adjs.ntype[ad2];
            const int I2 = adjs.natom[ad2];
            const int iat2 = ucell.itia2iat(T2, I2);
            const std::unordered_map<int, std::vector<double>>& nlm2 = nlm_tot[iat0][ad2];
            const ModuleBase::Vector3<int>& R_index2 = adjs.box[ad2];
            ModuleBase::Vector3<int> R_vector(R_index2[0] - R_index1[0],
                                              R_index2[1] - R_index1[1],
                                              R_index2[2] - R_index1[2]);
            body(iat1, iat2, R_vector, nlm1, nlm2);
        }
    }
}

/**
 * @brief accumulate the real-space HR contributions of one Hubbard atom
 *        (iat0) from the precomputed pot_onsite:
 *        HR(I,J,R) += <phi_I|chi_m> pot_onsite(m,m') <chi_m'|phi_{J,R}>
 *
 * @tparam TR           real-space Hamiltonian scalar type: double for
 *                      NSPIN=1,2 and std::complex<double> for NSPIN=4
 * @param ucell       [in] unit cell (atom index maps, npol)
 * @param hR          [in,out] real-space Hamiltonian container; matching
 *                      atom-pair blocks are updated in place
 * @param nlm_tot     [in] <phi|alpha^I> overlap table for all atoms
 * @param iat0        [in] global atom index of the Hubbard atom
 * @param adjs        [in] adjacent atom info of the Hubbard atom
 * @param pv          [in] parallel-orbitals descriptor providing local index maps
 * @param pot_onsite  [in] onsite potential matrix packed in npol*npol spin blocks
 */
template <typename TR>
void accumulate_hr_for_iat0(const UnitCell& ucell,
                            hamilt::HContainer<TR>* hR,
                            const NlmTot& nlm_tot,
                            const int iat0,
                            const AdjacentAtomInfo& adjs,
                            const Parallel_Orbitals& pv,
                            const std::vector<TR>& pot_onsite)
{
    for_adj_pair(ucell, iat0, adjs, nlm_tot,
        [&](int iat1, int iat2,
            const ModuleBase::Vector3<int>& R_vector,
            const std::unordered_map<int, std::vector<double>>& nlm1,
            const std::unordered_map<int, std::vector<double>>& nlm2)
        {
            hamilt::BaseMatrix<TR>* tmp = hR->find_matrix(iat1, iat2, R_vector[0], R_vector[1], R_vector[2]);
            if (tmp != nullptr)
            {
#ifdef _OPENMP
#pragma omp critical(dftu_hr_update)
#endif
                {
                    cal_hr_ijr<TR>(iat1, iat2, ucell.get_npol(), pv,
                                   nlm1, nlm2, pot_onsite, tmp->get_pointer());
                }
            }
        });
}

/**
 * @brief compute the occupation matrix of one Hubbard atom (iat0) from the
 *        real-space density matrix:
 *        occ(m,m') = sum_R DMR(I,J,R) * <phi_0|chi_m(I)> * <chi_m'(J)|phi_R>
 *        then MPI-Allreduce it and store it into the Plus_U occupation matrix.
 *
 * @param ucell         [in] unit cell (atom index maps, npol)
 * @param dftu          [in,out] Plus_U state receiving the occupation matrix
 * @param iat0          [in] global atom index of the Hubbard atom
 * @param target_L      [in] angular momentum channel of the correlated shell
 * @param current_spin  [in] active spin channel (0 for nspin=1/4)
 * @param nspin         [in] number of spin components
 * @param adjs          [in] adjacent atom info of the Hubbard atom
 * @param pv            [in] parallel-orbitals descriptor providing local index maps
 * @param nlm_tot       [in] <phi|alpha^I> overlap table for all atoms
 * @param dmR_current   [in] real-space density matrix of the active spin
 * @param occ           [out] flattened occupation matrix; overwritten
 */
inline void compute_occ_from_dmr(const UnitCell& ucell,
                                 Plus_U_Base& dftu,
                                 const int iat0,
                                 const int target_L,
                                 const int current_spin,
                                 const int nspin,
                                 const AdjacentAtomInfo& adjs,
                                 const Parallel_Orbitals& pv,
                                 const NlmTot& nlm_tot,
                                 const hamilt::HContainer<double>& dmR_current,
                                 std::vector<double>& occ)
{
    for_adj_pair(ucell, iat0, adjs, nlm_tot,
        [&](int iat1, int iat2,
            const ModuleBase::Vector3<int>& R_vector,
            const std::unordered_map<int, std::vector<double>>& nlm1,
            const std::unordered_map<int, std::vector<double>>& nlm2)
        {
            const hamilt::BaseMatrix<double>* tmp
                = dmR_current.find_matrix(iat1, iat2, R_vector[0], R_vector[1], R_vector[2]);
            if (tmp != nullptr)
            {
                cal_occ_ijr(iat1, iat2, ucell.get_npol(), pv,
                            nlm1, nlm2, tmp->get_pointer(), occ);
            }
        });
    Parallel_Reduce::reduce_all(occ.data(), occ.size());
    if (nspin == 1)
    {
        for (double& v : occ) { v *= 0.5; }
    }
    dftu.occmat().set_flat(iat0, target_L, current_spin, occ);
}

} // namespace DFTU_LCAO

#endif
