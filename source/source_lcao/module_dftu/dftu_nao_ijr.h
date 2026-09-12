#ifndef DFTU_LCAO_IJR_H
#define DFTU_LCAO_IJR_H

#include "source_basis/module_ao/parallel_orbitals.h"

#include <cassert>
#include <cmath>
#include <complex>
#include <unordered_map>
#include <vector>

namespace DFTU_LCAO
{

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
    auto row_indexes = pv.get_indexes_row(iat1);
    auto col_indexes = pv.get_indexes_col(iat2);
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
    auto row_indexes = pv.get_indexes_row(iat1);
    auto col_indexes = pv.get_indexes_col(iat2);
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

} // namespace DFTU_LCAO

#endif
