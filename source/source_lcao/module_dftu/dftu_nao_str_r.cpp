/// @file dftu_nao_str_r.cpp
/// @brief DFT+U stress calculation in real space (r-space) - implementation
///
/// See dftu_nao_str_r.h for the mathematical formula and detailed documentation.

#include "dftu_nao_str_r.h"
#include "dftu_nao_op.h"
#include "source_base/timer.h"

namespace hamilt
{

template <typename TK, typename TR>
void cal_str_IJR_nao_r(const DFTU<OperatorLCAO<TK, TR>>* dftu_op,
                      const int& iat1,
                      const int& iat2,
                      const Parallel_Orbitals* pv,
                      const std::unordered_map<int, std::vector<double>>& nlm1_all,
                      const std::unordered_map<int, std::vector<double>>& nlm2_all,
                      const std::vector<double>& pot_onsite_in,
                      const hamilt::BaseMatrix<double>** dmR_pointer,
                      const int nspin,
                      const ModuleBase::Vector3<double>& dis1,
                      const ModuleBase::Vector3<double>& dis2,
                      double* stress)
{
    // npol is the number of polarizations,
    // 1 for non-magnetic (one Hamiltonian matrix only has spin-up or spin-down),
    // 2 for magnetic (one Hamiltonian matrix has both spin-up and spin-down)
    const int npol = dftu_op->get_ucell()->get_npol();

    // ---------------------------------------------
    // calculate the Nonlocal matrix for each pair of orbitals
    // ---------------------------------------------
    auto row_indexes = pv->get_indexes_row(iat1);
    auto col_indexes = pv->get_indexes_col(iat2);
    const int m_size = int(sqrt(pot_onsite_in.size() / nspin));
    const int m_size2 = m_size * m_size;

    // step_trace = 0 for NSPIN=1,2; ={0, 1, local_col, local_col+1} for NSPIN=4
    std::vector<int> step_trace(npol * npol, 0);

    if (npol == 2)
    {
        step_trace[1] = 1;
        step_trace[2] = col_indexes.size();
        step_trace[3] = col_indexes.size() + 1;
    }

    // calculate the local matrix
    for (int is = 0; is < nspin; is++)
    {
        const int is0 = nspin == 2 ? is : 0;
        const int step_is = nspin == 4 ? is : 0;
        const double* dm_pointer = dmR_pointer[is0]->get_pointer();
        for (int iw1l = 0; iw1l < row_indexes.size(); iw1l += npol)
        {
            const std::vector<double>& nlm1 = nlm1_all.find(row_indexes[iw1l])->second;
            for (int iw2l = 0; iw2l < col_indexes.size(); iw2l += npol)
            {
                const std::vector<double>& nlm2 = nlm2_all.find(col_indexes[iw2l])->second;
#ifdef __DEBUG
                assert(nlm1.size() == nlm2.size());
#endif
                for (int m1 = 0; m1 < m_size; m1++)
                {
                    for (int m2 = 0; m2 < m_size; m2++)
                    {
                        double tmp = pot_onsite_in[m1 * m_size + m2 + is * m_size2]
                                     * dm_pointer[step_trace[step_is]];
                        // Voigt notation: stress[0]=xx, stress[1]=xy, stress[2]=xz,
                        //                 stress[3]=yy, stress[4]=yz, stress[5]=zz
                        stress[0] += tmp * (nlm1[m1 + m_size] * dis1.x * nlm2[m2]
                                     + nlm1[m1] * nlm2[m2 + m_size] * dis2.x);
                        stress[1] += tmp * (nlm1[m1 + m_size] * dis1.y * nlm2[m2]
                                     + nlm1[m1] * nlm2[m2 + m_size] * dis2.y);
                        stress[2] += tmp * (nlm1[m1 + m_size] * dis1.z * nlm2[m2]
                                     + nlm1[m1] * nlm2[m2 + m_size] * dis2.z);

                        stress[3] += tmp * (nlm1[m1 + m_size * 2] * dis1.y * nlm2[m2]
                                        + nlm1[m1] * nlm2[m2 + m_size * 2] * dis2.y);
                        stress[4] += tmp * (nlm1[m1 + m_size * 2] * dis1.z * nlm2[m2]
                                        + nlm1[m1] * nlm2[m2 + m_size * 2] * dis2.z);
                        stress[5] += tmp * (nlm1[m1 + m_size * 3] * dis1.z * nlm2[m2]
                                        + nlm1[m1] * nlm2[m2 + m_size * 3] * dis2.z);
                    }
                }
                dm_pointer += npol;
            }
            dm_pointer += (npol - 1) * col_indexes.size();
        }
    }
}

// explicit template instantiation
template void cal_str_IJR_nao_r<double, double>(
    const DFTU<OperatorLCAO<double, double>>* dftu_op,
    const int& iat1, const int& iat2,
    const Parallel_Orbitals* pv,
    const std::unordered_map<int, std::vector<double>>& nlm1_all,
    const std::unordered_map<int, std::vector<double>>& nlm2_all,
    const std::vector<double>& pot_onsite_in,
    const hamilt::BaseMatrix<double>** dmR_pointer,
    const int nspin,
    const ModuleBase::Vector3<double>& dis1,
    const ModuleBase::Vector3<double>& dis2,
    double* stress);

template void cal_str_IJR_nao_r<std::complex<double>, double>(
    const DFTU<OperatorLCAO<std::complex<double>, double>>* dftu_op,
    const int& iat1, const int& iat2,
    const Parallel_Orbitals* pv,
    const std::unordered_map<int, std::vector<double>>& nlm1_all,
    const std::unordered_map<int, std::vector<double>>& nlm2_all,
    const std::vector<double>& pot_onsite_in,
    const hamilt::BaseMatrix<double>** dmR_pointer,
    const int nspin,
    const ModuleBase::Vector3<double>& dis1,
    const ModuleBase::Vector3<double>& dis2,
    double* stress);

template void cal_str_IJR_nao_r<std::complex<double>, std::complex<double>>(
    const DFTU<OperatorLCAO<std::complex<double>, std::complex<double>>>* dftu_op,
    const int& iat1, const int& iat2,
    const Parallel_Orbitals* pv,
    const std::unordered_map<int, std::vector<double>>& nlm1_all,
    const std::unordered_map<int, std::vector<double>>& nlm2_all,
    const std::vector<double>& pot_onsite_in,
    const hamilt::BaseMatrix<double>** dmR_pointer,
    const int nspin,
    const ModuleBase::Vector3<double>& dis1,
    const ModuleBase::Vector3<double>& dis2,
    double* stress);

} // namespace hamilt
