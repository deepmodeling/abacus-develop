/// @file dftu_nao_op_kernel.cpp
/// @brief Per-pair kernel functions for the DFTU LCAO operator, split from
///        dftu_nao_op.cpp to keep the main file under the 500-line limit.
///        These are template member-function definitions that are explicitly
///        instantiated at the end of this translation unit.
#include "dftu_nao_op.h"

#include "source_base/timer.h"
#include "source_base/tool_title.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "source_estate/module_dm/density_matrix.h"
#include "source_lcao/module_operator_lcao/operator_lcao.h"
#include "source_base/parallel_reduce.h"

// Include the free function implementations for force/stress in real space
#include "dftu_nao_fs_r.h"

// The real-space atom-pair kernels cal_hr_ijr<TR>() and cal_occ_ijr() live
// in dftu_nao_ijr.h as free functions of namespace DFTU_LCAO.

// transfer_pot_onsite (generic: identity copy)
template <typename TK, typename TR>
void hamilt::DFTU<hamilt::OperatorLCAO<TK, TR>>::transfer_pot_onsite(std::vector<double>& pot_onsite_tmp, std::vector<TR>& pot_onsite)
{
#ifdef __DEBUG
    assert(pot_onsite.size() == pot_onsite_tmp.size());
#endif
    for (int i = 0; i < pot_onsite_tmp.size(); i++)
    {
        pot_onsite[i] = pot_onsite_tmp[i];
    }
}

// transfer_pot_onsite (noncollinear specialization: Pauli-to-spinor)
template <>
void hamilt::DFTU<hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>>::transfer_pot_onsite(
    std::vector<double>& pot_onsite_tmp,
    std::vector<std::complex<double>>& pot_onsite)
{
#ifdef __DEBUG
    assert(pot_onsite.size() == pot_onsite_tmp.size());
#endif

    // Pauli-to-spinor conversion for DFT+U potential:
    // V = V_0*I + V_x*sigma_x + V_y*sigma_y + V_z*sigma_z
    const int m_size = int(sqrt(pot_onsite.size()) / 2);
    const int m_size2 = m_size * m_size;
    pot_onsite.resize(pot_onsite_tmp.size());
    for (int m1 = 0; m1 < m_size; m1++)
    {
        for (int m2 = 0; m2 < m_size; m2++)
        {
            int index[4];
            index[0] = m1 * m_size + m2;
            index[1] = m1 * m_size + m2 + m_size2;
            index[2] = m2 * m_size + m1 + m_size2 * 2;
            index[3] = m2 * m_size + m1 + m_size2 * 3;
            pot_onsite[index[0]] = 0.5 * (pot_onsite_tmp[index[0]] + pot_onsite_tmp[index[3]]);
            pot_onsite[index[3]] = 0.5 * (pot_onsite_tmp[index[0]] - pot_onsite_tmp[index[3]]);
            pot_onsite[index[1]] = 0.5 * (pot_onsite_tmp[index[1]] - std::complex<double>(0.0, 1.0) * pot_onsite_tmp[index[2]]);
            pot_onsite[index[2]] = 0.5 * (pot_onsite_tmp[index[1]] + std::complex<double>(0.0, 1.0) * pot_onsite_tmp[index[2]]);
        }
    }
}

// cal_pot_onsite: pot = U * (1/2*delta - occ), energy = U * 1/2 * occ * occ
template <typename TK, typename TR>
void hamilt::DFTU<hamilt::OperatorLCAO<TK, TR>>::cal_pot_onsite(const std::vector<double>& occ,
                                                            const int m_size,
                                                            const double u_value,
                                                            double* pot_onsite,
                                                            double& eu)
{
    int spin_fold = occ.size() / m_size / m_size;
    if (spin_fold < 4) {
        for (int is = 0; is < spin_fold; ++is)
        {
            int start = is * m_size * m_size;
            for (int m1 = 0; m1 < m_size; m1++)
            {
                for (int m2 = 0; m2 < m_size; m2++)
                {
                    pot_onsite[start + m1 * m_size + m2] = u_value * (0.5 * (m1 == m2) - occ[start + m2 * m_size + m1]);
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
                pot_onsite[m1 * m_size + m2] = u_value * (1.0 * (m1 == m2) - occ[m2 * m_size + m1]);
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
                    pot_onsite[start + m1 * m_size + m2] = u_value * (0 - occ[start + m2 * m_size + m1]);
                    eu += u_value * 0.25 * occ[start + m2 * m_size + m1] * occ[start + m1 * m_size + m2];
                }
            }
        }
    }
}

// explicit template instantiation (matches dftu_nao_op.cpp)
template class hamilt::DFTU<hamilt::OperatorLCAO<double, double>>;
template class hamilt::DFTU<hamilt::OperatorLCAO<std::complex<double>, double>>;
template class hamilt::DFTU<hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>>;
