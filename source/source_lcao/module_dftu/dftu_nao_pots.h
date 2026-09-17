#ifndef DFTU_LCAO_POTS_H
#define DFTU_LCAO_POTS_H

#include <cmath>
#include <complex>
#include <vector>

class Plus_U_Base;
class UnitCell;
class Parallel_Orbitals;

namespace DFTU_LCAO {

/**
 * @brief On-site potential and Hubbard energy for one correlated shell:
 *   pot_onsite(m,m') = U_eff * (0.5 * delta_{m,m'} - occ(m,m'))
 *   EU = (U_eff / 2) * sum_{m,m'} occ(m,m') * (delta_{m,m'} - occ(m',m))
 *
 * Non-template core extracted from DFTU_onsite<OperatorLCAO<TK,TR>>::cal_pot_onsite.
 *
 * @param occ         flattened occupation matrix, size m_size*m_size*spin_fold
 * @param m_size      number of magnetic quantum states (2*l+1)
 * @param u_value     effective U parameter
 * @param pot_onsite  output potential, same size as occ
 * @param eu          output Hubbard energy contribution
 */
void cal_pot_onsite(const std::vector<double>& occ, const int m_size, const double u_value,
                    double* pot_onsite, double& eu);

/**
 * @brief transfer pot_onsite format from Pauli matrix storage to the target
 *        matrix-element type used by the Hamiltonian.
 *
 * Generic overload performs an identity copy. The complex specialization
 * converts Pauli blocks to spinor blocks for the non-collinear case.
 *
 * @param pot_onsite_tmp input onsite potential in Pauli-block storage
 * @param pot_onsite     output onsite potential in target storage
 */
template <typename TR>
void transfer_pot_onsite(const std::vector<double>& pot_onsite_tmp,
                         std::vector<TR>& pot_onsite)
{
    pot_onsite.resize(pot_onsite_tmp.size());
    for (int i = 0; i < static_cast<int>(pot_onsite_tmp.size()); i++)
    {
        pot_onsite[i] = static_cast<TR>(pot_onsite_tmp[i]);
    }
}

template <>
inline void transfer_pot_onsite<std::complex<double>>(const std::vector<double>& pot_onsite_tmp,
                                                      std::vector<std::complex<double>>& pot_onsite)
{
    const int m_size = int(sqrt(pot_onsite_tmp.size()) / 2);
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

/**
 * @brief one-body effective onsite potential element for a given (m0,m1) pair.
 *
 * Dispatches on UForm; only dud_fll (Dudarev simplified formalism with FLL
 * double counting) is currently implemented, the other forms return 0.
 *
 * @param dftu        Plus_U state providing U/J values and occupation matrices
 * @param T           atom type
 * @param iat         global atom index
 * @param L           angular momentum
 * @param N           radial index
 * @param spin        spin channel
 * @param m0          first magnetic quantum index (packed with polarization)
 * @param m1          second magnetic quantum index (packed with polarization)
 * @param new_occ_mat if true use occ_mat, otherwise use occ_mat_save
 * @return            onsite potential matrix element
 */
double get_onsite_pot(const Plus_U_Base& dftu,
                      const int T,
                      const int iat,
                      const int L,
                      const int N,
                      const int spin,
                      const int m0,
                      const int m1,
                      const bool new_occ_mat);

/**
 * @brief Calculate onsite effective potential matrix in the local orbital basis.
 *
 * Fills pot_onsite (length pv->nloc) with the onsite potential elements
 * projected onto the local orbital indices.
 *
 * @param dftu        Plus_U state
 * @param ucell       unit cell
 * @param pv          parallel orbitals descriptor
 * @param spin        spin channel
 * @param new_occ_mat if true use occ_mat, otherwise use occ_mat_save
 * @param pot_onsite  output buffer (length pv->nloc)
 */
template <typename T>
void cal_pot_onsite(const Plus_U_Base& dftu,
                    const UnitCell& ucell,
                    const Parallel_Orbitals* pv,
                    const int spin,
                    const bool new_occ_mat,
                    T* pot_onsite);

// Explicit instantiations
extern template void cal_pot_onsite<double>(const Plus_U_Base& dftu,
                                            const UnitCell& ucell,
                                            const Parallel_Orbitals* pv,
                                            const int spin,
                                            const bool new_occ_mat,
                                            double* pot_onsite);
extern template void cal_pot_onsite<std::complex<double>>(const Plus_U_Base& dftu,
                                                          const UnitCell& ucell,
                                                          const Parallel_Orbitals* pv,
                                                          const int spin,
                                                          const bool new_occ_mat,
                                                          std::complex<double>* pot_onsite);

/// DFT+U effective potential in the LCAO basis, with V = pot_onsite (the
/// on-site Hubbard correction potential in the full basis) and S the overlap
/// matrix. Returns the symmetrized k-space U-term potential:
///   pot_uterm = (V*S + (V*S)^T) / 2

/// @brief Compute the LCAO-basis U-term effective potential matrix.
///
/// @tparam T           matrix element type (double or std::complex<double>)
/// @param dftu         Plus_U state
/// @param ucell        unit cell
/// @param pv           parallel orbitals descriptor
/// @param spin         spin channel (isk[ik] from caller)
/// @param pot_uterm    output buffer (length pv->nloc)
/// @param sk           overlap matrix in k-space (length pv->nloc)
template <typename T>
void cal_pot_uterm(Plus_U_Base& dftu,
                   const UnitCell& ucell,
                   const Parallel_Orbitals* pv,
                   const int spin,
                   T* pot_uterm,
                   const T* sk);

// Explicit instantiations
extern template void cal_pot_uterm<double>(Plus_U_Base& dftu,
                                           const UnitCell& ucell,
                                           const Parallel_Orbitals* pv,
                                           const int spin,
                                           double* pot_uterm,
                                           const double* sk);
extern template void cal_pot_uterm<std::complex<double>>(Plus_U_Base& dftu,
                                                         const UnitCell& ucell,
                                                         const Parallel_Orbitals* pv,
                                                         const int spin,
                                                         std::complex<double>* pot_uterm,
                                                         const std::complex<double>* sk);

} // namespace DFTU_LCAO

#endif
