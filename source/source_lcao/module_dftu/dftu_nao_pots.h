#ifndef DFTU_LCAO_POTS_H
#define DFTU_LCAO_POTS_H

#include <complex>

class Plus_U;
class UnitCell;
class Parallel_Orbitals;

#ifdef __LCAO
namespace DFTU_LCAO {

/**
 * @brief one-body effective onsite potential element for a given (m0,m1) pair.
 *
 * Dispatches on cal_type; only case 3 (simplified formalism with FLL double
 * counting) is currently implemented, the other cases return 0.
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
double get_onsite_pot(const Plus_U& dftu,
                      const int T,
                      const int iat,
                      const int L,
                      const int N,
                      const int spin,
                      const int m0,
                      const int m1,
                      const bool new_occ_mat);

/**
 * @brief onsite effective potential matrix (complex) in the local orbital basis.
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
 * @param npol        number of polarizations
 */
void pot_onsite_complex(const Plus_U& dftu,
                        const UnitCell& ucell,
                        const Parallel_Orbitals* pv,
                        const int spin,
                        const bool new_occ_mat,
                        std::complex<double>* pot_onsite,
                        const int npol);

/**
 * @brief onsite effective potential matrix (real) in the local orbital basis.
 *
 * Real-valued counterpart of pot_onsite_complex.
 */
void pot_onsite_real(const Plus_U& dftu,
                     const UnitCell& ucell,
                     const Parallel_Orbitals* pv,
                     const int spin,
                     const bool new_occ_mat,
                     double* pot_onsite,
                     const int npol);

} // namespace DFTU_LCAO
#endif

#endif
