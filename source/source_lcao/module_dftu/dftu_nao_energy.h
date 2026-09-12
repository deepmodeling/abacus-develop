#ifndef DFTU_LCAO_ENERGY_H
#define DFTU_LCAO_ENERGY_H

class Plus_U_Base;
class UnitCell;

namespace DFTU_LCAO {

/**
 * @brief DFT+U energy correction with the double-counting term subtracted.
 *
 * Computes energy_u from occ_mat and the onsite potential, then writes the
 * result back to dftu via set_energy.
 *
 * @param dftu  Plus_U_Base state (mutable: set_energy is called at the end)
 * @param ucell unit cell
 * @param nspin number of spin channels (1, 2, or 4); sourced by the caller
 *        from the input parameter to avoid a PARAM read here
 * @param npol Pauli-component count (1 for collinear, 2 for noncollinear);
 *        sourced by the caller from ucell.get_npol()
 */
void cal_energy_correction(Plus_U_Base& dftu,
                           const UnitCell& ucell,
                           int nspin,
                           int npol);

/**
 * @brief Accumulate the DFT+U energy term (0.5 * U * (n - n^2)) for one
 *        (T, iat, l, n=0) channel in the collinear case (nspin=1 or 2).
 *        Returns the per-atom contribution to energy_u.
 */
double calc_energy_u_collinear(const Plus_U_Base& dftu,
                               int T,
                               int iat,
                               int l,
                               int n);

/**
 * @brief Accumulate the DFT+U energy term for one (T, iat, l, n=0) channel
 *        in the noncollinear case (nspin=4). Returns the per-atom
 *        contribution to energy_u.
 */
double calc_energy_u_noncollinear(const Plus_U_Base& dftu,
                                 int T,
                                 int iat,
                                 int l,
                                 int n,
                                 int npol);

/**
 * @brief Accumulate the double-counting correction energy_dc for one
 *        (T, iat, l, n=0) channel by summing onsite_pot * occ over the
 *        (m1, ipol1, m2, ipol2) grid. Dispatches on nspin to choose the
 *        spin loop count. Returns the per-atom contribution to energy_dc.
 */
double calc_energy_dc_block(const Plus_U_Base& dftu,
                            int T,
                            int iat,
                            int l,
                            int n,
                            int nspin,
                            int npol);

} // namespace DFTU_LCAO

#endif
