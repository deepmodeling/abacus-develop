#ifndef DFTU_LCAO_ENERGY_H
#define DFTU_LCAO_ENERGY_H

class Plus_U_Base;
class UnitCell;

namespace DFTU_LCAO {

/**
 * @brief DFT+U energy correction with the double-counting term subtracted:
 *        E = E_U - E_dc, with
 *          E_U  = (U_eff / 2) * sum_{m,m'} occ(m,m') * (delta_{m,m'} - occ(m',m))
 *          E_dc = sum_{m,m'} onsite_pot(m,m') * occ(m',m)
 *
 * @param dftu  Plus_U_Base state (mutable: set_energy is called at the end)
 * @param ucell unit cell
 * @param nspin number of spin channels (1, 2, or 4); sourced by the caller
 *        from the input parameter to avoid a PARAM read here
 */
void cal_energy_correction(Plus_U_Base& dftu,
                           const UnitCell& ucell,
                           int nspin);

/**
 * @brief Accumulate the DFT+U energy term (U_eff / 2) * (n - n^2) for one
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
                                 int n);

/**
 * @brief Accumulate the double-counting correction energy_dc for one
 *        (T, iat, l, n=0) channel:
 *        E_dc = sum_{m1,ipol1,m2,ipol2} onsite_pot(m1,ipol1;m2,ipol2) * occ(m2,ipol2;m1,ipol1)
 *        Returns the per-atom contribution to energy_dc.
 */
double calc_energy_dc_block(const Plus_U_Base& dftu,
                            int T,
                            int iat,
                            int l,
                            int n,
                            int nspin);

} // namespace DFTU_LCAO

#endif
