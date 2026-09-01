#ifndef DFTU_LCAO_ENERGY_H
#define DFTU_LCAO_ENERGY_H

class Plus_U;
class UnitCell;

#ifdef __LCAO
namespace DFTU_LCAO {

/**
 * @brief DFT+U energy correction with the double-counting term subtracted.
 *
 * Computes energy_u from occ_mat and the onsite potential, then writes the
 * result back to dftu via set_energy. The spin channel count is read from
 * the global input rather than a Plus_U member.
 *
 * @param dftu  Plus_U state (mutable: set_energy is called at the end)
 * @param ucell unit cell
 */
void cal_energy_correction(Plus_U& dftu, const UnitCell& ucell);

} // namespace DFTU_LCAO
#endif

#endif
