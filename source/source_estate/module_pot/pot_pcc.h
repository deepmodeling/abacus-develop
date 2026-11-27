#ifndef POT_PCC_H
#define POT_PCC_H

#include "pot_base.h"
#include "source_cell/cell_geometry.h"
#include "source_estate/pcc_0d.h"
#include "source_estate/pcc_2d.h"

namespace elecstate
{

/**
 * @brief Self-consistent vacuum point-counter-charge (PCC) correction.
 *
 * Owns the result of its last density update. The correction removes the
 * spurious interaction of the total charge with its periodic images by
 * expanding that charge to second order in its multipoles about the
 * mass-weighted ionic center (charge q, dipole d, second moment Q). The
 * charge seen here is ionic + electronic only; any solvent polarization charge
 * must be corrected by its own component.
 * In Hartree atomic units, for a unit positive test charge at r (or at the
 * normal coordinate z for slabs):
 *
 *   0D (cubic cell of edge L, volume V = L^3, Madelung constant alpha of the
 *   simple cubic lattice):
 *     v(r) = alpha q / L - 2 pi / (3 V) * (q r^2 - 2 d.r + Q)
 *
 *   2D (periodic area A, normal period L):
 *     v(z) = -pi q L / (3 A) - 2 pi / (A L) * (q z^2 - 2 d_z z + Q_zz)
 *
 * The correction energy is E = 1/2 * integral rho_total(r) v(r) dr and the
 * forces follow from the gradient of v at the ionic positions. The kernels are
 * implemented in pcc_0d.cpp and pcc_2d.cpp; see pcc_2d.h for the gauge of the
 * charged-slab constant.
 *
 * References:
 *   - O. Andreussi and N. Marzari, Phys. Rev. B 90, 245101 (2014)
 *     (parabolic 0D and 2D point-counter-charge corrections).
 *   - I. Dabo, B. Kozinsky, N. E. Singh-Miller and N. Marzari,
 *     Phys. Rev. B 77, 115139 (2008) (density-countercharge corrections).
 *   - M. J. Rutter, Electron. Struct. 3, 015002 (2021) (charged slabs: the
 *     zero-mean periodic potential of a charged plane fixes the 2D constant).
 *
 * The implementation follows the PCC corrections of the ENVIRON library
 * (www.quantum-environ.org).
 */
class PotPcc : public PotBase
{
  public:
    enum class Dimension { molecule, slab };

    /// electron_count is the expected number of electrons (INPUT nelec); a
    /// grid charge that differs from it by more than 1e-6 e is reported in the
    /// warning log, and the correction uses the grid charge.
    PotPcc(const ModulePW::PW_Basis* basis, double electron_count);

    PotPcc(const ModulePW::PW_Basis* basis, Dimension dimension, int open_axis, double electron_count);
    static void validate_kpoints(const std::vector<ModuleBase::Vector3<double>>& points, int count, int open_axis);

    void cal_v_eff(const Charge* charge, const UnitCell* cell, ModuleBase::matrix& potential) override;
    double get_energy() const override;
    const std::vector<double>& electron_potential() const;
    void add_force(const UnitCell& cell, ModuleBase::matrix& force) const;

  private:
    void prepare_ions(const UnitCell& cell, ChargeMoments& ionic_moments);
    ChargeMoments collect_electrons(const Charge& charge,
                                    const UnitCell& cell,
                                    std::vector<ModuleBase::Vector3<double>>& positions) const;

    ModuleBase::Vector3<double> relative_position(const ModuleBase::Vector3<double>& position) const;
    double correction_energy() const;
    double correction_potential(const ModuleBase::Vector3<double>& position) const;
    ModuleBase::Vector3<double> correction_force(double charge, const ModuleBase::Vector3<double>& position) const;

    const Dimension dimension_;
    const int open_axis_;
    const double electron_count_;
    unitcell::SlabCell slab_;
    Pcc2dParameters slab_parameters_;
    unitcell::OrthogonalCell geometry_;
    Pcc0dParameters parameters_;
    ChargeMoments moments_;
    std::vector<ModuleBase::Vector3<double>> ionic_positions_;
    std::vector<double> ionic_charges_;
    std::vector<double> electron_potential_;
    double energy_rydberg_ = 0.0;
};

} // namespace elecstate

#endif
