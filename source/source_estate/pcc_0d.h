#ifndef PCC_0D_H
#define PCC_0D_H

#include "charge_moments.h"

namespace unitcell
{
struct OrthogonalCell;
}

namespace elecstate
{

struct Pcc0dParameters
{
    double length = 0.0; ///< Cubic cell edge in Bohr
    double madelung = 2.837297479480619;
};

/// Return false unless the OrthogonalCell has equal edges within
/// relative_tolerance. No cell geometry is stored here.
bool make_pcc_0d_parameters(const unitcell::OrthogonalCell& cell,
                            double relative_tolerance,
                            Pcc0dParameters& parameters);

/// Numerical kernels require validated parameters and moments about a common
/// origin. Coordinates are relative to that origin. Energies/potentials are in
/// Hartree; gradients/forces are in Hartree/Bohr. The potential is for unit
/// positive charge; electron potential conversion belongs to the host adapter.
double pcc_0d_potential(const ChargeMoments& moments,
                        const ModuleBase::Vector3<double>& position,
                        const Pcc0dParameters& parameters);

ModuleBase::Vector3<double> pcc_0d_gradient(const ChargeMoments& moments,
                                          const ModuleBase::Vector3<double>& position,
                                          const Pcc0dParameters& parameters);

double pcc_0d_bilinear_energy(const ChargeMoments& left,
                             const ChargeMoments& right,
                             const Pcc0dParameters& parameters);

double pcc_0d_energy(const ChargeMoments& moments, const Pcc0dParameters& parameters);

ModuleBase::Vector3<double> pcc_0d_force(const ChargeMoments& moments,
                                       double ionic_charge,
                                       const ModuleBase::Vector3<double>& position,
                                       const Pcc0dParameters& parameters);

} // namespace elecstate

#endif
