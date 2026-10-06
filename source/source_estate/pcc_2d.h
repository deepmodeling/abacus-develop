#ifndef PCC_2D_H
#define PCC_2D_H

#include "charge_moments.h"

namespace unitcell
{
struct SlabCell;
}

namespace elecstate
{
struct Pcc2dParameters
{
    double area = 0.0; ///< Periodic area in Bohr^2
    double length = 0.0; ///< Normal period in Bohr
};

bool make_pcc_2d_parameters(const unitcell::SlabCell& cell,
                            Pcc2dParameters& parameters,
                            std::string& error);

/// Hartree kernels. Parameters must be validated; moments use normal-projected
/// positions r = normal * u, so second_moment is the normal second moment.
/// The open planar kernel is zero on the source plane; its periodic counterpart
/// has zero mean. This gauge fixes the charged-slab constant to -pi*L/(3*A);
/// see M. J. Rutter, Electron. Struct. 3, 015002 (2021).
double pcc_2d_potential(const ChargeMoments& moments,
                        double coordinate,
                        const ModuleBase::Vector3<double>& normal,
                        const Pcc2dParameters& parameters);

ModuleBase::Vector3<double> pcc_2d_gradient(const ChargeMoments& moments,
                                           double coordinate,
                                           const ModuleBase::Vector3<double>& normal,
                                           const Pcc2dParameters& parameters);

double pcc_2d_energy(const ChargeMoments& moments, const Pcc2dParameters& parameters);

ModuleBase::Vector3<double> pcc_2d_force(const ChargeMoments& moments,
                                        double ionic_charge,
                                        double coordinate,
                                        const ModuleBase::Vector3<double>& normal,
                                        const Pcc2dParameters& parameters);
} // namespace elecstate

#endif
