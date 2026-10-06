#ifndef CHARGE_MOMENTS_H
#define CHARGE_MOMENTS_H

#include "source_base/vector3.h"

namespace elecstate
{

/// Moments about one common origin. Charges carry their physical sign.
struct ChargeMoments
{
    double charge = 0.0;
    ModuleBase::Vector3<double> dipole;
    double second_moment = 0.0; ///< Trace of the second moment, in e Bohr^2
};

/// Local moments only; no MPI or cell/model dependencies. Use weight=1 for
/// point charges and the volume element for a signed grid density.
/// Zero samples are permitted.
ChargeMoments charge_moments(const double* charges,
                             const ModuleBase::Vector3<double>* relative_positions,
                             int count,
                             double weight);

ChargeMoments add_charge_moments(const ChargeMoments& left, const ChargeMoments& right);

} // namespace elecstate

#endif
