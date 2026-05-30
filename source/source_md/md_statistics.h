#ifndef MD_STATISTICS_H
#define MD_STATISTICS_H

#include "source_base/matrix.h"

/**
 * @brief Kinetic energy and temperature statistics result - pure data structure
 *
 * Replaces the kinetic energy write-back via reference parameter
 * in current_temp(kinetic, natom, frozen_freedom, allmass, vel).
 */
struct MDKineticState
{
    double kinetic     = 0.0; ///< kinetic energy (Hartree)
    double temperature = 0.0; ///< temperature (Hartree)

    /// Convenience conversion to Kelvin
    double temperature_kelvin(double hartree_to_k) const
    {
        return temperature * hartree_to_k;
    }
};

/**
 * @brief Stress statistics result - pure data structure
 *
 * Replaces the implicit write-back of both virial and stress
 * through reference parameters in compute_stress().
 * Separates the ionic kinetic contribution tensor t_vector
 * from the total stress.
 */
struct MDStressState
{
    ModuleBase::matrix t_vector; ///< ionic kinetic contribution tensor (3x3)
    ModuleBase::matrix stress;   ///< total stress = virial + t_vector/omega (3x3)
};

/**
 * @brief FIRE optimizer projection statistics - pure data structure
 *
 * Replaces the scattered P, sumforce, normvel local variables
 * in FIRE::check_fire().
 */
struct FIREProjection
{
    double power         = 0.0; ///< P = sum v_i · f_i
    double force_norm    = 0.0; ///< |f| = sqrt(sum |f_i|^2)
    double velocity_norm = 0.0; ///< |v| = sqrt(sum |v_i|^2)
    double max_force     = 0.0; ///< max |f_i| component
};

#endif // MD_STATISTICS_H
