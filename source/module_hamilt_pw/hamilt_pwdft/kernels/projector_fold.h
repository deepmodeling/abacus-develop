#ifndef HAMILT_PROJECTOR_FOLD_H
#define HAMILT_PROJECTOR_FOLD_H

#include "module_base/vector3.h"
#include <complex>
#include <memory>

namespace ModulePW { class PW_Basis_K; }

namespace hamilt
{

/// Factor CUDA complex-double projector products using atomic phase periodicity.
/// Geometry and the reciprocal basis remain fixed for the object's lifetime.
class ProjectorFold
{
  public:
    /// Return a plan when a bounded binary grid covers enough identical atoms.
    /// Positions are host arrays of length atoms; no coordinate is rounded.
    static std::unique_ptr<ProjectorFold> create(int atoms, int channels,
        const ModuleBase::Vector3<double>* direct_positions,
        const ModuleBase::Vector3<double>* cartesian_positions,
        const ModulePW::PW_Basis_K& basis);

    ~ProjectorFold();

    /// Prepare the current k point from the complete device projector table,
    /// with shape (projector_leading, atoms * channels), column-major.
    void prepare(int ik, const std::complex<double>* projectors, int projector_leading);

    /// Form P^H input for spin-component columns. Device input/output matrices
    /// have leading dimensions leading and atoms * channels, respectively.
    /// False leaves output untouched and requests the ordinary product.
    bool project(int columns, int leading, const std::complex<double>* input,
                 std::complex<double>* output);

    /// Add P coefficients^T to output. Device coefficients have shape
    /// (columns, atoms * channels), column-major; output has leading dimension
    /// leading. Padding is unchanged. False leaves output untouched.
    bool add(int columns, int leading, const std::complex<double>* coefficients,
             std::complex<double>* output);

  private:
    struct Impl;
    explicit ProjectorFold(std::unique_ptr<Impl> implementation);
    std::unique_ptr<Impl> implementation_;
};

}
#endif
