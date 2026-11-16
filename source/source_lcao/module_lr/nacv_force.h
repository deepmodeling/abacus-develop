#ifndef ABACUS_LR_NACV_FORCE_H
#define ABACUS_LR_NACV_FORCE_H
#include "lr_force.h"
namespace LR
{
template <typename T> struct NacvInputs;
/**
 * @brief Own geometry-fixed force intermediates for one NAC request batch.
 * @details The reference and parallel layouts are borrowed for the batch lifetime.
 * Potentials and derivative matrices are built lazily, including on empty MPI ranks.
 */
class NacvForce final : public LR_Force<double>
{
  public:
    /**
     * @param[in] data Geometry-local reference, density and integration dependencies.
     * @param[in] settings Resolved force switches, with no global configuration reads.
     */
    NacvForce(const NacvInputs<double> &data, const ForceSettings &settings);
    /** @return Immutable geometry intermediates, prepared on their first use. */
    const ForceGeometry &geometry();

  private:
    const NacvInputs<double> &data_;
    std::vector<hamilt::HContainer<double>> overlap_;
    std::vector<hamilt::HContainer<double>> kinetic_;
    std::unique_ptr<elecstate::Potential> local_;
    std::unique_ptr<elecstate::Potential> ground_;
    std::unique_ptr<ForceGeometry> geometry_;
    const ForceGeometry *geometry_terms() override;
    void prepare_geometry();
};
} // namespace LR
#endif
