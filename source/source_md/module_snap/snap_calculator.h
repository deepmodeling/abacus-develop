#ifndef SNAP_CALCULATOR_H
#define SNAP_CALCULATOR_H

#include "snap_descriptor.h"
#include "snap_model.h"

#include <array>
#include <vector>

namespace ModuleMD
{
namespace SNAP
{

struct EnergyResult
{
    double energy = 0.0;
    std::vector<double> descriptor;
};

struct ForceResult
{
    EnergyResult energy;
    std::vector<double> center_force;
    std::vector<std::vector<double> > neighbor_forces;
};

class Calculator
{
  public:
    explicit Calculator(const ModelData& model);

    EnergyResult compute(const std::vector<Neighbor>& neighbors, std::size_t central_element) const;

    ForceResult compute_with_analytic_force(const std::vector<Neighbor>& neighbors,
                                             std::size_t central_element) const;

    ForceResult compute_with_finite_difference_force(const std::vector<Neighbor>& neighbors,
                                                     std::size_t central_element,
                                                     double displacement) const;

  private:
    ModelData model_;
    DescriptorEvaluator descriptor_evaluator_;
    std::size_t descriptor_size_;
};

class CellAccumulator
{
  public:
    explicit CellAccumulator(std::size_t atom_count);

    void add(std::size_t center_index,
             const std::vector<std::size_t>& neighbor_indices,
             const std::vector<Neighbor>& neighbors,
             const ForceResult& result);

    double energy() const;
    const std::vector<std::array<double, 3> >& forces() const;
    const std::array<double, 9>& virial() const;

  private:
    double energy_ = 0.0;
    std::vector<std::array<double, 3> > forces_;
    std::array<double, 9> virial_;
};

} // namespace SNAP
} // namespace ModuleMD

#endif // SNAP_CALCULATOR_H
