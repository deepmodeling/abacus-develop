#ifndef SNAP_CALC_H
#define SNAP_CALC_H

#include "desc.h"
#include "model.h"

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

    ForceResult eval_force(const std::vector<Neighbor>& neighbors,
                                             std::size_t central_element) const;

    ForceResult eval_force_fd(const std::vector<Neighbor>& neighbors,
                                                     std::size_t central_element,
                                                     double displacement) const;

    std::size_t descriptor_size() const
    {
        return desc_size_;
    }

  private:
    ModelData model_;
    DescEval desc_eval_;
    std::size_t desc_size_;
};

class CellAccumulator
{
  public:
    explicit CellAccumulator(std::size_t atom_count);

    void add(std::size_t center_index,
             const std::vector<std::size_t>& neigh_indices,
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

#endif // SNAP_CALC_H
