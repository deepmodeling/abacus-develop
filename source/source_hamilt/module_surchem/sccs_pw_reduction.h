#ifndef SCCS_PW_REDUCTION_H
#define SCCS_PW_REDUCTION_H

#include "sccs_charge.h"
#include "sccs_poisson.h"

namespace ModuleSccs
{

class PoolPolarizationReduction : public PolarizationReduction
{
  public:
    explicit PoolPolarizationReduction(int process_count);

    void reduce_residual(double& square_sum,
                         double& maximum,
                         double& point_count) const override;
    void reduce_sum(double& value) const override;

  private:
    int process_count_ = 0;
};

class PoolChargeReduction : public ChargeReduction
{
  public:
    void reduce_sum(double& value) const override;
    void reduce_sum(double* values, int count) const override;
};

} // namespace ModuleSccs

#endif
