#ifndef SURCHEM_CHARGE_REDUCTION_H
#define SURCHEM_CHARGE_REDUCTION_H

namespace ModuleSurchem
{

// Sums local real-space integrals over the processes that share one grid.
class ChargeReduction
{
  public:
    virtual ~ChargeReduction() = default;

    virtual void reduce_sum(double& value) const = 0;
    virtual void reduce_sum(double* values, int count) const = 0;
};

// Single-process grid: only checks that the local values are finite.
class SerialChargeReduction : public ChargeReduction
{
  public:
    void reduce_sum(double& value) const override;
    void reduce_sum(double* values, int count) const override;
};

// Grid distributed over the processes of one pool.
class PoolChargeReduction : public ChargeReduction
{
  public:
    void reduce_sum(double& value) const override;
    void reduce_sum(double* values, int count) const override;
};

} // namespace ModuleSurchem

#endif
