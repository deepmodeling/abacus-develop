#ifndef HSOLVER_PRECISION_STRATEGY_H
#define HSOLVER_PRECISION_STRATEGY_H

/**
 * @file precision_strategy.h
 * @brief Precision selection strategy - template-based precision-agnostic solver wrapper
 *
 * Provides precision-agnostic solver interfaces with runtime precision configuration.
 * Separates precision selection logic from solver implementation via the strategy pattern.
 *
 * Usage:
 *   auto solver = make_precision_solver<DiagoCG>(PrecisionMode::kMixed, ...);
 *   solver.diag(...);
 */

#include "source_hsolver/precision_mode.h"
#include "source_hsolver/diago_cg.h"
#include <memory>
#include <stdexcept>
#include <string>

namespace hsolver
{

/**
 * @brief Base class for precision selection strategy
 *
 * @tparam SolverT Solver type (e.g., DiagoCG, DiagoDavid)
 * @tparam T Data type (e.g., double, complex<double>)
 * @tparam Device Device type
 */
template <template <typename, typename> class SolverT, typename T, typename Device = base_device::DEVICE_CPU>
class PrecisionStrategy
{
  public:
    using Real = typename GetTypeReal<T>::type;

    virtual ~PrecisionStrategy() = default;

    /**
     * @brief Get the current precision mode
     */
    virtual PrecisionMode get_mode() const = 0;

    /**
     * @brief Get string description of the current precision mode
     */
    virtual std::string get_mode_string() const
    {
        return precision_mode_to_string(get_mode());
    }

    /**
     * @brief Check whether the strategy is suitable for the given problem size
     *
     * For very small matrices (dim < 50), falls back to double precision.
     *
     * @param dim Matrix dimension
     * @return Recommended precision mode
     */
    static PrecisionMode recommend_mode(int dim)
    {
        if (dim < 50)
        {
            // Small matrix: double precision is more stable, performance difference is negligible
            return PrecisionMode::kDouble;
        }
        else if (dim < 200)
        {
            // Medium matrix: balanced mixed precision
            return PrecisionMode::kMixed;
        }
        else
        {
            // Large matrix: mixed precision provides clear benefit
            return PrecisionMode::kMixed;
        }
    }

    /**
     * @brief Auto-select the optimal precision mode
     *
     * Selects the best precision mode based on matrix dimension and user preference.
     *
     * @param mode_str User-specified precision mode ("auto", "double", "float", "mixed")
     * @param dim Matrix dimension
     * @return Final selected precision mode
     */
    static PrecisionMode auto_select_mode(const std::string& mode_str, int dim)
    {
        if (mode_str == "auto" || mode_str.empty())
        {
            return recommend_mode(dim);
        }
        return parse_precision_mode(mode_str);
    }
};

/**
 * @brief Double precision strategy
 */
template <template <typename, typename> class SolverT, typename T, typename Device = base_device::DEVICE_CPU>
class DoublePrecisionStrategy : public PrecisionStrategy<SolverT, T, Device>
{
  public:
    PrecisionMode get_mode() const override
    {
        return PrecisionMode::kDouble;
    }
};

/**
 * @brief Mixed precision strategy
 */
template <template <typename, typename> class SolverT, typename T, typename Device = base_device::DEVICE_CPU>
class MixedPrecisionStrategy : public PrecisionStrategy<SolverT, T, Device>
{
  public:
    PrecisionMode get_mode() const override
    {
        return PrecisionMode::kMixed;
    }
};

/**
 * @brief Float precision strategy (for fast prototyping and non-critical calculations)
 */
template <template <typename, typename> class SolverT, typename T, typename Device = base_device::DEVICE_CPU>
class FloatPrecisionStrategy : public PrecisionStrategy<SolverT, T, Device>
{
  public:
    PrecisionMode get_mode() const override
    {
        return PrecisionMode::kFloat;
    }
};

/**
 * @brief Precision strategy factory
 *
 * Creates the corresponding strategy object based on PrecisionMode.
 */
template <template <typename, typename> class SolverT, typename T, typename Device = base_device::DEVICE_CPU>
class PrecisionStrategyFactory
{
  public:
    static std::unique_ptr<PrecisionStrategy<SolverT, T, Device>> create(PrecisionMode mode)
    {
        switch (mode)
        {
            case PrecisionMode::kFloat:
                return std::make_unique<FloatPrecisionStrategy<SolverT, T, Device>>();
            case PrecisionMode::kMixed:
                return std::make_unique<MixedPrecisionStrategy<SolverT, T, Device>>();
            case PrecisionMode::kDouble:
            default:
                return std::make_unique<DoublePrecisionStrategy<SolverT, T, Device>>();
        }
    }

    /**
     * @brief Create strategy from string
     */
    static std::unique_ptr<PrecisionStrategy<SolverT, T, Device>> create_from_string(const std::string& mode_str)
    {
        return create(parse_precision_mode(mode_str));
    }
};

} // namespace hsolver

#endif // HSOLVER_PRECISION_STRATEGY_H
