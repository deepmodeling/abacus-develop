#ifndef HSOLVER_PRECISION_STRATEGY_H
#define HSOLVER_PRECISION_STRATEGY_H

/**
 * @file precision_strategy.h
 * @brief 精度选择策略 - 模板化的精度无关求解器包装
 *
 * 提供精度无关的求解器接口，支持运行时精度配置。
 * 通过策略模式分离精度选择逻辑和求解器实现。
 *
 * 使用方法:
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
 * @brief 精度选择策略基类
 *
 * @tparam SolverT 求解器类型 (如 DiagoCG, DiagoDavid)
 * @tparam T 数据类型 (double, complex<double> 等)
 * @tparam Device 设备类型
 */
template <template <typename, typename> class SolverT, typename T, typename Device = base_device::DEVICE_CPU>
class PrecisionStrategy
{
  public:
    using Real = typename GetTypeReal<T>::type;

    virtual ~PrecisionStrategy() = default;

    /**
     * @brief 获取当前精度模式
     */
    virtual PrecisionMode get_mode() const = 0;

    /**
     * @brief 获取精度模式对应的字符串描述
     */
    virtual std::string get_mode_string() const
    {
        return precision_mode_to_string(get_mode());
    }

    /**
     * @brief 检查是否适应当前问题规模
     *
     * 对于极小矩阵(dim < 50)，自动回退到双精度。
     *
     * @param dim 矩阵维度
     * @return 推荐的精度模式
     */
    static PrecisionMode recommend_mode(int dim)
    {
        if (dim < 50)
        {
            // 小矩阵：双精度更稳定，且性能差异不大
            return PrecisionMode::kDouble;
        }
        else if (dim < 200)
        {
            // 中等矩阵：混合精度平衡
            return PrecisionMode::kMixed;
        }
        else
        {
            // 大矩阵：混合精度收益明显
            return PrecisionMode::kMixed;
        }
    }

    /**
     * @brief 自动选择精度模式
     *
     * 根据矩阵维度和用户偏好自动选择最优精度模式。
     *
     * @param mode_str 用户指定的精度模式 ("auto", "double", "float", "mixed")
     * @param dim 矩阵维度
     * @return 最终选择的精度模式
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
 * @brief 双精度策略
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
 * @brief 混合精度策略
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
 * @brief 纯单精度策略 (用于快速原型和非关键计算)
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
 * @brief 精度策略工厂
 *
 * 根据 PrecisionMode 创建对应的策略对象。
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
     * @brief 从字符串创建策略
     */
    static std::unique_ptr<PrecisionStrategy<SolverT, T, Device>> create_from_string(const std::string& mode_str)
    {
        return create(parse_precision_mode(mode_str));
    }
};

} // namespace hsolver

#endif // HSOLVER_PRECISION_STRATEGY_H
