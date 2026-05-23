#ifndef HSOLVER_PRECISION_MODE_H
#define HSOLVER_PRECISION_MODE_H

#include <string>

namespace hsolver
{

/**
 * @brief Precision mode for diagonalization solvers.
 *
 * Controls the numerical precision used in iterative eigensolvers:
 * - kDouble: Pure double precision (default, highest accuracy)
 * - kFloat:  Pure single precision (fastest, for non-critical calculations)
 * - kMixed:  Mixed precision (Float iteration + Double refinement, recommended)
 */
enum class PrecisionMode
{
    kDouble = 0,  ///< Pure double precision
    kFloat  = 1,  ///< Pure single precision
    kMixed  = 2   ///< Mixed precision (float iteration + double refinement)
};

} // namespace hsolver

/**
 * @brief Parse precision mode from string.
 * @param mode_str "double", "float", "mixed", "single", or "auto"
 * @return Corresponding PrecisionMode enum value.
 */
inline hsolver::PrecisionMode parse_precision_mode(const std::string& mode_str)
{
    if (mode_str == "float" || mode_str == "single")
    {
        return hsolver::PrecisionMode::kFloat;
    }
    else if (mode_str == "mixed" || mode_str == "auto")
    {
        return hsolver::PrecisionMode::kMixed;
    }
    else
    {
        return hsolver::PrecisionMode::kDouble;
    }
}

/**
 * @brief Convert precision mode to string representation.
 */
inline std::string precision_mode_to_string(hsolver::PrecisionMode mode)
{
    switch (mode)
    {
        case hsolver::PrecisionMode::kFloat:  return "float";
        case hsolver::PrecisionMode::kMixed:  return "mixed";
        case hsolver::PrecisionMode::kDouble:
        default:                               return "double";
    }
}

#endif // HSOLVER_PRECISION_MODE_H
