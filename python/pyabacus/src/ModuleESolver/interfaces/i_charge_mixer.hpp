/**
 * @file i_charge_mixer.hpp
 * @brief Abstract interface for charge density mixer
 *
 * Defines the interface for mixing charge densities during SCF iterations.
 */

#ifndef PYABACUS_ESOLVER_I_CHARGE_MIXER_HPP
#define PYABACUS_ESOLVER_I_CHARGE_MIXER_HPP

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <string>

namespace py = pybind11;

namespace pyabacus {
namespace esolver {

/**
 * @brief Mixing method types
 */
enum class MixingMethod
{
    Plain,    ///< Simple linear mixing
    Pulay,    ///< Pulay mixing (DIIS)
    Broyden,  ///< Broyden mixing
    Anderson  ///< Anderson mixing
};

/**
 * @brief Configuration for charge mixing
 */
struct MixingConfig
{
    MixingMethod method = MixingMethod::Pulay;
    double beta = 0.7;           ///< Mixing parameter
    int ndim = 8;                ///< Number of history steps for Pulay/Broyden
    double gg0 = 1.0;            ///< Kerker mixing parameter
    bool mix_gg0 = false;        ///< Whether to use Kerker mixing
    bool mix_rho = true;         ///< Mix charge density (vs potential)
};

/**
 * @brief Abstract interface for charge density mixer
 *
 * This interface defines the contract for mixing charge densities
 * during SCF iterations to achieve convergence.
 */
class IChargeMixer
{
public:
    virtual ~IChargeMixer() = default;

    // ==================== Core Mixing Operations ====================

    /**
     * @brief Mix input and output charge densities
     * @param rho_in Input charge density (from previous iteration)
     * @param rho_out Output charge density (from current iteration)
     * @return Mixed charge density for next iteration
     */
    virtual py::array_t<double> mix(const py::array_t<double>& rho_in,
                                    const py::array_t<double>& rho_out) = 0;

    /**
     * @brief Reset mixer state (clear history)
     */
    virtual void reset() = 0;

    // ==================== State Queries ====================

    /**
     * @brief Get charge density difference from last mixing
     * @return drho value
     */
    virtual double get_drho() const = 0;

    /**
     * @brief Get current iteration count
     * @return Number of mixing iterations performed
     */
    virtual int get_iteration() const = 0;

    // ==================== Configuration ====================

    /**
     * @brief Set mixing configuration
     * @param config Mixing configuration
     */
    virtual void set_config(const MixingConfig& config) = 0;

    /**
     * @brief Get current mixing configuration
     * @return Current configuration
     */
    virtual MixingConfig get_config() const = 0;

    /**
     * @brief Set mixing parameter (beta)
     * @param beta Mixing parameter (0 < beta <= 1)
     */
    virtual void set_mixing_beta(double beta) = 0;

    /**
     * @brief Get mixing parameter
     * @return Current beta value
     */
    virtual double get_mixing_beta() const = 0;

    /**
     * @brief Set mixing method
     * @param method Mixing method to use
     */
    virtual void set_mixing_method(MixingMethod method) = 0;

    /**
     * @brief Get current mixing method
     * @return Current method
     */
    virtual MixingMethod get_mixing_method() const = 0;
};

/**
 * @brief Convert MixingMethod enum to string
 */
inline std::string mixing_method_to_string(MixingMethod method)
{
    switch (method)
    {
        case MixingMethod::Plain: return "plain";
        case MixingMethod::Pulay: return "pulay";
        case MixingMethod::Broyden: return "broyden";
        case MixingMethod::Anderson: return "anderson";
        default: return "unknown";
    }
}

/**
 * @brief Convert string to MixingMethod enum
 */
inline MixingMethod string_to_mixing_method(const std::string& str)
{
    if (str == "plain") return MixingMethod::Plain;
    if (str == "pulay") return MixingMethod::Pulay;
    if (str == "broyden") return MixingMethod::Broyden;
    if (str == "anderson") return MixingMethod::Anderson;
    return MixingMethod::Pulay; // default
}

} // namespace esolver
} // namespace pyabacus

#endif // PYABACUS_ESOLVER_I_CHARGE_MIXER_HPP
