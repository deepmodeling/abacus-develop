#ifndef ABACUS_LR_NACV_XC_NAMES_H
#define ABACUS_LR_NACV_XC_NAMES_H

#include <algorithm>
#include <cctype>
#include <set>
#include <string>

namespace LR
{
namespace NacvXC
{
/**
 * @brief Resolve built-in ground-state aliases before selecting a response kernel.
 * @param[in] name Case-insensitive functional name from INPUT or a pseudopotential.
 * @return Lowercase canonical name; unknown names are only converted to lowercase.
 * @details This maps names, not physical approximations: NONE does not become RPA.
 */
inline std::string canonical_xc_name(const std::string &name)
{
    std::string canonical = name;
    std::transform(
        canonical.begin(),
        canonical.end(),
        canonical.begin(),
        [](unsigned char character)
        {
            return static_cast<char>(std::tolower(character));
        });
    if (canonical == "pz" || canonical == "slapznogxnogc")
    {
        return "lda";
    }
    if (canonical == "slapwpbxpbc")
    {
        return "pbe";
    }
    if (canonical == "nox+noc")
    {
        return "none";
    }
    return canonical;
}

/**
 * @brief Enumerate canonical LR kernels carrying exact exchange.
 * @return Immutable lowercase names; weights/range separation come from Coulomb settings.
 * @details Membership does not select the Coulomb operator. The general kernel
 * [alpha + beta erfc(mu r)]/r comes from coulomb_param, assembled from the KS
 * functional and the explicit EXX parameters. Its overall scale is supplied by
 * hybrid_alpha = max(abs(alpha), abs(beta)); this list must not replace those
 * settings with a functional-name-dependent prefactor.
 */
inline const std::set<std::string> &hybrid_xc_list()
{
    static const std::set<std::string> names = {"hf",
        "hse",
        "pbe0",
        "b3lyp",
        "cam_pbeh",
        "lc_pbe",
        "lc_wpbe",
        "lrc_wpbe",
        "lrc_wpbeh"};
    return names;
}

/**
 * @brief Identify implemented LR kernels with a semilocal XC remainder.
 * @param[in] name Canonical lowercase kernel name, not an arbitrary KS name.
 * @return True for implemented LDA/GGA kernels and hybrids other than pure HF.
 * @details LibXC supplies each hybrid's weighted semilocal remainder after its
 * external parameters are set. Membership does not certify analytic NAC accuracy.
 */
inline bool has_local_xc(const std::string &name)
{
    if (name == "lda" || name == "pwlda" || name == "pbe")
    {
        return true;
    }
    return name != "hf" && hybrid_xc_list().count(name);
}
} // namespace NacvXC
} // namespace LR

#endif
