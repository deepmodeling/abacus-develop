#ifndef TD_OBS_PW_H
#define TD_OBS_PW_H

#include <array>
#include <vector>

namespace pw
{
/** @brief Independent requests; all ranks must use the same choices. */
struct CurrentRequest
{
    bool charge_total = false;
    bool charge_per_k = false;
    bool spin_current = false;
};

/** @brief Globally reduced currents; unrequested arrays stay zero or empty. */
struct TDVelResultPW
{
    std::array<double, 3> total{}; ///< Electron charge current per volume, in Hartree atomic units.
    std::vector<double> per_k;     ///< Same units, [global k including spin][Cartesian direction].
    std::array<double, 9> spin{};  ///< Spin angular momentum flux, [flow][spin], hbar/(Bohr^2*atomic time).
};

/** @brief Endpoint observables replicated on all ranks after the appropriate reductions. */
struct TDSpinResultPW
{
    std::array<double, 3> mag{};  ///< ABACUS Pauli moment in mu_B; no electron minus sign.
    double mag_abs = 0.0;         ///< Integral of the magnetization magnitude, in mu_B.
    double electrons = 0.0;       ///< Valence electron count; excludes the NLCC core density.
    std::vector<double> atom_mag; ///< [atom][Cartesian direction], in mu_B.
    std::array<double, 3> xc{};   ///< Rate of spin angular momentum S=M/2, hbar/atomic time.
    std::array<double, 3> soc{};  ///< Same units as xc.
    std::vector<double> atom_soc; ///< Pseudopotential atom decomposition, [atom][direction].
};
} // namespace pw
#endif
