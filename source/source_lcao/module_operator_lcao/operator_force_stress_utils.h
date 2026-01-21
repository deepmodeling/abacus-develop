#ifndef OPERATOR_FORCE_STRESS_UTILS_H
#define OPERATOR_FORCE_STRESS_UTILS_H

#include "source_base/matrix.h"
#include "source_cell/unitcell.h"
#include <vector>

namespace OperatorForceStress {

/**
 * @brief Setup step_trace array for handling npol (spin polarization)
 *
 * For npol=1 (non-spin-polarized): step_trace = {0}
 * For npol=2 (spin-polarized): step_trace = {0, 1, col_size, col_size+1}
 *
 * @param npol Number of spin polarizations (1 or 2)
 * @param col_size Number of columns in the density matrix block
 * @param step_trace Output vector to store step trace values
 */
inline void setup_step_trace(int npol, int col_size, std::vector<int>& step_trace)
{
    step_trace.resize(npol * npol, 0);
    if (npol == 2)
    {
        step_trace[1] = 1;
        step_trace[2] = col_size;
        step_trace[3] = col_size + 1;
    }
}

/**
 * @brief Structure to hold orbital quantum numbers
 */
struct OrbitalQuantumNumbers
{
    int L;  ///< Angular momentum quantum number
    int N;  ///< Principal quantum number
    int m;  ///< Magnetic quantum number (internal indexing)
    int M;  ///< Magnetic quantum number (standard convention)
};

/**
 * @brief Extract orbital quantum numbers from atom and orbital index
 *
 * @param atom Atom object containing orbital information
 * @param iw Orbital index within the atom
 * @return OrbitalQuantumNumbers structure with L, N, m, M values
 */
inline OrbitalQuantumNumbers get_orbital_qn(const Atom& atom, int iw)
{
    OrbitalQuantumNumbers qn;
    qn.L = atom.iw2l[iw];
    qn.N = atom.iw2n[iw];
    qn.m = atom.iw2m[iw];
    qn.M = (qn.m % 2 == 0) ? -qn.m / 2 : (qn.m + 1) / 2;
    return qn;
}

/**
 * @brief Helper function to extract real part from complex or real values
 *
 * Template specialization handles both std::complex<T> and double types
 */
template <typename T>
inline double get_real_part(const T& val)
{
    return val.real();
}

template <>
inline double get_real_part<double>(const double& val)
{
    return val;
}

/**
 * @brief Rearrange stress from 6-component vector to 3x3 matrix format
 *
 * Input format: [xx, xy, xz, yy, yz, zz]
 * Output format: 3x3 matrix with proper indexing
 *
 * @param stress Matrix to rearrange (must be at least 3x3)
 */
inline void rearrange_stress_matrix(ModuleBase::matrix& stress)
{
    stress.c[8] = stress.c[5]; // stress(2,2)
    stress.c[7] = stress.c[4]; // stress(2,1)
    stress.c[6] = stress.c[2]; // stress(2,0)
    stress.c[5] = stress.c[4]; // stress(1,2)
    stress.c[4] = stress.c[3]; // stress(1,1)
    stress.c[3] = stress.c[1]; // stress(1,0)
}

/**
 * @brief Finalize force and stress calculations with MPI reduction and post-processing
 *
 * Performs:
 * 1. MPI reduction of force and stress across all processes
 * 2. Apply force_factor (typically 2.0 for Hermitian matrices)
 * 3. Apply stress weight (lat0/omega) and stress_factor
 * 4. Rearrange stress matrix to 3x3 format
 *
 * @param cal_force Whether force calculation is enabled
 * @param cal_stress Whether stress calculation is enabled
 * @param ucell Unit cell containing lattice parameters
 * @param stress_tmp Temporary 6-component stress vector
 * @param force Force matrix to finalize
 * @param stress Stress matrix to finalize
 * @param force_factor Multiplicative factor for force (default: 2.0 for Hermitian)
 * @param stress_factor Multiplicative factor for stress (default: 2.0 for Hermitian)
 */
void finalize_force_stress(
    bool cal_force,
    bool cal_stress,
    const UnitCell* ucell,
    const std::vector<double>& stress_tmp,
    ModuleBase::matrix& force,
    ModuleBase::matrix& stress,
    double force_factor = 2.0,
    double stress_factor = 2.0);

} // namespace OperatorForceStress

#endif // OPERATOR_FORCE_STRESS_UTILS_H
