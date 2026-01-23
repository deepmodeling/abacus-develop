#ifndef OPERATOR_FORCE_STRESS_UTILS_HPP
#define OPERATOR_FORCE_STRESS_UTILS_HPP

#include "operator_force_stress_utils.h"
#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "source_lcao/module_hcontainer/hcontainer.h"
#include "source_basis/module_ao/parallel_orbitals.h"

namespace OperatorForceStress {

/**
 * @brief Template function for calculating force and stress from 2-center integrals
 *
 * This template unifies the force/stress calculation pattern for operators that use
 * 2-center integrals (e.g., overlap, kinetic energy). The sign conventions for force
 * and stress are controlled by template parameters to avoid runtime overhead.
 *
 * @tparam TK Type for k-space matrices (double or std::complex<double>)
 * @tparam TR Type for real-space matrices (typically double)
 * @tparam IntegralFunc Functor type for calculating integrals
 * @tparam ForceSign Sign convention for force (+1 or -1)
 * @tparam StressSign Sign convention for stress (+1 or -1)
 *
 * @param cal_force Whether to calculate forces
 * @param cal_stress Whether to calculate stress
 * @param dmR Density matrix in real space
 * @param ucell Unit cell containing atomic structure
 * @param gridD Grid driver for finding adjacent atoms
 * @param orb_cutoff Orbital cutoff radii for each atom type
 * @param paraV Parallel orbital distribution information
 * @param integral_calculator Functor that calculates integral and its derivatives
 * @param force Output force matrix (natom x 3)
 * @param stress Output stress matrix (3 x 3)
 */
template <typename TK, typename TR, typename IntegralFunc, int ForceSign, int StressSign>
void cal_force_stress_2center(
    const bool cal_force,
    const bool cal_stress,
    const hamilt::HContainer<double>* dmR,
    const UnitCell* ucell,
    const Grid_Driver* gridD,
    const std::vector<double>& orb_cutoff,
    const Parallel_Orbitals* paraV,
    IntegralFunc& integral_calculator,
    ModuleBase::matrix& force,
    ModuleBase::matrix& stress)
{
    const int npol = ucell->get_npol();
    std::vector<double> stress_tmp(6, 0);

    if (cal_force)
    {
        force.zero_out();
    }

    // Loop over all atom pairs and calculate force/stress contributions
    #pragma omp parallel
    {
        std::vector<double> stress_local(6, 0);
        ModuleBase::matrix force_local(force.nr, force.nc);

        #pragma omp for schedule(dynamic)
        for (int iat1 = 0; iat1 < ucell->nat; iat1++)
        {
            auto tau1 = ucell->get_tau(iat1);
            int T1 = 0, I1 = 0;
            ucell->iat2iait(iat1, &I1, &T1);
            Atom& atom1 = ucell->atoms[T1];

            // Find adjacent atoms
            AdjacentAtomInfo adjs;
            gridD->Find_atom(*ucell, tau1, T1, I1, &adjs);

            for (int ad = 0; ad < adjs.adj_num + 1; ++ad)
            {
                const int T2 = adjs.ntype[ad];
                const int I2 = adjs.natom[ad];
                const int iat2 = ucell->itia2iat(T2, I2);
                const ModuleBase::Vector3<int>& R_index = adjs.box[ad];

                // Check cutoff
                ModuleBase::Vector3<double> dtau = ucell->cal_dtau(iat1, iat2, R_index);
                if (dtau.norm() * ucell->lat0 >= orb_cutoff[T1] + orb_cutoff[T2])
                {
                    continue;
                }

                // Find density matrix for this atom pair
                const hamilt::BaseMatrix<double>* dm_matrix = dmR->find_matrix(iat1, iat2, R_index[0], R_index[1], R_index[2]);
                if (dm_matrix == nullptr)
                {
                    continue;
                }

                // Calculate force and stress for this atom pair
                double* force_tmp1 = (cal_force) ? &force_local(iat1, 0) : nullptr;
                double* force_tmp2 = (cal_force) ? &force_local(iat2, 0) : nullptr;

                Atom& atom2 = ucell->atoms[T2];
                auto row_indexes = paraV->get_indexes_row(iat1);
                auto col_indexes = paraV->get_indexes_col(iat2);

                if (row_indexes.size() == 0 || col_indexes.size() == 0)
                {
                    continue;
                }

                const double* dm_pointer = dm_matrix->get_pointer();
                double olm[4] = {0, 0, 0, 0}; // value, dx, dy, dz

                // step_trace = 0 for npol=1; ={0, 1, col_size, col_size+1} for npol=2
                std::vector<int> step_trace(npol * npol, 0);
                if (npol == 2)
                {
                    step_trace[1] = 1;
                    step_trace[2] = col_indexes.size();
                    step_trace[3] = col_indexes.size() + 1;
                }

                // Loop over orbital pairs
                for (int iw1l = 0; iw1l < row_indexes.size(); iw1l += npol)
                {
                    const int iw1 = row_indexes[iw1l] / npol;
                    const int L1 = atom1.iw2l[iw1];
                    const int N1 = atom1.iw2n[iw1];
                    const int m1 = atom1.iw2m[iw1];
                    const int M1 = (m1 % 2 == 0) ? -m1 / 2 : (m1 + 1) / 2;

                    for (int iw2l = 0; iw2l < col_indexes.size(); iw2l += npol)
                    {
                        const int iw2 = col_indexes[iw2l] / npol;
                        const int L2 = atom2.iw2l[iw2];
                        const int N2 = atom2.iw2n[iw2];
                        const int m2 = atom2.iw2m[iw2];
                        const int M2 = (m2 % 2 == 0) ? -m2 / 2 : (m2 + 1) / 2;

                        // Calculate integral and its gradient using provided functor
                        integral_calculator(T1, L1, N1, M1, T2, L2, N2, M2, dtau, olm);

                        // only charge should be considered
                        double dm_current = get_real_part(dm_pointer[0]);

                        // Calculate force contribution with compile-time sign
                        if (cal_force)
                        {
                            // Factor of 2 for Hermitian matrix will be applied later
                            for (int i = 0; i < 3; i++)
                            {
                                force_tmp1[i] += ForceSign * dm_current * olm[i + 1];
                                force_tmp2[i] -= ForceSign * dm_current * olm[i + 1];
                            }
                        }

                        // Calculate stress contribution with compile-time sign
                        if (cal_stress)
                        {
                            stress_local[0] += StressSign * dm_current * olm[1] * dtau.x; // xx
                            stress_local[1] += StressSign * dm_current * olm[1] * dtau.y; // xy
                            stress_local[2] += StressSign * dm_current * olm[1] * dtau.z; // xz
                            stress_local[3] += StressSign * dm_current * olm[2] * dtau.y; // yy
                            stress_local[4] += StressSign * dm_current * olm[2] * dtau.z; // yz
                            stress_local[5] += StressSign * dm_current * olm[3] * dtau.z; // zz
                        }

                        dm_pointer += npol;
                    }
                    dm_pointer += (npol - 1) * col_indexes.size();
                }
            }
        }

        #pragma omp critical
        {
            if (cal_force)
            {
                force += force_local;
            }
            if (cal_stress)
            {
                for (int i = 0; i < 6; i++)
                {
                    stress_tmp[i] += stress_local[i];
                }
            }
        }
    }

    // Finalize with MPI reduction and post-processing
    finalize_force_stress(cal_force, cal_stress, ucell, stress_tmp, force, stress, 1.0, 1.0);
}

} // namespace OperatorForceStress

#endif // OPERATOR_FORCE_STRESS_UTILS_HPP
