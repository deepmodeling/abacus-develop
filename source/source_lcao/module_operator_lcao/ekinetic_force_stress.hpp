#pragma once
#include "ekinetic_new.h"
#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"

namespace hamilt
{

// Helper function to get real part
template <typename T>
inline double get_real(const T& val)
{
    return val.real();
}

template <>
inline double get_real<double>(const double& val)
{
    return val;
}

template <typename TK, typename TR>
void EkineticNew<OperatorLCAO<TK, TR>>::cal_force_stress(const bool cal_force,
                                                          const bool cal_stress,
                                                          const HContainer<double>* dmR,
                                                          ModuleBase::matrix& force,
                                                          ModuleBase::matrix& stress)
{
    ModuleBase::TITLE("EkineticNew", "cal_force_stress");
    ModuleBase::timer::tick("EkineticNew", "cal_force_stress");

    const Parallel_Orbitals* paraV = dmR->get_paraV();
    const int npol = this->ucell->get_npol();
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
        for (int iat1 = 0; iat1 < this->ucell->nat; iat1++)
        {
            auto tau1 = ucell->get_tau(iat1);
            int T1 = 0, I1 = 0;
            ucell->iat2iait(iat1, &I1, &T1);
            Atom& atom1 = this->ucell->atoms[T1];

            // Find adjacent atoms
            AdjacentAtomInfo adjs;
            this->gridD->Find_atom(*ucell, tau1, T1, I1, &adjs);

            for (int ad = 0; ad < adjs.adj_num + 1; ++ad)
            {
                const int T2 = adjs.ntype[ad];
                const int I2 = adjs.natom[ad];
                const int iat2 = ucell->itia2iat(T2, I2);
                const ModuleBase::Vector3<int>& R_index = adjs.box[ad];

                // Check cutoff
                ModuleBase::Vector3<double> dtau = this->ucell->cal_dtau(iat1, iat2, R_index);
                if (dtau.norm() * this->ucell->lat0 >= orb_cutoff_[T1] + orb_cutoff_[T2])
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

                Atom& atom2 = this->ucell->atoms[T2];
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
                if (npol == 2) {
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

                        // Calculate kinetic integral and its gradient
                        intor_->calculate(T1, L1, N1, M1, T2, L2, N2, M2,
                                        dtau * this->ucell->lat0,
                                        &olm[0], &olm[1]);

                        // only charge should be considered
                        double dm_current = get_real(dm_pointer[0]);

                        // Calculate force contribution
                        if (cal_force)
                        {
                            // F = -sum(dm * dT/dr)
                            // Factor of 2 for Hermitian matrix will be applied later
                            for (int i = 0; i < 3; i++)
                            {
                                force_tmp1[i] += dm_current * olm[i + 1];
                                force_tmp2[i] -= dm_current * olm[i + 1];
                            }
                        }

                        // Calculate stress contribution
                        if (cal_stress)
                        {
                            // stress_ij = sum(dm * dT/dr_i * r_j)
                            stress_local[0] -= dm_current * olm[1] * dtau.x; // xx
                            stress_local[1] -= dm_current * olm[1] * dtau.y; // xy
                            stress_local[2] -= dm_current * olm[1] * dtau.z; // xz
                            stress_local[3] -= dm_current * olm[2] * dtau.y; // yy
                            stress_local[4] -= dm_current * olm[2] * dtau.z; // yz
                            stress_local[5] -= dm_current * olm[3] * dtau.z; // zz
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

    if (cal_force)
    {
#ifdef __MPI
        Parallel_Reduce::reduce_all(force.c, force.nr * force.nc);
#endif
    }

    if (cal_stress)
    {
#ifdef __MPI
        Parallel_Reduce::reduce_all(stress_tmp.data(), 6);
#endif
        const double weight = this->ucell->lat0 / this->ucell->omega;
        for (int i = 0; i < 6; i++)
        {
            stress.c[i] = stress_tmp[i] * weight;
        }
        // Rearrange to 3x3 matrix format
        stress.c[8] = stress.c[5]; // stress(2,2)
        stress.c[7] = stress.c[4]; // stress(2,1)
        stress.c[6] = stress.c[2]; // stress(2,0)
        stress.c[5] = stress.c[4]; // stress(1,2)
        stress.c[4] = stress.c[3]; // stress(1,1)
        stress.c[3] = stress.c[1]; // stress(1,0)
    }

    ModuleBase::timer::tick("EkineticNew", "cal_force_stress");
}

// Dummy implementations for cal_force_IJR and cal_stress_IJR
// These are not used in the simplified approach above
template <typename TK, typename TR>
void EkineticNew<OperatorLCAO<TK, TR>>::cal_force_IJR(
    const int& iat1,
    const int& iat2,
    const Parallel_Orbitals* paraV,
    const std::unordered_map<int, std::vector<double>>& nlm1_all,
    const std::unordered_map<int, std::vector<double>>& nlm2_all,
    const hamilt::BaseMatrix<TR>* dmR_pointer,
    double* force1,
    double* force2)
{
    // Not used in current implementation
}

template <typename TK, typename TR>
void EkineticNew<OperatorLCAO<TK, TR>>::cal_stress_IJR(
    const int& iat1,
    const int& iat2,
    const Parallel_Orbitals* paraV,
    const std::unordered_map<int, std::vector<double>>& nlm1_all,
    const std::unordered_map<int, std::vector<double>>& nlm2_all,
    const hamilt::BaseMatrix<TR>* dmR_pointer,
    const ModuleBase::Vector3<double>& dis1,
    const ModuleBase::Vector3<double>& dis2,
    double* stress)
{
    // Not used in current implementation
}

} // namespace hamilt
