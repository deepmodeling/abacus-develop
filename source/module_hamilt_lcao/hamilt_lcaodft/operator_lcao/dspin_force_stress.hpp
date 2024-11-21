#pragma once
#include "dspin_lcao.h"
#include "module_base/parallel_reduce.h"
#include "module_base/timer.h"

namespace hamilt
{

template <typename TK, typename TR>
void DeltaSpin<OperatorLCAO<TK, TR>>::cal_force_stress(const bool cal_force,
                                                  const bool cal_stress,
                                                  const HContainer<double>* dmR,
                                                  ModuleBase::matrix& force,
                                                  ModuleBase::matrix& stress)
{
    ModuleBase::TITLE("DeltaSpin", "cal_force_stress");

    // begin the calculation of force and stress
    ModuleBase::timer::tick("DeltaSpin", "cal_force_stress");

    spinconstrain::SpinConstrain<TK>& sc = spinconstrain::SpinConstrain<TK>::getScInstance();
    auto& constrain = sc.get_constrain();
    this->cal_constraint_atom_list(constrain);
    auto& lambda = sc.get_sc_lambda();

    const Parallel_Orbitals* paraV = dmR->get_paraV();
    const int npol = this->ucell->get_npol();
    std::vector<double> stress_tmp(6, 0);
    if (cal_force)
    {
        force.zero_out();
    }
    // 1. calculate <psi|beta> for each pair of atoms
    // loop over all on-site atoms
    int atom_index = 0;
    for (int iat0 = 0; iat0 < this->ucell->nat; iat0++)
    {
        if(!this->constraint_atom_list[iat0])
        {
            continue;
        }
      
        // skip the atoms without plus-U
        auto tau0 = ucell->get_tau(iat0);
        int T0, I0;
        ucell->iat2iait(iat0, &I0, &T0);
        
        // first step: find the adjacent atoms and filter the real adjacent atoms
        AdjacentAtomInfo adjs;
        this->gridD->Find_atom(*ucell, tau0, T0, I0, &adjs);

        std::vector<bool> is_adj(adjs.adj_num + 1, false);
        for (int ad = 0; ad < adjs.adj_num + 1; ++ad)
        {
            const int T1 = adjs.ntype[ad];
            const int I1 = adjs.natom[ad];
            const int iat1 = ucell->itia2iat(T1, I1);
            const ModuleBase::Vector3<int>& R_index1 = adjs.box[ad];
            // choose the real adjacent atoms
            // Note: the distance of atoms should less than the cutoff radius, 
            // When equal, the theoretical value of matrix element is zero, 
            // but the calculated value is not zero due to the numerical error, which would lead to result changes.
            if (this->ucell->cal_dtau(iat0, iat1, R_index1).norm() * this->ucell->lat0
                < this->orb_cutoff_[T1] + PARAM.inp.onsite_radius)
            {
                is_adj[ad] = true;
            }
        }
        filter_adjs(is_adj, adjs);
        const int max_l_plus_1 = this->ucell->atoms[T0].nwl + 1;
        const int length = max_l_plus_1 * max_l_plus_1;
        std::vector<std::unordered_map<int, std::vector<double>>> nlm_iat0(adjs.adj_num + 1);

        for (int ad = 0; ad < adjs.adj_num + 1; ++ad)
        {
            const int T1 = adjs.ntype[ad];
            const int I1 = adjs.natom[ad];
            const int iat1 = ucell->itia2iat(T1, I1);
            const ModuleBase::Vector3<double>& tau1 = adjs.adjacent_tau[ad];
            const Atom* atom1 = &ucell->atoms[T1];

            auto all_indexes = paraV->get_indexes_row(iat1);
            auto col_indexes = paraV->get_indexes_col(iat1);
            // insert col_indexes into all_indexes to get universal set with no repeat elements
            all_indexes.insert(all_indexes.end(), col_indexes.begin(), col_indexes.end());
            std::sort(all_indexes.begin(), all_indexes.end());
            all_indexes.erase(std::unique(all_indexes.begin(), all_indexes.end()), all_indexes.end());
            for (int iw1l = 0; iw1l < all_indexes.size(); iw1l += npol)
            {
                const int iw1 = all_indexes[iw1l] / npol;
                std::vector<std::vector<double>> nlm;
                // nlm is a vector of vectors, but size of outer vector is only 1 here
                // If we are calculating force, we need also to store the gradient
                // and size of outer vector is then 4
                // inner loop : all projectors (L0,M0)
                int L1 = atom1->iw2l[iw1];
                int N1 = atom1->iw2n[iw1];
                int m1 = atom1->iw2m[iw1];

                // convert m (0,1,...2l) to M (-l, -l+1, ..., l-1, l)
                int M1 = (m1 % 2 == 0) ? -m1 / 2 : (m1 + 1) / 2;

                ModuleBase::Vector3<double> dtau = tau0 - tau1;
                intor_->snap(T1, L1, N1, M1, T0, dtau * this->ucell->lat0, 1 /*cal_deri*/, nlm);
                // select the elements of nlm with target_L
                std::vector<double> nlm_target(length * 4);
                // select the elements of nlm with target_L (0, 1, 2, 3 ...)
                int target_L = 0, index=0;
                for(int iw =0;iw < this->ucell->atoms[T0].nw; iw++)
                {
                    const int L0 = this->ucell->atoms[T0].iw2l[iw];
                    // only the first zeta of each l-orbital is needed
                    if(L0 == target_L)
                    {
                        for(int m = 0; m < 2*L0+1; m++)
                        {
                            for (int n = 0; n < 4; n++) // value, deri_x, deri_y, deri_z
                            {
                                nlm_target[index + n * length] = nlm[n][iw + m];
                            }
                            index++;
                        }
                        target_L++;
                    }
                }
                nlm_iat0[ad].insert({all_indexes[iw1l], nlm_target});
            }
        }      

        // second iteration to calculate force and stress
        for (int ad1 = 0; ad1 < adjs.adj_num + 1; ++ad1)
        {
            const int T1 = adjs.ntype[ad1];
            const int I1 = adjs.natom[ad1];
            const int iat1 = ucell->itia2iat(T1, I1);
            double* force_tmp1 = (cal_force) ? &force(iat1, 0) : nullptr;
            double* force_tmp2 = (cal_force) ? &force(iat0, 0) : nullptr;
            ModuleBase::Vector3<int>& R_index1 = adjs.box[ad1];
            ModuleBase::Vector3<double> dis1 = adjs.adjacent_tau[ad1] - tau0;
            for (int ad2 = 0; ad2 < adjs.adj_num + 1; ++ad2)
            {
                const int T2 = adjs.ntype[ad2];
                const int I2 = adjs.natom[ad2];
                const int iat2 = ucell->itia2iat(T2, I2);
                ModuleBase::Vector3<int>& R_index2 = adjs.box[ad2];
                ModuleBase::Vector3<double> dis2 = adjs.adjacent_tau[ad2] - tau0;
                ModuleBase::Vector3<int> R_vector(R_index2[0] - R_index1[0],
                                                  R_index2[1] - R_index1[1],
                                                  R_index2[2] - R_index1[2]);
                const hamilt::BaseMatrix<double>* tmp = dmR->find_matrix(iat1, iat2, R_vector[0], R_vector[1], R_vector[2]);
                int row_size = paraV->get_row_size();
                int col_size = paraV->get_col_size();
                if(row_size == 0 || col_size == 0)
                {
                    continue;
                }
                // if not found , skip this pair of atoms
                if (tmp != nullptr)
                {
                    // calculate force
                    if (cal_force) {
                        this->cal_force_IJR(iat1,
                                            iat2,
                                            paraV,
                                            nlm_iat0[ad1],
                                            nlm_iat0[ad2],
                                            tmp,
                                            lambda[iat0],
                                            this->nspin,
                                            force_tmp1,
                                            force_tmp2);
                    }

                    // calculate stress
                    if (cal_stress) {
                        this->cal_stress_IJR(iat1,
                                             iat2,
                                             paraV,
                                             nlm_iat0[ad1],
                                             nlm_iat0[ad2],
                                             tmp,
                                             lambda[iat0],
                                             this->nspin,
                                             dis1,
                                             dis2,
                                             stress_tmp.data());
                    }
                }
            }
        }
    }

    if (cal_force)
    {
#ifdef __MPI
        // sum up the occupation matrix
        Parallel_Reduce::reduce_all(force.c, force.nr * force.nc);
#endif
        for (int i = 0; i < force.nr * force.nc; i++)
        {
            force.c[i] *= 2.0;
        }
    }

    // stress renormalization
    if (cal_stress)
    {
#ifdef __MPI
        // sum up the occupation matrix
        Parallel_Reduce::reduce_all(stress_tmp.data(), 6);
#endif
        const double weight = this->ucell->lat0 / this->ucell->omega;
        for (int i = 0; i < 6; i++)
        {
            stress.c[i] = stress_tmp[i] * weight;
        }
        stress.c[8] = stress.c[5]; // stress(2,2)
        stress.c[7] = stress.c[4]; // stress(2,1)
        stress.c[6] = stress.c[2]; // stress(2,0)
        stress.c[5] = stress.c[4]; // stress(1,2)
        stress.c[4] = stress.c[3]; // stress(1,1)
        stress.c[3] = stress.c[1]; // stress(1,0)
    }

    ModuleBase::timer::tick("DeltaSpin", "cal_force_stress");
}

template <typename TK, typename TR>
void DeltaSpin<OperatorLCAO<TK, TR>>::cal_force_IJR(const int& iat1,
                                               const int& iat2,
                                               const Parallel_Orbitals* paraV,
                                               const std::unordered_map<int, std::vector<double>>& nlm1_all,
                                               const std::unordered_map<int, std::vector<double>>& nlm2_all,
                                               const hamilt::BaseMatrix<double>* dmR_pointer,
                                               const ModuleBase::Vector3<double>& lambda,
                                               const int nspin,
                                               double* force1,
                                               double* force2)
{
    // npol is the number of polarizations,
    // 1 for non-magnetic (one Hamiltonian matrix only has spin-up or spin-down),
    // 2 for magnetic (one Hamiltonian matrix has both spin-up and spin-down)
    const int npol = this->ucell->get_npol();
    // ---------------------------------------------
    // calculate the Nonlocal matrix for each pair of orbitals
    // ---------------------------------------------
    auto row_indexes = paraV->get_indexes_row(iat1);
    auto col_indexes = paraV->get_indexes_col(iat2);
    // step_trace = 0 for NSPIN=2; ={0, 1, local_col, local_col+1} for NSPIN=4
    std::vector<int> step_trace(nspin, 0);
    if (nspin == 4) {
        step_trace[1] = 1;
        step_trace[2] = col_indexes.size();
        step_trace[3] = col_indexes.size() + 1;
    }
    double tmp[3];
    // calculate the local matrix
    for (int is = 1; is < nspin; is++)
    {
        const double lambda_tmp = nspin==2?lambda[2]:lambda[is-1];
        const double* dm_pointer = dmR_pointer->get_pointer();
        for (int iw1l = 0; iw1l < row_indexes.size(); iw1l += npol)
        {
            const std::vector<double>& nlm1 = nlm1_all.find(row_indexes[iw1l])->second;
            for (int iw2l = 0; iw2l < col_indexes.size(); iw2l += npol)
            {
                const std::vector<double>& nlm2 = nlm2_all.find(col_indexes[iw2l])->second;
#ifdef __DEBUG
                assert(nlm1.size() == nlm2.size());
#endif
                const int length = nlm1.size() / 4;
                const int lmax = sqrt(length);
                int index = 0;
                for(int l = 0; l<lmax; l++)
                {
                    for (int m = 0; m < 2*l+1; m++)
                    {
                        index = l*l + m;
                        tmp[0] = lambda_tmp * nlm1[index + length] * nlm2[index] * dm_pointer[step_trace[is]];
                        tmp[1] = lambda_tmp * nlm1[index + length * 2] * nlm2[index] * dm_pointer[step_trace[is]];
                        tmp[2] = lambda_tmp * nlm1[index + length * 3] * nlm2[index] * dm_pointer[step_trace[is]];
                        // force1 = - VU * <d phi_{I,R1}/d R1|chi_m> * <chi_m'|phi_{J,R2}>
                        // force2 = - VU * <phi_{I,R1}|d chi_m/d R0> * <chi_m'|phi_{J,R2>}
                        force1[0] += tmp[0];
                        force1[1] += tmp[1];
                        force1[2] += tmp[2];
                        force2[0] -= tmp[0];
                        force2[1] -= tmp[1];
                        force2[2] -= tmp[2];
                    }
                }
                dm_pointer += npol;
            }
            dm_pointer += (npol - 1) * col_indexes.size();
        }
    }
}

template <typename TK, typename TR>
void DeltaSpin<OperatorLCAO<TK, TR>>::cal_stress_IJR(const int& iat1,
                                                const int& iat2,
                                                const Parallel_Orbitals* paraV,
                                                const std::unordered_map<int, std::vector<double>>& nlm1_all,
                                                const std::unordered_map<int, std::vector<double>>& nlm2_all,
                                                const hamilt::BaseMatrix<double>* dmR_pointer,
                                                const ModuleBase::Vector3<double>& lambda,
                                                const int nspin,
                                                const ModuleBase::Vector3<double>& dis1,
                                                const ModuleBase::Vector3<double>& dis2,
                                                double* stress)
{
    // npol is the number of polarizations,
    // 1 for non-magnetic (one Hamiltonian matrix only has spin-up or spin-down),
    // 2 for magnetic (one Hamiltonian matrix has both spin-up and spin-down)
    const int npol = this->ucell->get_npol();
    // ---------------------------------------------
    // calculate the Nonlocal matrix for each pair of orbitals
    // ---------------------------------------------
    auto row_indexes = paraV->get_indexes_row(iat1);
    auto col_indexes = paraV->get_indexes_col(iat2);
    // step_trace = 0 for NSPIN=2; ={0, 1, local_col, local_col+1} for NSPIN=4
    std::vector<int> step_trace(nspin, 0);
    if (nspin == 4) {
        step_trace[1] = 1;
        step_trace[2] = col_indexes.size();
        step_trace[3] = col_indexes.size() + 1;
    }
    // calculate the local matrix
    for (int is = 1; is < nspin; is++)
    {
        const double lambda_tmp = nspin==2?lambda[2]:lambda[is-1];
        const double* dm_pointer = dmR_pointer->get_pointer();
        for (int iw1l = 0; iw1l < row_indexes.size(); iw1l += npol)
        {
            const std::vector<double>& nlm1 = nlm1_all.find(row_indexes[iw1l])->second;
            for (int iw2l = 0; iw2l < col_indexes.size(); iw2l += npol)
            {
                const std::vector<double>& nlm2 = nlm2_all.find(col_indexes[iw2l])->second;
#ifdef __DEBUG
                assert(nlm1.size() == nlm2.size());
#endif
                const int length = nlm1.size() / 4;
                const int lmax = sqrt(length);
                double tmp = lambda_tmp * dm_pointer[step_trace[is]];
                int index = 0;
                for(int l = 0; l<lmax; l++)
                {
                    for (int m = 0; m < 2*l+1; m++)
                    {
                        index = l*l + m;
                        stress[0]
                            += tmp * (nlm1[index + length] * dis1.x * nlm2[index] + nlm1[index] * nlm2[index + length] * dis2.x);
                        stress[1]
                            += tmp * (nlm1[index + length] * dis1.y * nlm2[index] + nlm1[index] * nlm2[index + length] * dis2.y);
                        stress[2]
                            += tmp * (nlm1[index + length] * dis1.z * nlm2[index] + nlm1[index] * nlm2[index + length] * dis2.z);
                        stress[3] += tmp
                                     * (nlm1[index + length * 2] * dis1.y * nlm2[index]
                                        + nlm1[index] * nlm2[index + length * 2] * dis2.y);
                        stress[4] += tmp
                                     * (nlm1[index + length * 2] * dis1.z * nlm2[index]
                                        + nlm1[index] * nlm2[index + length * 2] * dis2.z);
                        stress[5] += tmp
                                     * (nlm1[index + length * 3] * dis1.z * nlm2[index]
                                        + nlm1[index] * nlm2[index + length * 3] * dis2.z);
                    }
                }
                dm_pointer += npol;
            }
            dm_pointer += (npol - 1) * col_indexes.size();
        }
    }
}

} // namespace hamilt
