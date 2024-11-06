#pragma once
#include "nonlocal_new.h"
#include "module_base/parallel_reduce.h"
#include "module_base/timer.h"

namespace hamilt
{

template <typename TK, typename TR>
void NonlocalNew<OperatorLCAO<TK, TR>>::cal_force_stress(const bool cal_force,
                                                  const bool cal_stress,
                                                  const HContainer<TR>* dmR,
                                                  ModuleBase::matrix& force,
                                                  ModuleBase::matrix& stress)
{
    ModuleBase::TITLE("NonlocalNew", "cal_force_stress");

    // begin the calculation of force and stress
    ModuleBase::timer::tick("NonlocalNew", "cal_force_stress");

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
        // skip the atoms without plus-U
        auto tau0 = ucell->get_tau(iat0);
        int I0 = 0;
        ucell->iat2iait(iat0, &I0, &this->current_type);
        
        // first step: find the adjacent atoms and filter the real adjacent atoms
        AdjacentAtomInfo adjs;
        this->gridD->Find_atom(*ucell, tau0, this->current_type, I0, &adjs);

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
                < orb_cutoff_[T1] + this->ucell->infoNL.Beta[this->current_type].get_rcut_max())
            {
                is_adj[ad] = true;
            }
        }
        filter_adjs(is_adj, adjs);

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
                intor_->snap(T1, L1, N1, M1, this->current_type, dtau * this->ucell->lat0, true /*cal_deri*/, nlm);
                // select the elements of nlm with target_L
                const int length = nlm[0].size();
                std::vector<double> nlm_target(length * 4);
                // rearrange the nlm_target to store the gradient
                for(int index =0;index < length; index++)
                {
                    for (int n = 0; n < 4; n++) // value, deri_x, deri_y, deri_z
                    {
                        nlm_target[index + n * length] = nlm[n][index];
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
                const hamilt::BaseMatrix<TR>* tmp = dmR->find_matrix(iat1, iat2, R_vector[0], R_vector[1], R_vector[2]);
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

    ModuleBase::timer::tick("NonlocalNew", "cal_force_stress");
}

template <>
void NonlocalNew<OperatorLCAO<std::complex<double>, std::complex<double>>>::cal_force_IJR(const int& iat1,
                                               const int& iat2,
                                               const Parallel_Orbitals* paraV,
                                               const std::unordered_map<int, std::vector<double>>& nlm1_all,
                                               const std::unordered_map<int, std::vector<double>>& nlm2_all,
                                               const hamilt::BaseMatrix<std::complex<double>>* dmR_pointer,
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
    std::vector<int> step_trace(npol * npol, 0);
    if (npol == 2) {
        step_trace[1] = 1;
        step_trace[2] = col_indexes.size();
        step_trace[3] = col_indexes.size() + 1;
    }
    // calculate the local matrix
    const std::complex<double>* tmp_d = nullptr;
    const std::complex<double>* dm_pointer = dmR_pointer->get_pointer();
    for (int iw1l = 0; iw1l < row_indexes.size(); iw1l += npol)
    {
        const std::vector<double>& nlm1 = nlm1_all.find(row_indexes[iw1l])->second;
        const int length = nlm1.size() / 4;
        for (int iw2l = 0; iw2l < col_indexes.size(); iw2l += npol)
        {
            const std::vector<double>& nlm2 = nlm2_all.find(col_indexes[iw2l])->second;
#ifdef __DEBUG
            assert(nlm1.size() == nlm2.size());
#endif
            std::vector<std::complex<double>> nlm_tmp(12, ModuleBase::ZERO);
            for (int is = 0; is < 4; ++is)
            {
                for (int no = 0; no < this->ucell->atoms[this->current_type].ncpp.non_zero_count_soc[is]; no++)
                {
                    const int p1 = this->ucell->atoms[this->current_type].ncpp.index1_soc[is][no];
                    const int p2 = this->ucell->atoms[this->current_type].ncpp.index2_soc[is][no];
                    this->ucell->atoms[this->current_type].ncpp.get_d(is, p1, p2, tmp_d);
                    nlm_tmp[is*3] +=  nlm1[p1 + length] * nlm2[p2] * (*tmp_d);
                    nlm_tmp[is*3+1] += nlm1[p1 + length * 2] * nlm2[p2] * (*tmp_d);
                    nlm_tmp[is*3+2] += nlm1[p1 + length * 3] * nlm2[p2] * (*tmp_d);
                }
            }
            // calculate the force, transfer nlm_tmp to pauli matrix
            for(int i = 0; i < 3; i++)
            {
                double tmp = (dm_pointer[step_trace[0]] * nlm_tmp[i] 
                        + dm_pointer[step_trace[1]] * nlm_tmp[i+3]
                        + dm_pointer[step_trace[2]] * nlm_tmp[i+6]
                        + dm_pointer[step_trace[3]] * nlm_tmp[i+9]).real();
                force1[i] += tmp;
                force2[i] -= tmp;
            }
            dm_pointer += npol;
        }
        dm_pointer += (npol - 1) * col_indexes.size();
    }
}

template <>
void NonlocalNew<OperatorLCAO<std::complex<double>, std::complex<double>>>::cal_stress_IJR(const int& iat1,
                                                const int& iat2,
                                                const Parallel_Orbitals* paraV,
                                                const std::unordered_map<int, std::vector<double>>& nlm1_all,
                                                const std::unordered_map<int, std::vector<double>>& nlm2_all,
                                                const hamilt::BaseMatrix<std::complex<double>>* dmR_pointer,
                                                const ModuleBase::Vector3<double>& dis1,
                                                const ModuleBase::Vector3<double>& dis2,
                                                double* stress)
{
    // npol is the number of polarizations,
    // 1 for non-magnetic (one Hamiltonian matrix only has spin-up or spin-down),
    // 2 for magnetic (one Hamiltonian matrix has both spin-up and spin-down)
    const int npol = this->ucell->get_npol();
    const int npol2 = npol * npol;
    // ---------------------------------------------
    // calculate the Nonlocal matrix for each pair of orbitals
    // ---------------------------------------------
    auto row_indexes = paraV->get_indexes_row(iat1);
    auto col_indexes = paraV->get_indexes_col(iat2);
    // step_trace = 0 for NSPIN=2; ={0, 1, local_col, local_col+1} for NSPIN=4
    std::vector<int> step_trace(npol2, 0);
    if (npol == 2) {
        step_trace[1] = 1;
        step_trace[2] = col_indexes.size();
        step_trace[3] = col_indexes.size() + 1;
    }
    // calculate the local matrix
    const std::complex<double>* tmp_d = nullptr;
    const std::complex<double>* dm_pointer = dmR_pointer->get_pointer();
    for (int iw1l = 0; iw1l < row_indexes.size(); iw1l += npol)
    {
        const std::vector<double>& nlm1 = nlm1_all.find(row_indexes[iw1l])->second;
        const int length = nlm1.size() / npol2;
        for (int iw2l = 0; iw2l < col_indexes.size(); iw2l += npol)
        {
            const std::vector<double>& nlm2 = nlm2_all.find(col_indexes[iw2l])->second;
#ifdef __DEBUG
            assert(nlm1.size() == nlm2.size());
#endif
            std::vector<std::complex<double>> nlm_tmp(npol2 * 6, ModuleBase::ZERO);
            for (int is = 0; is < 4; ++is)
            {
                for (int no = 0; no < this->ucell->atoms[this->current_type].ncpp.non_zero_count_soc[is]; no++)
                {
                    const int p1 = this->ucell->atoms[this->current_type].ncpp.index1_soc[is][no];
                    const int p2 = this->ucell->atoms[this->current_type].ncpp.index2_soc[is][no];
                    this->ucell->atoms[this->current_type].ncpp.get_d(is, p1, p2, tmp_d);
                    nlm_tmp[is*6] += (nlm1[p1 + length] * dis1.x * nlm2[p2] + nlm1[p1] * nlm2[p2 + length] * dis2.x) * (*tmp_d);
                    nlm_tmp[is*6+1] += (nlm1[p1 + length] * dis1.y * nlm2[p2] + nlm1[p1] * nlm2[p2 + length] * dis2.y) * (*tmp_d);
                    nlm_tmp[is*6+2] += (nlm1[p1 + length] * dis1.z * nlm2[p2] + nlm1[p1] * nlm2[p2 + length] * dis2.z) * (*tmp_d);
                    nlm_tmp[is*6+3] += (nlm1[p1 + length * 2] * dis1.y * nlm2[p2] + nlm1[p1] * nlm2[p2 + length * 2] * dis2.y) * (*tmp_d);
                    nlm_tmp[is*6+4] += (nlm1[p1 + length * 2] * dis1.z * nlm2[p2] + nlm1[p1] * nlm2[p2 + length * 2] * dis2.z) * (*tmp_d);
                    nlm_tmp[is*6+5] += (nlm1[p1 + length * 3] * dis1.z * nlm2[p2] + nlm1[p1] * nlm2[p2 + length * 3] * dis2.z) * (*tmp_d);
                }
            }
            // calculate the force, transfer nlm_tmp to pauli matrix
            for(int i = 0; i < 6; i++)
            {
                stress[i] += (dm_pointer[step_trace[0]] * nlm_tmp[i] 
                        + dm_pointer[step_trace[1]] * nlm_tmp[i+6]
                        + dm_pointer[step_trace[2]] * nlm_tmp[i+12]
                        + dm_pointer[step_trace[3]] * nlm_tmp[i+18]).real();
            }
            dm_pointer += npol;
        }
        dm_pointer += (npol - 1) * col_indexes.size();
    }
}

template <typename TK, typename TR>
void NonlocalNew<OperatorLCAO<TK, TR>>::cal_force_IJR(const int& iat1,
                                               const int& iat2,
                                               const Parallel_Orbitals* paraV,
                                               const std::unordered_map<int, std::vector<double>>& nlm1_all,
                                               const std::unordered_map<int, std::vector<double>>& nlm2_all,
                                               const hamilt::BaseMatrix<TR>* dmR_pointer,
                                               double* force1,
                                               double* force2)
{
    // ---------------------------------------------
    // calculate the Nonlocal matrix for each pair of orbitals
    // ---------------------------------------------
    auto row_indexes = paraV->get_indexes_row(iat1);
    auto col_indexes = paraV->get_indexes_col(iat2);
    // calculate the local matrix
    const double* tmp_d = nullptr;
    const double* dm_pointer = dmR_pointer->get_pointer();
    for (int iw1l = 0; iw1l < row_indexes.size(); iw1l++)
    {
        const std::vector<double>& nlm1 = nlm1_all.find(row_indexes[iw1l])->second;
        const int length = nlm1.size() / 4;
        for (int iw2l = 0; iw2l < col_indexes.size(); iw2l++)
        {
            const std::vector<double>& nlm2 = nlm2_all.find(col_indexes[iw2l])->second;
#ifdef __DEBUG
            assert(nlm1.size() == nlm2.size());
#endif
            std::vector<double> nlm_tmp(3, 0.0);
            for (int no = 0; no < this->ucell->atoms[this->current_type].ncpp.non_zero_count_soc[0]; no++)
            {
                const int p1 = this->ucell->atoms[this->current_type].ncpp.index1_soc[0][no];
                const int p2 = this->ucell->atoms[this->current_type].ncpp.index2_soc[0][no];
                this->ucell->atoms[this->current_type].ncpp.get_d(0, p1, p2, tmp_d);
                nlm_tmp[0] += nlm1[p1 + length] * nlm2[p2] * (*tmp_d);
                nlm_tmp[1] += nlm1[p1 + length * 2] * nlm2[p2] * (*tmp_d);
                nlm_tmp[2] += nlm1[p1 + length * 3] * nlm2[p2] * (*tmp_d);
            }
            // calculate the force, transfer nlm_tmp to pauli matrix
            for(int i = 0; i < 3; i++)
            {
                force1[i] += dm_pointer[0] * nlm_tmp[i];
                force2[i] -= dm_pointer[0] * nlm_tmp[i];
            }
            dm_pointer++;
        }
    }
}

template <typename TK, typename TR>
void NonlocalNew<OperatorLCAO<TK, TR>>::cal_stress_IJR(const int& iat1,
                                                const int& iat2,
                                                const Parallel_Orbitals* paraV,
                                                const std::unordered_map<int, std::vector<double>>& nlm1_all,
                                                const std::unordered_map<int, std::vector<double>>& nlm2_all,
                                                const hamilt::BaseMatrix<TR>* dmR_pointer,
                                                const ModuleBase::Vector3<double>& dis1,
                                                const ModuleBase::Vector3<double>& dis2,
                                                double* stress)
{
    // ---------------------------------------------
    // calculate the Nonlocal matrix for each pair of orbitals
    // ---------------------------------------------
    auto row_indexes = paraV->get_indexes_row(iat1);
    auto col_indexes = paraV->get_indexes_col(iat2);
    // calculate the local matrix
    const double* tmp_d = nullptr;
    const double* dm_pointer = dmR_pointer->get_pointer();
    for (int iw1l = 0; iw1l < row_indexes.size(); iw1l++)
    {
        const std::vector<double>& nlm1 = nlm1_all.find(row_indexes[iw1l])->second;
        const int length = nlm1.size() / 4;
        for (int iw2l = 0; iw2l < col_indexes.size(); iw2l++)
        {
            const std::vector<double>& nlm2 = nlm2_all.find(col_indexes[iw2l])->second;
#ifdef __DEBUG
            assert(nlm1.size() == nlm2.size());
#endif
            std::vector<double> nlm_tmp(6, 0.0);
            for (int no = 0; no < this->ucell->atoms[this->current_type].ncpp.non_zero_count_soc[0]; no++)
            {
                const int p1 = this->ucell->atoms[this->current_type].ncpp.index1_soc[0][no];
                const int p2 = this->ucell->atoms[this->current_type].ncpp.index2_soc[0][no];
                this->ucell->atoms[this->current_type].ncpp.get_d(0, p1, p2, tmp_d);
                nlm_tmp[0] += (nlm1[p1 + length] * dis1.x * nlm2[p2] + nlm1[p1] * nlm2[p2 + length] * dis2.x) * (*tmp_d);
                nlm_tmp[1] += (nlm1[p1 + length] * dis1.y * nlm2[p2] + nlm1[p1] * nlm2[p2 + length] * dis2.y) * (*tmp_d);
                nlm_tmp[2] += (nlm1[p1 + length] * dis1.z * nlm2[p2] + nlm1[p1] * nlm2[p2 + length] * dis2.z) * (*tmp_d);
                nlm_tmp[3] += (nlm1[p1 + length * 2] * dis1.y * nlm2[p2] + nlm1[p1] * nlm2[p2 + length * 2] * dis2.y) * (*tmp_d);
                nlm_tmp[4] += (nlm1[p1 + length * 2] * dis1.z * nlm2[p2] + nlm1[p1] * nlm2[p2 + length * 2] * dis2.z) * (*tmp_d);
                nlm_tmp[5] += (nlm1[p1 + length * 3] * dis1.z * nlm2[p2] + nlm1[p1] * nlm2[p2 + length * 3] * dis2.z) * (*tmp_d);
            }
            // calculate the force, transfer nlm_tmp to pauli matrix
            for(int i = 0; i < 6; i++)
            {
                stress[i] += dm_pointer[0] * nlm_tmp[i];
            }
            dm_pointer++;
        }
    }
}

} // namespace hamilt
