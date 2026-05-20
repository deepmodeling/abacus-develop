#pragma once
#include "nonlocal.h"
#include "operator_force_stress_utils.h"
#include "source_base/timer.h"

namespace hamilt
{

template <typename TK, typename TR>
void Nonlocal<OperatorLCAO<TK, TR>>::cal_dH(hamilt::HContainer<double>* dhR_x,
                                                  hamilt::HContainer<double>* dhR_y,
                                                  hamilt::HContainer<double>* dhR_z)
{
    ModuleBase::TITLE("Nonlocal", "cal_dH");
    ModuleBase::timer::start("Nonlocal", "cal_dH");

    const Parallel_Orbitals* paraV = dhR_x->get_paraV();
    const int npol = this->ucell->get_npol();

    const int nat = this->ucell->nat;

    for (int iat0 = 0; iat0 < nat; iat0++)
    {
        auto tau0 = this->ucell->get_tau(iat0);
        int I0 = 0, T0 = 0;
        this->ucell->iat2iait(iat0, &I0, &T0);

        AdjacentAtomInfo adjs;
        this->gridD->Find_atom(*this->ucell, tau0, T0, I0, &adjs);

        std::vector<bool> is_adj(adjs.adj_num + 1, false);
        for (int ad = 0; ad < adjs.adj_num + 1; ++ad)
        {
            const int T1 = adjs.ntype[ad];
            const int I1 = adjs.natom[ad];
            const int iat1 = this->ucell->itia2iat(T1, I1);
            const ModuleBase::Vector3<int>& R_index1 = adjs.box[ad];
            if (this->ucell->cal_dtau(iat0, iat1, R_index1).norm() * this->ucell->lat0
                < this->orb_cutoff_[T1] + this->ucell->infoNL.Beta[T0].get_rcut_max())
            {
                is_adj[ad] = true;
            }
        }

        for (int ad1 = 0; ad1 < adjs.adj_num + 1; ++ad1)
        {
            if (!is_adj[ad1])
                continue;
            const int T1 = adjs.ntype[ad1];
            const int I1 = adjs.natom[ad1];
            const int iat1 = this->ucell->itia2iat(T1, I1);
            const ModuleBase::Vector3<int>& R_index1 = adjs.box[ad1];

            for (int ad2 = 0; ad2 < adjs.adj_num + 1; ++ad2)
            {
                if (!is_adj[ad2])
                    continue;
                const int T2 = adjs.ntype[ad2];
                const int I2 = adjs.natom[ad2];
                const int iat2 = this->ucell->itia2iat(T2, I2);
                const ModuleBase::Vector3<int>& R_index2 = adjs.box[ad2];

                if (paraV->get_row_size(iat1) <= 0 || paraV->get_col_size(iat2) <= 0)
                {
                    continue;
                }

                ModuleBase::Vector3<int> dR(R_index2.x - R_index1.x, R_index2.y - R_index1.y, R_index2.z - R_index1.z);

                hamilt::AtomPair<double> ap(iat1, iat2, dR.x, dR.y, dR.z, paraV);
                dhR_x->insert_pair(ap);
                dhR_y->insert_pair(ap);
                dhR_z->insert_pair(ap);
            }
        }
    }

    dhR_x->allocate(nullptr, true);
    dhR_y->allocate(nullptr, true);
    dhR_z->allocate(nullptr, true);

#pragma omp parallel
    {
#pragma omp for schedule(dynamic)
        for (int iat0 = 0; iat0 < nat; iat0++)
        {
            auto tau0 = this->ucell->get_tau(iat0);
            int I0 = 0, T0 = 0;
            this->ucell->iat2iait(iat0, &I0, &T0);

            AdjacentAtomInfo adjs;
            this->gridD->Find_atom(*this->ucell, tau0, T0, I0, &adjs);

            std::vector<bool> is_adj(adjs.adj_num + 1, false);
            for (int ad = 0; ad < adjs.adj_num + 1; ++ad)
            {
                const int T1 = adjs.ntype[ad];
                const int I1 = adjs.natom[ad];
                const int iat1 = this->ucell->itia2iat(T1, I1);
                const ModuleBase::Vector3<int>& R_index1 = adjs.box[ad];
                if (this->ucell->cal_dtau(iat0, iat1, R_index1).norm() * this->ucell->lat0
                    < this->orb_cutoff_[T1] + this->ucell->infoNL.Beta[T0].get_rcut_max())
                {
                    is_adj[ad] = true;
                }
            }

            std::vector<std::unordered_map<int, std::vector<double>>> nlm_iat0(adjs.adj_num + 1);

            for (int ad = 0; ad < adjs.adj_num + 1; ++ad)
            {
                if (!is_adj[ad])
                    continue;

                const int T1 = adjs.ntype[ad];
                const int I1 = adjs.natom[ad];
                const int iat1 = this->ucell->itia2iat(T1, I1);
                const ModuleBase::Vector3<double>& tau1 = adjs.adjacent_tau[ad];
                const Atom* atom1 = &this->ucell->atoms[T1];

                auto all_indexes = paraV->get_indexes_row(iat1);
                auto col_indexes = paraV->get_indexes_col(iat1);
                all_indexes.insert(all_indexes.end(), col_indexes.begin(), col_indexes.end());
                std::sort(all_indexes.begin(), all_indexes.end());
                all_indexes.erase(std::unique(all_indexes.begin(), all_indexes.end()), all_indexes.end());

                for (size_t iw1l = 0; iw1l < all_indexes.size(); iw1l += npol)
                {
                    const int iw1 = all_indexes[iw1l] / npol;
                    std::vector<std::vector<double>> nlm;

                    OperatorForceStress::OrbitalQuantumNumbers qn1 = OperatorForceStress::get_orbital_qn(*atom1, iw1);

                    ModuleBase::Vector3<double> dtau_at = tau0 - tau1;
                    this->intor_->snap(T1, qn1.L, qn1.N, qn1.M, T0, dtau_at * this->ucell->lat0, true, nlm);

                    const size_t length = nlm[0].size();
                    std::vector<double> nlm_target(length * 4);
                    for (size_t index = 0; index < length; index++)
                    {
                        for (int n = 0; n < 4; n++)
                            nlm_target[index + n * length] = nlm[n][index];
                    }
                    nlm_iat0[ad].insert({all_indexes[iw1l], nlm_target});
                }
            }

            for (int ad1 = 0; ad1 < adjs.adj_num + 1; ++ad1)
            {
                if (!is_adj[ad1])
                    continue;
                const int T1 = adjs.ntype[ad1];
                const int I1 = adjs.natom[ad1];
                const int iat1 = this->ucell->itia2iat(T1, I1);
                const ModuleBase::Vector3<int>& R_index1 = adjs.box[ad1];

                for (int ad2 = 0; ad2 < adjs.adj_num + 1; ++ad2)
                {
                    if (!is_adj[ad2])
                        continue;
                    const int T2 = adjs.ntype[ad2];
                    const int I2 = adjs.natom[ad2];
                    const int iat2 = this->ucell->itia2iat(T2, I2);
                    const ModuleBase::Vector3<int>& R_index2 = adjs.box[ad2];

                    ModuleBase::Vector3<int> dR(R_index2.x - R_index1.x,
                                                R_index2.y - R_index1.y,
                                                R_index2.z - R_index1.z);

                    hamilt::BaseMatrix<double>* mtx_x = dhR_x->find_matrix(iat1, iat2, dR.x, dR.y, dR.z);
                    hamilt::BaseMatrix<double>* mtx_y = dhR_y->find_matrix(iat1, iat2, dR.x, dR.y, dR.z);
                    hamilt::BaseMatrix<double>* mtx_z = dhR_z->find_matrix(iat1, iat2, dR.x, dR.y, dR.z);

                    if (!mtx_x || !mtx_y || !mtx_z)
                        continue;

                    double* ptr_x = mtx_x->get_pointer();
                    double* ptr_y = mtx_y->get_pointer();
                    double* ptr_z = mtx_z->get_pointer();
                    const int row_sz = mtx_x->get_row_size();
                    const int col_sz = mtx_x->get_col_size();

                    auto& nlm1_all = nlm_iat0[ad1];
                    auto& nlm2_all = nlm_iat0[ad2];

                    auto row_indexes = paraV->get_indexes_row(iat1);
                    auto col_indexes = paraV->get_indexes_col(iat2);

                    for (size_t iw1l = 0; iw1l < row_indexes.size(); iw1l++)
                    {
                        auto it1 = nlm1_all.find(row_indexes[iw1l]);
                        if (it1 == nlm1_all.end())
                            continue;
                        const std::vector<double>& nlm1 = it1->second;
                        const size_t length = nlm1.size() / 4;
                        const int iw1_row = paraV->global2local_row(row_indexes[iw1l]);

                        for (size_t iw2l = 0; iw2l < col_indexes.size(); iw2l++)
                        {
                            auto it2 = nlm2_all.find(col_indexes[iw2l]);
                            if (it2 == nlm2_all.end())
                                continue;
                            const std::vector<double>& nlm2 = it2->second;
                            const int iw2_col = paraV->global2local_col(col_indexes[iw2l]);

                            std::vector<double> nlm_tmp(3, 0.0);

                            for (int no = 0; no < this->ucell->atoms[T0].ncpp.non_zero_count_soc[0]; no++)
                            {
                                const int p1 = this->ucell->atoms[T0].ncpp.index1_soc[0][no];
                                const int p2 = this->ucell->atoms[T0].ncpp.index2_soc[0][no];
                                const double* tmp_d = nullptr;
                                this->ucell->atoms[T0].ncpp.get_d(0, p1, p2, tmp_d);
                                nlm_tmp[0] += nlm1[p1 + length] * nlm2[p2] * (*tmp_d);
                                nlm_tmp[1] += nlm1[p1 + length * 2] * nlm2[p2] * (*tmp_d);
                                nlm_tmp[2] += nlm1[p1 + length * 3] * nlm2[p2] * (*tmp_d);
                            }

                            int idx = iw1_row * col_sz + iw2_col;
                            ptr_x[idx] += nlm_tmp[0];
                            ptr_y[idx] += nlm_tmp[1];
                            ptr_z[idx] += nlm_tmp[2];
                        }
                    }
                }
            }
        }
    }

    ModuleBase::timer::end("Nonlocal", "cal_dH");
}

} // namespace hamilt
