#pragma once
#include "ekinetic.h"
#include "operator_force_stress_utils.hpp"
#include "source_base/timer.h"

namespace hamilt
{

template <typename TK, typename TR>
void EKinetic<OperatorLCAO<TK, TR>>::cal_dH(hamilt::HContainer<double>* dhR_x,
                                                  hamilt::HContainer<double>* dhR_y,
                                                  hamilt::HContainer<double>* dhR_z)
{
    ModuleBase::TITLE("EKinetic", "cal_dH");
    ModuleBase::timer::start("EKinetic", "cal_dH");

    const Parallel_Orbitals* paraV = dhR_x->get_paraV();
    const int npol = this->ucell->get_npol();

    const int nat = this->ucell->nat;

    for (int iat1 = 0; iat1 < nat; iat1++)
    {
        auto tau1 = this->ucell->get_tau(iat1);
        int T1 = 0, I1 = 0;
        this->ucell->iat2iait(iat1, &I1, &T1);

        AdjacentAtomInfo adjs;
        this->gridD->Find_atom(*this->ucell, tau1, T1, I1, &adjs);

        for (int ad = 0; ad < adjs.adj_num + 1; ++ad)
        {
            const int T2 = adjs.ntype[ad];
            const int I2 = adjs.natom[ad];
            const int iat2 = this->ucell->itia2iat(T2, I2);
            const ModuleBase::Vector3<int>& R_index = adjs.box[ad];

            ModuleBase::Vector3<double> dtau = this->ucell->cal_dtau(iat1, iat2, R_index);
            if (dtau.norm() * this->ucell->lat0 >= this->orb_cutoff_[T1] + this->orb_cutoff_[T2])
            {
                continue;
            }

            if (paraV->get_row_size(iat1) <= 0 || paraV->get_col_size(iat2) <= 0)
            {
                continue;
            }

            hamilt::AtomPair<double> ap(iat1, iat2, R_index.x, R_index.y, R_index.z, paraV);
            dhR_x->insert_pair(ap);
            dhR_y->insert_pair(ap);
            dhR_z->insert_pair(ap);
        }
    }

    dhR_x->allocate(nullptr, true);
    dhR_y->allocate(nullptr, true);
    dhR_z->allocate(nullptr, true);

#pragma omp parallel
    {
#pragma omp for schedule(dynamic)
        for (int iat1 = 0; iat1 < nat; iat1++)
        {
            auto tau1 = this->ucell->get_tau(iat1);
            int T1 = 0, I1 = 0;
            this->ucell->iat2iait(iat1, &I1, &T1);
            const Atom& atom1 = this->ucell->atoms[T1];

            AdjacentAtomInfo adjs;
            this->gridD->Find_atom(*this->ucell, tau1, T1, I1, &adjs);

            for (int ad = 0; ad < adjs.adj_num + 1; ++ad)
            {
                const int T2 = adjs.ntype[ad];
                const int I2 = adjs.natom[ad];
                const int iat2 = this->ucell->itia2iat(T2, I2);
                const ModuleBase::Vector3<int>& R_index = adjs.box[ad];

                ModuleBase::Vector3<double> dtau = this->ucell->cal_dtau(iat1, iat2, R_index);
                if (dtau.norm() * this->ucell->lat0 >= this->orb_cutoff_[T1] + this->orb_cutoff_[T2])
                {
                    continue;
                }

                hamilt::BaseMatrix<double>* mtx_x = dhR_x->find_matrix(iat1, iat2, R_index);
                hamilt::BaseMatrix<double>* mtx_y = dhR_y->find_matrix(iat1, iat2, R_index);
                hamilt::BaseMatrix<double>* mtx_z = dhR_z->find_matrix(iat1, iat2, R_index);

                if (!mtx_x || !mtx_y || !mtx_z)
                {
                    continue;
                }

                double* ptr_x = mtx_x->get_pointer();
                double* ptr_y = mtx_y->get_pointer();
                double* ptr_z = mtx_z->get_pointer();
                const int row_size = mtx_x->get_row_size();
                const int col_size = mtx_x->get_col_size();

                const Atom& atom2 = this->ucell->atoms[T2];

                auto row_indexes = paraV->get_indexes_row(iat1);
                auto col_indexes = paraV->get_indexes_col(iat2);

                if (row_indexes.size() == 0 || col_indexes.size() == 0)
                {
                    continue;
                }

                double olm[4] = {0, 0, 0, 0};
                double olm_rev[4] = {0, 0, 0, 0};

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

                        const ModuleBase::Vector3<double> dtau_scaled = dtau * this->ucell->lat0;

                        this->intor_->calculate(T1, L1, N1, M1, T2, L2, N2, M2, dtau_scaled, nullptr, olm);

                        const ModuleBase::Vector3<double> dtau_rev = (-1.0) * dtau_scaled;
                        this->intor_->calculate(T2, L2, N2, M2, T1, L1, N1, M1, dtau_rev, nullptr, olm_rev);

                        const int idx = (iw1l / npol) * col_size + (iw2l / npol);

                        ptr_x[idx] += olm[1] + olm_rev[1];
                        ptr_y[idx] += olm[2] + olm_rev[2];
                        ptr_z[idx] += olm[3] + olm_rev[3];
                    }
                }
            }
        }
    }

    ModuleBase::timer::end("EKinetic", "cal_dH");
}

} // namespace hamilt
