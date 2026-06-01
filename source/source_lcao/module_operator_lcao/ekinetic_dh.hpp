#pragma once
#include "ekinetic.h"
#include "operator_force_stress_utils.hpp"
#include "source_base/timer.h"

namespace hamilt
{

template <typename TK, typename TR>
void EKinetic<OperatorLCAO<TK, TR>>::cal_dH(std::vector<hamilt::HContainer<double>*>& dhR_x,
    std::vector<hamilt::HContainer<double>*>& dhR_y,
    std::vector<hamilt::HContainer<double>*>& dhR_z)
{
    ModuleBase::TITLE("EKinetic", "cal_dH");
    ModuleBase::timer::start("EKinetic", "cal_dH");

    const int nat = this->ucell->nat;
    assert(static_cast<int>(dhR_x.size()) == nat);
    const Parallel_Orbitals* paraV = dhR_x[0]->get_paraV();
    const int npol = this->ucell->get_npol();

    // Pass 1: build the same atom-pair structure in each per-atom-I container
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
            for (int iat = 0; iat < nat; ++iat)
            {
                dhR_x[iat]->insert_pair(ap);
                dhR_y[iat]->insert_pair(ap);
                dhR_z[iat]->insert_pair(ap);
            }
        }
    }

    for (int iat = 0; iat < nat; ++iat)
    {
        dhR_x[iat]->allocate(nullptr, true);
        dhR_y[iat]->allocate(nullptr, true);
        dhR_z[iat]->allocate(nullptr, true);
    }

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

                // d<phi_U|T|phi_V>/dtau_I is nonzero only for I in {U=iat1, V=iat2}:
                //   olm     = <phi_U|T|grad phi_V> -> d/dtau_V  -> container iat2
                //   olm_rev = <grad phi_U|T|phi_V> -> d/dtau_U  -> container iat1
                hamilt::BaseMatrix<double>* mtxU_x = dhR_x[iat1]->find_matrix(iat1, iat2, R_index);
                hamilt::BaseMatrix<double>* mtxU_y = dhR_y[iat1]->find_matrix(iat1, iat2, R_index);
                hamilt::BaseMatrix<double>* mtxU_z = dhR_z[iat1]->find_matrix(iat1, iat2, R_index);
                hamilt::BaseMatrix<double>* mtxV_x = dhR_x[iat2]->find_matrix(iat1, iat2, R_index);
                hamilt::BaseMatrix<double>* mtxV_y = dhR_y[iat2]->find_matrix(iat1, iat2, R_index);
                hamilt::BaseMatrix<double>* mtxV_z = dhR_z[iat2]->find_matrix(iat1, iat2, R_index);

                if (!mtxU_x || !mtxU_y || !mtxU_z || !mtxV_x || !mtxV_y || !mtxV_z)
                {
                    continue;
                }

                double* ptrU_x = mtxU_x->get_pointer();
                double* ptrU_y = mtxU_y->get_pointer();
                double* ptrU_z = mtxU_z->get_pointer();
                double* ptrV_x = mtxV_x->get_pointer();
                double* ptrV_y = mtxV_y->get_pointer();
                double* ptrV_z = mtxV_z->get_pointer();
                const int col_size = mtxU_x->get_col_size();

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

                        // calculate() writes the spatial gradient (d/dx,d/dy,d/dz) into
                        // grad_out[0..2] (here olm[0..2]); olm[3] is unused.
                        // d<phi|T|phi>/dtau_I = -<grad phi|T|phi>  (dtau = -grad); sign
                        // confirmed against the finite-difference reference.
                        // d/dtau_V (I = iat2)
                        ptrV_x[idx] -= olm[0];
                        ptrV_y[idx] -= olm[1];
                        ptrV_z[idx] -= olm[2];
                        // d/dtau_U (I = iat1)
                        ptrU_x[idx] -= olm_rev[0];
                        ptrU_y[idx] -= olm_rev[1];
                        ptrU_z[idx] -= olm_rev[2];
                    }
                }
            }
        }
    }

    ModuleBase::timer::end("EKinetic", "cal_dH");
}

} // namespace hamilt
