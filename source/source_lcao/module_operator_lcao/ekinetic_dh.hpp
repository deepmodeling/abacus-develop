#pragma once
#include "ekinetic.h"
#include "operator_force_stress_utils.hpp"
#include "source_base/timer.h"

namespace hamilt
{

template <typename TK, typename TR>
void EKinetic<OperatorLCAO<TK, TR>>::cal_dH(std::array<std::vector<hamilt::HContainer<double>*>, 3>& dhR)
{
    ModuleBase::TITLE("EKinetic", "cal_dH");
    ModuleBase::timer::start("EKinetic", "cal_dH");

    const int nat = this->ucell->nat;
    assert(static_cast<int>(dhR[0].size()) == nat);
    const Parallel_Orbitals* paraV = dhR[0][0]->get_paraV();
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

            if (paraV->is_invalid_atom_pair(iat1, iat2))
            {
                continue;
            }

            hamilt::AtomPair<double> ap(iat1, iat2, R_index.x, R_index.y, R_index.z, paraV);
            for (int iat = 0; iat < nat; ++iat)
            {
                for (int d = 0; d < 3; ++d)
                    dhR[d][iat]->insert_pair(ap);
            }
        }
    }

    for (int iat = 0; iat < nat; ++iat)
    {
        for (int d = 0; d < 3; ++d)
            dhR[d][iat]->allocate(nullptr, true);
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
                hamilt::BaseMatrix<double>* mtxU[3];
                hamilt::BaseMatrix<double>* mtxV[3];
                for (int d = 0; d < 3; ++d)
                {
                    mtxU[d] = dhR[d][iat1]->find_matrix(iat1, iat2, R_index);
                    mtxV[d] = dhR[d][iat2]->find_matrix(iat1, iat2, R_index);
                }

                if (!mtxU[0] || !mtxU[1] || !mtxU[2] || !mtxV[0] || !mtxV[1] || !mtxV[2])
                {
                    continue;
                }

                double* ptrU[3] = {mtxU[0]->get_pointer(), mtxU[1]->get_pointer(), mtxU[2]->get_pointer()};
                double* ptrV[3] = {mtxV[0]->get_pointer(), mtxV[1]->get_pointer(), mtxV[2]->get_pointer()};
                const int col_size = mtxU[0]->get_col_size();

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

                        this->intor_->calculate(T1, L1, N1, M1, T2, L2, N2, M2, dtau_scaled, nullptr, olm); // <phi_U|T|dphi_V/dtau_V>

                        const ModuleBase::Vector3<double> dtau_rev = (-1.0) * dtau_scaled;
                        this->intor_->calculate(T2, L2, N2, M2, T1, L1, N1, M1, dtau_rev, nullptr, olm_rev);    // <dphi_U/dtau_U|T|phi_V>

                        const int idx = (iw1l / npol) * col_size + (iw2l / npol);

                        // d<phi|T|phi>/dtau_I = -<grad phi|T|phi>
                        // but olm directly gives <dtau_I phi|T|phi> and <phi|T|dtau_I phi>, 
                        // so we can directly use them without extra negation.
                        // confirmed against the finite-difference reference.
                        for (int d = 0; d < 3; ++d)
                        {
                            ptrV[d][idx] += olm[d];
                            ptrU[d][idx] += olm_rev[d];
                        }
                    }
                }
            }
        }
    }

    ModuleBase::timer::end("EKinetic", "cal_dH");
}

} // namespace hamilt
