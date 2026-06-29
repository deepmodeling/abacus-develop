#pragma once
#include "td_pot_hybrid.h"
#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_base/libm/libm.h"

namespace hamilt
{
template <typename TK, typename TR>
void TD_pot_hybrid<OperatorLCAO<TK, TR>>::cal_force_stress(const bool cal_force,
                                                            const HContainer<TR>* dmR,
                                                            ModuleBase::matrix& force)
{
    if(!cal_force)
    {
        return;
    }
    Et = elecstate::H_TDDFT_pw::Et;
    const Parallel_Orbitals* paraV = dmR->get_paraV();
    #pragma omp parallel
    {
    ModuleBase::matrix force_local(force.nr, force.nc);
    #pragma omp for schedule(dynamic)
    for (int iat1 = 0; iat1 < ucell->nat; iat1++)
    {
        auto tau1 = ucell->get_tau(iat1);
        int T1, I1;
        ucell->iat2iait(iat1, &I1, &T1);
        AdjacentAtomInfo adjs;
        this->gridD->Find_atom(*ucell, tau1, T1, I1, &adjs);
        std::vector<bool> is_adj(adjs.adj_num + 1, false);
        for (int ad1 = 0; ad1 < adjs.adj_num + 1; ++ad1)
        {
            const int T2 = adjs.ntype[ad1];
            const int I2 = adjs.natom[ad1];
            const int iat2 = ucell->itia2iat(T2, I2);
            if (paraV->get_nrow_atom(iat1) <= 0 || paraV->get_ncol_atom(iat2) <= 0)
            {
                continue;
            }
            const ModuleBase::Vector3<int>& R_index2 = adjs.box[ad1];
            // choose the real adjacent atoms
            // Note: the distance of atoms should less than the cutoff radius,
            // When equal, the theoretical value of matrix element is zero,
            // but the calculated value is not zero due to the numerical error, which would lead to result changes.
            if (this->ucell->cal_dtau(iat1, iat2, R_index2).norm() * this->ucell->lat0
                < orb_cutoff_[T1] + orb_cutoff_[T2])
            {
                is_adj[ad1] = true;
            }
        }
        filter_adjs(is_adj, adjs);

        for (int ad = 0; ad < adjs.adj_num + 1; ++ad)
        {
            const int T2 = adjs.ntype[ad];
            const int I2 = adjs.natom[ad];
            const int iat2 = ucell->itia2iat(T2, I2);
            double* force_tmp1 = (cal_force) ? &force_local(iat2, 0) : nullptr;
            // force_tmp1[0]=0;
            // force_tmp1[1]=0;
            // force_tmp1[2]=0;
            double* force_tmp2 = (cal_force) ? &force_local(iat1, 0) : nullptr;
            const ModuleBase::Vector3<int>& R_index2 = adjs.box[ad];
            ModuleBase::Vector3<double> dtau = this->ucell->cal_dtau(iat1, iat2, R_index2);
            ModuleBase::Vector3<double> dR = this->ucell->cal_dtau(iat1, iat1, R_index2);

            const BaseMatrix<TR>* tmp = dmR->find_matrix(iat1, iat2, R_index2);
            if (tmp != nullptr)
            {
                this->cal_force_IJR(iat1, iat2, paraV, dtau, dR, tmp->get_pointer(), force_tmp1, force_tmp2);
            }
            else
            {
                ModuleBase::WARNING_QUIT("TD_pot_hybrid::calculate_HR", "R_index not found in HR");
            }
            // std::cout<<"iat1: "<<iat1<< " ad: " <<ad<<std::endl;
            // std::cout<<"T1: "<<T1<< " I1: " << I1 <<" T2: "<<T2<<" I2: "<<I2<<" R_ind_Z: "<<R_index2.z<<std::endl;
            // std::cout<<"force_z: "<<force_tmp1[2]<<std::endl;
            // std::cout<<"dtau: "<<dtau[0]<<" "<<dtau[1]<<" "<<dtau[2]<<std::endl;
        }
    }
    #pragma omp critical
    {
        if(cal_force)
        {
            force += force_local;
        }
    }
    }

    if (cal_force)
    {
#ifdef __MPI
        // sum up the occupation matrix
        Parallel_Reduce::reduce_all(force.c, force.nr * force.nc);
#endif
        // for (int i = 0; i < force.nr * force.nc; i++)
        // {
        //     force.c[i] *= 2.0;
        // }
    }
}
template <>
void TD_pot_hybrid<OperatorLCAO<std::complex<double>, std::complex<double>>>::cal_force_IJR(const int& iat1,
                                                                            const int& iat2,
                                                                            const Parallel_Orbitals* paraV,
                                                                            const ModuleBase::Vector3<double>& dtau,
                                                                            const ModuleBase::Vector3<double>& dR,
                                                                            std::complex<double>* dmR_pointer,
                                                                            double* force1,
                                                                            double* force2)
{
    ModuleBase::WARNING_QUIT("TD_pot_hybrid", "Force calculation for noncollinear spin-polarized systems is not yet supported by the hybrid gauge.");
}
template <typename TK, typename TR>
void TD_pot_hybrid<OperatorLCAO<TK, TR>>::cal_force_IJR(const int& iat1,
                                                        const int& iat2,
                                                        const Parallel_Orbitals* paraV,
                                                        const ModuleBase::Vector3<double>& dtau,
                                                        const ModuleBase::Vector3<double>& dR,
                                                        TR* dmR_pointer,
                                                        double* force1,
                                                        double* force2)
{
    // ---------------------------------------------
    // get info of orbitals of atom1 and atom2 from ucell
    // ---------------------------------------------
    int T1, I1;
    this->ucell->iat2iait(iat1, &I1, &T1);
    int T2, I2;
    this->ucell->iat2iait(iat2, &I2, &T2);
    Atom& atom1 = this->ucell->atoms[T1];
    Atom& atom2 = this->ucell->atoms[T2];

    const int* iw2l1 = atom1.iw2l.data();
    const int* iw2n1 = atom1.iw2n.data();
    const int* iw2m1 = atom1.iw2m.data();
    const int* iw2l2 = atom2.iw2l.data();
    const int* iw2n2 = atom2.iw2n.data();
    const int* iw2m2 = atom2.iw2m.data();

    // ---------------------------------------------
    // calculate the Ekinetic matrix for each pair of orbitals
    // ---------------------------------------------
    double olm[3] = {0, 0, 0};
    auto row_indexes = paraV->get_indexes_row(iat1);
    auto col_indexes = paraV->get_indexes_col(iat2);
    const int step_trace = col_indexes.size() + 1;
    ModuleBase::Vector3<double> shift = {0.0, 0.0, 1.0};
    const ModuleBase::Vector3<double>& tau1 = this->ucell->get_tau(iat1);
    const ModuleBase::Vector3<double> tau2 = tau1 + dtau;
    for (int iw1l = 0; iw1l < row_indexes.size(); iw1l++)
    {
        const int iw1 = row_indexes[iw1l];
        const int L1 = iw2l1[iw1];
        const int N1 = iw2n1[iw1];
        const int m1 = iw2m1[iw1];

        for (int iw2l = 0; iw2l < col_indexes.size(); iw2l++)
        {
            const int iw2 = col_indexes[iw2l];
            const int L2 = iw2l2[iw2];
            const int N2 = iw2n2[iw2];
            const int m2 = iw2m2[iw2];
            
            ModuleBase::Vector3<double> tmp_grad = r_calculator->get_psi_r_gradpsi(tau1 * this->ucell->lat0, T1, L1, m1, N1, tau2 * this->ucell->lat0, T2, L2, m2, N2, Et, dR * this->ucell->lat0);
            ModuleBase::Vector3<double> tmp_grad1 = r_calculator->get_psi_r_gradpsi(tau2 * this->ucell->lat0, T2, L2, m2, N2, tau1 * this->ucell->lat0, T1, L1, m1, N1, Et, dR * this->ucell->lat0);

            for(int i = 0; i < 3; i++)
            {
                force1[i] -= dmR_pointer[0] * tmp_grad[i];
                force2[i] -= dmR_pointer[0] * tmp_grad1[i];
            }
            dmR_pointer++;
        }
    }
}
}// namespace hamilt