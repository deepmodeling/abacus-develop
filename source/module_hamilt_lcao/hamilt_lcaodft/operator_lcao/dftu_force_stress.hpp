#pragma once
#include "dftu_new.h"
#include "module_base/parallel_reduce.h"

namespace hamilt
{

template <typename TK, typename TR>
void DFTUNew<OperatorLCAO<TK, TR>>::cal_force_stress(
  const bool cal_force, 
  const bool cal_stress, 
  ModuleBase::matrix& force, 
  ModuleBase::matrix& stress)
{
    ModuleBase::TITLE("DFTUNew", "cal_force_stress");
#ifdef __DEBUG
    assert(this->dm_in_dftu != nullptr);
#endif    
    ModuleBase::timer::tick("DFTUNew", "cal_force_stress");

    const Parallel_Orbitals* paraV = this->dm_in_dftu->get_DMR_pointer(1)->get_atom_pair(0).get_paraV();
    const int npol = this->ucell->get_npol();
    std::vector<double> stress_tmp(6, 0);
    if (cal_force)
    {
        force.zero_out();
    }
    // 1. calculate <psi|beta> for each pair of atoms
    // loop over all on-site atoms
    for (int iat0 = 0; iat0 < this->ucell->nat; iat0++)
    {
        // skip the atoms without plus-U
        auto tau0 = ucell->get_tau(iat0);
        int T0, I0;
        ucell->iat2iait(iat0, &I0, &T0);
        const int target_L = this->dftu->orbital_corr[T0];
        if(target_L == -1) continue;
        const int tlp1 = 2 * target_L + 1;
        AdjacentAtomInfo& adjs = this->adjs_all[iat0];

        std::vector<std::unordered_map<int, std::vector<double>>> nlm_tot;
        nlm_tot.resize(adjs.adj_num + 1);

        for (int ad = 0; ad < adjs.adj_num + 1; ++ad)
        {
            const int T1 = adjs.ntype[ad];
            const int I1 = adjs.natom[ad];
            const int iat1 = ucell->itia2iat(T1, I1);
            const ModuleBase::Vector3<double>& tau1 = adjs.adjacent_tau[ad];
            const Atom* atom1 = &ucell->atoms[T1];

            const ORB_gen_tables& uot = ORB_gen_tables::get_const_instance();
            const LCAO_Orbitals& orb = LCAO_Orbitals::get_const_instance();
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
#ifdef USE_NEW_TWO_CENTER
                //=================================================================
                //          new two-center integral (temporary)
                //=================================================================
                int L1 = atom1->iw2l[ iw1 ];
                int N1 = atom1->iw2n[ iw1 ];
                int m1 = atom1->iw2m[ iw1 ];

                // convert m (0,1,...2l) to M (-l, -l+1, ..., l-1, l)
                int M1 = (m1 % 2 == 0) ? -m1/2 : (m1+1)/2;

                ModuleBase::Vector3<double> dtau = tau0 - tau1;
                uot.two_center_bundle->overlap_orb_onsite->snap(
                        T1, L1, N1, M1, T0, dtau * this->ucell->lat0, 1 /*cal_deri*/, nlm);
#else
                ModuleBase::WARNING_QUIT("DFTUNew", "old two center integral method not implemented");
#endif
                // select the elements of nlm with target_L
                std::vector<double> nlm_target(tlp1 * 4);
                for(int iw =0;iw < this->ucell->atoms[T0].nw; iw++)
                {
                    const int L0 = this->ucell->atoms[T0].iw2l[iw];
                    if(L0 == target_L)
                    {
                        for(int m = 0; m < tlp1; m++)//-l, -l+1, ..., l-1, l
                        {
                            for(int n = 0; n < 4; n++)// value, deri_x, deri_y, deri_z
                            {
                                nlm_target[m + n * tlp1] = nlm[n][iw+m];
                                //if(dtau.norm2 == 0.0) std::cout<<__FILE__<<__LINE__<<" "<<m<<" "<<n<<" "<<(m+n*tlp1)<<" "<<iw+m<<" "<<nlm[n][iw+m]<<" "<<nlm_target[m + n * tlp1] << std::endl;
                            }
                        }
                        break;
                    }
                }
                nlm_tot[ad].insert({all_indexes[iw1l], nlm_target});
            }
        }
        //first iteration to calculate occupation matrix
        std::vector<double> occ(tlp1 * tlp1 * GlobalV::NSPIN, 0);
        for(int i=0;i<occ.size();i++)
        {
            const int is = i / (tlp1 * tlp1);
            const int ii = i % (tlp1 * tlp1);
            occ[i] = this->dftu->locale[iat0][target_L][0][is].c[ii];
        }
        hamilt::HContainer<double>* dmR_tmp[GlobalV::NSPIN];
        dmR_tmp[0] = this->dm_in_dftu->get_DMR_pointer(1);
        if(GlobalV::NSPIN==2) dmR_tmp[1] = this->dm_in_dftu->get_DMR_pointer(2);

        
        //calculate VU
        const double u_value = this->dftu->U[T0];
        std::vector<double> VU(occ.size());
        double eu_tmp = 0;
        this->cal_v_of_u(occ, tlp1, u_value, &VU[0], eu_tmp);

        // second iteration to calculate force and stress
        // calculate Force for atom J
        //     DMR_{I,J,R'-R} * <phi_{I,R}|chi_m> U*(1/2*delta(m, m')-occ(m, m')) 
        //     \frac{\partial <chi_m'|phi_{J,R'}>}{\partial \tau_J} for each pair of <IJR> atoms
        // calculate Stress for strain tensor \varepsilon_{\alpha\beta}
        //     -1/Omega * DMR_{I,J,R'-R} * [ \frac{\partial <phi_{I,R}|chi_m>}{\partial \tau_{J,\alpha}}\tau_{J,\beta} 
        //     U*(1/2*delta(m, m')-occ(m, m'))<chi_m'|phi_{J,R'}>  
        //   + <phi_{I,R}|chi_m> U*(1/2*delta(m, m')-occ(m, m'))
        //     \frac{\partial <chi_m'|phi_{J,R'}>}{\partial \tau_{J,\alpha}}\tau_{J,\beta}] for each pair of <IJR> atoms
        for (int ad1 = 0; ad1 < adjs.adj_num + 1; ++ad1)
        {
            const int T1 = adjs.ntype[ad1];
            const int I1 = adjs.natom[ad1];
            const int iat1 = ucell->itia2iat(T1, I1);
            double* force_tmp1 = (cal_force)? &force(iat1, 0) : nullptr;
            double* force_tmp2 = (cal_force)? &force(iat0, 0) : nullptr;
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
                const hamilt::BaseMatrix<double>* tmp[GlobalV::NSPIN];
                tmp[0] = dmR_tmp[0]->find_matrix(iat1, iat2, R_vector[0], R_vector[1], R_vector[2]);
                if(GlobalV::NSPIN == 2)
                {
                    tmp[1] = dmR_tmp[1]->find_matrix(iat1, iat2, R_vector[0], R_vector[1], R_vector[2]);
                }
                // if not found , skip this pair of atoms
                if (tmp[0] != nullptr)
                {
                    // calculate force
                    if (cal_force) this->cal_force_IJR(iat1, iat2, T0, paraV, nlm_tot[ad1], nlm_tot[ad2], VU, tmp, GlobalV::NSPIN, force_tmp1, force_tmp2);

                    // calculate stress
                    if (cal_stress) this->cal_stress_IJR(iat1, iat2, T0, paraV, nlm_tot[ad1], nlm_tot[ad2], VU, tmp, GlobalV::NSPIN, dis1, dis2, stress_tmp.data());
                }
            }
        }
    }

    if(cal_force)
    {
#ifdef __MPI
        // sum up the occupation matrix
        Parallel_Reduce::reduce_all(force.c, force.nr*force.nc);
#endif
        for(int i=0;i<force.nr*force.nc;i++)
        {
            force.c[i] *= 2.0;
        }
    }

    // stress renormalization
    if(cal_stress)
    {
#ifdef __MPI
        // sum up the occupation matrix
        Parallel_Reduce::reduce_all(stress_tmp.data(), 6);
#endif
        const double weight = this->ucell->lat0 / this->ucell->omega;
        for(int i=0;i<6;i++)
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

    ModuleBase::timer::tick("DFTUNew", "cal_force_stress");
}

template <typename TK, typename TR>
void DFTUNew<OperatorLCAO<TK, TR>>::cal_force_IJR(
    const int& iat1,
    const int& iat2,
    const int& T0,
    const Parallel_Orbitals* paraV,
    const std::unordered_map<int, std::vector<double>>& nlm1_all,
    const std::unordered_map<int, std::vector<double>>& nlm2_all,
    const std::vector<double>& vu_in,
    const hamilt::BaseMatrix<double>** dmR_pointer,
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
    const int m_size = int(sqrt(vu_in.size()/nspin));
    const int m_size2 = m_size * m_size;
    // step_trace = 0 for NSPIN=1,2; ={0, 1, local_col, local_col+1} for NSPIN=4
    std::vector<int> step_trace(npol, 0);
    if(npol == 2) step_trace[1] = col_indexes.size() + 1;
    double tmp[3];
    // calculate the local matrix
    for (int is = 0; is < nspin; is++)
    {
        double* dm_pointer = dmR_pointer[is]->get_pointer();
        for (int iw1l = 0; iw1l < row_indexes.size(); iw1l += npol)
        {
            const std::vector<double>& nlm1 = nlm1_all.find(row_indexes[iw1l])->second;
            for (int iw2l = 0; iw2l < col_indexes.size(); iw2l += npol)
            {
                const std::vector<double>& nlm2 = nlm2_all.find(col_indexes[iw2l])->second;
    #ifdef __DEBUG
                assert(nlm1.size() == nlm2.size());
    #endif
                for (int m1 = 0; m1 < m_size; m1++)
                {
                    for(int m2 = 0; m2 < m_size; m2++)
                    {
                        tmp[0] = vu_in[m1 * m_size + m2 + is*m_size2] * nlm1[m1 + m_size] * nlm2[m2] * dm_pointer[0];
                        tmp[1] = vu_in[m1 * m_size + m2 + is*m_size2] * nlm1[m1 + m_size*2] * nlm2[m2] * dm_pointer[0];
                        tmp[2] = vu_in[m1 * m_size + m2 + is*m_size2] * nlm1[m1 + m_size*3] * nlm2[m2] * dm_pointer[0];
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
                dm_pointer++;
            }
            dm_pointer += (npol - 1) * col_indexes.size();
        }
    }
}

template <typename TK, typename TR>
void DFTUNew<OperatorLCAO<TK, TR>>::cal_stress_IJR(
    const int& iat1,
    const int& iat2,
    const int& T0,
    const Parallel_Orbitals* paraV,
    const std::unordered_map<int, std::vector<double>>& nlm1_all,
    const std::unordered_map<int, std::vector<double>>& nlm2_all,
    const std::vector<double>& vu_in,
    const hamilt::BaseMatrix<double>** dmR_pointer,
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
    const int m_size = int(sqrt(vu_in.size()/nspin));
    const int m_size2 = m_size * m_size;
    // step_trace = 0 for NSPIN=1,2; ={0, 1, local_col, local_col+1} for NSPIN=4
    std::vector<int> step_trace(npol, 0);
    if(npol == 2) step_trace[1] = col_indexes.size() + 1;
    // calculate the local matrix
    for (int is = 0; is < nspin; is++)
    {
        double* dm_pointer = dmR_pointer[is]->get_pointer();
        for (int iw1l = 0; iw1l < row_indexes.size(); iw1l += npol)
        {
            const std::vector<double>& nlm1 = nlm1_all.find(row_indexes[iw1l])->second;
            for (int iw2l = 0; iw2l < col_indexes.size(); iw2l += npol)
            {
                const std::vector<double>& nlm2 = nlm2_all.find(col_indexes[iw2l])->second;
    #ifdef __DEBUG
                assert(nlm1.size() == nlm2.size());
    #endif
                for (int m1 = 0; m1 < m_size; m1++)
                {
                    for(int m2 = 0; m2 < m_size; m2++)
                    {
                        double tmp = vu_in[m1 * m_size + m2 + is*m_size2] * dm_pointer[0];
                        //std::cout<<__FILE__<<__LINE__<<" "<<tmp<<" "<<m1<<" "<<m2<<" "<<nlm1[m1 + m_size * 2]<<" "<<nlm2[m2 + m_size * 2]<<" "<<dis1.y<<" "<<dis2.y<<std::endl;
                        stress[0] += tmp * (nlm1[m1 + m_size] * dis1.x * nlm2[m2] + nlm1[m1] * nlm2[m2 + m_size] * dis2.x);
                        stress[1] += tmp * (nlm1[m1 + m_size] * dis1.y * nlm2[m2] + nlm1[m1] * nlm2[m2 + m_size] * dis2.y);
                        stress[2] += tmp * (nlm1[m1 + m_size] * dis1.z * nlm2[m2] + nlm1[m1] * nlm2[m2 + m_size] * dis2.z);
                        stress[3] += tmp * (nlm1[m1 + m_size*2] * dis1.y * nlm2[m2] + nlm1[m1] * nlm2[m2 + m_size*2] * dis2.y);
                        stress[4] += tmp * (nlm1[m1 + m_size*2] * dis1.z * nlm2[m2] + nlm1[m1] * nlm2[m2 + m_size*2] * dis2.z);
                        stress[5] += tmp * (nlm1[m1 + m_size*3] * dis1.z * nlm2[m2] + nlm1[m1] * nlm2[m2 + m_size*3] * dis2.z);
                    }
                }
                dm_pointer++;
            }
            dm_pointer += (npol - 1) * col_indexes.size();
        }
    }
}

} // namespace hamilt