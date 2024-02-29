#include "dftu_new.h"

#include "module_basis/module_ao/ORB_gen_tables.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_hamilt_lcao/hamilt_lcaodft/operator_lcao/operator_lcao.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer_funcs.h"
#include "module_base/timer.h"
#include "module_base/tool_title.h"
#ifdef _OPENMP
#include <unordered_set>
#endif
#include "module_base/parallel_reduce.h"

template <typename TK, typename TR>
const elecstate::DensityMatrix<TK, double>* hamilt::DFTUNew<hamilt::OperatorLCAO<TK, TR>>::dm_in_dftu = nullptr; 

template <typename TK, typename TR>
hamilt::DFTUNew<hamilt::OperatorLCAO<TK, TR>>::DFTUNew(
    LCAO_Matrix* LM_in,
    const std::vector<ModuleBase::Vector3<double>>& kvec_d_in,
    hamilt::HContainer<TR>* hR_in,
    std::vector<TK>* hK_in,
    const UnitCell* ucell_in,
    Grid_Driver* GridD_in,
    ModuleDFTU::DFTU* dftu_in,
    const Parallel_Orbitals* paraV)
    : hamilt::OperatorLCAO<TK, TR>(LM_in, kvec_d_in, hR_in, hK_in)
{
    this->cal_type = lcao_dftu;
    this->ucell = ucell_in;
    this->dftu = dftu_in;
#ifdef __DEBUG
    assert(this->ucell != nullptr);
#endif
    // initialize HR to allocate sparse Nonlocal matrix memory
    this->initialize_HR(GridD_in, paraV);
}

// destructor
template <typename TK, typename TR>
hamilt::DFTUNew<hamilt::OperatorLCAO<TK, TR>>::~DFTUNew()
{
}

// initialize_HR()
template <typename TK, typename TR>
void hamilt::DFTUNew<hamilt::OperatorLCAO<TK, TR>>::initialize_HR(Grid_Driver* GridD,
                                                                      const Parallel_Orbitals* paraV)
{
    ModuleBase::TITLE("DFTUNew", "initialize_HR");
    ModuleBase::timer::tick("DFTUNew", "initialize_HR");

    this->adjs_all.clear();
    this->adjs_all.reserve(this->ucell->nat);
    for (int iat0 = 0; iat0 < ucell->nat; iat0++)
    {
        auto tau0 = ucell->get_tau(iat0);
        int T0, I0;
        ucell->iat2iait(iat0, &I0, &T0);
        AdjacentAtomInfo adjs;
        GridD->Find_atom(*ucell, tau0, T0, I0, &adjs);
        std::vector<bool> is_adj(adjs.adj_num + 1, false);
        for (int ad1 = 0; ad1 < adjs.adj_num + 1; ++ad1)
        {
            const int T1 = adjs.ntype[ad1];
            const int I1 = adjs.natom[ad1];
            const int iat1 = ucell->itia2iat(T1, I1);
            const ModuleBase::Vector3<double>& tau1 = adjs.adjacent_tau[ad1];
            const ModuleBase::Vector3<int>& R_index1 = adjs.box[ad1];
            // choose the real adjacent atoms
            const LCAO_Orbitals& orb = LCAO_Orbitals::get_const_instance();
            // Note: the distance of atoms should less than the cutoff radius, 
            // When equal, the theoretical value of matrix element is zero, 
            // but the calculated value is not zero due to the numerical error, which would lead to result changes.
            if (this->ucell->cal_dtau(iat0, iat1, R_index1).norm() * this->ucell->lat0
                < orb.Phi[T1].getRcut() + GlobalV::onsite_radius)
            {
                is_adj[ad1] = true;
            }
        }
        filter_adjs(is_adj, adjs);
        this->adjs_all.push_back(adjs);
        
        /*for (int ad1 = 0; ad1 < adjs.adj_num + 1; ++ad1)
        {
            const int T1 = adjs.ntype[ad1];
            const int I1 = adjs.natom[ad1];
            const int iat1 = ucell->itia2iat(T1, I1);
            const ModuleBase::Vector3<int>& R_index1 = adjs.box[ad1];
            for (int ad2 = 0; ad2 < adjs.adj_num + 1; ++ad2)
            {
                const int T2 = adjs.ntype[ad2];
                const int I2 = adjs.natom[ad2];
                const int iat2 = ucell->itia2iat(T2, I2);
                ModuleBase::Vector3<int>& R_index2 = adjs.box[ad2];
                if (paraV->get_col_size(iat2) <= 0 || paraV->get_row_size(iat1) <= 0)
                {
                    continue;
                }
                hamilt::AtomPair<TR> tmp(iat1,
                                         iat2,
                                         R_index2.x - R_index1.x,
                                         R_index2.y - R_index1.y,
                                         R_index2.z - R_index1.z,
                                         paraV);
                this->hR->insert_pair(tmp);
            }
        }*/
    }
    // allocate the memory of BaseMatrix in HR, and set the new values to zero
    //this->hR->allocate(nullptr, true);

    ModuleBase::timer::tick("DFTUNew", "initialize_HR");
}

template <typename TK, typename TR>
void hamilt::DFTUNew<hamilt::OperatorLCAO<TK, TR>>::calculate_HR()
{
    ModuleBase::TITLE("DFTUNew", "calculate_HR");
    if(this->dm_in_dftu == nullptr && this->dftu->initialed_locale == false)
    {// skip the calculation if dm_in_dftu is nullptr
        return;
    }
    else
    {
        //will update this->dftu->locale and this->dftu->EU
        if(GlobalV::CURRENT_SPIN == 0) this->dftu->EU = 0.0;
    }
    ModuleBase::timer::tick("DFTUNew", "calculate_HR");

    const Parallel_Orbitals* paraV = this->hR->get_atom_pair(0).get_paraV();
    const int npol = this->ucell->get_npol();
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
        const int tlp1 = 2*target_L+1;
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
                        T1, L1, N1, M1, T0, dtau * this->ucell->lat0, 0 /*cal_deri*/, nlm);
#else
                ModuleBase::WARNING_QUIT("DFTUNew", "old two center integral method not implemented");
#endif
                // select the elements of nlm with target_L
                std::vector<double> nlm_target(tlp1);
                for(int iw =0;iw < this->ucell->atoms[T0].nw; iw++)
                {
                    const int L0 = this->ucell->atoms[T0].iw2l[iw];
                    if(L0 == target_L)
                    {
                        for(int m = 0; m < 2*L0+1; m++)
                        {
                            nlm_target[m] = nlm[0][iw+m];
                        }
                        break;
                    }
                }
                nlm_tot[ad].insert({all_indexes[iw1l], nlm_target});
            }
        }
        //first iteration to calculate occupation matrix
        const int spin_fold = (GlobalV::NSPIN == 4) ? 4 : 1;
        std::vector<double> occ(tlp1 * tlp1 * spin_fold, 0.0);
        if(this->dftu->initialed_locale == false)
        {
            hamilt::HContainer<double>* dmR_current = this->dm_in_dftu->get_DMR_pointer(GlobalV::CURRENT_SPIN+1);
            for (int ad1 = 0; ad1 < adjs.adj_num + 1; ++ad1)
            {
                const int T1 = adjs.ntype[ad1];
                const int I1 = adjs.natom[ad1];
                const int iat1 = ucell->itia2iat(T1, I1);
                ModuleBase::Vector3<int>& R_index1 = adjs.box[ad1];
                for (int ad2 = 0; ad2 < adjs.adj_num + 1; ++ad2)
                {
                    const int T2 = adjs.ntype[ad2];
                    const int I2 = adjs.natom[ad2];
                    const int iat2 = ucell->itia2iat(T2, I2);
                    ModuleBase::Vector3<int>& R_index2 = adjs.box[ad2];
                    ModuleBase::Vector3<int> R_vector(R_index2[0] - R_index1[0],
                                                    R_index2[1] - R_index1[1],
                                                    R_index2[2] - R_index1[2]);
                    const hamilt::BaseMatrix<double>* tmp = dmR_current->find_matrix(iat1, iat2, R_vector[0], R_vector[1], R_vector[2]);
                    // if not found , skip this pair of atoms
                    if (tmp != nullptr)
                    {
                        this->cal_occupations(iat1, iat2, T0, paraV, nlm_tot[ad1], nlm_tot[ad2], tmp->get_pointer(), occ);
                    }
                }
            }
    #ifdef __MPI
            // sum up the occupation matrix
            Parallel_Reduce::reduce_all(occ.data(), occ.size());
    #endif
            // save occ to dftu
            for(int i=0;i<occ.size();i++)
            {
                if(GlobalV::NSPIN==1) occ[i] *= 0.5;
                this->dftu->locale[iat0][target_L][0][GlobalV::CURRENT_SPIN].c[i] = occ[i];
            }
        }
        else // use readin locale to calculate occupation matrix
        {
            for(int i=0;i<occ.size();i++)
            {
                occ[i] = this->dftu->locale[iat0][target_L][0][GlobalV::CURRENT_SPIN].c[i];
            }
            // set initialed_locale to false to avoid using readin locale in next iteration
        }
        
        //calculate VU
        const double u_value = this->dftu->U[T0];
        std::vector<double> VU_tmp(occ.size());
        this->cal_v_of_u(occ, tlp1, u_value, VU_tmp.data(), this->dftu->EU);
        // transfer occ from pauli matrix format to normal format
        std::vector<TR> VU(occ.size());
        this->transfer_vu(VU_tmp, VU);

        // second iteration to calculate Hamiltonian matrix
        // calculate <psi_I|beta_m> U*(1/2*delta(m, m')-occ(m, m')) <beta_m'|psi_{J,R}> for each pair of <IJR> atoms
// 2. calculate <psi_I|beta>D<beta|psi_{J,R}> for each pair of <IJR> atoms
        for (int ad1 = 0; ad1 < adjs.adj_num + 1; ++ad1)
        {
            const int T1 = adjs.ntype[ad1];
            const int I1 = adjs.natom[ad1];
            const int iat1 = ucell->itia2iat(T1, I1);
            ModuleBase::Vector3<int>& R_index1 = adjs.box[ad1];
            for (int ad2 = 0; ad2 < adjs.adj_num + 1; ++ad2)
            {
                const int T2 = adjs.ntype[ad2];
                const int I2 = adjs.natom[ad2];
                const int iat2 = ucell->itia2iat(T2, I2);
                ModuleBase::Vector3<int>& R_index2 = adjs.box[ad2];
                ModuleBase::Vector3<int> R_vector(R_index2[0] - R_index1[0],
                                                  R_index2[1] - R_index1[1],
                                                  R_index2[2] - R_index1[2]);
                hamilt::BaseMatrix<TR>* tmp = this->hR->find_matrix(iat1, iat2, R_vector[0], R_vector[1], R_vector[2]);
                // if not found , skip this pair of atoms
                if (tmp != nullptr)
                {
                    this->cal_HR_IJR(iat1, iat2, T0, paraV, nlm_tot[ad1], nlm_tot[ad2], VU, tmp->get_pointer());
                }
            }
        }
    }

    //energy correction for NSPIN=1
    if(GlobalV::NSPIN==1) this->dftu->EU *= 2.0;
    // for readin onsite_dm, set initialed_locale to false to avoid using readin locale in next iteration
    if(GlobalV::CURRENT_SPIN == GlobalV::NSPIN-1 || GlobalV::NSPIN==4) this->dftu->initialed_locale = false;

    ModuleBase::timer::tick("DFTUNew", "calculate_HR");
}

// cal_HR_IJR()
template <typename TK, typename TR>
void hamilt::DFTUNew<hamilt::OperatorLCAO<TK, TR>>::cal_HR_IJR(
    const int& iat1,
    const int& iat2,
    const int& T0,
    const Parallel_Orbitals* paraV,
    const std::unordered_map<int, std::vector<double>>& nlm1_all,
    const std::unordered_map<int, std::vector<double>>& nlm2_all,
    const std::vector<TR>& VU,
    TR* data_pointer)
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
    const int m_size = int(sqrt(VU.size())/npol);
    // step_trace = 0 for NSPIN=1,2; ={0, 1, local_col, local_col+1} for NSPIN=4
    std::vector<int> step_trace(npol*npol, 0);
    for (int is = 0; is < npol; is++)
    {
        for (int is2 = 0; is2 < npol; is2++)
        {
            step_trace[is * npol + is2] = paraV->get_col_size(iat2) * is + is2;
        }
    }
    // calculate the local matrix
    const TR* tmp_d = nullptr;
    for (int iw1l = 0; iw1l < row_indexes.size(); iw1l += npol)
    {
        const std::vector<double>& nlm1 = nlm1_all.find(row_indexes[iw1l])->second;
        for (int iw2l = 0; iw2l < col_indexes.size(); iw2l += npol)
        {
            const std::vector<double>& nlm2 = nlm2_all.find(col_indexes[iw2l])->second;
#ifdef __DEBUG
            assert(nlm1.size() == nlm2.size());
#endif
            for (int is = 0; is < npol*npol; ++is)
            {
                int start = is * m_size * m_size;
                TR nlm_tmp = TR(0);
                for (int m1 = 0; m1 < m_size; m1++)
                {
                    for(int m2 = 0; m2 < m_size; m2++)
                    {
                        nlm_tmp += nlm1[m1] * nlm2[m2] * VU[m1 * m_size + m2 + start];
                    }
                }
                data_pointer[step_trace[is]] += nlm_tmp;
            }
            data_pointer += npol;
        }
        data_pointer += (npol - 1) * col_indexes.size();
    }
}

template <typename TK, typename TR>
void hamilt::DFTUNew<hamilt::OperatorLCAO<TK, TR>>::cal_occupations(
    const int& iat1,
    const int& iat2,
    const int& T0,
    const Parallel_Orbitals* paraV,
    const std::unordered_map<int, std::vector<double>>& nlm1_all,
    const std::unordered_map<int, std::vector<double>>& nlm2_all,
    const double* dm_pointer,
    std::vector<double>& occ)
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
    const int m_size = int(sqrt(occ.size())/npol);
    const int m_size2 = m_size * m_size;
#ifdef __DEBUG
    assert(m_size * m_size == occ.size());
#endif
    // step_trace = 0 for NSPIN=1,2; ={0, 1, local_col, local_col+1} for NSPIN=4
    std::vector<int> step_trace(npol * npol, 0);
    for (int is = 0; is < npol; is++)
    {
        for (int is2 = 0; is2 < npol; is2++)
        {
            step_trace[is * npol + is2] = paraV->get_col_size(iat2) * is + is2;
        }
    }
    // calculate the local matrix
    for (int iw1l = 0; iw1l < row_indexes.size(); iw1l += npol)
    {
        const std::vector<double>& nlm1 = nlm1_all.find(row_indexes[iw1l])->second;
        for (int iw2l = 0; iw2l < col_indexes.size(); iw2l += npol)
        {
            const std::vector<double>& nlm2 = nlm2_all.find(col_indexes[iw2l])->second;
#ifdef __DEBUG
            assert(nlm1.size() == nlm2.size());
#endif
            for (int is1 = 0; is1 < npol; ++is1)
            {
                for(int is2 = 0;is2 < npol; ++is2)
                {
                    for (int m1 = 0; m1 < m_size; ++m1)
                    {
                        for(int m2 = 0; m2 < m_size; ++m2)
                        {
                            occ[m1 * m_size + m2 + (is1 * npol + is2) * m_size2] += 
                                nlm1[m1] * nlm2[m2] * dm_pointer[step_trace[is1 * npol + is2]];
                        }
                    }
                }
            }
            dm_pointer += npol;
        }
        dm_pointer += (npol - 1) * col_indexes.size();
    }
}

template<typename TK, typename TR>
void hamilt::DFTUNew<hamilt::OperatorLCAO<TK, TR>>::transfer_vu(std::vector<double>& vu_tmp, std::vector<TR>& vu)
{
#ifdef __DEBUG
    assert(vu.size() == vu_tmp.size());
#endif
    for(int i=0;i<vu_tmp.size();i++)
    {
        vu[i] = vu_tmp[i];
    }
}

template<>
void hamilt::DFTUNew<hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>>::transfer_vu(std::vector<double>& vu_tmp, std::vector<std::complex<double>>& vu)
{
#ifdef __DEBUG
    assert(vu.size() == vu_tmp.size());
#endif
    
    // TR == std::complex<double> transfer from double to std::complex<double>
    const int m_size = int(sqrt(vu.size())/2);
    const int m_size2 = m_size * m_size;
    vu.resize(vu_tmp.size());
    for(int m1=0;m1<m_size;m1++)
    {
        for(int m2=0;m2<m_size;m2++)
        {
            int index[4];
            index[0] = m1*m_size+m2;
            index[1] = m1*m_size+m2 + m_size2;
            index[2] = m2*m_size+m1 + m_size2 * 2;
            index[3] = m2*m_size+m1 + m_size2 * 3;
            vu[index[0]] = 0.5 * (vu_tmp[index[0]] + vu_tmp[index[3]]);
            vu[index[3]] = 0.5 * (vu_tmp[index[0]] - vu_tmp[index[3]]);
            //vu should be complex<double> type, but here we use double type for test
            vu[index[1]] = 0.5 * (vu_tmp[index[1]] + std::complex<double>(0.0, 1.0) * vu_tmp[index[2]]);
            vu[index[2]] = 0.5 * (vu_tmp[index[1]] - std::complex<double>(0.0, 1.0) * vu_tmp[index[2]]);
        }
    }
}

template <typename TK, typename TR>
void hamilt::DFTUNew<hamilt::OperatorLCAO<TK, TR>>::cal_v_of_u(
    const std::vector<double>& occ,
    const int m_size,
    const double u_value,
    double* VU,
    double& EU)
{
    // calculate the local matrix
    int spin_fold = occ.size() / m_size / m_size;
    if(spin_fold < 4)
    for(int is=0;is<spin_fold;++is)
    {
        int start = is * m_size * m_size;
        for (int m1 = 0; m1 < m_size; m1++)
        {
            for(int m2 = 0; m2 < m_size; m2++)
            {
                VU[start + m1 * m_size + m2] = u_value * (0.5 * (m1 == m2) - occ[start + m2 * m_size + m1]);
                EU += u_value * 0.5 * occ[start + m2 * m_size + m1] * occ[start + m1 * m_size + m2];
            }
        }
    }
    else
    {
        for (int m1 = 0; m1 < m_size; m1++)
        {
            for(int m2 = 0; m2 < m_size; m2++)
            {
                VU[m1 * m_size + m2] = u_value * (1.0 * (m1 == m2) - occ[m2 * m_size + m1]);
                EU += u_value * 0.25 * occ[m2 * m_size + m1] * occ[m1 * m_size + m2];
            }
        }
        for(int is=1;is<spin_fold;++is)
        {
            int start = is * m_size * m_size;
            for (int m1 = 0; m1 < m_size; m1++)
            {
                for(int m2 = 0; m2 < m_size; m2++)
                {
                    VU[start + m1 * m_size + m2] = u_value * (0 - occ[start + m2 * m_size + m1]);
                    EU += u_value * 0.25 * occ[start + m2 * m_size + m1] * occ[start + m1 * m_size + m2];
                }
            }
        }
    }
}

// contributeHR()
template <typename TK, typename TR>
void hamilt::DFTUNew<hamilt::OperatorLCAO<TK, TR>>::contributeHR()
{
    ModuleBase::TITLE("DFTUNew", "contributeHR");
    ModuleBase::timer::tick("DFTUNew", "contributeHR");
    this->calculate_HR();
    ModuleBase::timer::tick("DFTUNew", "contributeHR");
    return;
}

#include "dftu_force_stress.hpp"

template class hamilt::DFTUNew<hamilt::OperatorLCAO<double, double>>;
template class hamilt::DFTUNew<hamilt::OperatorLCAO<std::complex<double>, double>>;
template class hamilt::DFTUNew<hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>>;