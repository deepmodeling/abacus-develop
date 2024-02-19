#include "td_nonlocal_lcao.h"

#include "module_basis/module_ao/ORB_gen_tables.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_hamilt_lcao/hamilt_lcaodft/operator_lcao/operator_lcao.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer_funcs.h"
#include "module_base/timer.h"
#include "module_base/tool_title.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#ifdef _OPENMP
#include <unordered_set>
#endif

template <typename TK, typename TR>
hamilt::TDNonlocal<hamilt::OperatorLCAO<TK, TR>>::TDNonlocal(
    LCAO_Matrix* LM_in,
    const std::vector<ModuleBase::Vector3<double>>& kvec_d_in,
    hamilt::HContainer<TR>* hR_in,
    std::vector<TK>* hK_in,
    const UnitCell* ucell_in,
    Grid_Driver* GridD_in,
    const Parallel_Orbitals* paraV)
    : hamilt::OperatorLCAO<TK, TR>(LM_in, kvec_d_in, hR_in, hK_in)
{
    this->cal_type = lcao_tddft_velocity;
    this->ucell = ucell_in;
    this->LM = LM_in;
    this->Grid = GridD_in;
#ifdef __DEBUG
    assert(this->ucell != nullptr);
#endif
    // initialize HR to get adjs info.
    this->init_td();
    this->initialize_HR(Grid,this->LM->ParaV);
}

// destructor
template <typename TK, typename TR>
hamilt::TDNonlocal<hamilt::OperatorLCAO<TK, TR>>::~TDNonlocal()
{
    if (this->allocated)
    {
        delete this->hR_tmp;
    }
}
template <typename TK, typename TR>
void hamilt::TDNonlocal<hamilt::OperatorLCAO<TK, TR>>::init_td(void)
{
    //calculate At in cartesian coorinates.
	double l_norm[3]={GlobalC::ucell.a1.norm() ,GlobalC::ucell.a2.norm() ,GlobalC::ucell.a3.norm()};
    double (&A)[3] = elecstate::H_TDDFT_pw::At;
	cart_At = -(GlobalC::ucell.a1*A[0]/l_norm[0] + GlobalC::ucell.a2*A[1]/l_norm[1] + GlobalC::ucell.a3*A[2]/l_norm[2]);
    std::cout << "cart_At: " << cart_At[0] << " " <<cart_At[1]<< " " << cart_At[2] << std::endl;
}
// initialize_HR()
template <typename TK, typename TR>
void hamilt::TDNonlocal<hamilt::OperatorLCAO<TK, TR>>::initialize_HR(Grid_Driver* GridD,
                                                                      const Parallel_Orbitals* paraV)
{
    if (elecstate::H_TDDFT_pw::stype != 1)
    {
        return;
    }
    ModuleBase::TITLE("TDNonlocal", "initialize_HR");
    ModuleBase::timer::tick("TDNonlocal", "initialize_HR");

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
                < orb.Phi[T1].getRcut() + this->ucell->infoNL.Beta[T0].get_rcut_max())
            {
                is_adj[ad1] = true;
            }
        }
        filter_adjs(is_adj, adjs);
        this->adjs_all.push_back(adjs);
    }

    ModuleBase::timer::tick("TDNonlocal", "initialize_HR");
}

// initialize_HR_tmp()
template <typename TK, typename TR>
void hamilt::TDNonlocal<hamilt::OperatorLCAO<TK, TR>>::initialize_HR_tmp(const Parallel_Orbitals* paraV)
{
    if (elecstate::H_TDDFT_pw::stype != 1)
    {
        return;
    }
    ModuleBase::TITLE("TDNonlocal", "initialize_HR_tmp");
    ModuleBase::timer::tick("TDNonlocal", "initialize_HR_tmp");
    
    for (int i = 0; i < this->hR->size_atom_pairs(); ++i)
    {
        hamilt::AtomPair<TR>& tmp = this->hR->get_atom_pair(i);
        for(int ir = 0;ir < tmp.get_R_size(); ++ir )
        {
            const int* R_index = tmp.get_R_index(ir);
            const int iat1 = tmp.get_atom_i();
            const int iat2 = tmp.get_atom_j(); 

            hamilt::AtomPair<std::complex<double>> tmp1(iat1, iat2, R_index[0], R_index[1], R_index[2], paraV);
            this->hR_tmp->insert_pair(tmp1);
        }
    }
    this->hR_tmp->allocate(nullptr,true);

    ModuleBase::timer::tick("TDNonlocal", "initialize_HR_tmp");
}

template <typename TK, typename TR>
void hamilt::TDNonlocal<hamilt::OperatorLCAO<TK, TR>>::calculate_HR()
{
    ModuleBase::TITLE("TDNonlocal", "calculate_HR");
    ModuleBase::timer::tick("TDNonlocal", "calculate_HR");

    const Parallel_Orbitals* paraV = this->hR_tmp->get_atom_pair(0).get_paraV();
    const int npol = this->ucell->get_npol();
    // 1. calculate <psi|beta> for each pair of atoms
#ifdef _OPENMP
#pragma omp parallel
{
    std::unordered_set<int> atom_row_list;
    #pragma omp for
    for (int iat0 = 0; iat0 < this->ucell->nat; iat0++)
    {
        atom_row_list.insert(iat0);
    }
#endif
    for (int iat0 = 0; iat0 < this->ucell->nat; iat0++)
    {
        auto tau0 = ucell->get_tau(iat0);
        int T0, I0;
        ucell->iat2iait(iat0, &I0, &T0);
        AdjacentAtomInfo& adjs = this->adjs_all[iat0];

        std::vector<std::unordered_map<int, std::vector<std::complex<double>>>> nlm_tot;
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
#ifdef _OPENMP
            if(atom_row_list.find(iat1) == atom_row_list.end())
            {
                all_indexes.clear();
            }
#endif
            auto col_indexes = paraV->get_indexes_col(iat1);
            // insert col_indexes into all_indexes to get universal set with no repeat elements
            all_indexes.insert(all_indexes.end(), col_indexes.begin(), col_indexes.end());
            std::sort(all_indexes.begin(), all_indexes.end());
            all_indexes.erase(std::unique(all_indexes.begin(), all_indexes.end()), all_indexes.end());
            for (int iw1l = 0; iw1l < all_indexes.size(); iw1l += npol)
            {
                const int iw1 = all_indexes[iw1l] / npol;
                std::vector<std::vector<std::complex<double>>> nlm;
                // nlm is a vector of vectors, but size of outer vector is only 1 here
                // If we are calculating force, we need also to store the gradient
                // and size of outer vector is then 4
                // inner loop : all projectors (L0,M0)
#ifdef USE_NEW_TWO_CENTER
                // snap_psibeta_half_tddft() are used to calculate <psi|exp(iAr)|beta> no matter whether to use new two-center or not for now.
                int L1 = atom1->iw2l[ iw1 ];
                int N1 = atom1->iw2n[ iw1 ];
                int m1 = atom1->iw2m[ iw1 ];

                // convert m (0,1,...2l) to M (-l, -l+1, ..., l-1, l)
                int M1 = (m1 % 2 == 0) ? -m1/2 : (m1+1)/2;

                ModuleBase::Vector3<double> dtau = tau0 - tau1;
                uot.snap_psibeta_half_tddft(orb,
                                            this->ucell->infoNL,
                                            nlm,
                                            tau1 * this->ucell->lat0,
                                            T1,
                                            atom1->iw2l[iw1],
                                            atom1->iw2m[iw1],
                                            atom1->iw2n[iw1],
                                            tau0 * this->ucell->lat0,
                                            T0,
                                            -cart_At,
                                            0);
#else
                uot.snap_psibeta_half_tddft(orb,
                                            this->ucell->infoNL,
                                            nlm,
                                            tau1 * this->ucell->lat0,
                                            T1,
                                            atom1->iw2l[iw1],
                                            atom1->iw2m[iw1],
                                            atom1->iw2n[iw1],
                                            tau0 * this->ucell->lat0,
                                            T0,
                                            -cart_At,
                                            0);
#endif
                nlm_tot[ad].insert({all_indexes[iw1l], nlm[0]});
            }
        }
// 2. calculate <psi_I|beta>D<beta|psi_{J,R}> for each pair of <IJR> atoms
        for (int ad1 = 0; ad1 < adjs.adj_num + 1; ++ad1)
        {
            const int T1 = adjs.ntype[ad1];
            const int I1 = adjs.natom[ad1];
            const int iat1 = ucell->itia2iat(T1, I1);
#ifdef _OPENMP
            if(atom_row_list.find(iat1) == atom_row_list.end())
            {
                continue;
            }
#endif
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
                hamilt::BaseMatrix<std::complex<double>>* tmp = this->hR_tmp->find_matrix(iat1, iat2, R_vector[0], R_vector[1], R_vector[2]);
                // if not found , skip this pair of atoms
                if (tmp != nullptr)
                {
                    this->cal_HR_IJR(iat1, iat2, T0, paraV, nlm_tot[ad1], nlm_tot[ad2], tmp->get_pointer());
                }
            }
        }
    }
#ifdef _OPENMP
}
#endif

    ModuleBase::timer::tick("TDNonlocal", "calculate_HR");
}

// cal_HR_IJR()
template <typename TK, typename TR>
void hamilt::TDNonlocal<hamilt::OperatorLCAO<TK, TR>>::cal_HR_IJR(
    const int& iat1,
    const int& iat2,
    const int& T0,
    const Parallel_Orbitals* paraV,
    const std::unordered_map<int, std::vector<std::complex<double>>>& nlm1_all,
    const std::unordered_map<int, std::vector<std::complex<double>>>& nlm2_all,
    std::complex<double>* data_pointer)
{

    // npol is the number of polarizations,
    // 1 for non-magnetic (one Hamiltonian matrix only has spin-up or spin-down),
    // 2 for magnetic (one Hamiltonian matrix has both spin-up and spin-down)
    const int npol = this->ucell->get_npol();
    // ---------------------------------------------
    // calculate the Nonlocal matrix for each pair of orbitals
    // ---------------------------------------------
    double olm[3] = {0, 0, 0};
    auto row_indexes = paraV->get_indexes_row(iat1);
    auto col_indexes = paraV->get_indexes_col(iat2);
    // step_trace = 0 for NSPIN=1,2; ={0, 1, local_col, local_col+1} for NSPIN=4
    std::vector<int> step_trace(npol * npol, 0);
    for (int is = 0; is < npol; is++)
    {
        for (int is2 = 0; is2 < npol; is2++)
        {
            step_trace[is + is2 * npol] = col_indexes.size() * is + is2;
        }
    }
    // calculate the local matrix
    const std::complex<double>* tmp_d = nullptr;
    for (int iw1l = 0; iw1l < row_indexes.size(); iw1l += npol)
    {
        const std::vector<std::complex<double>>& nlm1 = nlm1_all.find(row_indexes[iw1l])->second;
        for (int iw2l = 0; iw2l < col_indexes.size(); iw2l += npol)
        {
            const std::vector<std::complex<double>>& nlm2 = nlm2_all.find(col_indexes[iw2l])->second;
#ifdef __DEBUG
            assert(nlm1.size() == nlm2.size());
#endif
            for (int is = 0; is < npol * npol; ++is)
            {
                std::complex<double> nlm_tmp = std::complex<double>{0,0};
                for (int no = 0; no < this->ucell->atoms[T0].ncpp.non_zero_count_soc[is]; no++)
                {
                    const int p1 = this->ucell->atoms[T0].ncpp.index1_soc[is][no];
                    const int p2 = this->ucell->atoms[T0].ncpp.index2_soc[is][no];
                    //this->ucell->atoms[T0].ncpp.get_d(is, p1, p2, tmp_d);
                    this->ucell->atoms[T0].ncpp.get_d(is, p2, p1, tmp_d);
                    nlm_tmp += nlm1[p1] * std::conj(nlm2[p2]) * (*tmp_d);
                }
                data_pointer[step_trace[is]] += nlm_tmp;
            }
            data_pointer += npol;
        }
        data_pointer += (npol - 1) * col_indexes.size();
    }
}

// set_hR_tmp()
template <typename TK, typename TR>
void hamilt::TDNonlocal<hamilt::OperatorLCAO<TK, TR>>::set_HR_fixed(void* hR_tmp_in)
{
    this->hR_tmp = static_cast<hamilt::HContainer<std::complex<double>>*>(hR_tmp_in);
    this->allocated = false;
}

// contributeHR()
template <typename TK, typename TR>
void hamilt::TDNonlocal<hamilt::OperatorLCAO<TK, TR>>::contributeHR()
{
    ModuleBase::TITLE("TDNonlocal", "contributeHR");
    ModuleBase::timer::tick("TDNonlocal", "contributeHR");
    if (elecstate::H_TDDFT_pw::stype != 1)
    {
        return;
    }
    if (!this->hR_tmp_done)
    {   
        if (this->hR_tmp == nullptr)
        {
            this->hR_tmp = new hamilt::HContainer<std::complex<double>>(this->LM->ParaV);
            // allocate memory for hR_tmp use the same memory as hR
            this->initialize_HR_tmp(this->LM->ParaV);
            this->allocated = true;
        }
        if(this->next_sub_op != nullptr)
        {
            // pass pointer of hR_tmp to the next node
            static_cast<OperatorLCAO<TK, TR>*>(this->next_sub_op)->set_HR_fixed(this->hR_tmp);
        }
        // calculate the values in hR_tmp
        this->calculate_HR();
        this->hR_tmp_done = true;
    }
    ModuleBase::timer::tick("TDNonlocal", "contributeHR");
    return;
}
template<typename TK, typename TR>
void hamilt::TDNonlocal<hamilt::OperatorLCAO<TK, TR>>::contributeHk(int ik)
{
    return;
}
template<>
void hamilt::TDNonlocal<hamilt::OperatorLCAO<std::complex<double>, double>>::contributeHk(int ik)
{
    if (GlobalV::ESOLVER_TYPE != "tddft" || elecstate::H_TDDFT_pw::stype != 1)
    {
        return;
    }
    else{        
        ModuleBase::TITLE("TDNonlocal", "contributeHk");
        ModuleBase::timer::tick("TDNonlocal", "contributeHk");
        //folding inside HR to HK
        if(ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER())
        {
            const int nrow = this->LM->ParaV->get_row_size();
            folding_HR(*this->hR_tmp, this->hK->data(), this->kvec_d[ik], nrow, 1);
        }
        else
        {
            const int ncol = this->LM->ParaV->get_col_size();
            folding_HR(*this->hR_tmp, this->hK->data(), this->kvec_d[ik], ncol, 0);
        }
        
        ModuleBase::timer::tick("TDNonlocal", "contributeHk");
    }
}
template class hamilt::TDNonlocal<hamilt::OperatorLCAO<double, double>>;
template class hamilt::TDNonlocal<hamilt::OperatorLCAO<std::complex<double>, double>>;
template class hamilt::TDNonlocal<hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>>;