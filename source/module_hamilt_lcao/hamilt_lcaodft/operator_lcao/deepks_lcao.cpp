#include "deepks_lcao.h"
#include "module_base/timer.h"
#include "module_base/tool_title.h"
#ifdef __DEEPKS
#include "module_hamilt_lcao/module_deepks/LCAO_deepks.h"
#endif
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer_funcs.h"
#ifdef _OPENMP
#include <unordered_set>
#endif

namespace hamilt
{

template class DeePKS<OperatorLCAO<double, double>>;

template class DeePKS<OperatorLCAO<std::complex<double>, double>>;

template class DeePKS<OperatorLCAO<std::complex<double>, std::complex<double>>>;

template <typename TK, typename TR>
DeePKS<OperatorLCAO<TK, TR>>::DeePKS(Local_Orbital_Charge* loc_in,
    LCAO_Matrix* LM_in,
    const std::vector<ModuleBase::Vector3<double>>& kvec_d_in,
    HContainer<TR>* hR_in,
    std::vector<TK>* hK_in,
    const UnitCell* ucell_in,
    Grid_Driver* GridD_in,
    const int& nks_in,
    elecstate::DensityMatrix<TK,double>* DM_in) : loc(loc_in), nks(nks_in), ucell(ucell_in), OperatorLCAO<TK, TR>(LM_in, kvec_d_in, hR_in, hK_in), DM(DM_in)
{
    this->cal_type = lcao_deepks;
#ifdef __DEEPKS
    this->initialize_HR(GridD_in, LM_in->ParaV);
#endif
}

template <typename TK, typename TR>
DeePKS<OperatorLCAO<TK, TR>>::~DeePKS()
{
    if (this->H_V_delta != nullptr)
    {
        delete this->H_V_delta;
    }
}

#ifdef __DEEPKS
// initialize_HR()
template <typename TK, typename TR>
void hamilt::DeePKS<hamilt::OperatorLCAO<TK, TR>>::initialize_HR(Grid_Driver* GridD,
                                                                      const Parallel_Orbitals* paraV)
{
    ModuleBase::TITLE("DeePKS", "initialize_HR");
    ModuleBase::timer::tick("DeePKS", "initialize_HR");

    //this->H_V_delta = new HContainer<TR>(paraV);
    if(std::is_same<TK, double>::value)
    {
        this->H_V_delta = new HContainer<TR>(paraV);
        this->H_V_delta->fix_gamma();
    }

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
                < orb.Phi[T1].getRcut() + orb.Alpha[0].getRcut())
            {
                is_adj[ad1] = true;
            }
        }
        filter_adjs(is_adj, adjs);
        this->adjs_all.push_back(adjs);
        for (int ad1 = 0; ad1 < adjs.adj_num + 1; ++ad1)
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
                if(std::is_same<TK, double>::value)
                {
                    this->H_V_delta->insert_pair(tmp);
                }
            }
        }
    }
    // allocate the memory of BaseMatrix in HR, and set the new values to zero
    if(std::is_same<TK, double>::value)
    {
        this->H_V_delta->allocate(true);
    }

    ModuleBase::timer::tick("DeePKS", "initialize_HR");
}
#endif

template<>
void DeePKS<OperatorLCAO<double, double>>::contributeHR()
{
    ModuleBase::TITLE("DeePKS", "contributeHR");
#ifdef __DEEPKS
    if(GlobalC::ld.get_hr_cal())
    {
        ModuleBase::timer::tick("DeePKS", "contributeHR");
        const Parallel_Orbitals* pv = this->LM->ParaV;
        GlobalC::ld.cal_projected_DM(this->DM->get_DMK_vector(),
            *this->ucell,
            GlobalC::ORB,
            GlobalC::GridD);
        GlobalC::ld.cal_descriptor();        
        GlobalC::ld.cal_gedm(this->ucell->nat);
        //GlobalC::ld.add_v_delta(*this->ucell,
        //    GlobalC::ORB,
        //    GlobalC::GridD);
        // recalculate the H_V_delta
        this->H_V_delta->set_zero();
        this->calculate_HR();

        GlobalC::ld.set_hr_cal(false);

        ModuleBase::timer::tick("DeePKS", "contributeHR");
    }
    // save H_V_delta to hR
    this->hR->add(*this->H_V_delta);
#endif
}

template<>
void DeePKS<OperatorLCAO<std::complex<double>, double>>::contributeHR()
{
#ifdef __DEEPKS
    ModuleBase::TITLE("DeePKS", "contributeHR");
    // if DM_K changed, HR of DeePKS need to refresh.
    // the judgement is based on the status of HR in GlobalC::ld
    // this operator should be informed that DM_K has changed and HR need to recalculate.
    if(GlobalC::ld.get_hr_cal())
    {
        ModuleBase::timer::tick("DeePKS", "contributeHR");

        GlobalC::ld.cal_projected_DM_k(this->DM->get_DMK_vector(),
            *this->ucell,
            GlobalC::ORB,
            GlobalC::GridD,
            this->nks,
            this->kvec_d);
        GlobalC::ld.cal_descriptor();
        // calculate dE/dD
        GlobalC::ld.cal_gedm(this->ucell->nat);

        // calculate H_V_deltaR from saved <alpha(0)|psi(R)>
        //GlobalC::ld.add_v_delta_k(*this->ucell, GlobalC::ORB, GlobalC::GridD, this->LM->ParaV->nnr);
        
        // recalculate the H_V_delta
        if(this->H_V_delta == nullptr)
        {
            this->H_V_delta = new hamilt::HContainer<double>(*this->hR);
        }
        this->H_V_delta->set_zero();
        this->calculate_HR();

        GlobalC::ld.set_hr_cal(false);
        
        ModuleBase::timer::tick("DeePKS", "contributeHR");
    } 
    // save H_V_delta to hR
    this->hR->add(*this->H_V_delta);
#endif
}
template<>
void DeePKS<OperatorLCAO<std::complex<double>, std::complex<double>>>::contributeHR()
{
#ifdef __DEEPKS
    ModuleBase::TITLE("DeePKS", "contributeHR");
    // if DM_K changed, HR of DeePKS need to refresh.
    // the judgement is based on the status of HR in GlobalC::ld
    // this operator should be informed that DM_K has changed and HR need to recalculate.
    if(GlobalC::ld.get_hr_cal())
    {
        ModuleBase::timer::tick("DeePKS", "contributeHR");

        GlobalC::ld.cal_projected_DM_k(this->DM->get_DMK_vector(),
            *this->ucell,
            GlobalC::ORB,
            GlobalC::GridD,
            this->nks,
            this->kvec_d);
        GlobalC::ld.cal_descriptor();
        // calculate dE/dD
        GlobalC::ld.cal_gedm(this->ucell->nat);

        // calculate H_V_deltaR from saved <alpha(0)|psi(R)>
        //GlobalC::ld
        //    .add_v_delta_k(*this->ucell, GlobalC::ORB, GlobalC::GridD, this->LM->ParaV->nnr);
        
        // recalculate the H_V_delta
        if(this->H_V_delta == nullptr)
        {
            this->H_V_delta = new hamilt::HContainer<std::complex<double>>(*this->hR);
        }
        this->H_V_delta->set_zero();
        this->calculate_HR();

        GlobalC::ld.set_hr_cal(false);
        
        ModuleBase::timer::tick("DeePKS", "contributeHR");
    } 
    // save H_V_delta to hR
    this->hR->add(*this->H_V_delta);

#endif
}

#ifdef __DEEPKS
template <typename TK, typename TR>
void hamilt::DeePKS<hamilt::OperatorLCAO<TK, TR>>::calculate_HR()
{
    ModuleBase::TITLE("DeePKS", "calculate_HR");
    if(this->H_V_delta->size_atom_pairs() == 0)
    {
        return;
    }
    ModuleBase::timer::tick("DeePKS", "calculate_HR");

    const Parallel_Orbitals* paraV = this->H_V_delta->get_atom_pair(0).get_paraV();
    const int npol = this->ucell->get_npol();

    const LCAO_Orbitals& orb = LCAO_Orbitals::get_const_instance();

    // 1. calculate <psi|alpha> for each pair of atoms
#ifdef _OPENMP
#pragma omp parallel
{
#endif
    // prepair the vector of Loop-L0-N0
    // calculate the length of Loop-L0-N0
    int L0_size = 0;
    std::vector<int> L0s;
    for (int L0 = 0; L0 <= orb.Alpha[0].getLmax(); ++L0)
    {
        L0_size += orb.Alpha[0].getNchi(L0);
        L0s.insert(L0s.end(), orb.Alpha[0].getNchi(L0), L0);
    }
    std::vector<const double*> gedms(L0_size);

    std::unordered_set<int> atom_row_list;
#ifdef _OPENMP
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

        //--------------------------------------------------
        // prepair the vector of Loop-L0-N0
        int index = 0;
        for (int L0 = 0; L0 <= orb.Alpha[0].getLmax();++L0)
        {
            for (int N0 = 0;N0 < orb.Alpha[0].getNchi(L0);++N0)
            {
                const int inl = GlobalC::ld.get_inl(T0, I0, L0, N0);
                gedms[index] = GlobalC::ld.get_gedms(inl);
                index++;
            }
        }
        //--------------------------------------------------

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
                uot.two_center_bundle->overlap_orb_alpha->snap(
                        T1, L1, N1, M1, 0, dtau * ucell->lat0, 0 /*calc_deri*/, nlm);
#else
                uot.snap_psialpha_half(
                        orb,
						nlm, 0, tau1, T1,
						atom1->iw2l[ iw1 ], // L1
						atom1->iw2m[ iw1 ], // m1
						atom1->iw2n[ iw1 ], // N1
						tau0, T0, I0);
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
                hamilt::BaseMatrix<TR>* tmp = this->H_V_delta->find_matrix(iat1, iat2, R_vector[0], R_vector[1], R_vector[2]);
                // if not found , skip this pair of atoms
                if (tmp != nullptr)
                {
                    this->cal_HR_IJR(iat1, iat2, T0, paraV, nlm_tot[ad1], nlm_tot[ad2], L0s.data(), gedms.data(), L0_size, tmp->get_pointer());
                }
            }
        }
    }
#ifdef _OPENMP
}
#endif
    ModuleBase::timer::tick("DeePKS", "calculate_HR");
}

// cal_HR_IJR()
template <typename TK, typename TR>
void hamilt::DeePKS<hamilt::OperatorLCAO<TK, TR>>::cal_HR_IJR(
    const int& iat1,
    const int& iat2,
    const int& T0,
    const Parallel_Orbitals* paraV,
    const std::unordered_map<int, std::vector<double>>& nlm1_all,
    const std::unordered_map<int, std::vector<double>>& nlm2_all,
    const int* L0s,
    const double** gedms,
    const int size_gedms,
    TR* data_pointer)
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
    vector<int> step_trace(2, 0);
    step_trace[1] = col_indexes.size() + 1;
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
            for (int is = 0; is < npol; ++is)
            {
                int ib=0;
                TR nlm_tmp = TR(0);
                for (int index = 0; index< size_gedms; ++index)
                {
                    const int L0 = L0s[index];
                    const int nm = 2*L0+1;
                    for (int m1 = 0;m1 < nm;++m1)
                    {
                        for (int m2 = 0; m2 < nm; ++m2)
                        {
                            nlm_tmp += gedms[index][m1*nm+m2]*nlm1[ib+m1]*nlm2[ib+m2];
                        }
                    }
                    ib+=nm;
                }
                data_pointer[step_trace[is]] += nlm_tmp;
            }
            data_pointer += npol;
        }
        data_pointer += (npol - 1) * col_indexes.size();
    }
}

inline void get_h_delta_k(int ik, double*& h_delta_k)
{
    h_delta_k = GlobalC::ld.H_V_delta.data();
    return;
}
inline void get_h_delta_k(int ik, std::complex<double>*& h_delta_k)
{
    h_delta_k = GlobalC::ld.H_V_delta_k[ik].data();
    return;
}

// contributeHk()
template<typename TK, typename TR>
void hamilt::DeePKS<hamilt::OperatorLCAO<TK, TR>>::contributeHk(int ik)
{
    ModuleBase::TITLE("DeePKS", "contributeHk");
    ModuleBase::timer::tick("DeePKS", "contributeHk");
    
    TK* h_delta_k = nullptr;
    get_h_delta_k(ik, h_delta_k);
    // set SK to zero and then calculate SK for each k vector
    ModuleBase::GlobalFunc::ZEROS(h_delta_k, this->hK->size());

    if(ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER())
    {
        const int nrow = this->LM->ParaV->get_row_size();
        hamilt::folding_HR(*this->H_V_delta, h_delta_k, this->kvec_d[ik], nrow, 1);
    }
    else
    {
        const int ncol = this->LM->ParaV->get_col_size();
        hamilt::folding_HR(*this->H_V_delta, h_delta_k, this->kvec_d[ik], ncol, 0);
    }
    ModuleBase::timer::tick("DeePKS", "contributeHk");
}

#endif

}