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
    bool pre_cal_nlm = false;
    if(ucell->nat<100) // less than 100 atom , cost memory for high performance
    { // pre calculate nlm in initialization
        this->nlm_tot.resize(ucell->nat);
        pre_cal_nlm = true;
    }
    else
    { // calculate nlm on the fly
        this->nlm_tot.resize(1);
    }
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
        if(pre_cal_nlm)
        {
            this->pre_calculate_nlm(iat0, nlm_tot[iat0]);
        }
    }
    // allocate the memory of BaseMatrix in HR, and set the new values to zero
    if(std::is_same<TK, double>::value)
    {
        // only gamma-only has full size of Hamiltonian of DeePKS now, 
        // multi-k keep same size of nonlocal operator, H_V_delta will be allocated by hR
        this->H_V_delta->allocate(nullptr, true);
        // expand hR with H_V_delta, only gamma-only case now
        this->hR->add(*this->H_V_delta);
        this->hR->allocate(nullptr, false);
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
        GlobalC::ld.cal_projected_DM(this->DM,
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

        GlobalC::ld.cal_projected_DM_k(this->DM,
            *this->ucell,
            GlobalC::ORB,
            GlobalC::GridD);
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

        GlobalC::ld.cal_projected_DM_k(this->DM,
            *this->ucell,
            GlobalC::ORB,
            GlobalC::GridD);
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
void hamilt::DeePKS<hamilt::OperatorLCAO<TK, TR>>::pre_calculate_nlm(const int iat0, std::vector<std::unordered_map<int, std::vector<double>>>& nlm_in)
{
    const Parallel_Orbitals* paraV = this->LM->ParaV;
    const int npol = this->ucell->get_npol();
    auto tau0 = ucell->get_tau(iat0);
    int T0, I0;
    ucell->iat2iait(iat0, &I0, &T0);
    AdjacentAtomInfo& adjs = this->adjs_all[iat0];
    nlm_in.resize(adjs.adj_num + 1);

    for (int ad = 0; ad < adjs.adj_num + 1; ++ad)
    {
        const int T1 = adjs.ntype[ad];
        const int I1 = adjs.natom[ad];
        const int iat1 = ucell->itia2iat(T1, I1);
        const ModuleBase::Vector3<double>& tau1 = adjs.adjacent_tau[ad];
        const Atom* atom1 = &ucell->atoms[T1];

        const ORB_gen_tables& uot = ORB_gen_tables::get_const_instance();
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
            nlm_in[ad].insert({all_indexes[iw1l], nlm[0]});
            if(npol == 2) nlm_in[ad].insert({all_indexes[iw1l+1], nlm[0]});
        }
    }
}

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
    for (int iat0 = 0; iat0 < this->ucell->nat; iat0++)
    {
        auto tau0 = ucell->get_tau(iat0);
        int T0, I0;
        ucell->iat2iait(iat0, &I0, &T0);
        AdjacentAtomInfo& adjs = this->adjs_all[iat0];

        //trace alpha orbital
        std::vector<int> trace_alpha_row;
        std::vector<int> trace_alpha_col;
        std::vector<double> gedms;
        int ib=0;
        for (int L0 = 0; L0 <= orb.Alpha[0].getLmax();++L0)
        {
            for (int N0 = 0;N0 < orb.Alpha[0].getNchi(L0);++N0)
            {
                const int inl = GlobalC::ld.get_inl(T0, I0, L0, N0);
                const double* pgedm = GlobalC::ld.get_gedms(inl);
                const int nm = 2*L0+1;
        
                for (int m1=0; m1<nm; ++m1) // m1 = 1 for s, 3 for p, 5 for d
                {
                    for (int m2=0; m2<nm; ++m2) // m1 = 1 for s, 3 for p, 5 for d
                    {
                        trace_alpha_row.push_back(ib+m1);
                        trace_alpha_col.push_back(ib+m2);
                        gedms.push_back(pgedm[m1*nm+m2]);
                    }
                }
                ib+=nm;
            }
        }
        const int trace_alpha_size = trace_alpha_row.size();
        //--------------------------------------------------

        // if nlm_tot is not calculated already, calculate it on the fly now
        int iat00 = iat0;
        if(nlm_tot.size() != this->ucell->nat)
        {
            iat00 = 0;
            nlm_tot[iat00].clear();
            this->pre_calculate_nlm(iat0, nlm_tot[iat00]);
        }
        std::vector<std::unordered_map<int, std::vector<double>>>& nlm_iat = nlm_tot[iat00];

// 2. calculate <psi_I|beta>D<beta|psi_{J,R}> for each pair of <IJR> atoms
        for (int ad1 = 0; ad1 < adjs.adj_num + 1; ++ad1)
        {
            const int T1 = adjs.ntype[ad1];
            const int I1 = adjs.natom[ad1];
            const int iat1 = ucell->itia2iat(T1, I1);
            ModuleBase::Vector3<int>& R_index1 = adjs.box[ad1];
            auto row_indexes = paraV->get_indexes_row(iat1);
            const int row_size = row_indexes.size();
            if(row_size == 0) continue;

            std::vector<double> s_1t(trace_alpha_size * row_size);
            for(int irow=0;irow<row_size;irow++)
            {
                const double* row_ptr = nlm_iat[ad1][row_indexes[irow]].data();
                double* ps1t = &s_1t[irow * trace_alpha_size];
                for(int i=0;i<trace_alpha_size;i++)
                {
                    ps1t[i] = row_ptr[trace_alpha_row[i]] * gedms[i];
                }
            }
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
                if (tmp == nullptr) continue;
                auto col_indexes = paraV->get_indexes_col(iat2);
                const int col_size = col_indexes.size();
                std::vector<double> hr_current(row_size * col_size, 0);
                std::vector<double> s_2t(trace_alpha_size * col_size);
                for(int icol=0;icol<col_size;icol++)
                {
                    const double* col_ptr = nlm_iat[ad2][col_indexes[icol]].data();
                    double* ps2t = &s_2t[icol * trace_alpha_size];
                    for(int i=0;i<trace_alpha_size;i++)
                    {
                        ps2t[i] = col_ptr[trace_alpha_col[i]];
                    }
                }
                /*for(int irow = 0;irow<row_size;irow++)
                {
                    for(int icol=0;icol<col_size;icol++)
                    {
                        for(int ialpha=0;ialpha<trace_alpha_size;ialpha++)
                        {
                            tmp->get_pointer()[irow*col_size+icol] += 
                                s_1t[irow*trace_alpha_size+ialpha] * s_2t[icol*trace_alpha_size+ialpha];
                        }
                    }
                }*/
                //dgemm for s_2t and s_1t to get HR_12
                constexpr char transa='T', transb='N';
                const double gemm_alpha = 1.0, gemm_beta = 1.0;
                dgemm_(
                    &transa, &transb, 
                    &col_size, 
                    &row_size,
                    &trace_alpha_size, 
                    &gemm_alpha, 
                    s_2t.data(), 
                    &trace_alpha_size,
                    s_1t.data(),  
                    &trace_alpha_size,     
                    &gemm_beta,      
                    hr_current.data(),    
                    &col_size);
                // add data of HR to target BaseMatrix
                this->cal_HR_IJR(hr_current.data(), row_size, col_size, tmp->get_pointer());
            }
        }
    }
    ModuleBase::timer::tick("DeePKS", "calculate_HR");
}

// cal_HR_IJR()
template <typename TK, typename TR>
void hamilt::DeePKS<hamilt::OperatorLCAO<TK, TR>>::cal_HR_IJR(
    const double* hr_in,
    const int& row_size,
    const int& col_size,
    TR* data_pointer)
{

    // npol is the number of polarizations,
    // 1 for non-magnetic (one Hamiltonian matrix only has spin-up or spin-down),
    // 2 for magnetic (one Hamiltonian matrix has both spin-up and spin-down)
    const int npol = this->ucell->get_npol();
    // step_trace = 0 for NSPIN=1,2; ={0, 1, local_col, local_col+1} for NSPIN=4
    vector<int> step_trace(2, 0);
    step_trace[1] = col_size + 1;
    // calculate the local matrix
    for (int iw1l = 0; iw1l < row_size; iw1l += npol)
    {
        for (int iw2l = 0; iw2l < col_size; iw2l += npol)
        {
            for (int is = 0; is < npol; ++is)
            {
                data_pointer[step_trace[is]] += TR(*hr_in);
            }
            data_pointer += npol;
            hr_in += npol;
        }
        data_pointer += (npol - 1) * col_size;
        hr_in += (npol - 1) * col_size;
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

template class DeePKS<OperatorLCAO<double, double>>;

template class DeePKS<OperatorLCAO<std::complex<double>, double>>;

template class DeePKS<OperatorLCAO<std::complex<double>, std::complex<double>>>;

}