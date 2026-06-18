#include "td_deepks_lcao.h"

#include "source_base/timer.h"
#include "source_base/tool_title.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "source_estate/module_pot/H_TDDFT_pw.h"
#include "source_io/module_parameter/parameter.h"
#include "source_lcao/module_operator_lcao/operator_lcao.h"
#include "source_lcao/module_rt/td_info.h"

#ifdef _OPENMP
#include <omp.h>
#include <unordered_set>
#endif

template <typename TK, typename TR>
cal_r_overlap_R hamilt::TDDeePKS<hamilt::OperatorLCAO<TK, TR>>::r_calculator;
template <typename TK, typename TR>
bool hamilt::TDDeePKS<hamilt::OperatorLCAO<TK, TR>>::init_done = false;

template <typename TK, typename TR>
hamilt::TDDeePKS<hamilt::OperatorLCAO<TK, TR>>::TDDeePKS(HS_Matrix_K<TK>* hsk_in,
                                                        const std::vector<ModuleBase::Vector3<double>>& kvec_d_in,
                                                        hamilt::HContainer<TR>* hR_in,
                                                        const UnitCell* ucell_in,
                                                        const LCAO_Orbitals& orb,
                                                        const Grid_Driver* GridD_in
#ifdef __MLALGO
                                                        ,
                                                        LCAO_Deepks<TK>* ld_in
#endif
                                                        )
    : hamilt::OperatorLCAO<TK, TR>(hsk_in, kvec_d_in, hR_in), orb_(orb)
{
    this->cal_type = calculation_type::lcao_tddft_periodic;
    this->ucell = ucell_in;
    this->Grid = GridD_in;
#ifdef __MLALGO
    this->ld = ld_in;
#endif
#ifdef __DEBUG
    assert(this->ucell != nullptr);
#endif
    // only meaningful in velocity gauge (td_stype=1)
    if (elecstate::H_TDDFT_pw::stype != 1)
    {
        return;
    }
    // initialize HR to get adjs info (same cutoff as DeePKS V_delta).
    this->initialize_HR(this->Grid);
    // initialize the r-weighted projection tables <phi|alpha>, <phi|r|alpha> once.
    if (!init_done)
    {
        r_calculator.init_alpha(*this->ucell, *this->hsk->get_pv(), orb_);
        init_done = true;
    }
}

template <typename TK, typename TR>
hamilt::TDDeePKS<hamilt::OperatorLCAO<TK, TR>>::~TDDeePKS()
{
}

// initialize_HR()
template <typename TK, typename TR>
void hamilt::TDDeePKS<hamilt::OperatorLCAO<TK, TR>>::initialize_HR(const Grid_Driver* GridD)
{
    if (elecstate::H_TDDFT_pw::stype != 1)
    {
        return;
    }
    ModuleBase::TITLE("TDDeePKS", "initialize_HR");
    ModuleBase::timer::start("TDDeePKS", "initialize_HR");

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
            const ModuleBase::Vector3<int>& R_index1 = adjs.box[ad1];
            // choose the real adjacent atoms (same cutoff as DeePKS V_delta:
            // Phi[T1].Rcut + Alpha[0].Rcut, centered on the descriptor atom iat0)
            if (this->ucell->cal_dtau(iat0, iat1, R_index1).norm() * this->ucell->lat0
                < orb_.Phi[T1].getRcut() + orb_.Alpha[0].getRcut())
            {
                is_adj[ad1] = true;
            }
        }
        filter_adjs(is_adj, adjs);
        this->adjs_all.push_back(adjs);
    }

    ModuleBase::timer::end("TDDeePKS", "initialize_HR");
}

template <typename TK, typename TR>
void hamilt::TDDeePKS<hamilt::OperatorLCAO<TK, TR>>::calculate_current()
{
#ifdef __MLALGO
    ModuleBase::TITLE("TDDeePKS", "calculate_current");

    // current_term is only allocated when out_current==1; nothing to do otherwise.
    if (TD_info::out_current != 1)
    {
        return;
    }
    // DeePKS rt-TDDFT current is currently only implemented for nspin=1.
    // nspin=2 needs the spin-resolved (gedm +/- gedm_mag) handling and nspin=4
    // needs the npol=2 spinor layout, neither of which is validated here.
    if (PARAM.inp.nspin != 1)
    {
        ModuleBase::WARNING_QUIT("TDDeePKS::calculate_current",
                                 "DeePKS rt-TDDFT current (out_current=1) only supports nspin=1 for now.");
    }

    ModuleBase::timer::start("TDDeePKS", "calculate_current");

    const Parallel_Orbitals* paraV = TD_info::td_vel_op->get_current_term_pointer(0)->get_atom_pair(0).get_paraV();
    const int npol = this->ucell->get_npol();

    for (int iat0 = 0; iat0 < this->ucell->nat; iat0++)
    {
        auto tau0 = ucell->get_tau(iat0);
        int T0 = 0;
        int I0 = 0;
        ucell->iat2iait(iat0, &I0, &T0);
        AdjacentAtomInfo& adjs = this->adjs_all[iat0];

        // --------------------------------------------------------------------
        // build the gedm triples (alpha_row, alpha_col, value) for this iat0.
        // The alpha index is ib + m running over (L0 -> N0 -> m), matching the
        // column order of phialpha and of get_psi_r_alpha's output (nlm).
        // Replicated from DeePKS::calculate_HR (deepks_lcao.cpp). nspin=1: no mag.
        // --------------------------------------------------------------------
        std::vector<int> trace_alpha_row;
        std::vector<int> trace_alpha_col;
        std::vector<double> gedms;
        if (!PARAM.inp.deepks_equiv)
        {
            int ib = 0;
            for (int L0 = 0; L0 <= orb_.Alpha[0].getLmax(); ++L0)
            {
                for (int N0 = 0; N0 < orb_.Alpha[0].getNchi(L0); ++N0)
                {
                    const int inl = this->ld->deepks_param.inl_index[T0](I0, L0, N0);
                    const double* pgedm = this->ld->gedm[inl];
                    const int nm = 2 * L0 + 1;
                    for (int m1 = 0; m1 < nm; ++m1)
                    {
                        for (int m2 = 0; m2 < nm; ++m2)
                        {
                            trace_alpha_row.push_back(ib + m1);
                            trace_alpha_col.push_back(ib + m2);
                            gedms.push_back(pgedm[m1 * nm + m2]);
                        }
                    }
                    ib += nm;
                }
            }
        }
        else
        {
            const double* pgedm = this->ld->gedm[iat0];
            int nproj = 0;
            for (int il = 0; il < this->ld->deepks_param.lmaxd + 1; il++)
            {
                nproj += (2 * il + 1) * orb_.Alpha[0].getNchi(il);
            }
            for (int iproj = 0; iproj < nproj; iproj++)
            {
                for (int jproj = 0; jproj < nproj; jproj++)
                {
                    trace_alpha_row.push_back(iproj);
                    trace_alpha_col.push_back(jproj);
                    gedms.push_back(pgedm[iproj * nproj + jproj]);
                }
            }
        }

        // --------------------------------------------------------------------
        // compute <phi|alpha> and <phi|r|alpha> (nlm[0..3]) for every orbital of
        // every neighbor atom of iat0.
        // --------------------------------------------------------------------
        std::vector<std::vector<std::unordered_map<int, std::vector<double>>>> nlm_tot;
        nlm_tot.resize(adjs.adj_num + 1);
        for (int i = 0; i < adjs.adj_num + 1; i++)
        {
            nlm_tot[i].resize(4);
        }
        for (int ad = 0; ad < adjs.adj_num + 1; ++ad)
        {
            const int T1 = adjs.ntype[ad];
            const int I1 = adjs.natom[ad];
            const int iat1 = ucell->itia2iat(T1, I1);
            const ModuleBase::Vector3<double>& tau1 = adjs.adjacent_tau[ad];
            const Atom* atom1 = &ucell->atoms[T1];
            auto all_indexes = paraV->get_indexes_row(iat1);
            auto col_indexes = paraV->get_indexes_col(iat1);
            all_indexes.insert(all_indexes.end(), col_indexes.begin(), col_indexes.end());
            std::sort(all_indexes.begin(), all_indexes.end());
            all_indexes.erase(std::unique(all_indexes.begin(), all_indexes.end()), all_indexes.end());
            for (int iw1l = 0; iw1l < all_indexes.size(); iw1l += npol)
            {
                const int iw1 = all_indexes[iw1l] / npol;
                std::vector<std::vector<double>> nlm;
                // R1 = orbital position (tau1), R2 = descriptor center (tau0)
                r_calculator.get_psi_r_alpha(orb_,
                                             nlm,
                                             tau1 * this->ucell->lat0,
                                             T1,
                                             atom1->iw2l[iw1],
                                             atom1->iw2m[iw1],
                                             atom1->iw2n[iw1],
                                             tau0 * this->ucell->lat0);
                for (int dir = 0; dir < 4; dir++)
                {
                    nlm_tot[ad][dir].insert({all_indexes[iw1l], nlm[dir]});
                }
            }
        }

        // --------------------------------------------------------------------
        // contract over alpha and accumulate -i[V_delta,r] into current_term
        // for each <iat1, iat2, R> atom pair.
        // --------------------------------------------------------------------
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
                const ModuleBase::Vector3<int>& R_index2 = adjs.box[ad2];
                const ModuleBase::Vector3<int> R_vector(R_index2[0] - R_index1[0],
                                                        R_index2[1] - R_index1[1],
                                                        R_index2[2] - R_index1[2]);
                std::complex<double>* tmp_c[3] = {nullptr, nullptr, nullptr};
                for (int i = 0; i < 3; i++)
                {
                    hamilt::BaseMatrix<std::complex<double>>* matrix_ptr
                        = TD_info::td_vel_op->get_current_term_pointer(i)->find_matrix(iat1,
                                                                                       iat2,
                                                                                       R_vector[0],
                                                                                       R_vector[1],
                                                                                       R_vector[2]);
                    if (matrix_ptr != nullptr)
                    {
                        tmp_c[i] = matrix_ptr->get_pointer();
                    }
                }
                // the pair must exist in current_term (it always does: DeePKS
                // expands hR with V_delta_R, and current_term spans hR). Skip if
                // not found, for safety.
                if (tmp_c[0] != nullptr)
                {
                    this->cal_current_IJR(iat1,
                                          iat2,
                                          paraV,
                                          nlm_tot[ad1],
                                          nlm_tot[ad2],
                                          trace_alpha_row,
                                          trace_alpha_col,
                                          gedms,
                                          tmp_c);
                }
            }
        }
    }
    ModuleBase::timer::end("TDDeePKS", "calculate_current");
#endif
}

template <typename TK, typename TR>
void hamilt::TDDeePKS<hamilt::OperatorLCAO<TK, TR>>::cal_current_IJR(
    const int& iat1,
    const int& iat2,
    const Parallel_Orbitals* paraV,
    const std::vector<std::unordered_map<int, std::vector<double>>>& nlm1_all,
    const std::vector<std::unordered_map<int, std::vector<double>>>& nlm2_all,
    const std::vector<int>& trace_alpha_row,
    const std::vector<int>& trace_alpha_col,
    const std::vector<double>& gedms,
    std::complex<double>** current_mat_p)
{
    const int npol = this->ucell->get_npol();
    auto row_indexes = paraV->get_indexes_row(iat1);
    auto col_indexes = paraV->get_indexes_col(iat2);
    // step_trace = 0 for NSPIN=1,2; ={0, 1, local_col, local_col+1} for NSPIN=4
    std::vector<int> step_trace(npol * npol, 0);
    for (int is = 0; is < npol; is++)
    {
        for (int is2 = 0; is2 < npol; is2++)
        {
            step_trace[is * npol + is2] = col_indexes.size() * is + is2;
        }
    }
    const int trace_alpha_size = trace_alpha_row.size();
    for (int iw1l = 0; iw1l < row_indexes.size(); iw1l += npol)
    {
        std::vector<const std::vector<double>*> nlm1;
        for (int dir = 0; dir < 4; dir++)
        {
            nlm1.push_back(&(nlm1_all[dir].find(row_indexes[iw1l])->second));
        }
        for (int iw2l = 0; iw2l < col_indexes.size(); iw2l += npol)
        {
            std::vector<const std::vector<double>*> nlm2;
            for (int dir = 0; dir < 4; dir++)
            {
                nlm2.push_back(&(nlm2_all[dir].find(col_indexes[iw2l])->second));
            }
            for (int is = 0; is < npol * npol; ++is)
            {
                for (int dir = 0; dir < 3; dir++)
                {
                    std::complex<double> nlm_r_tmp = std::complex<double>{0, 0};
                    std::complex<double> imag_unit = std::complex<double>{0, 1};
                    for (int t = 0; t < trace_alpha_size; t++)
                    {
                        const int a1 = trace_alpha_row[t];
                        const int a2 = trace_alpha_col[t];
                        // <phi_i|r|alpha> gedm <alpha|phi_j> - <phi_i|alpha> gedm <alpha|r|phi_j>
                        nlm_r_tmp += (nlm1[dir + 1]->at(a1) * nlm2[0]->at(a2)
                                      - nlm1[0]->at(a1) * nlm2[dir + 1]->at(a2))
                                     * gedms[t];
                    }
                    // -i[r,V_delta], factor 2.0 due to the unit transformation
                    // (identical convention to the Vnl term in cal_vcomm_r_IJR)
                    current_mat_p[dir][step_trace[is]] -= imag_unit * nlm_r_tmp / 2.0;
                }
            }
            for (int dir = 0; dir < 3; dir++)
            {
                current_mat_p[dir] += npol;
            }
        }
        for (int dir = 0; dir < 3; dir++)
        {
            current_mat_p[dir] += (npol - 1) * col_indexes.size();
        }
    }
}

template <typename TK, typename TR>
void hamilt::TDDeePKS<hamilt::OperatorLCAO<TK, TR>>::contributeHR()
{
    ModuleBase::TITLE("TDDeePKS", "contributeHR");

    if (elecstate::H_TDDFT_pw::stype != 1)
    {
        return;
    }

    ModuleBase::timer::start("TDDeePKS", "contributeHR");

    // recompute whenever the TD operators recompute (same gating as TDNonlocal).
    // NOTE: we do NOT reset TD_info::evolve_once here; TDNonlocal (which runs after
    // this operator in the sub-chain) is responsible for resetting it.
    if (!this->current_done || TD_info::evolve_once)
    {
        // forward the shared hR_tmp pointer to the next sub-operator (TDNonlocal),
        // so it accumulates Vnl into the same buffer created by TDEkinetic.
        // This operator itself does not write hR_tmp.
        if (this->next_sub_op != nullptr)
        {
            static_cast<OperatorLCAO<TK, TR>*>(this->next_sub_op)->set_HR_fixed(this->hR_tmp);
        }
        this->calculate_current();
        this->current_done = true;
    }

    ModuleBase::timer::end("TDDeePKS", "contributeHR");
    return;
}

template <typename TK, typename TR>
void hamilt::TDDeePKS<hamilt::OperatorLCAO<TK, TR>>::contributeHk(int ik)
{
    return;
}

template <typename TK, typename TR>
void hamilt::TDDeePKS<hamilt::OperatorLCAO<TK, TR>>::set_HR_fixed(void* hR_tmp_in)
{
    this->hR_tmp = static_cast<hamilt::HContainer<std::complex<double>>*>(hR_tmp_in);
}

template class hamilt::TDDeePKS<hamilt::OperatorLCAO<std::complex<double>, double>>;
template class hamilt::TDDeePKS<hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>>;
