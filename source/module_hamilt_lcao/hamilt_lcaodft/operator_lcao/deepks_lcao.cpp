#include "deepks_lcao.h"

#include "source_base/timer.h"
#include "source_base/tool_title.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_hamilt_lcao/module_deepks/LCAO_deepks.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer_funcs.h"
#include "module_parameter/parameter.h"
#ifdef _OPENMP
#include <unordered_set>
#endif

namespace hamilt
{

template <typename TK, typename TR>
DeePKS<OperatorLCAO<TK, TR>>::DeePKS(HS_Matrix_K<TK>* hsk_in,
                                     const std::vector<ModuleBase::Vector3<double>>& kvec_d_in,
                                     HContainer<TR>* hR_in,
                                     const UnitCell* ucell_in,
                                     const Grid_Driver* GridD_in,
                                     const TwoCenterIntegrator* intor_orb_alpha,
                                     const LCAO_Orbitals* ptr_orb,
                                     const int& nks_in,
                                     elecstate::DensityMatrix<TK, double>* DM_in
#ifdef __MLALGO
                                     ,
                                     LCAO_Deepks<TK>* ld_in
#endif
                                     )
    : OperatorLCAO<TK, TR>(hsk_in, kvec_d_in, hR_in), DM(DM_in), ucell(ucell_in), intor_orb_alpha_(intor_orb_alpha),
      ptr_orb_(ptr_orb), nks(nks_in)
{
    this->cal_type = calculation_type::lcao_deepks;
    this->gd = GridD_in;
#ifdef __MLALGO
    this->ld = ld_in;
    this->initialize_HR(GridD_in);
#endif
}

template <typename TK, typename TR>
DeePKS<OperatorLCAO<TK, TR>>::~DeePKS()
{
    if (this->V_delta_R != nullptr)
    {
        delete this->V_delta_R;
    }
}

#ifdef __MLALGO
// initialize_HR()
template <typename TK, typename TR>
void hamilt::DeePKS<hamilt::OperatorLCAO<TK, TR>>::initialize_HR(const Grid_Driver* GridD)
{
    ModuleBase::TITLE("DeePKS", "initialize_HR");
    ModuleBase::timer::tick("DeePKS", "initialize_HR");

    auto* paraV = this->hR->get_paraV(); // get parallel orbitals from HR
    // TODO: if paraV is nullptr, AtomPair can not use paraV for constructor, I will repair it in the future.

    this->V_delta_R = new HContainer<TR>(paraV);
    if (std::is_same<TK, double>::value)
    {
        // this->V_delta_R = new HContainer<TR>(paraV);
        this->V_delta_R->fix_gamma();
    }

    this->adjs_all.clear();
    this->adjs_all.reserve(this->ucell->nat);

    for (int iat0 = 0; iat0 < ucell->nat; iat0++)
    {
        auto tau0 = ucell->get_tau(iat0);
        int T0 = 0;
        int I0 = 0;
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
            // Note: the distance of atoms should less than the cutoff radius,
            // When equal, the theoretical value of matrix element is zero,
            // but the calculated value is not zero due to the numerical error, which would lead to result changes.
            if (this->ucell->cal_dtau(iat0, iat1, R_index1).norm() * this->ucell->lat0
                < ptr_orb_->Phi[T1].getRcut() + ptr_orb_->Alpha[0].getRcut())
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
                // if (std::is_same<TK, double>::value)
                // {
                this->V_delta_R->insert_pair(tmp);
                // }
            }
        }
    }
    // allocate the memory of BaseMatrix in HR, and set the new values to zero
    // if (std::is_same<TK, double>::value)
    // {
    this->V_delta_R->allocate(nullptr, true);
    // expand hR with V_delta_R
    // update : for computational rigor, gamma-only and multi-k cases both have full size of Hamiltonian of DeePKS now
    this->hR->add(*this->V_delta_R);
    this->hR->allocate(nullptr, false);
    // }

    ModuleBase::timer::tick("DeePKS", "initialize_HR");
}
#endif

template <typename TK, typename TR>
void hamilt::DeePKS<hamilt::OperatorLCAO<TK, TR>>::contributeHR()
{
#ifdef __MLALGO
    ModuleBase::TITLE("DeePKS", "contributeHR");
    // if DM changed, HR of DeePKS need to refresh.
    // the judgement is based on the status of HR in ld
    // this operator should be informed that DM has changed and HR need to recalculate.
    if (this->ld->get_hr_cal())
    {
        ModuleBase::timer::tick("DeePKS", "contributeHR");

        const int inlmax = ptr_orb_->Alpha[0].getTotal_nchi() * this->ucell->nat;

        DeePKS_domain::cal_pdm<TK>(this->ld->init_pdm,
                                   inlmax,
                                   this->ld->lmaxd,
                                   this->ld->inl2l,
                                   this->ld->inl_index,
                                   this->kvec_d,
                                   this->ld->dm_r,
                                   this->ld->phialpha,
                                   *this->ucell,
                                   *ptr_orb_,
                                   *(this->gd),
                                   *(this->hR->get_paraV()),
                                   this->ld->pdm);

        std::vector<torch::Tensor> descriptor;
        DeePKS_domain::cal_descriptor(this->ucell->nat,
                                      inlmax,
                                      this->ld->inl2l,
                                      this->ld->pdm,
                                      descriptor,
                                      this->ld->des_per_atom);
        if (PARAM.inp.deepks_equiv)
        {
            DeePKS_domain::cal_edelta_gedm_equiv(this->ucell->nat,
                                                 this->ld->lmaxd,
                                                 this->ld->nmaxd,
                                                 inlmax,
                                                 this->ld->des_per_atom,
                                                 this->ld->inl2l,
                                                 descriptor,
                                                 this->ld->gedm,
                                                 this->ld->E_delta,
                                                 GlobalV::MY_RANK);
        }
        else
        {
            DeePKS_domain::cal_edelta_gedm(this->ucell->nat,
                                           inlmax,
                                           this->ld->des_per_atom,
                                           this->ld->inl2l,
                                           descriptor,
                                           this->ld->pdm,
                                           this->ld->model_deepks,
                                           this->ld->gedm,
                                           this->ld->E_delta);
        }

        // // recalculate the V_delta_R
        // if (this->V_delta_R == nullptr)
        // {
        //     this->V_delta_R = new hamilt::HContainer<std::complex<double>>(*this->hR);
        // }
        this->V_delta_R->set_zero();
        this->calculate_HR();

        this->ld->set_hr_cal(false);

        ModuleBase::timer::tick("DeePKS", "contributeHR");
    }
    // save V_delta_R to hR
    this->hR->add(*this->V_delta_R);
#endif
}

#ifdef __MLALGO

template <typename TK, typename TR>
void hamilt::DeePKS<hamilt::OperatorLCAO<TK, TR>>::calculate_HR()
{
    ModuleBase::TITLE("DeePKS", "calculate_HR");
    ModuleBase::timer::tick("DeePKS", "calculate_HR");
    if (this->V_delta_R->size_atom_pairs() == 0)
    {
        return;
    }

    const Parallel_Orbitals* paraV = this->V_delta_R->get_paraV();
    const int npol = this->ucell->get_npol();

    #pragma omp parallel for schedule(dynamic)
    for (int iat0 = 0; iat0 < this->ucell->nat; iat0++)
    {
        auto tau0 = ucell->get_tau(iat0);
        int T0 = 0;
        int I0 = 0;
        ucell->iat2iait(iat0, &I0, &T0);
        AdjacentAtomInfo& adjs = this->adjs_all[iat0];

        // trace alpha orbital
        std::vector<int> trace_alpha_row;
        std::vector<int> trace_alpha_col;
        std::vector<double> gedms;
        if (!PARAM.inp.deepks_equiv)
        {
            int ib = 0;
            for (int L0 = 0; L0 <= ptr_orb_->Alpha[0].getLmax(); ++L0)
            {
                for (int N0 = 0; N0 < ptr_orb_->Alpha[0].getNchi(L0); ++N0)
                {
                    const int inl = this->ld->inl_index[T0](I0, L0, N0);
                    const double* pgedm = this->ld->gedm[inl];
                    const int nm = 2 * L0 + 1;

                    for (int m1 = 0; m1 < nm; ++m1) // m1 = 1 for s, 3 for p, 5 for d
                    {
                        for (int m2 = 0; m2 < nm; ++m2) // m1 = 1 for s, 3 for p, 5 for d
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
            for (int il = 0; il < this->ld->lmaxd + 1; il++)
            {
                nproj += (2 * il + 1) * ptr_orb_->Alpha[0].getNchi(il);
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
        const int trace_alpha_size = trace_alpha_row.size();

        for (int ad1 = 0; ad1 < adjs.adj_num + 1; ++ad1)
        {
            const int T1 = adjs.ntype[ad1];
            const int I1 = adjs.natom[ad1];
            const int iat1 = ucell->itia2iat(T1, I1);
            ModuleBase::Vector3<int>& R_index1 = adjs.box[ad1];
            auto row_indexes = paraV->get_indexes_row(iat1);
            const int row_size = row_indexes.size();
            if (row_size == 0)
            {
                continue;
            }
            ModuleBase::Vector3<int> dR1(adjs.box[ad1].x, adjs.box[ad1].y, adjs.box[ad1].z);
            if (this->ld->phialpha[0]->find_matrix(iat0, iat1, dR1.x, dR1.y, dR1.z) == nullptr)
            {
                continue;
            }

            std::vector<double> s_1t(trace_alpha_size * row_size);
            for (int irow = 0; irow < row_size; irow++)
            {
                const hamilt::BaseMatrix<double>* overlap_1 = this->ld->phialpha[0]->find_matrix(iat0, iat1, dR1);
                const double* row_ptr = overlap_1->get_pointer() + row_indexes[irow] * overlap_1->get_col_size();
                double* ps1t = &s_1t[irow * trace_alpha_size];
                for (int i = 0; i < trace_alpha_size; i++)
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
                hamilt::BaseMatrix<TR>* tmp
                    = this->V_delta_R->find_matrix(iat1, iat2, R_vector[0], R_vector[1], R_vector[2]);
                // if not found , skip this pair of atoms
                if (tmp == nullptr)
                {
                    continue;
                }
                auto col_indexes = paraV->get_indexes_col(iat2);
                const int col_size = col_indexes.size();

                if (col_size == 0)
                {
                    continue;
                }
                ModuleBase::Vector3<int> dR2(adjs.box[ad2].x, adjs.box[ad2].y, adjs.box[ad2].z);
                if (this->ld->phialpha[0]->find_matrix(iat0, iat2, dR2.x, dR2.y, dR2.z) == nullptr)
                {
                    continue;
                }

                std::vector<double> hr_current(row_size * col_size, 0);
                std::vector<double> s_2t(trace_alpha_size * col_size);
                for (int icol = 0; icol < col_size; icol++)
                {
                    const hamilt::BaseMatrix<double>* overlap_2 = this->ld->phialpha[0]->find_matrix(iat0, iat2, dR2);
                    const double* col_ptr = overlap_2->get_pointer() + col_indexes[icol] * overlap_2->get_col_size();
                    double* ps2t = &s_2t[icol * trace_alpha_size];
                    for (int i = 0; i < trace_alpha_size; i++)
                    {
                        ps2t[i] = col_ptr[trace_alpha_col[i]];
                    }
                }
                // dgemm for s_2t and s_1t to get HR_12
                constexpr char transa = 'T', transb = 'N';
                const double gemm_alpha = 1.0, gemm_beta = 1.0;

                dgemm_(&transa,
                       &transb,
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
            #pragma omp critical
            {
                this->cal_HR_IJR(hr_current.data(), row_size, col_size, tmp->get_pointer());
            }
            }
        }
    }
    ModuleBase::timer::tick("DeePKS", "calculate_HR");
}

// cal_HR_IJR()
template <typename TK, typename TR>
void hamilt::DeePKS<hamilt::OperatorLCAO<TK, TR>>::cal_HR_IJR(const double* hr_in,
                                                              const int& row_size,
                                                              const int& col_size,
                                                              TR* data_pointer)
{

    //! npol is the number of polarizations,
    //! 1 for non-magnetic (one Hamiltonian matrix only has spin-up or spin-down),
    //! 2 for magnetic (one Hamiltonian matrix has both spin-up and spin-down)
    const int npol = this->ucell->get_npol();

    // step_trace = 0 for NSPIN=1,2; ={0, 1, local_col, local_col+1} for NSPIN=4
    vector<int> step_trace(2, 0);
    step_trace[1] = col_size + 1;

    //! calculate the local matrix
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

// contributeHk()
template <typename TK, typename TR>
void hamilt::DeePKS<hamilt::OperatorLCAO<TK, TR>>::contributeHk(int ik)
{
    ModuleBase::TITLE("DeePKS", "contributeHk");
    ModuleBase::timer::tick("DeePKS", "contributeHk");

    TK* h_delta_k = this->ld->V_delta[ik].data();
    // set SK to zero and then calculate SK for each k vector
    ModuleBase::GlobalFunc::ZEROS(h_delta_k, this->hsk->get_size());

    if (ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER(PARAM.inp.ks_solver))
    {
        const int nrow = this->hsk->get_pv()->get_row_size();
        hamilt::folding_HR(*this->V_delta_R, h_delta_k, this->kvec_d[ik], nrow, 1);
    }
    else
    {
        const int ncol = this->hsk->get_pv()->get_col_size();
        hamilt::folding_HR(*this->V_delta_R, h_delta_k, this->kvec_d[ik], ncol, 0);
    }
    ModuleBase::timer::tick("DeePKS", "contributeHk");
}

#endif

template class DeePKS<OperatorLCAO<double, double>>;
template class DeePKS<OperatorLCAO<std::complex<double>, double>>;
template class DeePKS<OperatorLCAO<std::complex<double>, std::complex<double>>>;

} // namespace hamilt
