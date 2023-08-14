#include "nonlocal_new.h"

#include "module_basis/module_ao/ORB_gen_tables.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_hamilt_lcao/hamilt_lcaodft/operator_lcao/operator_lcao.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer_funcs.h"

template <typename TK, typename TR>
hamilt::NonlocalNew<hamilt::OperatorLCAO<TK>, TR>::NonlocalNew(
    LCAO_Matrix* LM_in,
    const std::vector<ModuleBase::Vector3<double>>& kvec_d_in,
    hamilt::HContainer<TR>* HR_in,
    TK* HK_pointer_in,
    const UnitCell* ucell_in,
    Grid_Driver* GridD_in,
    const Parallel_Orbitals* paraV)
    : hamilt::OperatorLCAO<TK>(LM_in, kvec_d_in)
{
    this->ucell = ucell_in;
    this->HR = HR_in;
    this->HK_pointer = HK_pointer_in;
#ifdef __DEBUG
    assert(this->ucell != nullptr);
    assert(this->HR != nullptr);
    assert(this->HK_pointer != nullptr);
#endif
    // initialize HR to allocate sparse Nonlocal matrix memory
    this->initialize_HR(GridD_in, paraV);
}

// initialize_HR()
template <typename TK, typename TR>
void hamilt::NonlocalNew<hamilt::OperatorLCAO<TK>, TR>::initialize_HR(Grid_Driver* GridD,
                                                                      const Parallel_Orbitals* paraV)
{
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
            if(this->ucell->cal_dtau(iat0, iat1, R_index1).norm() * this->ucell->lat0 <= orb.Phi[T1].getRcut() + this->ucell->infoNL.Beta[T0].get_rcut_max())
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
                if(paraV->get_col_size(iat2)<=0 || paraV->get_row_size(iat1) <= 0)
                {
                    continue;
                }
                hamilt::AtomPair<TR> tmp(iat1,
                                         iat2,
                                         R_index2.x - R_index1.x,
                                         R_index2.y - R_index1.y,
                                         R_index2.z - R_index1.z,
                                         paraV);
                HR->insert_pair(tmp);
            }
        }
    }
    HR->allocate(1);
}

template <typename TK, typename TR>
void hamilt::NonlocalNew<hamilt::OperatorLCAO<TK>, TR>::calculate_HR()
{
    const Parallel_Orbitals* paraV = this->HR->get_atom_pair(0).get_paraV();
    const int npol = this->ucell->get_npol();
    // 1. calculate <psi|beta> for each pair of atoms
    for (int iat0 = 0; iat0 < this->ucell->nat; iat0++)
    {
        auto tau0 = ucell->get_tau(iat0);
        int T0, I0;
        ucell->iat2iait(iat0, &I0, &T0);
        AdjacentAtomInfo& adjs = this->adjs_all[iat0];

        std::vector<std::unordered_map<int, std::vector<double>>> nlm_tot;
        nlm_tot.resize(adjs.adj_num + 1);

#ifdef _OPENMP
#pragma omp parallel for
#endif
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
                uot.snap_psibeta_half(orb,
                                      this->ucell->infoNL,
                                      nlm,
                                      tau1,
                                      T1,
                                      atom1->iw2l[iw1], // L1
                                      atom1->iw2m[iw1], // m1
                                      atom1->iw2n[iw1], // N1
                                      tau0,
                                      T0,
                                      0 /*cal_deri*/); // R0,T0
                nlm_tot[ad].insert({all_indexes[iw1l], nlm[0]});
            }
        }
// 2. calculate <psi_I|beta>D<beta|psi_{J,R}> for each pair of <IJR> atoms
#ifdef _OPENMP
#pragma omp parallel for
#endif
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
                hamilt::BaseMatrix<TR>* tmp = this->HR->find_matrix(iat1, iat2, R_vector[0], R_vector[1], R_vector[2]);
                if (tmp != nullptr)
                {
                    this->cal_HR_IJR(iat1, iat2, T0, paraV, nlm_tot[ad1], nlm_tot[ad2], tmp->get_pointer());
                }
            }
        }
    }
}

// cal_HR_IJR()
template <typename TK, typename TR>
void hamilt::NonlocalNew<hamilt::OperatorLCAO<TK>, TR>::cal_HR_IJR(const int& iat1,
                                                                   const int& iat2,
                                                                   const int& T0,
                                                                   const Parallel_Orbitals* paraV,
                                                                   const std::unordered_map<int, std::vector<double>> & nlm1_all,
                                                                   const std::unordered_map<int, std::vector<double>> & nlm2_all,
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
    std::vector<int> step_trace(npol * npol, 0);
    for (int is = 0; is < npol; is++)
    {
        for (int is2 = 0; is2 < npol; is2++)
        {
            step_trace[is * npol + is2] = col_indexes.size() * is + is2;
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
            for (int is = 0; is < npol * npol; ++is)
            {
                TR nlm_tmp = TR(0);
                for (int no = 0; no < this->ucell->atoms[T0].ncpp.non_zero_count_soc[is]; no++)
                {
                    const int p1 = this->ucell->atoms[T0].ncpp.index1_soc[is][no];
                    const int p2 = this->ucell->atoms[T0].ncpp.index2_soc[is][no];
                    this->ucell->atoms[T0].ncpp.get_d(is, p2, p1, tmp_d);
                    nlm_tmp += nlm1[p1] * nlm2[p2] * (*tmp_d);
                }
                data_pointer[step_trace[is]] += nlm_tmp;
            }
            data_pointer += npol;
        }
        data_pointer += (npol-1) * col_indexes.size();
    }
}

// contributeHR()
template <typename TK, typename TR>
void hamilt::NonlocalNew<hamilt::OperatorLCAO<TK>, TR>::contributeHR()
{
    if (this->HR_fixed_done)
    {
        return;
    }
    this->calculate_HR();
    this->HR_fixed_done = true;
}

// contributeHk()
template <typename TK, typename TR>
void hamilt::NonlocalNew<hamilt::OperatorLCAO<TK>, TR>::contributeHk(int ik)
{
    const int ncol = this->HR->get_atom_pair(0).get_paraV()->get_col_size();
    hamilt::folding_HR(*this->HR, this->HK_pointer, this->kvec_d[ik], ncol, 0);
}

template class hamilt::NonlocalNew<hamilt::OperatorLCAO<double>, double>;
template class hamilt::NonlocalNew<hamilt::OperatorLCAO<std::complex<double>>, double>;
template class hamilt::NonlocalNew<hamilt::OperatorLCAO<std::complex<double>>, std::complex<double>>;