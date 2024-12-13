#include "td_nonlocal_lcao.h"

#include "module_parameter/parameter.h"
#include "module_base/timer.h"
#include "module_base/tool_title.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_hamilt_lcao/hamilt_lcaodft/operator_lcao/operator_lcao.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer_funcs.h"
#include "module_hamilt_lcao/module_tddft/snap_psibeta_half_tddft.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#ifdef _OPENMP
#include <unordered_set>
#include <omp.h>
#endif

template <typename TK, typename TR>
hamilt::TDNonlocal<hamilt::OperatorLCAO<TK, TR>>::TDNonlocal(HS_Matrix_K<TK>* hsk_in,
                                                             const std::vector<ModuleBase::Vector3<double>>& kvec_d_in,
                                                             hamilt::HContainer<TR>* hR_in,
                                                             const UnitCell* ucell_in,
                                                             const LCAO_Orbitals& orb,
                                                             const Grid_Driver* GridD_in)
    : hamilt::OperatorLCAO<TK, TR>(hsk_in, kvec_d_in, hR_in), orb_(orb)
{
    this->cal_type = calculation_type::lcao_tddft_velocity;
    this->ucell = ucell_in;
    this->Grid = GridD_in;
#ifdef __DEBUG
    assert(this->ucell != nullptr);
#endif
    // initialize HR to get adjs info.
    this->init_td();
    this->initialize_HR(Grid);
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
void hamilt::TDNonlocal<hamilt::OperatorLCAO<TK, TR>>::init_td()
{
    // calculate At in cartesian coorinates.
    this->cart_At = TD_Velocity::td_vel_op->cart_At;
}
// initialize_HR()
template <typename TK, typename TR>
void hamilt::TDNonlocal<hamilt::OperatorLCAO<TK, TR>>::initialize_HR(const Grid_Driver* GridD)
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
            // Note: the distance of atoms should less than the cutoff radius,
            // When equal, the theoretical value of matrix element is zero,
            // but the calculated value is not zero due to the numerical error, which would lead to result changes.
            if (this->ucell->cal_dtau(iat0, iat1, R_index1).norm() * this->ucell->lat0
                < orb_.Phi[T1].getRcut() + this->ucell->infoNL.Beta[T0].get_rcut_max())
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
        for (int ir = 0; ir < tmp.get_R_size(); ++ir)
        {
            const ModuleBase::Vector3<int> R_index = tmp.get_R_index(ir);
            const int iat1 = tmp.get_atom_i();
            const int iat2 = tmp.get_atom_j();

            hamilt::AtomPair<std::complex<double>> tmp1(iat1, iat2, R_index, paraV);
            this->hR_tmp->insert_pair(tmp1);
        }
    }
    this->hR_tmp->allocate(nullptr, true);

    ModuleBase::timer::tick("TDNonlocal", "initialize_HR_tmp");
}

template <typename TK, typename TR>
void hamilt::TDNonlocal<hamilt::OperatorLCAO<TK, TR>>::calculate_HR()
{
    ModuleBase::TITLE("TDNonlocal", "calculate_HR");
    ModuleBase::timer::tick("TDNonlocal", "calculate_HR");

    const Parallel_Orbitals* paraV = this->hR_tmp->get_atom_pair(0).get_paraV();
    const int npol = this->ucell->get_npol();
    const int nlm_dim = TD_Velocity::out_current ? 4 : 1;
    // 1. calculate <psi|beta> for each pair of atoms

    for (int iat0 = 0; iat0 < this->ucell->nat; iat0++)
    {
        const auto tau0 = ucell->get_tau(iat0);
        int T0, I0;
        ucell->iat2iait(iat0, &I0, &T0);
        const AdjacentAtomInfo& adjs = this->adjs_all[iat0];
        std::vector<std::vector<std::unordered_map<int, std::vector<std::complex<double>>>>> nlm_tot;
        nlm_tot.resize(adjs.adj_num + 1);
        for (int i = 0; i < adjs.adj_num + 1; i++)
        {
            nlm_tot[i].resize(nlm_dim);
        }

        #pragma omp parallel
        {
            #pragma omp for schedule(dynamic)
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
                    std::vector<std::vector<std::complex<double>>> nlm;
                    // nlm is a vector of vectors, but size of outer vector is only 1 when out_current is false
                    // and size of outer vector is 4 when out_current is true (3 for <psi|r_i * exp(-iAr)|beta>, 1 for
                    // <psi|exp(-iAr)|beta>) inner loop : all projectors (L0,M0)

                    // snap_psibeta_half_tddft() are used to calculate <psi|exp(-iAr)|beta>
                    // and <psi|rexp(-iAr)|beta> as well if current are needed
                    module_tddft::snap_psibeta_half_tddft(orb_,
                                                          this->ucell->infoNL,
                                                          nlm,
                                                          tau1 * this->ucell->lat0,
                                                          T1,
                                                          atom1->iw2l[iw1],
                                                          atom1->iw2m[iw1],
                                                          atom1->iw2n[iw1],
                                                          tau0 * this->ucell->lat0,
                                                          T0,
                                                          cart_At,
                                                          TD_Velocity::out_current);
                    for (int dir = 0; dir < nlm_dim; dir++)
                    {
                        nlm_tot[ad][dir].insert({all_indexes[iw1l], nlm[dir]});
                    }
                }
            }

#ifdef _OPENMP
            // record the iat number of the adjacent atoms
            std::set<int> ad_atom_set;
            for (int ad = 0; ad < adjs.adj_num + 1; ++ad)
            {
                const int T1 = adjs.ntype[ad];
                const int I1 = adjs.natom[ad];
                const int iat1 = ucell->itia2iat(T1, I1);
                ad_atom_set.insert(iat1);
            }

            // split the ad_atom_set into num_threads parts
            const int num_threads = omp_get_num_threads();
            const int thread_id = omp_get_thread_num();
            std::set<int> ad_atom_set_thread;
            int i = 0;
            for(const auto iat1 : ad_atom_set)
            {
                if (i % num_threads == thread_id)
                {
                    ad_atom_set_thread.insert(iat1);
                }
                i++;
            }
#endif

            // 2. calculate <psi_I|beta>D<beta|psi_{J,R}> for each pair of <IJR> atoms
            for (int ad1 = 0; ad1 < adjs.adj_num + 1; ++ad1)
            {
                const int T1 = adjs.ntype[ad1];
                const int I1 = adjs.natom[ad1];
                const int iat1 = ucell->itia2iat(T1, I1);

#ifdef _OPENMP
                if (ad_atom_set_thread.find(iat1) == ad_atom_set_thread.end())
                {
                    continue;
                }
#endif
                
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
                    hamilt::BaseMatrix<std::complex<double>>* tmp
                        = this->hR_tmp->find_matrix(iat1, iat2, R_vector[0], R_vector[1], R_vector[2]);
                    // if not found , skip this pair of atoms
                    if (tmp != nullptr)
                    {
                        if (TD_Velocity::out_current)
                        {
                            std::complex<double>* tmp_c[3] = {nullptr, nullptr, nullptr};
                            for (int i = 0; i < 3; i++)
                            {
                                tmp_c[i] = TD_Velocity::td_vel_op->get_current_term_pointer(i)
                                                ->find_matrix(iat1, iat2, R_vector[0], R_vector[1], R_vector[2])
                                                ->get_pointer();
                            }
                            this->cal_HR_IJR(iat1,
                                             iat2,
                                             T0,
                                             paraV,
                                             nlm_tot[ad1],
                                             nlm_tot[ad2],
                                             tmp->get_pointer(),
                                             tmp_c);
                        }
                        else
                        {
                            this->cal_HR_IJR(iat1,
                                             iat2,
                                             T0,
                                             paraV,
                                             nlm_tot[ad1],
                                             nlm_tot[ad2],
                                             tmp->get_pointer(),
                                             nullptr);
                        }
                    }
                }
            }
        }
    }

    ModuleBase::timer::tick("TDNonlocal", "calculate_HR");
}

// cal_HR_IJR()
template <typename TK, typename TR>
void hamilt::TDNonlocal<hamilt::OperatorLCAO<TK, TR>>::cal_HR_IJR(
    const int& iat1,
    const int& iat2,
    const int& T0,
    const Parallel_Orbitals* paraV,
    const std::vector<std::unordered_map<int, std::vector<std::complex<double>>>>& nlm1_all,
    const std::vector<std::unordered_map<int, std::vector<std::complex<double>>>>& nlm2_all,
    std::complex<double>* data_pointer,
    std::complex<double>** data_pointer_c)
{
    const int nlm_dim = TD_Velocity::out_current ? 4 : 1;
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
    const std::complex<double>* tmp_d = nullptr;
    for (int iw1l = 0; iw1l < row_indexes.size(); iw1l += npol)
    {
        // const std::vector<std::complex<double>>* nlm1 = &(nlm1_all[0].find(row_indexes[iw1l])->second);
        std::vector<const std::vector<std::complex<double>>*> nlm1;
        for (int dir = 0; dir < nlm_dim; dir++)
        {
            nlm1.push_back(&(nlm1_all[dir].find(row_indexes[iw1l])->second));
        }

        for (int iw2l = 0; iw2l < col_indexes.size(); iw2l += npol)
        {
            std::vector<const std::vector<std::complex<double>>*> nlm2;
            for (int dir = 0; dir < nlm_dim; dir++)
            {
                nlm2.push_back(&(nlm2_all[dir].find(col_indexes[iw2l])->second));
            }
#ifdef __DEBUG
            assert(nlm1.size() == nlm2.size());
#endif
            for (int is = 0; is < npol * npol; ++is)
            {
                std::complex<double> nlm_tmp = std::complex<double>{0, 0};
                for (int no = 0; no < this->ucell->atoms[T0].ncpp.non_zero_count_soc[is]; no++)
                {
                    const int p1 = this->ucell->atoms[T0].ncpp.index1_soc[is][no];
                    const int p2 = this->ucell->atoms[T0].ncpp.index2_soc[is][no];
                    this->ucell->atoms[T0].ncpp.get_d(is, p1, p2, tmp_d);
                    nlm_tmp += nlm1[0]->at(p1) * std::conj(nlm2[0]->at(p2)) * (*tmp_d);
                }
                data_pointer[step_trace[is]] += nlm_tmp;
                if (data_pointer_c != nullptr)
                {
                    for (int dir = 0; dir < 3; dir++)
                    {
                        std::complex<double> nlm_r_tmp = std::complex<double>{0, 0};
                        std::complex<double> imag_unit = std::complex<double>{0, 1};
                        for (int no = 0; no < this->ucell->atoms[T0].ncpp.non_zero_count_soc[is]; no++)
                        {
                            const int p1 = this->ucell->atoms[T0].ncpp.index1_soc[is][no];
                            const int p2 = this->ucell->atoms[T0].ncpp.index2_soc[is][no];
                            this->ucell->atoms[T0].ncpp.get_d(is, p1, p2, tmp_d);
                            //<psi|rexp(-iAr)|beta><beta|exp(iAr)|psi>-<psi|exp(-iAr)|beta><beta|rexp(iAr)|psi>
                            // multiply d in the end
                            nlm_r_tmp += (nlm1[dir + 1]->at(p1) * std::conj(nlm2[0]->at(p2))
                                          - nlm1[0]->at(p1) * std::conj(nlm2[dir + 1]->at(p2)))
                                         * (*tmp_d);
                        }
                        // -i[r,Vnl], 2.0 due to the unit transformation
                        data_pointer_c[dir][step_trace[is]] -= imag_unit * nlm_r_tmp / 2.0;
                    }
                }
            }
            data_pointer += npol;
            if (data_pointer_c != nullptr)
            {
                for (int dir = 0; dir < 3; dir++)
                {
                    data_pointer_c[dir] += npol;
                }
            }
        }
        data_pointer += (npol - 1) * col_indexes.size();
        if (data_pointer_c != nullptr)
        {
            for (int dir = 0; dir < 3; dir++)
            {
                data_pointer_c[dir] += (npol - 1) * col_indexes.size();
            }
        }
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
            this->hR_tmp = new hamilt::HContainer<std::complex<double>>(this->hsk->get_pv());
            // allocate memory for hR_tmp use the same memory as hR
            this->initialize_HR_tmp(this->hsk->get_pv());
            this->allocated = true;
        }
        if (this->next_sub_op != nullptr)
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
template <typename TK, typename TR>
void hamilt::TDNonlocal<hamilt::OperatorLCAO<TK, TR>>::contributeHk(int ik)
{
    return;
}
template <>
void hamilt::TDNonlocal<hamilt::OperatorLCAO<std::complex<double>, double>>::contributeHk(int ik)
{
    if (TD_Velocity::tddft_velocity == false)
    {
        return;
    }
    else
    {
        ModuleBase::TITLE("TDNonlocal", "contributeHk");
        ModuleBase::timer::tick("TDNonlocal", "contributeHk");
        // folding inside HR to HK
        if (ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER(PARAM.inp.ks_solver))
        {
            const int nrow = this->hsk->get_pv()->get_row_size();
            folding_HR(*this->hR_tmp, this->hsk->get_hk(), this->kvec_d[ik], nrow, 1);
        }
        else
        {
            const int ncol = this->hsk->get_pv()->get_col_size();
            folding_HR(*this->hR_tmp, this->hsk->get_hk(), this->kvec_d[ik], ncol, 0);
        }

        ModuleBase::timer::tick("TDNonlocal", "contributeHk");
    }
}
template class hamilt::TDNonlocal<hamilt::OperatorLCAO<double, double>>;
template class hamilt::TDNonlocal<hamilt::OperatorLCAO<std::complex<double>, double>>;
template class hamilt::TDNonlocal<hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>>;
