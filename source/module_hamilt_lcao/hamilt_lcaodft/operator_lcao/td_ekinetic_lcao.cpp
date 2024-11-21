#include "td_ekinetic_lcao.h"

#include "module_parameter/parameter.h"
#include "module_base/global_variable.h"
#include "module_base/libm/libm.h"
#include "module_base/timer.h"
#include "module_base/tool_title.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_elecstate/potentials/H_TDDFT_pw.h"
#include "module_hamilt_lcao/hamilt_lcaodft/center2_orb-orb11.h"
#include "module_hamilt_lcao/hamilt_lcaodft/spar_hsr.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer_funcs.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

namespace hamilt
{
template <typename TK, typename TR>
TDEkinetic<OperatorLCAO<TK, TR>>::TDEkinetic(HS_Matrix_K<TK>* hsk_in,
                                             hamilt::HContainer<TR>* hR_in,
                                             const K_Vectors* kv_in,
                                             const UnitCell* ucell_in,
                                             const std::vector<double>& orb_cutoff,
                                             Grid_Driver* GridD_in,
                                             const TwoCenterIntegrator* intor)
    : OperatorLCAO<TK, TR>(hsk_in, kv_in->kvec_d, hR_in), orb_cutoff_(orb_cutoff), kv(kv_in), intor_(intor)
{
    this->ucell = ucell_in;
    this->cal_type = calculation_type::lcao_tddft_velocity;
    this->Grid = GridD_in;
    this->init_td();
    // initialize HR to get adjs info.
    this->initialize_HR(Grid);
}
template <typename TK, typename TR>
TDEkinetic<OperatorLCAO<TK, TR>>::~TDEkinetic()
{
    if (this->hR_tmp != nullptr)
    {
        delete this->hR_tmp;
    }
    TD_Velocity::td_vel_op = nullptr;
}
// term A^2*S
template <typename TK, typename TR>
void TDEkinetic<OperatorLCAO<TK, TR>>::td_ekinetic_scalar(std::complex<double>* Hloc,const TR& overlap, int nnr)
{
    return;
}
// term A^2*S
template <>
void TDEkinetic<OperatorLCAO<std::complex<double>, double>>::td_ekinetic_scalar(std::complex<double>* Hloc,
                                                                                const double& overlap,
                                                                                int nnr)
{
    // the correction term A^2/2. From Hatree to Ry, it needs to be multiplied by 2.0
    std::complex<double> tmp = {cart_At.norm2() * overlap, 0};
    Hloc[nnr] += tmp;
    return;
}
// term A dot ∇
template <typename TK, typename TR>
void TDEkinetic<OperatorLCAO<TK, TR>>::td_ekinetic_grad(std::complex<double>* Hloc,
                                                        int nnr,
                                                        ModuleBase::Vector3<double> grad_overlap)
{
    // the correction term -iA dot ∇r
    //∇ refer to the integral ∫𝜙(𝑟)𝜕/𝜕𝑟𝜙(𝑟−𝑅)𝑑𝑟,but abacus only provide the integral of ∫𝜙(𝑟)𝜕/𝜕R𝜙(𝑟−𝑅)𝑑𝑟. An extra
    //minus must be counted in. The final term is iA dot ∇R. From Hatree to Ry, it needs to be multiplied by 2.0
    std::complex<double> tmp = {0, grad_overlap * cart_At};
    Hloc[nnr] += tmp * 2.0;
    return;
}

template <typename TK, typename TR>
void TDEkinetic<OperatorLCAO<TK, TR>>::calculate_HR()
{
    ModuleBase::TITLE("TDEkinetic", "calculate_HR");
    if (this->hR_tmp == nullptr || this->hR_tmp->size_atom_pairs() <= 0)
    {
        ModuleBase::WARNING_QUIT("TDEkinetic::calculate_HR", "hR_tmp is nullptr or empty");
    }
    ModuleBase::timer::tick("TDEkinetic", "calculate_HR");

    const Parallel_Orbitals* paraV = this->hR_tmp->get_atom_pair(0).get_paraV();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int iat1 = 0; iat1 < this->ucell->nat; iat1++)
    {
        auto tau1 = ucell->get_tau(iat1);
        int T1, I1;
        ucell->iat2iait(iat1, &I1, &T1);
        AdjacentAtomInfo& adjs = this->adjs_all[iat1];
        for (int ad = 0; ad < adjs.adj_num + 1; ++ad)
        {
            const int T2 = adjs.ntype[ad];
            const int I2 = adjs.natom[ad];
            const int iat2 = ucell->itia2iat(T2, I2);
            const ModuleBase::Vector3<int>& R_index2 = adjs.box[ad];
            ModuleBase::Vector3<double> dtau = this->ucell->cal_dtau(iat1, iat2, R_index2);

            hamilt::BaseMatrix<std::complex<double>>* tmp = this->hR_tmp->find_matrix(iat1, iat2, R_index2);
            if (tmp != nullptr)
            {
                if (TD_Velocity::out_current)
                {
                    std::complex<double>* tmp_c[3] = {nullptr, nullptr, nullptr};
                    for (int i = 0; i < 3; i++)
                    {
                        tmp_c[i] = td_velocity.get_current_term_pointer(i)->find_matrix(iat1, iat2, R_index2)->get_pointer();
                    }
                    this->cal_HR_IJR(iat1, iat2, paraV, dtau, tmp->get_pointer(), tmp_c);
                }
                else
                {
                    this->cal_HR_IJR(iat1, iat2, paraV, dtau, tmp->get_pointer(), nullptr);
                }
            }
            else
            {
                ModuleBase::WARNING_QUIT("TDEkinetic::calculate_HR", "R_index not found in HR");
            }
        }
    }
    ModuleBase::timer::tick("TDEkinetic", "calculate_HR");
}

template <typename TK, typename TR>
void TDEkinetic<OperatorLCAO<TK, TR>>::cal_HR_IJR(const int& iat1,
                                                  const int& iat2,
                                                  const Parallel_Orbitals* paraV,
                                                  const ModuleBase::Vector3<double>& dtau,
                                                  std::complex<double>* hr_mat_p,
                                                  std::complex<double>** current_mat_p)
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

    // npol is the number of polarizations,
    // 1 for non-magnetic (one Hamiltonian matrix only has spin-up or spin-down),
    // 2 for magnetic (one Hamiltonian matrix has both spin-up and spin-down)
    const int npol = this->ucell->get_npol();

    const int* iw2l1 = atom1.iw2l.data();
    const int* iw2n1 = atom1.iw2n.data();
    const int* iw2m1 = atom1.iw2m.data();
    const int* iw2l2 = atom2.iw2l.data();
    const int* iw2n2 = atom2.iw2n.data();
    const int* iw2m2 = atom2.iw2m.data();
    // ---------------------------------------------
    // get tau1 (in cell <0,0,0>) and tau2 (in cell R)
    // in principle, only dtau is needed in this function
    // snap_psipsi should be refactored to use dtau directly
    // ---------------------------------------------
    const ModuleBase::Vector3<double>& tau1 = this->ucell->get_tau(iat1);
    const ModuleBase::Vector3<double> tau2 = tau1 + dtau;
    // ---------------------------------------------
    // calculate the Ekinetic matrix for each pair of orbitals
    // ---------------------------------------------
    double grad[3] = {0, 0, 0};
    double overlap = 0;
    auto row_indexes = paraV->get_indexes_row(iat1);
    auto col_indexes = paraV->get_indexes_col(iat2);
    const int step_trace = col_indexes.size() + 1;
    for (int iw1l = 0; iw1l < row_indexes.size(); iw1l += npol)
    {
        const int iw1 = row_indexes[iw1l] / npol;
        const int L1 = iw2l1[iw1];
        const int N1 = iw2n1[iw1];
        const int m1 = iw2m1[iw1];

        // convert m (0,1,...2l) to M (-l, -l+1, ..., l-1, l)
        int M1 = (m1 % 2 == 0) ? -m1 / 2 : (m1 + 1) / 2;

        for (int iw2l = 0; iw2l < col_indexes.size(); iw2l += npol)
        {
            const int iw2 = col_indexes[iw2l] / npol;
            const int L2 = iw2l2[iw2];
            const int N2 = iw2n2[iw2];
            const int m2 = iw2m2[iw2];

            // convert m (0,1,...2l) to M (-l, -l+1, ..., l-1, l)
            int M2 = (m2 % 2 == 0) ? -m2 / 2 : (m2 + 1) / 2;

            // calculate <psi|∇R|psi>, which equals to -<psi|∇r|psi>.
            intor_->calculate(T1, L1, N1, M1, T2, L2, N2, M2, dtau * this->ucell->lat0, &overlap, grad);
            ModuleBase::Vector3<double> grad_overlap(grad[0], grad[1], grad[2]);

            for (int ipol = 0; ipol < npol; ipol++)
            {
                // key change
                td_ekinetic_scalar(hr_mat_p, overlap, ipol * step_trace);
                td_ekinetic_grad(hr_mat_p, ipol * step_trace, grad_overlap);
            }
            hr_mat_p += npol;
            // current grad part
            if (current_mat_p != nullptr)
            {
                for (int dir = 0; dir < 3; dir++)
                {
                    for (int ipol = 0; ipol < npol; ipol++)
                    {
                        // part of Momentum operator, -i∇r,used to calculate the current
                        // here is actually i∇R
                        current_mat_p[dir][ipol * step_trace] += std::complex<double>(0, grad_overlap[dir]);
                        // part of Momentum operator, eA,used to calculate the current
                        current_mat_p[dir][ipol * step_trace] += std::complex<double>(overlap * cart_At[dir], 0);
                    }
                    current_mat_p[dir] += npol;
                }
            }
        }
        hr_mat_p += (npol - 1) * col_indexes.size();
        if (current_mat_p != nullptr)
        {
            for (int dir = 0; dir < 3; dir++) {
                current_mat_p[dir] += (npol - 1) * col_indexes.size();
}
        }
    }
}
// init two center integrals and vector potential for td_ekintic term
template <typename TK, typename TR>
void TDEkinetic<OperatorLCAO<TK, TR>>::init_td()
{
    TD_Velocity::td_vel_op = &td_velocity;
    // calculate At in cartesian coorinates.
    td_velocity.cal_cart_At(elecstate::H_TDDFT_pw::At);
    this->cart_At = td_velocity.cart_At;
    std::cout << "cart_At: " << cart_At[0] << " " << cart_At[1] << " " << cart_At[2] << std::endl;
}

template <typename TK, typename TR>
void hamilt::TDEkinetic<hamilt::OperatorLCAO<TK, TR>>::set_HR_fixed(void* hR_tmp_in)
{
    this->hR_tmp = static_cast<hamilt::HContainer<std::complex<double>>*>(hR_tmp_in);
    this->allocated = false;
}
template <typename TK, typename TR>
void TDEkinetic<OperatorLCAO<TK, TR>>::initialize_HR(Grid_Driver* GridD)
{
    if (elecstate::H_TDDFT_pw::stype != 1)
    {
        return;
    }
    ModuleBase::TITLE("TDEkinetic", "initialize_HR");
    ModuleBase::timer::tick("TDEkinetic", "initialize_HR");

    auto* paraV = this->hR->get_paraV();// get parallel orbitals from HR
    // TODO: if paraV is nullptr, AtomPair can not use paraV for constructor, I will repair it in the future.

    this->adjs_all.clear();
    this->adjs_all.reserve(this->ucell->nat);
    for (int iat1 = 0; iat1 < ucell->nat; iat1++)
    {
        auto tau1 = ucell->get_tau(iat1);
        int T1, I1;
        ucell->iat2iait(iat1, &I1, &T1);
        AdjacentAtomInfo adjs;
        GridD->Find_atom(*ucell, tau1, T1, I1, &adjs);
        std::vector<bool> is_adj(adjs.adj_num + 1, false);
        for (int ad1 = 0; ad1 < adjs.adj_num + 1; ++ad1)
        {
            const int T2 = adjs.ntype[ad1];
            const int I2 = adjs.natom[ad1];
            const int iat2 = ucell->itia2iat(T2, I2);
            if (paraV->get_row_size(iat1) <= 0 || paraV->get_col_size(iat2) <= 0)
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
        this->adjs_all.push_back(adjs);
    }
    ModuleBase::timer::tick("TDEkinetic", "initialize_HR");
}
template <typename TK, typename TR>
void TDEkinetic<OperatorLCAO<TK, TR>>::initialize_HR_tmp()
{
    if (elecstate::H_TDDFT_pw::stype != 1)
    {
        return;
    }
    ModuleBase::TITLE("TDEkinetic", "initialize_HR_tmp");
    ModuleBase::timer::tick("TDEkinetic", "initialize_HR_tmp");

    auto* paraV = this->hR->get_paraV();// get parallel orbitals from HR
    // TODO: if paraV is nullptr, AtomPair can not use paraV for constructor, I will repair it in the future.
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

    ModuleBase::timer::tick("TDEkinetic", "initialize_HR_tmp");
}

template <typename TK, typename TR>
void TDEkinetic<OperatorLCAO<TK, TR>>::contributeHR()
{
    // const Parallel_Orbitals* paraV = this->hR->get_atom_pair(0).get_paraV();
    ModuleBase::TITLE("TDEkinetic", "contributeHR");
    ModuleBase::timer::tick("TDEkinetic", "contributeHR");
    // skip if not TDDFT velocity gauge
    if (elecstate::H_TDDFT_pw::stype != 1)
    {
        return;
    }
    if (!this->hR_tmp_done)
    {
        const Parallel_Orbitals* paraV = this->hR->get_atom_pair(0).get_paraV();
        // if this Operator is the first node of the sub_chain, then hR_tmp is nullptr
        if (this->hR_tmp == nullptr)
        {
            this->hR_tmp = new hamilt::HContainer<std::complex<double>>(this->hR->get_paraV());
            // allocate memory for hR_tmp use the same memory as hR
            this->initialize_HR_tmp();
            this->allocated = true;
        }
        if (this->next_sub_op != nullptr)
        {
            // pass pointer of hR_tmp to the next node
            static_cast<OperatorLCAO<TK, TR>*>(this->next_sub_op)->set_HR_fixed(this->hR_tmp);
        }
        // initialize current term if needed
        if (TD_Velocity::out_current)
        {
            td_velocity.initialize_current_term(this->hR_tmp, paraV);
        }
        // calculate the values in hR_tmp
        this->calculate_HR();
        this->hR_tmp_done = true;
    }

    ModuleBase::timer::tick("TDEkinetic", "contributeHR");
    return;
}

template <typename TK, typename TR>
void TDEkinetic<OperatorLCAO<TK, TR>>::contributeHk(int ik)
{
    return;
}
template <>
void TDEkinetic<OperatorLCAO<std::complex<double>, double>>::contributeHk(int ik)
{
    if (TD_Velocity::tddft_velocity == false)
    {
        return;
    }
    else
    {
        ModuleBase::TITLE("TDEkinetic", "contributeHk");
        ModuleBase::timer::tick("TDEkinetic", "contributeHk");
        const Parallel_Orbitals* paraV = this->hR_tmp->get_atom_pair(0).get_paraV();
        // save HR data for output
        int spin_tot = PARAM.inp.nspin;
        if (spin_tot == 4)
            ;
        else if (!output_hR_done && TD_Velocity::out_mat_R)
        {
            for (int spin_now = 0; spin_now < spin_tot; spin_now++)
            {
                sparse_format::cal_HContainer_cd(*(paraV),
                                                 spin_now,
                                                 1e-10,
                                                 *hR_tmp,
                                                 td_velocity.HR_sparse_td_vel[spin_now]);
            }
            output_hR_done = true;
        }
        // folding inside HR to HK
        if (ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER(PARAM.inp.ks_solver))
        {
            const int nrow = paraV->get_row_size();
            hamilt::folding_HR(*this->hR_tmp, this->hsk->get_hk(), this->kvec_d[ik], nrow, 1);
        }
        else
        {
            const int ncol = paraV->get_col_size();
            hamilt::folding_HR(*this->hR_tmp, this->hsk->get_hk(), this->kvec_d[ik], ncol, 0);
        }

        ModuleBase::timer::tick("TDEkinetic", "contributeHk");
    }
}

template class TDEkinetic<hamilt::OperatorLCAO<double, double>>;
template class TDEkinetic<hamilt::OperatorLCAO<std::complex<double>, double>>;
template class TDEkinetic<hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>>;

} // namespace hamilt
