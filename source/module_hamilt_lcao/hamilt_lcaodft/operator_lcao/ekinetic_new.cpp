#include "ekinetic_new.h"

#include "module_basis/module_ao/ORB_gen_tables.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_hamilt_lcao/hamilt_lcaodft/operator_lcao/operator_lcao.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer_funcs.h"
#include "module_base/timer.h"
#include "module_base/tool_title.h"

// Constructor
template <typename TK, typename TR>
hamilt::EkineticNew<hamilt::OperatorLCAO<TK, TR>>::EkineticNew(
    LCAO_Matrix* LM_in,
    const std::vector<ModuleBase::Vector3<double>>& kvec_d_in,
    hamilt::HContainer<TR>* hR_in,
    std::vector<TK>* hK_in,
    const UnitCell* ucell_in,
    Grid_Driver* GridD_in,
    const Parallel_Orbitals* paraV)
    : hamilt::OperatorLCAO<TK, TR>(LM_in, kvec_d_in, hR_in, hK_in)
{
    this->cal_type = lcao_fixed;
    this->ucell = ucell_in;
#ifdef __DEBUG
    assert(this->ucell != nullptr);
#endif
    // initialize HR to allocate sparse Ekinetic matrix memory
    this->initialize_HR(GridD_in, paraV);
}

// destructor
template <typename TK, typename TR>
hamilt::EkineticNew<hamilt::OperatorLCAO<TK, TR>>::~EkineticNew()
{
    if (this->allocated)
    {
        delete this->HR_fixed;
    }
}

// initialize_HR()
template <typename TK, typename TR>
void hamilt::EkineticNew<hamilt::OperatorLCAO<TK, TR>>::initialize_HR(Grid_Driver* GridD,
                                                                      const Parallel_Orbitals* paraV)
{
    ModuleBase::TITLE("EkineticNew", "initialize_HR");
    ModuleBase::timer::tick("EkineticNew", "initialize_HR");

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
            const LCAO_Orbitals& orb = LCAO_Orbitals::get_const_instance();
            // Note: the distance of atoms should less than the cutoff radius, 
            // When equal, the theoretical value of matrix element is zero, 
            // but the calculated value is not zero due to the numerical error, which would lead to result changes.
            if (this->ucell->cal_dtau(iat1, iat2, R_index2).norm() * this->ucell->lat0
                < orb.Phi[T1].getRcut() + orb.Phi[T2].getRcut())
            {
                is_adj[ad1] = true;
            }
        }
        filter_adjs(is_adj, adjs);
        this->adjs_all.push_back(adjs);
        for (int ad = 0; ad < adjs.adj_num + 1; ++ad)
        {
            const int T2 = adjs.ntype[ad];
            const int I2 = adjs.natom[ad];
            int iat2 = ucell->itia2iat(T2, I2);
            ModuleBase::Vector3<int>& R_index = adjs.box[ad];
            hamilt::AtomPair<TR> tmp(iat1, iat2, R_index.x, R_index.y, R_index.z, paraV);
            this->hR->insert_pair(tmp);
        }
    }
    // allocate the memory of BaseMatrix in HR, and set the new values to zero
    this->hR->allocate(true);

    ModuleBase::timer::tick("EkineticNew", "initialize_HR");
}

template <typename TK, typename TR>
void hamilt::EkineticNew<hamilt::OperatorLCAO<TK, TR>>::calculate_HR()
{
    ModuleBase::TITLE("EkineticNew", "calculate_HR");
    if(this->HR_fixed==nullptr || this->HR_fixed->size_atom_pairs()<=0)
    {
        ModuleBase::WARNING_QUIT("hamilt::EkineticNew::calculate_HR", "HR_fixed is nullptr or empty");
    }
    ModuleBase::timer::tick("EkineticNew", "calculate_HR");

    const Parallel_Orbitals* paraV = this->HR_fixed->get_atom_pair(0).get_paraV();
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

            hamilt::BaseMatrix<TR>* tmp = this->HR_fixed->find_matrix(iat1, iat2, R_index2.x, R_index2.y, R_index2.z);
            if (tmp != nullptr)
            {
                this->cal_HR_IJR(iat1, iat2, paraV, dtau, tmp->get_pointer());
            }
            else
            {
                ModuleBase::WARNING_QUIT("hamilt::EkineticNew::calculate_HR", "R_index not found in HR");
            }
        }
    }

    ModuleBase::timer::tick("EkineticNew", "calculate_HR");
}

// cal_HR_IJR()
template <typename TK, typename TR>
void hamilt::EkineticNew<hamilt::OperatorLCAO<TK, TR>>::cal_HR_IJR(const int& iat1,
                                                                   const int& iat2,
                                                                   const Parallel_Orbitals* paraV,
                                                                   const ModuleBase::Vector3<double>& dtau,
                                                                   TR* data_pointer)
{
    const ORB_gen_tables& uot = ORB_gen_tables::get_const_instance();
    const LCAO_Orbitals& orb = LCAO_Orbitals::get_const_instance();
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

    const int* iw2l1 = atom1.iw2l;
    const int* iw2n1 = atom1.iw2n;
    const int* iw2m1 = atom1.iw2m;
    const int* iw2l2 = atom2.iw2l;
    const int* iw2n2 = atom2.iw2n;
    const int* iw2m2 = atom2.iw2m;
#ifndef USE_NEW_TWO_CENTER
    // ---------------------------------------------
    // get tau1 (in cell <0,0,0>) and tau2 (in cell R)
    // in principle, only dtau is needed in this function
    // snap_psipsi should be refactored to use dtau directly
    // ---------------------------------------------
    const ModuleBase::Vector3<double>& tau1 = this->ucell->get_tau(iat1);
    const ModuleBase::Vector3<double> tau2 = tau1 + dtau;
#endif
    // ---------------------------------------------
    // calculate the Ekinetic matrix for each pair of orbitals
    // ---------------------------------------------
    double olm[3] = {0, 0, 0};
    auto row_indexes = paraV->get_indexes_row(iat1);
    auto col_indexes = paraV->get_indexes_col(iat2);
    const int step_trace = col_indexes.size() + 1;
    for (int iw1l = 0; iw1l < row_indexes.size(); iw1l += npol)
    {
        const int iw1 = row_indexes[iw1l] / npol;
        const int L1 = iw2l1[iw1];
        const int N1 = iw2n1[iw1];
        const int m1 = iw2m1[iw1];
#ifdef USE_NEW_TWO_CENTER
        int M1 = (m1 % 2 == 0) ? -m1/2 : (m1+1)/2;
#endif
        for (int iw2l = 0; iw2l < col_indexes.size(); iw2l += npol)
        {
            const int iw2 = col_indexes[iw2l] / npol;
            const int L2 = iw2l2[iw2];
            const int N2 = iw2n2[iw2];
            const int m2 = iw2m2[iw2];
#ifdef USE_NEW_TWO_CENTER
            //=================================================================
            //          new two-center integral (temporary)
            //=================================================================
            // convert m (0,1,...2l) to M (-l, -l+1, ..., l-1, l)
            int M2 = (m2 % 2 == 0) ? -m2/2 : (m2+1)/2;
            uot.two_center_bundle->kinetic_orb->calculate(T1, L1, N1, M1,
                    T2, L2, N2, M2, dtau * this->ucell->lat0, olm);
#else
            uot.snap_psipsi(orb, // orbitals
                            olm,
                            0,
                            'T', // olm, job of derivation, dtype of Operator
                            tau1,
                            T1,
                            L1,
                            m1,
                            N1, // info of atom1
                            tau2,
                            T2,
                            L2,
                            m2,
                            N2 // info of atom2
            );
#endif
            for (int ipol = 0; ipol < npol; ipol++)
            {
                data_pointer[ipol * step_trace] += olm[0];
            }
            data_pointer += npol;
        }
        data_pointer += (npol - 1) * col_indexes.size();
    }
}

// set_HR_fixed()
template <typename TK, typename TR>
void hamilt::EkineticNew<hamilt::OperatorLCAO<TK, TR>>::set_HR_fixed(void* HR_fixed_in)
{
    this->HR_fixed = static_cast<hamilt::HContainer<TR>*>(HR_fixed_in);
    this->allocated = false;
}

// contributeHR()
template <typename TK, typename TR>
void hamilt::EkineticNew<hamilt::OperatorLCAO<TK, TR>>::contributeHR()
{
    ModuleBase::TITLE("EkineticNew", "contributeHR");
    ModuleBase::timer::tick("EkineticNew", "contributeHR");

    if (!this->HR_fixed_done)
    {
        // if this Operator is the first node of the sub_chain, then HR_fixed is nullptr
        if (this->HR_fixed == nullptr)
        {
            this->HR_fixed = new hamilt::HContainer<TR>(*this->hR);
            this->HR_fixed->set_zero();
            this->allocated = true;
        }
        if(this->next_sub_op != nullptr)
        {
            // pass pointer of HR_fixed to the next node
            static_cast<OperatorLCAO<TK, TR>*>(this->next_sub_op)->set_HR_fixed(this->HR_fixed);
        }
        // calculate the values in HR_fixed
        this->calculate_HR();
        this->HR_fixed_done = true;
    }
    // last node of sub-chain, add HR_fixed into HR
    if(this->next_sub_op == nullptr)
    {
        this->hR->add(*(this->HR_fixed));
    }

    ModuleBase::timer::tick("EkineticNew", "contributeHR");
    return;
}

template class hamilt::EkineticNew<hamilt::OperatorLCAO<double, double>>;
template class hamilt::EkineticNew<hamilt::OperatorLCAO<std::complex<double>, double>>;
template class hamilt::EkineticNew<hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>>;