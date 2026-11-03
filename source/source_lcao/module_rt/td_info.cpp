#include "td_info.h"

#include "source_base/atom_basis_view.h"
#include "source_base/basis_flattened_indexing.h"
#include "source_base/block_cyclic_distribution_information.h"
#include "source_base/libm/libm.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "source_io/module_parameter/parameter.h"

bool TD_info::out_mat_R = false;
int TD_info::out_current = 0;
bool TD_info::out_current_k = false;
bool TD_info::evolve_once = false;

TD_info::TD_info(const UnitCell* ucell_in,
                 const BasisFlattenedIndexing& basis,
                 const BlockCyclicDistributionInformation& layout,
                 const LCAO_Orbitals& orb,
                 const elecstate::TDFieldManager& field_manager)
    : basis_(basis), layout_(layout), field_manager_(field_manager)
{
    if (out_current == 2 || field_manager_.gauge() == 2)
    {
        const bool cal_force = PARAM.inp.cal_force;
        const int nlocal = PARAM.globalv.nlocal;
        r_calculator.init(*ucell_in, layout, orb, cal_force, nlocal);
    }
    return;
}

AtomBasisView TD_info::row_basis(const int iat) const
{
    return AtomBasisView::rows(basis_, layout_, iat);
}

AtomBasisView TD_info::col_basis(const int iat) const
{
    return AtomBasisView::cols(basis_, layout_, iat);
}
TD_info::~TD_info()
{
    if (field_manager_.gauge() == 1)
    {
        this->destroy_HS_R_td_sparse();
    }
    for (int dir = 0; dir < 3; dir++)
    {
        if (this->current_term[dir] != nullptr)
        {
            delete this->current_term[dir];
        }
    }
    for (int dir = 0; dir < 3; dir++)
    {
        if (this->grad_overlap[dir] != nullptr)
        {
            delete this->grad_overlap[dir];
        }
    }
}

const ModuleBase::Vector3<double>& TD_info::endpoint_A() const
{
    return field_manager_.A_right_ha();
}

void TD_info::update_field()
{
    // update hybrid gauge phase
    if (field_manager_.gauge() == 2)
    {
        for (const auto& phase_pair: phase_hybrid)
        {
            const ModuleBase::Vector3<int>& r_index = phase_pair.first;
            ModuleBase::Vector3<double> dR = double(r_index.x) * a1 + double(r_index.y) * a2 + double(r_index.z) * a3;
            const double arg_td = endpoint_A() * dR * lat0;
            double sinp, cosp;
            ModuleBase::libm::sincos(arg_td, &sinp, &cosp);
            phase_hybrid[r_index] = std::complex<double>(cosp, sinp);
        }
    }
}

template <typename TR>
void TD_info::initialize_phase_hybrid(const UnitCell& ucell, const hamilt::HContainer<TR>* hR)
{
    this->a1 = ucell.a1;
    this->a2 = ucell.a2;
    this->a3 = ucell.a3;
    this->lat0 = ucell.lat0;
    for (int i = 0; i < hR->size_atom_pairs(); ++i)
    {
        hamilt::AtomPair<TR>& tmp = hR->get_atom_pair(i);
        for (int ir = 0; ir < tmp.get_R_size(); ++ir)
        {
            const ModuleBase::Vector3<int> r_index = tmp.get_R_index(ir);
            if (phase_hybrid.count(r_index))
                continue;

            ModuleBase::Vector3<double> dR = double(r_index.x) * a1 + double(r_index.y) * a2 + double(r_index.z) * a3;
            const double arg_td = endpoint_A() * dR * lat0;
            double sinp, cosp;
            ModuleBase::libm::sincos(arg_td, &sinp, &cosp);
            phase_hybrid[r_index] = std::complex<double>(cosp, sinp);
        }
    }
}
void TD_info::initialize_current_term(const hamilt::HContainer<std::complex<double>>* HR)
{
    ModuleBase::TITLE("TD_info", "initialize_current_term");
    ModuleBase::timer::start("TD_info", "initialize_current_term");

    for (int dir = 0; dir < 3; dir++)
    {
        delete this->current_term[dir];
        this->current_term[dir] = new hamilt::HContainer<std::complex<double>>(*HR);
    }

    ModuleBase::timer::end("TD_info", "initialize_current_term");
}

void TD_info::destroy_HS_R_td_sparse(void)
{
    std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, std::complex<double>>>> empty_HR_sparse_td_vel_up;
    std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, std::complex<double>>>> empty_HR_sparse_td_vel_down;
    HR_sparse_td_vel[0].swap(empty_HR_sparse_td_vel_up);
    HR_sparse_td_vel[1].swap(empty_HR_sparse_td_vel_down);
}

void TD_info::calculate_grad_overlap(const UnitCell& ucell,
                                     const Grid_Driver& GridD,
                                     const std::vector<double>& orb_cutoff,
                                     const TwoCenterIntegrator* intor)
{
    ModuleBase::TITLE("TD_info", "calculate_grad_overlap");
    ModuleBase::timer::start("TD_info", "calculate_grad_overlap");
    for (int dir = 0; dir < 3; dir++)
    {
        if (this->grad_overlap[dir] != nullptr)
        {
            delete this->grad_overlap[dir];
        }
        this->grad_overlap[dir] =
            new hamilt::HContainer<double>(basis_, layout_, static_cast<double*>(nullptr));
    }
    for (int iat1 = 0; iat1 < ucell.nat; iat1++)
    {
        auto tau1 = ucell.get_tau(iat1);
        int T1 = 0;
        int I1 = 0;
        ucell.iat2iait(iat1, &I1, &T1);
        AdjacentAtomInfo adjs;
        GridD.Find_atom(ucell, tau1, T1, I1, &adjs);
        for (int ad = 0; ad < adjs.adj_num + 1; ++ad)
        {
            const int T2 = adjs.ntype[ad];
            const int I2 = adjs.natom[ad];
            int iat2 = ucell.itia2iat(T2, I2);
            const ModuleBase::Vector3<int>& R_index = adjs.box[ad];
            // choose the real adjacent atoms
            // Note: the distance of atoms should less than the cutoff radius,
            // When equal, the theoretical value of matrix element is zero,
            // but the calculated value is not zero due to the numerical error, which would lead to result changes.
            if (ucell.cal_dtau(iat1, iat2, R_index).norm() * ucell.lat0 >= orb_cutoff[T1] + orb_cutoff[T2])
            {
                continue;
            }
            for (int dir = 0; dir < 3; dir++)
            {
                this->grad_overlap[dir]->insert_pair(iat1, iat2, R_index);
            }
        }
    }
    // allocate the memory of BaseMatrix in grad_overlap, and set the new values to zero
    for (int dir = 0; dir < 3; dir++)
    {
        this->grad_overlap[dir]->allocate(nullptr, true);
    }

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int iap = 0; iap < this->grad_overlap[0]->size_atom_pairs(); ++iap)
    {
        hamilt::AtomPair<double>& tmp = this->grad_overlap[0]->get_atom_pair(iap);
        const int iat1 = tmp.get_atom_i();
        const int iat2 = tmp.get_atom_j();

        for (int iR = 0; iR < tmp.get_R_size(); ++iR)
        {
            const ModuleBase::Vector3<int> R_index = tmp.get_R_index(iR);
            auto dtau = ucell.cal_dtau(iat1, iat2, R_index);

            double* p_data[3] = {nullptr, nullptr, nullptr};
            for (int i = 0; i < 3; i++)
            {
                p_data[i] = this->grad_overlap[i]->get_atom_pair(iap).get_pointer(iR);
            }
            // ---------------------------------------------
            // get info of orbitals of atom1 and atom2 from ucell
            // ---------------------------------------------
            int T1 = 0;
            int I1 = 0;
            ucell.iat2iait(iat1, &I1, &T1);
            int T2 = 0;
            int I2 = 0;
            ucell.iat2iait(iat2, &I2, &T2);
            Atom& atom1 = ucell.atoms[T1];
            Atom& atom2 = ucell.atoms[T2];

            // npol is the number of polarizations,
            // 1 for non-magnetic (one Hamiltonian matrix only has spin-up or spin-down),
            // 2 for magnetic (one Hamiltonian matrix has both spin-up and spin-down)
            const int npol = ucell.get_npol();

            const int* iw2l1 = atom1.iw2l.data();
            const int* iw2n1 = atom1.iw2n.data();
            const int* iw2m1 = atom1.iw2m.data();
            const int* iw2l2 = atom2.iw2l.data();
            const int* iw2n2 = atom2.iw2n.data();
            const int* iw2m2 = atom2.iw2m.data();

            // ---------------------------------------------
            // calculate the overlap matrix for each pair of orbitals
            // ---------------------------------------------
            double olm[3] = {0, 0, 0};
            const AtomBasisView rows = this->row_basis(iat1);
            const AtomBasisView cols = this->col_basis(iat2);
            const int step_trace = cols.nloc() + 1;
            for (int iw1l = 0; iw1l < rows.nloc(); iw1l += npol)
            {
                const int iw1 = rows[iw1l] / npol;
                const int L1 = iw2l1[iw1];
                const int N1 = iw2n1[iw1];
                const int m1 = iw2m1[iw1];

                // convert m (0,1,...2l) to M (-l, -l+1, ..., l-1, l)
                int M1 = (m1 % 2 == 0) ? -m1 / 2 : (m1 + 1) / 2;

                for (int iw2l = 0; iw2l < cols.nloc(); iw2l += npol)
                {
                    const int iw2 = cols[iw2l] / npol;
                    const int L2 = iw2l2[iw2];
                    const int N2 = iw2n2[iw2];
                    const int m2 = iw2m2[iw2];

                    // convert m (0,1,...2l) to M (-l, -l+1, ..., l-1, l)
                    int M2 = (m2 % 2 == 0) ? -m2 / 2 : (m2 + 1) / 2;
                    intor->calculate(T1, L1, N1, M1, T2, L2, N2, M2, dtau * ucell.lat0, nullptr, olm);
                    for (int dir = 0; dir < 3; dir++)
                    {
                        for (int ipol = 0; ipol < npol; ipol++)
                        {
                            p_data[dir][ipol * step_trace] += olm[dir];
                        }
                        p_data[dir] += npol;
                    }
                }
                for (int dir = 0; dir < 3; dir++)
                {
                    p_data[dir] += (npol - 1) * cols.nloc();
                }
            }
        }
    }
    ModuleBase::timer::end("TD_info", "calculate_grad_overlap");
}
template void TD_info::initialize_phase_hybrid<std::complex<double>>(const UnitCell& ucell,
                                                                     const hamilt::HContainer<std::complex<double>>* hR);
template void TD_info::initialize_phase_hybrid<double>(const UnitCell& ucell, const hamilt::HContainer<double>* hR);
