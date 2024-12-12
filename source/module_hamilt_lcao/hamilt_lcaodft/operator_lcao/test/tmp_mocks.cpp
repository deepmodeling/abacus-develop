
#include "module_cell/unitcell.h"

// constructor of Atom
Atom::Atom() {}
Atom::~Atom() {}

Atom_pseudo::Atom_pseudo() {}
Atom_pseudo::~Atom_pseudo() {}

Magnetism::Magnetism() {}
Magnetism::~Magnetism() {}

InfoNonlocal::InfoNonlocal() {}
InfoNonlocal::~InfoNonlocal() {}

pseudo::pseudo() {}
pseudo::~pseudo() {}

// constructor of UnitCell
UnitCell::UnitCell() {}
UnitCell::~UnitCell() {}

void UnitCell::set_iat2iwt(const int& npol_in) {
    this->iat2iwt.resize(this->nat);
    this->npol = npol_in;
    int iat = 0;
    int iwt = 0;
    for (int it = 0; it < this->ntype; it++) {
        for (int ia = 0; ia < atoms[it].na; ia++) {
            this->iat2iwt[iat] = iwt;
            iwt += atoms[it].nw * this->npol;
            ++iat;
        }
    }
    return;
}

// mock of OperatorLCAO
#include "module_hamilt_lcao/hamilt_lcaodft/operator_lcao/operator_lcao.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer_funcs.h"

/* template<typename FPTYPE, typename Device>
hamilt::Operator<FPTYPE, Device>::Operator(){}

template<typename FPTYPE, typename Device>
hamilt::Operator<FPTYPE, Device>::~Operator(){}

template<typename FPTYPE, typename Device>
typename hamilt::Operator<FPTYPE, Device>::hpsi_info hamilt::Operator<FPTYPE,
Device>::hPsi(hpsi_info&) const
{
    return hpsi_info(nullptr, 0, nullptr);
}

template<typename FPTYPE, typename Device>
void hamilt::Operator<FPTYPE, Device>::init(const int ik_in)
{
    return;
}

template<typename FPTYPE, typename Device>
void hamilt::Operator<FPTYPE, Device>::add(Operator* next)
{
    return;
}

template<typename FPTYPE, typename Device>
FPTYPE* hamilt::Operator<FPTYPE, Device>::get_hpsi(const hpsi_info& info) const
{
    return nullptr;
}

template class hamilt::Operator<double, base_device::DEVICE_CPU>;
template class hamilt::Operator<std::complex<double>,
base_device::DEVICE_CPU>;*/

// mock of OperatorLCAO
template <typename TK, typename TR>
void hamilt::OperatorLCAO<TK, TR>::init(const int ik_in) {
    if (!this->hr_done) {
        OperatorLCAO<TK, TR>* last = this;
        while (last != nullptr) {
            last->contributeHR();
            last = dynamic_cast<OperatorLCAO<TK, TR>*>(last->next_sub_op);
        }
        this->hr_done = true;
    }
    this->contributeHk(ik_in);
    return;
}
template <typename TK, typename TR>
void hamilt::OperatorLCAO<TK, TR>::contributeHk(int ik) {
    if (!this->is_first_node) {
        return;
    } else {
        const int ncol = this->hR->get_atom_pair(0).get_paraV()->get_col_size();
        hamilt::folding_HR(*this->hR, this->hsk->get_hk(), this->kvec_d[ik], ncol, 0);
    }
}
template <typename TK, typename TR>
void hamilt::OperatorLCAO<TK, TR>::get_hs_pointers() {
    return;
}
template class hamilt::OperatorLCAO<double, double>;
template class hamilt::OperatorLCAO<std::complex<double>, double>;
template class hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>;

// mock of TwoCenterIntegrator and LCAO_Orbitals
#include "module_basis/module_nao/two_center_integrator.h"
TwoCenterIntegrator::TwoCenterIntegrator() {}

void TwoCenterIntegrator::tabulate(const RadialCollection& bra,
                                   const RadialCollection& ket,
                                   const char op,
                                   const int nr,
                                   const double cutoff) {}

void TwoCenterIntegrator::calculate(
    const int itype1,
    const int l1,
    const int izeta1,
    const int m1,
    const int itype2,
    const int l2,
    const int izeta2,
    const int m2,
    const ModuleBase::Vector3<double>& vR, // vR = R2 - R1
    double* out,
    double* grad_out) const {
    out[0] = 1.0;
}

void TwoCenterIntegrator::snap(
    const int itype1,
    const int l1,
    const int izeta1,
    const int m1,
    const int itype2,
    const ModuleBase::Vector3<double>& vR, // vR = R2 - R1
    const bool deriv,
    std::vector<std::vector<double>>& out) const {
    out.resize(1);
    for (int i = 0; i < out.size(); ++i) {
        out[i].resize(5, 1.0);
    }
}

#include "module_basis/module_ao/ORB_read.h"
LCAO_Orbitals::LCAO_Orbitals() { this->Phi = new Numerical_Orbital[1]; }
LCAO_Orbitals::~LCAO_Orbitals() { delete[] Phi; }

#include "module_cell/module_neighbor/sltk_grid_driver.h"
// mock find_atom() function
void Grid_Driver::Find_atom(const UnitCell& ucell,
                            const ModuleBase::Vector3<double>& tau,
                            const int& T,
                            const int& I,
                            AdjacentAtomInfo* adjs) {
    adjs->adj_num = ucell.nat - 1;
    adjs->adjacent_tau.resize(ucell.nat);
    adjs->ntype.resize(ucell.nat, 0);
    adjs->natom.resize(ucell.nat);
    adjs->box.resize(ucell.nat);
    for (int iat = 0; iat < ucell.nat; iat++) {
        adjs->natom[iat] = iat;
        adjs->box[iat].x = 1;
        adjs->box[iat].y = 1;
        adjs->box[iat].z = 1;
        adjs->adjacent_tau[iat] = ucell.get_tau(iat);
    }
}
Grid::Grid(const int& test_grid_in) : test_grid(test_grid_in) {}
Grid::~Grid() {}
Grid_Driver::Grid_Driver(const int& test_d_in,
                         const int& test_grid_in)
    : Grid(test_grid_in), test_deconstructor(test_d_in) {
}
Grid_Driver::~Grid_Driver() {}

// filter_adjs delete not adjacent atoms in adjs
void filter_adjs(const std::vector<bool>& is_adj, AdjacentAtomInfo& adjs) {
    const int size = adjs.adj_num + 1;
    for (int i = size - 1; i >= 0; --i) {
        if (!is_adj[i]) {
            adjs.adj_num--;
            adjs.ntype.erase(adjs.ntype.begin() + i);
            adjs.natom.erase(adjs.natom.begin() + i);
            adjs.adjacent_tau.erase(adjs.adjacent_tau.begin() + i);
            adjs.box.erase(adjs.box.begin() + i);
        }
    }
}

Numerical_Nonlocal::Numerical_Nonlocal() { this->rcut_max = 1.0; }
Numerical_Nonlocal::~Numerical_Nonlocal() {}

Numerical_Orbital::Numerical_Orbital() { this->rcut = 1.0; }
Numerical_Orbital::~Numerical_Orbital() {}

void Numerical_Orbital::set_orbital_info(const int&,
                                         const std::string&,
                                         const int&,
                                         const int*,
                                         const int&) {}
