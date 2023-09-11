
#include "module_cell/unitcell.h"

// constructor of Atom
Atom::Atom()
{
}
Atom::~Atom()
{
}

Atom_pseudo::Atom_pseudo()
{
}
Atom_pseudo::~Atom_pseudo()
{
}

Magnetism::Magnetism()
{
}
Magnetism::~Magnetism()
{
}

InfoNonlocal::InfoNonlocal()
{
}
InfoNonlocal::~InfoNonlocal()
{
}

pseudo_nc::pseudo_nc()
{
}
pseudo_nc::~pseudo_nc()
{
}

// constructor of UnitCell
UnitCell::UnitCell()
{
}
UnitCell::~UnitCell()
{
}

void UnitCell::set_iat2iwt(const int& npol_in)
{
	this->iat2iwt.resize(this->nat);
	this->npol = npol_in;
	int iat=0;
	int iwt=0;
	for(int it = 0;it < this->ntype; it++)
	{
		for(int ia=0; ia<atoms[it].na; ia++)
		{
			this->iat2iwt[iat] = iwt;
			iwt += atoms[it].nw * this->npol;
			++iat;
		}	
	}
	return;
}

// mock of OperatorLCAO
#include "module_hamilt_lcao/hamilt_lcaodft/operator_lcao/operator_lcao.h"

/* template<typename FPTYPE, typename Device>
hamilt::Operator<FPTYPE, Device>::Operator(){}

template<typename FPTYPE, typename Device>
hamilt::Operator<FPTYPE, Device>::~Operator(){}

template<typename FPTYPE, typename Device>
typename hamilt::Operator<FPTYPE, Device>::hpsi_info hamilt::Operator<FPTYPE, Device>::hPsi(hpsi_info&) const 
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

template class hamilt::Operator<double, psi::DEVICE_CPU>;
template class hamilt::Operator<std::complex<double>, psi::DEVICE_CPU>;*/

// mock of OperatorLCAO
template<typename TK>
void hamilt::OperatorLCAO<TK>::init(const int ik_in)
{
    if(!this->hr_done)
    {
        OperatorLCAO<TK>* last = this;
        while(last != nullptr)
        {
            last->contributeHR();
            last = dynamic_cast<OperatorLCAO<TK>*>(last->next_sub_op);
        }
        this->hr_done = true;
    }
    this->contributeHk(ik_in);
    return;
}
template<typename TK>
void hamilt::OperatorLCAO<TK>::get_hs_pointers()
{
    return;
}
template class hamilt::OperatorLCAO<double>;
template class hamilt::OperatorLCAO<std::complex<double>>;

// mock of ORB_gen_tables and LCAO_Orbitals
#include "module_basis/module_ao/ORB_gen_tables.h"
const ORB_gen_tables& ORB_gen_tables::get_const_instance()
{
    static ORB_gen_tables instance;
    return instance;
}
ORB_gen_tables::ORB_gen_tables() {}
ORB_gen_tables::~ORB_gen_tables() {}
ORB_gaunt_table::ORB_gaunt_table() {}
ORB_gaunt_table::~ORB_gaunt_table() {}
ORB_table_phi::ORB_table_phi() {}
ORB_table_phi::~ORB_table_phi() {}
ORB_table_alpha::ORB_table_alpha() {}
ORB_table_alpha::~ORB_table_alpha() {}
ORB_table_beta::ORB_table_beta() {}
ORB_table_beta::~ORB_table_beta() {}
// mock of snap_psipsi
void ORB_gen_tables::snap_psipsi(
    const LCAO_Orbitals &orb,
    double olm[],
    const int &job, ///<[in]0 for matrix element of either S or T, 1 for its derivatives
    const char &dtype, ///<[in] derivative type, 'S' for overlap, 'T' for kinetic energy, 'D' for descriptor in deepks
    const ModuleBase::Vector3<double> &R1,
    const int &I1,
    const int &l1,
    const int &m1,
    const int &n1,
    const ModuleBase::Vector3<double> &R2,
    const int &I2,
    const int &l2,
    const int &m2,
    const int &n2,
    bool cal_syns,
    double dmax)const
{
    if(dtype == 'S')
    {
        olm[0] = 1.0;
    }
    else if(dtype == 'T')
    {
        olm[0] = 2.0;
    }
    else if(dtype == 'D')
    {
        olm[0] = 3.0;
    }
}

void ORB_gen_tables::snap_psibeta_half(
	const LCAO_Orbitals &orb,
	const InfoNonlocal &infoNL_,
	std::vector<std::vector<double>> &nlm,
	const ModuleBase::Vector3<double> &R1,
	const int &T1,
	const int &L1,
	const int &m1,
	const int &N1,
	const ModuleBase::Vector3<double> &R0, // The projector.
	const int &T0,
	const bool &calc_deri)const // mohan add 2021-04-25)
{
    nlm.resize(1);
    for(int i = 0; i < nlm.size(); ++i)
    {
        nlm[i].resize(5, 1.0);
    }
}

#include "module_basis/module_ao/ORB_read.h"
const LCAO_Orbitals& LCAO_Orbitals::get_const_instance()
{
    static LCAO_Orbitals instance;
    return instance;
}
LCAO_Orbitals::LCAO_Orbitals() {this->Phi = new Numerical_Orbital[1];}
LCAO_Orbitals::~LCAO_Orbitals() { delete[] Phi; }

#include "module_cell/module_neighbor/sltk_grid_driver.h"
// mock find_atom() function
void Grid_Driver::Find_atom(const UnitCell& ucell,
                            const ModuleBase::Vector3<double>& tau,
                            const int& T,
                            const int& I,
                            AdjacentAtomInfo* adjs)
{
    adjs->adj_num = ucell.nat - 1 ;
    adjs->adjacent_tau.resize(ucell.nat);
    adjs->ntype.resize(ucell.nat, 0);
    adjs->natom.resize(ucell.nat);
    adjs->box.resize(ucell.nat);
    for(int iat = 0;iat<ucell.nat;iat++)
    {
        adjs->natom[iat] = iat;
        adjs->box[iat].x = 1;
        adjs->box[iat].y = 1;
        adjs->box[iat].z = 1;
        adjs->adjacent_tau[iat] = ucell.get_tau(iat);
    }
}
Grid::Grid(const int &test_grid_in):test_grid(test_grid_in)
{}
Grid::~Grid(){}
Grid_Driver::Grid_Driver(const int &test_d_in, 
		const int &test_gd_in, 
		const int &test_grid_in) :Grid(test_grid_in), test_deconstructor(test_d_in), test_grid_driver(test_gd_in) {}
Grid_Driver::~Grid_Driver() {}

// filter_adjs delete not adjacent atoms in adjs
void filter_adjs(const std::vector<bool>& is_adj, AdjacentAtomInfo& adjs)
{
	const int size = adjs.adj_num+1;
	for(int i = size-1; i >= 0; --i)
	{
		if(!is_adj[i])
		{
			adjs.adj_num--;
			adjs.ntype.erase(adjs.ntype.begin()+i);
			adjs.natom.erase(adjs.natom.begin()+i);
			adjs.adjacent_tau.erase(adjs.adjacent_tau.begin()+i);
			adjs.box.erase(adjs.box.begin()+i);
		}
	}
}

Numerical_Nonlocal::Numerical_Nonlocal()
{
    this->rcut_max = 1.0;
}
Numerical_Nonlocal::~Numerical_Nonlocal()
{
}

Numerical_Orbital::Numerical_Orbital()
{
    this->rcut = 1.0;
}
Numerical_Orbital::~Numerical_Orbital()
{
}
