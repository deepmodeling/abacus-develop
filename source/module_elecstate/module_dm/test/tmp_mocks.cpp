#include "module_cell/klist.h"

K_Vectors::K_Vectors()
{
}

K_Vectors::~K_Vectors()
{
}

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