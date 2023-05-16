#include "../sltk_atom_arrange.h"
#include<iostream>
#include<string>

namespace ModuleBase
{
	void TITLE(const std::string &class_name,const std::string &function_name,const bool disable)
{
	
}
	void TITLE(std::ofstream &ofs,const std::string &class_name,const std::string &function_name,const bool disable){}
	
	class timer
	{
	private:
	
	public:
		timer();
		~timer();
	void tick(const std::string &class_name,const std::string &name);
	};
		timer::timer(){}
		timer::~timer(){}
	void timer::tick(const std::string &class_name,const std::string &name){};
	
	void Matrix3::Identity(void)
	{
		e11 = 1;e12 = 0;e13 = 0;
		e21 = 0;e22 = 1;e23 = 0;
		e31 = 0;e32 = 0;e33 = 1;
	}
	IntArray::~IntArray(){}
	IntArray::IntArray(int, int){}
}
//constructions and destructions
Magnetism::Magnetism(){}
Magnetism::~Magnetism(){}
UnitCell::~UnitCell(){}
#ifdef __LCAO
InfoNonlocal::InfoNonlocal(){}
InfoNonlocal::~InfoNonlocal(){}
#endif
UnitCell::UnitCell(){}
Grid_Driver::Grid_Driver(
	const int &test_d_in, 
	const int &test_gd_in, 
	const int &test_grid_in)
:test_deconstructor(test_d_in),
test_grid_driver(test_gd_in),
Grid(test_grid_in){}
Grid_Driver::~Grid_Driver(){}
AtomLink::AtomLink(const FAtom &atom, AtomLink* const pointNext)
: fatom(atom), next_p(pointNext) {}
Grid::Grid(const int &test_grid_in):test_grid(test_grid_in)
{
	init_cell_flag = false;
	this->atomlink = new AtomLink[1];
}
Grid::~Grid()
{
	delete[] atomlink;
	this->delete_Cell();
	
}
CellSet::CellSet()
: address(NullPtr), length(0)
{
	in_grid[0] = 0;
	in_grid[1] = 0;
	in_grid[2] = 0;
}
FAtom::FAtom()
{
	d_x = 0.0;	
	d_y = 0.0;	
	d_z = 0.0;	
	as = nullptr;
	type = 0;		
	natom = 0;
}
FAtom::~FAtom(){}
Atom_input::Atom_input
(
	std::ofstream &ofs_in,
	const UnitCell &ucell,
	const int amount,
	const int ntype,
	const bool boundary_in,
	const double radius_in,
	const int &test_atom_in
	)
:d_amount(amount),
d_amount_expand(amount),
periodic_boundary(boundary_in),
radius(radius_in),
expand_flag(false),
glayerX(1),
glayerX_minus(0),
glayerY(1),
glayerY_minus(0),
glayerZ(1),
glayerZ_minus(0),
x_min_expand(0),
y_min_expand(0),
z_min_expand(0),
x_max_expand(0),
y_max_expand(0),
z_max_expand(0),
d_current(0),
test_atom_input(test_atom_in),
type(0),
natom(-1)
{}
Atom_input::~Atom_input(){}


void Grid::delete_vector(const Atom_input &input)
{
	std::cout<<"test"<<"\n";
}
void Grid::init(std::basic_ofstream<char, std::char_traits<char> >&, UnitCell const&, Atom_input const&){};
void Grid_Driver::Find_atom(UnitCell const&, ModuleBase::Vector3<double> const&, int const&, int const&, AdjacentAtomInfo*){};

