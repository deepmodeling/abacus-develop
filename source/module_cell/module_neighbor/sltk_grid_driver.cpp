#include "sltk_grid_driver.h"
#include "module_parameter/parameter.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/timer.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace GlobalC
{
Grid_Driver GridD;
}
Grid_Driver::Grid_Driver(
	const int &test_d_in, 
	const int &test_grid_in)
:test_deconstructor(test_d_in),
Grid(test_grid_in)
{
	//	ModuleBase::TITLE("Grid_Driver","Grid_Driver");
}

Grid_Driver::~Grid_Driver()
{
}


void Grid_Driver::Find_atom(
	const UnitCell &ucell, 
	const ModuleBase::Vector3<double> &cartesian_pos, 
	const int &ntype, 
	const int &nnumber,
	AdjacentAtomInfo *adjs)
{
	ModuleBase::timer::tick("Grid_Driver","Find_atom");
//	std::cout << "lenght in Find atom = " << atomlink[offset].fatom.getAdjacentSet()->getLength() << std::endl;

	// store result in member adj_info when parameter adjs is NULL
	AdjacentAtomInfo* local_adjs = adjs == nullptr ? &this->adj_info : adjs;
	local_adjs->clear();
	const std::vector<FAtom *> & all_atom = Cell[this->true_cell_x][this->true_cell_y][this->true_cell_z].atom_map[ntype][nnumber].getAdjacent();
	//std::cout << "ntype = "<< ntype << "  atom size = " << all_atom.size() << std::endl;

	ModuleBase::Vector3<double> vec1(ucell.latvec.e11, ucell.latvec.e12, ucell.latvec.e13);
	ModuleBase::Vector3<double> vec2(ucell.latvec.e21, ucell.latvec.e22, ucell.latvec.e23);
	ModuleBase::Vector3<double> vec3(ucell.latvec.e31, ucell.latvec.e32, ucell.latvec.e33);

	for(const FAtom * atom : all_atom)
	{
		// std::cout << "atom type = " << atom.getType() << " number = " << atom.getNatom() << " box = " << atom.getCellX() << " " << atom.getCellY() << " " << atom.getCellZ() 
		// << " tau = " << atom.x() << " " << atom.y() << " " << atom.z() << std::endl;
		local_adjs->ntype.push_back(atom->getType());
		local_adjs->natom.push_back(atom->getNatom());
		local_adjs->box.push_back(ModuleBase::Vector3<int>(atom->getCellX(), atom->getCellY(), atom->getCellZ()));
		if (expand_flag)
		{
			local_adjs->adjacent_tau.push_back(ModuleBase::Vector3<double>(atom->x(), atom->y(), atom->z()));
		}
		else
		{
			local_adjs->adjacent_tau.push_back(Calculate_adjacent_site(atom->x(), atom->y(), atom->z(),
			                  vec1[0], vec2[0], vec3[0],
			                  vec1[1], vec2[1], vec3[1],
			                  vec1[2], vec2[2], vec3[2],
			                  atom->getCellX(), atom->getCellY(), atom->getCellZ()));
		}//end if expand_flag
		local_adjs->adj_num++;
	}
	local_adjs->ntype.push_back(ntype);
	local_adjs->natom.push_back(nnumber);
	local_adjs->box.push_back(ModuleBase::Vector3<int>(0, 0, 0));
	local_adjs->adjacent_tau.push_back(ModuleBase::Vector3<double>(cartesian_pos.x, cartesian_pos.y, cartesian_pos.z));


	ModuleBase::timer::tick("Grid_Driver","Find_atom");
	return;
}

ModuleBase::Vector3<double> Grid_Driver::Calculate_adjacent_site(const double x, const double y, const double z,
																	const double &box11, const double &box12, const double &box13,
																	const double &box21, const double &box22, const double &box23,
																	const double &box31, const double &box32, const double &box33,
																	const short box_x, const short box_y, const short box_z) const
{
	ModuleBase::Vector3<double> adjacent_site(0, 0, 0);
	adjacent_site.x = x + box_x * box11 + box_y * box12 + box_z * box13;
	adjacent_site.y = y + box_x * box21 + box_y * box22 + box_z * box23;
	adjacent_site.z = z + box_x * box31 + box_y * box32 + box_z * box33;

	return adjacent_site;
}

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
			adjs.adjacent_tau.erase(adjs.adjacent_tau.begin()+i);//info of adjacent_tau is not used in future
			adjs.box.erase(adjs.box.begin()+i);
		}
	}
}
