#include "sltk_grid_driver.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/timer.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace GlobalC
{
Grid_Driver GridD(GlobalV::test_deconstructor, GlobalV::test_grid_driver, GlobalV::test_grid);
}

Grid_Driver::Grid_Driver(
	const int &test_d_in, 
	const int &test_gd_in, 
	const int &test_grid_in)
:test_deconstructor(test_d_in),
test_grid_driver(test_gd_in),
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
	//if (test_grid_driver) ModuleBase::TITLE(GlobalV::ofs_running, "Grid_Driver", "Find_atom");
 #ifdef _OPENMP
    int in_parallel = omp_in_parallel();
 #else
    int in_parallel = 0;
#endif
    // to avoid writing conflict, disable timing in omp parallel block
	if (!in_parallel) ModuleBase::timer::tick("Grid_Driver","Find_atom");

	if (test_grid_driver > 1)
	{
//		GlobalV::ofs_running << " Atom cartesian_pos = " << cartesian_pos.x
//		<< " " << cartesian_pos.y
//		<< " " << cartesian_pos.z << std::endl;
	}

//----------------------------------------------------------
// CALL MEMBER FUNCTION :
// NAME : Locate_offset (Using Hash method to get the position)
// NAME : Find_adjacent_Atom ( find the adjacent information)
//----------------------------------------------------------
	//const int offset = this->Locate_offset(cartesian_pos);
	const int offset = this->Locate_offset(ucell, cartesian_pos, ntype, nnumber);

//	std::cout << "lenght in Find atom = " << atomlink[offset].fatom.getAdjacentSet()->getLength() << std::endl;

	// store result in member adj_info when parameter adjs is NULL
	AdjacentAtomInfo* local_adjs = adjs == nullptr ? &this->adj_info : adjs;
	this->Find_adjacent_atom(offset, this->atomlink[offset].fatom.getAdjacentSet(), *local_adjs);

	if (!in_parallel) ModuleBase::timer::tick("Grid_Driver","Find_atom");
	return;
}


int Grid_Driver::Locate_offset(
	const UnitCell &ucell, 
	const ModuleBase::Vector3<double> &cartesian_pos, 
	const int &ntype, 
	const int &nnumber)const
{
//	if (test_grid_driver) ModuleBase::TITLE(GlobalV::ofs_running, "Grid_Driver", "Locate_offset");

//----------------------------------------------------------
// EXPLAIN : Create an AtomLink object
//----------------------------------------------------------
	AtomLink temp;

	temp.fatom.setX(cartesian_pos.x);
	temp.fatom.setY(cartesian_pos.y);
	temp.fatom.setZ(cartesian_pos.z);
	temp.fatom.setType(ntype);
	temp.fatom.setNatom(nnumber);

//----------------------------------------------------------
// EXPLAIN : Find the Hash number of this atom position
//----------------------------------------------------------
	AtomLink* Search = this->getHashCode(ucell, temp.fatom);

//----------------------------------------------------------
// EXPLAIN : If we don't get the index for one Hash try,
// we need to search in array.
//----------------------------------------------------------
	for (; Search < this->cordon_p ; Search = Search->next_p)
	{
		if (Search->fatom.x() == cartesian_pos.x &&
		        Search->fatom.y() == cartesian_pos.y &&
		        Search->fatom.z() == cartesian_pos.z)
		{
			const int offset = Search - this->atomlink;
			return offset;
		}
	}

	// Peize Lin update 2019-05-01
	// throw std::runtime_error("Locate_Atom wrong. "+TO_STRING(__FILE__)+" line "+TO_STRING(__LINE__));
	// mohan update 2021-06-21
	ModuleBase::WARNING_QUIT("Locate_Atom", "something wrong!");

	return 0; //meaningless mohan add 2021-06-21

}

void Grid_Driver::Find_adjacent_atom(const int offset, std::shared_ptr<AdjacentSet> as, AdjacentAtomInfo &adjs) const
{
//	if (test_grid_driver) ModuleBase::TITLE(GlobalV::ofs_running, "Grid_Driver", "Find_adjacent_atom");

	// alias of variables
	auto &adj_num = adjs.adj_num;
	auto &ntype = adjs.ntype;
	auto &natom = adjs.natom;
	auto &adjacent_tau = adjs.adjacent_tau;
	auto &box = adjs.box;

//----------------------------------------------------------
// CALL OTHER CLASS MEMBER FUNCTION :
// NAME : getLength(get the adjacent number of this atom)
//----------------------------------------------------------
	adj_num = as->getLength();

	if (test_grid_driver > 1) 
	{
//		OUT(GlobalV::ofs_running,"adj_num",this->adj_num);
	}	
	
	//std::cout << "\n length = " << adj_num << std::endl;
	//BLOCK_HERE("Find_adjacent_atom");

//----------------------------------------------------------
// add center atom at the end of these arrays
// ywcui add 2009-03-30
//----------------------------------------------------------
	ntype.resize(adj_num+1);
	natom.resize(adj_num+1);
	adjacent_tau.resize(adj_num+1);
	box.resize(adj_num+1);

	// the last one is the atom itself.
	ntype[adj_num] = this->atomlink[offset].fatom.getType();
	natom[adj_num] = this->atomlink[offset].fatom.getNatom();
	adjacent_tau[adj_num].x = this->atomlink[offset].fatom.x();
	adjacent_tau[adj_num].y = this->atomlink[offset].fatom.y();
	adjacent_tau[adj_num].z = this->atomlink[offset].fatom.z();
	// mohan add 2010-07-01
	box[adj_num].x = 0;
	box[adj_num].y = 0;
	box[adj_num].z = 0;

	if (test_grid_driver > 1)
	{
/*
		GlobalV::ofs_running << " "
		<< std::setw(15) << "Box"
		<< std::setw(8) << "Offset"
		<< std::setw(5) << "Type"
		<< std::setw(8) << "Natom"
		<< std::setw(10) << "X"
		<< std::setw(10) << "Y"
		<< std::setw(10) << "Z"
		<< std::setw(10) << "Distance" << std::endl;
*/
	}


	for (int i = 0;i < adj_num;i++)
	{
//----------------------------------------------------------
// LOCAL VARIABLE :
// offset ( the position in this->atomlink , this is
// exactly what adjacentSet record!)
//----------------------------------------------------------
		const int offset_i = as->offset[i];
		AdjacentSet::getBox(as->box[i], box[i].x, box[i].y, box[i].z);
//----------------------------------------------------------
// USE MEMBER VARIABEL :
// ntype : get the adjacent atom type index
// natom : get the adjacent atom index in this type
//----------------------------------------------------------
		ntype[i] = this->atomlink[offset_i].fatom.getType();
		natom[i] = this->atomlink[offset_i].fatom.getNatom();

		if (test_grid_driver > 1)
		{
			//if(i<100)
			{
/*
				GlobalV::ofs_running << " "
				<< std::setw(5) << box[i].x
				<< std::setw(5) << box[i].y
				<< std::setw(5) << box[i].z
				<< std::setw(8) << offset_i
				<< std::setw(5) << ntype[i]
				<< std::setw(8) << natom[i];
*/
			}
		}

		if (expand_flag)
		{
//----------------------------------------------------------
// EXPLAIN : If expand case, all adjacent site is recorded
// in store_x, store_y, store_z, so we don't need to consider
// box_x, box_y, box_z
// CALL MEMBER FUNCTION :
// NAME : Distance
//----------------------------------------------------------
			adjacent_tau[i].x = this->atomlink[offset_i].fatom.x();
			adjacent_tau[i].y = this->atomlink[offset_i].fatom.y();
			adjacent_tau[i].z = this->atomlink[offset_i].fatom.z();

			if (test_grid_driver > 1)
			{
				//if(i<100)
				{
					const double distance = Distance(this->atomlink[offset], this->atomlink[offset_i]);
/*
					GlobalV::ofs_running << std::setw(10) << this->atomlink[offset_i].fatom.x()
					<< std::setw(10) << this->atomlink[offset_i].fatom.y()
					<< std::setw(10) << this->atomlink[offset_i].fatom.z()
					<< std::setw(10) << distance << std::endl;
*/
					assert(distance <= this->sradius);
				}
			}
		}
		else
		{
//----------------------------------------------------------
// EXPLAIN : If no expand case, first calculate adjacent
// site according to box_x,box_y,box_z, in this case,
// adjacent atoms will be in one of 27 grid.
// CALL MEMBER FUNCTION :
// NAME : Calculate_adjacent_site
// NAME : Distance
//----------------------------------------------------------
			adjacent_tau[i] = Calculate_adjacent_site(offset_i,
			                  vec1[0], vec2[0], vec3[0],
			                  vec1[1], vec2[1], vec3[1],
			                  vec1[2], vec2[2], vec3[2],
			                  box[i].x, box[i].y, box[i].z);

			if (test_grid_driver > 1)
			{
				if(i<100)
				{
					const double distance = this->Distance(this->atomlink[offset], adjacent_tau[i]);
/*
					GlobalV::ofs_running << std::setw(10) << adjacent_tau[i].x
					<< std::setw(10) << adjacent_tau[i].y
					<< std::setw(10) << adjacent_tau[i].z
					<< std::setw(10) << distance;
*/
					assert(distance <= this->sradius);
				}
			}
		}//end if expand_flag
	}// end adjacent number

	return;
}

//==========================================================
// For expand case
//==========================================================
double Grid_Driver::Distance(const AtomLink& a1, const AtomLink& a2)const
{
	const double dx = a1.fatom.x() - a2.fatom.x();
	const double dy = a1.fatom.y() - a2.fatom.y();
	const double dz = a1.fatom.z() - a2.fatom.z();
	return sqrt(dx*dx + dy*dy + dz*dz);
}

//==========================================================
// For not_expand case
//==========================================================
double Grid_Driver::Distance(const AtomLink& a1, const ModuleBase::Vector3<double> &adjacent_site)const
{
	const double dx = a1.fatom.x() - adjacent_site.x;
	const double dy = a1.fatom.y() - adjacent_site.y;
	const double dz = a1.fatom.z() - adjacent_site.z;
	return sqrt(dx*dx + dy*dy + dz*dz);
}


ModuleBase::Vector3<double> Grid_Driver::Calculate_adjacent_site
(
    const short offset,
    const double &box11, const double &box12, const double &box13,
    const double &box21, const double &box22, const double &box23,
    const double &box31, const double &box32, const double &box33,
    const short box_x, const short box_y, const short box_z
)const
{
	ModuleBase::Vector3<double> adjacent_site(0, 0, 0);
	adjacent_site.x = this->atomlink[offset].fatom.x() +
	                  box_x * box11 + box_y * box12 + box_z * box13;

	adjacent_site.y = this->atomlink[offset].fatom.y() +
	                  box_x * box21 + box_y * box22 + box_z * box23;

	adjacent_site.z = this->atomlink[offset].fatom.z() +
	                  box_x * box31 + box_y * box32 + box_z * box33;

	if (test_grid_driver > 3)
	{
/*
		GlobalV::ofs_running << " Offset_i cartesian :" 
		<< " " << this->atomlink[offset].fatom.x()
		<< " " << this->atomlink[offset].fatom.y()
		<< " " << this->atomlink[offset].fatom.z() 
		<< std::endl;
		GlobalV::ofs_running << " box1 : " << box11 << " " << box12 << " " << box13 << std::endl;
		GlobalV::ofs_running << " box2 : " << box21 << " " << box22 << " " << box23 << std::endl;
		GlobalV::ofs_running << " box3 : " << box31 << " " << box32 << " " << box33 << std::endl;
		GlobalV::ofs_running << " adjacent site:" 
		<< " " << adjacent_site.x
		<< " " << adjacent_site.y
		<< " " << adjacent_site.z 
		<< std::endl;
*/
	}

	return adjacent_site;
}

AdjacentAtomInfo Grid_Driver::get_adjs(const UnitCell& ucell_in, const size_t &iat)
{
    const int it = ucell_in.iat2it[iat];
    const int ia = ucell_in.iat2ia[iat];
    const ModuleBase::Vector3<double> &tau = ucell_in.atoms[it].tau[ia];

    AdjacentAtomInfo adjs;
    this->Find_atom(ucell_in, tau, it, ia, &adjs);
    return adjs;
}

std::vector<AdjacentAtomInfo> Grid_Driver::get_adjs(const UnitCell& ucell_in)
{
    std::vector<AdjacentAtomInfo> adjs(ucell_in.nat);
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(size_t iat = 0; iat < ucell_in.nat; iat++)
    {
        adjs[iat] = Grid_Driver::get_adjs(ucell_in, iat);
    }
    return adjs;
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
