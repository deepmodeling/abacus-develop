#include "sltk_grid.h"
#include "sltk_atom_input.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/memory.h"
#include "module_base/timer.h"

//==================
// Class CellSet
//==================
CellSet::CellSet()
{
	in_grid[0] = 0;
	in_grid[1] = 0;
	in_grid[2] = 0;
}


const char* const Grid::ERROR[3] =
{
	"Function Grid<Input>::buildCache: exception std::logic_error\n"
	"\tThe max length of input must be a positive number!",
	"Function Grid<Input>::buildCache: exception std::out_of_range\n"
	"\tLogic error! The atom amount is above the maxAmount of input!",
	"Function Grid<Input>::foldHashTable::AtomLinkPointStack::push: exception std::out_of_range\n"
	"\tLogic error! The atom amount in one grid must be less then the value of "
	"MAX_ATOM_IN_ONE_GRID!"
};

Grid::Grid(const int &test_grid_in):test_grid(test_grid_in)
{
//	ModuleBase::TITLE("Grid","Grid");
//----------------------------------------------------------
// EXPLAIN : init_cell_flag (use this flag in case memory
// leak)
//----------------------------------------------------------
	init_cell_flag = false;
}

Grid::~Grid()
{
	this->delete_Cell();
}

void Grid::init(std::ofstream &ofs_in,
				const UnitCell &ucell, 
				const Atom_input &input)
{
	ModuleBase::TITLE("SLTK_Grid", "init");

	this->setMemberVariables(ofs_in, input);
	this->Build_Hash_Table(ucell, input);
	this->setBoundaryAdjacent(ofs_in, input);

	return;
}

//==========================================================
// MEMBER FUNCTION :
// NAME : setMemberVariables(read in data from Atom_input)
//==========================================================
void Grid::setMemberVariables(
	std::ofstream &ofs_in, //  output data to ofs
	const Atom_input &input)
{
	ModuleBase::TITLE("SLTK_Grid", "setMemberVariables");

	this->delete_Cell();
	// mohan add 2010-09-05
	AdjacentSet::call_times = 0;

	this->natom = input.getAmount();
	if(test_grid)ModuleBase::GlobalFunc::OUT(ofs_in, "natom", natom);

	this->pbc = input.getBoundary();
	if(test_grid)ModuleBase::GlobalFunc::OUT(ofs_in, "PeriodicBoundary", this->pbc);

	this->sradius = input.getRadius();
	if(test_grid)ModuleBase::GlobalFunc::OUT(ofs_in, "Radius(unit:lat0)", sradius);

	for (int i = 0;i < 3;i++)
	{
		this->vec1[i] = input.vec1[i];
		this->vec2[i] = input.vec2[i];
		this->vec3[i] = input.vec3[i];
	}

	this->lat_now = input.getLatNow();
	if(test_grid)ModuleBase::GlobalFunc::OUT(ofs_in,"lat0(unit:Bohr)", lat_now);

	this->expand_flag = input.getExpandFlag();
	if(test_grid)ModuleBase::GlobalFunc::OUT(ofs_in,"Expand_flag", expand_flag);
	
	// output std::vector
	if(test_grid)ModuleBase::GlobalFunc::OUT(ofs_in,"Vec1",vec1[0],vec1[1],vec1[2]);
	if(test_grid)ModuleBase::GlobalFunc::OUT(ofs_in,"Vec2",vec2[0],vec2[1],vec2[2]);
	if(test_grid)ModuleBase::GlobalFunc::OUT(ofs_in,"Vec3",vec3[0],vec3[1],vec3[2]);

	// output grid length
	this->grid_length[0] = input.Clength0();
	this->grid_length[1] = input.Clength1();
	this->grid_length[2] = input.Clength2();

	if(test_grid)ModuleBase::GlobalFunc::OUT(ofs_in,"Grid_length",grid_length[0],grid_length[1],grid_length[2]);
//----------------------------------------------------------
// EXPLAIN : (d_minX,d_minY,d_minZ)minimal value of
// x[] ,y[] , z[]
//----------------------------------------------------------
	this->d_minX = input.minX();
	this->d_minY = input.minY();
	this->d_minZ = input.minZ();
	if(test_grid)ModuleBase::GlobalFunc::OUT(ofs_in,"MinCoordinate",d_minX,d_minY,d_minZ);
//----------------------------------------------------------
//layer: grid layer after expand
//----------------------------------------------------------
	this->cell_x_length = input.getCellXLength();
	this->cell_y_length = input.getCellYLength();
	this->cell_z_length = input.getCellZLength();
	if(test_grid)ModuleBase::GlobalFunc::OUT(ofs_in,"CellLength(unit: lat0)",cell_x_length,cell_y_length,cell_z_length);
//----------------------------------------------------------
// set dx, dy, dz
//----------------------------------------------------------
	this->dx = input.getCellX();
	this->dy = input.getCellY();
	this->dz = input.getCellZ();
	if(test_grid)ModuleBase::GlobalFunc::OUT(ofs_in,"CellNumber",dx,dy,dz);


	Cell = new CellSet**[dx];
	for (int i = 0;i < dx;i++)
	{
		Cell[i] = new CellSet*[dy];

		for (int j = 0;j < dy;j++)
		{
			Cell[i][j] = new CellSet[dz];
		}
	}
	this->init_cell_flag = true;

	this->true_cell_x = input.getGrid_layerX_minus();
	this->true_cell_y = input.getGrid_layerY_minus();
	this->true_cell_z = input.getGrid_layerZ_minus();

	return;
}

void Grid::setBoundaryAdjacent(
	std::ofstream &ofs_in,
	const Atom_input &input)
{
	if (test_grid) ModuleBase::TITLE(ofs_in, "Grid", "setBoundaryAdjacent");

	if (expand_flag)
	{
		const int i = input.getGrid_layerX_minus();
		const int j = input.getGrid_layerY_minus();
		const int k = input.getGrid_layerZ_minus();
		this->Construct_Adjacent_expand(i, j, k);
	}
	else
	{
		this->Construct_Adjacent_begin();
	}

	if(test_grid)ModuleBase::GlobalFunc::OUT(ofs_in,"Adjacent_set call times",AdjacentSet::call_times);

	//if (test_grid)ModuleBase::GlobalFunc::DONE(ofs_in, "Construct_Adjacent");

	return;
}

#include "module_base/mathzone.h"
void Grid::In_Which_Cell(const UnitCell &ucell, int &a, int &b, int &c, const FAtom &atom)const
{
	if (expand_flag)
	{
		// EXPLAIN : In expand grid case,
		// the input cell is exactly the same as input file.
		a = atom.getCellX() - this->d_minX;
		b = atom.getCellY() - this->d_minY;
		c = atom.getCellZ() - this->d_minZ;
	}
	else
	{
//----------------------------------------------------------
// EXPLAIN : Not expand case , the cell is 'cubic',
// the three dimension length :
// cell_x_length = |radius|
// cell_y_length = |radius|
// cell_z_length = |radius|
//
// So we don't need crystal coordinate to locate the atom.
// We use cartesian coordinate directly.
//----------------------------------------------------------
		a = static_cast<int>(std::floor((atom.x() - this->d_minX) / this->cell_x_length));
		b = static_cast<int>(std::floor((atom.y() - this->d_minY) / this->cell_y_length));
		c = static_cast<int>(std::floor((atom.z() - this->d_minZ) / this->cell_z_length));

		/*
		std::cout<<"\n"<<std::setw(12)<<atom.x()
		<<std::setw(12)<<atom.y()
		<<std::setw(12)<<atom.z()
		<<std::setw(6)<<a
		<<std::setw(6)<<b
		<<std::setw(6)<<c
		<<std::setw(12)<<d_minX
		<<std::setw(12)<<d_minY
		<<std::setw(12)<<d_minZ;
		*/
	}

	return;
}

void Grid::Build_Hash_Table(const UnitCell &ucell, const Atom_input &input)
{
    ModuleBase::timer::tick("Grid", "Build_Hash_Table");

	for (int i = 0; i < input.getAmount(); ++i)
	{
		const FAtom &atom = input.getFakeAtom(i);
		all_atom_map[{atom.x(), atom.y(), atom.z()}] = atom;

		int a, b, c;
		this->In_Which_Cell(ucell, a, b, c, atom);
		this->Cell[a][b][c].atom_map[{atom.x(), atom.y(), atom.z()}] = atom;
	}
    ModuleBase::timer::tick("Grid", "Build_Hash_Table");
}

void Grid::Construct_Adjacent_expand(
	const int true_i, 
	const int true_j, 
	const int true_k)
{
//	if (test_grid)ModuleBase::TITLE(ofs_running, "Grid", "Construct_Adjacent_expand");

//----------------------------------------------------------
// EXPlAIN : In expand grid case, use
// AdjacentSet::index_expand() to record the grid number,
// We use formula (i*dy*dz + j*dz + k) to store the
// displacement of cell.
// Of course , an alternative operatiion is to store the
// (i,j,k),but we want to use memory as small as possible
// for this storage.
//----------------------------------------------------------
    ModuleBase::timer::tick("Grid", "Construct_Adjacent_expand");

	AdjacentSet::setExpandFlag(this->expand_flag);

	AdjacentSet::setDx(this->dx);

	AdjacentSet::setDy(this->dy);

	AdjacentSet::setDz(this->dz);

	// mohan add 2009-10-20
	AdjacentSet::setTrueX(true_i);
	
	AdjacentSet::setTrueY(true_j);
	
	AdjacentSet::setTrueZ(true_k);
	

	AdjacentSet::setCenter(true_i * dy * dz + true_j * dz + true_k);
	

//	if(test_grid)OUT(ofs_running,"GridCenter",true_i,true_j,true_k);
//	if(test_grid)OUT(ofs_running,"GridDim",dx,dy,dz);

//-----------------------------------------------------------
// EXPLAIN : (true_i,true_j,true_k) is the cell we want
// to found AdjacentSet.And other cell save the displacement
// of center_grid in 'in_grid'
//-----------------------------------------------------------
	for (int i = 0;i < this->dx;i++)
	{
		for (int j = 0;j < this->dy;j++)
		{
			for (int k = 0;k < this->dz;k++)
			{
				this->Cell[i][j][k].in_grid[0] = i - true_i;
				this->Cell[i][j][k].in_grid[1] = j - true_j;
				this->Cell[i][j][k].in_grid[2] = k - true_k;
			}
		}
	}

//----------------------------------------------------------
// EXPLAIN : Only construct AdjacentSet for 'true' cell.
//----------------------------------------------------------
    for (auto& pair : this->Cell[true_i][true_j][true_k].atom_map)
	{
		FAtom& fatom = pair.second;

		if (this->pbc)
		{
			Construct_Adjacent_expand_periodic(true_i, true_j, true_k, fatom);
			// std::cout << "fatom1 = " << fatom.getNatom() << "  " << fatom.getAdjacent().size() << std::endl;
		}
		else
		{
			ModuleBase::WARNING_QUIT("Construct_Adjacent_expand", "\n Expand case, must use periodic boundary.");
		}
	}
    ModuleBase::timer::tick("Grid", "Construct_Adjacent_expand");
}

void Grid::Construct_Adjacent_expand_periodic(const int true_i, 
												const int true_j, 
												const int true_k, 
												FAtom& fatom)
{
//	if (test_grid)ModuleBase::TITLE(ofs_running, "Grid", "Construct_Adjacent_expand_periodic");
    ModuleBase::timer::tick("Grid", "Construct_Adjacent_expand_periodic");

	for (int i = 0;i < this->dx;i++)
	{
		for (int j = 0;j < this->dy;j++)
		{
			for (int k = 0;k < this->dz;k++)
			{
			    for (auto& pair : this->Cell[i][j][k].atom_map)
				{
					FAtom& fatom2 = pair.second;
					Construct_Adjacent_final(true_i, true_j, true_k, fatom, i, j, k, fatom2);
				}
			}
		}
	}
    ModuleBase::timer::tick("Grid", "Construct_Adjacent_expand_periodic");

}

void Grid::Construct_Adjacent_begin(void)
{
//	if (test_grid)ModuleBase::TITLE(ofs_running, "Grid", "Construct_Adjacent_begin");

//----------------------------------------------------------
// EXPLAIN : Searching in all cells in this grid
//----------------------------------------------------------

	for (int i = 0;i < this->dx;i++)
	{
		for (int j = 0;j < this->dy;j++)
		{
			for (int k = 0;k < this->dz;k++)
			{
//----------------------------------------------------------
// EXPLAIN : Cell length == Number of atoms in this cell
//----------------------------------------------------------
			    for (auto& pair : this->Cell[i][j][k].atom_map)
				{
					FAtom& fatom = pair.second;
					//pbc: periodic boundary condition
					if (this->pbc)
					{
						Construct_Adjacent_periodic(i, j, k, fatom);
					}
					else
					{
						Construct_Adjacent_nature(i, j, k, fatom);
					}
					
				}//ia
			}//k
		}//j
	}//i

	return;
}

void Grid::Construct_Adjacent_nature
(
    const int i,
    const int j,
    const int k,
    FAtom & fatom1
)
{
//	if(test_grid)ModuleBase::TITLE(ofs_running,"Grid","Construct_Adjacent_nature");
	for (int i2 = i - 1;i2 <= i + 1;i2++)
	{
		if (i2<dx && i2 >= 0)
		{
			for (int j2 = j - 1;j2 <= j + 1;j2++)
			{
				if (j2<dy && j2 >= 0)
				{
					for (int k2 = k - 1;k2 <= k + 1;k2++)
					{
						if (k2<dz && k2 >= 0)
						{
							for (auto& pair : this->Cell[i2][j2][k2].atom_map)
							{
								FAtom& fatom2 = pair.second;
								Construct_Adjacent_final(i, j, k, fatom1, i2, j2, k2, fatom2);
							}//ia2
						}
					}//k2
				}
			}//j2
		}
	}//2

	return;
}

void Grid::Construct_Adjacent_periodic
(
    const int i,
    const int j,
    const int k,
    FAtom & fatom1
)
{
//	if(test_grid)ModuleBase::TITLE(ofs_running,"Grid","Construct_Adjacent_periodic");
	bool first_i = true;

	for (int i2 = i - 1;i2 <= i + 1;i2++)
	{
		bool first_j = true;

		for (int j2 = j - 1;j2 <= j + 1;j2++)
		{
			bool first_k = true;

			for (int k2 = k - 1;k2 <= k + 1;k2++)
			{
				int temp_i = i2;
				int temp_j = j2;
				int temp_k = k2;

				int g0 = 0;
				int g1 = 0;
				int g2 = 0;

				if (i2 < 0)
				{
					g0 = -1;

					if (first_i)
					{
						if (dx >= 2)
						{
							i2--;
							temp_i--;
						}

						first_i = false;
					}

					i2 += dx;
				}
				else if (i2 >= dx)
				{
					g0 = 1;
					i2 -= dx;
				}

				if (j2 < 0)
				{
					g1 = -1;

					if (first_j)
					{
						if (dy >= 2)
						{
							j2--;
							temp_j--;
						}

						first_j = false;
					}

					j2 += dy;
				}
				else if (j2 >= dy)
				{
					g1 = 1;
					j2 -= dy;
				}

				if (k2 < 0)
				{
					g2 = -1;

					if (first_k)
					{
						if (dz >= 2)
						{
							k2--;
							temp_k--;
						}

						first_k = false;
					}

					k2 += dz;
				}
				else if (k2 >= dz)
				{
					g2 = 1;
					k2 -= dz;
				}

				Cell[i2][j2][k2].in_grid[0] = g0;

				Cell[i2][j2][k2].in_grid[1] = g1;
				Cell[i2][j2][k2].in_grid[2] = g2;

				for (auto& pair : this->Cell[i2][j2][k2].atom_map)
				{
					FAtom& fatom2 = pair.second;
					Construct_Adjacent_final(i, j, k, fatom1, i2, j2, k2, fatom2);
				}//ia2

				i2 = temp_i;

				j2 = temp_j;

				k2 = temp_k;//resume i2 j2 k2
			}//k2
		}//j2
	}//i2

	return;
}

void Grid::Construct_Adjacent_final(const int i, const int j, const int k, FAtom & fatom1,
 									const int i2, const int j2, const int k2, FAtom & fatom2)
{
	ModuleBase::timer::tick("Grid_Driver","Construct_Adjacent_final");
//----------------------------------------------------------
// EXPLAIN : 		expand_case				not_expand_case
// (i,j,k,ia) 		only the 'true' cell	only the 'true' grid
// (i2,j2,k2,ia2) 	all atoms in grid		all atoms in 27*cell
//----------------------------------------------------------
// (suitable for small cell periodic condition)
// Expand_Case : many 'pseudo' cells, only one true cell,
// one grid(true grid).
// Advantage : only the atoms in 'true' cell need to construct
// AdjacentSet.
// Disadvantage : must search all atoms in true grid to construct
// AdjacentSet.
//
// (suitable for large cell periodic/nature condition,here
// we discuss periodic case,once you known this case, nature
// boundary is easy to understand)
// Not_Expand_Case : 27 'pseudo' grid,only one true grid,
// many true cells.
// Advantage : (the disadvantage above is the advantage now)
// only need to search 27*cells to construct AdjacentSet
// for each cell.
// Disadvantage : (the advantave mentioned above)
// need to construct adjacent for each cell.
//----------------------------------------------------------
	const double x  = fatom1.x();
	const double y  = fatom1.y();
	const double z  = fatom1.z();
	double x2 = fatom2.x();
	double y2 = fatom2.y();
	double z2 = fatom2.z();
//----------------------------------------------------------
// EXPLAIN : in different case 'in_grid' has different
// meaning.
//----------------------------------------------------------
// NAME : 			expand_case		 |  not_expand_case
// in_which_grid	'not available'	 |  one of 27 adjacent grid
// in_which_cell	one of all cells |  'not available'
//----------------------------------------------------------
// The solution here is we save these datas in one structrue
// named : 'in_grid'
//----------------------------------------------------------
	const int b0 = Cell[i2][j2][k2].in_grid[0];
	const int b1 = Cell[i2][j2][k2].in_grid[1];
	const int b2 = Cell[i2][j2][k2].in_grid[2];

	if (!expand_flag)
	{
		x2 = x2 + b0 * vec1[0] + b1 * vec2[0] + b2 * vec3[0];
		y2 = y2 + b0 * vec1[1] + b1 * vec2[1] + b2 * vec3[1];
		z2 = z2 + b0 * vec1[2] + b1 * vec2[2] + b2 * vec3[2];
	}

//----------------------------------------------------------
// EXPlAIN : Calculate distance between two atoms.
//----------------------------------------------------------
	double delta_x = x - x2;
	double delta_y = y - y2;
	double delta_z = z - z2;

	double dr = sqrt(delta_x * delta_x + delta_y * delta_y + delta_z * delta_z);

	ModuleBase::timer::tick("Grid_Driver","addAdjacent");

	if (dr != 0.0 && dr <= this->sradius)
	{
		fatom1.addAdjacent(fatom2);
	}
	ModuleBase::timer::tick("Grid_Driver","addAdjacent");

	ModuleBase::timer::tick("Grid_Driver","Construct_Adjacent_final");
}
//2015-05-07
void Grid::delete_vector(const Atom_input &input)
{
	if (expand_flag)
	{
		const int i = input.getGrid_layerX_minus();
		const int j = input.getGrid_layerY_minus();
		const int k = input.getGrid_layerZ_minus();
		for (auto &pair : Cell[i][j][k].atom_map)
		{
			if (this->pbc)
			{
				pair.second.clearAdjacent();
			}
		}
	}
}
