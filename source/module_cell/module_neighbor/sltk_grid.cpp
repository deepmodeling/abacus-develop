#include "sltk_grid.h"

#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/memory.h"
#include "module_base/timer.h"
#include "sltk_atom_input.h"

//==================
// Class CellSet
//==================
CellSet::CellSet()
{
    in_grid[0] = 0;
    in_grid[1] = 0;
    in_grid[2] = 0;
}

Grid::Grid(const int& test_grid_in) : test_grid(test_grid_in)
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

void Grid::init(std::ofstream& ofs_in, const UnitCell& ucell, const Atom_input& input)
{
    ModuleBase::TITLE("SLTK_Grid", "init");

    this->setMemberVariables(ofs_in, input);
    this->Build_Hash_Table(ucell, input);
    this->setBoundaryAdjacent(ofs_in, input);
}

//==========================================================
// MEMBER FUNCTION :
// NAME : setMemberVariables(read in data from Atom_input)
//==========================================================
void Grid::setMemberVariables(std::ofstream& ofs_in, //  output data to ofs
                              const Atom_input& input)
{
    ModuleBase::TITLE("SLTK_Grid", "setMemberVariables");

    this->delete_Cell();
    // mohan add 2010-09-05
    // AdjacentSet::call_times = 0;

    this->pbc = input.getBoundary();
    this->sradius2 = input.getRadius() * input.getRadius();
    this->sradius = input.getRadius();
    this->expand_flag = input.getExpandFlag();

    if (test_grid)
    {
        ModuleBase::GlobalFunc::OUT(ofs_in, "PeriodicBoundary", this->pbc);
        ModuleBase::GlobalFunc::OUT(ofs_in, "Radius(unit:lat0)", sradius);
        ModuleBase::GlobalFunc::OUT(ofs_in, "Expand_flag", expand_flag);
    }

    //----------------------------------------------------------
    // EXPLAIN : (d_minX,d_minY,d_minZ)minimal value of
    // x[] ,y[] , z[]
    //----------------------------------------------------------
    this->d_minX = input.minX();
    this->d_minY = input.minY();
    this->d_minZ = input.minZ();
    if (test_grid)
    {
        ModuleBase::GlobalFunc::OUT(ofs_in, "MinCoordinate", d_minX, d_minY, d_minZ);
    }
    //----------------------------------------------------------
    // set dx, dy, dz
    //----------------------------------------------------------
    this->cell_nx = input.getCell_nX();
    this->cell_ny = input.getCell_nY();
    this->cell_nz = input.getCell_nZ();
    if (test_grid)
    {
        ModuleBase::GlobalFunc::OUT(ofs_in, "CellNumber", cell_nx, cell_ny, cell_nz);
    }

    Cell.resize(cell_nx);
    for (int i = 0; i < cell_nx; i++)
    {
        Cell[i].resize(cell_ny);
        for (int j = 0; j < cell_ny; j++)
        {
            Cell[i][j].resize(cell_nz);
        }
    }
    this->init_cell_flag = true;

    this->true_cell_x = input.getGrid_layerX_minus();
    this->true_cell_y = input.getGrid_layerY_minus();
    this->true_cell_z = input.getGrid_layerZ_minus();
}

void Grid::setBoundaryAdjacent(std::ofstream& ofs_in, const Atom_input& input)
{
    if (expand_flag)
    {
        this->Construct_Adjacent_expand(true_cell_x, true_cell_y, true_cell_z);
    }
    else
    {
        this->Construct_Adjacent_begin();
    }
}

void Grid::Build_Hash_Table(const UnitCell& ucell, const Atom_input& input)
{
    ModuleBase::timer::tick("Grid", "Build_Hash_Table");

    // TODO in case expand == false, the following code is over malloc
    for (int i = 0; i < cell_nx; i++)
    {
        for (int j = 0; j < cell_ny; j++)
        {
            for (int k = 0; k < cell_nz; k++)
            {
                Cell[i][j][k].atom_map.resize(ucell.ntype);
                for (int it = 0; it < ucell.ntype; ++it)
                {
                    Cell[i][j][k].atom_map[it].resize(ucell.atoms[it].na);
                }
            }
        }
    }
    ModuleBase::Vector3<double> vec1(ucell.latvec.e11, ucell.latvec.e12, ucell.latvec.e13);
    ModuleBase::Vector3<double> vec2(ucell.latvec.e21, ucell.latvec.e22, ucell.latvec.e23);
    ModuleBase::Vector3<double> vec3(ucell.latvec.e31, ucell.latvec.e32, ucell.latvec.e33);

    for (int ix = -input.getGrid_layerX_minus(); ix < input.getGrid_layerX(); ix++)
    {
        for (int iy = -input.getGrid_layerY_minus(); iy < input.getGrid_layerY(); iy++)
        {
            for (int iz = -input.getGrid_layerZ_minus(); iz < input.getGrid_layerZ(); iz++)
            {
                for (int i = 0; i < ucell.ntype; i++)
                {
                    for (int j = 0; j < ucell.atoms[i].na; j++)
                    {
                        double x = ucell.atoms[i].tau[j].x + vec1[0] * ix + vec2[0] * iy + vec3[0] * iz;
                        double y = ucell.atoms[i].tau[j].y + vec1[1] * ix + vec2[1] * iy + vec3[1] * iz;
                        double z = ucell.atoms[i].tau[j].z + vec1[2] * ix + vec2[2] * iy + vec3[2] * iz;
                        FAtom atom(x, y, z, i, j, ix, iy, iz);
                        int a, b, c;
                        if (expand_flag)
                        {
                            // EXPLAIN : In expand grid case,
                            // the input cell is exactly the same as input file.
                            a = atom.getCellX() + true_cell_x;
                            b = atom.getCellY() + true_cell_y;
                            c = atom.getCellZ() + true_cell_z;
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
                            a = static_cast<int>(std::floor((atom.x() - this->d_minX) / this->sradius));
                            b = static_cast<int>(std::floor((atom.y() - this->d_minY) / this->sradius));
                            c = static_cast<int>(std::floor((atom.z() - this->d_minZ) / this->sradius));
                        }

                        this->Cell[a][b][c].atom_map[atom.getType()][atom.getNatom()] = atom;
                    }
                }
            }
        }
    }
    ModuleBase::timer::tick("Grid", "Build_Hash_Table");
}

void Grid::Construct_Adjacent_expand(const int true_i, const int true_j, const int true_k)
{
    ModuleBase::timer::tick("Grid", "Construct_Adjacent_expand");

    //-----------------------------------------------------------
    // EXPLAIN : (true_i,true_j,true_k) is the cell we want
    // to found AdjacentSet.And other cell save the displacement
    // of center_grid in 'in_grid'
    //-----------------------------------------------------------
    for (int i = 0; i < this->cell_nx; i++)
    {
        for (int j = 0; j < this->cell_ny; j++)
        {
            for (int k = 0; k < this->cell_nz; k++)
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
    for (auto& atom_vector: this->Cell[true_i][true_j][true_k].atom_map)
    {
        for (auto& fatom: atom_vector)
        {
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
    }
    ModuleBase::timer::tick("Grid", "Construct_Adjacent_expand");
}

void Grid::Construct_Adjacent_expand_periodic(const int true_i, const int true_j, const int true_k, FAtom& fatom)
{
    //	if (test_grid)ModuleBase::TITLE(ofs_running, "Grid", "Construct_Adjacent_expand_periodic");
    ModuleBase::timer::tick("Grid", "Construct_Adjacent_expand_periodic");

    for (int i = 0; i < this->cell_nx; i++)
    {
        for (int j = 0; j < this->cell_ny; j++)
        {
            for (int k = 0; k < this->cell_nz; k++)
            {
                for (auto& atom_vector: this->Cell[i][j][k].atom_map)
                {
                    for (auto& fatom2: atom_vector)
                    {
                        Construct_Adjacent_final(true_i, true_j, true_k, fatom, i, j, k, fatom2);
                    }
                }
            }
        }
    }
    ModuleBase::timer::tick("Grid", "Construct_Adjacent_expand_periodic");
}

void Grid::Construct_Adjacent_begin()
{
    //	if (test_grid)ModuleBase::TITLE(ofs_running, "Grid", "Construct_Adjacent_begin");

    //----------------------------------------------------------
    // EXPLAIN : Searching in all cells in this grid
    //----------------------------------------------------------

    for (int i = 0; i < this->cell_nx; i++)
    {
        for (int j = 0; j < this->cell_ny; j++)
        {
            for (int k = 0; k < this->cell_nz; k++)
            {
                //----------------------------------------------------------
                // EXPLAIN : Cell length == Number of atoms in this cell
                //----------------------------------------------------------
                for (auto& atom_vector: this->Cell[i][j][k].atom_map)
                {
                    for (auto& fatom2: atom_vector)
                    {
                        // pbc: periodic boundary condition
                        if (this->pbc)
                        {
                            Construct_Adjacent_periodic(i, j, k, fatom2);
                        }
                        else
                        {
                            Construct_Adjacent_nature(i, j, k, fatom2);
                        }
                    }

                } // ia
            } // k
        } // j
    } // i

    return;
}

void Grid::Construct_Adjacent_nature(const int i, const int j, const int k, FAtom& fatom1)
{
    //	if(test_grid)ModuleBase::TITLE(ofs_running,"Grid","Construct_Adjacent_nature");
    for (int i2 = i - 1; i2 <= i + 1; i2++)
    {
        if (i2 < cell_nx && i2 >= 0)
        {
            for (int j2 = j - 1; j2 <= j + 1; j2++)
            {
                if (j2 < cell_ny && j2 >= 0)
                {
                    for (int k2 = k - 1; k2 <= k + 1; k2++)
                    {
                        if (k2 < cell_nz && k2 >= 0)
                        {
                            for (auto& atom_vector: this->Cell[i2][j2][k2].atom_map)
                            {
                                for (auto& fatom2: atom_vector)
                                {
                                    Construct_Adjacent_final(i, j, k, fatom1, i2, j2, k2, fatom2);
                                } // ia2
                            }
                        }
                    } // k2
                }
            } // j2
        }
    } // 2

    return;
}

void Grid::Construct_Adjacent_periodic(const int i, const int j, const int k, FAtom& fatom1)
{
    //	if(test_grid)ModuleBase::TITLE(ofs_running,"Grid","Construct_Adjacent_periodic");
    bool first_i = true;

    for (int i2 = i - 1; i2 <= i + 1; i2++)
    {
        bool first_j = true;

        for (int j2 = j - 1; j2 <= j + 1; j2++)
        {
            bool first_k = true;

            for (int k2 = k - 1; k2 <= k + 1; k2++)
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
                        if (cell_nx >= 2)
                        {
                            i2--;
                            temp_i--;
                        }

                        first_i = false;
                    }

                    i2 += cell_nx;
                }
                else if (i2 >= cell_nx)
                {
                    g0 = 1;
                    i2 -= cell_nx;
                }

                if (j2 < 0)
                {
                    g1 = -1;

                    if (first_j)
                    {
                        if (cell_ny >= 2)
                        {
                            j2--;
                            temp_j--;
                        }

                        first_j = false;
                    }

                    j2 += cell_ny;
                }
                else if (j2 >= cell_ny)
                {
                    g1 = 1;
                    j2 -= cell_ny;
                }

                if (k2 < 0)
                {
                    g2 = -1;

                    if (first_k)
                    {
                        if (cell_nz >= 2)
                        {
                            k2--;
                            temp_k--;
                        }

                        first_k = false;
                    }

                    k2 += cell_nz;
                }
                else if (k2 >= cell_nz)
                {
                    g2 = 1;
                    k2 -= cell_nz;
                }

                Cell[i2][j2][k2].in_grid[0] = g0;

                Cell[i2][j2][k2].in_grid[1] = g1;
                Cell[i2][j2][k2].in_grid[2] = g2;

                for (auto& atom_vector: this->Cell[i2][j2][k2].atom_map)
                {
                    for (auto& fatom2: atom_vector)
                    {
                        Construct_Adjacent_final(i, j, k, fatom1, i2, j2, k2, fatom2);
                    } // ia2
                }

                i2 = temp_i;

                j2 = temp_j;

                k2 = temp_k; // resume i2 j2 k2
            } // k2
        } // j2
    } // i2

    return;
}

void Grid::Construct_Adjacent_final(const int i,
                                    const int j,
                                    const int k,
                                    FAtom& fatom1,
                                    const int i2,
                                    const int j2,
                                    const int k2,
                                    FAtom& fatom2)
{
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
    const double x = fatom1.x();
    const double y = fatom1.y();
    const double z = fatom1.z();
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

    //----------------------------------------------------------
    // EXPlAIN : Calculate distance between two atoms.
    //----------------------------------------------------------
    double delta_x = x - x2;
    double delta_y = y - y2;
    double delta_z = z - z2;

    double dr = delta_x * delta_x + delta_y * delta_y + delta_z * delta_z;

    if (dr != 0.0 && dr <= this->sradius2)
    {
        fatom1.addAdjacent(fatom2);
    }
}
// 2015-05-07
void Grid::delete_vector(int i, int j, int k)
{
    if (expand_flag)
    {
        if (this->pbc)
        {
            this->Cell[i][j][k].atom_map.clear();
        }
    }
}
