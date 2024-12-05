#include "sltk_grid.h"

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

Grid::Grid(const int& test_grid_in) : test_grid(test_grid_in)
{
}

Grid::~Grid()
{
}

void Grid::init(std::ofstream& ofs_in, const UnitCell& ucell, const double radius, bool pbc_in)
{
    ModuleBase::TITLE("SLTK_Grid", "init");
    this->pbc = pbc_in;
    this->sradius2 = radius * radius;
    this->sradius = radius;

    this->x_min = ucell.atoms[0].tau[0].x;
    this->y_min = ucell.atoms[0].tau[0].y;
    this->z_min = ucell.atoms[0].tau[0].z;
    this->x_max = ucell.atoms[0].tau[0].x;
    this->y_max = ucell.atoms[0].tau[0].y;
    this->z_max = ucell.atoms[0].tau[0].z;


    // calculate min & max value
    for (int i = 0; i < ucell.ntype; i++)
    {
        for (int j = 0; j < ucell.atoms[i].na; j++)
        {
            x_min = std::min(x_min, ucell.atoms[i].tau[j].x);
            x_max = std::max(x_max, ucell.atoms[i].tau[j].x);
            y_min = std::min(y_min, ucell.atoms[i].tau[j].y);
            y_max = std::max(y_max, ucell.atoms[i].tau[j].y);
            z_min = std::min(z_min, ucell.atoms[i].tau[j].z);
            z_max = std::max(z_max, ucell.atoms[i].tau[j].z);
        }
    }
    this->Check_Expand_Condition(ucell);
    this->setMemberVariables(ofs_in);
    this->Build_Hash_Table(ucell);
    this->setBoundaryAdjacent(ofs_in);
}
void Grid::Check_Expand_Condition(const UnitCell& ucell)
{
    //	ModuleBase::TITLE(GlobalV::ofs_running, "Atom_input", "Check_Expand_Condition");

    if (!pbc)
    {
        return;
    }
    double a23_1 = ucell.latvec.e22 * ucell.latvec.e33 - ucell.latvec.e23 * ucell.latvec.e32;
    double a23_2 = ucell.latvec.e21 * ucell.latvec.e33 - ucell.latvec.e23 * ucell.latvec.e31;
    double a23_3 = ucell.latvec.e21 * ucell.latvec.e32 - ucell.latvec.e22 * ucell.latvec.e31;
    double a23_norm = sqrt(a23_1 * a23_1 + a23_2 * a23_2 + a23_3 * a23_3);
    double extend_v = a23_norm * sradius;
    double extend_d1 = extend_v / ucell.omega * ucell.lat0 * ucell.lat0 * ucell.lat0;
    int extend_d11 = static_cast<int>(extend_d1);
    // 2016-09-05, LiuXh
    if (extend_d1 - extend_d11 > 0.0)
    {
        extend_d11 += 1;
    }

    double a31_1 = ucell.latvec.e32 * ucell.latvec.e13 - ucell.latvec.e33 * ucell.latvec.e12;
    double a31_2 = ucell.latvec.e31 * ucell.latvec.e13 - ucell.latvec.e33 * ucell.latvec.e11;
    double a31_3 = ucell.latvec.e31 * ucell.latvec.e12 - ucell.latvec.e32 * ucell.latvec.e11;
    double a31_norm = sqrt(a31_1 * a31_1 + a31_2 * a31_2 + a31_3 * a31_3);
    double extend_d2 = a31_norm * sradius / ucell.omega * ucell.lat0 * ucell.lat0 * ucell.lat0;
    int extend_d22 = static_cast<int>(extend_d2);
    // 2016-09-05, LiuXh
    if (extend_d2 - extend_d22 > 0.0)
    {
        extend_d22 += 1;
    }

    double a12_1 = ucell.latvec.e12 * ucell.latvec.e23 - ucell.latvec.e13 * ucell.latvec.e22;
    double a12_2 = ucell.latvec.e11 * ucell.latvec.e23 - ucell.latvec.e13 * ucell.latvec.e21;
    double a12_3 = ucell.latvec.e11 * ucell.latvec.e22 - ucell.latvec.e12 * ucell.latvec.e21;
    double a12_norm = sqrt(a12_1 * a12_1 + a12_2 * a12_2 + a12_3 * a12_3);
    double extend_d3 = a12_norm * sradius / ucell.omega * ucell.lat0 * ucell.lat0 * ucell.lat0;
    int extend_d33 = static_cast<int>(extend_d3);
    // 2016-09-05, LiuXh
    if (extend_d3 - extend_d33 > 0.0)
    {
        extend_d33 += 1;
    }

    glayerX = extend_d11 + 1;
    glayerY = extend_d22 + 1;
    glayerZ = extend_d33 + 1;
    glayerX_minus = extend_d11;
    glayerY_minus = extend_d22;
    glayerZ_minus = extend_d33;
}
//==========================================================
// MEMBER FUNCTION :
// NAME : setMemberVariables(read in data from Atom_input)
//==========================================================
void Grid::setMemberVariables(std::ofstream& ofs_in)
{
    ModuleBase::TITLE("SLTK_Grid", "setMemberVariables");

    // mohan add 2010-09-05
    // AdjacentSet::call_times = 0;


    if (test_grid)
    {
        ModuleBase::GlobalFunc::OUT(ofs_in, "PeriodicBoundary", this->pbc);
        ModuleBase::GlobalFunc::OUT(ofs_in, "Radius(unit:lat0)", sradius);
    }

    //----------------------------------------------------------
    // EXPLAIN : (d_minX,d_minY,d_minZ)minimal value of
    // x[] ,y[] , z[]
    //----------------------------------------------------------
    //----------------------------------------------------------
    // set dx, dy, dz
    //----------------------------------------------------------
    this->cell_nx = glayerX + glayerX_minus;
    this->cell_ny = glayerY + glayerY_minus;
    this->cell_nz = glayerZ + glayerZ_minus;
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

    this->true_cell_x = glayerX_minus;
    this->true_cell_y = glayerX_minus;
    this->true_cell_z = glayerX_minus;
}

void Grid::setBoundaryAdjacent(std::ofstream& ofs_in)
{
    this->Construct_Adjacent_expand(true_cell_x, true_cell_y, true_cell_z);
}

void Grid::Build_Hash_Table(const UnitCell& ucell)
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

    for (int ix = -glayerX_minus; ix < glayerX; ix++)
    {
        for (int iy = -glayerY_minus; iy < glayerY; iy++)
        {
            for (int iz = -glayerZ_minus; iz < glayerZ; iz++)
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
                        a = atom.getCellX() + true_cell_x;
                        b = atom.getCellY() + true_cell_y;
                        c = atom.getCellZ() + true_cell_z;

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

    for (auto& atom_vector: this->Cell[true_i][true_j][true_k].atom_map)
    {
        for (auto& fatom: atom_vector)
        {
            Construct_Adjacent_expand_periodic(true_i, true_j, true_k, fatom);
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


void Grid::Construct_Adjacent_final(const int i,
                                    const int j,
                                    const int k,
                                    FAtom& fatom1,
                                    const int i2,
                                    const int j2,
                                    const int k2,
                                    FAtom& fatom2)
{
    const double x = fatom1.x();
    const double y = fatom1.y();
    const double z = fatom1.z();
    double x2 = fatom2.x();
    double y2 = fatom2.y();
    double z2 = fatom2.z();

    double delta_x = x - x2;
    double delta_y = y - y2;
    double delta_z = z - z2;

    double dr = delta_x * delta_x + delta_y * delta_y + delta_z * delta_z;

    if (dr != 0.0 && dr <= this->sradius2)
    {
        fatom1.addAdjacent(fatom2);
    }
}

