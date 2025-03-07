#include "sltk_grid.h"

#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/memory.h"
#include "module_base/timer.h"

Grid::Grid(const int& test_grid_in) : test_grid(test_grid_in)
{
}

Grid::~Grid()
{
    this->clear_atoms();
}

void Grid::init(std::ofstream& ofs_in, const UnitCell& ucell, const double radius_in, const bool boundary)
{
    ModuleBase::TITLE("SLTK_Grid", "init");
    ModuleBase::timer::tick("atom_arrange", "grid_d.init");
    this->pbc = boundary;
    this->sradius2 = radius_in * radius_in;
    this->sradius = radius_in;

    ModuleBase::GlobalFunc::OUT(ofs_in, "PeriodicBoundary", this->pbc);
    ModuleBase::GlobalFunc::OUT(ofs_in, "Radius(unit:lat0)", sradius);

    this->Check_Expand_Condition(ucell);
    ModuleBase::GlobalFunc::OUT(ofs_in, "glayer", glayerX, glayerY, glayerZ);
    ModuleBase::GlobalFunc::OUT(ofs_in, "glayer_minus", glayerX_minus, glayerY_minus, glayerZ_minus);

    this->setMemberVariables(ofs_in, ucell);
    this->Construct_Adjacent(ucell);
    ModuleBase::timer::tick("atom_arrange", "grid_d.init");
}

void Grid::Check_Expand_Condition(const UnitCell& ucell)
{
    //	ModuleBase::TITLE(GlobalV::ofs_running, "Atom_input", "Check_Expand_Condition");

    if (!pbc)
    {
        return;
    }

    /*2016-07-19, LiuXh
        // the unit of extent_1DX,Y,Z is lat0.
        // means still how far can be included now.
        double extent_1DX = glayerX * clength0 - dmaxX;
        while (radius > extent_1DX)
        {
            glayerX++;
            extent_1DX = glayerX * clength0 - dmaxX;
        }
        double extent_1DY = glayerY * clength1 - dmaxY;
        while (radius > extent_1DY)
        {
            glayerY++;
            extent_1DY = glayerY * clength1 - dmaxY;
        }
        double extent_1DZ = glayerZ * clength2 - dmaxZ;
        while (radius > extent_1DZ)
        {
            glayerZ++;
            extent_1DZ = glayerZ * clength2 - dmaxZ;
        }

        // in case the cell is not retangle.
        // mohan added 2009-10-23
        // if this is not added, it's a serious bug.
        glayerX++;
        glayerY++;
        glayerZ++;
        if(test_atom_input)
        {
            GlobalV::ofs_running << " Extend distance from the (maxX,maxY,maxZ) direct position in this unitcell: " <<
    std::endl;
        }

        if(test_atom_input)OUT(GlobalV::ofs_running,"ExtentDim+",extent_1DX,extent_1DY,extent_1DZ);

        double extent_1DX_minus = glayerX_minus * clength0 + dminX;
        while (radius > extent_1DX_minus)
        {
            glayerX_minus++;
            extent_1DX_minus = glayerX_minus * clength0 + dminX;
        }
        double extent_1DY_minus = glayerY_minus * clength1 + dminY;
        while (radius > extent_1DY_minus)
        {
            glayerY_minus++;
            extent_1DY_minus = glayerY_minus * clength1 + dminY;
        }
        double extent_1DZ_minus = glayerZ_minus * clength2 + dminZ;
        while (radius > extent_1DZ_minus)
        {
            glayerZ_minus++;
            extent_1DZ_minus = glayerZ_minus * clength2 + dminZ;
        }

        // in case the cell is not retangle.
        // mohan added 2009-10-23
        // if this is not added, it's a serious bug.
        glayerX_minus++;
        glayerY_minus++;
        glayerZ_minus++;

        //glayerX_minus++;
        //glayerY_minus++;
        //glayerZ_minus++;
    2016-07-19, LiuXh*/
    // Begin, 2016-07-19, LiuXh
    double a23_1 = ucell.latvec.e22 * ucell.latvec.e33 - ucell.latvec.e23 * ucell.latvec.e32;
    double a23_2 = ucell.latvec.e21 * ucell.latvec.e33 - ucell.latvec.e23 * ucell.latvec.e31;
    double a23_3 = ucell.latvec.e21 * ucell.latvec.e32 - ucell.latvec.e22 * ucell.latvec.e31;
    double a23_norm = sqrt(a23_1 * a23_1 + a23_2 * a23_2 + a23_3 * a23_3);
    double extend_v = a23_norm * sradius;
    double extend_d1 = extend_v / ucell.omega * ucell.lat0 * ucell.lat0 * ucell.lat0;
    int extend_d11 = std::ceil(extend_d1);

    double a31_1 = ucell.latvec.e32 * ucell.latvec.e13 - ucell.latvec.e33 * ucell.latvec.e12;
    double a31_2 = ucell.latvec.e31 * ucell.latvec.e13 - ucell.latvec.e33 * ucell.latvec.e11;
    double a31_3 = ucell.latvec.e31 * ucell.latvec.e12 - ucell.latvec.e32 * ucell.latvec.e11;
    double a31_norm = sqrt(a31_1 * a31_1 + a31_2 * a31_2 + a31_3 * a31_3);
    double extend_d2 = a31_norm * sradius / ucell.omega * ucell.lat0 * ucell.lat0 * ucell.lat0;
    int extend_d22 = std::ceil(extend_d2);

    double a12_1 = ucell.latvec.e12 * ucell.latvec.e23 - ucell.latvec.e13 * ucell.latvec.e22;
    double a12_2 = ucell.latvec.e11 * ucell.latvec.e23 - ucell.latvec.e13 * ucell.latvec.e21;
    double a12_3 = ucell.latvec.e11 * ucell.latvec.e22 - ucell.latvec.e12 * ucell.latvec.e21;
    double a12_norm = sqrt(a12_1 * a12_1 + a12_2 * a12_2 + a12_3 * a12_3);
    double extend_d3 = a12_norm * sradius / ucell.omega * ucell.lat0 * ucell.lat0 * ucell.lat0;
    int extend_d33 = std::ceil(extend_d3);
    // 2016-09-05, LiuXh

    glayerX = extend_d11 + 1;
    glayerY = extend_d22 + 1;
    glayerZ = extend_d33 + 1;
    glayerX_minus = extend_d11;
    glayerY_minus = extend_d22;
    glayerZ_minus = extend_d33;
    // End, 2016-09-05, LiuXh

}


void Grid::setMemberVariables(std::ofstream& ofs_in, //  output data to ofs
                              const UnitCell& ucell)
{
    ModuleBase::TITLE("SLTK_Grid", "setMemberVariables");

    this->clear_atoms();

    // random selection, in order to estimate again.
    for (int it = 0; it < ucell.ntype; it++)
    {
        if (ucell.atoms[it].na > 0)
        {
            this->x_min = ucell.atoms[it].tau[0].x;
            this->y_min = ucell.atoms[it].tau[0].y;
            this->z_min = ucell.atoms[it].tau[0].z;
            this->x_max = ucell.atoms[it].tau[0].x;
            this->y_max = ucell.atoms[it].tau[0].y;
            this->z_max = ucell.atoms[it].tau[0].z;
            break;
        }
    }

    ModuleBase::Vector3<double> vec1(ucell.latvec.e11, ucell.latvec.e12, ucell.latvec.e13);
    ModuleBase::Vector3<double> vec2(ucell.latvec.e21, ucell.latvec.e22, ucell.latvec.e23);
    ModuleBase::Vector3<double> vec3(ucell.latvec.e31, ucell.latvec.e32, ucell.latvec.e33);

    // calculate min & max value
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
                        x_min = std::min(x_min, x);
                        x_max = std::max(x_max, x);
                        y_min = std::min(y_min, y);
                        y_max = std::max(y_max, y);
                        z_min = std::min(z_min, z);
                        z_max = std::max(z_max, z);
                    }
                }
            }
        }
    }
    ModuleBase::GlobalFunc::OUT(ofs_in, "Find the coordinate range of the input atom(unit:lat0).");
    ModuleBase::GlobalFunc::OUT(ofs_in, "min_tau", x_min, y_min, z_min);
    ModuleBase::GlobalFunc::OUT(ofs_in, "max_tau", x_max, y_max, z_max);

    this->box_edge_length = sradius + 0.1; // To avoid edge cases, the size of the box is slightly increased.

/*  warning box algorithm   
    this->box_nx = std::ceil((this->x_max - this->x_min) / box_edge_length) + 1;
    this->box_ny = std::ceil((this->y_max - this->y_min) / box_edge_length) + 1;
    this->box_nz = std::ceil((this->z_max - this->z_min) / box_edge_length) + 1;
    ModuleBase::GlobalFunc::OUT(ofs_in, "BoxNumber", box_nx, box_ny, box_nz);

    atoms_in_box.resize(this->box_nx);
    for (int i = 0; i < this->box_nx; i++)
    {
        atoms_in_box[i].resize(this->box_ny);
        for (int j = 0; j < this->box_ny; j++)
        {
            atoms_in_box[i][j].resize(this->box_nz);
        }
    }
 */
    this->box_nx = glayerX + glayerX_minus;
    this->box_ny = glayerY + glayerY_minus;
    this->box_nz = glayerZ + glayerZ_minus;
    ModuleBase::GlobalFunc::OUT(ofs_in, "BoxNumber", box_nx, box_ny, box_nz);

    atoms_in_box.resize(this->box_nx);
    for (int i = 0; i < this->box_nx; i++)
    {
        atoms_in_box[i].resize(this->box_ny);
        for (int j = 0; j < this->box_ny; j++)
        {
            atoms_in_box[i][j].resize(this->box_nz);
        }
    }
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
                        int box_i_x, box_i_y, box_i_z;
                        //this->getBox(box_i_x, box_i_y, box_i_z, x, y, z);
                        box_i_x = ix + glayerX_minus;
                        box_i_y = iy + glayerY_minus;
                        box_i_z = iz + glayerZ_minus;
                        this->atoms_in_box[box_i_x][box_i_y][box_i_z].push_back(atom);
                    }
                }
            }
        }
    }
    
    this->all_adj_info.resize(ucell.ntype);
    for (int i = 0; i < ucell.ntype; i++)
    {
        this->all_adj_info[i].resize(ucell.atoms[i].na);
    }
}

void Grid::Construct_Adjacent(const UnitCell& ucell)
{
    ModuleBase::timer::tick("Grid", "Construct_Adjacent_expand");

    for  (int i_type = 0; i_type < ucell.ntype; i_type++)
    {
        for (int j_atom = 0; j_atom < ucell.atoms[i_type].na; j_atom++)
        {

            FAtom atom(ucell.atoms[i_type].tau[j_atom].x,
                     ucell.atoms[i_type].tau[j_atom].y,
                     ucell.atoms[i_type].tau[j_atom].z,
                     i_type,
                     j_atom,
                     0, 0 ,0);

            this->Construct_Adjacent_near_box(atom);
        }
    }
    ModuleBase::timer::tick("Grid", "Construct_Adjacent_expand");
}

void Grid::Construct_Adjacent_near_box(const FAtom& fatom)
{
    //	if (test_grid)ModuleBase::TITLE(ofs_running, "Grid", "Construct_Adjacent_expand_periodic");
    ModuleBase::timer::tick("Grid", "Construct_Adjacent_expand_periodic");
    int box_i_x, box_i_y, box_i_z;
    this->getBox(box_i_x, box_i_y, box_i_z, fatom.x, fatom.y, fatom.z);

/*     for (int box_i_x_adj = std::max(box_i_x - 1, 0); box_i_x_adj <= std::min(box_i_x + 1, box_nx - 1); box_i_x_adj++)
    {
        for (int box_i_y_adj = std::max(box_i_y - 1, 0); box_i_y_adj <= std::min(box_i_y + 1, box_ny - 1); box_i_y_adj++)
        {
            for (int box_i_z_adj = std::max(box_i_z - 1, 0); box_i_z_adj <= std::min(box_i_z + 1, box_nz - 1); box_i_z_adj++)
            {
 */             
    for (int box_i_x_adj = 0; box_i_x_adj < glayerX + glayerX_minus; box_i_x_adj++)
    {
        for (int box_i_y_adj = 0; box_i_y_adj < glayerY + glayerY_minus; box_i_y_adj++)
        {
            for (int box_i_z_adj = 0; box_i_z_adj < glayerZ + glayerZ_minus; box_i_z_adj++)
            {
                for (auto &fatom2 : this->atoms_in_box[box_i_x_adj][box_i_y_adj][box_i_z_adj])
                {
                    this->Construct_Adjacent_final(fatom, &fatom2);
                }
            }
        }
    }
    ModuleBase::timer::tick("Grid", "Construct_Adjacent_expand_periodic");
}

void Grid::Construct_Adjacent_final(const FAtom& fatom1,
                                    FAtom* fatom2)
{
    double delta_x = fatom1.x - fatom2->x;
    double delta_y = fatom1.y - fatom2->y;
    double delta_z = fatom1.z - fatom2->z;

    double dr = delta_x * delta_x + delta_y * delta_y + delta_z * delta_z;


    // 20241204 zhanghaochong
    // dr == 0 means the same atom
    // the atom itself is neighbour atom, but the order itself must on last in the list.
    // so we will add itself on find atom function, and skip here.
    // I dont know why, but if we add self here, test 701_LJ_MD_Anderson will assert
    if (dr != 0.0 && dr <= this->sradius2)
    {
        all_adj_info[fatom1.type][fatom1.natom].push_back(fatom2);
    }
}
