#include "sltk_atom_input.h"

#include "module_base/memory.h"
#include "module_parameter/parameter.h"
#include "sltk_grid.h"

//==========================================================
// define constructor and deconstructor
//==========================================================
Atom_input::Atom_input(std::ofstream& ofs_in,
                       const UnitCell& ucell,
                       const int amount,
                       const int ntype,
                       const bool boundary_in,
                       const double radius_in,
                       const int& test_atom_in)
    : periodic_boundary(boundary_in), radius(radius_in), expand_flag(false), glayerX(1), glayerX_minus(0), glayerY(1),
      glayerY_minus(0), glayerZ(1), glayerZ_minus(0), test_atom_input(test_atom_in)
{
    ModuleBase::TITLE("Atom_input", "Atom_input");

    if (test_atom_input)
    {
        ModuleBase::GlobalFunc::OUT(ofs_in, "ntype", ntype);
        ModuleBase::GlobalFunc::OUT(ofs_in, "Amount(atom number)", amount);
        ModuleBase::GlobalFunc::OUT(ofs_in, "Periodic_boundary", periodic_boundary);
        ModuleBase::GlobalFunc::OUT(ofs_in, "Searching radius(lat0)", radius);
    }

    if (radius < 0.0)
    {
        ModuleBase::WARNING_QUIT("atom_arrange::init", " search radius < 0,forbidden");
    }
    // random selection, in order to estimate again.
    this->x_min = ucell.atoms[0].tau[0].x;
    this->y_min = ucell.atoms[0].tau[0].y;
    this->z_min = ucell.atoms[0].tau[0].z;
    this->x_max = ucell.atoms[0].tau[0].x;
    this->y_max = ucell.atoms[0].tau[0].y;
    this->z_max = ucell.atoms[0].tau[0].z;

    // calculate min & max value
    for (int i = 0; i < ntype; i++)
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

    if (test_atom_input)
    {
        ModuleBase::GlobalFunc::OUT(ofs_in, "Find the coordinate range of the input atom(unit:lat0).");
        ModuleBase::GlobalFunc::OUT(ofs_in, "min_tau", x_min, y_min, z_min);
        ModuleBase::GlobalFunc::OUT(ofs_in, "max_tau", x_max, y_max, z_max);
    }

    //----------------------------------------------------------
    // CALL MEMBER FUNCTION :
    // NAME : Check_Expand_Condition(check if swe need to
    // expand grid,and generate 6 MEMBER VARIABLE(number of
    // layers for 6 dimension)
    // initial value for "glayerX,Y,Z" : 1
    // (if > 2 ,expand flag = 1)
    // initial value for "glayerX,Y,Z_minus" : 0
    // ( if > 1 ,expand flag = 1)
    //----------------------------------------------------------

    this->Check_Expand_Condition(ucell);

    if (test_atom_input)
    {
        ModuleBase::GlobalFunc::OUT(ofs_in, "glayer+", glayerX, glayerY, glayerZ);
        ModuleBase::GlobalFunc::OUT(ofs_in, "glayer-", glayerX_minus, glayerY_minus, glayerZ_minus);
        ModuleBase::GlobalFunc::OUT(ofs_in, "expand_flag", expand_flag);
    }

    //----------------------------------------------------------
    // CALL MEMBER FUNCTION :
    // NAME : calculate_cells
    // Calculate how many cells we need in each direction.
    //----------------------------------------------------------
    this->calculate_cells();
    if (test_atom_input)
    {
        ModuleBase::GlobalFunc::OUT(ofs_in, "CellDim", cell_nx, cell_ny, cell_nz);
    }
    return;
}

Atom_input::~Atom_input()
{
}

//============================================
// !!!! May still have bug, be very careful!!
// should use the same algorithm to generate
// dxe, dye, dze in grid_meshcell.cpp.
//============================================
void Atom_input::Check_Expand_Condition(const UnitCell& ucell)
{
    //	ModuleBase::TITLE(GlobalV::ofs_running, "Atom_input", "Check_Expand_Condition");

    if (!periodic_boundary)
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
    double extend_v = a23_norm * radius;
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
    double extend_d2 = a31_norm * radius / ucell.omega * ucell.lat0 * ucell.lat0 * ucell.lat0;
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
    double extend_d3 = a12_norm * radius / ucell.omega * ucell.lat0 * ucell.lat0 * ucell.lat0;
    int extend_d33 = static_cast<int>(extend_d3);
    // 2016-09-05, LiuXh
    if (extend_d3 - extend_d33 > 0.0)
    {
        extend_d33 += 1;
    }

    glayerX = extend_d11 + 1;
    glayerY = extend_d22 + 1;
    glayerZ = extend_d33 + 1;
    // Begin, 2016-09-05, LiuXh
    // glayerX_minus = extend_d11 +1;
    // glayerY_minus = extend_d22 +1;
    // glayerZ_minus = extend_d33 +1;
    glayerX_minus = extend_d11;
    glayerY_minus = extend_d22;
    glayerZ_minus = extend_d33;
    // End, 2016-09-05, LiuXh

    if (glayerX == 1)
    {
        glayerX++;
    }
    if (glayerY == 1)
    {
        glayerY++;
    }
    if (glayerZ == 1)
    {
        glayerZ++;
    }
    if (glayerX_minus == 1)
    {
        glayerX_minus++;
    }
    if (glayerY_minus == 1)
    {
        glayerY_minus++;
    }
    if (glayerZ_minus == 1)
    {
        glayerZ_minus++;
    }
    // End, 2016-07-19, LiuXh
    /*
        if(test_atom_input)
        {
            GlobalV::ofs_running << " Extend distance from the (minX,minY,minZ) direct position in this unitcell: " <<
       std::endl;
        }

        if(test_atom_input)OUT(GlobalV::ofs_running,"ExtentDim-",extent_1DX_minus,extent_1DY_minus,extent_1DZ_minus);
    */
    //----------------------------------------------------------
    // EXPLAIN : if extent don't satisfty the searching
    // requiment, we must expand one more layer
    //----------------------------------------------------------

    if (glayerX > 2 || glayerY > 2 || glayerZ > 2)
    {
        this->expand_flag = true;
    }
    else if (glayerX_minus > 1 || glayerX_minus > 1 || glayerX_minus > 1)
    {
        this->expand_flag = true;
    }
    else
    {
        this->expand_flag = false;
    }
    return;
}

void Atom_input::calculate_cells()
{
    ModuleBase::TITLE("Atom_input", "calculate_cells");
    //----------------------------------------------------------
    // EXPLAIN :
    // Expand_Case : Simple , we already know the cell numbres,
    // all the trouble is only to construct adjacentset using all
    // the cells.
    // Not_Expand_Case : Using searching radius to construct
    // the cells ,  trouble here,but is the convenience of searching
    // time , we then only need to search 27-adjacent cell for each cell.
    //----------------------------------------------------------
    if (expand_flag)
    {
        cell_nx = glayerX + glayerX_minus;
        cell_ny = glayerY + glayerY_minus;
        cell_nz = glayerZ + glayerZ_minus;
    }
    else
    {
        // maybe a bug, if we don't use direct
        // coordinates, mohan note 2011-04-14
        double real_nx, real_ny, real_nz;
        real_nx = (x_max - x_min) / radius;
        real_ny = (y_max - y_min) / radius;
        real_nz = (z_max - z_min) / radius;
        cell_nx = static_cast<int>(real_nx) + 1;
        cell_ny = static_cast<int>(real_ny) + 1;
        cell_nz = static_cast<int>(real_nz) + 1;
    }

    //================
    // Wrong !
    //================
    //	if(int_nx != real_nx) this->cell_nx++;
    //	if(int_ny != real_ny) this->cell_ny++;
    //	if(int_nz != real_nz) this->cell_nz++;
    //=======================================
    // Not need because if int_nx = real_nx,
    // the position belong to the next cell
    //=======================================
    return;
}
