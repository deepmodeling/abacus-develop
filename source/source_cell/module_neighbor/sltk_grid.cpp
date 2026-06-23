#include "sltk_grid.h"

#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_base/memory.h"
#include "source_base/timer.h"

#include <algorithm>
#include <cmath>

namespace
{
int clamp_box_index(const int index, const int size)
{
    if (index < 0)
    {
        return 0;
    }
    if (index >= size)
    {
        return size - 1;
    }
    return index;
}
}

Grid::Grid(const int& test_grid_in) : test_grid(test_grid_in)
{
}

Grid::~Grid()
{
    this->clear_atoms();
}

void Grid::init(std::ofstream& ofs_in,
                const UnitCell& ucell,
                const double radius_in,
                const bool boundary,
                const NeighborBuildMode build_mode,
                const NeighborSearchMode search_mode,
                const NeighborReferenceMode reference_mode)
{
    ModuleBase::TITLE("Grid", "init");
    ModuleBase::timer::start("Grid", "init");
    this->pbc = boundary;
    this->sradius2 = radius_in * radius_in;
    this->sradius = radius_in;

//    ModuleBase::GlobalFunc::OUT(ofs_in, "PeriodicBoundary", this->pbc);
    ModuleBase::GlobalFunc::OUT(ofs_in, "Radius (unit: lattice constant)", sradius);

    this->Check_Expand_Condition(ucell);
    ModuleBase::GlobalFunc::OUT(ofs_in, "Max number of cells", glayerX, glayerY, glayerZ);
    ModuleBase::GlobalFunc::OUT(ofs_in, "Min number of cells", glayerX_minus, glayerY_minus, glayerZ_minus);

    this->setMemberVariables(ofs_in, ucell);
    this->Build_AtomPack_Search_Path(ucell, search_mode, reference_mode);
    if (build_mode == NeighborBuildMode::AtomPackAndLegacy)
    {
        this->Build_Legacy_Search_Path(ucell);
    }
    ModuleBase::timer::end("Grid", "init");
}

void Grid::Check_Expand_Condition(const UnitCell& ucell)
{
    //	ModuleBase::TITLE(GlobalV::ofs_running, "Atom_input", "Check_Expand_Condition");

    if (!pbc)
    {
        glayerX = 1;
        glayerY = 1;
        glayerZ = 1;
        glayerX_minus = 0;
        glayerY_minus = 0;
        glayerZ_minus = 0;
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

    const int ix_begin = this->pbc ? -glayerX_minus : 0;
    const int ix_end = this->pbc ? glayerX : 1;
    const int iy_begin = this->pbc ? -glayerY_minus : 0;
    const int iy_end = this->pbc ? glayerY : 1;
    const int iz_begin = this->pbc ? -glayerZ_minus : 0;
    const int iz_end = this->pbc ? glayerZ : 1;

    // calculate min & max value
    for (int ix = ix_begin; ix < ix_end; ix++)
    {
        for (int iy = iy_begin; iy < iy_end; iy++)
        {
            for (int iz = iz_begin; iz < iz_end; iz++)
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

//    ofs_in << " RANGE OF ATOMIC COORDINATES (unit: lat0)" << std::endl;
    ModuleBase::GlobalFunc::OUT(ofs_in, "Min coordinates of atoms", x_min, y_min, z_min);
    ModuleBase::GlobalFunc::OUT(ofs_in, "Max coordinates of atoms", x_max, y_max, z_max);

    this->box_edge_length = sradius + 0.1; // To avoid edge cases, the size of the box is slightly increased.

    this->box_nx = std::max(1, static_cast<int>(std::floor((this->x_max - this->x_min) / box_edge_length)) + 1);
    this->box_ny = std::max(1, static_cast<int>(std::floor((this->y_max - this->y_min) / box_edge_length)) + 1);
    this->box_nz = std::max(1, static_cast<int>(std::floor((this->z_max - this->z_min) / box_edge_length)) + 1);
    ModuleBase::GlobalFunc::OUT(ofs_in, "Number of needed cells", box_nx, box_ny, box_nz);

}

void Grid::Build_Legacy_Search_Path(const UnitCell& ucell)
{
    atoms_in_box.clear();
    all_adj_info.clear();

    ModuleBase::Vector3<double> vec1(ucell.latvec.e11, ucell.latvec.e12, ucell.latvec.e13);
    ModuleBase::Vector3<double> vec2(ucell.latvec.e21, ucell.latvec.e22, ucell.latvec.e23);
    ModuleBase::Vector3<double> vec3(ucell.latvec.e31, ucell.latvec.e32, ucell.latvec.e33);

    const int ix_begin = this->pbc ? -glayerX_minus : 0;
    const int ix_end = this->pbc ? glayerX : 1;
    const int iy_begin = this->pbc ? -glayerY_minus : 0;
    const int iy_end = this->pbc ? glayerY : 1;
    const int iz_begin = this->pbc ? -glayerZ_minus : 0;
    const int iz_end = this->pbc ? glayerZ : 1;

    atoms_in_box.resize(this->box_nx);
    for (int i = 0; i < this->box_nx; i++)
    {
        atoms_in_box[i].resize(this->box_ny);
        for (int j = 0; j < this->box_ny; j++)
        {
            atoms_in_box[i][j].resize(this->box_nz);
        }
    }
    for (int ix = ix_begin; ix < ix_end; ix++)
    {
        for (int iy = iy_begin; iy < iy_end; iy++)
        {
            for (int iz = iz_begin; iz < iz_end; iz++)
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
                        this->getBox(box_i_x, box_i_y, box_i_z, x, y, z);
                        box_i_x = clamp_box_index(box_i_x, this->box_nx);
                        box_i_y = clamp_box_index(box_i_y, this->box_ny);
                        box_i_z = clamp_box_index(box_i_z, this->box_nz);
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

    this->Construct_Adjacent(ucell);
}

void Grid::Construct_Adjacent(const UnitCell& ucell)
{
    ModuleBase::timer::start("Grid", "constru_adj");

    for (int i_type = 0; i_type < ucell.ntype; i_type++)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
        for (int j_atom = 0; j_atom < ucell.atoms[i_type].na; j_atom++)
        {
            FAtom atom(ucell.atoms[i_type].tau[j_atom].x,
                       ucell.atoms[i_type].tau[j_atom].y,
                       ucell.atoms[i_type].tau[j_atom].z,
                       i_type,
                       j_atom,
                       0, 0, 0);

            this->Construct_Adjacent_near_box(atom);
        }
    }
    ModuleBase::timer::end("Grid", "constru_adj");
}

void Grid::Construct_Adjacent_near_box(const FAtom& fatom)
{
    int box_i_x=0;
    int box_i_y=0;
    int box_i_z=0;
    this->getBox(box_i_x, box_i_y, box_i_z, fatom.x, fatom.y, fatom.z);
    box_i_x = clamp_box_index(box_i_x, this->box_nx);
    box_i_y = clamp_box_index(box_i_y, this->box_ny);
    box_i_z = clamp_box_index(box_i_z, this->box_nz);

    const int search_layer = std::max(1, static_cast<int>(std::ceil(this->sradius / this->box_edge_length)));
    const int box_x_begin = std::max(0, box_i_x - search_layer);
    const int box_x_end = std::min(this->box_nx - 1, box_i_x + search_layer);
    const int box_y_begin = std::max(0, box_i_y - search_layer);
    const int box_y_end = std::min(this->box_ny - 1, box_i_y + search_layer);
    const int box_z_begin = std::max(0, box_i_z - search_layer);
    const int box_z_end = std::min(this->box_nz - 1, box_i_z + search_layer);

    for (int box_i_x_adj = box_x_begin; box_i_x_adj <= box_x_end; box_i_x_adj++)
    {
        for (int box_i_y_adj = box_y_begin; box_i_y_adj <= box_y_end; box_i_y_adj++)
        {
            for (int box_i_z_adj = box_z_begin; box_i_z_adj <= box_z_end; box_i_z_adj++)
            {
                for (auto &fatom2 : this->atoms_in_box[box_i_x_adj][box_i_y_adj][box_i_z_adj])
                {
                    this->Construct_Adjacent_final(fatom, &fatom2);
                }
            }
        }
    }
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

void Grid::Build_AtomPack_Search_Path(const UnitCell& ucell,
                                      const NeighborSearchMode search_mode,
                                      const NeighborReferenceMode reference_mode)
{
    // Build the Phase 2.1 integer-indexed path. The production query list is
    // neighbor_pairs; neighbor_pairs_27 is populated only when tests explicitly
    // request a full-search reference.
    atom_pack = ModuleNeighbor::build_atom_pack_from_unitcell(ucell, sradius, pbc);
    grid_storage = ModuleNeighbor::build_grid_storage_from_atom_pack(atom_pack, box_edge_length);
    if (search_mode == NeighborSearchMode::Full27)
    {
        neighbor_pairs = ModuleNeighbor::build_neighbor_pairs_27(atom_pack, grid_storage, sradius);
    }
    else
    {
        neighbor_pairs = ModuleNeighbor::build_neighbor_pairs_14(atom_pack, grid_storage, sradius);
    }

    neighbor_pairs_27.clear();
    if (reference_mode == NeighborReferenceMode::Full27)
    {
        neighbor_pairs_27 = ModuleNeighbor::build_neighbor_pairs_27(atom_pack, grid_storage, sradius);
    }

    // Convert the global pair list to per-center lookup indices so Find_atom()
    // can keep its existing single-atom query interface without scanning all
    // neighbor pairs each time.
    neighbor_pair_indices.clear();
    neighbor_pair_indices.resize(ucell.ntype);
    for (int it = 0; it < ucell.ntype; ++it)
    {
        neighbor_pair_indices[it].resize(ucell.atoms[it].na);
    }

    for (int pair_index = 0; pair_index < static_cast<int>(neighbor_pairs.size()); ++pair_index)
    {
        const ModuleNeighbor::NeighborPair& pair = neighbor_pairs[pair_index];
        neighbor_pair_indices[pair.center_type][pair.center_natom].push_back(pair_index);
    }
}
