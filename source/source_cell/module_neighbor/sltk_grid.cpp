#include "sltk_grid.h"

#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_base/timer.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>

namespace
{
constexpr double minimum_cartesian_bin_edge = 0.1;
constexpr std::size_t max_cartesian_bin_count = 1000000;

int count_bins_for_range(const double coordinate_range, const double box_edge_length)
{
    if (coordinate_range <= 0.0 || box_edge_length <= 0.0)
    {
        return 1;
    }

    const double bin_count = std::floor(coordinate_range / box_edge_length) + 1.0;
    if (!std::isfinite(bin_count)
        || bin_count > static_cast<double>(std::numeric_limits<int>::max()))
    {
        return std::numeric_limits<int>::max();
    }
    return std::max(1, static_cast<int>(bin_count));
}

std::size_t saturated_bin_count_product(const int nx, const int ny, const int nz)
{
    if (nx <= 0 || ny <= 0 || nz <= 0)
    {
        return 0;
    }

    const std::size_t sx = static_cast<std::size_t>(nx);
    const std::size_t sy = static_cast<std::size_t>(ny);
    const std::size_t sz = static_cast<std::size_t>(nz);
    const std::size_t max_size = std::numeric_limits<std::size_t>::max();
    if (sx > max_size / sy)
    {
        return max_size;
    }
    const std::size_t xy = sx * sy;
    if (xy > max_size / sz)
    {
        return max_size;
    }
    return xy * sz;
}

void set_box_counts_for_edge(const double range_x,
                             const double range_y,
                             const double range_z,
                             const double box_edge_length,
                             int& box_nx,
                             int& box_ny,
                             int& box_nz)
{
    box_nx = count_bins_for_range(range_x, box_edge_length);
    box_ny = count_bins_for_range(range_y, box_edge_length);
    box_nz = count_bins_for_range(range_z, box_edge_length);
}

double estimate_bin_count_for_edge(const double range_x,
                                   const double range_y,
                                   const double range_z,
                                   const double box_edge_length)
{
    const double nx = std::floor(range_x / box_edge_length) + 1.0;
    const double ny = std::floor(range_y / box_edge_length) + 1.0;
    const double nz = std::floor(range_z / box_edge_length) + 1.0;
    return nx * ny * nz;
}

double coarsen_edge_for_bin_cap(const double initial_edge,
                                const double range_x,
                                const double range_y,
                                const double range_z)
{
    double edge = initial_edge;
    const double max_range = std::max(range_x, std::max(range_y, range_z));

    for (int iteration = 0; iteration < 16; ++iteration)
    {
        const double estimated_count = estimate_bin_count_for_edge(range_x, range_y, range_z, edge);
        if (std::isfinite(estimated_count)
            && estimated_count <= static_cast<double>(max_cartesian_bin_count))
        {
            return edge;
        }

        if (!std::isfinite(estimated_count))
        {
            return std::max(edge, max_range);
        }

        const double scale = std::max(2.0,
                                      std::cbrt(estimated_count
                                                / static_cast<double>(max_cartesian_bin_count)));
        edge *= scale;
    }

    return std::max(edge, max_range);
}

double squared_distance(const FAtom& lhs, const FAtom& rhs)
{
    const double delta_x = lhs.x - rhs.x;
    const double delta_y = lhs.y - rhs.y;
    const double delta_z = lhs.z - rhs.z;
    return delta_x * delta_x + delta_y * delta_y + delta_z * delta_z;
}

bool is_neighbor_by_current_rule(const double distance2, const double radius2)
{
    const bool is_not_zero_distance = distance2 != 0.0;
    const bool is_within_search_radius = distance2 <= radius2;
    return is_not_zero_distance && is_within_search_radius;
}
} // namespace

Grid::Grid(const int& test_grid_in) : test_grid(test_grid_in)
{
}

Grid::~Grid()
{
    this->clear_atoms();
}

void Grid::init(std::ofstream& ofs_in, const UnitCell& ucell, const double radius_in, const bool boundary)
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
    this->Construct_Adjacent(ucell);
    ModuleBase::timer::end("Grid", "init");
}

void Grid::Check_Expand_Condition(const UnitCell& ucell)
{
    //	ModuleBase::TITLE(GlobalV::ofs_running, "Atom_input", "Check_Expand_Condition");

    glayerX = 0;
    glayerX_minus = 0;
    glayerY = 0;
    glayerY_minus = 0;
    glayerZ = 0;
    glayerZ_minus = 0;

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

//    ofs_in << " RANGE OF ATOMIC COORDINATES (unit: lat0)" << std::endl;
    ModuleBase::GlobalFunc::OUT(ofs_in, "Min coordinates of atoms", x_min, y_min, z_min);
    ModuleBase::GlobalFunc::OUT(ofs_in, "Max coordinates of atoms", x_max, y_max, z_max);

    // Keep each Cartesian bin at least as wide as the search radius so a 27-bin scan is sufficient.
    // Coarsen very sparse AABB grids to avoid allocating millions of empty bins.
    this->box_edge_length = std::max(sradius, minimum_cartesian_bin_edge);

    const bool has_expanded_images = (glayerX + glayerX_minus) > 0 && (glayerY + glayerY_minus) > 0
                                     && (glayerZ + glayerZ_minus) > 0;
    if (has_expanded_images)
    {
        const double range_x = std::max(0.0, this->x_max - this->x_min);
        const double range_y = std::max(0.0, this->y_max - this->y_min);
        const double range_z = std::max(0.0, this->z_max - this->z_min);

        this->box_edge_length = coarsen_edge_for_bin_cap(this->box_edge_length, range_x, range_y, range_z);
        set_box_counts_for_edge(range_x, range_y, range_z, this->box_edge_length, box_nx, box_ny, box_nz);
        while (saturated_bin_count_product(box_nx, box_ny, box_nz) > max_cartesian_bin_count)
        {
            this->box_edge_length *= 2.0;
            set_box_counts_for_edge(range_x, range_y, range_z, this->box_edge_length, box_nx, box_ny, box_nz);
        }
    }
    else
    {
        this->box_nx = 0;
        this->box_ny = 0;
        this->box_nz = 0;
    }
    ModuleBase::GlobalFunc::OUT(ofs_in, "Number of needed cells", box_nx, box_ny, box_nz);

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
                        int box_i_x = 0;
                        int box_i_y = 0;
                        int box_i_z = 0;
                        this->getBox(box_i_x, box_i_y, box_i_z, x, y, z);
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
    ModuleBase::timer::start("Grid", "constru_adj");
    for (int i_type = 0; i_type < ucell.ntype; i_type++)
    {
        for (int j_atom = 0; j_atom < ucell.atoms[i_type].na; j_atom++)
        {
            FAtom atom(ucell.atoms[i_type].tau[j_atom].x,
                       ucell.atoms[i_type].tau[j_atom].y,
                       ucell.atoms[i_type].tau[j_atom].z,
                       i_type,
                       j_atom,
                       0,
                       0,
                       0);

            std::vector<FAtom*>& adjacent_atoms = this->all_adj_info[i_type][j_atom];
            this->Construct_Adjacent_near_box(atom, adjacent_atoms);
        }
    }
    ModuleBase::timer::end("Grid", "constru_adj");
}

void Grid::Construct_Adjacent_near_box(const FAtom& fatom, std::vector<FAtom*>& adjacent_atoms)
{
    if (box_nx <= 0 || box_ny <= 0 || box_nz <= 0)
    {
        return;
    }

    int box_i_x=0;
    int box_i_y=0;
    int box_i_z=0;
    this->getBox(box_i_x, box_i_y, box_i_z, fatom.x, fatom.y, fatom.z);

    const int box_x_begin = std::max(0, box_i_x - 1);
    const int box_x_end = std::min(box_nx - 1, box_i_x + 1);
    const int box_y_begin = std::max(0, box_i_y - 1);
    const int box_y_end = std::min(box_ny - 1, box_i_y + 1);
    const int box_z_begin = std::max(0, box_i_z - 1);
    const int box_z_end = std::min(box_nz - 1, box_i_z + 1);

    for (int box_i_x_adj = box_x_begin; box_i_x_adj <= box_x_end; box_i_x_adj++)
    {
        for (int box_i_y_adj = box_y_begin; box_i_y_adj <= box_y_end; box_i_y_adj++)
        {
            for (int box_i_z_adj = box_z_begin; box_i_z_adj <= box_z_end; box_i_z_adj++)
            {
                for (auto& fatom2 : this->atoms_in_box[box_i_x_adj][box_i_y_adj][box_i_z_adj])
                {
                    this->Construct_Adjacent_final(fatom, &fatom2, adjacent_atoms);
                }
            }
        }
    }
}

void Grid::Construct_Adjacent_final(const FAtom& fatom1,
                                    FAtom* fatom2,
                                    std::vector<FAtom*>& adjacent_atoms) const
{
    const double distance2 = squared_distance(fatom1, *fatom2);

    // 20241204 zhanghaochong
    // dr == 0 means the same atom
    // the atom itself is neighbour atom, but the order itself must on last in the list.
    // so we will add itself on find atom function, and skip here.
    // I dont know why, but if we add self here, test 701_LJ_MD_Anderson will assert
    if (is_neighbor_by_current_rule(distance2, this->sradius2))
    {
        adjacent_atoms.push_back(fatom2);
    }
}
