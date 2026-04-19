#include "source_cell/module_neighbor_search/neighbor_search.h"
#include <cmath>
#include <algorithm>
#include <limits>


InputAtoms NeighborSearch::ucell_to_input_atoms(const UnitCell& ucell)
{
    InputAtoms input_atoms;
    int atom_count = 0;

    input_atoms.x_low = ucell.atoms[0].tau[0].x;
    input_atoms.x_high = ucell.atoms[0].tau[0].x;
    input_atoms.y_low = ucell.atoms[0].tau[0].y;
    input_atoms.y_high = ucell.atoms[0].tau[0].y;
    input_atoms.z_low = ucell.atoms[0].tau[0].z;
    input_atoms.z_high = ucell.atoms[0].tau[0].z;

    for (int i = 0; i < ucell.ntype; i++)
    {
        for (int j = 0; j < ucell.atoms[i].na; j++)
        {
            NeighborAtom atom(
                ucell.atoms[i].tau[j].x,
                ucell.atoms[i].tau[j].y,
                ucell.atoms[i].tau[j].z,
                i,
                j,
                atom_count
            );
            input_atoms.InputAtom.push_back(atom);

            input_atoms.x_low = std::min(input_atoms.x_low, atom.position_x);
            input_atoms.x_high = std::max(input_atoms.x_high, atom.position_x);
            input_atoms.y_low = std::min(input_atoms.y_low, atom.position_y);
            input_atoms.y_high = std::max(input_atoms.y_high, atom.position_y);
            input_atoms.z_low = std::min(input_atoms.z_low, atom.position_z);
            input_atoms.z_high = std::max(input_atoms.z_high, atom.position_z);

            atom_count++;
        }
    }

    input_atoms.n_atoms = atom_count;
    return input_atoms;
}

void NeighborSearch::init(const UnitCell& ucell, double sr, int mpi_rank)
{
    search_radius = sr / ucell.lat0;
    Check_Expand_Condition(ucell);
    setMemberVariables(ucell);
    InputAtoms atoms = ucell_to_input_atoms(ucell);

    int mpi_size = 1;
    int nx, ny, nz;
    decompose(mpi_size, nx, ny, nz);

    z = mpi_rank / (nx * ny);
    y = (mpi_rank % (nx * ny)) / ny;
    x = mpi_rank % (nx * ny) % ny;

    wide_x = (atoms.x_high - atoms.x_low) / nx;
    wide_y = (atoms.y_high - atoms.y_low) / ny;
    wide_z = (atoms.z_high - atoms.z_low) / nz;

    for (int i = 0; i < all_atoms.size(); i++)
    {
        bool in_x = (all_atoms[i].position_x >= atoms.x_low + x * wide_x) &&
                    (all_atoms[i].position_x <= atoms.x_low + (x + 1) * wide_x);
        bool in_y = (all_atoms[i].position_y >= atoms.y_low + y * wide_y) &&
                    (all_atoms[i].position_y <= atoms.y_low + (y + 1) * wide_y);
        bool in_z = (all_atoms[i].position_z >= atoms.z_low + z * wide_z) &&
                    (all_atoms[i].position_z <= atoms.z_low + (z + 1) * wide_z);

        if (in_x && in_y && in_z)
        {
            all_atoms[i].isghost = false;
            inside_atoms.push_back(all_atoms[i]);
        }
        else if (distance(
            all_atoms[i].position_x,
            all_atoms[i].position_y,
            all_atoms[i].position_z,
            atoms.x_low,
            atoms.y_low,
            atoms.z_low) <= search_radius * search_radius)
        {
            all_atoms[i].isghost = true;
            ghost_atoms.push_back(all_atoms[i]);
        }
    }

    neighbor_list.initialize(ucell.nat, 10000000);
}

void NeighborSearch::build_neighbors()
{
    bin_manager.init_bins(search_radius, inside_atoms, ghost_atoms);
    bin_manager.do_binning(inside_atoms, ghost_atoms);
    bin_manager.build_atom_neighbors(neighbor_list, inside_atoms);
}

void NeighborSearch::Check_Expand_Condition(const UnitCell& ucell)
{
    double a23_1 = ucell.latvec.e22 * ucell.latvec.e33 - ucell.latvec.e23 * ucell.latvec.e32;
    double a23_2 = ucell.latvec.e21 * ucell.latvec.e33 - ucell.latvec.e23 * ucell.latvec.e31;
    double a23_3 = ucell.latvec.e21 * ucell.latvec.e32 - ucell.latvec.e22 * ucell.latvec.e31;
    double a23_norm = sqrt(a23_1 * a23_1 + a23_2 * a23_2 + a23_3 * a23_3);
    double extend_v = a23_norm * search_radius;
    double extend_d1 = extend_v / ucell.omega * ucell.lat0 * ucell.lat0 * ucell.lat0;
    int extend_d11 = std::ceil(extend_d1);

    double a31_1 = ucell.latvec.e32 * ucell.latvec.e13 - ucell.latvec.e33 * ucell.latvec.e12;
    double a31_2 = ucell.latvec.e31 * ucell.latvec.e13 - ucell.latvec.e33 * ucell.latvec.e11;
    double a31_3 = ucell.latvec.e31 * ucell.latvec.e12 - ucell.latvec.e32 * ucell.latvec.e11;
    double a31_norm = sqrt(a31_1 * a31_1 + a31_2 * a31_2 + a31_3 * a31_3);
    double extend_d2 = a31_norm * search_radius / ucell.omega * ucell.lat0 * ucell.lat0 * ucell.lat0;
    int extend_d22 = std::ceil(extend_d2);

    double a12_1 = ucell.latvec.e12 * ucell.latvec.e23 - ucell.latvec.e13 * ucell.latvec.e22;
    double a12_2 = ucell.latvec.e11 * ucell.latvec.e23 - ucell.latvec.e13 * ucell.latvec.e21;
    double a12_3 = ucell.latvec.e11 * ucell.latvec.e22 - ucell.latvec.e12 * ucell.latvec.e21;
    double a12_norm = sqrt(a12_1 * a12_1 + a12_2 * a12_2 + a12_3 * a12_3);
    double extend_d3 = a12_norm * search_radius / ucell.omega * ucell.lat0 * ucell.lat0 * ucell.lat0;
    int extend_d33 = std::ceil(extend_d3);

    glayerX = extend_d11 + 1;
    glayerY = extend_d22 + 1;
    glayerZ = extend_d33 + 1;
    glayerX_minus = extend_d11;
    glayerY_minus = extend_d22;
    glayerZ_minus = extend_d33;
}

void NeighborSearch::setMemberVariables(const UnitCell& ucell)
{
    all_atoms.clear();

    ModuleBase::Vector3<double> vec1(ucell.latvec.e11, ucell.latvec.e12, ucell.latvec.e13);
    ModuleBase::Vector3<double> vec2(ucell.latvec.e21, ucell.latvec.e22, ucell.latvec.e23);
    ModuleBase::Vector3<double> vec3(ucell.latvec.e31, ucell.latvec.e32, ucell.latvec.e33);

    int atom_count = 0;

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

                        NeighborAtom atom(x, y, z, i, j, atom_count);
                        all_atoms.push_back(atom);
                        atom_count++;
                    }
                }
            }
        }
    }
}

double NeighborSearch::distance(
    double position_x,
    double position_y,
    double position_z,
    double x_low,
    double y_low,
    double z_low)
{
    double dx = std::max(0.0, std::max(x_low + x * wide_x - position_x, position_x - (x_low + (x + 1) * wide_x)));
    double dy = std::max(0.0, std::max(y_low + y * wide_y - position_y, position_y - (y_low + (y + 1) * wide_y)));
    double dz = std::max(0.0, std::max(z_low + z * wide_z - position_z, position_z - (z_low + (z + 1) * wide_z)));
    return dx * dx + dy * dy + dz * dz;
}

void NeighborSearch::decompose(int mpi_size, int &nx, int &ny, int &nz)
{
    nx = 1;
    ny = 1;
    nz = mpi_size;

    int cube = cbrt(mpi_size);
    for (int i = cube; i >= 1; i--)
    {
        if (mpi_size % i == 0)
        {
            nx = i;
            ny = mpi_size / i;
            break;
        }
    }

    int sq = sqrt(ny);
    for (int i = sq; i >= 1; i--)
    {
        if (ny % i == 0)
        {
            nz = ny / i;
            ny = i;
            break;
        }
    }
}