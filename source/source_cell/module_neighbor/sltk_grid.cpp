#include "sltk_grid.h"

#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_base/timer.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <tuple>

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

int box_count_from_range(const double minimum,
                         const double maximum,
                         const double edge_length)
{
    const double count = std::floor((maximum - minimum) / edge_length) + 1.0;
    if (!std::isfinite(count) || count < 1.0
        || count > static_cast<double>(std::numeric_limits<int>::max()))
    {
        throw std::overflow_error("Grid box count exceeds the supported range.");
    }
    return static_cast<int>(count);
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
    if (!std::isfinite(radius_in) || radius_in < 0.0)
    {
        throw std::invalid_argument("Grid search radius must be finite and non-negative.");
    }
    ModuleBase::timer::start("Grid", "init");
    this->pbc = boundary;
    this->sradius2 = radius_in * radius_in;
    this->sradius = radius_in;
    this->query_radius2 = this->sradius2;
    this->query_radius = this->sradius;

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

void Grid::init_mpi(std::ofstream& ofs_in,
                    const UnitCell& ucell,
                    const double radius_in,
                    const bool boundary,
                    ModuleNeighbor::NeighborMpiComm communicator,
                    ModuleNeighbor::MpiGhostExchangeStats* stats)
{
    ModuleBase::TITLE("Grid", "init_mpi");
    if (!std::isfinite(radius_in) || radius_in < 0.0)
    {
        throw std::invalid_argument("Grid MPI search radius must be finite and non-negative.");
    }
    ModuleBase::timer::start("Grid", "init_mpi");

    this->clear_atoms();
    this->pbc = boundary;
    this->sradius = radius_in;
    this->sradius2 = radius_in * radius_in;
    this->query_radius = radius_in;
    this->query_radius2 = radius_in * radius_in;
    this->box_edge_length = radius_in + 0.1;

    const std::array<double, 3> origin{{0.0, 0.0, 0.0}};
    const std::array<std::array<double, 3>, 3> lattice_vectors{{
        std::array<double, 3>{{ucell.latvec.e11, ucell.latvec.e12, ucell.latvec.e13}},
        std::array<double, 3>{{ucell.latvec.e21, ucell.latvec.e22, ucell.latvec.e23}},
        std::array<double, 3>{{ucell.latvec.e31, ucell.latvec.e32, ucell.latvec.e33}},
    }};
    mpi_domain_.initialize_lattice(communicator, origin, lattice_vectors, radius_in, boundary);

    std::vector<ModuleNeighbor::MpiAtomRecord> all_records;
    all_records.reserve(ucell.nat);
    local_center_mask.resize(ucell.ntype);
    int global_index = 0;
    for (int type = 0; type < ucell.ntype; ++type)
    {
        local_center_mask[type].assign(ucell.atoms[type].na, false);
        for (int natom = 0; natom < ucell.atoms[type].na; ++natom)
        {
            const ModuleBase::Vector3<double>& tau = ucell.atoms[type].tau[natom];
            all_records.push_back(ModuleNeighbor::MpiAtomRecord(tau.x,
                                                                 tau.y,
                                                                 tau.z,
                                                                 global_index,
                                                                 std::array<int, 3>{{0, 0, 0}},
                                                                 false,
                                                                 type,
                                                                 natom));
            ++global_index;
        }
    }

    const std::vector<int> local_indices = mpi_domain_.select_local_atoms(all_records);
    std::vector<ModuleNeighbor::MpiAtomRecord> local_records;
    local_records.reserve(local_indices.size());
    for (const int index: local_indices)
    {
        local_records.push_back(all_records[index]);
        local_center_mask[all_records[index].type][all_records[index].natom] = true;
    }

    const std::vector<ModuleNeighbor::MpiAtomRecord> ghosts
        = mpi_domain_.exchange_ghost_atoms(local_records, stats);
    const std::size_t initial_record_count = local_records.size() + ghosts.size();
    if (initial_record_count
        > static_cast<std::size_t>(std::numeric_limits<int>::max()))
    {
        throw std::overflow_error("Grid MPI atom count exceeds the supported integer index range.");
    }
    atom_pack.reserve(static_cast<int>(initial_record_count));
    for (const ModuleNeighbor::MpiAtomRecord& record: local_records)
    {
        atom_pack.append_mpi_record(record);
    }
    for (const ModuleNeighbor::MpiAtomRecord& record: ghosts)
    {
        atom_pack.append_mpi_record(record);
    }

    // Expand both local and received records with the UnitCell image layers.
    // This covers self-rank images and combinations between a communicated
    // shift and another periodic direction in multi-dimensional topologies.
    Append_Periodic_Images(ucell, local_records);
    Append_Periodic_Images(ucell, ghosts);

    if (!atom_pack.empty())
    {
        grid_storage = ModuleNeighbor::build_grid_storage_from_atom_pack(atom_pack, box_edge_length);
        // MPI uses the same Half14 path as replicated Grid. Full27 remains
        // available through the replicated regression interface.
        neighbor_pairs = ModuleNeighbor::build_neighbor_pairs_14(atom_pack, grid_storage, sradius);
        x_min = grid_storage.x_min;
        y_min = grid_storage.y_min;
        z_min = grid_storage.z_min;
        x_max = grid_storage.x_max;
        y_max = grid_storage.y_max;
        z_max = grid_storage.z_max;
        box_nx = grid_storage.box_nx;
        box_ny = grid_storage.box_ny;
        box_nz = grid_storage.box_nz;
    }
    else
    {
        grid_storage.clear();
        neighbor_pairs.clear();
        x_min = y_min = z_min = 0.0;
        x_max = y_max = z_max = 0.0;
        box_nx = box_ny = box_nz = 0;
    }
    neighbor_pairs_27.clear();
    Build_Paged_Neighbor_List(ucell);
    mpi_mode_ = true;

    ModuleBase::GlobalFunc::OUT(ofs_in, "MPI local atoms", static_cast<int>(local_records.size()));
    ModuleBase::GlobalFunc::OUT(ofs_in, "MPI ghost atoms", static_cast<int>(ghosts.size()));
    ModuleBase::timer::end("Grid", "init_mpi");
}

bool Grid::is_local_center(const int type, const int natom) const
{
    return type >= 0 && type < static_cast<int>(local_center_mask.size()) && natom >= 0
           && natom < static_cast<int>(local_center_mask[type].size()) && local_center_mask[type][natom];
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

    this->box_nx = box_count_from_range(this->x_min, this->x_max, box_edge_length);
    this->box_ny = box_count_from_range(this->y_min, this->y_max, box_edge_length);
    this->box_nz = box_count_from_range(this->z_min, this->z_max, box_edge_length);
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
    // neighbor_pairs is the active query list. Build the full 27-box list only
    // when the caller requests an independent correctness reference.
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

    local_center_mask.resize(ucell.ntype);
    for (int type = 0; type < ucell.ntype; ++type)
    {
        local_center_mask[type].assign(ucell.atoms[type].na, true);
    }
    Build_Paged_Neighbor_List(ucell);
}

void Grid::Build_Paged_Neighbor_List(const UnitCell& ucell)
{
    std::vector<int> atom_counts(ucell.ntype, 0);
    for (int type = 0; type < ucell.ntype; ++type)
    {
        atom_counts[type] = ucell.atoms[type].na;
    }
    paged_neighbor_list.build(neighbor_pairs, atom_counts);
}

void Grid::Append_Periodic_Images(
    const UnitCell& ucell,
    const std::vector<ModuleNeighbor::MpiAtomRecord>& records)
{
    if (!pbc || records.empty())
    {
        return;
    }

    const std::vector<std::array<int, 3>> image_shifts
        = ModuleNeighbor::build_periodic_image_shifts(ucell, sradius, true);

    const std::array<std::array<double, 3>, 3> lattice_vectors{{
        std::array<double, 3>{{ucell.latvec.e11, ucell.latvec.e12, ucell.latvec.e13}},
        std::array<double, 3>{{ucell.latvec.e21, ucell.latvec.e22, ucell.latvec.e23}},
        std::array<double, 3>{{ucell.latvec.e31, ucell.latvec.e32, ucell.latvec.e33}},
    }};
    std::set<std::tuple<int, int, int, int, int, int>> existing;
    for (int index = 0; index < atom_pack.size(); ++index)
    {
        existing.insert(std::make_tuple(atom_pack.global_index[index],
                                        atom_pack.type[index],
                                        atom_pack.natom[index],
                                        atom_pack.cell_x[index],
                                        atom_pack.cell_y[index],
                                        atom_pack.cell_z[index]));
    }

    for (const ModuleNeighbor::MpiAtomRecord& record: records)
    {
        for (const std::array<int, 3>& extra_shift: image_shifts)
        {
            if (extra_shift[0] == 0 && extra_shift[1] == 0 && extra_shift[2] == 0)
            {
                continue;
            }
            ModuleNeighbor::MpiAtomRecord image = record;
            for (int axis = 0; axis < 3; ++axis)
            {
                image.x += static_cast<double>(extra_shift[axis]) * lattice_vectors[axis][0];
                image.y += static_cast<double>(extra_shift[axis]) * lattice_vectors[axis][1];
                image.z += static_cast<double>(extra_shift[axis]) * lattice_vectors[axis][2];
                image.pbc_shift[axis] += extra_shift[axis];
            }
            image.is_ghost = true;
            const auto key = std::make_tuple(image.global_index,
                                             image.type,
                                             image.natom,
                                             image.pbc_shift[0],
                                             image.pbc_shift[1],
                                             image.pbc_shift[2]);
            if (existing.insert(key).second)
            {
                atom_pack.append_mpi_record(image);
            }
        }
    }
}
