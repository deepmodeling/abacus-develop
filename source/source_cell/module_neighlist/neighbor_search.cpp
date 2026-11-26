#include "source_cell/module_neighlist/neighbor_search.h"
#include "source_cell/mdcell.h"
#include "source_cell/unitcell.h"

#include <cmath>
#include <algorithm>
#include <limits>
#include <cassert>
#include <array>
#include <stdexcept>
#include <unordered_map>
#include <vector>

// ========== Getter methods ==========

double NeighborSearch::get_search_radius() const {
    return search_radius_;
}

const std::vector<NeighborAtom>& NeighborSearch::get_all_atoms() const {
    return all_atoms_;
}

const std::vector<NeighborAtom>& NeighborSearch::get_inside_atoms() const {
    return inside_atoms_;
}

const std::vector<NeighborAtom>& NeighborSearch::get_ghost_atoms() const {
    return ghost_atoms_;
}

NeighborList& NeighborSearch::get_neighbor_list() {
    return neighbor_list_;
}

const NeighborList& NeighborSearch::get_neighbor_list() const {
    return neighbor_list_;
}

// ========== Main public interface ==========

void NeighborSearch::init_from_mdcell_(const MDCell& cell, double sr)
{
    inside_atoms_.clear();
    ghost_atoms_.clear();
    all_atoms_.clear();
    bin_manager_.clear();

    search_radius_ = sr / cell.lat0();

    const int owned_count = cell.nowned_atoms();
    const int ghost_count = cell.nghost();
    const int total_atoms = owned_count + ghost_count;

    all_atoms_.reserve(total_atoms);
    inside_atoms_.reserve(owned_count);
    ghost_atoms_.reserve(ghost_count);

    int atom_id = 0;
    for (int iat = 0; iat < owned_count; ++iat)
    {
        const LocalAtom& local = cell.owned_atoms()[iat];
        NeighborAtom atom(local.cart.x,
                          local.cart.y,
                          local.cart.z,
                          local.type,
                          local.type_index,
                          atom_id,
                          local.owner_rank);
        all_atoms_.push_back(atom);
        inside_atoms_.push_back(atom);
        ++atom_id;
    }

    for (int iat = 0; iat < ghost_count; ++iat)
    {
        const LocalAtom& local = cell.ghost_atoms()[iat];
        NeighborAtom atom(local.cart.x,
                          local.cart.y,
                          local.cart.z,
                          local.type,
                          local.type_index,
                          atom_id,
                          local.owner_rank);
        all_atoms_.push_back(atom);
        ghost_atoms_.push_back(atom);
        ++atom_id;
    }

    const long long requested_page_size = static_cast<long long>(total_atoms) * neighbor_reserve_factor;
    const int page_size = static_cast<int>(std::min<long long>(std::numeric_limits<int>::max(),
                                                               std::max(1LL, requested_page_size)));
    neighbor_list_.initialize(owned_count, page_size);
    candidate_neighbor_list_.initialize(owned_count, page_size);
}

void NeighborSearch::init_from_unitcell_(const UnitCell& ucell, double sr)
{
    search_radius_ = sr / ucell.lat0;

    // clear possible residual data from previous runs
    inside_atoms_.clear();
    ghost_atoms_.clear();
    all_atoms_.clear();
    bin_manager_.clear();

    int atom_id = 0;
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
                atom_id
            );
            inside_atoms_.push_back(atom);
            all_atoms_.push_back(atom);
            ++atom_id;
        }
    }

    int glayerX ;
    int glayerY ;
    int glayerZ ;

    int glayerX_minus ;
    int glayerY_minus ;
    int glayerZ_minus ;

    check_expand_condition(ucell, glayerX_minus, glayerX, glayerY_minus, glayerY, glayerZ_minus, glayerZ);
    set_member_variables(ucell, glayerX_minus, glayerX, glayerY_minus, glayerY, glayerZ_minus, glayerZ, atom_id);
    const long long requested_page_size = static_cast<long long>(atom_id) * neighbor_reserve_factor;
    const int page_size = static_cast<int>(std::min<long long>(std::numeric_limits<int>::max(),
                                                               std::max(1LL, requested_page_size)));
    neighbor_list_.initialize(ucell.nat, page_size);
    candidate_neighbor_list_.initialize(ucell.nat, page_size);
}

void NeighborSearch::init(BaseCell& cell, double sr)
{
    if (cell.kind() == BaseCell::Kind::mdcell)
    {
        MDCell& mdcell = static_cast<MDCell&>(cell);
        init_from_mdcell_(mdcell, sr);
        return;
    }

    assert(cell.kind() == BaseCell::Kind::unitcell);
    UnitCell& ucell = static_cast<UnitCell&>(cell);
    init_from_unitcell_(ucell, sr);
}

void NeighborSearch::build_neighbors()
{
    bin_manager_.init_bins(search_radius_, all_atoms_);
    bin_manager_.do_binning(all_atoms_);
    bin_manager_.build_atom_neighbors(candidate_neighbor_list_, inside_atoms_, all_atoms_);
    bin_manager_.build_atom_neighbors(neighbor_list_, inside_atoms_, all_atoms_);
}

void NeighborSearch::refresh_mdcell(const MDCell& cell, double cutoff)
{
    const int owned_count = cell.nowned_atoms();
    const int ghost_count = cell.nghost();
    const int expected = owned_count + ghost_count;
    if (static_cast<int>(all_atoms_.size()) != expected ||
        static_cast<int>(inside_atoms_.size()) != owned_count)
    {
        throw std::runtime_error("MDCell neighbor layout changed without rebuilding the candidate list.");
    }
    for (int i = 0; i < owned_count; ++i)
    {
        const LocalAtom& atom = cell.owned_atoms()[i];
        all_atoms_[i].position_x = atom.cart.x;
        all_atoms_[i].position_y = atom.cart.y;
        all_atoms_[i].position_z = atom.cart.z;
        inside_atoms_[i].position_x = atom.cart.x;
        inside_atoms_[i].position_y = atom.cart.y;
        inside_atoms_[i].position_z = atom.cart.z;
    }
    const int offset = owned_count;
    for (int i = 0; i < ghost_count; ++i)
    {
        const LocalAtom& atom = cell.ghost_atoms()[i];
        all_atoms_[offset + i].position_x = atom.cart.x;
        all_atoms_[offset + i].position_y = atom.cart.y;
        all_atoms_[offset + i].position_z = atom.cart.z;
        ghost_atoms_[i].position_x = atom.cart.x;
        ghost_atoms_[i].position_y = atom.cart.y;
        ghost_atoms_[i].position_z = atom.cart.z;
    }
    filter_candidate_neighbors_(cutoff, cell.lat0());
}

void NeighborSearch::filter_candidate_neighbors_(double cutoff, double lat0)
{
    const double cutoff2 = cutoff * cutoff;
    neighbor_list_.reset();
    std::vector<int> active;
    for (int i = 0; i < candidate_neighbor_list_.get_ncentral_atoms(); ++i)
    {
        active.clear();
        const NeighborAtom& center = all_atoms_[i];
        for (int j = 0; j < candidate_neighbor_list_.get_numneigh(i); ++j)
        {
            const int index = candidate_neighbor_list_.get_firstneigh(i)[j];
            const NeighborAtom& neighbor = all_atoms_[index];
            const double dx = (center.position_x - neighbor.position_x) * lat0;
            const double dy = (center.position_y - neighbor.position_y) * lat0;
            const double dz = (center.position_z - neighbor.position_z) * lat0;
            if (dx * dx + dy * dy + dz * dz < cutoff2)
            {
                active.push_back(index);
            }
        }
        neighbor_list_.set_neighbors(i, active);
    }
}


// ========== Internal methods ==========

double NeighborSearch::cross_product_norm(double a1, double a2, double a3,
                                          double b1, double b2, double b3)
{
    double c1 = a2 * b3 - a3 * b2;
    double c2 = a3 * b1 - a1 * b3;
    double c3 = a1 * b2 - a2 * b1;
    return sqrt(c1 * c1 + c2 * c2 + c3 * c3);
}

void NeighborSearch::check_expand_condition(const UnitCell& ucell, int& glayerX_minus, int& glayerX, int& glayerY_minus, int& glayerY, int& glayerZ_minus, int& glayerZ)
{
    const ModuleBase::Matrix3& lat = ucell.latvec;
    const double omega = ucell.omega;
    const double lat0 = ucell.lat0;
    const double lat0_cubed = lat0 * lat0 * lat0;

    double a23_norm = cross_product_norm(lat.e21, lat.e22, lat.e23, lat.e31, lat.e32, lat.e33);
    int extend_d11 = std::ceil(a23_norm * search_radius_ / omega * lat0_cubed);

    double a31_norm = cross_product_norm(lat.e31, lat.e32, lat.e33, lat.e11, lat.e12, lat.e13);
    int extend_d22 = std::ceil(a31_norm * search_radius_ / omega * lat0_cubed);

    double a12_norm = cross_product_norm(lat.e11, lat.e12, lat.e13, lat.e21, lat.e22, lat.e23);
    int extend_d33 = std::ceil(a12_norm * search_radius_ / omega * lat0_cubed);

    glayerX = extend_d11 + positive_layer_offset;
    glayerY = extend_d22 + positive_layer_offset;
    glayerZ = extend_d33 + positive_layer_offset;
    glayerX_minus = extend_d11;
    glayerY_minus = extend_d22;
    glayerZ_minus = extend_d33;
}

void NeighborSearch::set_member_variables(const UnitCell& ucell, int glayerX_minus, int glayerX, int glayerY_minus,
                                          int glayerY, int glayerZ_minus, int glayerZ, int& atom_id)
{
    ModuleBase::Vector3<double> vec1(ucell.latvec.e11, ucell.latvec.e12, ucell.latvec.e13);
    ModuleBase::Vector3<double> vec2(ucell.latvec.e21, ucell.latvec.e22, ucell.latvec.e23);
    ModuleBase::Vector3<double> vec3(ucell.latvec.e31, ucell.latvec.e32, ucell.latvec.e33);

    for (int ix = -glayerX_minus; ix < glayerX; ix++)
    {
        for (int iy = -glayerY_minus; iy < glayerY; iy++)
        {
            for (int iz = -glayerZ_minus; iz < glayerZ; iz++)
            {
                if(ix==0 && iy==0 && iz==0)
                {
                    continue;
                }
                for (int i = 0; i < ucell.ntype; i++)
                {
                    for (int j = 0; j < ucell.atoms[i].na; j++)
                    {
                        double atom_x = ucell.atoms[i].tau[j].x + vec1[0] * ix + vec2[0] * iy + vec3[0] * iz;
                        double atom_y = ucell.atoms[i].tau[j].y + vec1[1] * ix + vec2[1] * iy + vec3[1] * iz;
                        double atom_z = ucell.atoms[i].tau[j].z + vec1[2] * ix + vec2[2] * iy + vec3[2] * iz;

                        NeighborAtom atom(atom_x, atom_y, atom_z, i, j, atom_id);
                        ghost_atoms_.push_back(atom);
                        all_atoms_.push_back(atom);
                        ++atom_id;
                    }
                }
            }
        }
    }
}
