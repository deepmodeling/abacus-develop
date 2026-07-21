#include "source_cell/module_neighlist/neighbor_search.h"
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

void NeighborSearch::init_distributed(const std::vector<LocalAtom>& owned_atoms,
                                      const std::vector<LocalAtom>& ghost_atoms,
                                      double sr,
                                      double lat0)
{
    inside_atoms_.clear();
    ghost_atoms_.clear();
    all_atoms_.clear();
    bin_manager_.clear();

    search_radius_ = sr / lat0;

    const std::size_t total_atoms = ModuleNeighList::checked_size_sum(owned_atoms.size(),
                                                                       ghost_atoms.size(),
                                                                       "NeighborSearch distributed atom count");
    if (total_atoms > static_cast<std::size_t>(std::numeric_limits<ModuleNeighList::LocalAtomIndex>::max()))
    {
        throw std::overflow_error("NeighborSearch distributed atom count exceeds local atom index range.");
    }

    all_atoms_.reserve(total_atoms);
    inside_atoms_.reserve(owned_atoms.size());
    ghost_atoms_.reserve(ghost_atoms.size());

    for (size_t iat = 0; iat < owned_atoms.size(); ++iat)
    {
        const LocalAtom& local = owned_atoms[iat];
        NeighborAtom atom(local.cart.x,
                          local.cart.y,
                          local.cart.z,
                          local.type,
                          local.type_index,
                          ModuleNeighList::checked_local_atom_index(all_atoms_.size(),
                                                                    "NeighborSearch owned atom id"),
                          local.global_id,
                          local.owner_rank);
        all_atoms_.push_back(atom);
        inside_atoms_.push_back(atom);
    }

    for (size_t iat = 0; iat < ghost_atoms.size(); ++iat)
    {
        const LocalAtom& local = ghost_atoms[iat];
        NeighborAtom atom(local.cart.x,
                          local.cart.y,
                          local.cart.z,
                          local.type,
                          local.type_index,
                          ModuleNeighList::checked_local_atom_index(all_atoms_.size(),
                                                                    "NeighborSearch ghost atom id"),
                          local.global_id,
                          local.owner_rank);
        all_atoms_.push_back(atom);
        ghost_atoms_.push_back(atom);
    }

    const std::size_t page_size = ModuleNeighList::checked_size_product(all_atoms_.size(),
                                                                        neighbor_reserve_factor,
                                                                        "NeighborSearch page size");
    neighbor_list_.initialize(inside_atoms_.size(), page_size);
}

void NeighborSearch::init(const AtomProvider& ucell, double sr)
{
    search_radius_ = sr / ucell.get_lat0();

    // clear possible residual data from previous runs
    inside_atoms_.clear();
    ghost_atoms_.clear();
    all_atoms_.clear();
    bin_manager_.clear();

    for (int i = 0; i < ucell.get_ntype(); i++)
    {
        for (int j = 0; j < ucell.get_na(i); j++)
        {
            const ModuleNeighList::LocalAtomIndex atom_count
                = ModuleNeighList::checked_local_atom_index(all_atoms_.size(),
                                                            "NeighborSearch atom id");
            NeighborAtom atom(
                ucell.get_tau(i,j).x,
                ucell.get_tau(i,j).y,
                ucell.get_tau(i,j).z,
                i,
                j,
                atom_count
            );
            inside_atoms_.push_back(atom);
            all_atoms_.push_back(atom);
        }
    }

    int glayerX ;
    int glayerY ;
    int glayerZ ;

    int glayerX_minus ;
    int glayerY_minus ;
    int glayerZ_minus ;

    check_expand_condition(ucell, glayerX_minus, glayerX, glayerY_minus, glayerY, glayerZ_minus, glayerZ);
    set_member_variables(ucell, glayerX_minus, glayerX, glayerY_minus, glayerY, glayerZ_minus, glayerZ);
    const std::size_t page_size = ModuleNeighList::checked_size_product(all_atoms_.size(),
                                                                        neighbor_reserve_factor,
                                                                        "NeighborSearch page size");
    neighbor_list_.initialize(inside_atoms_.size(), page_size);
}

void NeighborSearch::build_neighbors()
{
    bin_manager_.init_bins(search_radius_, all_atoms_);
    bin_manager_.do_binning(all_atoms_);
    bin_manager_.build_atom_neighbors(neighbor_list_, inside_atoms_, all_atoms_);
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

void NeighborSearch::check_expand_condition(const AtomProvider& ucell, int& glayerX_minus, int& glayerX, int& glayerY_minus, int& glayerY, int& glayerZ_minus, int& glayerZ)
{
    const auto& lat = ucell.get_latvec();
    const double omega = ucell.get_omega();
    const double lat0 = ucell.get_lat0();
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

void NeighborSearch::set_member_variables(const AtomProvider& ucell, int glayerX_minus, int glayerX, int glayerY_minus, int glayerY, int glayerZ_minus, int glayerZ)
{
    ModuleBase::Vector3<double> vec1(ucell.get_latvec().e11, ucell.get_latvec().e12, ucell.get_latvec().e13);
    ModuleBase::Vector3<double> vec2(ucell.get_latvec().e21, ucell.get_latvec().e22, ucell.get_latvec().e23);
    ModuleBase::Vector3<double> vec3(ucell.get_latvec().e31, ucell.get_latvec().e32, ucell.get_latvec().e33);

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
                for (int i = 0; i < ucell.get_ntype(); i++)
                {
                    for (int j = 0; j < ucell.get_na(i); j++)
                    {
                        double atom_x = ucell.get_tau(i,j).x + vec1[0] * ix + vec2[0] * iy + vec3[0] * iz;
                        double atom_y = ucell.get_tau(i,j).y + vec1[1] * ix + vec2[1] * iy + vec3[1] * iz;
                        double atom_z = ucell.get_tau(i,j).z + vec1[2] * ix + vec2[2] * iy + vec3[2] * iz;

                        const ModuleNeighList::LocalAtomIndex atom_count
                            = ModuleNeighList::checked_local_atom_index(all_atoms_.size(),
                                                                        "NeighborSearch atom id");
                        NeighborAtom atom(atom_x, atom_y, atom_z, i, j, atom_count);
                        ghost_atoms_.push_back(atom);
                        all_atoms_.push_back(atom);
                    }
                }
            }
        }
    }
}
