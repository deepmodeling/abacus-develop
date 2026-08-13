#include "source_cell/module_neighlist/neighbor_search.h"
#include <cmath>
#include <algorithm>
#include <cstdint>
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
            const ModuleBase::Vector3<double> position = ucell.get_tau(i, j);
            const ModuleNeighList::LocalAtomIndex atom_count
                = ModuleNeighList::checked_local_atom_index(all_atoms_.size(),
                                                            "NeighborSearch atom id");
            NeighborAtom atom(
                position.x,
                position.y,
                position.z,
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
    if (glayerX_minus < 0 || glayerX < 0 || glayerY_minus < 0 || glayerY < 0
        || glayerZ_minus < 0 || glayerZ < 0)
    {
        throw std::invalid_argument("NeighborSearch periodic image layers must be non-negative.");
    }

    const ModuleBase::Matrix3& lattice = ucell.get_latvec();
    const ModuleBase::Vector3<double> vec1(lattice.e11, lattice.e12, lattice.e13);
    const ModuleBase::Vector3<double> vec2(lattice.e21, lattice.e22, lattice.e23);
    const ModuleBase::Vector3<double> vec3(lattice.e31, lattice.e32, lattice.e33);

    const std::size_t image_count_x
        = ModuleNeighList::checked_size_sum(static_cast<std::size_t>(glayerX_minus),
                                            static_cast<std::size_t>(glayerX),
                                            "NeighborSearch x image count");
    const std::size_t image_count_y
        = ModuleNeighList::checked_size_sum(static_cast<std::size_t>(glayerY_minus),
                                            static_cast<std::size_t>(glayerY),
                                            "NeighborSearch y image count");
    const std::size_t image_count_z
        = ModuleNeighList::checked_size_sum(static_cast<std::size_t>(glayerZ_minus),
                                            static_cast<std::size_t>(glayerZ),
                                            "NeighborSearch z image count");
    const std::size_t image_count_yz
        = ModuleNeighList::checked_size_product(image_count_y,
                                                image_count_z,
                                                "NeighborSearch yz image count");
    const std::size_t image_count
        = ModuleNeighList::checked_size_product(image_count_x,
                                                image_count_yz,
                                                "NeighborSearch periodic image count");
    if (image_count == 0)
    {
        throw std::invalid_argument("NeighborSearch periodic image range must contain the primary cell.");
    }

    const std::size_t origin_xy
        = ModuleNeighList::checked_size_sum(
            ModuleNeighList::checked_size_product(static_cast<std::size_t>(glayerX_minus),
                                                  image_count_y,
                                                  "NeighborSearch primary image index"),
            static_cast<std::size_t>(glayerY_minus),
            "NeighborSearch primary image index");
    const std::size_t origin_image
        = ModuleNeighList::checked_size_sum(
            ModuleNeighList::checked_size_product(origin_xy,
                                                  image_count_z,
                                                  "NeighborSearch primary image index"),
            static_cast<std::size_t>(glayerZ_minus),
            "NeighborSearch primary image index");

    const std::size_t base_atom_count = inside_atoms_.size();
    const std::size_t ghost_image_count = image_count - 1;
    const std::size_t generated_atom_count
        = ModuleNeighList::checked_size_product(ghost_image_count,
                                                base_atom_count,
                                                "NeighborSearch generated atom count");
    const std::size_t final_atom_count
        = ModuleNeighList::checked_size_sum(base_atom_count,
                                            generated_atom_count,
                                            "NeighborSearch total atom count");
    if (final_atom_count > static_cast<std::size_t>(std::numeric_limits<ModuleNeighList::LocalAtomIndex>::max()))
    {
        throw std::overflow_error("NeighborSearch total atom count exceeds local atom index range.");
    }
    if (generated_atom_count == 0)
    {
        return;
    }

    const NeighborAtom placeholder(0.0, 0.0, 0.0, 0, 0, 0);
    all_atoms_.insert(all_atoms_.end(), generated_atom_count, placeholder);
    ghost_atoms_.assign(generated_atom_count, placeholder);

    // Thread creation costs more than it saves for small unit cells.
    const std::size_t parallel_threshold = 100000;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(generated_atom_count >= parallel_threshold)
#endif
    for (std::int64_t generated_index = 0;
         generated_index < static_cast<std::int64_t>(generated_atom_count);
         ++generated_index)
    {
        const std::size_t output_index = static_cast<std::size_t>(generated_index);
        const std::size_t ghost_image = output_index / base_atom_count;
        const std::size_t base_atom = output_index % base_atom_count;
        const std::size_t image = ghost_image < origin_image ? ghost_image : ghost_image + 1;
        const std::size_t image_x = image / image_count_yz;
        const std::size_t image_yz = image % image_count_yz;
        const std::size_t image_y = image_yz / image_count_z;
        const std::size_t image_z = image_yz % image_count_z;
        const double ix = static_cast<double>(image_x) - static_cast<double>(glayerX_minus);
        const double iy = static_cast<double>(image_y) - static_cast<double>(glayerY_minus);
        const double iz = static_cast<double>(image_z) - static_cast<double>(glayerZ_minus);

        const NeighborAtom& source = inside_atoms_[base_atom];
        const ModuleNeighList::LocalAtomIndex atom_id
            = static_cast<ModuleNeighList::LocalAtomIndex>(base_atom_count + output_index);
        const NeighborAtom atom(source.position_x + vec1[0] * ix + vec2[0] * iy + vec3[0] * iz,
                                source.position_y + vec1[1] * ix + vec2[1] * iy + vec3[1] * iz,
                                source.position_z + vec1[2] * ix + vec2[2] * iy + vec3[2] * iz,
                                source.atom_type,
                                source.atom_index,
                                atom_id);
        ghost_atoms_[output_index] = atom;
        all_atoms_[base_atom_count + output_index] = atom;
    }
}
