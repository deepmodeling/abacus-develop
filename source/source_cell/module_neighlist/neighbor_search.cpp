#include "source_cell/module_neighlist/neighbor_search.h"
#include <cmath>
#include <algorithm>
#include <limits>
#include <cassert>
#include "source_base/timer.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace
{
const int neighbor_setup_openmp_threshold = 100000;
}

// ========== Getter methods ==========

double NeighborSearch::get_search_radius() const {
    return search_radius_;
}

int NeighborSearch::get_x() const {
    return x_;
}

int NeighborSearch::get_y() const {
    return y_;
}

int NeighborSearch::get_z() const {
    return z_;
}

double NeighborSearch::get_wide_x() const {
    return wide_x_;
}

double NeighborSearch::get_wide_y() const {
    return wide_y_;
}

double NeighborSearch::get_wide_z() const {
    return wide_z_;
}

int NeighborSearch::get_glayerX() const {
    return glayerX_;
}

int NeighborSearch::get_glayerY() const {
    return glayerY_;
}

int NeighborSearch::get_glayerZ() const {
    return glayerZ_;
}

int NeighborSearch::get_glayerX_minus() const {
    return glayerX_minus_;
}

int NeighborSearch::get_glayerY_minus() const {
    return glayerY_minus_;
}

int NeighborSearch::get_glayerZ_minus() const {
    return glayerZ_minus_;
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

// ========== Setter methods ==========

void NeighborSearch::set_search_radius(double sr) {
    search_radius_ = sr;
}

void NeighborSearch::set_position(int x, int y, int z) {
    x_ = x;
    y_ = y;
    z_ = z;
}

void NeighborSearch::set_width(double wx, double wy, double wz) {
    wide_x_ = wx;
    wide_y_ = wy;
    wide_z_ = wz;
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

InputAtoms NeighborSearch::ucell_to_input_atoms(const AtomProvider& ucell)
{
    InputAtoms input_atoms;
    int atom_count = 0;
    assert(ucell.get_natom() > 0);

    input_atoms.x_low = input_atoms.y_low = input_atoms.z_low = std::numeric_limits<double>::max();
    input_atoms.x_high = input_atoms.y_high = input_atoms.z_high = std::numeric_limits<double>::lowest();

    for (int i = 0; i < ucell.get_ntype(); i++)
    {
        for (int j = 0; j < ucell.get_na(i); j++)
        {
            NeighborAtom atom(
                ucell.get_tau(i,j).x,
                ucell.get_tau(i,j).y,
                ucell.get_tau(i,j).z,
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

void NeighborSearch::check_expand_condition(const AtomProvider& ucell)
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

    glayerX_ = extend_d11 + positive_layer_offset;
    glayerY_ = extend_d22 + positive_layer_offset;
    glayerZ_ = extend_d33 + positive_layer_offset;
    glayerX_minus_ = extend_d11;
    glayerY_minus_ = extend_d22;
    glayerZ_minus_ = extend_d33;
}

void NeighborSearch::set_member_variables(const AtomProvider& ucell)
{
    ModuleBase::timer::start("NeighborSearch", "set_member_variables");
    all_atoms_.clear();

    ModuleBase::Vector3<double> vec1(ucell.get_latvec().e11, ucell.get_latvec().e12, ucell.get_latvec().e13);
    ModuleBase::Vector3<double> vec2(ucell.get_latvec().e21, ucell.get_latvec().e22, ucell.get_latvec().e23);
    ModuleBase::Vector3<double> vec3(ucell.get_latvec().e31, ucell.get_latvec().e32, ucell.get_latvec().e33);

    std::vector<NeighborAtom> base_atoms;
    base_atoms.reserve(ucell.get_natom());

    for (int i = 0; i < ucell.get_ntype(); i++)
    {
        for (int j = 0; j < ucell.get_na(i); j++)
        {
            base_atoms.push_back(NeighborAtom(
                ucell.get_tau(i,j).x,
                ucell.get_tau(i,j).y,
                ucell.get_tau(i,j).z,
                i,
                j,
                static_cast<int>(base_atoms.size())
            ));
        }
    }

    const int atoms_per_image = static_cast<int>(base_atoms.size());
    const int nimage_x = glayerX_ + glayerX_minus_;
    const int nimage_y = glayerY_ + glayerY_minus_;
    const int nimage_z = glayerZ_ + glayerZ_minus_;
    const int nimages = nimage_x * nimage_y * nimage_z;
    const int total_atoms = nimages * atoms_per_image;

#ifdef _OPENMP
    const bool use_parallel = total_atoms >= neighbor_setup_openmp_threshold && omp_get_max_threads() > 1;
    if (use_parallel)
    {
        all_atoms_.assign(total_atoms, NeighborAtom(0.0, 0.0, 0.0, 0, 0, 0));

#pragma omp parallel for schedule(static)
        for (int image = 0; image < nimages; image++)
        {
            const int image_z = image % nimage_z;
            const int image_y = (image / nimage_z) % nimage_y;
            const int image_x = image / (nimage_y * nimage_z);
            const int ix = image_x - glayerX_minus_;
            const int iy = image_y - glayerY_minus_;
            const int iz = image_z - glayerZ_minus_;
            const bool is_inside = ix == 0 && iy == 0 && iz == 0;
            const double shift_x = vec1[0] * ix + vec2[0] * iy + vec3[0] * iz;
            const double shift_y = vec1[1] * ix + vec2[1] * iy + vec3[1] * iz;
            const double shift_z = vec1[2] * ix + vec2[2] * iy + vec3[2] * iz;
            const int atom_offset = image * atoms_per_image;

            for (int atom = 0; atom < atoms_per_image; atom++)
            {
                const NeighborAtom& base_atom = base_atoms[atom];
                const int atom_id = atom_offset + atom;

                all_atoms_[atom_id] = NeighborAtom(
                    base_atom.position_x + shift_x,
                    base_atom.position_y + shift_y,
                    base_atom.position_z + shift_z,
                    base_atom.atom_type,
                    base_atom.atom_index,
                    atom_id
                );
                all_atoms_[atom_id].is_inside = is_inside;
            }
        }
        ModuleBase::timer::end("NeighborSearch", "set_member_variables");
        return;
    }
#endif

    all_atoms_.reserve(total_atoms);
    for (int image = 0; image < nimages; image++)
    {
        const int image_z = image % nimage_z;
        const int image_y = (image / nimage_z) % nimage_y;
        const int image_x = image / (nimage_y * nimage_z);
        const int ix = image_x - glayerX_minus_;
        const int iy = image_y - glayerY_minus_;
        const int iz = image_z - glayerZ_minus_;
        const bool is_inside = ix == 0 && iy == 0 && iz == 0;
        const double shift_x = vec1[0] * ix + vec2[0] * iy + vec3[0] * iz;
        const double shift_y = vec1[1] * ix + vec2[1] * iy + vec3[1] * iz;
        const double shift_z = vec1[2] * ix + vec2[2] * iy + vec3[2] * iz;
        const int atom_offset = image * atoms_per_image;

        for (int atom = 0; atom < atoms_per_image; atom++)
        {
            const NeighborAtom& base_atom = base_atoms[atom];
            const int atom_id = atom_offset + atom;

            all_atoms_.push_back(NeighborAtom(
                base_atom.position_x + shift_x,
                base_atom.position_y + shift_y,
                base_atom.position_z + shift_z,
                base_atom.atom_type,
                base_atom.atom_index,
                atom_id
            ));
            all_atoms_.back().is_inside = is_inside;
        }
    }
    ModuleBase::timer::end("NeighborSearch", "set_member_variables");
}

// ========== Main public interface ==========

void NeighborSearch::init(const AtomProvider& ucell, double sr, int mpi_rank)
{
    ModuleBase::timer::start("NeighborSearch", "init");
    // clear possible residual data from previous runs
    inside_atoms_.clear();
    ghost_atoms_.clear();
    all_atoms_.clear();
    // clear any existing bin manager state
    bin_manager_.clear();

    search_radius_ = sr / ucell.get_lat0();
    check_expand_condition(ucell);
    set_member_variables(ucell);
    InputAtoms atoms = ucell_to_input_atoms(ucell);

    int mpi_size = 1;
    int nx, ny, nz;
    decompose(mpi_size, nx, ny, nz);

    z_ = mpi_rank / (nx * ny);
    y_ = (mpi_rank % (nx * ny)) / nx;
    x_ = mpi_rank % (nx * ny) % nx;

    wide_x_ = (atoms.x_high - atoms.x_low) / nx;
    wide_y_ = (atoms.y_high - atoms.y_low) / ny;
    wide_z_ = (atoms.z_high - atoms.z_low) / nz;
    assert(wide_x_ >= 0);
    assert(wide_y_ >= 0);
    assert(wide_z_ >= 0);

    auto classify_atom = [&](const NeighborAtom& atom)
    {
        int in_x;
        int in_y;
        int in_z;

        if(wide_x_ < coord_tolerance)
        {
            if(std::abs(atom.position_x - atoms.x_low) < coord_tolerance)
            {
                in_x = x_;
            }
            else
            {
                in_x = std::numeric_limits<int>::max();
            }
        }
        else
        {
            in_x = std::min(
                static_cast<int>(std::floor((atom.position_x - atoms.x_low) / wide_x_)),
                nx - 1
            );
        }
        if(wide_y_ < coord_tolerance)
        {
            if(std::abs(atom.position_y - atoms.y_low) < coord_tolerance)
            {
                in_y = y_;
            }
            else
            {
                in_y = std::numeric_limits<int>::max();
            }
        }
        else
        {
            in_y = std::min(
                static_cast<int>(std::floor((atom.position_y - atoms.y_low) / wide_y_)),
                ny - 1
            );
        }
        if(wide_z_ < coord_tolerance)
        {
            if(std::abs(atom.position_z - atoms.z_low) < coord_tolerance)
            {
                in_z = z_;
            }
            else
            {
                in_z = std::numeric_limits<int>::max();
            }
        }
        else
        {
            in_z = std::min(
                static_cast<int>(std::floor((atom.position_z - atoms.z_low) / wide_z_)),
                nz - 1
            );
        }

        if (in_x == x_ && in_y == y_ && in_z == z_ &&
            atom.position_x <= atoms.x_high &&
            atom.position_y <= atoms.y_high &&
            atom.position_z <= atoms.z_high &&
            atom.is_inside)
        {
            return 1;
        }
        else if (distance(
            atom.position_x,
            atom.position_y,
            atom.position_z,
            atoms.x_low,
            atoms.y_low,
            atoms.z_low) <= search_radius_ * search_radius_)
        {
            return 2;
        }

        return 0;
    };

#ifdef _OPENMP
    const bool use_parallel_classification =
        static_cast<int>(all_atoms_.size()) >= neighbor_setup_openmp_threshold && omp_get_max_threads() > 1;
    if (use_parallel_classification)
    {
        std::vector<unsigned char> categories(all_atoms_.size(), 0);
        int ninside = 0;
        int nghost = 0;

#pragma omp parallel for schedule(static) reduction(+:ninside, nghost)
        for (int i = 0; i < static_cast<int>(all_atoms_.size()); i++)
        {
            const int category = classify_atom(all_atoms_[i]);
            categories[i] = static_cast<unsigned char>(category);
            ninside += category == 1;
            nghost += category == 2;
        }

        inside_atoms_.reserve(ninside);
        ghost_atoms_.reserve(nghost);

        int inside_index = 0;
        int ghost_index = 0;
        for (size_t i = 0; i < all_atoms_.size(); i++)
        {
            if (categories[i] == 1)
            {
                inside_atoms_.push_back(all_atoms_[i]);
                inside_index++;
            }
            else if (categories[i] == 2)
            {
                ghost_atoms_.push_back(all_atoms_[i]);
                ghost_index++;
            }
        }
        assert(inside_index == ninside);
        assert(ghost_index == nghost);
    }
    else
#endif
    {
        inside_atoms_.reserve(ucell.get_natom());
        ghost_atoms_.reserve(all_atoms_.size());
        for (size_t i = 0; i < all_atoms_.size(); i++)
        {
            const int category = classify_atom(all_atoms_[i]);
            if (category == 1)
            {
                inside_atoms_.push_back(all_atoms_[i]);
            }
            else if (category == 2)
            {
                ghost_atoms_.push_back(all_atoms_[i]);
            }
        }
    }

    neighbor_list_.initialize(inside_atoms_.size(), all_atoms_.size() * neighbor_reserve_factor);
    ModuleBase::timer::end("NeighborSearch", "init");
}

void NeighborSearch::build_neighbors()
{
    ModuleBase::timer::start("NeighborSearch", "build_neighbors");
    bin_manager_.init_bins(search_radius_, inside_atoms_, ghost_atoms_);
    bin_manager_.do_binning(inside_atoms_, ghost_atoms_);
    bin_manager_.build_atom_neighbors(neighbor_list_, inside_atoms_);
    ModuleBase::timer::end("NeighborSearch", "build_neighbors");
}

// ========== Utility methods ==========

double NeighborSearch::distance(
    double position_x,
    double position_y,
    double position_z,
    double x_low,
    double y_low,
    double z_low)
{
    double dx = std::max(0.0, std::max(x_low + x_ * wide_x_ - position_x, position_x - (x_low + (x_ + 1) * wide_x_)));
    double dy = std::max(0.0, std::max(y_low + y_ * wide_y_ - position_y, position_y - (y_low + (y_ + 1) * wide_y_)));
    double dz = std::max(0.0, std::max(z_low + z_ * wide_z_ - position_z, position_z - (z_low + (z_ + 1) * wide_z_)));
    return dx * dx + dy * dy + dz * dz;
}

void NeighborSearch::decompose(int mpi_size, int &nx, int &ny, int &nz)
{
    nx = 1;
    ny = 1;
    nz = mpi_size;

    int cube = static_cast<int>(cbrt(mpi_size));
    for (int i = cube; i >= 1; i--)
    {
        if (mpi_size % i == 0)
        {
            nx = i;
            ny = mpi_size / i;
            break;
        }
    }

    int sq = static_cast<int>(sqrt(ny));
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
