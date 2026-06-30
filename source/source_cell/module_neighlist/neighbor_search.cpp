#include "source_cell/module_neighlist/neighbor_search.h"
#include <cmath>
#include <algorithm>
#include <limits>
#include <cassert>
#include <array>
#include <unordered_map>
#include <vector>

namespace
{
struct OriginalAtom
{
    OriginalAtom(const std::array<double, 3>& frac_in, int atom_type_in, int atom_index_in)
        : frac(frac_in), atom_type(atom_type_in), atom_index(atom_index_in)
    {
    }

    std::array<double, 3> frac;
    int atom_type = 0;
    int atom_index = 0;
};

struct PeriodicInterval
{
    PeriodicInterval(double lo_in, double hi_in, int shift_in) : lo(lo_in), hi(hi_in), shift(shift_in)
    {
    }

    double lo = 0.0;
    double hi = 0.0;
    int shift = 0;
};

struct FractionalDomain
{
    FractionalDomain(const std::array<double, 3>& lo_in, const std::array<double, 3>& hi_in)
        : lo(lo_in), hi(hi_in)
    {
    }

    std::array<double, 3> lo;
    std::array<double, 3> hi;
};

double dot_product(const ModuleBase::Vector3<double>& a, const ModuleBase::Vector3<double>& b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

ModuleBase::Vector3<double> cross_product(const ModuleBase::Vector3<double>& a,
                                          const ModuleBase::Vector3<double>& b)
{
    return ModuleBase::Vector3<double>(a.y * b.z - a.z * b.y,
                                       a.z * b.x - a.x * b.z,
                                       a.x * b.y - a.y * b.x);
}

double norm(const ModuleBase::Vector3<double>& v)
{
    return std::sqrt(dot_product(v, v));
}

double wrap_fractional(double value)
{
    value -= std::floor(value);
    if (value >= 1.0 - 1.0e-12)
    {
        return 0.0;
    }
    if (value < 1.0e-12)
    {
        return 0.0;
    }
    return value;
}

int clamp_index(int value, int low, int high)
{
    return std::min(std::max(value, low), high);
}

int fractional_domain_index(double frac, int n)
{
    return clamp_index(static_cast<int>(std::floor(frac * n)), 0, n - 1);
}

long long bin_key(int ix, int iy, int iz, const std::array<int, 3>& nbin)
{
    return (static_cast<long long>(ix) * nbin[1] + iy) * nbin[2] + iz;
}

std::vector<PeriodicInterval> split_periodic_interval(double lo, double hi)
{
    std::vector<PeriodicInterval> intervals;
    if (hi <= lo)
    {
        return intervals;
    }

    const int first_shift = static_cast<int>(std::floor(lo));
    const int last_shift = static_cast<int>(std::ceil(hi)) - 1;
    for (int shift = first_shift; shift <= last_shift; ++shift)
    {
        const double local_lo = std::max(0.0, lo - shift);
        const double local_hi = std::min(1.0, hi - shift);
        if (local_lo < local_hi)
        {
            intervals.push_back(PeriodicInterval(local_lo, local_hi, shift));
        }
    }
    return intervals;
}

bool inside_interval(double value, double lo, double hi)
{
    return value >= lo && value < hi;
}

bool inside_block(const OriginalAtom& atom,
                  const PeriodicInterval& bx,
                  const PeriodicInterval& by,
                  const PeriodicInterval& bz)
{
    return inside_interval(atom.frac[0], bx.lo, bx.hi) &&
           inside_interval(atom.frac[1], by.lo, by.hi) &&
           inside_interval(atom.frac[2], bz.lo, bz.hi);
}
} // namespace

constexpr double NeighborSearch::coord_tolerance;

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
    all_atoms_.clear();

    ModuleBase::Vector3<double> vec1(ucell.get_latvec().e11, ucell.get_latvec().e12, ucell.get_latvec().e13);
    ModuleBase::Vector3<double> vec2(ucell.get_latvec().e21, ucell.get_latvec().e22, ucell.get_latvec().e23);
    ModuleBase::Vector3<double> vec3(ucell.get_latvec().e31, ucell.get_latvec().e32, ucell.get_latvec().e33);

    int atom_count = 0;

    for (int ix = -glayerX_minus_; ix < glayerX_; ix++)
    {
        for (int iy = -glayerY_minus_; iy < glayerY_; iy++)
        {
            for (int iz = -glayerZ_minus_; iz < glayerZ_; iz++)
            {
                for (int i = 0; i < ucell.get_ntype(); i++)
                {
                    for (int j = 0; j < ucell.get_na(i); j++)
                    {
                        double atom_x = ucell.get_tau(i,j).x + vec1[0] * ix + vec2[0] * iy + vec3[0] * iz;
                        double atom_y = ucell.get_tau(i,j).y + vec1[1] * ix + vec2[1] * iy + vec3[1] * iz;
                        double atom_z = ucell.get_tau(i,j).z + vec1[2] * ix + vec2[2] * iy + vec3[2] * iz;

                        NeighborAtom atom(atom_x, atom_y, atom_z, i, j, atom_count);
                        if(ix==0 && iy==0 && iz==0)
                        {
                            atom.is_inside = true;
                        }
                        else
                        {
                            atom.is_inside = false;
                        }
                        all_atoms_.push_back(atom);
                        atom_count++;
                    }
                }
            }
        }
    }
}

void NeighborSearch::set_local_member_variables(const AtomProvider& ucell,
                                                const InputAtoms& atoms,
                                                int nx,
                                                int ny,
                                                int nz)
{
    all_atoms_.clear();
    inside_atoms_.clear();
    ghost_atoms_.clear();

    ModuleBase::Vector3<double> vec1(ucell.get_latvec().e11, ucell.get_latvec().e12, ucell.get_latvec().e13);
    ModuleBase::Vector3<double> vec2(ucell.get_latvec().e21, ucell.get_latvec().e22, ucell.get_latvec().e23);
    ModuleBase::Vector3<double> vec3(ucell.get_latvec().e31, ucell.get_latvec().e32, ucell.get_latvec().e33);

    const auto domain_index = [](double position, double low, double wide, int n) {
        if (wide < coord_tolerance)
        {
            return std::abs(position - low) < coord_tolerance ? 0 : std::numeric_limits<int>::max();
        }
        return std::min(std::max(static_cast<int>(std::floor((position - low) / wide)), 0), n - 1);
    };

    const double search_radius2 = search_radius_ * search_radius_;
    for (int ix = -glayerX_minus_; ix < glayerX_; ix++)
    {
        for (int iy = -glayerY_minus_; iy < glayerY_; iy++)
        {
            for (int iz = -glayerZ_minus_; iz < glayerZ_; iz++)
            {
                const bool central_image = (ix == 0 && iy == 0 && iz == 0);
                for (int it = 0; it < ucell.get_ntype(); it++)
                {
                    for (int ia = 0; ia < ucell.get_na(it); ia++)
                    {
                        double atom_x = ucell.get_tau(it, ia).x + vec1[0] * ix + vec2[0] * iy + vec3[0] * iz;
                        double atom_y = ucell.get_tau(it, ia).y + vec1[1] * ix + vec2[1] * iy + vec3[1] * iz;
                        double atom_z = ucell.get_tau(it, ia).z + vec1[2] * ix + vec2[2] * iy + vec3[2] * iz;

                        const int in_x = domain_index(atom_x, atoms.x_low, wide_x_, nx);
                        const int in_y = domain_index(atom_y, atoms.y_low, wide_y_, ny);
                        const int in_z = domain_index(atom_z, atoms.z_low, wide_z_, nz);

                        const bool owned = central_image &&
                                           in_x == x_ &&
                                           in_y == y_ &&
                                           in_z == z_ &&
                                           atom_x <= atoms.x_high &&
                                           atom_y <= atoms.y_high &&
                                           atom_z <= atoms.z_high;
                        const bool ghost = !owned &&
                                           distance(atom_x, atom_y, atom_z, atoms.x_low, atoms.y_low, atoms.z_low)
                                               <= search_radius2;

                        if (!owned && !ghost)
                        {
                            continue;
                        }

                        NeighborAtom atom(atom_x, atom_y, atom_z, it, ia, static_cast<int>(all_atoms_.size()));
                        atom.is_inside = owned;
                        all_atoms_.push_back(atom);
                        if (owned)
                        {
                            inside_atoms_.push_back(atom);
                        }
                        else
                        {
                            ghost_atoms_.push_back(atom);
                        }
                    }
                }
            }
        }
    }
}

void NeighborSearch::set_local_member_variables_by_halo(const AtomProvider& ucell, int nx, int ny, int nz)
{
    all_atoms_.clear();
    inside_atoms_.clear();
    ghost_atoms_.clear();

    const ModuleBase::Matrix3& lat = ucell.get_latvec();
    const ModuleBase::Matrix3 inv_lat = lat.Inverse();

    const ModuleBase::Vector3<double> a1(lat.e11, lat.e12, lat.e13);
    const ModuleBase::Vector3<double> a2(lat.e21, lat.e22, lat.e23);
    const ModuleBase::Vector3<double> a3(lat.e31, lat.e32, lat.e33);

    const ModuleBase::Vector3<double> a2xa3 = cross_product(a2, a3);
    const ModuleBase::Vector3<double> a3xa1 = cross_product(a3, a1);
    const ModuleBase::Vector3<double> a1xa2 = cross_product(a1, a2);

    const double volume = std::abs(dot_product(a1, a2xa3));
    assert(volume > coord_tolerance);

    const std::array<double, 3> heights = {{
        volume / norm(a2xa3),
        volume / norm(a3xa1),
        volume / norm(a1xa2)
    }};

    std::array<double, 3> margin = {{
        search_radius_ / heights[0] + coord_tolerance,
        search_radius_ / heights[1] + coord_tolerance,
        search_radius_ / heights[2] + coord_tolerance
    }};

    const std::array<double, 3> domain_lo = {{
        static_cast<double>(x_) / nx,
        static_cast<double>(y_) / ny,
        static_cast<double>(z_) / nz
    }};
    const std::array<double, 3> domain_hi = {{
        static_cast<double>(x_ + 1) / nx,
        static_cast<double>(y_ + 1) / ny,
        static_cast<double>(z_ + 1) / nz
    }};
    const FractionalDomain domain(domain_lo, domain_hi);

    const std::array<double, 3> halo_lo = {{
        domain.lo[0] - margin[0],
        domain.lo[1] - margin[1],
        domain.lo[2] - margin[2]
    }};
    const std::array<double, 3> halo_hi = {{
        domain.hi[0] + margin[0],
        domain.hi[1] + margin[1],
        domain.hi[2] + margin[2]
    }};
    const FractionalDomain halo(halo_lo, halo_hi);

    std::vector<OriginalAtom> original_atoms;
    original_atoms.reserve(ucell.get_natom());
    for (int it = 0; it < ucell.get_ntype(); ++it)
    {
        for (int ia = 0; ia < ucell.get_na(it); ++ia)
        {
            const ModuleBase::Vector3<double> cart = ucell.get_tau(it, ia);
            const ModuleBase::Vector3<double> frac = cart * inv_lat;
            const std::array<double, 3> wrapped_frac = {{
                wrap_fractional(frac.x),
                wrap_fractional(frac.y),
                wrap_fractional(frac.z)
            }};
            original_atoms.push_back(OriginalAtom(wrapped_frac, it, ia));
        }
    }

    std::array<int, 3> nbin;
    for (int idim = 0; idim < 3; ++idim)
    {
        nbin[idim] = std::max(1, static_cast<int>(std::ceil(1.0 / std::max(margin[idim], coord_tolerance))));
    }

    std::unordered_map<long long, std::vector<int>> bins;
    bins.reserve(original_atoms.size());
    for (int iat = 0; iat < static_cast<int>(original_atoms.size()); ++iat)
    {
        const OriginalAtom& atom = original_atoms[iat];
        const int ix = clamp_index(static_cast<int>(std::floor(atom.frac[0] * nbin[0])), 0, nbin[0] - 1);
        const int iy = clamp_index(static_cast<int>(std::floor(atom.frac[1] * nbin[1])), 0, nbin[1] - 1);
        const int iz = clamp_index(static_cast<int>(std::floor(atom.frac[2] * nbin[2])), 0, nbin[2] - 1);
        bins[bin_key(ix, iy, iz, nbin)].push_back(iat);
    }

    const std::vector<PeriodicInterval> intervals_x = split_periodic_interval(halo.lo[0], halo.hi[0]);
    const std::vector<PeriodicInterval> intervals_y = split_periodic_interval(halo.lo[1], halo.hi[1]);
    const std::vector<PeriodicInterval> intervals_z = split_periodic_interval(halo.lo[2], halo.hi[2]);

    for (const PeriodicInterval& bx : intervals_x)
    {
        const int ix_begin = clamp_index(static_cast<int>(std::floor(bx.lo * nbin[0])), 0, nbin[0] - 1);
        const int ix_end = clamp_index(static_cast<int>(std::ceil(bx.hi * nbin[0])) - 1, 0, nbin[0] - 1);
        for (const PeriodicInterval& by : intervals_y)
        {
            const int iy_begin = clamp_index(static_cast<int>(std::floor(by.lo * nbin[1])), 0, nbin[1] - 1);
            const int iy_end = clamp_index(static_cast<int>(std::ceil(by.hi * nbin[1])) - 1, 0, nbin[1] - 1);
            for (const PeriodicInterval& bz : intervals_z)
            {
                const int iz_begin = clamp_index(static_cast<int>(std::floor(bz.lo * nbin[2])), 0, nbin[2] - 1);
                const int iz_end = clamp_index(static_cast<int>(std::ceil(bz.hi * nbin[2])) - 1, 0, nbin[2] - 1);

                for (int ix = ix_begin; ix <= ix_end; ++ix)
                {
                    for (int iy = iy_begin; iy <= iy_end; ++iy)
                    {
                        for (int iz = iz_begin; iz <= iz_end; ++iz)
                        {
                            const auto bin_iter = bins.find(bin_key(ix, iy, iz, nbin));
                            if (bin_iter == bins.end())
                            {
                                continue;
                            }

                            for (const int atom_id : bin_iter->second)
                            {
                                const OriginalAtom& original = original_atoms[atom_id];
                                if (!inside_block(original, bx, by, bz))
                                {
                                    continue;
                                }

                                const bool central_image = bx.shift == 0 && by.shift == 0 && bz.shift == 0;
                                const bool owned = central_image &&
                                                   fractional_domain_index(original.frac[0], nx) == x_ &&
                                                   fractional_domain_index(original.frac[1], ny) == y_ &&
                                                   fractional_domain_index(original.frac[2], nz) == z_;

                                const ModuleBase::Vector3<double> frac_image(original.frac[0] + bx.shift,
                                                                             original.frac[1] + by.shift,
                                                                             original.frac[2] + bz.shift);
                                const ModuleBase::Vector3<double> cart_image = frac_image * lat;

                                NeighborAtom atom(cart_image.x,
                                                  cart_image.y,
                                                  cart_image.z,
                                                  original.atom_type,
                                                  original.atom_index,
                                                  static_cast<int>(all_atoms_.size()));
                                atom.is_inside = owned;
                                all_atoms_.push_back(atom);
                                if (owned)
                                {
                                    inside_atoms_.push_back(atom);
                                }
                                else
                                {
                                    ghost_atoms_.push_back(atom);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

// ========== Main public interface ==========

void NeighborSearch::init(const AtomProvider& ucell, double sr, int mpi_rank)
{
    this->init(ucell, sr, mpi_rank, 1);
}

void NeighborSearch::init(const AtomProvider& ucell, double sr, int mpi_rank, int mpi_size)
{
    // clear possible residual data from previous runs
    inside_atoms_.clear();
    ghost_atoms_.clear();
    all_atoms_.clear();
    // clear any existing bin manager state
    bin_manager_.clear();

    search_radius_ = sr / ucell.get_lat0();
    check_expand_condition(ucell);
    InputAtoms atoms = ucell_to_input_atoms(ucell);

    assert(mpi_size > 0);
    assert(mpi_rank >= 0);

    const double span_x = atoms.x_high - atoms.x_low;
    const double span_y = atoms.y_high - atoms.y_low;
    const double span_z = atoms.z_high - atoms.z_low;

    int nx, ny, nz;
    decompose(mpi_size, span_x, span_y, span_z, nx, ny, nz);

    const int active_size = nx * ny * nz;
    assert(active_size > 0);
    assert(active_size <= mpi_size);

    wide_x_ = span_x / nx;
    wide_y_ = span_y / ny;
    wide_z_ = span_z / nz;
    assert(wide_x_ >= 0);
    assert(wide_y_ >= 0);
    assert(wide_z_ >= 0);

    if (mpi_rank >= active_size)
    {
        x_ = -1;
        y_ = -1;
        z_ = -1;
        neighbor_list_.initialize(0, neighbor_reserve_factor);
        return;
    }

    z_ = mpi_rank / (nx * ny);
    y_ = (mpi_rank % (nx * ny)) / nx;
    x_ = mpi_rank % nx;

    if (mpi_size > 1)
    {
        set_local_member_variables_by_halo(ucell, nx, ny, nz);
    }
    else
    {
        set_member_variables(ucell);

        int in_x, in_y, in_z;
        const auto domain_index = [](double position, double low, double wide, int n) {
            if (wide < coord_tolerance)
            {
                return std::abs(position - low) < coord_tolerance ? 0 : std::numeric_limits<int>::max();
            }
            return std::min(std::max(static_cast<int>(std::floor((position - low) / wide)), 0), n - 1);
        };

        for (size_t i = 0; i < all_atoms_.size(); i++)
        {
            in_x = domain_index(all_atoms_[i].position_x, atoms.x_low, wide_x_, nx);
            in_y = domain_index(all_atoms_[i].position_y, atoms.y_low, wide_y_, ny);
            in_z = domain_index(all_atoms_[i].position_z, atoms.z_low, wide_z_, nz);

            if (in_x == x_ && in_y == y_ && in_z == z_ &&
                all_atoms_[i].position_x <= atoms.x_high &&
                all_atoms_[i].position_y <= atoms.y_high &&
                all_atoms_[i].position_z <= atoms.z_high &&
                all_atoms_[i].is_inside)
            {
                inside_atoms_.push_back(all_atoms_[i]);
            }
            else if (distance(
                all_atoms_[i].position_x,
                all_atoms_[i].position_y,
                all_atoms_[i].position_z,
                atoms.x_low,
                atoms.y_low,
                atoms.z_low) <= search_radius_ * search_radius_)
            {
                ghost_atoms_.push_back(all_atoms_[i]);
            }
        }
    }

    neighbor_list_.initialize(inside_atoms_.size(), std::max(1, static_cast<int>(all_atoms_.size()) * neighbor_reserve_factor));
}

void NeighborSearch::build_neighbors()
{
    bin_manager_.init_bins(search_radius_, inside_atoms_, ghost_atoms_);
    bin_manager_.do_binning(inside_atoms_, ghost_atoms_);
    bin_manager_.build_atom_neighbors(neighbor_list_, inside_atoms_);
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
    assert(mpi_size > 0);
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

void NeighborSearch::decompose(int mpi_size, double span_x, double span_y, double span_z, int& nx, int& ny, int& nz)
{
    assert(mpi_size > 0);

    nx = 1;
    ny = 1;
    nz = 1;

    span_x = std::max(0.0, span_x);
    span_y = std::max(0.0, span_y);
    span_z = std::max(0.0, span_z);

    const bool can_split_x = span_x > coord_tolerance;
    const bool can_split_y = span_y > coord_tolerance;
    const bool can_split_z = span_z > coord_tolerance;
    if (!can_split_x && !can_split_y && !can_split_z)
    {
        return;
    }

    std::vector<int> factors;
    int remaining = mpi_size;
    for (int factor = 2; factor * factor <= remaining; ++factor)
    {
        while (remaining % factor == 0)
        {
            factors.push_back(factor);
            remaining /= factor;
        }
    }
    if (remaining > 1)
    {
        factors.push_back(remaining);
    }
    std::sort(factors.rbegin(), factors.rend());

    for (const int factor : factors)
    {
        int* best_dim = nullptr;
        double best_score = -1.0;

        const auto try_dimension = [&](bool can_split, double span, int& dim) {
            if (!can_split)
            {
                return;
            }
            const double score = span / dim;
            if (score > best_score)
            {
                best_score = score;
                best_dim = &dim;
            }
        };

        try_dimension(can_split_x, span_x, nx);
        try_dimension(can_split_y, span_y, ny);
        try_dimension(can_split_z, span_z, nz);
        assert(best_dim != nullptr);
        *best_dim *= factor;
    }
}
