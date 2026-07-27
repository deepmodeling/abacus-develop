#ifndef STERNHEIMER_FD_HAMILTONIAN_H
#define STERNHEIMER_FD_HAMILTONIAN_H

#include "source_lcao/module_ri/sternheimer_fd_nonlocal_projector.h"
#include "source_lcao/module_ri/sternheimer_kq.h"

#include <array>
#include <algorithm>
#include <cmath>
#include <complex>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <vector>

namespace ModuleRI
{

class SternheimerFDHamiltonian
{
  public:
    using Complex = std::complex<double>;
    using Vector = std::vector<Complex>;
    using Matrix = std::vector<Vector>;

    struct Eigenpairs
    {
        std::vector<double> eigenvalues;
        std::vector<Vector> eigenvectors;
    };

    struct Grid
    {
        int nx = 0;
        int ny = 0;
        int nz = 0;
        double hx = 0.0;
        double hy = 0.0;
        double hz = 0.0;
        bool periodic = true;
        SternheimerReducedKPoint kpoint{0.0, 0.0, 0.0};
        std::array<std::array<double, 3>, 3> lattice_vectors{};

        int size() const;
    };

    // kinetic_prefactor = 0.5 gives the Hartree convention -1/2 laplacian.
    // ABACUS PW eigenvalues are in Ry, so ABACUS bindings pass 1.0.
    SternheimerFDHamiltonian(Grid grid,
                             std::vector<double> local_potential,
                             double kinetic_prefactor = 0.5,
                             std::shared_ptr<const SternheimerFDNonlocalProjector> nonlocal_projector = nullptr);

    const Grid& grid() const;
    const std::vector<double>& local_potential() const;
    double kinetic_prefactor() const;
    const SternheimerReducedKPoint& kpoint() const;
    const SternheimerFDNonlocalProjector* nonlocal_projector() const;

    void apply(const Vector& psi, Vector& hpsi) const;
    void apply(const Vector& psi, Vector& hpsi, int* threads_used) const;

    Matrix dense_matrix(int max_size = 4096) const;

    // Diagonalizes the zero-order Hermitian FD Hamiltonian on small debug grids.
    // The shifted Sternheimer operator H - eps + i omega is non-Hermitian and
    // must be handled as a linear operator, not by this routine.
    Eigenpairs diagonalize_dense(int max_size = 4096) const;

    int index(int ix, int iy, int iz) const;

  private:
    struct ShiftedGridPoint
    {
        int index = -1;
        Complex phase{1.0, 0.0};
    };

    ShiftedGridPoint shifted_grid_point(int ix, int iy, int iz) const;

    Grid grid_;
    std::vector<double> local_potential_;
    double kinetic_prefactor_ = 0.5;
    std::shared_ptr<const SternheimerFDNonlocalProjector> nonlocal_projector_;
};

using SternheimerFDLatticeVectors = std::array<std::array<double, 3>, 3>;
using SternheimerFDReducedRotation = std::array<std::array<double, 3>, 3>;

inline SternheimerFDLatticeVectors sternheimer_fd_grid_lattice_vectors(
    const SternheimerFDHamiltonian::Grid& grid)
{
    bool has_explicit_lattice = false;
    for (const auto& vector: grid.lattice_vectors)
    {
        for (const double component: vector)
        {
            has_explicit_lattice = has_explicit_lattice || component != 0.0;
        }
    }
    if (has_explicit_lattice)
    {
        return grid.lattice_vectors;
    }
    return {{{grid.nx * grid.hx, 0.0, 0.0},
             {0.0, grid.ny * grid.hy, 0.0},
             {0.0, 0.0, grid.nz * grid.hz}}};
}

inline SternheimerFDLatticeVectors sternheimer_fd_grid_dual_vectors(
    const SternheimerFDHamiltonian::Grid& grid)
{
    const SternheimerFDLatticeVectors lattice = sternheimer_fd_grid_lattice_vectors(grid);
    const auto cross = [](const std::array<double, 3>& lhs, const std::array<double, 3>& rhs) {
        return std::array<double, 3>{lhs[1] * rhs[2] - lhs[2] * rhs[1],
                                     lhs[2] * rhs[0] - lhs[0] * rhs[2],
                                     lhs[0] * rhs[1] - lhs[1] * rhs[0]};
    };
    const auto dot = [](const std::array<double, 3>& lhs, const std::array<double, 3>& rhs) {
        return lhs[0] * rhs[0] + lhs[1] * rhs[1] + lhs[2] * rhs[2];
    };
    SternheimerFDLatticeVectors dual{{cross(lattice[1], lattice[2]),
                                      cross(lattice[2], lattice[0]),
                                      cross(lattice[0], lattice[1])}};
    const double determinant = dot(lattice[0], dual[0]);
    for (auto& vector: dual)
    {
        for (double& component: vector)
        {
            component /= determinant;
        }
    }
    return dual;
}

inline std::vector<int> sternheimer_fd_second_order_stencil_symmetry_indices(
    const SternheimerFDHamiltonian::Grid& grid,
    const std::vector<SternheimerFDReducedRotation>& rotations,
    const double tolerance = 1.0e-10)
{
    if (grid.nx <= 0 || grid.ny <= 0 || grid.nz <= 0 || tolerance <= 0.0)
    {
        throw std::invalid_argument("Invalid grid or tolerance for Sternheimer FD stencil symmetry.");
    }

    using Offset = std::array<int, 3>;
    std::map<Offset, double> stencil;
    const SternheimerFDLatticeVectors dual = sternheimer_fd_grid_dual_vectors(grid);
    const std::array<double, 3> dimensions{
        static_cast<double>(grid.nx), static_cast<double>(grid.ny), static_cast<double>(grid.nz)};
    std::array<std::array<double, 3>, 3> coefficients{};
    for (int left = 0; left != 3; ++left)
    {
        for (int right = 0; right != 3; ++right)
        {
            for (int component = 0; component != 3; ++component)
            {
                coefficients[left][right]
                    += dual[left][component] * dual[right][component]
                       * dimensions[left] * dimensions[right];
            }
        }
    }
    const auto add = [&stencil](const Offset& offset, const double weight) {
        stencil[offset] += weight;
    };
    for (int direction = 0; direction != 3; ++direction)
    {
        Offset positive{0, 0, 0};
        Offset negative{0, 0, 0};
        positive[direction] = 1;
        negative[direction] = -1;
        add(positive, coefficients[direction][direction]);
        add(negative, coefficients[direction][direction]);
    }
    for (int left = 0; left != 3; ++left)
    {
        for (int right = left + 1; right != 3; ++right)
        {
            for (const int left_sign: {-1, 1})
            {
                for (const int right_sign: {-1, 1})
                {
                    Offset offset{0, 0, 0};
                    offset[left] = left_sign;
                    offset[right] = right_sign;
                    add(offset,
                        0.5 * coefficients[left][right]
                            * static_cast<double>(left_sign * right_sign));
                }
            }
        }
    }
    double weight_scale = 0.0;
    for (const auto& entry: stencil)
    {
        const double weight = entry.second;
        weight_scale = std::max(weight_scale, std::abs(weight));
    }
    for (auto iter = stencil.begin(); iter != stencil.end();)
    {
        if (std::abs(iter->second) <= tolerance * std::max(1.0, weight_scale))
        {
            iter = stencil.erase(iter);
        }
        else
        {
            ++iter;
        }
    }

    std::vector<int> preserving;
    for (std::size_t isym = 0; isym != rotations.size(); ++isym)
    {
        bool valid = true;
        std::set<Offset> mapped_offsets;
        for (const auto& entry: stencil)
        {
            const Offset& offset = entry.first;
            const double weight = entry.second;
            Offset mapped{};
            for (int target = 0; target != 3; ++target)
            {
                double mapped_index = 0.0;
                for (int source = 0; source != 3; ++source)
                {
                    mapped_index += static_cast<double>(offset[source])
                                    / dimensions[source]
                                    * rotations[isym][source][target]
                                    * dimensions[target];
                }
                mapped[target] = static_cast<int>(std::llround(mapped_index));
                if (std::abs(mapped_index - static_cast<double>(mapped[target])) > tolerance)
                {
                    valid = false;
                    break;
                }
            }
            if (!valid)
            {
                break;
            }
            const auto mapped_iter = stencil.find(mapped);
            if (mapped_iter == stencil.end()
                || std::abs(mapped_iter->second - weight)
                       > tolerance
                             * std::max(std::max(1.0, weight_scale),
                                        std::max(std::abs(weight),
                                                 std::abs(mapped_iter->second)))
                || !mapped_offsets.insert(mapped).second)
            {
                valid = false;
                break;
            }
        }
        if (valid && mapped_offsets.size() == stencil.size())
        {
            preserving.push_back(static_cast<int>(isym));
        }
    }
    return preserving;
}

inline std::array<double, 3> sternheimer_fd_grid_cartesian_position(
    const SternheimerFDHamiltonian::Grid& grid,
    const int ix,
    const int iy,
    const int iz)
{
    const SternheimerFDLatticeVectors lattice = sternheimer_fd_grid_lattice_vectors(grid);
    const std::array<double, 3> reduced{static_cast<double>(ix) / grid.nx,
                                        static_cast<double>(iy) / grid.ny,
                                        static_cast<double>(iz) / grid.nz};
    std::array<double, 3> cartesian{};
    for (int direction = 0; direction != 3; ++direction)
    {
        for (int component = 0; component != 3; ++component)
        {
            cartesian[component] += reduced[direction] * lattice[direction][component];
        }
    }
    return cartesian;
}

inline std::array<double, 3> sternheimer_fd_grid_lattice_translation(
    const SternheimerFDHamiltonian::Grid& grid,
    const std::array<int, 3>& translation)
{
    const SternheimerFDLatticeVectors lattice = sternheimer_fd_grid_lattice_vectors(grid);
    std::array<double, 3> cartesian{};
    for (int direction = 0; direction != 3; ++direction)
    {
        for (int component = 0; component != 3; ++component)
        {
            cartesian[component] += static_cast<double>(translation[direction]) * lattice[direction][component];
        }
    }
    return cartesian;
}

} // namespace ModuleRI

#endif
