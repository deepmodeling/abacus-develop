#include "source_dftb/sk_matrix.h"

#include <cmath>
#include <stdexcept>
#include <string>

namespace ModuleDFTB
{
namespace
{
std::size_t matrix_index(std::size_t dimension, std::size_t row, std::size_t col)
{
    if (row >= dimension || col >= dimension)
    {
        throw std::out_of_range("DFTB Bloch matrix index is out of range");
    }
    return row * dimension + col;
}

bool canonical_pair(const DftbPairImage& pair)
{
    if (pair.atom_a < pair.atom_b) return true;
    if (pair.atom_a > pair.atom_b) return false;
    for (const double component : pair.translation_bohr)
    {
        if (component > 1.0e-12) return true;
        if (component < -1.0e-12) return false;
    }
    return false;
}

const DftbPairParameters& find_parameters(std::size_t species_a,
                                          std::size_t species_b,
                                          const std::vector<DftbPairParameters>& parameters)
{
    for (const auto& pair : parameters)
    {
        if (pair.species_a == species_a && pair.species_b == species_b)
        {
            if (pair.ab == nullptr || pair.ba == nullptr)
            {
                throw std::runtime_error("DFTB pair parameters contain a null directed SKF file");
            }
            return pair;
        }
    }
    throw std::runtime_error("Missing directed SKF pair for species indices "
                             + std::to_string(species_a) + "-" + std::to_string(species_b));
}
} // namespace

std::complex<double>& DftbBlochMatrices::h(std::size_t row, std::size_t col)
{
    return hamiltonian[matrix_index(dimension, row, col)];
}

std::complex<double>& DftbBlochMatrices::s(std::size_t row, std::size_t col)
{
    return overlap[matrix_index(dimension, row, col)];
}

const std::complex<double>& DftbBlochMatrices::h(std::size_t row, std::size_t col) const
{
    return hamiltonian[matrix_index(dimension, row, col)];
}

const std::complex<double>& DftbBlochMatrices::s(std::size_t row, std::size_t col) const
{
    return overlap[matrix_index(dimension, row, col)];
}

double calculate_repulsive_energy_hartree(const std::vector<DftbSpAtom>& atoms,
                                           const std::vector<DftbPairImage>& pair_images,
                                           const std::vector<DftbPairParameters>& parameters)
{
    double energy = 0.0;
    for (const auto& pair : pair_images)
    {
        if (pair.atom_a >= atoms.size() || pair.atom_b >= atoms.size())
        {
            throw std::out_of_range("DFTB repulsive pair references an invalid atom index");
        }
        if (!canonical_pair(pair))
        {
            throw std::invalid_argument("DFTB repulsive pair is not in the canonical half-list orientation");
        }
        const auto& atom_a = atoms[pair.atom_a];
        const auto& atom_b = atoms[pair.atom_b];
        const auto& sk = find_parameters(atom_a.species, atom_b.species, parameters);
        if (!sk.ab->has_repulsive_spline)
        {
            throw std::runtime_error("Native repulsive energy currently requires an SKF spline representation");
        }
        double distance_squared = 0.0;
        for (std::size_t axis = 0; axis < 3; ++axis)
        {
            const double displacement = atom_b.position_bohr[axis] + pair.translation_bohr[axis]
                                       - atom_a.position_bohr[axis];
            distance_squared += displacement * displacement;
        }
        const double distance = std::sqrt(distance_squared);
        if (!(distance > 0.0) || !std::isfinite(distance))
        {
            throw std::invalid_argument("DFTB repulsive pair must have finite, positive length");
        }
        energy += sk.ab->repulsive.evaluate(distance).energy_hartree;
    }
    return energy;
}

DftbBlochMatrices assemble_sp_bloch_matrices(const std::vector<DftbSpAtom>& atoms,
                                               const std::vector<DftbPairImage>& pair_images,
                                               const std::vector<DftbPairParameters>& parameters,
                                               const std::array<double, 3>& k_bohr_inverse)
{
    constexpr std::size_t orbitals_per_atom = 4;
    DftbBlochMatrices matrices;
    matrices.dimension = orbitals_per_atom * atoms.size();
    matrices.hamiltonian.resize(matrices.dimension * matrices.dimension);
    matrices.overlap.resize(matrices.dimension * matrices.dimension);

    for (std::size_t atom_index = 0; atom_index < atoms.size(); ++atom_index)
    {
        const auto& atom = atoms[atom_index];
        if (atom.homonuclear_data == nullptr || !atom.homonuclear_data->has_atomic_data)
        {
            throw std::runtime_error("Each DFTB s+p atom requires homonuclear SKF atomic data");
        }
        const std::size_t base = orbitals_per_atom * atom_index;
        matrices.h(base, base) = atom.homonuclear_data->onsite_hartree[0];
        matrices.s(base, base) = 1.0;
        for (std::size_t p = 1; p < orbitals_per_atom; ++p)
        {
            matrices.h(base + p, base + p) = atom.homonuclear_data->onsite_hartree[1];
            matrices.s(base + p, base + p) = 1.0;
        }
    }

    for (const auto& sk : parameters)
    {
        if (sk.ab == nullptr || sk.ba == nullptr)
        {
            throw std::runtime_error("DFTB pair parameters contain a null directed SKF file");
        }
        if (sk.ab->homonuclear != sk.ba->homonuclear)
        {
            throw std::runtime_error("DFTB pair parameters mix homonuclear and heteronuclear SKF files");
        }
        if (sk.ab->homonuclear)
        {
            if (sk.ab->filename != sk.ba->filename)
            {
                throw std::runtime_error("Homonuclear pair parameters must reference one SKF file");
            }
        }
        else
        {
            validate_pair_directions(*sk.ab, *sk.ba);
        }
    }

    for (const auto& pair : pair_images)
    {
        if (pair.atom_a >= atoms.size() || pair.atom_b >= atoms.size())
        {
            throw std::out_of_range("DFTB pair image references an invalid atom index");
        }
        if (!canonical_pair(pair))
        {
            throw std::invalid_argument("DFTB pair image is not in the canonical half-list orientation");
        }
        const auto& atom_a = atoms[pair.atom_a];
        const auto& atom_b = atoms[pair.atom_b];
        const auto& sk = find_parameters(atom_a.species, atom_b.species, parameters);
        const std::array<double, 3> displacement{{
            atom_b.position_bohr[0] + pair.translation_bohr[0] - atom_a.position_bohr[0],
            atom_b.position_bohr[1] + pair.translation_bohr[1] - atom_a.position_bohr[1],
            atom_b.position_bohr[2] + pair.translation_bohr[2] - atom_a.position_bohr[2]}};
        const auto block = evaluate_sp_pair_block(*sk.ab, *sk.ba, displacement);
        double phase_angle = 0.0;
        for (std::size_t axis = 0; axis < 3; ++axis)
        {
            phase_angle += k_bohr_inverse[axis] * pair.translation_bohr[axis];
        }
        const std::complex<double> phase(std::cos(phase_angle), std::sin(phase_angle));
        const std::size_t base_a = orbitals_per_atom * pair.atom_a;
        const std::size_t base_b = orbitals_per_atom * pair.atom_b;
        for (std::size_t row_b = 0; row_b < orbitals_per_atom; ++row_b)
        {
            for (std::size_t col_a = 0; col_a < orbitals_per_atom; ++col_a)
            {
                const std::size_t block_index = row_b * orbitals_per_atom + col_a;
                const auto h = block.hamiltonian[block_index] * phase;
                const auto s = block.overlap[block_index] * phase;
                matrices.h(base_b + row_b, base_a + col_a) += h;
                matrices.s(base_b + row_b, base_a + col_a) += s;
                matrices.h(base_a + col_a, base_b + row_b) += std::conj(h);
                matrices.s(base_a + col_a, base_b + row_b) += std::conj(s);
            }
        }
    }
    return matrices;
}

} // namespace ModuleDFTB
