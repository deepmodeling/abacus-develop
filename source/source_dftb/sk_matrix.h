#ifndef SOURCE_DFTB_SK_MATRIX_H
#define SOURCE_DFTB_SK_MATRIX_H

#include "source_dftb/sk_block.h"

#include <array>
#include <complex>
#include <cstddef>
#include <vector>

namespace ModuleDFTB
{

struct DftbSpAtom
{
    std::size_t species = 0;
    std::array<double, 3> position_bohr{};
    const SkfData* homonuclear_data = nullptr;
};

struct DftbPairParameters
{
    std::size_t species_a = 0;
    std::size_t species_b = 0;
    const SkfData* ab = nullptr;
    const SkfData* ba = nullptr;
};

/** One half-list pair: atom B is shifted from its cell by translation_bohr. */
struct DftbPairImage
{
    std::size_t atom_a = 0;
    std::size_t atom_b = 0;
    std::array<double, 3> translation_bohr{};
};

struct DftbBlochMatrices
{
    std::size_t dimension = 0;
    std::vector<std::complex<double>> hamiltonian;
    std::vector<std::complex<double>> overlap;

    std::complex<double>& h(std::size_t row, std::size_t col);
    std::complex<double>& s(std::size_t row, std::size_t col);
    const std::complex<double>& h(std::size_t row, std::size_t col) const;
    const std::complex<double>& s(std::size_t row, std::size_t col) const;
};

/**
 * Assemble dense H(k), S(k) for one s+p shell per atom. Pair images must be a
 * canonical half list: atom_a < atom_b, or for a self-image the first nonzero
 * translation component is positive. The routine inserts the Hermitian
 * counterpart with the conjugate Bloch phase. Positions, translations, and k
 * are in Bohr, Bohr, and inverse Bohr, respectively.
 */
double calculate_repulsive_energy_hartree(const std::vector<DftbSpAtom>& atoms,
                                           const std::vector<DftbPairImage>& pair_images,
                                           const std::vector<DftbPairParameters>& parameters);

DftbBlochMatrices assemble_sp_bloch_matrices(const std::vector<DftbSpAtom>& atoms,
                                               const std::vector<DftbPairImage>& pair_images,
                                               const std::vector<DftbPairParameters>& parameters,
                                               const std::array<double, 3>& k_bohr_inverse);

} // namespace ModuleDFTB

#endif // SOURCE_DFTB_SK_MATRIX_H
