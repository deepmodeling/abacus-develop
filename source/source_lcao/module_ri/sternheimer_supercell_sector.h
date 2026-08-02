#ifndef STERNHEIMER_SUPERCELL_SECTOR_H
#define STERNHEIMER_SUPERCELL_SECTOR_H

#include "source_lcao/module_ri/sternheimer_kq.h"

#include <array>
#include <complex>
#include <vector>

namespace ModuleRI
{

struct SternheimerSupercellSector
{
    std::vector<double> eigenvalues;
    std::vector<std::vector<std::complex<double>>> coefficients;
    double max_orthonormality_error = 0.0;
    double max_full_space_residual = 0.0;
};

struct SternheimerSupercellKPointSector
{
    SternheimerReducedKPoint kpoint{0.0, 0.0, 0.0};
    SternheimerSupercellSector sector;
};

SternheimerSupercellSector recover_sternheimer_supercell_sector(
    const std::vector<double>& eigenvalues,
    const std::vector<std::vector<std::complex<double>>>& eigenvectors,
    const std::array<int, 3>& repeats,
    const SternheimerReducedKPoint& primitive_kpoint);

std::vector<SternheimerSupercellKPointSector> recover_all_sternheimer_supercell_sectors(
    const std::vector<double>& eigenvalues,
    const std::vector<std::vector<std::complex<double>>>& eigenvectors,
    const std::array<int, 3>& repeats);

} // namespace ModuleRI

#endif
