#ifndef STERNHEIMER_SUPERCELL_PERTURBATION_H
#define STERNHEIMER_SUPERCELL_PERTURBATION_H

#include "source_lcao/module_ri/sternheimer_abfs_perturbation.h"

#include <array>
#include <string>
#include <vector>

namespace ModuleRI
{

struct SternheimerSupercellTranslationSum
{
    std::array<int, 3> repeats{1, 1, 1};
    SternheimerReducedKPoint primitive_qpoint{0.0, 0.0, 0.0};
    int atoms_per_primitive = 0;
    int basis_atom = -1;
    int channel_within_atom = -1;
};

SternheimerSupercellTranslationSum parse_sternheimer_supercell_translation_sum(
    const std::string& specification);

SternheimerABFBlochGridChannel combine_sternheimer_supercell_translation_channel(
    const std::vector<SternheimerABFBlochGridChannel>& channels,
    const SternheimerSupercellTranslationSum& config);

std::vector<SternheimerABFBlochGridChannel>
combine_all_sternheimer_supercell_translation_channels(
    const std::vector<SternheimerABFBlochGridChannel>& channels,
    const SternheimerSupercellTranslationSum& config);

int sternheimer_supercell_primitive_cell_count(const SternheimerSupercellTranslationSum& config);

} // namespace ModuleRI

#endif
