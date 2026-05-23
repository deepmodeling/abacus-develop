#include "dpmd_potential.h"

#include <stdexcept>

namespace ModuleMD
{

void DPMDPotential::init(const std::string& model_file)
{
    model_file_ = model_file;

#ifndef USE_DEEPMD
    throw std::runtime_error(
        "DPMDPotential requires deepmd-kit, but USE_DEEPMD is not enabled.");
#else
    dp_model_ = std::make_unique<deepmd::DeepPot>();
    dp_model_->init(model_file);
#endif
}

PotentialResult DPMDPotential::compute(
    const std::vector<std::array<double, 3>>& positions,
    const std::vector<int>& atom_types,
    const std::array<double, 9>& cell,
    bool pbc)
{
#ifndef USE_DEEPMD
    throw std::runtime_error(
        "DPMDPotential::compute called without deepmd-kit support.");
#else
    const int natom = static_cast<int>(positions.size());

    if (static_cast<int>(atom_types.size()) != natom)
    {
        throw std::runtime_error("DPMDPotential: atom_types size mismatch.");
    }

    std::vector<double> coord;
    coord.reserve(3 * natom);

    for (int i = 0; i < natom; ++i)
    {
        coord.push_back(positions[i][0]);
        coord.push_back(positions[i][1]);
        coord.push_back(positions[i][2]);
    }

    std::vector<double> box(cell.begin(), cell.end());

    double energy = 0.0;
    std::vector<double> force(3 * natom, 0.0);
    std::vector<double> virial(9, 0.0);

    // deepmd-kit DeepPot C++ API:
    // compute energy, force, virial from coord, atom types, and box.
    dp_model_->compute(
        energy,
        force,
        virial,
        coord,
        atom_types,
        box);

    PotentialResult result;
    result.energy = energy;
    result.force.resize(natom);

    for (int i = 0; i < natom; ++i)
    {
        result.force[i][0] = force[3 * i + 0];
        result.force[i][1] = force[3 * i + 1];
        result.force[i][2] = force[3 * i + 2];
    }

    for (int i = 0; i < 9; ++i)
    {
        result.virial[i] = virial[i];
    }

    return result;
#endif
}

} // namespace ModuleMD