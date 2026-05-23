#include "nep_potential.h"

// 这里换成你实际引入的 NEP 头文件
#include "nep/nep_cpu.h"

#include <stdexcept>

namespace ModuleMD
{

void NEPPotential::init(const std::string& model_file)
{
    model_file_ = model_file;

    nep_model_ = std::make_unique<NEPCalculator>();
    nep_model_->load_model(model_file);
}

PotentialResult NEPPotential::compute(
    const std::vector<std::array<double, 3>>& positions,
    const std::vector<int>& atom_types,
    const std::array<double, 9>& cell,
    bool pbc)
{
    const int natom = static_cast<int>(positions.size());

    if (static_cast<int>(atom_types.size()) != natom)
    {
        throw std::runtime_error("NEPPotential: atom_types size mismatch.");
    }

    PotentialResult result;
    result.force.resize(natom);

    // 下面是假设接口，你需要根据实际 NEP 源码改名。
    nep_model_->compute(
        natom,
        positions,
        atom_types,
        cell,
        pbc,
        result.energy,
        result.force,
        result.virial);

    return result;
}

} // namespace ModuleMD