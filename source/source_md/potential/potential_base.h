#pragma once

#include <string>
#include <vector>
#include <array>

namespace ModuleMD
{

struct PotentialResult
{
    double energy = 0.0;

    // force[i][0], force[i][1], force[i][2]
    std::vector<std::array<double, 3>> force;

    // virial: 3x3, row-major
    std::array<double, 9> virial = {0.0};
};

class PotentialBase
{
public:
    virtual ~PotentialBase() = default;

    virtual void init(const std::string& model_file) = 0;

    virtual PotentialResult compute(
        const std::vector<std::array<double, 3>>& positions,
        const std::vector<int>& atom_types,
        const std::array<double, 9>& cell,
        bool pbc) = 0;

    virtual std::string name() const = 0;
};

} // namespace ModuleMD