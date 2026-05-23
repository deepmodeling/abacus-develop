#pragma once

#include "../potential_base.h"

#include <memory>
#include <string>
#include <vector>
#include <array>

namespace ModuleMD
{

class NEPCalculator; // 你可以在 nep/nep3.h 中定义或适配这个类

class NEPPotential : public PotentialBase
{
public:
    NEPPotential() = default;
    ~NEPPotential() override = default;

    void init(const std::string& model_file) override;

    PotentialResult compute(
        const std::vector<std::array<double, 3>>& positions,
        const std::vector<int>& atom_types,
        const std::array<double, 9>& cell,
        bool pbc) override;

    std::string name() const override
    {
        return "NEP";
    }

private:
    std::unique_ptr<NEPCalculator> nep_model_;
    std::string model_file_;
};

} // namespace ModuleMD