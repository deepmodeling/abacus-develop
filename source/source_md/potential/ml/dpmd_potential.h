#pragma once

#include "../potential_base.h"

#ifdef USE_DEEPMD
#include "deepmd/DeepPot.h"
#endif

#include <memory>
#include <string>
#include <vector>
#include <array>

namespace ModuleMD
{

class DPMDPotential : public PotentialBase
{
public:
    DPMDPotential() = default;
    ~DPMDPotential() override = default;

    void init(const std::string& model_file) override;

    PotentialResult compute(
        const std::vector<std::array<double, 3>>& positions,
        const std::vector<int>& atom_types,
        const std::array<double, 9>& cell,
        bool pbc) override;

    std::string name() const override
    {
        return "DPMD";
    }

private:
#ifdef USE_DEEPMD
    std::unique_ptr<deepmd::DeepPot> dp_model_;
#endif

    std::string model_file_;
};

} // namespace ModuleMD