#pragma once

#include "potential_base.h"

#include <memory>
#include <string>

namespace ModuleMD
{

class PotentialFactory
{
public:
    static std::unique_ptr<PotentialBase> create(
        const std::string& potential_type,
        const std::string& model_file);
};

} // namespace ModuleMD