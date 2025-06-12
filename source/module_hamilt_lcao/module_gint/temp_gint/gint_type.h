#pragma once

#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include "source_base/vector3.h"
#include "source_base/matrix3.h"

namespace ModuleGint
{
    using Matrix3 = ModuleBase::Matrix3;
    using Vec3d = ModuleBase::Vector3<double>;
    using Vec3i = ModuleBase::Vector3<int>;

    template <typename T>
    using HContainer = hamilt::HContainer<T>;
}