#include "source_base/constants.h"
#include "source_lcao/module_rt/td_folding.h"

#include <gtest/gtest.h>

namespace
{
TEST(TdFoldingTest, BlochAndHybridPhases)
{
    using Complex = std::complex<double>;
    using Vector = ModuleBase::Vector3<double>;
    Parallel_Orbitals orbitals;
    orbitals.para_2d().set_serial(1, 1);
    const int atom_start = 0;
    orbitals.set_atomic_trace(&atom_start, 1, 1);
    const ModuleBase::Vector3<int> r(1, -1, 2);
    hamilt::HContainer<Complex> hr(&orbitals);
    hamilt::AtomPair<Complex> pair(
        0,
        0,
        r,
        &orbitals);
    hr.insert_pair(pair);
    hr.allocate(nullptr, true);
    const Complex value(1.0, 0.2);
    hr.get_atom_pair(0).get_pointer(0)[0] = value;
    const Vector k(0.25, -0.125, 0.375);
    const Vector direct(r.x, r.y, r.z);
    const std::vector<Vector> fields = {{0.0, 0.0, 0.0}, {0.07, -0.11, 0.05}};
    for (const auto &field : fields)
    {
        const double field_angle = field * direct;
        const Complex field_exponent(0.0, field_angle);
        const std::map<ModuleBase::Vector3<int>, Complex> phases = {{r, std::exp(field_exponent)}};
        Complex folded(0.0, 0.0);
        module_rt::folding_HR_td(
            hr,
            &folded,
            k,
            field,
            phases,
            1,
            1);
        const double angle = ModuleBase::TWO_PI * (k * direct) + field_angle;
        const Complex exponent(0.0, angle);
        const Complex expected = value * std::exp(exponent);
        const double error = std::abs(folded - expected);
        EXPECT_LT(error, 1.0e-12);
    }
}
} // namespace
