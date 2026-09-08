#include "../symm_rotation.h"
#include "source_io/module_parameter/parameter.h"
#include "gtest/gtest.h"

class TestParameters
{
  public:
    TestParameters(Parameter& parameters, const int nspin)
        : parameters_(parameters), original_nspin_(parameters.inp.nspin)
    {
        parameters_.input.nspin = nspin;
    }
    ~TestParameters() { parameters_.input.nspin = original_nspin_; }

  private:
    Parameter& parameters_;
    const int original_nspin_;
};

// K-point generation is outside this test: use explicit stars, but provide
// the virtual symbols needed by the existing lightweight rotation test target.
void ModuleCell::ReciprocalGrid::renew(const int&)
{
    ADD_FAILURE() << "Unexpected k-point generation";
}
void K_Vectors::renew(const int&)
{
    ADD_FAILURE() << "Unexpected k-point generation";
}
void K_Vectors::reduce_by_symmetry(const UnitCell&, const ModuleSymmetry::Symmetry&,
                                  bool, std::string&, bool&, const int, std::ofstream&)
{
    ADD_FAILURE() << "Unexpected k-point reduction";
}

namespace
{
using Complex = std::complex<double>;

// Independent dense product in the stored (transposed density) convention.
std::vector<Complex> rotate_reference(const std::vector<Complex>& density,
                                      const std::vector<Complex>& rotation,
                                      const int n)
{
    std::vector<Complex> result(n * n, 0.0);
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            for (int a = 0; a < n; ++a)
            {
                for (int b = 0; b < n; ++b)
                {
                    result[i + j * n] += rotation[a + i * n] * density[a + b * n]
                                         * std::conj(rotation[b + j * n]);
                }
            }
        }
    }
    return result;
}

void check_little_group_restoration(const int nspin)
{
    const TestParameters parameters(PARAM, nspin);
    const int channels = nspin == 2 ? 2 : 1;
    const int n = 4;
    Parallel_2D pv;
    pv.init(n, n, 1, MPI_COMM_WORLD);
    ModuleSymmetry::Symmetry_rotation rotation;
    std::vector<Complex> identity(n * n, 0.0);
    std::vector<Complex> little(n * n, 0.0);
    std::vector<Complex> representative(n * n, 0.0);
    std::vector<Complex> alternate(n * n, 0.0);
    const int sign[n] = {1, 1, -1, -1};
    for (int i = 0; i < n; ++i)
    {
        identity[i + i * n] = 1.0;
        little[i + i * n] = sign[i];
        const int row = (i + 1) % n;
        representative[row + i * n] = std::polar(1.0, 0.3 * i);
        alternate[row + i * n] = double(sign[row]) * representative[row + i * n];
    }
    auto local = [&pv, n](const std::vector<Complex>& dense) {
        std::vector<Complex> result(pv.get_local_size());
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < n; ++j)
            {
                if (pv.in_this_processor(i, j))
                {
                    result[pv.global2local_row(i) + pv.global2local_col(j) * pv.get_row_size()]
                        = dense[i + j * n];
                }
            }
        }
        return result;
    };
    rotation.set_density_rotations_for_testing(
        {{{0, local(identity)}, {1, local(little)}, {2, local(representative)}, {3, local(alternate)}}},
        {{0, 1}}, 4);
    K_Vectors kv;
    kv.set_nkstot(channels);
    kv.set_nkstot_nospin(2);
    kv.kstars = {{{0, {0.25, 0.0, 0.0}}, {2, {0.0, 0.25, 0.0}}}};
    std::vector<std::vector<Complex>> inputs;
    std::vector<std::vector<Complex>> expected;
    for (int spin = 0; spin < channels; ++spin)
    {
        std::vector<Complex> density(n * n);
        std::vector<Complex> projected(n * n);
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < n; ++j)
            {
                const Complex value((spin + 1) * (2.0 + i + j), 0.2 * (i - j));
                density[i + j * n] = value;
                // For this C2 little group, averaging removes exactly the odd blocks.
                projected[i + j * n] = sign[i] == sign[j] ? 0.5 * value : Complex(0.0);
            }
        }
        inputs.push_back(local(density));
        expected.push_back(local(projected));
        expected.push_back(local(rotate_reference(projected, representative, n)));
    }
    const auto restored = rotation.restore_dm(kv, inputs, pv);
    ASSERT_EQ(restored.size(), expected.size());
    for (size_t k = 0; k < expected.size(); ++k)
    {
        for (size_t i = 0; i < expected[k].size(); ++i)
        {
            EXPECT_NEAR(std::abs(restored[k][i] - expected[k][i]), 0.0, 1e-12);
        }
    }
    // A different representative of the same star must give the same density.
    kv.kstars = {{{0, {0.25, 0.0, 0.0}}, {3, {0.0, 0.25, 0.0}}}};
    const auto changed_representative = rotation.restore_dm(kv, inputs, pv);
    for (size_t k = 0; k < expected.size(); ++k)
    {
        for (size_t i = 0; i < expected[k].size(); ++i)
        {
            EXPECT_NEAR(std::abs(changed_representative[k][i] - restored[k][i]), 0.0, 1e-12);
        }
    }
}
} // namespace

TEST(SymmetryDensityRestoration, LittleGroupAndStarWeight)
{
    check_little_group_restoration(1);
}

TEST(SymmetryDensityRestoration, IndependentSpinChannels)
{
    check_little_group_restoration(2);
}

TEST(SymmetryDensityRestoration, SpinorDensity)
{
    check_little_group_restoration(4);
}
