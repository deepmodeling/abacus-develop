#include "source_lcao/module_deltaspin/lambda_loop_helper.h"
#include "source_lcao/module_deltaspin/spin_constrain.h"

#include "gtest/gtest.h"

namespace spinconstrain
{

TEST(LambdaLoopHelper, GradientDecayUsesAtomsOfEachType)
{
    SpinConstrain<std::complex<double>>& sc = SpinConstrain<std::complex<double>>::getScInstance();
    const std::map<int, int> atom_counts = {{0, 1}, {1, 2}};
    const std::vector<ModuleBase::Vector3<double>> spin(3, {0.0, 0.0, 0.0});
    const std::vector<ModuleBase::Vector3<double>> new_spin
        = {{0.0, 0.0, 10.0}, {0.0, 0.0, 0.1}, {0.0, 0.0, 0.2}};
    const std::vector<ModuleBase::Vector3<double>> old_lambda(3, {0.0, 0.0, 0.0});
    const std::vector<ModuleBase::Vector3<double>> new_lambda(3, {0.0, 0.0, 1.0});
    const std::vector<ModuleBase::Vector3<int>> constrain(3, {0, 0, 1});
    const std::vector<double> decay_grad = {0.0, 0.9};

    sc.set_atomCounts(atom_counts);
    sc.set_constrain(constrain.data(), 3);
    sc.set_decay_grad(decay_grad.data(), 2);

    EXPECT_TRUE(check_gradient_decay(sc, new_spin, spin, new_lambda, old_lambda, false, std::cout));
}

} // namespace spinconstrain
