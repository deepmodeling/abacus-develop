#include "source_pw/module_rt/unittests/td_test_utils.h"

#include <cmath>

namespace td_test
{
TEST_F(Observables, CurrentsMatchPauliAndHamiltonianDerivative)
{
    const ModuleBase::Vector3<double> a(0.13, -0.07, 0.09);
    pw::TDVelPW<double, Device> observer;
    observer.initialize(cell_, basis_, kv_, pp_);
    pw::CurrentRequest request;
    request.charge_total = true;
    request.charge_per_k = true;
    request.spin_current = true;
    const pw::TDVelResultPW actual = observer.calculate(*psi_, state_, a, request);
    hamilt::Velocity<double, Device> velocity(&basis_, kv_.isk.data(), &pp_, &cell_, true, nullptr, 0, 0);
    velocity.init(0, a);
    const int n = psi_->get_nbasis();
    std::array<double, 3> charge{};
    std::array<double, 9> spin{};
    const double epsilon = 1e-4;
    for (int band = 0; band < psi_->get_nbands(); ++band)
    {
        const Complex* data = psi_->get_pointer() + band * n;
        const Vector x(data, data + n);
        Vector vx(3 * n);
        velocity.act(psi_.get(), 2, x.data(), vx.data(), false);
        for (int direction = 0; direction < 3; ++direction)
        {
            const Vector component(vx.begin() + direction * n, vx.begin() + (direction + 1) * n);
            ModuleBase::Vector3<double> plus = a;
            ModuleBase::Vector3<double> minus = a;
            plus[direction] += epsilon;
            minus[direction] -= epsilon;
            const Vector hp = hamiltonian(plus, x);
            const Vector hm = hamiltonian(minus, x);
            double error = 0;
            double norm = 0;
            for (int i = 0; i < n; ++i)
            {
                const Complex derivative = (hp[i] - hm[i]) / (2 * epsilon);
                error += std::norm(component[i] - derivative);
                norm += std::norm(derivative);
            }
            EXPECT_LE(std::sqrt(error), 1e-10 + 1e-6 * std::sqrt(norm));
            charge[direction] -= state_.wg(0, band) * inner(x, component).real() / cell_.omega;
            for (int axis = 0; axis < 3; ++axis)
            {
                const Vector sx = pauli(x, axis);
                Vector vsx(3 * n);
                velocity.act(psi_.get(), 2, sx.data(), vsx.data(), false);
                const Vector v_sigma(vsx.begin() + direction * n, vsx.begin() + (direction + 1) * n);
                const Vector sigma_v = pauli(component, axis);
                const Complex anticommutator = inner(x, v_sigma) + inner(x, sigma_v);
                spin[3 * direction + axis] += 0.25 * state_.wg(0, band) * anticommutator.real() / cell_.omega;
            }
        }
    }
    double signal = 0;
    for (int i = 0; i < 3; ++i)
    {
        EXPECT_NEAR(actual.total[i], charge[i], 1e-10);
        EXPECT_NEAR(actual.per_k[i], charge[i], 1e-10);
    }
    for (int i = 0; i < 9; ++i)
    {
        EXPECT_NEAR(actual.spin[i], spin[i], 1e-10);
        signal += std::abs(spin[i]);
    }
    EXPECT_GT(signal, 1e-8);
}
} // namespace td_test
