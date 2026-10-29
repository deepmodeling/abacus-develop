#include "module_hamilt_general/operator.h"
#include <mpi.h>
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <vector>

using Complex = std::complex<double>;
using CPU = base_device::DEVICE_CPU;
using Operator = hamilt::Operator<Complex, CPU>;

class ScaleOperator : public Operator
{
public:
    explicit ScaleOperator(double scale, hamilt::calculation_type type) : scale_(scale)
    {
        this->cal_type = type;
    }
    void act(int columns, int leading, int polarizations, const Complex* input,
             Complex* output, int) const override
    {
        for (size_t i = 0; i < size_t(columns) * leading / polarizations; ++i) { output[i] += scale_ * input[i]; }
    }
private:
    double scale_;
};

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    double maximum = 0;
    for (int polarizations : {1, 2})
    {
        const int bands = 3, leading = 7 * polarizations;
        int plane_waves[] = {7};
        std::vector<Complex> input(bands * leading), output(input.size());
        for (size_t i = 0; i < input.size(); ++i) { input[i] = Complex(0.25 * i, -0.125 * i); }
        psi::Psi<Complex, CPU> wavefunctions(input.data(), 1, bands, leading, plane_waves);
        wavefunctions.npol = polarizations;
        std::unique_ptr<Operator> head(new ScaleOperator(2, hamilt::calculation_type::pw_ekinetic));
        auto* second = new ScaleOperator(3, hamilt::calculation_type::pw_nonlocal);
        head->add(second);
        head->init(0);
        int handled = 0;
        Operator::Action replacement = [&](const Operator* op, int columns, int stride, int npol,
                                            const Complex* values, Complex* target, int npw) {
            if (op != second) { return false; }
            if (columns != bands * polarizations || stride != leading || npol != polarizations || npw != 7)
            {
                std::cerr << "Invalid operator-action dimensions\n";
                std::exit(2);
            }
            ++handled;
            for (size_t i = 0; i < input.size(); ++i) { target[i] += 7.0 * values[i]; }
            return true;
        };
        for (int call = 0; call < 3; ++call)
        {
            Operator::hpsi_info info(&wavefunctions, psi::Range(true, 0, 0, bands - 1), output.data());
            head->hPsi(info, call == 1 ? replacement : Operator::Action{});
            const double scale = call == 1 ? 9 : 5;
            for (size_t i = 0; i < input.size(); ++i)
            {
                const double difference = std::abs(output[i] - scale * input[i]);
                if (!std::isfinite(difference)) { return 4; }
                maximum = std::max(maximum, difference);
            }
        }
        if (handled != 1) { return 3; }
    }
    MPI_Finalize();
    std::cout << "operator_action_error=" << maximum << std::endl;
    return maximum == 0 ? 0 : 4;
}
