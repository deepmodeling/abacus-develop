#include "module_hamilt_pw/hamilt_pwdft/kernels/projector_fold.h"
#include "module_basis/module_pw/pw_basis_k.h"
#include "module_hsolver/kernels/math_kernel_op.h"
#include <mpi.h>
#include <cuda_runtime.h>
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <random>
#include <vector>

using Complex = std::complex<double>;
using GPU = base_device::DEVICE_GPU;
using CPU = base_device::DEVICE_CPU;

void require(bool condition, const char* message)
{
    if (!condition) { std::cerr << message << '\n'; std::exit(1); }
}
void cuda_check(cudaError_t status) { require(status == cudaSuccess, cudaGetErrorString(status)); }

class Buffer
{
public:
    Complex* data = nullptr;
    size_t size;
    explicit Buffer(size_t count) : size(count) { cuda_check(cudaMalloc(&data, size * sizeof(Complex))); }
    ~Buffer() { cuda_check(cudaFree(data)); }
    void upload(const std::vector<Complex>& values)
    {
        require(values.size() == size, "Buffer upload shape mismatch");
        cuda_check(cudaMemcpy(data, values.data(), size * sizeof(Complex), cudaMemcpyHostToDevice));
    }
    std::vector<Complex> download() const
    {
        std::vector<Complex> result(size);
        cuda_check(cudaMemcpy(result.data(), data, size * sizeof(Complex), cudaMemcpyDeviceToHost));
        return result;
    }
};

void run(int scenario, double& projection_error, double& addition_error)
{
    ModulePW::PW_Basis_K basis("gpu", "double");
    basis.initmpi(1, 0, MPI_COMM_WORLD);
    const ModuleBase::Matrix3 lattice(1, 0.2, 0.05, 0.1, 1.1, 0.15, 0, 0.05, 0.9);
    basis.initgrids(6.0, lattice, 28, 24, 30);
    ModuleBase::Vector3<double> kpoints[2] = {{0.25, -0.125, 0.375}, {-0.2, 0.15, -0.1}};
    basis.initparameters(false, 50.0, 2, kpoints, 1, true);
    basis.setuptransform();
    basis.collect_local_pw();
    const int atoms = 8, channels = 5, projectors = atoms * channels;
    std::vector<ModuleBase::Vector3<double>> direct(atoms), cartesian(atoms);
    for (int atom = 0; atom < atoms; ++atom)
    {
        direct[atom] = {0.5 * (atom / 4), 0.5 * ((atom / 2) % 2), 0.5 * (atom % 2)};
    }
    if (scenario == 1 || scenario == 2) { direct[0].x += 0.01234567; }
    if (scenario == 2) { direct[3].y += 0.02345671; }
    if (scenario == 3)
    {
        for (int atom = 0; atom < atoms; ++atom) { direct[atom].z += 0.00314159 * (atom + 1); }
    }
    if (scenario == 4)
    {
        const int coordinates[8][3] = {{0,0,0}, {1,3,5}, {2,6,1}, {3,1,4},
                                       {4,5,2}, {5,2,7}, {6,7,3}, {7,4,6}};
        for (int atom = 0; atom < atoms; ++atom)
        {
            direct[atom] = {coordinates[atom][0] / 8.0, coordinates[atom][1] / 8.0, coordinates[atom][2] / 8.0};
        }
    }
    for (int atom = 0; atom < atoms; ++atom)
    {
        const auto& r = direct[atom];
        cartesian[atom] = {r.x * lattice.e11 + r.y * lattice.e21 + r.z * lattice.e31,
                           r.x * lattice.e12 + r.y * lattice.e22 + r.z * lattice.e32,
                           r.x * lattice.e13 + r.y * lattice.e23 + r.z * lattice.e33};
    }
    auto folded = hamilt::ProjectorFold::create(atoms, channels, direct.data(), cartesian.data(), basis);
    if (scenario == 3)
    {
        require(folded == nullptr, "Irregular geometry unexpectedly selected folding");
        std::cout << "geometry_fallback=pass\n";
        return;
    }
    require(folded != nullptr, "Regular geometry did not select folding");
    const int projector_leading = basis.npwk_max + 3, leading = basis.npwk_max + 7;
    const Complex padding(17.25, -23.5), one(1, 0), zero(0, 0);
    std::mt19937 generator(3107 + scenario);
    std::uniform_real_distribution<double> uniform(-0.5, 0.5);
    for (int ik : {0, 1, 0})
    {
        const int rows = basis.npwk[ik];
        std::vector<Complex> matrix(size_t(projector_leading) * projectors, padding);
        for (int atom = 0; atom < atoms; ++atom)
        for (int channel = 0; channel < channels; ++channel)
        for (int g = 0; g < rows; ++g)
        {
            const auto momentum = basis.getgdirect(ik, g) + kpoints[ik];
            const double phase = -2 * std::acos(-1.0) * (momentum * direct[atom]);
            const Complex radial(std::sin(0.13 * (channel + 1) * (g + 1)), std::cos(0.17 * (channel + 2) * (g + 1)));
            matrix[size_t(atom * channels + channel) * projector_leading + g]
                = radial * std::polar(1.0, phase) / std::sqrt(double(rows));
        }
        Buffer device_matrix(matrix.size());
        device_matrix.upload(matrix);
        folded->prepare(ik, device_matrix.data, projector_leading);
        if (scenario == 4)
        {
            std::vector<Complex> preserved(size_t(projectors) * 2, padding);
            std::vector<Complex> untouched_input(size_t(leading) * 2, padding);
            Buffer input(untouched_input.size()), output(preserved.size());
            input.upload(untouched_input);
            output.upload(preserved);
            require(!folded->project(2, leading, input.data, output.data), "Uneconomical residue grid selected folding");
            require(!folded->add(2, leading, output.data, input.data), "Uneconomical residue grid selected expansion");
            require(output.download() == preserved && input.download() == untouched_input, "Fallback modified an output");
            continue;
        }
        for (int fields : {1, 2, 31, 127, 128, 129, 259, 7})
        {
            std::vector<Complex> input(size_t(leading) * fields, padding);
            std::vector<Complex> coefficients(size_t(projectors) * fields), initial(input.size(), padding);
            for (int field = 0; field < fields; ++field)
            for (int g = 0; g < rows; ++g)
            {
                input[size_t(field) * leading + g] = Complex(uniform(generator), uniform(generator)) / std::sqrt(double(rows));
                initial[size_t(field) * leading + g] = Complex(uniform(generator), uniform(generator));
            }
            for (auto& value : coefficients) { value = Complex(uniform(generator), uniform(generator)); }
            Buffer device_input(input.size()), device_coefficients(coefficients.size());
            Buffer projected(size_t(projectors) * fields), output(initial.size());
            device_input.upload(input); device_coefficients.upload(coefficients); output.upload(initial);
            require(folded->project(fields, leading, device_input.data, projected.data), "Unexpected contraction fallback");
            require(folded->add(fields, leading, device_coefficients.data, output.data), "Unexpected expansion fallback");
            std::vector<Complex> expected_projection(projected.size), expected_output(initial);
            hsolver::gemm_op<Complex, CPU>()(nullptr, 'C', 'N', projectors, fields, rows,
                &one, matrix.data(), projector_leading, input.data(), leading, &zero, expected_projection.data(), projectors);
            hsolver::gemm_op<Complex, CPU>()(nullptr, 'N', 'T', rows, fields, projectors,
                &one, matrix.data(), projector_leading, coefficients.data(), fields, &one, expected_output.data(), leading);
            const auto actual_projection = projected.download(), actual_output = output.download();
            for (size_t index = 0; index < actual_projection.size(); ++index)
            {
                const double difference = std::abs(actual_projection[index] - expected_projection[index]);
                require(std::isfinite(difference), "Nonfinite folded projection");
                projection_error = std::max(projection_error, difference);
            }
            for (size_t index = 0; index < actual_output.size(); ++index)
            {
                const double difference = std::abs(actual_output[index] - expected_output[index]);
                require(std::isfinite(difference), "Nonfinite folded expansion");
                addition_error = std::max(addition_error, difference);
                if (index % leading >= size_t(rows)) { require(actual_output[index] == padding, "Output padding changed"); }
            }
            require(device_matrix.download() == matrix && device_input.download() == input
                    && device_coefficients.download() == coefficients, "A folded input changed");
        }
    }
    std::cout << "scenario=" << scenario << " projection_error=" << projection_error
              << " addition_error=" << addition_error << std::endl;
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    hsolver::createGpuBlasHandle();
    double projection_error = 0, addition_error = 0;
    for (int scenario : {0, 1, 2, 3, 4}) { run(scenario, projection_error, addition_error); }
    hsolver::destoryBLAShandle();
    MPI_Finalize();
    return projection_error < 1e-11 && addition_error < 1e-11 ? 0 : 2;
}
