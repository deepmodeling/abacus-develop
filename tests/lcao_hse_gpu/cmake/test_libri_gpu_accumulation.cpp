#include <RI/global/MPI_Wrapper.h>
#include <RI/ri/LRI.h>
#include <array>
#include <cmath>
#include <complex>
#include <iostream>
#include <map>
#include <mpi.h>
#include <stdexcept>
#include <vector>

namespace
{
template <typename T>
RI::Tensor<T> make_tensor(const RI::Shape_Vector& shape, const T imaginary)
{
    RI::Tensor<T> tensor(shape);
    for (std::size_t i = 0; i < tensor.data->size(); ++i)
    {
        (*tensor.data)[i] = T(0.5 + 0.07 * i) + imaginary * T(0.1 + 0.03 * i);
    }
    return tensor;
}

template <typename T>
void check_accumulation(const char* scalar_name, const T imaginary)
{
    using LRI = RI::LRI<int, int, 3, T>;
    using Tensors = std::map<int, std::map<typename LRI::TAC, RI::Tensor<T>>>;
    const typename LRI::TAC key = {0, {0, 0, 0}};
    const std::map<int, std::array<double, 3>> atoms = {{0, {0.0, 0.0, 0.0}}};
    const std::array<std::array<double, 3>, 3> lattice = {{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}};
    const std::vector<RI::Label::ab_ab> labels(RI::Label::array_ab_ab.begin(), RI::Label::array_ab_ab.end());
    LRI lri;
    lri.set_parallel(MPI_COMM_WORLD, atoms, lattice, {1, 1, 1}, labels);
    for (const auto label: RI::Label::array_ab)
    {
        const bool coefficient = label == RI::Label::ab::a || label == RI::Label::ab::b;
        const RI::Shape_Vector shape = coefficient ? RI::Shape_Vector{2, 2, 2} : RI::Shape_Vector{2, 2};
        Tensors tensors;
        tensors[0][key] = make_tensor(shape, imaginary);
        lri.set_tensors_map2(tensors, {label});
    }

    const RI::Tensor<T> seed = make_tensor<T>({2, 2}, imaginary);
    for (const auto label: labels)
    {
        Tensors reference;
        lri.cal_mode = RI::LRI_Cal_Mode::CPU;
        lri.cal_loop3({label}, reference, 1.0);
        const RI::Tensor<T> contribution = reference.at(0).at(key).copy();
        const double norm = contribution.norm(2);
        if (!std::isfinite(norm) || norm < 1e-8)
        {
            throw std::runtime_error("Empty or invalid CPU reference contraction");
        }

        for (const auto mode: {RI::LRI_Cal_Mode::CPU, RI::LRI_Cal_Mode::GPU})
        {
            lri.cal_mode = mode;
            Tensors result;
            // Tensor copies share storage unless copy() is used explicitly.
            result[0][key] = seed.copy();
            double accumulated_factor = 0.0;
            for (const double factor: {-1.0, 1.0, 1.0, -1.0})
            {
                lri.cal_loop3({label}, result, factor);
                accumulated_factor += factor;
                const RI::Tensor<T>& actual = result.at(0).at(key);
                if (result.size() != 1 || result.at(0).size() != 1 || actual.shape.size() != 2 || actual.shape[0] != 2
                    || actual.shape[1] != 2)
                {
                    throw std::runtime_error("Unexpected contraction output shape");
                }
                for (std::size_t i = 0; i < seed.data->size(); ++i)
                {
                    const T expected = (*seed.data)[i] + accumulated_factor * (*contribution.data)[i];
                    const double error = std::abs((*actual.data)[i] - expected);
                    const double tolerance = 1e-11 * (1.0 + std::abs(expected));
                    if (!std::isfinite(error) || error > tolerance)
                    {
                        std::cerr << scalar_name << " label=" << static_cast<int>(label)
                                  << " mode=" << static_cast<int>(mode) << " factor=" << factor << " element=" << i
                                  << " error=" << error << '\n';
                        throw std::runtime_error("Signed accumulation differs from the CPU reference");
                    }
                }
            }
        }
    }
    std::cout << scalar_name << ": all 18 contractions preserve signed accumulation\n";
}
} // namespace

int main(int argc, char** argv)
{
    // This standalone dependency test has no ABACUS runtime to initialize MPI.
    int provided = 0;
    if (MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided) != MPI_SUCCESS)
    {
        return 1;
    }
    int status = 0;
    try
    {
        const int ranks = RI::MPI_Wrapper::mpi_get_size(MPI_COMM_WORLD);
        if (ranks != 1 || provided < MPI_THREAD_MULTIPLE)
        {
            throw std::runtime_error("This test requires one MPI rank and MPI_THREAD_MULTIPLE");
        }
        check_accumulation<double>("real", 0.0);
        const std::complex<double> imaginary(0.0, 1.0);
        check_accumulation<std::complex<double>>("complex", imaginary);
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        status = 1;
    }
    MPI_Finalize();
    return status;
}
