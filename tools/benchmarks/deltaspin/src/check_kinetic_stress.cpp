#include "module_hamilt_pw/hamilt_pwdft/stress_func.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include <ATen/core/tensor.h>
#include <ATen/core/tensor_map.h>
#include <mpi.h>
#include <omp.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int ranks = 0, rank = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    POOL_WORLD = MPI_COMM_WORLD;
    GlobalV::NPROC_IN_POOL = ranks;
    GlobalV::RANK_IN_POOL = rank;
    GlobalV::NBANDS = 7;
    ModuleSymmetry::Symmetry::symm_flag = -1;
    GlobalC::ucell.ntype = 1;
    GlobalC::ucell.nat = 1;
    GlobalC::ucell.atoms = new Atom[1];
    GlobalC::ucell.atoms[0].na = 1;
    GlobalC::ucell.atoms[0].nw = 1;
    fftw_init_threads();
    fftw_plan_with_nthreads(omp_get_max_threads());
    double max_error = 0.0, phase_error = 0.0;
    for (int npol : {1, 2})
    {
        GlobalV::NPOL = npol;
        GlobalC::ucell.set_iat2iwt(npol);
        for (double height : {0.0, 0.7})
        {
            ModulePW::PW_Basis_K basis("gpu", "double");
            basis.initmpi(ranks, rank, MPI_COMM_WORLD);
            basis.initgrids(8.0, ModuleBase::Matrix3(1, 0.1, 0, 0, 1, 0.2, 0, 0, 1), 18, 20, 22);
            ModuleBase::Vector3<double> kvec[2] = {{0.25, 0.125, -0.1}, {-0.2, 0.3, 0.25}};
            basis.initparameters(false, 10.0, 2, kvec);
            basis.setuptransform();
            basis.collect_local_pw(1.3, height, 0.4);
            GlobalC::ucell.lat0 = basis.lat0;
            GlobalC::ucell.omega = basis.omega;
            K_Vectors kpoints;
            kpoints.set_nks(2);
            kpoints.ngk.assign(basis.npwk, basis.npwk + 2);
            ModuleBase::matrix weights(2, GlobalV::NBANDS);
            for (int ik = 0; ik < 2; ++ik)
            {
                for (int band = 0; band < GlobalV::NBANDS; ++band)
                {
                    weights(ik, band) = band % 3 == 0 ? 0.0 : 0.5 / (band + 1);
                }
            }
            psi::Psi<std::complex<double>> host(2, GlobalV::NBANDS, npol * basis.npwk_max, basis.npwk);
            std::mt19937 generator(97 + rank);
            std::uniform_real_distribution<double> uniform(-0.1, 0.1);
            for (size_t i = 0; i < host.size(); ++i) { host.get_pointer()[i] = {uniform(generator), uniform(generator)}; }
            ModuleBase::matrix original(3, 3);
            for (int phase = 0; phase < 2; ++phase)
            {
                if (phase)
                {
                    for (size_t i = 0; i < host.size(); ++i)
                    {
                        const auto value = host.get_pointer()[i];
                        host.get_pointer()[i] = {-value.imag(), value.real()};
                    }
                }
                psi::Psi<std::complex<double>, base_device::DEVICE_GPU> device(host);
                Stress_Func<double, base_device::DEVICE_CPU> reference;
                ModuleBase::matrix expected(3, 3), actual(3, 3);
                reference.stress_kin(expected, weights, &GlobalC::ucell.symm, &kpoints, &basis, &host);
                auto device_weights = ct::TensorMap(weights.c, ct::DataType::DT_DOUBLE, ct::DeviceType::CpuDevice,
                                                   {2, GlobalV::NBANDS}).to_device<ct::DEVICE_GPU>();
                const int max_blocks = hamilt::kinetic_stress_blocks_gpu(basis.npwk_max);
                ct::Tensor partial(ct::DataType::DT_DOUBLE, ct::DeviceType::GpuDevice, {max_blocks, 6});
                double sums[6] = {};
                for (int ik = 0; ik < 2; ++ik)
                {
                    hamilt::kinetic_stress_gpu(basis.npwk[ik], basis.npwk_max, GlobalV::NBANDS, npol, basis.tpiba,
                                               basis.erf_ecut, basis.erf_height, basis.erf_sigma,
                                               basis.get_gcar_data<double>() + 3 * ik * basis.npwk_max,
                                               basis.get_kvec_c_data<double>() + 3 * ik,
                                               device_weights.data<double>() + ik * GlobalV::NBANDS,
                                               &device(ik, 0, 0), partial.data<double>());
                    auto host_partial = partial.to_device<ct::DEVICE_CPU>();
                    for (int block = 0; block < hamilt::kinetic_stress_blocks_gpu(basis.npwk[ik]); ++block)
                    {
                        for (int component = 0; component < 6; ++component)
                        {
                            sums[component] += host_partial.data<double>()[6 * block + component];
                        }
                    }
                }
                int component = 0;
                for (int row = 0; row < 3; ++row)
                {
                    for (int column = 0; column <= row; ++column)
                    {
                        actual(row, column) = actual(column, row) = 2.0 * sums[component++] / basis.omega;
                    }
                }
                Parallel_Reduce::reduce_all(actual.c, 9);
                for (int index = 0; index < 9; ++index)
                {
                    max_error = std::max(max_error, std::abs(expected.c[index] - actual.c[index]));
                    if (phase) { phase_error = std::max(phase_error, std::abs(original.c[index] - actual.c[index])); }
                    else { original.c[index] = actual.c[index]; }
                }
            }
        }
    }
    if (rank == 0) { std::cout << "stress_error=" << max_error << " phase_invariance_error=" << phase_error << '\n'; }
    fftw_cleanup_threads();
    MPI_Finalize();
    return max_error < 1e-12 && phase_error < 1e-12 ? 0 : 2;
}
